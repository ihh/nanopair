#include <math.h>
#include "basecaller.h"

unsigned long basecall_index_wrapper (Basecall_viterbi_matrix* matrix, int state, int n_event) {
  int states, n_events;
  n_events = matrix->events->n_events;
  states = matrix->model->states;
  return Basecall_index(state,n_event);
}

Basecall_viterbi_matrix* new_basecall_viterbi_matrix (Seq_event_pair_model* model, Fast5_event_array* events) {
  Basecall_viterbi_matrix* matrix = SafeMalloc (sizeof (Basecall_viterbi_matrix));
  matrix->model = model;
  matrix->events = events;

  int n_events = events->n_events;
  int states = model->states;
  matrix->matrix_cells = states * (n_events + 1);

  matrix->matchEventLogLike = SafeMalloc (matrix->matrix_cells * sizeof(long double));
  matrix->noDelete = SafeMalloc (states * sizeof(long double));
  matrix->shortDelete = SafeMalloc (states * sizeof(long double));
  matrix->matchEventYes = SafeMalloc (states * sizeof(long double));
  matrix->matchEventNo = SafeMalloc (states * sizeof(long double));

  matrix->logKmerProb = SafeMalloc (states * sizeof(double));
  matrix->logKmerConditionalProb = SafeMalloc (states * sizeof(double));
  for (int prefix = 0; prefix < states; prefix += AlphabetSize) {
    double norm = 0;
    for (int state = prefix; state < prefix + AlphabetSize; ++state)
      norm += model->kmerProb[state];
    for (int state = prefix; state < prefix + AlphabetSize; ++state) {
      matrix->logKmerProb[state] = log (model->kmerProb[state]);
      matrix->logKmerConditionalProb[state] = log (model->kmerProb[state] / norm);
    }
  }

  matrix->vitStart = SafeMalloc ((n_events + 1) * sizeof(double));
  matrix->vitMatch = SafeMalloc (matrix->matrix_cells * sizeof(double));

  return matrix;
}

void delete_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix) {
  SafeFree (matrix->matchEventLogLike);
  SafeFree (matrix->matchEventYes);
  SafeFree (matrix->matchEventNo);
  SafeFree (matrix->shortDelete);
  SafeFree (matrix->noDelete);
  SafeFree (matrix->logKmerProb);
  SafeFree (matrix->logKmerConditionalProb);
  SafeFree (matrix->vitStart);
  SafeFree (matrix->vitMatch);
  SafeFree (matrix);
}

void fill_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix) {
  Seq_event_pair_model* model = matrix->model;
  Fast5_event_array* events = matrix->events;
  int n_events = events->n_events;
  int states = model->states;
  
  /* calculate logs of transition probabilities & emit precisions */

  matrix->longDelete = log(model->emitProb) + log(model->pBeginDelete) + log(model->pExtendDelete);
  matrix->emitNo = log(1. - model->emitProb);
  
  /* calculate emit densities & state-dependent transition probabilities */
  for (int state = 0; state < states; ++state) {
    matrix->matchEventYes[state] = log (model->pMatchEvent[state]);
    matrix->matchEventNo[state] = log (1. - model->pMatchEvent[state]);

    matrix->noDelete[state] = log(model->emitProb) + log(1. - model->pBeginDelete) + log(1. - model->pMatchSkip[state]);
    matrix->shortDelete[state] = log(model->emitProb) + log(model->pBeginDelete * (1. - model->pExtendDelete)
							    + (1. - model->pBeginDelete) * model->pMatchSkip[state]);

    double mean = model->matchMean[state];
    double precision = model->matchPrecision[state];
    double logPrecision = log (precision);
    double logTick = log (model->pMatchTick[state]);
    double logNoTick = log (1. - model->pMatchTick[state]);

    matrix->matchEventLogLike[Basecall_index(state,0)] = -INFINITY;
    for (int n_event = 1; n_event <= n_events; ++n_event) {
      Fast5_event* event = &events->event[n_event - 1];
      matrix->matchEventLogLike[Basecall_index(state,n_event)]
	= log_event_density (event,
			     mean,
			     precision,
			     logPrecision,
			     logTick,
			     logNoTick);
    }
  }

  matrix->vitStart[0] = 0.;
  for (int state = 0; state < states; ++state)
    matrix->vitMatch[Basecall_index(state,0)] = -INFINITY;
  
  int firstBaseMultiplier = pow (AlphabetSize, model->order - 1);
  long double longDel = matrix->longDelete;
  for (int n_event = 1; n_event <= n_events; ++n_event) {
    long double st = -INFINITY;
    for (int state = 0; state < states; ++state) {
      long double m = max_func (matrix->vitStart[n_event - 1]
				+ matrix->logKmerProb[state],   /* Start -> Emit (x^N y) */
				matrix->vitMatch[Basecall_index(state,n_event-1)]
				+ matrix->matchEventYes[state]);  /* Emit -> Emit (y) */
      int prefix = state / AlphabetSize;
      long double logCondProb = matrix->logKmerConditionalProb[state];
      long double noDel = matrix->noDelete[state];
      for (int prevBase = 0; prevBase < AlphabetSize; ++prevBase) {
	int prevState = prefix + prevBase * firstBaseMultiplier;
	m = max_func (m,
		      matrix->vitMatch[Basecall_index(prevState,n_event-1)]
		      + matrix->matchEventNo[prevState]
		      + noDel
		      + logCondProb);   /* Emit -> Emit (xy) */
	long double logPrevCondProb = matrix->logKmerConditionalProb[prevState];
	long double shortDel = matrix->shortDelete[prevState];
	int prePrefix = prevState / AlphabetSize;
	for (int prevPrevBase = 0; prevPrevBase < AlphabetSize; ++prevPrevBase) {
	  int prevPrevState = prePrefix + prevPrevBase * firstBaseMultiplier;
	  m = max_func (m,
			matrix->vitMatch[Basecall_index(prevPrevState,n_event-1)]
			+ matrix->matchEventNo[prevPrevState]
			+ logPrevCondProb
			+ shortDel
			+ logCondProb);   /* Emit -> Emit (xxy) */
	}
      }
      int idx = Basecall_index(state,n_event);
      m += matrix->matchEventLogLike[idx];
      matrix->vitMatch[idx] = m;
      st = max_func (st, m + longDel);  /* Emit -> Start */
    }
    matrix->vitStart[n_event] = st;
  }

  long double end = -INFINITY;
  for (int state = 0; state < states; ++state)
    end = max_func (end,
		    matrix->vitMatch[Basecall_index(state,n_events)]
		    + matrix->emitNo);
  matrix->vitResult = end;
}

char* get_basecall_viterbi_matrix_traceback (Basecall_viterbi_matrix* matrix) {
  /* TODO: write me */
  return NULL;
}
