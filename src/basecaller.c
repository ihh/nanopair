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
  for (int prefix = 0; prefix < states; prefix += 4) {
    double norm = 0;
    for (int state = prefix; state < prefix + 4; ++state)
      norm += model->kmerProb[state];
    for (int state = prefix; state < prefix + 4; ++state) {
      matrix->logKmerProb[state] = log (model->kmerProb[state]);
      matrix->logKmerConditionalProb[state] = log (model->kmerProb[state] / norm);
    }
  }

  matrix->vitStart = SafeMalloc ((n_events + 1) * sizeof(double));
  matrix->vitMatch = SafeMalloc (matrix->matrix_cells * sizeof(double));

  /* calculate logs of transition probabilities & emit precisions */

  /* TODO: WRITE ME */
  
  /* calculate emit densities & state-dependent transition probabilities */
  for (int state = 0; state < states; ++state) {
    matrix->matchEventYes[state] = log (model->pMatchEvent[state]);
    matrix->matchEventNo[state] = log (1. - model->pMatchEvent[state]);

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
  int n_events = matrix->events->n_events;
  /*
  int states = matrix->model->states;
  */

  matrix->vitStart[0] = 0.;
  for (int n_event = 1; n_event <= n_events; ++n_event) {
    /*
    Fast5_event* event = &matrix->events->event[n_event - 1];
    */
    /* TODO: write me */
  }
}

char* get_basecall_viterbi_matrix_traceback (Basecall_viterbi_matrix* matrix) {
  /* TODO: write me */
  return NULL;
}
