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

  matrix->nullEmitDensity = SafeMalloc ((n_events + 1) * sizeof(long double));

  matrix->matchEmitDensity = SafeMalloc (matrix->matrix_cells * sizeof(long double));
  matrix->matchEmitYes = SafeMalloc (states * sizeof(long double));
  matrix->matchEmitNo = SafeMalloc (states * sizeof(long double));

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
  matrix->nullEmitYes = log (model->pNullEmit);
  matrix->nullEmitNo = log (1. - model->pNullEmit);

  matrix->startEmitYes = log (model->pStartEmit);
  matrix->startEmitNo = log (1. - model->pStartEmit);

  matrix->beginDeleteYes = log (model->pBeginDelete);
  matrix->beginDeleteNo = log (1. - model->pBeginDelete);

  matrix->extendDeleteYes = log (model->pExtendDelete);
  matrix->extendDeleteNo = log (1. - model->pExtendDelete);

  double logNullPrecision = log (model->nullPrecision);

  /* calculate emit densities & state-dependent transition probabilities */
  matrix->nullEmitDensity[0] = -INFINITY;
  matrix->nullModel = matrix->nullEmitNo;
  for (int n_event = 0; n_event < n_events; ++n_event) {
    Fast5_event* event = &events->event[n_event];
    double loglike = log_event_density (event,
					model->nullMean,
					model->nullPrecision,
					logNullPrecision);
    matrix->nullEmitDensity[n_event + 1] = loglike;
    matrix->nullModel += loglike + matrix->nullEmitYes * event->ticks;
  }

  for (int state = 0; state < states; ++state) {
    matrix->matchEmitYes[state] = log (model->pMatchEmit[state]);
    matrix->matchEmitNo[state] = log (1. - model->pMatchEmit[state]);

    double mean = model->matchMean[state];
    double precision = model->matchPrecision[state];
    double logPrecision = log (precision);

    matrix->matchEmitDensity[Basecall_index(state,0)] = -INFINITY;
    for (int n_event = 1; n_event <= n_events; ++n_event) {
      Fast5_event* event = &events->event[n_event - 1];
      matrix->matchEmitDensity[Basecall_index(state,n_event)]
	= log_event_density (event,
			     mean,
			     precision,
			     logPrecision);
    }
  }

  return matrix;
}

void delete_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix) {
  SafeFree (matrix->nullEmitDensity);
  SafeFree (matrix->matchEmitDensity);
  SafeFree (matrix->matchEmitYes);
  SafeFree (matrix->matchEmitNo);
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
    Fast5_event* event = &matrix->events->event[n_event - 1];
    matrix->vitStart[n_event]
      = matrix->vitStart[n_event - 1]
      + matrix->startEmitYes * event->ticks
      + matrix->nullEmitDensity[n_event];  /* Start -> Start (output) */
  }
  
  /* TODO: write me */
}

char* get_basecall_viterbi_matrix_traceback (Basecall_viterbi_matrix* matrix) {
  /* TODO: write me */
  return NULL;
}
