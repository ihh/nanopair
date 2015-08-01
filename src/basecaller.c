#include <math.h>
#include "basecaller.h"

/* helper for finding index of max item in a set indexed by 2-tuples (for Viterbi traceback) */
void update_max2 (long double *current_max,
		  int* current_max_idx1,
		  int* current_max_idx2,
		  long double candidate_max,
		  int candidate_max_idx1,
		  int candidate_max_idx2);


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
  matrix->matchEventYes = SafeMalloc (states * sizeof(long double));
  matrix->matchEventNo = SafeMalloc (states * sizeof(long double));
  matrix->matchSkipYes = SafeMalloc (states * sizeof(long double));
  matrix->matchSkipNo = SafeMalloc (states * sizeof(long double));
  
  matrix->logKmerProb = SafeMalloc (states * sizeof(long double));
  matrix->logKmerConditionalProb = SafeMalloc (states * sizeof(long double));
  for (int prefix = 0; prefix < states; prefix += AlphabetSize) {
    double norm = 0;
    for (int state = prefix; state < prefix + AlphabetSize; ++state)
      norm += model->kmerProb[state];
    for (int state = prefix; state < prefix + AlphabetSize; ++state) {
      matrix->logKmerProb[state] = log (model->kmerProb[state]);
      matrix->logKmerConditionalProb[state] = log (model->kmerProb[state] / norm);
    }
  }

  matrix->vitStart = SafeMalloc ((n_events + 1) * sizeof(long double));
  matrix->vitMatch = SafeMalloc (matrix->matrix_cells * sizeof(long double));

  return matrix;
}

void delete_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix) {
  SafeFree (matrix->matchEventLogLike);
  SafeFree (matrix->matchEventYes);
  SafeFree (matrix->matchEventNo);
  SafeFree (matrix->matchSkipYes);
  SafeFree (matrix->matchSkipNo);
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
  
  Logger* logger = model->logger;

  /* calculate logs of transition probabilities & emit precisions */

  matrix->beginDeleteYes = log(model->pBeginDelete);
  matrix->beginDeleteNo = log(1. - model->pBeginDelete);
  matrix->extendDeleteYes = log(model->pExtendDelete);
  matrix->extendDeleteNo = log(1. - model->pExtendDelete);
  matrix->emitYes = log(model->emitProb);
  matrix->emitNo = log(1. - model->emitProb);

  matrix->longDelete = matrix->emitYes * model->order + matrix->beginDeleteYes + matrix->extendDeleteYes * (model->order - 1) + matrix->extendDeleteNo;

  long double longDel = matrix->longDelete;

  /* calculate emit densities & state-dependent transition probabilities */
  for (int state = 0; state < states; ++state) {
    matrix->matchEventYes[state] = log (model->pMatchEvent[state]);
    matrix->matchEventNo[state] = log (1. - model->pMatchEvent[state]);

    matrix->matchSkipYes[state] = log (model->pMatchSkip[state]);
    matrix->matchSkipNo[state] = log (1. - model->pMatchSkip[state]);

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

  if (LogThisAt(2))
    init_progress ("Viterbi basecaller");

  for (int n_event = 1; n_event <= n_events; ++n_event) {

    if (LogThisAt(2))
      log_progress (n_event / (double) n_events, "event %d/%d", n_event, n_events);

    long double st = -INFINITY;
    for (int state = 0; state < states; ++state) {
      long double m = max_func (matrix->vitStart[n_event - 1]
				+ matrix->logKmerProb[state],   /* Start -> Emit (x^N y) */
				matrix->vitMatch[Basecall_index(state,n_event-1)]
				+ matrix->matchEventYes[state]);  /* Emit -> Emit (y) */
      int prefix = state / AlphabetSize;
      long double logCondProb = matrix->logKmerConditionalProb[state];
      for (int prevBase = 0; prevBase < AlphabetSize; ++prevBase) {
	int prevState = prefix + prevBase * firstBaseMultiplier;
	m = max_func (m,
		      matrix->vitMatch[Basecall_index(prevState,n_event-1)]
		      + matrix->matchEventNo[prevState]
		      + matrix->emitYes
		      + matrix->beginDeleteNo
		      + matrix->matchSkipNo[state]
		      + logCondProb);   /* Emit -> Emit (xy) */
	long double logPrevCondProb = matrix->logKmerConditionalProb[prevState];
	int prePrefix = prevState / AlphabetSize;
	for (int prevPrevBase = 0; prevPrevBase < AlphabetSize; ++prevPrevBase) {
	  int prevPrevState = prePrefix + prevPrevBase * firstBaseMultiplier;
	  m = max_func (m,
			matrix->vitMatch[Basecall_index(prevPrevState,n_event-1)]
			+ matrix->matchEventNo[prevPrevState]
			+ matrix->emitYes*2
			+ matrix->beginDeleteNo
			+ matrix->matchSkipYes[prevState]
			+ matrix->matchSkipNo[state]
			+ logPrevCondProb
			+ logCondProb);   /* Emit -> Emit (xxy) */
	}
      }
      unsigned long idx = Basecall_index(state,n_event);
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

  if (LogThisAt(2))
    fprintf (stderr, "Viterbi log-likelihood is %Lg\n", matrix->vitResult);
}

void update_max2 (long double *current_max,
		  int* current_max_idx1,
		  int* current_max_idx2,
		  long double candidate_max,
		  int candidate_max_idx1,
		  int candidate_max_idx2) {
  if (candidate_max > *current_max) {
    *current_max = candidate_max;
    *current_max_idx1 = candidate_max_idx1;
    if (current_max_idx2)
      *current_max_idx2 = candidate_max_idx2;
  }
}

char* get_basecall_viterbi_matrix_traceback (Basecall_viterbi_matrix* matrix) {
  const int Start = -1;
  const int NoState = -2;
  enum { None, StartEmit, EmitEmit0, EmitEmit1, EmitEmit2, EmitStart } bestTrans;
  List *emittedBase;

  Seq_event_pair_model* model = matrix->model;
  Fast5_event_array* events = matrix->events;
  int n_events = events->n_events;
  int states = model->states;

  int firstBaseMultiplier = pow (AlphabetSize, model->order - 1);
  long double longDel = matrix->longDelete;

  int state, bestPrevState, endState = -1;
  long double loglike = -INFINITY;
  for (state = 0; state < states; ++state)
    update_max2 (&loglike,
		 &endState, NULL,
		 matrix->vitMatch[Basecall_index(state,n_events)]
		 + matrix->emitNo,  /* Match -> End */
		 state, -1);
  Assert (endState >= 0, "Traceback failed");

  int n_event = n_events;
  emittedBase = newList (IntCopy, IntDelete, IntPrint);

  state = endState;
  
  while (!(n_event == 0 && state == Start)) {
    unsigned long idx = Basecall_index(state,n_event);
    loglike = -INFINITY;
    bestTrans = None;
    bestPrevState = NoState;
    
    if (state >= 0) {
      update_max2 (&loglike,
		   (int*) &bestTrans, &bestPrevState,
		   matrix->vitStart[n_event - 1]
		   + matrix->logKmerProb[state],
		   StartEmit, Start);  /* Start -> Emit (x^N y) */

      update_max2 (&loglike,
		   (int*) &bestTrans, &bestPrevState,
		   matrix->vitMatch[Basecall_index(state,n_event-1)]
		   + matrix->matchEventYes[state],
		   EmitEmit0, state);  /* Emit -> Emit (y) */

      int prefix = state / AlphabetSize;
      long double logCondProb = matrix->logKmerConditionalProb[state];
      for (int prevBase = 0; prevBase < AlphabetSize; ++prevBase) {
	int prevState = prefix + prevBase * firstBaseMultiplier;
	update_max2 (&loglike,
		     (int*) &bestTrans, &bestPrevState,
		     matrix->vitMatch[Basecall_index(prevState,n_event-1)]
		     + matrix->matchEventNo[prevState]
		     + matrix->emitYes
		     + matrix->beginDeleteNo
		     + matrix->matchSkipNo[state]
		     + logCondProb,
		     EmitEmit1, prevState);   /* Emit -> Emit (xy) */

	long double logPrevCondProb = matrix->logKmerConditionalProb[prevState];
	int prePrefix = prevState / AlphabetSize;
	for (int prevPrevBase = 0; prevPrevBase < AlphabetSize; ++prevPrevBase) {
	  int prevPrevState = prePrefix + prevPrevBase * firstBaseMultiplier;
	  update_max2 (&loglike,
		       (int*) &bestTrans, &bestPrevState,
		       matrix->vitMatch[Basecall_index(prevPrevState,n_event-1)]
		       + matrix->matchEventNo[prevPrevState]
		       + matrix->emitYes*2
		       + matrix->beginDeleteNo
		       + matrix->matchSkipYes[prevState]
		       + matrix->matchSkipNo[state]
		       + logPrevCondProb
		       + logCondProb,
		       EmitEmit2, prevPrevState);   /* Emit -> Emit (xxy) */
	}
      }
      Assert (loglike + matrix->matchEventLogLike[idx] == matrix->vitMatch[idx], "Traceback error");

    } else if (state == Start) {
      for (int prevState = 0; prevState < states; ++prevState)
	update_max2 (&loglike,
		     (int*) &bestTrans, &bestPrevState,
		     matrix->vitMatch[Basecall_index(prevState,n_event)],
		     EmitStart, prevState);  /* Emit -> Start */
      Assert (loglike + longDel == matrix->vitStart[n_event], "Traceback error");

    } else {
      Abort ("Unknown traceback state");
    }

    switch (bestTrans) {

    case StartEmit:
      Assert (state >= 0, "oops");
      --n_event;
      int mul = 1;
      for (int pos = 0; pos < model->order; ++pos) {
	ListInsertBefore (emittedBase, emittedBase->head, IntNew((state / mul) % AlphabetSize));
	mul *= AlphabetSize;
      }
      break;

    case EmitEmit0:
      Assert (state >= 0, "oops");
      --n_event;
      break;

    case EmitEmit1:
      Assert (state >= 0, "oops");
      --n_event;
      ListInsertBefore (emittedBase, emittedBase->head, IntNew(state % AlphabetSize));
      break;

    case EmitEmit2:
      Assert (state >= 0, "oops");
      --n_event;
      ListInsertBefore (emittedBase, emittedBase->head, IntNew(state % AlphabetSize));
      ListInsertBefore (emittedBase, emittedBase->head, IntNew((state / AlphabetSize) % AlphabetSize));
      break;

    case EmitStart:
      Assert (state == Start, "oops");
      break;

    case None:
    default:
      Abort ("Unknown traceback transition");
      break;
    }

    state = bestPrevState;
  }
  
  char *seq = SafeMalloc ((ListSize(emittedBase) + 1) * sizeof(char));
  int pos = 0;
  for (ListNode* node = emittedBase->head; node; node = node->next)
    seq[pos++] = token2base (*(int*)node->value);
  seq[pos] = '\0';
  
  deleteList (emittedBase);
  
  return seq;
}

char* basecall_fast5_event_array (Seq_event_pair_model* model, Fast5_event_array* events) {
  Basecall_viterbi_matrix* matrix = new_basecall_viterbi_matrix (model, events);
  fill_basecall_viterbi_matrix (matrix);
  char *seq = get_basecall_viterbi_matrix_traceback (matrix);
  delete_basecall_viterbi_matrix (matrix);
  return seq;
}

