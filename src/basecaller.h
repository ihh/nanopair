#ifndef BASECALLER_INCLUDED
#define BASECALLER_INCLUDED

#include "seqevtpair.h"

typedef struct Basecall_viterbi_matrix {
  /* model, data, dimensions */
  Seq_event_pair_model *model;
  Fast5_event_array *events;
  unsigned long matrix_cells;
  /* precomputed scores */
  long double *nullEmitDensity, *matchEmitDensity, *matchEmitYes, *matchEmitNo;
  long double nullEmitYes, nullEmitNo;
  long double startEmitYes, startEmitNo;
  long double beginDeleteYes, beginDeleteNo;
  long double extendDeleteYes, extendDeleteNo;
  long double *logKmerProb, *logKmerConditionalProb;
  /* matrix */
  long double *vitStart, *vitMatch, vitTotal;
  long double nullModel;
} Basecall_viterbi_matrix;

/* Basecall_index macros assume local variable 'states' is equal to matrix->model->states */
#define Unsafe_basecall_index(state,n_event) (((unsigned long)n_event) * ((unsigned long)states) + state)
#ifdef SEQEVTPAIR_DEBUG
#define Basecall_index(state,n_event) ((n_event < 0 || n_event > n_events || state < 0 || state >= states) ? (Abort("Index (%d,%d) out of bounds for dimensions (%d,%d)",state,n_event,states,n_events), 0) : Unsafe_basecall_index(state,n_event))
#else /* SEQEVTPAIR_DEBUG */
#define Basecall_index(seqpos,n_event) Unsafe_basecall_index(seqpos,n_event)
#endif /* SEQEVTPAIR_DEBUG */
unsigned long basecall_index_wrapper (Basecall_viterbi_matrix* matrix, int state, int n_event);

Basecall_viterbi_matrix* new_basecall_viterbi_matrix (Seq_event_pair_model* model, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix);

void fill_basecall_viterbi_matrix (Basecall_viterbi_matrix* matrix);
char* get_basecall_viterbi_matrix_traceback (Basecall_viterbi_matrix* matrix);

#endif /* BASECALLER_INCLUDED */
