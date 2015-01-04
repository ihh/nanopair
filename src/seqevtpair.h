#ifndef SEQEVTPAIR_INCLUDED
#define SEQEVTPAIR_INCLUDED

#include "fast5events.h"
#include "kseqcontainer.h"
#include "xmlparser.h"

/* Parameters */

typedef struct Seq_event_pair_model {
  double *pEmitCurrent;
  double pBeginDelete, pExtendDelete;
  int order, states;
  double *currentMean, *currentPrecision;
  double pStartEmitCurrent, startCurrentMean, startCurrentPrecision;
} Seq_event_pair_model;

Seq_event_pair_model* new_seq_event_pair_model (int order);
void delete_seq_event_pair_model (Seq_event_pair_model* model);

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml);
xmlChar* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model);

/* Forward-backward matrix */

typedef struct Seq_event_pair_fb_matrix {
  /* parameters */
  Seq_event_pair_model* model;
  /* data */
  int seqlen;
  int* dsq;
  Fast5_event_array* events;
  /* dynamic programming matrices */
  long double *fwdStart, *fwdMatch, *fwdDelete;
  long double *backStart, *backMatch, *backDelete;
  long double fwdLikelihood;
} Seq_event_pair_fb_matrix;

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, int* dsq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

/* Expected counts */

typedef struct Seq_event_pair_counts {
  long double *nEmitCurrent_yes, *nEmitCurrent_no;
  long double nBeginDelete_yes, nBeginDelete_no;
  long double nExtendDelete_yes, nExtendDelete_no;
  int order, states;
  long double *currentMoment0, *currentMoment1, *currentMoment2;
  long double startCurrentMoment0, startCurrentMoment1, startCurrentMoment2;
} Seq_event_pair_counts;

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model);
void delete_seq_event_pair_counts (Seq_event_pair_counts* counts);

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts);
void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);

/* Single Baum-Welch iteration */

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
void update_seq_event_pair_model (Seq_event_pair_model* matrix, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Fast5_event_array* events, double minimum_fractional_log_likelihood_increase);

#endif /* SEQEVTPAIR_INCLUDED */
