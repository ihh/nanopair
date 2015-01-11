#ifndef SEQEVTPAIR_INCLUDED
#define SEQEVTPAIR_INCLUDED

#include "fast5events.h"
#include "kseqcontainer.h"
#include "xmlparser.h"

/* Parameters */

typedef struct Seq_event_pair_model {
  double *pMatchEmit;
  double pBeginDelete, pExtendDelete;
  int order, states;
  double *matchMean, *matchPrecision;
  double pStartEmit, startMean, startPrecision;
} Seq_event_pair_model;

Seq_event_pair_model* new_seq_event_pair_model (int order);
void delete_seq_event_pair_model (Seq_event_pair_model* model);

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml);
xmlChar* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model);

/* Data + precomputed log-likelihoods for DP */

typedef struct Seq_event_pair_data {
  /* parameters */
  Seq_event_pair_model* model;
  /* data */
  int seqlen, matrix_cells;
  char *seq;  /* not owned */
  Fast5_event_array* events;  /* not owned */
  int *state;  /* owned */
  /* precalculated emit & transition scores */
  long double *startEmitDensity, *matchEmitDensity, *matchEmitYes, *matchEmitNo;
  long double startEmitYes, startEmitNo;
  long double beginDeleteYes, beginDeleteNo;
  long double extendDeleteYes, extendDeleteNo;
} Seq_event_pair_data;

Seq_event_pair_data* new_seq_event_pair_data (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_data (Seq_event_pair_data* matrix);

void precalc_seq_event_pair_data (Seq_event_pair_data* data);

/* Seq_event_pair_index macro assumes seqlen & order are the same as in Seq_event_pair_fb_data */
#define Seq_event_pair_index(seqpos,n_event) (n_event*(seqlen-order+1) + seqpos - order)

/* Forward-backward matrix */

typedef struct Seq_event_pair_fb_matrix {
  /* data */
  Seq_event_pair_data* data;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *fwdStart, *fwdMatch, *fwdDelete, fwdEnd;
  long double *backStart, *backMatch, *backDelete;
} Seq_event_pair_fb_matrix;

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

/* DP helper functions */

double log_gaussian_density (double x, double mean, double precision, double log_precision);
double log_event_density (Fast5_event* event, double mean, double precision, double log_precision);

long double log_sum_exp (long double a, long double b);  /* returns log(exp(a) + exp(b)) */

/* Expected counts */

typedef struct Seq_event_pair_counts {
  long double nStartEmitYes, nStartEmitNo;
  long double *nMatchEmitYes, *nMatchEmitNo;
  long double nBeginDeleteYes, nBeginDeleteNo;
  long double nExtendDeleteYes, nExtendDeleteNo;
  int order, states;
  long double *matchMoment0, *matchMoment1, *matchMoment2;
  long double startMoment0, startMoment1, startMoment2;
} Seq_event_pair_counts;

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model);
void delete_seq_event_pair_counts (Seq_event_pair_counts* counts);

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts);
void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);

/* Single Baum-Welch iteration */

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* matrix, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Vector* event_arrays, double minimum_fractional_log_likelihood_increase);

/* Alignment */

typedef struct Seq_event_pair_alignment {
  Fast5_event_array *events;  /* not owned */
  char *seqname, *seq;  /* not owned */
  int *n_events_at_pos;  /* owned */
  int first_seqpos, seqpos_len, n_first_event;
} Seq_event_pair_alignment;

Seq_event_pair_alignment* new_seq_event_pair_alignment (Fast5_event_array *events, char *seq, int first_seqpos, int seqpos_len);
void delete_seq_event_pair_alignment (Seq_event_pair_alignment* align);

void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, FILE* out);

/* Viterbi matrix */

typedef struct Seq_event_pair_viterbi_matrix {
  /* data */
  Seq_event_pair_data* data;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *vitStart, *vitMatch, *vitDelete, vitEnd;
} Seq_event_pair_viterbi_matrix;

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);

Seq_event_pair_alignment* fill_seq_event_pair_viterbi_matrix_and_get_trace (Seq_event_pair_viterbi_matrix* matrix);

#endif /* SEQEVTPAIR_INCLUDED */
