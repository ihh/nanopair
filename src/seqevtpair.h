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
  double pStartEmit;
  double pNullEmit, nullMean, nullPrecision;
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
  long double *nullEmitDensity, *matchEmitDensity, *matchEmitYes, *matchEmitNo;
  long double nullEmitYes, nullEmitNo;
  long double startEmitYes, startEmitNo;
  long double beginDeleteYes, beginDeleteNo;
  long double extendDeleteYes, extendDeleteNo;
  /* null model log-likelihood */
  long double nullModel;
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
  long double nNullEmitYes, nNullEmitNo;
  long double *nMatchEmitYes, *nMatchEmitNo;
  long double nBeginDeleteYes, nBeginDeleteNo;
  long double nExtendDeleteYes, nExtendDeleteNo;
  int order, states;
  long double *matchMoment0, *matchMoment1, *matchMoment2;
  long double nullMoment0, nullMoment1, nullMoment2;
} Seq_event_pair_counts;

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model);
void delete_seq_event_pair_counts (Seq_event_pair_counts* counts);

Seq_event_pair_counts* new_seq_event_pair_counts_minimal_prior (Seq_event_pair_model* model);  /* Laplace +1's */

void reset_seq_event_null_counts (Seq_event_pair_counts* counts);
void reset_seq_event_pair_counts (Seq_event_pair_counts* counts);

void add_weighted_seq_event_pair_counts (Seq_event_pair_counts* counts, Seq_event_pair_counts* inc, long double weight);

void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);  /* NB does not touch delete counts */
void inc_seq_event_null_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);

void optimize_seq_event_null_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);
void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

/* Single Baum-Welch iteration */

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
long double inc_seq_event_pair_counts_via_fb (Seq_event_pair_model* model, Seq_event_pair_counts* counts, int seqlen, char *seq, Fast5_event_array* events);  /* wraps creation & destruction of FB matrix, returns log-likelihood ratio (forward/null) */

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays);

/* Alignment */

typedef struct Seq_event_pair_alignment {
  Fast5_event_array *events;  /* not owned */
  const char *seqname, *seq;  /* not owned */
  int *events_at_pos;  /* owned. NB: events_at_pos[n] refers to seq[start_seqpos + n] */
  /* end_seqpos = index of last aligned character + 1
     thus, number of aligned bases = end_seqpos - start_seqpos */
  int seqlen, start_seqpos, end_seqpos, start_n_event;
  long double log_likelihood_ratio;
} Seq_event_pair_alignment;

Seq_event_pair_alignment* new_seq_event_pair_alignment (Fast5_event_array *events, char *seq, int seqlen);
void delete_seq_event_pair_alignment (Seq_event_pair_alignment* align);

void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, int strand, FILE* out);

/* Viterbi matrix */

typedef struct Seq_event_pair_viterbi_matrix {
  /* data */
  Seq_event_pair_data* data;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *vitStart, *vitMatch, *vitDelete, vitEnd;
} Seq_event_pair_viterbi_matrix;

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);

void fill_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);
Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix);

/* wrappers */

void print_seq_evt_pair_alignments_as_gff_cigar (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold);

#endif /* SEQEVTPAIR_INCLUDED */
