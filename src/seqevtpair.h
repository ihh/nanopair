#ifndef SEQEVTPAIR_INCLUDED
#define SEQEVTPAIR_INCLUDED

#include "fast5events.h"
#include "kseqcontainer.h"
#include "xmlparser.h"

/* #define the following for various levels of debugging info/code
   (if defined at all) adds out-of-bounds guards to DP matrix accessors, logs progress messages
   (if defined to 10 or greater) dumps DP matrices to filename SEQEVTMATRIX_FILENAME
 */
#define SEQEVTPAIR_DEBUG 10
#define SEQEVTMATRIX_FILENAME "npmatrix"

/* state encoding functions */
int base2token (char base);
char token2base (int token);

void encode_state_identifier (int state, int order, char* state_id);
int decode_state_identifier (int order, char* state_id);

/* Parameters */

typedef struct Seq_event_pair_model {
  double *pMatchEmit;
  double pSkip, pBeginDelete, pExtendDelete;
  int order, states;
  double *matchMean, *matchPrecision;
  double pStartEmit;
  double pNullEmit, nullMean, nullPrecision;
  double *kmerProb;
} Seq_event_pair_model;

Seq_event_pair_model* new_seq_event_pair_model (int order);
void delete_seq_event_pair_model (Seq_event_pair_model* model);

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml);
xmlChar* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model);

xmlChar* make_squiggle_svg (Fast5_event_array *events, Seq_event_pair_model* model);


/* Data + precomputed log-likelihoods for DP */

typedef struct Seq_event_pair_data {
  /* parameters */
  Seq_event_pair_model* model;
  /* data */
  int seqlen;
  unsigned long matrix_cells;
  char *seq;  /* not owned */
  Fast5_event_array* events;  /* not owned */
  int *state;  /* owned */
  /* precalculated emit & transition scores */
  long double *nullEmitDensity, *matchEmitDensity, *matchEmitYes, *matchEmitNo;
  long double nullEmitYes, nullEmitNo;
  long double startEmitYes, startEmitNo;
  long double skipYes, skipNo;
  long double beginDeleteYes, beginDeleteNo;
  long double extendDeleteYes, extendDeleteNo;
  /* null model log-likelihood */
  long double nullModel;
} Seq_event_pair_data;

Seq_event_pair_data* new_seq_event_pair_data (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_data (Seq_event_pair_data* matrix);

void precalc_seq_event_pair_data (Seq_event_pair_data* data);

/* Seq_event_pair_index macro assumes seqlen & order are the same as in Seq_event_pair_fb_data */
#define Unsafe_seq_event_pair_index(seqpos,n_event) (((unsigned long)n_event) * ((unsigned long)(seqlen-order+1)) + (unsigned long)(seqpos - order))
#ifdef SEQEVTPAIR_DEBUG
#define Seq_event_pair_index(seqpos,n_event) ((n_event < 0 || n_event > n_events || seqpos < 0 || seqpos > seqlen) ? (Abort("Index (%d,%d) out of bounds for dimensions (%d,%d)",seqpos,n_event,seqlen,n_events), 0) : Unsafe_seq_event_pair_index(seqpos,n_event))
#else /* SEQEVTPAIR_DEBUG */
#define Seq_event_pair_index(seqpos,n_event) Unsafe_seq_event_pair_index(seqpos,n_event)
#endif /* SEQEVTPAIR_DEBUG */
unsigned long seq_event_pair_index_wrapper (Seq_event_pair_data* data, int seqpos, int n_event);

/* Forward-backward matrix */

typedef struct Seq_event_pair_fb_matrix {
  /* data */
  Seq_event_pair_data* data;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *fwdStart, *fwdMatch, *fwdDelete, fwdTotal;
  long double *backStart, *backMatch, *backDelete, backTotal;
} Seq_event_pair_fb_matrix;

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

/* DP helper functions */

double log_gaussian_density (double x, double mean, double precision, double log_precision);
double log_event_density (Fast5_event* event, double mean, double precision, double log_precision);

void dump_seq_event_pair_matrix (FILE* file, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxDelete, long double mxTotal);
void dump_seq_event_pair_matrix_to_file (const char* filename, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxDelete, long double mxTotal);

/* Expected counts */

typedef struct Seq_event_pair_counts {
  long double nStartEmitYes, nStartEmitNo;
  long double nNullEmitYes, nNullEmitNo;
  long double *nMatchEmitYes, *nMatchEmitNo;
  long double nSkipYes, nSkipNo;
  long double nBeginDeleteYes, nBeginDeleteNo;
  long double nExtendDeleteYes, nExtendDeleteNo;
  int order, states;
  long double *matchMoment0, *matchMoment1, *matchMoment2;
  long double nullMoment0, nullMoment1, nullMoment2;
  long double loglike;
} Seq_event_pair_counts;

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model);
void delete_seq_event_pair_counts (Seq_event_pair_counts* counts);

xmlChar* convert_seq_event_pair_counts_to_xml_string (Seq_event_pair_counts* counts);

Seq_event_pair_counts* new_seq_event_pair_counts_minimal_prior (Seq_event_pair_model* model);  /* Laplace +1's */

void reset_seq_event_null_counts (Seq_event_pair_counts* counts);
void reset_seq_event_pair_counts (Seq_event_pair_counts* counts);

void add_weighted_seq_event_pair_counts (Seq_event_pair_counts* counts, Seq_event_pair_counts** inc, int n_inc);  /* weights each set of counts in inc[] proportionally to its log-likelihood */

void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);  /* NB does not touch delete counts or loglike */
void inc_seq_event_null_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events);

void optimize_seq_event_null_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);
void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

void optimize_seq_event_model_for_events (Seq_event_pair_model* model, Vector* event_arrays);
int init_seq_event_model_from_fast5 (Seq_event_pair_model* model, const char* filename);  /* returns nonzero for failure */

/* Single Baum-Welch iteration */

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
void inc_seq_event_pair_counts_via_fb (Seq_event_pair_model* model, Seq_event_pair_counts* counts, int seqlen, char *seq, Fast5_event_array* events);  /* wraps creation & destruction of FB matrix */
Seq_event_pair_counts* get_seq_event_pair_counts (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays, int both_strands);  /* wraps iteration over sequences and reads */

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays, int both_strands);
void fit_seq_event_null_model (Seq_event_pair_model* model, Vector* event_arrays);
void fit_seq_event_kmer_model (Seq_event_pair_model* model, Kseq_container* seqs);

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

typedef void (*WriteSeqEventPairAlignmentFunction) (Seq_event_pair_alignment*, Seq_event_pair_model*, int, FILE*);
void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, Seq_event_pair_model* model, int strand, FILE* out);
void write_seq_event_pair_alignment_as_stockholm (Seq_event_pair_alignment* align, Seq_event_pair_model* model, int strand, FILE* out);

/* Viterbi matrix */

typedef struct Seq_event_pair_viterbi_matrix {
  /* data */
  Seq_event_pair_data* data;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *vitStart, *vitMatch, *vitDelete, vitTotal;
} Seq_event_pair_viterbi_matrix;

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);

void fill_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);
Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix);

/* wrappers */

void print_seq_evt_pair_alignments_as_gff_cigar (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, int both_strands, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold);
void print_seq_evt_pair_alignments_as_stockholm (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, int both_strands, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold);
void print_seq_evt_pair_alignments_generic (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, int both_strands, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold, WriteSeqEventPairAlignmentFunction write_func);

#endif /* SEQEVTPAIR_INCLUDED */
