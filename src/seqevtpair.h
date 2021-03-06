#ifndef SEQEVTPAIR_INCLUDED
#define SEQEVTPAIR_INCLUDED

#include "fast5events.h"
#include "kseqcontainer.h"
#include "xmlparser.h"
#include "logger.h"

/* state encoding functions */
#define AlphabetSize 4

int base2token (char base);
char token2base (int token);

void encode_state_identifier (int state, int order, char* state_id);
int decode_state_identifier (int order, char* state_id);

/* Algorithm configuration */
typedef struct Seq_event_pair_config {
  int seq_evt_pair_EM_max_iterations;
  double seq_evt_pair_EM_min_fractional_loglike_increment;
  const char *debug_matrix_filename;
  int both_strands, all_vs_all;
} Seq_event_pair_config;

void init_seq_event_pair_config (Seq_event_pair_config *config);
int parse_seq_event_pair_config_general (int* argcPtr, char*** argvPtr, Seq_event_pair_config *config);
int parse_seq_event_pair_config_training (int* argcPtr, char*** argvPtr, Seq_event_pair_config *config);

/* Parameters */

typedef struct Seq_event_pair_model {
  /* transducer */
  double *pMatchSkip, *pMatchEvent, *pMatchTick;
  double pBeginDelete, pExtendDelete;
  int order, states;
  double *matchMean, *matchPrecision;
  double pStartEvent;
  double pNullEvent, pNullTick, nullMean, nullPrecision;
  /* generator */
  double *kmerProb;
  double emitProb;
  /* prior */
  StringDoubleMap *prior;
  /* configuration */
  Logger *logger;  /* not owned */
  Seq_event_pair_config config;
} Seq_event_pair_model;

StringDoubleMap* new_seq_event_pair_model_default_prior();

Seq_event_pair_model* new_seq_event_pair_model (int order);
Seq_event_pair_model* new_seq_event_pair_model_from_prior (int order, StringDoubleMap *prior);
void delete_seq_event_pair_model (Seq_event_pair_model* model);

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml);
xmlChar* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model);

xmlChar* make_squiggle_svg (Fast5_event_array *events, Seq_event_pair_model* model);

double get_seq_event_pair_pseudocount (Seq_event_pair_model* model, const char* param);
double get_seq_event_prior_mode (Seq_event_pair_model* model, const char* param);
double get_seq_event_log_beta_prior (Seq_event_pair_model* model, const char* param, double x);
#define BetaMode(YesCount,NoCount) (1. / (1. + ((NoCount) / (YesCount))))

/* Representations of model states, transitions, paths */

typedef enum { PairMatchState = 0, PairSkipState = 1, PairDeleteState = 2, PairStartState = 3, PairEndState = 4 } Seq_event_pair_state;

typedef enum { UndefinedTransition, StartStartOutTransition, StartMatchOutTransition, MatchMatchOutTransition, MatchSkipInTransition, SkipSkipInTransition, MatchMatchInOutTransition, SkipMatchInOutTransition, MatchDeleteInTransition, DeleteDeleteInTransition, DeleteMatchInOutTransition, MatchEndTransition } Seq_event_pair_transition;

typedef struct Labeled_seq_event_pair_transition {
  Seq_event_pair_model* model;
  Seq_event_pair_transition trans;
  int src_kmer, dest_kmer;
  Fast5_event* emission;
} Labeled_seq_event_pair_transition;

Labeled_seq_event_pair_transition* new_labeled_seq_event_pair_transition (Seq_event_pair_model* model, Seq_event_pair_transition trans, int src_kmer, int dest_kmer, Fast5_event* emission);

void* copy_labeled_seq_event_pair_transition (void*);
void delete_labeled_seq_event_pair_transition (void*);
void print_labeled_seq_event_pair_transition (FILE*, void*);

typedef List Labeled_seq_event_pair_path;
Labeled_seq_event_pair_path* new_labeled_seq_event_pair_path();
void delete_labeled_seq_event_pair_path (Labeled_seq_event_pair_path*);
void print_labeled_seq_event_pair_path (FILE*, Labeled_seq_event_pair_path*);

int seq_event_pair_transition_absorb_count (Seq_event_pair_transition trans);
int seq_event_pair_transition_emit_count (Seq_event_pair_transition trans);

long double seq_event_pair_transition_trans_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int src_kmer_state, int dest_kmer_state);
long double seq_event_pair_transition_emit_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int kmer_state, Fast5_event* emission);
long double seq_event_pair_transition_generator_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int dest_kmer_state);

long double labeled_seq_event_pair_transition_conditional_loglike (Labeled_seq_event_pair_transition*);
long double labeled_seq_event_pair_path_conditional_loglike (Labeled_seq_event_pair_path*);

long double labeled_seq_event_pair_transition_joint_loglike (Labeled_seq_event_pair_transition*);
long double labeled_seq_event_pair_path_joint_loglike (Labeled_seq_event_pair_path*);

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
  long double *nullEventLogLike, *matchEventLogLike, *matchSkipYes, *matchSkipNo, *matchEventYes, *matchEventNo;
  long double nullEventYes, nullEventNo;
  long double startEventYes, startEventNo;
  long double beginDeleteYes, beginDeleteNo;
  long double extendDeleteYes, extendDeleteNo;
  /* null model log-likelihood */
  long double nullModel;
  /* prior log-likelihoods (currently probability parameters only; ignores prior on Gaussian params) */
  long double logPrior, logNullPrior;
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
  long double *fwdStart, *fwdMatch, *fwdSkip, *fwdDelete, fwdResult;
  long double *backStart, *backMatch, *backSkip, *backDelete, backResult;
} Seq_event_pair_fb_matrix;

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

/* DP helper functions */

double log_event_density (Fast5_event* event, double mean, double precision, double log_precision, double log_pTick, double log_pNoTick);

void dump_seq_event_pair_matrix (FILE* file, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxSkip, long double *mxDelete, long double mxResult);
void dump_seq_event_pair_matrix_to_file (const char* filename, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxSkip, long double *mxDelete, long double mxResult);

/* Expected counts */

typedef struct Seq_event_pair_counts {
  long double nStartEventYes, nStartEventNo;
  long double nNullEventYes, nNullEventNo;
  long double *nMatchEventYes, *nMatchEventNo;
  long double *nMatchSkipYes, *nMatchSkipNo;
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
Seq_event_pair_counts* get_seq_event_pair_counts (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays);  /* wraps iteration over sequences and reads */

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays);
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
  Labeled_seq_event_pair_path* basecall_path;  /* equivalent path through generator (basecalling) HMM */
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
  long double *vitStart, *vitMatch, *vitSkip, *vitDelete, vitResult;
} Seq_event_pair_viterbi_matrix;

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);

void fill_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix);
Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix);

/* wrappers */

void print_seq_evt_pair_alignments_as_gff_cigar (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold);
void print_seq_evt_pair_alignments_as_stockholm (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold);
void print_seq_evt_pair_alignments_generic (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold, WriteSeqEventPairAlignmentFunction write_func);

#endif /* SEQEVTPAIR_INCLUDED */
