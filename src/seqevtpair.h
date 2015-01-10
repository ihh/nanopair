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

/* Forward-backward matrix */

typedef struct Seq_event_pair_fb_matrix {
  /* parameters */
  Seq_event_pair_model* model;
  /* data */
  int seqlen;
  char *seq;
  int *state;
  Fast5_event_array* events;
  /* dynamic programming matrices: entries are all in log-probability space */
  long double *fwdStart, *fwdMatch, *fwdDelete, fwdEnd;
  long double *backStart, *backMatch, *backDelete;
  /* precalculated emit & transition scores */
  long double *startEmitDensity, *matchEmitDensity, *matchEmitYes, *matchEmitNo;
} Seq_event_pair_fb_matrix;

/* Seq_event_pair_index macro assumes seqlen & order are the same as in Seq_event_pair_fb_matrix */
#define Seq_event_pair_index(seqpos,n_event) (n_event*(seqlen-order+1) + seqpos - order)

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

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

/* accum_count(back_src,fwd_src,trans,back_dest,matrix,count,event,moment0,moment1,moment2)
   increments *count1 and *count2 by weight = exp(fwd_src + trans + back_dest - matrix->fwdEnd)
   also increments (moment0,moment1,moment2) by weight-ed event->(ticks,sumticks_cur,sumticks_cur_sq)
   returns log_sum_exp(back_src,trans + back_dest)
 */
long double accum_count (long double back_src,
			 long double fwd_src,
			 long double trans,
			 long double back_dest,
			 Seq_event_pair_fb_matrix *matrix,
			 long double *count1,
			 long double *count2,
			 Fast5_event *event,
			 long double *moment0,
			 long double *moment1,
			 long double *moment2);

/* Single Baum-Welch iteration */

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* matrix, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

/* Full Baum-Welch wrapper */

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Vector* event_arrays, double minimum_fractional_log_likelihood_increase);

#endif /* SEQEVTPAIR_INCLUDED */
