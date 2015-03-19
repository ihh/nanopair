#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include "seqevtpair.h"
#include "xmlutil.h"
#include "xmlkeywords.h"
#include "kseqcontainer.h"
#include "logsumexp.h"

/* logging stuff */
#define SEQEVTPAIR_PRECALC_WARN_PERIOD 1000000
#define SEQEVTPAIR_DP_WARN_PERIOD 500

/* log(sqrt(2*pi)) */
static const double log_sqrt2pi = 1.83787706640935;

/* path we check for model data in HDF5 file */
const char* model_path = "/Analyses/Basecall_2D_000/BaseCalled_template/Model";

/* names of various member fields in HDF5 file */
#define FAST5_MODEL_KMER "kmer"
#define FAST5_MODEL_LEVEL_MEAN "level_mean"
#define FAST5_MODEL_LEVEL_STDV "level_stdv"

/* default names */
static const char* unknown_seqname = "UnknownSequence";
static const char* gff3_source = "nanopair";
static const char* gff3_feature = "Match";
static const char* gff3_gap_attribute = "Gap";

const double seq_evt_pair_EM_max_iterations = 0.0001;
const double seq_evt_pair_EM_min_fractional_loglike_increment = 0.0001;

int base2token (char base);
char token2base (int token);

void encode_state_identifier (int state, int order, char* state_id);
int decode_state_identifier (int order, char* state_id);

/* accum_count(back_src,fwd_src,trans,back_dest,matrix,count,event,moment0,moment1,moment2)
   increments *count1 and *count2 by (weight = exp(fwd_src + trans + back_dest - matrix->fwdEnd)) * event->ticks
   also increments (moment0,moment1,moment2) by weight * event->(ticks,sumticks_cur,sumticks_cur_sq)
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

/* helper for finding index of max item in a list */
void update_max (long double *current_max, int* current_max_idx, long double candidate_max, int candidate_max_idx);

/* Metrichor_state_iterator
   Used to populate the emission parameters of a Seq_event_pair_model from a Metrichor model description in a FAST5 file */
typedef struct Metrichor_state_iterator {
  Seq_event_pair_model *model;
  size_t kmer_offset, level_mean_offset, level_sd_offset;
} Metrichor_state_iterator;

herr_t populate_seq_event_model_emit_params (void *elem, hid_t type_id, unsigned ndim, 
					     const hsize_t *point, void *operator_data)
{
  Metrichor_state_iterator *iter = (Metrichor_state_iterator*) operator_data;
  int state = decode_state_identifier (iter->model->order, (char*) elem + iter->kmer_offset);
  iter->model->matchMean[state] = *((double*) (elem + iter->level_mean_offset));
  double sd = *((double*) (elem + iter->level_sd_offset));
  iter->model->matchPrecision[state] = 1 / (sd * sd);
  return 0;
}

/* main function bodies */
Seq_event_pair_model* new_seq_event_pair_model (int order) {
  Seq_event_pair_model *model;
  model = SafeMalloc (sizeof (Seq_event_pair_model));
  model->order = order;
  model->states = pow(4,order);
  model->pMatchEmit = SafeMalloc (model->states * sizeof(double));
  model->matchMean = SafeMalloc (model->states * sizeof(double));
  model->matchPrecision = SafeMalloc (model->states * sizeof(double));
  return model;
}

void delete_seq_event_pair_model (Seq_event_pair_model* model) {
  SafeFree (model->pMatchEmit);
  SafeFree (model->matchMean);
  SafeFree (model->matchPrecision);
  SafeFree (model);
}

int base2token (char base) {
  return tokenize (toupper(base), dna_alphabet);
}

char token2base (int tok) {
  return tok < 0 || tok >= 4 ? 'N' : dna_alphabet[tok];
}

void encode_state_identifier (int state, int order, char* state_id) {
  int k;
  for (k = 0; k < order; ++k, state = state < 0 ? state : state / 4)
    state_id[order - k - 1] = token2base (state < 0 ? -1 : (state % 4));
  state_id[order] = '\0';
}

int decode_state_identifier (int order, char* state_id) {
  int k, token, chartok, mul;
  for (token = 0, mul = 1, k = 0; k < order; ++k, mul *= 4) {
    chartok = base2token (state_id[order - k - 1]);
    if (chartok < 0)
      return -1;
    token += mul * chartok;
  }
  return token;
}

double emitProbToMeanLength (double p) { return p / (1. - p); }
double meanLengthToEmitProb (double l) { return l / (1. + l); }

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml) {
  xmlNode *modelNode, *statesNode, *stateNode, *deleteNode, *startNode, *nullNode;
  Seq_event_pair_model *model;
  int state;
  modelNode = xmlTreeFromString (xml);
  model = new_seq_event_pair_model (CHILDINT(modelNode,ORDER));

  deleteNode = CHILD(modelNode,DELETE);
  model->pBeginDelete = CHILDFLOAT(deleteNode,BEGIN);
  model->pExtendDelete = CHILDFLOAT(deleteNode,EXTEND);

  statesNode = CHILD(modelNode,STATES);
  for (stateNode = statesNode->children; stateNode; stateNode = stateNode->next)
    if (MATCHES(stateNode,STATE)) {
      state = decode_state_identifier (model->order, (char*) CHILDSTRING(stateNode,KMER));
      model->pMatchEmit[state] = meanLengthToEmitProb (CHILDFLOAT(stateNode,WAIT));
      model->matchMean[state] = CHILDFLOAT(stateNode,MEAN);
      model->matchPrecision[state] = 1 / pow (CHILDFLOAT(stateNode,STDV), 2);
    }

  startNode = CHILD(modelNode,START);
  model->pStartEmit = meanLengthToEmitProb (CHILDFLOAT(startNode,WAIT));

  nullNode = CHILD(modelNode,NULLMODEL);
  model->pNullEmit = meanLengthToEmitProb (CHILDFLOAT(nullNode,WAIT));
  model->nullMean = CHILDFLOAT(nullNode,MEAN);
  model->nullPrecision = 1 / pow (CHILDFLOAT(nullNode,STDV), 2);

  deleteXmlTree (modelNode);
  return model;
}

xmlChar* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model) {
  xmlTextWriterPtr writer;
  char* id;
  int state;

  id = SafeMalloc ((model->order + 1) * sizeof(char));

  writer = newXmlTextWriter();
  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(MODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(ORDER), "%d", model->order);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(DELETE));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(BEGIN), "%g", model->pBeginDelete);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EXTEND), "%g", model->pExtendDelete);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATES));
  for (state = 0; state < model->states; ++state) {
    xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATE));
    encode_state_identifier (state, model->order, id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(KMER), "%s", id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(WAIT), "%g", emitProbToMeanLength (model->pMatchEmit[state]));
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->matchMean[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(STDV), "%g", 1 / sqrt(model->matchPrecision[state]));
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(START));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(WAIT), "%g", emitProbToMeanLength (model->pStartEmit));
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(NULLMODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(WAIT), "%g", emitProbToMeanLength (model->pNullEmit));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->nullMean);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(STDV), "%g", 1 / sqrt(model->nullPrecision));
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);

  SafeFree (id);

  return deleteXmlTextWriterLeavingText (writer);
}

/* Pairwise alignment data */

Seq_event_pair_data* new_seq_event_pair_data (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_data* data;
  int n_events, seqpos;
  unsigned long matrix_cells;

  Assert (seqlen >= model->order, "Sequence length is %d, which is less than model order (%d)", seqlen, model->order);

  data = SafeMalloc (sizeof (Seq_event_pair_data));
  data->model = model;
  data->seqlen = seqlen;
  data->seq = seq;
  data->events = events;

  n_events = events->n_events;

  matrix_cells = ((unsigned long) (n_events + 1)) * (unsigned long) (seqlen - model->order + 1);
  Warn ("Allocating DP scratch space with %d*%d = %Lu cells", n_events + 1, seqlen - model->order + 1, matrix_cells);

  data->matrix_cells = matrix_cells;

  data->nullEmitDensity = SafeMalloc ((n_events + 1) * sizeof(long double));

  data->matchEmitDensity = SafeMalloc (matrix_cells * sizeof(long double));
  data->matchEmitYes = SafeMalloc ((seqlen + 1) * sizeof(long double));
  data->matchEmitNo = SafeMalloc ((seqlen + 1) * sizeof(long double));

  data->state = SafeMalloc ((seqlen + 1) * sizeof(int));

  /* calculate states */
  for (seqpos = 0; seqpos < model->order; ++seqpos)
    data->state[seqpos] = -1;
  for (seqpos = model->order; seqpos <= seqlen; ++seqpos)
    data->state[seqpos] = decode_state_identifier (model->order, data->seq + seqpos - model->order);

  return data;
}

void delete_seq_event_pair_data (Seq_event_pair_data* data) {
  SafeFree (data->state);
  SafeFree (data->nullEmitDensity);
  SafeFree (data->matchEmitDensity);
  SafeFree (data->matchEmitYes);
  SafeFree (data->matchEmitNo);
  SafeFree (data);
}

void precalc_seq_event_pair_data (Seq_event_pair_data* data) {
  Seq_event_pair_model* model;
  int seqlen, n_events, order, seqpos, n_event, state;
  long double logNullPrecision, loglike;
  double mean, precision, logPrecision;
  Fast5_event* event;

  model = data->model;
  order = model->order;
  seqlen = data->seqlen;
  n_events = data->events->n_events;

  /* calculate logs of transition probabilities & emit precisions */
  data->nullEmitYes = log (model->pNullEmit);
  data->nullEmitNo = log (1. - model->pNullEmit);

  data->startEmitYes = log (model->pStartEmit);
  data->startEmitNo = log (1. - model->pStartEmit);

  data->beginDeleteYes = log (model->pBeginDelete);
  data->beginDeleteNo = log (1. - model->pBeginDelete);

  data->extendDeleteYes = log (model->pExtendDelete);
  data->extendDeleteNo = log (1. - model->pExtendDelete);

  logNullPrecision = log (model->nullPrecision);

  /* calculate emit densities & state-dependent transition probabilities */
  data->nullEmitDensity[0] = -INFINITY;
  data->nullModel = data->nullEmitNo;
  for (n_event = 0; n_event < n_events; ++n_event) {
    event = &data->events->event[n_event];
    loglike = log_event_density (event,
				 model->nullMean,
				 model->nullPrecision,
				 logNullPrecision);
    data->nullEmitDensity[n_event + 1] = loglike;
    data->nullModel += loglike + data->nullEmitYes * event->ticks;
  }

  for (seqpos = 0; seqpos < order; ++seqpos) {
    data->matchEmitYes[seqpos] = -INFINITY;
    data->matchEmitNo[seqpos] = -INFINITY;
  }

  for (seqpos = order; seqpos <= seqlen; ++seqpos) {

#ifdef SEQEVTPAIR_DEBUG
    if (seqpos % SEQEVTPAIR_PRECALC_WARN_PERIOD == 0)
      Warn ("Precalculating likelihoods at sequence position %d", seqpos);
#endif /* SEQEVTPAIR_DEBUG */

    state = data->state[seqpos];

    data->matchEmitYes[seqpos] = state < 0 ? data->nullEmitYes : log (model->pMatchEmit[state]);
    data->matchEmitNo[seqpos] = state < 0 ? data->nullEmitNo : log (1. - model->pMatchEmit[state]);

    mean = state < 0 ? model->nullMean : model->matchMean[state];
    precision = state < 0 ? model->nullPrecision : model->matchPrecision[state];
    logPrecision = log (precision);

    data->matchEmitDensity[Seq_event_pair_index(seqpos,0)] = -INFINITY;
    for (n_event = 1; n_event <= n_events; ++n_event) {
      event = &data->events->event[n_event - 1];
      data->matchEmitDensity[Seq_event_pair_index(seqpos,n_event)]
	= log_event_density (event,
			     mean,
			     precision,
			     logPrecision);
    }
  }
}

unsigned long seq_event_pair_index_wrapper (Seq_event_pair_data* data, int seqpos, int n_event) {
  int seqlen, n_events, order;
  seqlen = data->seqlen;
  order = data->model->order;
  n_events = data->events->n_events;
  return Seq_event_pair_index(seqpos,n_event);
}

/* Forward-backward matrix */

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_fb_matrix* mx;
  int n_events;
  unsigned long matrix_cells;

  mx = SafeMalloc (sizeof (Seq_event_pair_fb_matrix));
  mx->data = new_seq_event_pair_data (model, seqlen, seq, events);

  n_events = mx->data->events->n_events;
  matrix_cells = mx->data->matrix_cells;

#ifdef SEQEVTPAIR_DEBUG
  Warn ("Allocating Forward-Backward matrix of size %d*%d (approx.)", n_events, seqlen);
#endif /* SEQEVTPAIR_DEBUG */

  mx->fwdStart = SafeMalloc ((n_events + 1) * sizeof(long double));
  mx->backStart = SafeMalloc ((n_events + 1) * sizeof(long double));

  mx->fwdMatch = SafeMalloc (matrix_cells * sizeof(long double));
  mx->fwdDelete = SafeMalloc (matrix_cells * sizeof(long double));
  mx->backMatch = SafeMalloc (matrix_cells * sizeof(long double));
  mx->backDelete = SafeMalloc (matrix_cells * sizeof(long double));

  return mx;
}

void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* mx) {
  SafeFree (mx->fwdStart);
  SafeFree (mx->backStart);
  SafeFree (mx->fwdMatch);
  SafeFree (mx->fwdDelete);
  SafeFree (mx->backMatch);
  SafeFree (mx->backDelete);
  SafeFree (mx->data);
  SafeFree (mx);
}

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model) {
  Seq_event_pair_counts* counts;
  counts = SafeMalloc (sizeof (Seq_event_pair_counts));
  counts->order = model->order;
  counts->states = model->states;
  counts->nMatchEmitYes = SafeMalloc (model->states * sizeof(long double));
  counts->nMatchEmitNo = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment0 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment1 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment2 = SafeMalloc (model->states * sizeof(long double));
  return counts;
}

void delete_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  SafeFree (counts->nMatchEmitYes);
  SafeFree (counts->nMatchEmitNo);
  SafeFree (counts->matchMoment0);
  SafeFree (counts->matchMoment1);
  SafeFree (counts->matchMoment2);
  SafeFree (counts);
}

void reset_seq_event_null_counts (Seq_event_pair_counts* counts) {
  counts->nNullEmitYes = 0;
  counts->nNullEmitNo = 0;
  counts->nullMoment0 = 0;
  counts->nullMoment1 = 0;
  counts->nullMoment2 = 0;
}

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  int state;
  for (state = 0; state < counts->states; ++state) {
    counts->nMatchEmitYes[state] = 0;
    counts->nMatchEmitNo[state] = 0;
    counts->matchMoment0[state] = 0;
    counts->matchMoment1[state] = 0;
    counts->matchMoment2[state] = 0;
  }
  counts->nStartEmitYes = 0;
  counts->nStartEmitNo = 0;
  counts->nBeginDeleteYes = 0;
  counts->nBeginDeleteNo = 0;
  counts->nExtendDeleteYes = 0;
  counts->nExtendDeleteNo = 0;
  counts->loglike = 0;
}

Seq_event_pair_counts* new_seq_event_pair_counts_minimal_prior (Seq_event_pair_model* model) {
  Seq_event_pair_counts* counts;
  int state;

  counts = new_seq_event_pair_counts (model);
  reset_seq_event_null_counts (counts);
  reset_seq_event_pair_counts (counts);

  for (state = 0; state < counts->states; ++state) {
    counts->nMatchEmitYes[state] += 1.;
    counts->nMatchEmitNo[state] += 1.;
    counts->matchMoment0[state] += 1.;
  }
  counts->nStartEmitYes += 1.;
  counts->nStartEmitNo += 1.;
  counts->nBeginDeleteYes += 1.;
  counts->nBeginDeleteNo += 1.;
  counts->nExtendDeleteYes += 1.;
  counts->nExtendDeleteNo += 1.;
  counts->nNullEmitYes += 1.;
  counts->nNullEmitNo += 1.;
  /* counts->nullMoment0 is untouched */

  return counts;
}

void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events) {
  int n_event, state;
  long moves, new_moves;
  Fast5_event *event;
  moves = 0;
  for (n_event = 0; n_event < events->n_events; ++n_event) {
    event = &events->event[n_event];
    new_moves = moves + event->move;
    if (new_moves < counts->order) {
      counts->nStartEmitYes += 1.;
      counts->nullMoment0 += event->ticks;
      counts->nullMoment1 += event->sumticks_cur;
      counts->nullMoment2 += event->sumticks_cur_sq;
    } else {
      state = decode_state_identifier (counts->order, event->model_state);
      if (moves < counts->order)
	counts->nStartEmitNo += 1.;
      counts->nMatchEmitNo[state] += event->move;
      counts->nMatchEmitYes[state] += event->ticks;
      counts->matchMoment0[state] += event->ticks;
      counts->matchMoment1[state] += event->sumticks_cur;
      counts->matchMoment2[state] += event->sumticks_cur_sq;
    }
    counts->nNullEmitYes += 1.;
    moves = new_moves;
  }
  counts->nNullEmitNo += 1.;
}

void add_weighted_seq_event_pair_counts (Seq_event_pair_counts* counts, Seq_event_pair_counts** inc, int n_inc) {
  int state, n;
  long double loglike_all = -INFINITY;
  double weight;
  for (n = 0; n < n_inc; ++n)
    loglike_all = log_sum_exp (loglike_all, inc[n]->loglike);
  for (n = 0; n < n_inc; ++n) {
    weight = exp (inc[n]->loglike - loglike_all);
    for (state = 0; state < counts->states; ++state) {
      counts->nMatchEmitYes[state] += weight * inc[n]->nMatchEmitYes[state];
      counts->nMatchEmitNo[state] += weight * inc[n]->nMatchEmitNo[state];
      counts->matchMoment0[state] += weight * inc[n]->matchMoment0[state];
      counts->matchMoment1[state] += weight * inc[n]->matchMoment1[state];
      counts->matchMoment2[state] += weight * inc[n]->matchMoment2[state];
    }
    counts->nStartEmitYes += weight * inc[n]->nStartEmitYes;
    counts->nStartEmitNo += weight * inc[n]->nStartEmitNo;
    counts->nBeginDeleteYes += weight * inc[n]->nBeginDeleteYes;
    counts->nBeginDeleteNo += weight * inc[n]->nBeginDeleteNo;
    counts->nExtendDeleteYes += weight * inc[n]->nExtendDeleteYes;
    counts->nExtendDeleteNo += weight * inc[n]->nExtendDeleteNo;
  }
  counts->loglike += loglike_all;
}

double log_gaussian_density (double x, double mean, double precision, double log_precision) {
  double xz;
  xz = x - mean;
  return log_precision/2. - log_sqrt2pi - precision*xz*xz/2.;
}

double log_event_density (Fast5_event* event, double mean, double precision, double log_precision) {
  return event->ticks * (log_precision/2. - log_sqrt2pi - precision*mean*mean/2)
    - precision*(event->sumticks_cur_sq/2. - event->sumticks_cur*mean);
}

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event, state;
  long double mat, del;
  unsigned long idx, inputIdx, outputIdx;
  Fast5_event* event;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;

  /* update data */
  precalc_seq_event_pair_data (data);

  /* fill forward */
  matrix->fwdStart[0] = 0.;
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    matrix->fwdStart[n_event]
      = matrix->fwdStart[n_event - 1]
      + data->startEmitYes * event->ticks
      + data->nullEmitDensity[n_event];  /* Start -> Start (output) */
  }

  matrix->fwdEnd = -INFINITY;

  for (n_event = 0; n_event <= n_events; ++n_event) {
#ifdef SEQEVTPAIR_DEBUG
    if (n_event % SEQEVTPAIR_DP_WARN_PERIOD == 0)
      Warn ("Filling forward matrix, event %d", n_event + 1);
#endif /* SEQEVTPAIR_DEBUG */

    for (seqpos = order; seqpos <= seqlen; ++seqpos) {
      idx = Seq_event_pair_index(seqpos,n_event);
      inputIdx = Seq_event_pair_index(seqpos-1,n_event);

      mat = matrix->fwdStart[n_event] + data->startEmitNo;   /* Start -> Match (input) */

      if (n_event > 0) {
	outputIdx = Seq_event_pair_index(seqpos,n_event-1);
	event = &data->events->event[n_event - 1];
	mat = log_sum_exp (mat,
			   matrix->fwdMatch[outputIdx]
			   + data->matchEmitYes[seqpos] * event->ticks
			   + data->matchEmitDensity[idx]);     /* Match -> Match (output) */
      }

      if (seqpos == order) {
	del = -INFINITY;
      } else {  /* seqpos > order */
	del = log_sum_exp
	  (matrix->fwdMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteYes,  /* Match -> Delete (input) */
	   matrix->fwdDelete[inputIdx] + data->extendDeleteYes);   /* Delete -> Delete (input) */

	mat = log_sum_exp
	  (mat,
	   matrix->fwdMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteNo);  /* Match -> Match (input) */

	mat = log_sum_exp (mat, del + data->extendDeleteNo);  /* Delete -> Match */
      }

      matrix->fwdMatch[idx] = mat;
      matrix->fwdDelete[idx] = del;

      matrix->fwdEnd = log_sum_exp
	(matrix->fwdEnd,
	 matrix->fwdMatch[Seq_event_pair_index(seqpos,n_events)] + data->matchEmitNo[seqpos]);  /* Match -> End (input) */
    }
  }

#if defined(SEQEVTPAIR_DEBUG) && SEQEVTPAIR_DEBUG >= 10
  dump_seq_event_pair_matrix_to_file (SEQEVTMATRIX_FILENAME, "Forward", data, matrix->fwdStart, matrix->fwdMatch, matrix->fwdDelete);
#endif /* SEQEVTPAIR_DEBUG >= 10 */

  /* fill backward & accumulate counts */
  for (n_event = n_events; n_event >= 0; --n_event) {
#ifdef SEQEVTPAIR_DEBUG
    if (n_event % SEQEVTPAIR_DP_WARN_PERIOD == 0)
      Warn ("Filling backward matrix, event %d", n_event + 1);
#endif /* SEQEVTPAIR_DEBUG */

    event = n_event < n_events ? &data->events->event[n_event] : NULL;

    matrix->backStart[n_event] = -INFINITY;

    for (seqpos = seqlen; seqpos >= order; --seqpos) {
      state = data->state[seqpos];
      idx = Seq_event_pair_index(seqpos,n_event);

      if (n_event == n_events) {
	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEmitNo[seqpos],
			   0.,
			   matrix,
			   state < 0 ? NULL : &counts->nMatchEmitNo[state],
			   NULL,
			   NULL, NULL, NULL, NULL);   /* Match -> End (input) */
      } else {  /* n_event < n_events */
	outputIdx = Seq_event_pair_index(seqpos,n_event+1);

	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEmitYes[seqpos] * event->ticks
			   + data->matchEmitDensity[outputIdx],
			   matrix->backMatch[outputIdx],
			   matrix,
			   state < 0 ? NULL : &counts->nMatchEmitYes[state],
			   NULL,
			   state < 0 ? NULL : event,
			   state < 0 ? NULL : &counts->matchMoment0[state],
			   state < 0 ? NULL : &counts->matchMoment1[state],
			   state < 0 ? NULL : &counts->matchMoment2[state]);  /* Match -> Match (output) */
      }

      if (seqpos == seqlen) {
	del = -INFINITY;
      } else {
	inputIdx = Seq_event_pair_index(seqpos+1,n_event);

	del = accum_count (-INFINITY,
			   matrix->fwdDelete[idx],
			   data->extendDeleteYes,
			   matrix->backDelete[inputIdx],
			   matrix,
			   &counts->nExtendDeleteYes, NULL,
			   NULL, NULL, NULL, NULL);   /* Delete -> Delete (input) */

	mat = accum_count (mat,
			   matrix->fwdMatch[idx],
			   data->matchEmitNo[seqpos] + data->beginDeleteYes,
			   matrix->backDelete[inputIdx],
			   matrix,
			   state < 0 ? NULL : &counts->nMatchEmitNo[state],
			   &counts->nBeginDeleteYes,
			   NULL, NULL, NULL, NULL);  /* Match -> Delete (input) */
	mat = accum_count (mat,
			   matrix->fwdMatch[idx],
			   data->matchEmitNo[seqpos] + data->beginDeleteNo,
			   matrix->backMatch[inputIdx],
			   matrix,
			   state < 0 ? NULL : &counts->nMatchEmitNo[state],
			   &counts->nBeginDeleteNo,
			   NULL, NULL, NULL, NULL);  /* Match -> Match (input) */
      }

      del = accum_count (del,
			 matrix->fwdDelete[idx],
			 data->extendDeleteNo,
			 mat,
			 matrix,
			 &counts->nExtendDeleteNo, NULL,
			 NULL, NULL, NULL, NULL);  /* Delete -> Match */

      matrix->backMatch[idx] = mat;
      matrix->backDelete[idx] = del;

      matrix->backStart[n_event] = accum_count (matrix->backStart[n_event],
						matrix->fwdStart[n_event],
						data->startEmitNo,
						mat,
						matrix,
						&counts->nStartEmitNo, NULL,
						NULL, NULL, NULL, NULL);   /* Start -> Match (input) */
    }

    if (n_event < n_events)
      matrix->backStart[n_event] = accum_count (matrix->backStart[n_event],
						matrix->fwdStart[n_event],
						data->startEmitYes * event->ticks
						+ data->nullEmitDensity[n_event + 1],
						matrix->backStart[n_event + 1],
						matrix,
						&counts->nStartEmitYes, NULL,
						NULL, NULL, NULL, NULL);  /* Start -> Start (output) */
  }

#if defined(SEQEVTPAIR_DEBUG) && SEQEVTPAIR_DEBUG >= 10
  dump_seq_event_pair_matrix_to_file (SEQEVTMATRIX_FILENAME, "Backward", data, matrix->backStart, matrix->backMatch, matrix->backDelete);
#endif /* SEQEVTPAIR_DEBUG >= 10 */
}

void inc_seq_event_null_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events) {
  Fast5_event* event;
  int n_event;
  for (n_event = 1; n_event <= events->n_events; ++n_event) {
    event = &events->event[n_event - 1];
    counts->nullMoment0 += event->ticks;
    counts->nullMoment1 += event->sumticks_cur;
    counts->nullMoment2 += event->sumticks_cur_sq;
  }
  counts->nNullEmitYes += events->n_events;
  counts->nNullEmitNo += 1.;
}

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
			 long double *moment2) {
  long double weight;
  weight = exp (fwd_src + trans + back_dest - matrix->fwdEnd);
  if (event) {
    *moment0 += weight * event->ticks;
    *moment1 += weight * event->sumticks_cur;
    *moment2 += weight * event->sumticks_cur_sq;
    weight *= event->ticks;
  }
  if (count1)
    *count1 += weight;
  if (count2)
    *count2 += weight;
  return log_sum_exp (back_src, trans + back_dest);
}

void optimize_seq_event_null_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior) {
  Seq_event_pair_counts* dummy_prior;

  if (prior == NULL) {
    dummy_prior = new_seq_event_pair_counts (model);
    reset_seq_event_pair_counts (dummy_prior);
    prior = dummy_prior;
  } else
    dummy_prior = NULL;

  model->pNullEmit = (counts->nNullEmitYes + prior->nNullEmitYes) / (counts->nNullEmitYes + prior->nNullEmitYes + counts->nNullEmitNo + prior->nNullEmitNo);
  model->nullMean = (counts->nullMoment1 + prior->nullMoment1) / (counts->nullMoment0 + prior->nullMoment0);
  model->nullPrecision = 1. / ((counts->nullMoment2 + prior->nullMoment2) / (counts->nullMoment0 + prior->nullMoment0) - model->nullMean * model->nullMean);

  if (dummy_prior != NULL)
    delete_seq_event_pair_counts (dummy_prior);
}

void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior) {
  int state;
  Seq_event_pair_counts* dummy_prior;

  if (prior == NULL) {
    dummy_prior = new_seq_event_pair_counts (model);
    reset_seq_event_pair_counts (dummy_prior);
    prior = dummy_prior;
  } else
    dummy_prior = NULL;

  for (state = 0; state < model->states; ++state) {
    model->pMatchEmit[state] = (counts->nMatchEmitYes[state] + prior->nMatchEmitYes[state]) / (counts->nMatchEmitYes[state] + prior->nMatchEmitYes[state] + counts->nMatchEmitNo[state] + prior->nMatchEmitNo[state]);
    model->matchMean[state] = (counts->matchMoment1[state] + prior->matchMoment1[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]);
    model->matchPrecision[state] = 1. / ((counts->matchMoment2[state] + prior->matchMoment2[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]) - model->matchMean[state] * model->matchMean[state]);
  }
  model->pBeginDelete = (counts->nBeginDeleteYes + prior->nBeginDeleteYes) / (counts->nBeginDeleteYes + prior->nBeginDeleteYes + counts->nBeginDeleteNo + prior->nBeginDeleteNo);
  model->pExtendDelete = (counts->nExtendDeleteYes + prior->nExtendDeleteYes) / (counts->nExtendDeleteYes + prior->nExtendDeleteYes + counts->nExtendDeleteNo + prior->nExtendDeleteNo);
  model->pStartEmit = (counts->nStartEmitYes + prior->nStartEmitYes) / (counts->nStartEmitYes + prior->nStartEmitYes + counts->nStartEmitNo + prior->nStartEmitNo);

  if (dummy_prior != NULL)
    delete_seq_event_pair_counts (dummy_prior);
}

void inc_seq_event_pair_counts_via_fb (Seq_event_pair_model* model, Seq_event_pair_counts* counts, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_fb_matrix* matrix;
  long double fwdEnd, nullModel;
  matrix = new_seq_event_pair_fb_matrix (model, seqlen, seq, events);
  fill_seq_event_pair_fb_matrix_and_inc_counts (matrix, counts);
  fwdEnd = matrix->fwdEnd;
  nullModel = matrix->data->nullModel;
  delete_seq_event_pair_fb_matrix (matrix);
  counts->loglike += fwdEnd - nullModel;
}

void optimize_seq_event_model_for_events (Seq_event_pair_model* model, Vector* event_arrays) {
  Seq_event_pair_counts *counts, *prior;
  void **events_iter;
  Fast5_event_array *events;

  prior = new_seq_event_pair_counts_minimal_prior (model);
  counts = new_seq_event_pair_counts (model);

  reset_seq_event_null_counts (counts);
  reset_seq_event_pair_counts (counts);
  for (events_iter = event_arrays->begin; events_iter != event_arrays->end; ++events_iter) {
    events = (Fast5_event_array*) *events_iter;
    inc_seq_event_null_counts_from_fast5 (counts, events);
    inc_seq_event_pair_counts_from_fast5 (counts, events);
  }
  optimize_seq_event_null_model_for_counts (model, counts, prior);
  optimize_seq_event_pair_model_for_counts (model, counts, prior);

  delete_seq_event_pair_counts (counts);
  delete_seq_event_pair_counts (prior);
}

int init_seq_event_model_from_fast5 (Seq_event_pair_model* model, const char* filename) {
  hid_t file_id, strtype_id;
  int ret = 0;

  /* open file with default properties */
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* check if opening file was succeful */
  if ( file_id < 0 )
    {
      Warn("Failed to open/init input file %s", filename);
      return -1;
    }

  /* get root group */
  hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if (root_id < 0)
    {
      Warn("Failed to open root group in file %s",filename);
      return -1;
    }

  /* see if path exists */
  if (H5LTpath_valid ( file_id, model_path, 1))
    {
      /* get dataset */
      hid_t model_id = H5Oopen(file_id, model_path, H5P_DEFAULT);
      if ( model_id < 0 )
	{
	  Warn("Failed to open dataset %s in file %s",model_path,filename);
	  ret = -1;
	}
      else
	{
	  /* get information about fields in model */
	  hid_t model_type_id = H5Dget_type( model_id );

	  Metrichor_state_iterator iter;

	  int kmer_idx = H5Tget_member_index( model_type_id, FAST5_MODEL_KMER );
	  int level_mean_idx = H5Tget_member_index( model_type_id, FAST5_MODEL_LEVEL_MEAN );
	  int level_sd_idx = H5Tget_member_index( model_type_id, FAST5_MODEL_LEVEL_STDV );

	  iter.kmer_offset = H5Tget_member_offset( model_type_id, kmer_idx );
	  iter.level_mean_offset = H5Tget_member_offset( model_type_id, level_mean_idx );
	  iter.level_sd_offset = H5Tget_member_offset( model_type_id, level_sd_idx );

	  /* check that the kmer strings are the same length as our model order */
	  strtype_id = H5Tget_member_type( model_type_id, kmer_idx );
	  int kmer_len = H5Tget_size (strtype_id);
	  Assert (model->order == kmer_len, "Length of kmers in file '%s' (%d) does not match model order supplied to %s (%d)", filename, kmer_len, __FUNCTION__, model->order);

	  /* get dimensions */
	  hid_t model_space_id = H5Dget_space( model_id );
	  hssize_t model_npoints = H5Sget_simple_extent_npoints( model_space_id );
	  size_t state_size = H5Tget_size( model_type_id );

	  /* read into memory buffer */
	  void *buf = SafeMalloc (model_npoints * state_size);
	  H5Dread( model_id, model_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );

	  /* convert */
	  iter.model = model;
	  H5Diterate( buf, model_type_id, model_space_id, populate_seq_event_model_emit_params, &iter );

	  /* free buffer */
	  SafeFree (buf);

	  /* close objects */
	  H5Tclose(strtype_id);
	  H5Sclose(model_space_id);
	  H5Tclose(model_type_id);
	  H5Dclose(model_id);
	}
    }
  else {
    Warn("Path %s not valid in file %s",model_path,filename);
    ret = -1;
  }

  /* close HDF5 resources */
  H5Gclose(root_id);
  H5Fclose(file_id);

  /* if we haven't screwed up yet, do some last-minute parameter fudging */
  if (ret == 0) {
    /* parameterize the null model by averaging moments of all the states */
    double m1 = 0, m2 = 0;
    int state;
    for (state = 0; state < model->states; ++state) {
      m1 += model->matchMean[state];
      m2 += model->matchMean[state] * model->matchMean[state] + 1 / model->matchPrecision[state];
    }
    m1 /= model->states;
    m2 /= model->states;
    model->nullMean = m1;
    model->nullPrecision = 1 / (m2 - m1 * m1);

    /* set all boolean probabilities to 0.5 */
    model->pBeginDelete = model->pExtendDelete = model->pStartEmit = model->pNullEmit = 0.5;
    for (state = 0; state < model->states; ++state)
      model->pMatchEmit[state] = 0.5;
  }

  /* return */
  return ret;
}

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays) {
  int iter;
  long double loglike, prev_loglike;
  Seq_event_pair_counts *counts, *prior;

  prior = new_seq_event_pair_counts_minimal_prior (model);

  prev_loglike = 0.;
  for (iter = 0; iter < seq_evt_pair_EM_max_iterations; ++iter) {
    counts = get_seq_event_pair_counts (model, seqs, event_arrays);
    loglike = counts->loglike;

    optimize_seq_event_pair_model_for_counts (model, counts, prior);
    delete_seq_event_pair_counts (counts);

#ifdef SEQEVTPAIR_DEBUG
    Warn ("Baum-Welch iteration %d: log-likelihood %Lg", iter + 1, loglike);
#endif /* SEQEVTPAIR_DEBUG */

    if (iter > 0 && prev_loglike != 0. && abs(prev_loglike) != INFINITY
	&& abs((loglike-prev_loglike)/prev_loglike) < seq_evt_pair_EM_min_fractional_loglike_increment)
      break;
    prev_loglike = loglike;
  }

  delete_seq_event_pair_counts (prior);
}

Seq_event_pair_counts* get_seq_event_pair_counts (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays) {
  int n_seq, *seqrev_len;
  char **rev, **seqrev;
  void **events_iter;
  Fast5_event_array *events;
  Seq_event_pair_counts *counts, **seq_counts;

  rev = SafeMalloc (seqs->n * sizeof(char*));
  for (n_seq = 0; n_seq < seqs->n; ++n_seq)
    rev[n_seq] = new_revcomp_seq (seqs->seq[n_seq], seqs->len[n_seq]);

  seqrev = SafeMalloc (2 * seqs->n * sizeof(char*));
  seqrev_len = SafeMalloc (2 * seqs->n * sizeof(int));
  for (n_seq = 0; n_seq < seqs->n; ++n_seq) {
    seqrev[2*n_seq] = seqs->seq[n_seq];
    seqrev[2*n_seq+1] = rev[n_seq];
    seqrev_len[2*n_seq] = seqrev_len[2*n_seq+1] = seqs->len[n_seq];
  }

  counts = new_seq_event_pair_counts (model);
  reset_seq_event_pair_counts (counts);

  seq_counts = SafeMalloc (2 * seqs->n * sizeof(Seq_event_pair_counts*));
  for (n_seq = 0; n_seq < 2 * seqs->n; ++n_seq)
    seq_counts[n_seq] = new_seq_event_pair_counts (model);

  for (events_iter = event_arrays->begin; events_iter != event_arrays->end; ++events_iter) {
    events = (Fast5_event_array*) *events_iter;
    for (n_seq = 0; n_seq < 2 * seqs->n; ++n_seq) {
      reset_seq_event_pair_counts (seq_counts[n_seq]);
      inc_seq_event_pair_counts_via_fb (model, seq_counts[n_seq], seqrev_len[n_seq], seqrev[n_seq], events);
#ifdef SEQEVTPAIR_DEBUG
      Warn ("FAST5 file \"%s\", sequence \"%s\" (%s strand): log-likelihood = %Lg", events->name == NULL ? "<none>" : events->name, seqs->name[n_seq/2], (n_seq % 2) ? "reverse" : "forward", seq_counts[n_seq]->loglike);
#endif /* SEQEVTPAIR_DEBUG */
    }
    add_weighted_seq_event_pair_counts (counts, seq_counts, 2 * seqs->n);
  }

  for (n_seq = 0; n_seq < seqs->n; ++n_seq)
    SafeFree (rev[n_seq]);
  SafeFree (rev);

  for (n_seq = 0; n_seq < 2 * seqs->n; ++n_seq)
    delete_seq_event_pair_counts (seq_counts[n_seq]);
  SafeFree (seq_counts);
  SafeFree (seqrev);
  SafeFree (seqrev_len);

  return counts;
}

Seq_event_pair_alignment* new_seq_event_pair_alignment (Fast5_event_array *events, char *seq, int seqlen) {
  Seq_event_pair_alignment* align;
  align = SafeMalloc (sizeof (Seq_event_pair_alignment));
  align->events = events;
  align->seqname = unknown_seqname;
  align->seq = seq;
  align->seqlen = seqlen;
  align->events_at_pos = NULL;
  align->start_seqpos = 0;
  align->end_seqpos = 0;
  align->start_n_event = events->n_events;
  return align;
}

void delete_seq_event_pair_alignment (Seq_event_pair_alignment* align) {
  SafeFreeOrNull (align->events_at_pos);
  SafeFree (align);
}

void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, int strand, FILE* out) {
  int n, deleted;
  StringVector *cigar;
  char buf[30];

  fprintf (out, "%s\t%s\t%s\t%d\t%d\t%Lg\t%s\t.\t%s=",
	   align->seqname, gff3_source, gff3_feature,
	   strand > 0 ? (align->start_seqpos + 1) : (align->seqlen - align->end_seqpos),
	   strand > 0 ? (align->end_seqpos + 1) : (align->seqlen - align->start_seqpos),
	   align->log_likelihood_ratio,
	   strand > 0 ? "+" : "-",
	   gff3_gap_attribute);

  cigar = newStringVector();
  if (align->start_n_event > 0) {
    sprintf (buf, "I%d ", align->start_n_event);
    StringVectorPushBack (cigar, buf);
  }

  deleted = 0;
  for (n = 0; n <= align->end_seqpos - align->start_seqpos; ++n) {
    ++deleted;
    if (align->events_at_pos[n] > 0) {
      sprintf (buf, "D%d ", deleted);
      StringVectorPushBack (cigar, buf);
      
      sprintf (buf, "I%d ", align->events_at_pos[n]);
      StringVectorPushBack (cigar, buf);

      deleted = 0;
    }
  }
  if (deleted) {
    sprintf (buf, "D%d ", deleted);
    StringVectorPushBack (cigar, buf);
  }

  if (strand > 0)
    for (n = 0; n < (int) VectorSize(cigar); ++n)
      fprintf (out, "%s", VectorGet(cigar,n));
  else
    for (n = ((int) VectorSize(cigar)) - 1; n >= 0; --n)
      fprintf (out, "%s", VectorGet(cigar,n));

  fprintf (out, "\n");

  deleteVector (cigar);
}

void append_seqrow_column_to_seqevt_alignment (StringVector *seqrow, char res) {
  char buf[] = "*";
  buf[0] = res;
  StringVectorPushBack (seqrow, buf);
}

void append_evtrow_column_to_seqevt_alignment (StringVector *seqrow, StringVector *evtrow, Fast5_event *evt) {
  char buf[] = "*", gap[] = "-";
  int len, n;
  if (evt == NULL) {
    while (StringVectorSize(seqrow) > StringVectorSize(evtrow))
      StringVectorPushBack (evtrow, gap);
  } else {
    len = strlen (evt->model_state);
    if (evt->move > 0) {
      for (n = 0; n < evt->move; ++n) {
	buf[0] = toupper (evt->model_state[len - evt->move + n]);
	StringVectorPushBack (evtrow, buf);
      }
    } else {
      buf[0] = tolower (evt->model_state[len-1]);
      StringVectorPushBack (evtrow, buf);
    }
    while (StringVectorSize(evtrow) > StringVectorSize(seqrow))
      StringVectorPushBack (seqrow, gap);
  }
}

void write_seq_event_pair_alignment_as_stockholm (Seq_event_pair_alignment* align, int strand, FILE* out) {
  StringVector *seqrow, *evtrow;
  int n_event, n, seqpos_offset, name_width;
  char namefmt[100], *revcomp;

  revcomp = strand > 0 ? NULL : new_revcomp_seq (align->seq, strlen(align->seq));

  seqrow = newStringVector();
  evtrow = newStringVector();

  n_event = 0;
  for (n = 0; n < align->start_n_event; ++n)
    append_evtrow_column_to_seqevt_alignment (seqrow, evtrow, &align->events->event[n_event++]);

  for (seqpos_offset = 0; seqpos_offset <= align->end_seqpos - align->start_seqpos; ++seqpos_offset) {
    append_seqrow_column_to_seqevt_alignment (seqrow, (strand > 0 ? align->seq : revcomp)[align->start_seqpos + seqpos_offset]);
    if (align->events_at_pos[seqpos_offset] == 0)
      append_evtrow_column_to_seqevt_alignment (seqrow, evtrow, NULL);
    else
      for (n = 0; n < align->events_at_pos[seqpos_offset]; ++n)
	append_evtrow_column_to_seqevt_alignment (seqrow, evtrow, &align->events->event[n_event++]);
  }

  name_width = 1 + MAX (strlen (align->seqname), strlen (align->events->name));
  sprintf (namefmt, "%%-%ds", name_width);

  fprintf (out, "# STOCKHOLM 1.0\n");
  fprintf (out, "#=GF SCORE %Lg\n", align->log_likelihood_ratio);
  fprintf (out, "#=GS %s STRAND %c\n", align->seqname, strand > 0 ? '+' : '-');
  fprintf (out, "#=GS %s START  %d\n", align->seqname, strand > 0 ? (align->start_seqpos + 1) : (align->seqlen - align->start_seqpos));
  fprintf (out, "#=GS %s END    %d\n", align->seqname, strand > 0 ? (align->end_seqpos + 1) : (align->seqlen - align->end_seqpos));

  fprintf (out, namefmt, align->seqname);
  for (n = 0; n < (int) StringVectorSize(seqrow); ++n)
    fprintf (out, "%s", StringVectorGet(seqrow,n));
  fprintf (out, "\n");

  fprintf (out, namefmt, align->events->name);
  for (n = 0; n < (int) StringVectorSize(evtrow); ++n)
    fprintf (out, "%s", StringVectorGet(evtrow,n));
  fprintf (out, "\n");

  fprintf (out, "//\n");

  deleteVector (seqrow);
  deleteVector (evtrow);
  SafeFreeOrNull (revcomp);
}

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_viterbi_matrix* mx;
  int n_events;
  unsigned long matrix_cells;

  mx = SafeMalloc (sizeof (Seq_event_pair_viterbi_matrix));
  mx->data = new_seq_event_pair_data (model, seqlen, seq, events);

  n_events = mx->data->events->n_events;
  matrix_cells = mx->data->matrix_cells;

  mx->vitStart = SafeMalloc ((n_events + 1) * sizeof(long double));
  mx->vitMatch = SafeMalloc (matrix_cells * sizeof(long double));
  mx->vitDelete = SafeMalloc (matrix_cells * sizeof(long double));

  return mx;
}

void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* mx) {
  SafeFree (mx->vitStart);
  SafeFree (mx->vitMatch);
  SafeFree (mx->vitDelete);
  SafeFree (mx->data);
  SafeFree (mx);
}

void update_max (long double *current_max, int *current_max_idx, long double candidate_max, int candidate_max_idx) {
  if (candidate_max > *current_max) {
    *current_max = candidate_max;
    *current_max_idx = candidate_max_idx;
  }
}

void fill_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* matrix) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event;
  long double mat, del;
  unsigned long idx, inputIdx, outputIdx;
  Fast5_event* event;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;

  /* update data */
  precalc_seq_event_pair_data (data);

  /* fill Viterbi */
  matrix->vitStart[0] = 0.;
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    matrix->vitStart[n_event]
      = matrix->vitStart[n_event - 1]
      + data->startEmitYes * event->ticks
      + data->nullEmitDensity[n_event];  /* Start -> Start (output) */
  }

  matrix->vitEnd = -INFINITY;

  for (seqpos = order; seqpos <= seqlen; ++seqpos) {
    for (n_event = 0; n_event <= n_events; ++n_event) {
      idx = Seq_event_pair_index(seqpos,n_event);

      mat = matrix->vitStart[n_event] + data->startEmitNo;   /* Start -> Match (input) */

      if (n_event > 0) {
	outputIdx = Seq_event_pair_index(seqpos,n_event-1);
	event = &data->events->event[n_event - 1];
	mat = max_func (mat,
			matrix->vitMatch[outputIdx]
			+ data->matchEmitYes[seqpos] * event->ticks
			+ data->matchEmitDensity[idx]);     /* Match -> Match (output) */
      }

      if (seqpos == order) {
	del = -INFINITY;
      } else {  /* seqpos > order */
	inputIdx = Seq_event_pair_index(seqpos-1,n_event);
	del = max_func
	  (matrix->vitMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteYes,  /* Match -> Delete (input) */
	   matrix->vitDelete[inputIdx] + data->extendDeleteYes);   /* Delete -> Delete (input) */

	mat = max_func
	  (mat,
	   matrix->vitMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteNo);  /* Match -> Match (input) */

	mat = max_func (mat, del + data->extendDeleteNo);  /* Delete -> Match */
      }

      matrix->vitMatch[idx] = mat;
      matrix->vitDelete[idx] = del;
    }

    matrix->vitEnd = max_func
      (matrix->vitEnd,
       matrix->vitMatch[Seq_event_pair_index(seqpos,n_events)] + data->matchEmitNo[seqpos]);  /* Match -> End (input) */
  }

#if defined(SEQEVTPAIR_DEBUG) && SEQEVTPAIR_DEBUG >= 10
  dump_seq_event_pair_matrix_to_file (SEQEVTMATRIX_FILENAME, "Viterbi", data, matrix->vitStart, matrix->vitMatch, matrix->vitDelete);
#endif /* SEQEVTPAIR_DEBUG >= 10 */
}

Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event, end_seqpos, start_seqpos, start_n_event, n, k;
  long double loglike;
  unsigned long idx, inputIdx, outputIdx;
  Fast5_event* event;
  enum { None, StartMatchIn, MatchMatchOut, MatchMatchIn, MatchDeleteIn, DeleteDeleteIn, DeleteMatch } trans;
  enum { Start, Match, Delete } state;
  Vector *events_emitted;
  Seq_event_pair_alignment *align;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;

  start_seqpos = seqlen;
  start_n_event = n_events;
  
  /* traceback from End state */
  end_seqpos = -1;
  loglike = -INFINITY;
  for (seqpos = seqlen; seqpos >= order; --seqpos)
    update_max (&loglike,
		&end_seqpos,
		matrix->vitMatch[Seq_event_pair_index(seqpos,n_events)] + data->matchEmitNo[seqpos],
		seqpos);
  Assert (end_seqpos >= 0, "Traceback failed");

  seqpos = end_seqpos;
  n_event = n_events;
  state = Match;

  events_emitted = newVector (IntCopy, IntDelete, IntPrint);
  VectorPushBack (events_emitted, IntNew(0));

  while (state != Start) {
    idx = Seq_event_pair_index(seqpos,n_event);
    if (seqpos > 0)
      inputIdx = Seq_event_pair_index(seqpos-1,n_event);
    if (n_event > 0)
      outputIdx = Seq_event_pair_index(seqpos,n_event-1);

    event = n_event > 0 ? &data->events->event[n_event - 1] : NULL;

    loglike = -INFINITY;
    trans = None;

    switch (state) {
    case Match:
      update_max (&loglike,
		  (int*) &trans,
		  matrix->vitStart[n_event] + data->startEmitNo,
		  StartMatchIn);   /* Start -> Match (input) */      

      if (n_event > 0)
	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitMatch[outputIdx]
		    + data->matchEmitYes[seqpos] * event->ticks
		    + data->matchEmitDensity[idx],
		    MatchMatchOut);    /* Match -> Match (output) */

      if (seqpos > order) {
	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteNo,
		    MatchMatchIn);    /* Match -> Match (input) */

	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitDelete[idx] + data->extendDeleteNo,
		    DeleteMatch);    /* Delete -> Match */
      }
      break;

    case Delete:
	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitMatch[inputIdx] + data->matchEmitNo[seqpos-1] + data->beginDeleteYes,
		    MatchDeleteIn);   /* Match -> Delete (input) */

	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitDelete[inputIdx] + data->extendDeleteYes,
		    DeleteDeleteIn);   /* Delete -> Delete (input) */

	break;

    default:
      Abort ("Unknown traceback state");
      break;
    }

    switch (trans) {

    case StartMatchIn:
      state = Start;
      start_n_event = n_event;
      start_seqpos = seqpos;
      break;

    case MatchMatchOut:
      Assert (state == Match, "oops");
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      break;

    case MatchMatchIn:
      Assert (state == Match, "oops");
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case MatchDeleteIn:
      state = Match;
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case DeleteDeleteIn:
      Assert (state == Delete, "oops");
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case DeleteMatch:
      state = Delete;
      break;

    case None:
    default:
      Abort ("Unknown traceback transition");
      break;
    }
  }

  align = new_seq_event_pair_alignment (data->events, data->seq, seqlen);
  align->log_likelihood_ratio = matrix->vitEnd - matrix->data->nullModel;
  align->start_seqpos = start_seqpos;
  align->end_seqpos = end_seqpos;
  align->start_n_event = start_n_event;
  align->events_at_pos = SafeMalloc (VectorSize(events_emitted) * sizeof(int));
  for (n = ((int) VectorSize(events_emitted)) - 1, k = 0; n >= 0; --n, ++k)
    align->events_at_pos[k] = *((int*) VectorGet (events_emitted, n));

  deleteVector (events_emitted);

  return align;
}

void print_seq_evt_pair_alignments_as_gff_cigar (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold) {
  print_seq_evt_pair_alignments_generic (model, seqlen, seq, seqname, events, out, log_odds_ratio_threshold, write_seq_event_pair_alignment_as_gff_cigar);
}

void print_seq_evt_pair_alignments_as_stockholm (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold) {
  print_seq_evt_pair_alignments_generic (model, seqlen, seq, seqname, events, out, log_odds_ratio_threshold, write_seq_event_pair_alignment_as_stockholm);
}
  
void print_seq_evt_pair_alignments_generic (Seq_event_pair_model* model, int seqlen, char *seq, char *seqname, Fast5_event_array* events, FILE *out, double log_odds_ratio_threshold, WriteSeqEventPairAlignmentFunction write_func) {
  Seq_event_pair_viterbi_matrix* matrix;
  Seq_event_pair_alignment* align;
  char *rev;
  int strand;

  rev = new_revcomp_seq (seq, seqlen);

  for (strand = +1; strand >= -1; strand -= 2) {
    matrix = new_seq_event_pair_viterbi_matrix (model, seqlen, strand > 0 ? seq : rev, events);
    fill_seq_event_pair_viterbi_matrix (matrix);

    align = get_seq_event_pair_viterbi_matrix_traceback (matrix);
    align->seqname = seqname;

    if (align->log_likelihood_ratio >= log_odds_ratio_threshold)
      (*write_func) (align, strand, out);

    delete_seq_event_pair_viterbi_matrix (matrix);
    delete_seq_event_pair_alignment (align);
  }

  SafeFree (rev);
}

void dump_seq_event_pair_matrix_to_file (const char* filename, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxDelete) {
  FILE *file;
  file = fopen (filename, "a");
  Assert (file != NULL, "Could not write DP matrix to file %s", filename);
  dump_seq_event_pair_matrix (file, algorithm, data, mxStart, mxMatch, mxDelete);
  fclose (file);
}

void dump_seq_event_pair_matrix (FILE* file, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double *mxDelete) {
  int n_event, seqpos, n_events, seqlen, order;
  Fast5_event *event;
  char *id;
  unsigned long idx;
  time_t rawtime;
  struct tm *rawtime_tm;

  order = data->model->order;
  n_events = data->events->n_events;
  seqlen = data->seqlen;

  id = SafeMalloc ((data->model->order + 1) * sizeof(char));

  time (&rawtime);
  rawtime_tm = localtime (&rawtime);
  fprintf (file, "%s matrix (%s)\n", algorithm, asctime(rawtime_tm));

  for (seqpos = order; seqpos <= seqlen; ++seqpos) {
    encode_state_identifier (data->state[seqpos], data->model->order, id);
    idx = Seq_event_pair_index(seqpos,0);
    fprintf (file, "Event %d, seqpos %d (base=%c,state=%s): match %Lg, delete %Lg\n", 0, seqpos, data->seq[seqpos-1], id, mxMatch[idx], mxDelete[idx]);
  }
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    fprintf (file, "Event %d (n=%g,sum=%g,sumsq=%g): start %Lg\n", n_event, event->ticks, event->sumticks_cur, event->sumticks_cur_sq, mxStart[n_event]);
    for (seqpos = order; seqpos <= seqlen; ++seqpos) {
      encode_state_identifier (data->state[seqpos], data->model->order, id);
      idx = Seq_event_pair_index(seqpos,n_event);
      fprintf (file, "Event %d (n=%g,sum=%g,sumsq=%g), seqpos %d (base=%c,state=%s): match %Lg, delete %Lg\n", n_event, event->ticks, event->sumticks_cur, event->sumticks_cur_sq, seqpos, data->seq[seqpos-1], id, mxMatch[idx], mxDelete[idx]);
    }
  }

  SafeFree (id);
}

void xmlTextWriterBooleanCount (xmlTextWriterPtr writer, const char* element, double yes, double no) {
  xmlTextWriterStartElement (writer, (xmlChar*) element);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(YES), "%g", yes);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(NO), "%g", no);
  xmlTextWriterEndElement (writer);
}

xmlChar* convert_seq_event_pair_counts_to_xml_string (Seq_event_pair_counts* counts) {
  xmlTextWriterPtr writer;
  char* id;
  int state;

  id = SafeMalloc ((counts->order + 1) * sizeof(char));

  writer = newXmlTextWriter();
  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(COUNTS));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(ORDER), "%d", counts->order);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(DELETE));
  xmlTextWriterBooleanCount (writer, XMLPREFIX(BEGIN), counts->nBeginDeleteYes, counts->nBeginDeleteNo);
  xmlTextWriterBooleanCount (writer, XMLPREFIX(EXTEND), counts->nExtendDeleteYes, counts->nExtendDeleteNo);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATES));
  for (state = 0; state < counts->states; ++state) {
    xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATE));
    encode_state_identifier (state, counts->order, id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(KMER), "%s", id);
    xmlTextWriterBooleanCount (writer, XMLPREFIX(WAIT), counts->nMatchEmitYes[state], counts->nMatchEmitNo[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M0), "%g", counts->matchMoment0[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M1), "%g", counts->matchMoment1[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M2), "%g", counts->matchMoment2[state]);
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(START));
  xmlTextWriterBooleanCount (writer, XMLPREFIX(WAIT), counts->nStartEmitYes, counts->nStartEmitNo);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(NULLMODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M0), "%g", counts->nullMoment0);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M1), "%g", counts->nullMoment1);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M2), "%g", counts->nullMoment2);
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);

  SafeFree (id);

  return deleteXmlTextWriterLeavingText (writer);
}
