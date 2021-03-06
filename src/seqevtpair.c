#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <float.h>

#include <hdf5.h>
#include <hdf5_hl.h>

#include <gsl/gsl_randist.h>

#include "seqevtpair.h"
#include "xmlutil.h"
#include "xmlkeywords.h"
#include "kseqcontainer.h"
#include "logsumexp.h"

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

/* helpers for finding index of max item in an indexed set (for Viterbi traceback) */
void update_max (long double *current_max, int* current_max_idx, long double candidate_max, int candidate_max_idx);

/* squiggle plot config */
const double pixels_per_tick = .02;
const double pixels_per_lev = 10;

/* normalize a model */
void normalize_model (Seq_event_pair_model* model);

/* accum_count(back_src,fwd_src,trans,back_dest,matrix,count,event,moment0,moment1,moment2)
   increments *count by (weight = exp(fwd_src + trans + back_dest - matrix->fwdResult))
   also increments (moment0,moment1,moment2) by weight * event->(ticks,sumticks_cur,sumticks_cur_sq)
   returns log_sum_exp(back_src,trans + back_dest)
 */
long double accum_count (long double back_src,
			 long double fwd_src,
			 long double trans,
			 long double back_dest,
			 Seq_event_pair_fb_matrix *matrix,
			 long double *count,
			 Fast5_event *event,
			 long double *moment0,
			 long double *moment1,
			 long double *moment2);

/* Metrichor_state_iterator
   Used to populate the emission parameters of a Seq_event_pair_model from a Metrichor model description in a FAST5 file */
typedef struct Metrichor_state_iterator {
  Seq_event_pair_model *model;
  size_t kmer_offset, level_mean_offset, level_sd_offset;
} Metrichor_state_iterator;

double read_metrichor_model_attribute (hid_t dataset, const char* attr_name)
{
  herr_t ret;
  hid_t attr;
  double attr_val;
  attr = H5Aopen_name(dataset,attr_name);
  ret = H5Aread(attr, H5T_NATIVE_DOUBLE, &attr_val);
  ret = H5Aclose(attr);
  return attr_val;
}

herr_t populate_seq_event_model_emit_params (void *elem, hid_t type_id, unsigned ndim, 
					     const hsize_t *point, void *operator_data)
{
  Metrichor_state_iterator *iter = (Metrichor_state_iterator*) operator_data;
  int state = decode_state_identifier (iter->model->order, (char*) elem + iter->kmer_offset);
  iter->model->matchMean[state] = *((double*) (elem + iter->level_mean_offset));
  double sd = MAX (DBL_MIN, *((double*) (elem + iter->level_sd_offset)));
  iter->model->matchPrecision[state] = 1 / (sd * sd);
  return 0;
}

/* main function bodies */
void init_seq_event_pair_config (Seq_event_pair_config *config) {
  /* EM convergence criteria */
  config->seq_evt_pair_EM_max_iterations = 100;
  config->seq_evt_pair_EM_min_fractional_loglike_increment = 0.0001;
  config->debug_matrix_filename = "npmatrix";
  config->both_strands = 0;
  config->all_vs_all = 0;
}

int parse_seq_event_pair_config_general (int* argcPtr, char*** argvPtr, Seq_event_pair_config *config) {
  if (*argcPtr > 0) {
    const char* arg = **argvPtr;
    if (strcmp (arg, "-matrixdumpfile") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      const char* val = (*argvPtr)[1];
      config->debug_matrix_filename = val;
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;

    } else if (strcmp (arg, "-bothstrands") == 0) {
      config->both_strands = 1;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    }

  }
  return 0;
}

int parse_seq_event_pair_config_training (int* argcPtr, char*** argvPtr, Seq_event_pair_config *config) {
  if (*argcPtr > 0) {
    const char* arg = **argvPtr;
    if (strcmp (arg, "-maxiter") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      const char* val = (*argvPtr)[1];
      config->seq_evt_pair_EM_max_iterations = atoi (val);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;

    } else if (strcmp (arg, "-mininc") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      const char* val = (*argvPtr)[1];
      config->seq_evt_pair_EM_min_fractional_loglike_increment = atof (val);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;

    } else if (strcmp (arg, "-allvsall") == 0) {
      config->all_vs_all = 1;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    }

  }

  return parse_seq_event_pair_config_general (argcPtr, argvPtr, config);
}

StringDoubleMap* new_seq_event_pair_model_default_prior() {
  StringDoubleMap *prior = newStringDoubleMap();
  StringDoubleMapSet (prior, "no_skip", 9);
  StringDoubleMapSet (prior, "no_delete", 99);
  return prior;
}

Seq_event_pair_model* new_seq_event_pair_model (int order) {
  StringDoubleMap* prior = new_seq_event_pair_model_default_prior();
  Seq_event_pair_model* model = new_seq_event_pair_model_from_prior (order, prior);
  deleteStringDoubleMap (prior);
  return model;
}

Seq_event_pair_model* new_seq_event_pair_model_from_prior (int order, StringDoubleMap *prior) {
  Seq_event_pair_model *model;
  model = SafeMalloc (sizeof (Seq_event_pair_model));
  model->order = order;
  model->states = pow(AlphabetSize,order);
  model->pMatchEvent = SafeMalloc (model->states * sizeof(double));
  model->pMatchTick = SafeMalloc (model->states * sizeof(double));
  model->pMatchSkip = SafeMalloc (model->states * sizeof(double));
  model->matchMean = SafeMalloc (model->states * sizeof(double));
  model->matchPrecision = SafeMalloc (model->states * sizeof(double));
  model->kmerProb = SafeMalloc (model->states * sizeof(double));

  model->prior = RBTreeDeepCopy (prior);
  
  model->pBeginDelete = get_seq_event_prior_mode(model,"delete");
  model->pExtendDelete = get_seq_event_prior_mode(model,"extend");
  model->pStartEvent = get_seq_event_prior_mode(model,"event");
  model->pNullEvent = get_seq_event_prior_mode(model,"event");
  model->pNullTick = get_seq_event_prior_mode(model,"tick");
  model->emitProb = get_seq_event_prior_mode(model,"emit");
  model->nullMean = 0;
  model->nullPrecision = 1;
  for (int state = 0; state < model->states; ++state) {
    model->pMatchEvent[state] = get_seq_event_prior_mode(model,"event");
    model->pMatchTick[state] = get_seq_event_prior_mode(model,"tick");
    model->pMatchSkip[state] = get_seq_event_prior_mode(model,"skip");
    model->matchMean[state] = 0;
    model->matchPrecision[state] = 1;
    model->kmerProb[state] = 1. / (double) model->states;
  }

  model->logger = NULL;
  init_seq_event_pair_config (&model->config);
  
  return model;
}

void delete_seq_event_pair_model (Seq_event_pair_model* model) {
  deleteStringDoubleMap (model->prior);
  SafeFree (model->kmerProb);
  SafeFree (model->pMatchEvent);
  SafeFree (model->pMatchTick);
  SafeFree (model->matchMean);
  SafeFree (model->matchPrecision);
  SafeFree (model);
}

double get_seq_event_pair_pseudocount (Seq_event_pair_model* model, const char* param) {
  StringDoubleMapNode *node = StringDoubleMapFind (model->prior, param);
  return node ? *((double*)node->value) : 1.;
}

double get_seq_event_prior_mode (Seq_event_pair_model* model, const char* param) {
  char *no_param = SafeMalloc ((strlen(param) + 4) * sizeof(char));
  sprintf (no_param, "no_%s", param);
  double p = BetaMode (get_seq_event_pair_pseudocount(model,param), get_seq_event_pair_pseudocount(model,no_param));
  SafeFree (no_param);
  return p;
}

double get_seq_event_log_beta_prior (Seq_event_pair_model* model, const char* param, double x) {
  char *no_param = SafeMalloc ((strlen(param) + 4) * sizeof(char));
  sprintf (no_param, "no_%s", param);
  double lp = log (gsl_ran_beta_pdf (x,
				     get_seq_event_pair_pseudocount(model,param) + 1,
				     get_seq_event_pair_pseudocount(model,no_param) + 1));
  SafeFree (no_param);
  return lp;
}

int base2token (char base) {
  return tokenize (toupper(base), dna_alphabet);
}

char token2base (int tok) {
  return tok < 0 || tok >= AlphabetSize ? 'N' : dna_alphabet[tok];
}

void encode_state_identifier (int state, int order, char* state_id) {
  int k;
  for (k = 0; k < order; ++k, state = state < 0 ? state : state / AlphabetSize)
    state_id[order - k - 1] = token2base (state < 0 ? -1 : (state % AlphabetSize));
  state_id[order] = '\0';
}

int decode_state_identifier (int order, char* state_id) {
  int k, token, chartok, mul;
  for (token = 0, mul = 1, k = 0; k < order; ++k, mul *= AlphabetSize) {
    char c = state_id[order - k - 1];
    chartok = base2token (c);
    Assert (chartok >= 0, "Unknown character in k-mer: %c", c);
    token += mul * chartok;
  }
  return token;
}

double emitProbToMeanLength (double p) { return p / (1. - p); }
double meanLengthToEventProb (double l) { return l / (1. + l); }

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml) {
  xmlNode *modelNode, *statesNode, *stateNode, *deleteNode, *startNode, *nullNode, *pseudoNode;
  Seq_event_pair_model *model;
  int state;
  modelNode = xmlTreeFromString (xml);
  model = new_seq_event_pair_model (CHILDINT(modelNode,ORDER));

  for (pseudoNode = modelNode->children; pseudoNode; pseudoNode = pseudoNode->next)
    if (MATCHES(pseudoNode,PSEUDO)) {
      const char* param = (const char*) CHILDSTRING(pseudoNode,PARAM);
      const double count = CHILDFLOAT(pseudoNode,COUNT);
      StringDoubleMapSet (model->prior, param, count);
    }

  deleteNode = CHILD(modelNode,DELETE);
  model->pBeginDelete = CHILDFLOAT(deleteNode,BEGIN);
  model->pExtendDelete = CHILDFLOAT(deleteNode,EXTEND);

  statesNode = CHILD(modelNode,STATES);
  for (stateNode = statesNode->children; stateNode; stateNode = stateNode->next)
    if (MATCHES(stateNode,STATE)) {
      state = decode_state_identifier (model->order, (char*) CHILDSTRING(stateNode,KMER));
      model->pMatchEvent[state] = meanLengthToEventProb (CHILDFLOAT(stateNode,EVENTS));
      model->pMatchTick[state] = meanLengthToEventProb (CHILDFLOAT(stateNode,TICKS));
      model->pMatchSkip[state] = CHILDFLOAT(stateNode,SKIP);
      model->matchMean[state] = CHILDFLOAT(stateNode,MEAN);
      model->matchPrecision[state] = 1 / MAX (DBL_MIN, pow (CHILDFLOAT(stateNode,STDV), 2));
      model->kmerProb[state] = MAX (DBL_MIN, CHILDFLOAT(stateNode,FREQ));
    }

  model->emitProb = meanLengthToEventProb (CHILDFLOAT(modelNode,SEQLEN));
  
  startNode = CHILD(modelNode,START);
  model->pStartEvent = meanLengthToEventProb (CHILDFLOAT(startNode,EVENTS));

  nullNode = CHILD(modelNode,NULLMODEL);
  model->pNullEvent = meanLengthToEventProb (CHILDFLOAT(nullNode,EVENTS));
  model->pNullTick = meanLengthToEventProb (CHILDFLOAT(nullNode,TICKS));
  model->nullMean = CHILDFLOAT(nullNode,MEAN);
  model->nullPrecision = 1 / pow (MAX (DBL_MIN, CHILDFLOAT(nullNode,STDV)), 2);

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

  for (StringDoubleMapNode *pseudoNode = RBTreeFirst(model->prior);
       !RBTreeIteratorFinished(model->prior,pseudoNode);
       pseudoNode = RBTreeSuccessor(model->prior,pseudoNode)) {
    xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(PSEUDO));
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(PARAM), "%s", (char*) pseudoNode->key);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(COUNT), "%g", *(double*)pseudoNode->value);
    xmlTextWriterEndElement (writer);
  }

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(DELETE));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(BEGIN), "%g", model->pBeginDelete);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EXTEND), "%g", model->pExtendDelete);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATES));
  for (state = 0; state < model->states; ++state) {
    xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(STATE));
    encode_state_identifier (state, model->order, id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(KMER), "%s", id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(SKIP), "%g", model->pMatchSkip[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EVENTS), "%g", emitProbToMeanLength (model->pMatchEvent[state]));
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(TICKS), "%g", emitProbToMeanLength (model->pMatchTick[state]));
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->matchMean[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(STDV), "%g", 1 / sqrt(model->matchPrecision[state]));
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(FREQ), "%g", model->kmerProb[state]);
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(SEQLEN), "%g", emitProbToMeanLength (model->emitProb));

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(START));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EVENTS), "%g", emitProbToMeanLength (model->pStartEvent));
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(NULLMODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EVENTS), "%g", emitProbToMeanLength (model->pNullEvent));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(TICKS), "%g", emitProbToMeanLength (model->pNullTick));
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

  Logger *logger = model->logger;

  Assert (seqlen >= model->order, "Sequence length is %d, which is less than model order (%d)", seqlen, model->order);

  data = SafeMalloc (sizeof (Seq_event_pair_data));
  data->model = model;
  data->seqlen = seqlen;
  data->seq = seq;
  data->events = events;

  n_events = events->n_events;

  matrix_cells = ((unsigned long) (n_events + 1)) * (unsigned long) (seqlen - model->order + 1);
  if (LogThisAt(1))
    Warn ("Allocating DP scratch space with %d*%d = %Lu cells", n_events + 1, seqlen - model->order + 1, matrix_cells);

  data->matrix_cells = matrix_cells;

  data->nullEventLogLike = SafeMalloc ((n_events + 1) * sizeof(long double));

  data->matchEventLogLike = SafeMalloc (matrix_cells * sizeof(long double));
  data->matchEventYes = SafeMalloc ((seqlen + 1) * sizeof(long double));
  data->matchEventNo = SafeMalloc ((seqlen + 1) * sizeof(long double));
  data->matchSkipYes = SafeMalloc ((seqlen + 1) * sizeof(long double));
  data->matchSkipNo = SafeMalloc ((seqlen + 1) * sizeof(long double));

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
  SafeFree (data->nullEventLogLike);
  SafeFree (data->matchEventLogLike);
  SafeFree (data->matchEventYes);
  SafeFree (data->matchEventNo);
  SafeFree (data->matchSkipYes);
  SafeFree (data->matchSkipNo);
  SafeFree (data);
}

void precalc_seq_event_pair_data (Seq_event_pair_data* data) {
  Seq_event_pair_model* model;
  Logger *logger;
  int seqlen, n_events, order, seqpos, n_event, state;
  long double logNullPrecision, logNullTick, logNoNullTick, loglike;
  double mean, precision, logPrecision, logTick, logNoTick;
  Fast5_event* event;

  model = data->model;
  order = model->order;
  seqlen = data->seqlen;
  n_events = data->events->n_events;

  logger = model->logger;

  /* calculate logs of transition probabilities & emit precisions */
  data->nullEventYes = log (model->pNullEvent);
  data->nullEventNo = log (1. - model->pNullEvent);

  data->startEventYes = log (model->pStartEvent);
  data->startEventNo = log (1. - model->pStartEvent);

  data->beginDeleteYes = log (model->pBeginDelete);
  data->beginDeleteNo = log (1. - model->pBeginDelete);

  data->extendDeleteYes = log (model->pExtendDelete);
  data->extendDeleteNo = log (1. - model->pExtendDelete);

  logNullPrecision = log (model->nullPrecision);
  logNullTick = log (model->pNullTick);
  logNoNullTick = log (1. - model->pNullTick);

  /* calculate emit densities & state-dependent transition probabilities */
  data->nullEventLogLike[0] = -INFINITY;
  data->nullModel = data->nullEventNo;
  for (n_event = 0; n_event < n_events; ++n_event) {
    event = &data->events->event[n_event];
    loglike = log_event_density (event,
				 model->nullMean,
				 model->nullPrecision,
				 logNullPrecision,
				 logNullTick,
				 logNoNullTick);
    data->nullEventLogLike[n_event + 1] = loglike;
    data->nullModel += loglike + data->nullEventYes;
  }

  for (seqpos = 0; seqpos < order; ++seqpos) {
    data->matchEventYes[seqpos] = -INFINITY;
    data->matchEventNo[seqpos] = -INFINITY;
  }

  if (LogThisAt(2))
    init_progress ("Pre-DP calculations");
  for (seqpos = order; seqpos <= seqlen; ++seqpos) {

    if (LogThisAt(1))
      log_progress ((seqpos - order) / (double) (seqlen - order), "base %d/%d", seqpos, seqlen);

    state = data->state[seqpos];

    data->matchEventYes[seqpos] = state < 0 ? data->nullEventYes : log (model->pMatchEvent[state]);
    data->matchEventNo[seqpos] = state < 0 ? data->nullEventNo : log (1. - model->pMatchEvent[state]);

    data->matchSkipYes[seqpos] = log (model->pMatchSkip[state]);
    data->matchSkipNo[seqpos] = log (1. - model->pMatchSkip[state]);

    mean = state < 0 ? model->nullMean : model->matchMean[state];
    precision = state < 0 ? model->nullPrecision : model->matchPrecision[state];
    logPrecision = log (precision);
    logTick = log (model->pMatchTick[state]);
    logNoTick = log (1. - model->pMatchTick[state]);
    
    data->matchEventLogLike[Seq_event_pair_index(seqpos,0)] = -INFINITY;
    for (n_event = 1; n_event <= n_events; ++n_event) {
      event = &data->events->event[n_event - 1];
      data->matchEventLogLike[Seq_event_pair_index(seqpos,n_event)]
	= log_event_density (event,
			     mean,
			     precision,
			     logPrecision,
			     logTick,
			     logNoTick);
    }
  }

  /* calculate parameter prior */
  data->logPrior = data->logNullPrior = 0;
  data->logNullPrior += get_seq_event_log_beta_prior (model, "event", model->pNullEvent);
  data->logPrior += get_seq_event_log_beta_prior (model, "tick", model->pNullTick);
  data->logPrior += get_seq_event_log_beta_prior (model, "event", model->pStartEvent);
  for (state = 0; state < model->states; ++state) {
    data->logPrior += get_seq_event_log_beta_prior (model, "event", model->pMatchEvent[state]);
    data->logPrior += get_seq_event_log_beta_prior (model, "tick", model->pMatchTick[state]);
    data->logPrior += get_seq_event_log_beta_prior (model, "skip", model->pMatchSkip[state]);
  }
  data->logPrior += get_seq_event_log_beta_prior (model, "delete", model->pBeginDelete);
  data->logPrior += get_seq_event_log_beta_prior (model, "extend", model->pExtendDelete);
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
  Logger *logger;
  Seq_event_pair_fb_matrix* mx;
  int n_events;
  unsigned long matrix_cells;

  logger = model->logger;

  mx = SafeMalloc (sizeof (Seq_event_pair_fb_matrix));
  mx->data = new_seq_event_pair_data (model, seqlen, seq, events);

  n_events = mx->data->events->n_events;
  matrix_cells = mx->data->matrix_cells;

  if (LogThisAt(1))
    Warn ("Allocating Forward-Backward matrix of size %d*%d (approx.)", n_events, seqlen);

  mx->fwdStart = SafeMalloc ((n_events + 1) * sizeof(long double));
  mx->backStart = SafeMalloc ((n_events + 1) * sizeof(long double));

  mx->fwdMatch = SafeMalloc (matrix_cells * sizeof(long double));
  mx->fwdSkip = SafeMalloc (matrix_cells * sizeof(long double));
  mx->fwdDelete = SafeMalloc (matrix_cells * sizeof(long double));
  mx->backMatch = SafeMalloc (matrix_cells * sizeof(long double));
  mx->backSkip = SafeMalloc (matrix_cells * sizeof(long double));
  mx->backDelete = SafeMalloc (matrix_cells * sizeof(long double));

  return mx;
}

void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* mx) {
  delete_seq_event_pair_data (mx->data);
  SafeFree (mx->fwdStart);
  SafeFree (mx->backStart);
  SafeFree (mx->fwdMatch);
  SafeFree (mx->fwdSkip);
  SafeFree (mx->fwdDelete);
  SafeFree (mx->backMatch);
  SafeFree (mx->backSkip);
  SafeFree (mx->backDelete);
  SafeFree (mx);
}

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model) {
  Seq_event_pair_counts* counts;
  counts = SafeMalloc (sizeof (Seq_event_pair_counts));
  counts->order = model->order;
  counts->states = model->states;
  counts->nMatchEventYes = SafeMalloc (model->states * sizeof(long double));
  counts->nMatchEventNo = SafeMalloc (model->states * sizeof(long double));
  counts->nMatchSkipYes = SafeMalloc (model->states * sizeof(long double));
  counts->nMatchSkipNo = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment0 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment1 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment2 = SafeMalloc (model->states * sizeof(long double));
  return counts;
}

void delete_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  SafeFree (counts->nMatchEventYes);
  SafeFree (counts->nMatchEventNo);
  SafeFree (counts->nMatchSkipYes);
  SafeFree (counts->nMatchSkipNo);
  SafeFree (counts->matchMoment0);
  SafeFree (counts->matchMoment1);
  SafeFree (counts->matchMoment2);
  SafeFree (counts);
}

void reset_seq_event_null_counts (Seq_event_pair_counts* counts) {
  counts->nNullEventYes = 0;
  counts->nNullEventNo = 0;
  counts->nullMoment0 = 0;
  counts->nullMoment1 = 0;
  counts->nullMoment2 = 0;
}

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  int state;
  for (state = 0; state < counts->states; ++state) {
    counts->nMatchEventYes[state] = 0;
    counts->nMatchEventNo[state] = 0;
    counts->nMatchSkipYes[state] = 0;
    counts->nMatchSkipNo[state] = 0;
    counts->matchMoment0[state] = 0;
    counts->matchMoment1[state] = 0;
    counts->matchMoment2[state] = 0;
  }
  counts->nStartEventYes = 0;
  counts->nStartEventNo = 0;
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

  double pseudoEventYes = get_seq_event_pair_pseudocount (model, "event");
  double pseudoEventNo = get_seq_event_pair_pseudocount (model, "no_event");
  double pseudoSkipYes = get_seq_event_pair_pseudocount (model, "skip");
  double pseudoSkipNo = get_seq_event_pair_pseudocount (model, "no_skip");
  for (state = 0; state < counts->states; ++state) {
    counts->nMatchEventYes[state] += pseudoEventYes;
    counts->nMatchEventNo[state] += pseudoEventNo;
    counts->nMatchSkipYes[state] += pseudoSkipYes;
    counts->nMatchSkipNo[state] += pseudoSkipNo;
    counts->matchMoment0[state] += 1.;
  }
  counts->nStartEventYes += pseudoEventYes;
  counts->nStartEventNo += pseudoEventNo;
  counts->nBeginDeleteYes += get_seq_event_pair_pseudocount (model, "delete");
  counts->nBeginDeleteNo += get_seq_event_pair_pseudocount (model, "no_delete");
  counts->nExtendDeleteYes += get_seq_event_pair_pseudocount (model, "extend");
  counts->nExtendDeleteNo += get_seq_event_pair_pseudocount (model, "no_extend");
  counts->nNullEventYes += pseudoEventYes;
  counts->nNullEventNo += pseudoEventNo;
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
      counts->nStartEventYes += 1.;
      counts->nullMoment0 += event->ticks;
      counts->nullMoment1 += event->sumticks_cur;
      counts->nullMoment2 += event->sumticks_cur_sq;
    } else {
      state = decode_state_identifier (counts->order, event->model_state);
      if (moves < counts->order)
	counts->nStartEventNo += 1.;
      counts->nMatchEventNo[state] += event->move;
      counts->nMatchEventYes[state] += event->ticks;
      counts->matchMoment0[state] += event->ticks;
      counts->matchMoment1[state] += event->sumticks_cur;
      counts->matchMoment2[state] += event->sumticks_cur_sq;
    }
    counts->nNullEventYes += 1.;
    moves = new_moves;
  }
  counts->nNullEventNo += 1.;
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
      counts->nMatchEventYes[state] += weight * inc[n]->nMatchEventYes[state];
      counts->nMatchEventNo[state] += weight * inc[n]->nMatchEventNo[state];
      counts->nMatchSkipYes[state] += weight * inc[n]->nMatchSkipYes[state];
      counts->nMatchSkipNo[state] += weight * inc[n]->nMatchSkipNo[state];
      counts->matchMoment0[state] += weight * inc[n]->matchMoment0[state];
      counts->matchMoment1[state] += weight * inc[n]->matchMoment1[state];
      counts->matchMoment2[state] += weight * inc[n]->matchMoment2[state];
    }
    counts->nStartEventYes += weight * inc[n]->nStartEventYes;
    counts->nStartEventNo += weight * inc[n]->nStartEventNo;
    counts->nBeginDeleteYes += weight * inc[n]->nBeginDeleteYes;
    counts->nBeginDeleteNo += weight * inc[n]->nBeginDeleteNo;
    counts->nExtendDeleteYes += weight * inc[n]->nExtendDeleteYes;
    counts->nExtendDeleteNo += weight * inc[n]->nExtendDeleteNo;
  }
  counts->loglike += loglike_all;
}

double log_event_density (Fast5_event* event, double mean, double precision, double log_precision, double log_pTick, double log_pNoTick) {
  return event->ticks * (log_pTick + log_precision/2. - log_sqrt2pi - precision*mean*mean/2)
    - precision*(event->sumticks_cur_sq/2. - event->sumticks_cur*mean)
    + log_pNoTick;
}

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event, state, nextState;
  long double mat, del, skp, st, count;
  unsigned long idx, inputIdx, outputIdx, ioIdx;
  Fast5_event *event;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;

  Logger *logger = model->logger;

  /* update data */
  precalc_seq_event_pair_data (data);

  /* fill forward */
  matrix->fwdStart[0] = 0.;
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    matrix->fwdStart[n_event]
      = matrix->fwdStart[n_event - 1]
      + data->startEventYes
      + data->nullEventLogLike[n_event];  /* Start -> Start (output) */
  }

  matrix->fwdResult = matrix->backResult = -INFINITY;

  if (LogThisAt(2))
    init_progress ("Forward matrix");
  for (seqpos = order; seqpos <= seqlen; ++seqpos) {

    if (LogThisAt(2))
      log_progress ((seqpos - order) / (double) (seqlen - order), "base %d/%d", seqpos, seqlen);

    for (n_event = 0; n_event <= n_events; ++n_event) {

      idx = Seq_event_pair_index(seqpos,n_event);
      mat = -INFINITY;
      
      if (n_event > 0) {
	outputIdx = Seq_event_pair_index(seqpos,n_event-1);
	event = &data->events->event[n_event - 1];
	mat = log_sum_exp (mat,
			   matrix->fwdMatch[outputIdx]
			   + data->matchEventYes[seqpos]
			   + data->matchEventLogLike[idx]);     /* Match -> Match (output) */

	mat = log_sum_exp (mat,
			   matrix->fwdStart[n_event-1]
			   + data->startEventNo
			   + data->matchEventLogLike[idx]);     /* Start -> Match (output) */
      } else
	event = NULL;

      if (seqpos == order) {
	skp = del = -INFINITY;
      } else {  /* seqpos > order */
	inputIdx = Seq_event_pair_index(seqpos-1,n_event);
	del = log_sum_exp
	  (matrix->fwdMatch[inputIdx]
	   + data->matchEventNo[seqpos-1]
	   + data->beginDeleteYes,     /* Match -> Delete (input) */

	   matrix->fwdDelete[inputIdx]
	   + data->extendDeleteYes);   /* Delete -> Delete (input) */

	skp = log_sum_exp (matrix->fwdMatch[inputIdx]
			   + data->matchEventNo[seqpos-1]
			   + data->beginDeleteNo
			   + data->matchSkipYes[seqpos],   /* Match -> Skip (input) */

			   matrix->fwdSkip[inputIdx]
			   + data->matchSkipYes[seqpos]);  /* Skip -> Skip (input) */

	if (n_event > 0) {
	  ioIdx = Seq_event_pair_index(seqpos-1,n_event-1);
	  mat = log_sum_exp (mat,
			     matrix->fwdMatch[ioIdx]
			     + data->matchEventNo[seqpos-1]
			     + data->beginDeleteNo
			     + data->matchSkipNo[seqpos]
			     + data->matchEventLogLike[idx]);     /* Match -> Match (input/output) */

	  mat = log_sum_exp (mat,
			     matrix->fwdSkip[ioIdx]
			     + data->matchSkipNo[seqpos]
			     + data->matchEventLogLike[idx]);     /* Skip -> Match (input/output) */

	  mat = log_sum_exp (mat,
			     matrix->fwdDelete[ioIdx]
			     + data->extendDeleteNo
			     + data->matchEventLogLike[idx]);     /* Delete -> Match (input/output) */
	}
      }

      matrix->fwdMatch[idx] = mat;
      matrix->fwdSkip[idx] = skp;
      matrix->fwdDelete[idx] = del;
    }

    matrix->fwdResult = log_sum_exp
      (matrix->fwdResult,
       matrix->fwdMatch[Seq_event_pair_index(seqpos,n_events)]
       + data->matchEventNo[seqpos]);  /* Match -> End (input) */
  }

  if (LogThisAt(2))
    fprintf (stderr, "Forward log-likelihood is %Lg (prior %Lg, null model %Lg, null prior %Lg)\n", matrix->fwdResult, matrix->data->logPrior, matrix->data->nullModel, matrix->data->logNullPrior);

  if (LogThisAt(10))
    dump_seq_event_pair_matrix_to_file (model->config.debug_matrix_filename, "Forward", data, matrix->fwdStart, matrix->fwdMatch, matrix->fwdSkip, matrix->fwdDelete, matrix->fwdResult);

  /* fill backward & accumulate counts */
  for (n_event = n_events; n_event >= 0; --n_event)
    matrix->backStart[n_event] = -INFINITY;

  if (LogThisAt(2))
    init_progress ("Backward matrix");
  for (seqpos = seqlen; seqpos >= order; --seqpos) {

    if (LogThisAt(2))
      log_progress ((seqlen - seqpos) / (double) (seqlen - order), "base %d/%d", seqpos, seqlen);

    state = data->state[seqpos];

    for (n_event = n_events; n_event >= 0; --n_event) {
      event = n_event < n_events ? &data->events->event[n_event] : NULL;
      idx = Seq_event_pair_index(seqpos,n_event);

      if (n_event == n_events) {
	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEventNo[seqpos],
			   0.,
			   matrix,
			   &count,
			   NULL, NULL, NULL, NULL);   /* Match -> End */
	counts->nMatchEventNo[state] += count;

      } else {  /* n_event < n_events */
	outputIdx = Seq_event_pair_index(seqpos,n_event+1);

	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEventYes[seqpos]
			   + data->matchEventLogLike[outputIdx],
			   matrix->backMatch[outputIdx],
			   matrix,
			   &count,
			   event,
			   &counts->matchMoment0[state],
			   &counts->matchMoment1[state],
			   &counts->matchMoment2[state]);  /* Match -> Match (output) */
	counts->nMatchEventYes[state] += count;

	matrix->backStart[n_event] = accum_count (matrix->backStart[n_event],
						  matrix->fwdStart[n_event],
						  data->startEventNo
						  + data->matchEventLogLike[outputIdx],
						  matrix->backMatch[outputIdx],
						  matrix,
						  &count,
						  NULL, NULL, NULL, NULL);   /* Start -> Match (output) */
	counts->nStartEventNo += count;
      }

      if (seqpos == seqlen) {
	skp = del = -INFINITY;
      } else {
	inputIdx = Seq_event_pair_index(seqpos+1,n_event);
	nextState = data->state[seqpos+1];

	del = accum_count (-INFINITY,
			   matrix->fwdDelete[idx],
			   data->extendDeleteYes,
			   matrix->backDelete[inputIdx],
			   matrix,
			   &count,
			   NULL, NULL, NULL, NULL);   /* Delete -> Delete (input) */
	counts->nExtendDeleteYes += count;
	
	mat = accum_count (mat,
			   matrix->fwdMatch[idx],
			   data->matchEventNo[seqpos]
			   + data->beginDeleteYes,
			   matrix->backDelete[inputIdx],
			   matrix,
			   &count,
			   NULL, NULL, NULL, NULL);  /* Match -> Delete (input) */
	counts->nMatchEventNo[state] += count;
	counts->nBeginDeleteYes += count;

	mat = accum_count (mat,
			   matrix->fwdMatch[idx],
			   data->matchEventNo[seqpos]
			   + data->beginDeleteNo
			   + data->matchSkipYes[seqpos+1],
			   matrix->backSkip[inputIdx],
			   matrix,
			   &count,
			   NULL, NULL, NULL, NULL);  /* Match -> Skip (input) */
	counts->nMatchEventNo[state] += count;
	counts->nBeginDeleteNo += count;
	counts->nMatchSkipYes[nextState] += count;

	skp = accum_count (-INFINITY,
			   matrix->fwdSkip[idx],
			   data->matchSkipYes[seqpos+1],
			   matrix->backSkip[inputIdx],
			   matrix,
			   &count,
			   NULL, NULL, NULL, NULL);  /* Skip -> Skip (input) */
	counts->nMatchSkipYes[nextState] += count;

	if (n_event < n_events) {
	  ioIdx = Seq_event_pair_index(seqpos+1,n_event+1);

	  mat = accum_count (mat,
			     matrix->fwdMatch[idx],
			     data->matchEventNo[seqpos]
			     + data->beginDeleteNo
			     + data->matchSkipNo[seqpos+1]
			     + data->matchEventLogLike[ioIdx],
			     matrix->backMatch[ioIdx],
			     matrix,
			     &count,
			     event,
			     &counts->matchMoment0[nextState],
			     &counts->matchMoment1[nextState],
			     &counts->matchMoment2[nextState]);  /* Match -> Match (input/output) */
	  counts->nMatchEventNo[state] += count;
	  counts->nBeginDeleteNo += count;
	  counts->nMatchSkipNo[nextState] += count;

	  skp = accum_count (skp,
			     matrix->fwdSkip[idx],
			     data->matchSkipNo[seqpos+1]
			     + data->matchEventLogLike[ioIdx],
			     matrix->backMatch[ioIdx],
			     matrix,
			     &count,
			     event,
			     &counts->matchMoment0[nextState],
			     &counts->matchMoment1[nextState],
			     &counts->matchMoment2[nextState]);  /* Skip -> Match (input/output) */
	  counts->nMatchSkipNo[nextState] += count;

	  del = accum_count (del,
			     matrix->fwdDelete[idx],
			     data->extendDeleteNo
			     + data->matchEventLogLike[ioIdx],
			     matrix->backMatch[ioIdx],
			     matrix,
			     &count,
			     NULL, NULL, NULL, NULL);            /* Delete -> Match (input/output) */
	  counts->nExtendDeleteNo += count;
	}
      }

      matrix->backMatch[idx] = mat;
      matrix->backSkip[idx] = skp;
      matrix->backDelete[idx] = del;
    }
  }

  for (n_event = n_events - 1; n_event >= 0; --n_event) {
    event = &data->events->event[n_event];

    st = accum_count (matrix->backStart[n_event],
		      matrix->fwdStart[n_event],
		      data->startEventYes
		      + data->nullEventLogLike[n_event + 1],
		      matrix->backStart[n_event + 1],
		      matrix,
		      &count,
		      NULL, NULL, NULL, NULL);  /* Start -> Start (output) */
    counts->nStartEventYes += count;

    matrix->backStart[n_event] = st;
  }
  matrix->backResult = matrix->backStart[0];

  if (LogThisAt(2))
    fprintf (stderr, "Backward log-likelihood is %Lg\n", matrix->backResult);

  if (LogThisAt(10))
    dump_seq_event_pair_matrix_to_file (model->config.debug_matrix_filename, "Backward", data, matrix->backStart, matrix->backMatch, matrix->backSkip, matrix->backDelete, matrix->backResult);
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
  counts->nNullEventYes += events->n_events;
  counts->nNullEventNo += 1.;
}

long double accum_count (long double back_src,
			 long double fwd_src,
			 long double trans,
			 long double back_dest,
			 Seq_event_pair_fb_matrix *matrix,
			 long double *count,
			 Fast5_event *event,
			 long double *moment0,
			 long double *moment1,
			 long double *moment2) {
  long double weight;
  weight = exp (fwd_src + trans + back_dest - matrix->fwdResult);
#if defined(NAN_DEBUG)  /* define NAN_DEBUG in util.h */
  if (isnan(weight)) {
    Abort ("NaN error in accum_count");
  }
#endif
  if (count)
    *count = weight;
  if (event) {
    *moment0 += weight * event->ticks;
    *moment1 += weight * event->sumticks_cur;
    *moment2 += weight * event->sumticks_cur_sq;
  }
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

  model->pNullEvent = BetaMode (counts->nNullEventYes + prior->nNullEventYes, counts->nNullEventNo + prior->nNullEventNo);
  model->pNullTick = BetaMode (counts->nullMoment0 + prior->nullMoment0, counts->nNullEventYes + prior->nNullEventYes);
  model->nullMean = (counts->nullMoment1 + prior->nullMoment1) / (counts->nullMoment0 + prior->nullMoment0);
  model->nullPrecision = 1. / MAX (DBL_MIN, ((counts->nullMoment2 + prior->nullMoment2) / (counts->nullMoment0 + prior->nullMoment0) - model->nullMean * model->nullMean));

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
    model->pMatchEvent[state] = BetaMode (counts->nMatchEventYes[state] + prior->nMatchEventYes[state], counts->nMatchEventNo[state] + prior->nMatchEventNo[state]);
    model->pMatchTick[state] = BetaMode (counts->matchMoment0[state] + prior->matchMoment0[state], counts->nMatchEventYes[state] + prior->nMatchEventYes[state]);
    model->pMatchSkip[state] = BetaMode (counts->nMatchSkipYes[state] + prior->nMatchSkipYes[state], counts->nMatchSkipNo[state] + prior->nMatchSkipNo[state]);
    model->matchMean[state] = (counts->matchMoment1[state] + prior->matchMoment1[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]);
    model->matchPrecision[state] = 1. / MAX (DBL_MIN, ((counts->matchMoment2[state] + prior->matchMoment2[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]) - model->matchMean[state] * model->matchMean[state]));
  }
  model->pBeginDelete = BetaMode (counts->nBeginDeleteYes + prior->nBeginDeleteYes, counts->nBeginDeleteNo + prior->nBeginDeleteNo);
  model->pExtendDelete = BetaMode (counts->nExtendDeleteYes + prior->nExtendDeleteYes, counts->nExtendDeleteNo + prior->nExtendDeleteNo);
  model->pStartEvent = BetaMode (counts->nStartEventYes + prior->nStartEventYes, counts->nStartEventNo + prior->nStartEventNo);

  if (dummy_prior != NULL)
    delete_seq_event_pair_counts (dummy_prior);
}

void inc_seq_event_pair_counts_via_fb (Seq_event_pair_model* model, Seq_event_pair_counts* counts, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_fb_matrix* matrix;
  long double fwdResult, nullModel, logPrior, logNullPrior;
  matrix = new_seq_event_pair_fb_matrix (model, seqlen, seq, events);
  fill_seq_event_pair_fb_matrix_and_inc_counts (matrix, counts);
  fwdResult = matrix->fwdResult;
  nullModel = matrix->data->nullModel;
  logPrior = matrix->data->logPrior;
  logNullPrior = matrix->data->logNullPrior;
  delete_seq_event_pair_fb_matrix (matrix);
  counts->loglike += fwdResult + logPrior - nullModel - logNullPrior;
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
      /* get model dataset */
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
	  int kmer_len = (int) H5Tget_size (strtype_id);
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

	  /* normalize */
	  normalize_model (model);

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
    model->nullPrecision = 1 / MAX (DBL_MIN, (m2 - m1 * m1));

    /* set all boolean probabilities to 0.5 */
    model->pBeginDelete = model->pExtendDelete = model->pStartEvent = model->pNullEvent = model->pNullTick = 0.5;
    for (state = 0; state < model->states; ++state) {
      model->pMatchEvent[state] = 0.5;
      model->pMatchTick[state] = 0.5;
    }
  }

  /* return */
  return ret;
}

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays) {
  int iter;
  long double loglike, prev_loglike;
  Seq_event_pair_counts *counts, *prior;

  Logger *logger = model->logger;

  prior = new_seq_event_pair_counts_minimal_prior (model);

  prev_loglike = 0.;
  for (iter = 0; iter < model->config.seq_evt_pair_EM_max_iterations; ++iter) {
    counts = get_seq_event_pair_counts (model, seqs, event_arrays);
    loglike = counts->loglike;

    optimize_seq_event_pair_model_for_counts (model, counts, prior);
    delete_seq_event_pair_counts (counts);

    if (LogThisAt(1))
      Warn ("Baum-Welch iteration %d: log-likelihood ratio score %Lg", iter + 1, loglike);

    if (iter > 0 && prev_loglike != 0. && fabsl(prev_loglike) != INFINITY
	&& ((loglike-prev_loglike)/fabsl(prev_loglike)) < model->config.seq_evt_pair_EM_min_fractional_loglike_increment)
      break;
    prev_loglike = loglike;
  }

  delete_seq_event_pair_counts (prior);
}

void fit_seq_event_null_model (Seq_event_pair_model* model, Vector* event_arrays) {
  Seq_event_pair_counts *counts = new_seq_event_pair_counts (model);
  Seq_event_pair_counts *prior = new_seq_event_pair_counts_minimal_prior (model);
  reset_seq_event_null_counts (counts);
  for (int n = 0; n < (int) VectorSize(event_arrays); ++n)
    inc_seq_event_null_counts_from_fast5 (counts, (Fast5_event_array*) VectorGet (event_arrays, n));
  optimize_seq_event_null_model_for_counts (model, counts, prior);
  delete_seq_event_pair_counts (counts);
  delete_seq_event_pair_counts (prior);
}

void fit_seq_event_kmer_model (Seq_event_pair_model* model, Kseq_container* seqs) {
  int *kmerCount = SafeCalloc (model->states, sizeof(int));
  int totalKmers = 0;
  for (int n = 0; n < seqs->n; ++n)
    for (int seqpos = 0; seqpos <= seqs->len[n] - model->order; ++seqpos) {
      ++kmerCount[decode_state_identifier (model->order, seqs->seq[n] + seqpos)];
      ++totalKmers;
    }
  /* add Laplace pseudocounts */
  for (int state = 0; state < model->states; ++state) {
    ++kmerCount[state];
    ++totalKmers;
  }
  /* convert to probabilities */
  for (int state = 0; state < model->states; ++state)
    model->kmerProb[state] = ((double) kmerCount[state]) / (double) totalKmers;
  model->emitProb = meanLengthToEventProb (((double) totalKmers) / (double) seqs->n);
  SafeFree (kmerCount);
}

Seq_event_pair_counts* get_seq_event_pair_counts (Seq_event_pair_model* model, Kseq_container* seqs, Vector* event_arrays) {
  int n_seq, *seqrev_len;
  char **rev, **seqrev;
  void **events_iter;
  Fast5_event_array *events;
  Seq_event_pair_counts *counts, **seq_counts, **seqrev_counts;

  Logger *logger = model->logger;

  if (!model->config.all_vs_all)
    Assert (seqs->n == (int) VectorSize(event_arrays), "Must have same number of reads as sequences; or, specify -allvsall");

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
  reset_seq_event_null_counts (counts);

  seqrev_counts = SafeMalloc (2 * seqs->n * sizeof(Seq_event_pair_counts*));
  seq_counts = SafeMalloc (seqs->n * sizeof(Seq_event_pair_counts*));
  for (n_seq = 0; n_seq < 2 * seqs->n; ++n_seq)
    seqrev_counts[n_seq] = new_seq_event_pair_counts (model);
  for (n_seq = 0; n_seq < seqs->n; ++n_seq)
    seq_counts[n_seq] = seqrev_counts[2*n_seq];

  for (events_iter = event_arrays->begin; events_iter != event_arrays->end; ++events_iter) {
    events = (Fast5_event_array*) *events_iter;
    int n_seq_begin = model->config.all_vs_all ? 0 : (2*(events_iter - event_arrays->begin));
    int n_seq_end = model->config.all_vs_all ? (2 * seqs->n) : (n_seq_begin + 2);
    int n_seq_step = model->config.both_strands ? 1 : 2;
    for (n_seq = n_seq_begin; n_seq < n_seq_end; n_seq += n_seq_step) {
      reset_seq_event_pair_counts (seqrev_counts[n_seq]);
      inc_seq_event_pair_counts_via_fb (model, seqrev_counts[n_seq], seqrev_len[n_seq], seqrev[n_seq], events);
      if (LogThisAt(2))
	Warn ("FAST5 file \"%s\", sequence \"%s\" (%s strand): log-likelihood = %Lg", events->name == NULL ? "<none>" : events->name, seqs->name[n_seq/2], (n_seq % 2) ? "reverse" : "forward", seqrev_counts[n_seq]->loglike);
    }
    add_weighted_seq_event_pair_counts (counts, model->config.both_strands ? seqrev_counts : seq_counts, (model->config.both_strands ? 2 : 1) * seqs->n);
  }

  for (n_seq = 0; n_seq < seqs->n; ++n_seq)
    SafeFree (rev[n_seq]);
  SafeFree (rev);

  for (n_seq = 0; n_seq < 2 * seqs->n; ++n_seq)
    delete_seq_event_pair_counts (seqrev_counts[n_seq]);
  SafeFree (seqrev_counts);
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
  align->basecall_path = new_labeled_seq_event_pair_path();
  return align;
}

void delete_seq_event_pair_alignment (Seq_event_pair_alignment* align) {
  delete_labeled_seq_event_pair_path (align->basecall_path);
  SafeFreeOrNull (align->events_at_pos);
  SafeFree (align);
}

void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, Seq_event_pair_model* model, int strand, FILE* out) {
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
    len = (int) strlen (evt->model_state);
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

void write_seq_event_pair_alignment_as_stockholm (Seq_event_pair_alignment* align, Seq_event_pair_model* model, int strand, FILE* out) {
  StringVector *seqrow, *evtrow;
  int n_event, n, seqpos_offset, name_width;
  char namefmt[100], *revcomp;

  Logger *logger = model->logger;

  revcomp = strand > 0 ? NULL : new_revcomp_seq (align->seq, (int) strlen(align->seq));

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

  name_width = 1 + MAX ((int) strlen (align->seqname), (int) strlen (align->events->name));
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

  if (LogThisAt(3)) {
    char state_id[20];
    fprintf (out, "#=GF event NumEvent\tEventMean\tEventSD\tModelState\tModelMean\tSeqPos\tSeqState\tStateMean\tStateSD\n");
    for (n_event = n = 0; n < align->start_n_event; ++n_event, ++n) {
      fprintf (out, "#=GR event %d\t%g\t%g\t%s\t%g - -\t%g\t%g\n", n_event, align->events->event[n_event].mean, align->events->event[n_event].stdv, align->events->event[n_event].model_state, align->events->event[n_event].model_level, model->nullMean, sqrt(1/model->nullPrecision));
    }
    for (seqpos_offset = 0; seqpos_offset <= align->end_seqpos - align->start_seqpos; ++seqpos_offset) {
      int state = decode_state_identifier (model->order, (char*) (align->seq + align->start_seqpos + seqpos_offset - model->order + 1));
      encode_state_identifier (state, model->order, state_id);
      for (n = 0; n < align->events_at_pos[seqpos_offset]; ++n_event, ++n) {
	fprintf (out, "#=GR event %d\t%g\t%g\t%s\t%g\t%d\t%s\t%g\t%g\n", n_event, align->events->event[n_event].mean, align->events->event[n_event].stdv, align->events->event[n_event].model_state, align->events->event[n_event].model_level, align->start_seqpos + seqpos_offset, state_id, model->matchMean[state], sqrt(1/model->matchPrecision[state]));
      }
    }
  }

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
  mx->vitSkip = SafeMalloc (matrix_cells * sizeof(long double));
  mx->vitDelete = SafeMalloc (matrix_cells * sizeof(long double));

  return mx;
}

void delete_seq_event_pair_viterbi_matrix (Seq_event_pair_viterbi_matrix* mx) {
  delete_seq_event_pair_data (mx->data);
  SafeFree (mx->vitStart);
  SafeFree (mx->vitMatch);
  SafeFree (mx->vitSkip);
  SafeFree (mx->vitDelete);
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
  long double mat, del, skp;
  unsigned long idx, inputIdx, outputIdx, ioIdx;
  Fast5_event* event;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;

  Logger *logger = model->logger;

  /* update data */
  precalc_seq_event_pair_data (data);

  /* fill Viterbi */
  matrix->vitStart[0] = 0.;
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    matrix->vitStart[n_event]
      = matrix->vitStart[n_event - 1]
      + data->startEventYes
      + data->nullEventLogLike[n_event];  /* Start -> Start (output) */
  }

  matrix->vitResult = -INFINITY;

  if (LogThisAt(2))
    init_progress ("Viterbi matrix");
  for (seqpos = order; seqpos <= seqlen; ++seqpos) {
    for (n_event = 0; n_event <= n_events; ++n_event) {

      if (LogThisAt(2))
	log_progress ((seqpos - order) / (double) (seqlen - order), "base %d/%d", seqpos, seqlen);

      idx = Seq_event_pair_index(seqpos,n_event);
      mat = -INFINITY;

      if (n_event > 0) {
	outputIdx = Seq_event_pair_index(seqpos,n_event-1);
	event = &data->events->event[n_event - 1];
	mat = max_func (mat,
			matrix->vitMatch[outputIdx]
			+ data->matchEventYes[seqpos]
			+ data->matchEventLogLike[idx]);     /* Match -> Match (output) */

	mat = max_func (mat,
			matrix->vitStart[n_event-1]
			+ data->startEventNo
			+ data->matchEventLogLike[idx]);     /* Start -> Match (output) */
      }

      if (seqpos == order) {
	skp = del = -INFINITY;
      } else {  /* seqpos > order */
	inputIdx = Seq_event_pair_index(seqpos-1,n_event);
	del = max_func
	  (matrix->vitMatch[inputIdx]
	   + data->matchEventNo[seqpos-1]
	   + data->beginDeleteYes,     /* Match -> Delete (input) */

	   matrix->vitDelete[inputIdx]
	   + data->extendDeleteYes);   /* Delete -> Delete (input) */

	skp = max_func
	  (matrix->vitMatch[inputIdx]
	   + data->matchEventNo[seqpos-1]
	   + data->beginDeleteNo
	   + data->matchSkipYes[seqpos],  /* Match -> Skip (input) */

	   matrix->vitSkip[inputIdx]
	   + data->matchSkipYes[seqpos]);  /* Skip -> Skip (input) */

	if (n_event > 0) {
	  ioIdx = Seq_event_pair_index(seqpos-1,n_event-1);
	  mat = max_func (mat,
			  matrix->vitMatch[ioIdx]
			  + data->matchEventNo[seqpos-1]
			  + data->beginDeleteNo
			  + data->matchSkipNo[seqpos]
			  + data->matchEventLogLike[idx]);     /* Match -> Match (input/output) */

	  mat = max_func (mat,
			  matrix->vitSkip[ioIdx]
			  + data->matchSkipNo[seqpos]
			  + data->matchEventLogLike[idx]);     /* Skip -> Match (input/output) */

	  mat = max_func (mat,
			  matrix->vitDelete[ioIdx]
			  + data->extendDeleteNo
			  + data->matchEventLogLike[idx]);     /* Delete -> Match (input/output) */
	}
      }

      matrix->vitMatch[idx] = mat;
      matrix->vitSkip[idx] = skp;
      matrix->vitDelete[idx] = del;
    }

    matrix->vitResult = max_func
      (matrix->vitResult,
       matrix->vitMatch[Seq_event_pair_index(seqpos,n_events)]
       + data->matchEventNo[seqpos]);  /* Match -> End (input) */
  }

  if (LogThisAt(2))
    fprintf (stderr, "Viterbi log-likelihood is %Lg (prior %Lg, null model %Lg, null prior %Lg)\n", matrix->vitResult, matrix->data->logPrior, matrix->data->nullModel, matrix->data->logNullPrior);

  if (LogThisAt(10))
    dump_seq_event_pair_matrix_to_file (model->config.debug_matrix_filename, "Viterbi", data, matrix->vitStart, matrix->vitMatch, matrix->vitSkip, matrix->vitDelete, matrix->vitResult);
}

Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, states, n_event, end_seqpos, start_seqpos, start_n_event, src_kmer, dest_kmer, n, k;
  long double loglike;
  unsigned long idx, inputIdx, outputIdx, ioIdx;
  Fast5_event* event;
  Seq_event_pair_transition trans;
  Seq_event_pair_state state;
  Vector *events_emitted;
  Seq_event_pair_alignment *align;

  data = matrix->data;
  seqlen = data->seqlen;
  n_events = data->events->n_events;
  model = data->model;
  order = model->order;
  states = model->states;

  Logger* logger = model->logger;
  
  start_seqpos = seqlen;
  start_n_event = n_events;
  
  /* traceback from End state */
  end_seqpos = -1;
  loglike = -INFINITY;
  for (seqpos = seqlen; seqpos >= order; --seqpos)
    update_max (&loglike,
		&end_seqpos,
		matrix->vitMatch[Seq_event_pair_index(seqpos,n_events)]
		+ data->matchEventNo[seqpos],  /* Match -> End */
		seqpos);
  Assert (end_seqpos >= 0, "Traceback failed");

  seqpos = end_seqpos;
  n_event = n_events;
  state = PairMatchState;
  inputIdx = outputIdx = ioIdx = -1;

  events_emitted = newVector (IntCopy, IntDelete, IntPrint);
  VectorPushBack (events_emitted, IntNew(0));

  align = new_seq_event_pair_alignment (data->events, data->seq, seqlen);
  int end_kmer = decode_state_identifier (order, data->seq + end_seqpos - order);
  ListInsertBefore (align->basecall_path, align->basecall_path->head, new_labeled_seq_event_pair_transition (model, MatchEndTransition, end_kmer, end_kmer, NULL));
  
  while (state != PairStartState) {
    idx = Seq_event_pair_index(seqpos,n_event);
    if (seqpos > 0)
      inputIdx = Seq_event_pair_index(seqpos-1,n_event);
    if (n_event > 0)
      outputIdx = Seq_event_pair_index(seqpos,n_event-1);
    if (seqpos > 0 && n_event > 0)
      ioIdx = Seq_event_pair_index(seqpos-1,n_event-1);

    event = n_event > 0 ? &data->events->event[n_event - 1] : NULL;

    loglike = -INFINITY;
    trans = UndefinedTransition;

    switch (state) {
    case PairMatchState:
      if (n_event > 0) {
	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitMatch[outputIdx]
		    + data->matchEventYes[seqpos]
		    + data->matchEventLogLike[idx],
		    MatchMatchOutTransition);    /* Match -> Match (output) */

	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitStart[n_event-1]
		    + data->startEventNo
		    + data->matchEventLogLike[idx],
		    StartMatchOutTransition);   /* Start -> Match (output) */      

	if (seqpos > order) {
	  update_max (&loglike,
		      (int*) &trans,
		      matrix->vitMatch[ioIdx]
		      + data->matchEventNo[seqpos-1]
		      + data->beginDeleteNo
		      + data->matchSkipNo[seqpos]
		      + data->matchEventLogLike[idx],
		      MatchMatchInOutTransition);    /* Match -> Match (input/output) */

	  update_max (&loglike,
		      (int*) &trans,
		      matrix->vitSkip[ioIdx]
		      + data->matchSkipNo[seqpos]
		      + data->matchEventLogLike[idx],
		      SkipMatchInOutTransition);    /* Skip -> Match (input/output) */

	  update_max (&loglike,
		      (int*) &trans,
		      matrix->vitDelete[ioIdx]
		      + data->extendDeleteNo
		      + data->matchEventLogLike[idx],
		      DeleteMatchInOutTransition);    /* Delete -> Match (input/output) */
	}
      } else
	Abort ("Traceback error");
      Assert (loglike == matrix->vitMatch[idx], "Traceback error");
      break;

    case PairSkipState:
      update_max (&loglike,
		  (int*) &trans,
		  matrix->vitMatch[inputIdx]
		  + data->matchEventNo[seqpos-1]
		  + data->beginDeleteNo
		  + data->matchSkipYes[seqpos],
		  MatchSkipInTransition);    /* Match -> Skip (input) */

      update_max (&loglike,
		  (int*) &trans,
		  matrix->vitSkip[inputIdx]
		  + data->matchSkipYes[seqpos],
		  SkipSkipInTransition);    /* Skip -> Skip (input) */

      Assert (loglike == matrix->vitSkip[idx], "Traceback error");
      break;

    case PairDeleteState:
	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitMatch[inputIdx]
		    + data->matchEventNo[seqpos-1]
		    + data->beginDeleteYes,
		    MatchDeleteInTransition);   /* Match -> Delete (input) */

	update_max (&loglike,
		    (int*) &trans,
		    matrix->vitDelete[inputIdx]
		    + data->extendDeleteYes,
		    DeleteDeleteInTransition);   /* Delete -> Delete (input) */

	Assert (loglike == matrix->vitDelete[idx], "Traceback error");
	break;

    default:
      Abort ("Unknown traceback state");
      break;
    }

    dest_kmer = decode_state_identifier (order, data->seq + seqpos - order);

    switch (trans) {

    case StartMatchOutTransition:
      state = PairStartState;
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      start_n_event = n_event;
      start_seqpos = seqpos;
      break;

    case MatchMatchOutTransition:
      Assert (state == PairMatchState, "oops");
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      break;

    case MatchSkipInTransition:
      Assert (state == PairSkipState, "oops");
      state = PairMatchState;
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case SkipSkipInTransition:
      Assert (state == PairSkipState, "oops");
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case MatchMatchInOutTransition:
      Assert (state == PairMatchState, "oops");
      --seqpos;
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case SkipMatchInOutTransition:
      Assert (state == PairMatchState, "oops");
      state = PairSkipState;
      --seqpos;
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case MatchDeleteInTransition:
      Assert (state == PairDeleteState, "oops");
      state = PairMatchState;
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case DeleteDeleteInTransition:
      Assert (state == PairDeleteState, "oops");
      --seqpos;
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case DeleteMatchInOutTransition:
      Assert (state == PairMatchState, "oops");
      state = PairDeleteState;
      --seqpos;
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      VectorPushBack (events_emitted, IntNew(0));
      break;

    case UndefinedTransition:
    default:
      Abort ("Unknown traceback transition");
      break;
    }

    src_kmer = decode_state_identifier (order, data->seq + seqpos - order);
    ListInsertBefore (align->basecall_path, align->basecall_path->head, new_labeled_seq_event_pair_transition (model, trans, src_kmer, dest_kmer, event));
  }

  align->log_likelihood_ratio = matrix->vitResult + matrix->data->logPrior - matrix->data->nullModel - matrix->data->logNullPrior;
  align->start_seqpos = start_seqpos - 1;
  align->end_seqpos = end_seqpos - 1;
  align->start_n_event = start_n_event;
  align->events_at_pos = SafeMalloc (VectorSize(events_emitted) * sizeof(int));
  for (n = ((int) VectorSize(events_emitted)) - 1, k = 0; n >= 0; --n, ++k)
    align->events_at_pos[k] = *((int*) VectorGet (events_emitted, n));

  for (n_event = start_n_event - 1; n_event >= 0; --n_event)
    ListInsertBefore (align->basecall_path, align->basecall_path->head, new_labeled_seq_event_pair_transition (model, StartStartOutTransition, src_kmer, src_kmer, &data->events->event[n_event]));
  
  deleteVector (events_emitted);

  if (LogThisAt(5)) {
    fprintf (stderr, "Traceback:\n");
    print_labeled_seq_event_pair_path (stderr, align->basecall_path);
  }
  
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

  for (strand = +1; strand >= (model->config.both_strands ? -1 : +1); strand -= 2) {
    matrix = new_seq_event_pair_viterbi_matrix (model, seqlen, strand > 0 ? seq : rev, events);
    fill_seq_event_pair_viterbi_matrix (matrix);

    align = get_seq_event_pair_viterbi_matrix_traceback (matrix);
    align->seqname = seqname;

    if (align->log_likelihood_ratio >= log_odds_ratio_threshold)
      (*write_func) (align, model, strand, out);

    delete_seq_event_pair_viterbi_matrix (matrix);
    delete_seq_event_pair_alignment (align);
  }

  SafeFree (rev);
}

void dump_seq_event_pair_matrix_to_file (const char* filename, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double* mxSkip, long double *mxDelete, long double mxResult) {
  FILE *file;
  file = fopen (filename, "a");
  Assert (file != NULL, "Could not write DP matrix to file %s", filename);
  dump_seq_event_pair_matrix (file, algorithm, data, mxStart, mxMatch, mxSkip, mxDelete, mxResult);
  fclose (file);
}

void dump_seq_event_pair_matrix (FILE* file, const char* algorithm, Seq_event_pair_data *data, long double *mxStart, long double *mxMatch, long double* mxSkip, long double *mxDelete, long double mxResult) {
  int n_event, seqpos, n_events, seqlen, order;
  Fast5_event *event, *next_event;
  char *id;
  unsigned long idx, out_idx, inout_idx;
  time_t rawtime;
  struct tm *rawtime_tm;

  order = data->model->order;
  n_events = data->events->n_events;
  seqlen = data->seqlen;

  id = SafeMalloc ((data->model->order + 1) * sizeof(char));

  time (&rawtime);
  rawtime_tm = localtime (&rawtime);
  fprintf (file, "%s matrix: %s\n", algorithm, asctime(rawtime_tm));
  for (n_event = 0; n_event <= n_events; ++n_event) {
    event = n_event > 0 ? &data->events->event[n_event - 1] : NULL;
    next_event = n_event < n_events ? &data->events->event[n_event] : NULL;
    fprintf (file, "(%d,*) Event %d", n_event, n_event);
    if (n_event > 0)
      fprintf (file, " (n=%g,sum=%g,sumsq=%g)", event->ticks, event->sumticks_cur, event->sumticks_cur_sq);
    fprintf (file, ", pre-seq\n  start %Lg (in>s:0", mxStart[n_event]);
    if (n_event < n_events)
      fprintf (file, ", out>s:%Lg", data->startEventYes + data->nullEventLogLike[n_event + 1]);
    fprintf (file, ")\n");
    for (seqpos = order; seqpos <= seqlen; ++seqpos) {
      encode_state_identifier (data->state[seqpos], data->model->order, id);
      idx = Seq_event_pair_index(seqpos,n_event);
      out_idx = n_event < n_events ? Seq_event_pair_index(seqpos,n_event+1) : -1;
      inout_idx = (n_event < n_events && seqpos < seqlen) ? Seq_event_pair_index(seqpos+1,n_event+1) : -1;
      fprintf (file, "(%d,%d) Event %d", n_event, seqpos, n_event);
      if (n_event > 0)
	fprintf (file, " (n=%g,sum=%g,sumsq=%g)", event->ticks, event->sumticks_cur, event->sumticks_cur_sq);
      fprintf (file, ", seqpos %d (base=%c,state=%s)\n  start %Lg", seqpos, data->seq[seqpos-1], id, mxStart[n_event]);
      if (seqpos < seqlen)
	fprintf (file, " (in>s:0%s", n_event < n_events ? ", " : ")");
      if (n_event < n_events)
	fprintf (file, "%sout>s:%Lg, out>m:%Lg)", seqpos < seqlen ? "" : " (", data->startEventYes + data->nullEventLogLike[n_event + 1], data->startEventNo + data->matchEventLogLike[out_idx]);
      fprintf (file, "\n  match %Lg (", mxMatch[idx]);
      if (seqpos < seqlen) {
	fprintf (file, "in>d:%Lg, in>sk:%Lg, ", data->matchEventNo[seqpos] + data->beginDeleteYes, data->matchEventNo[seqpos] + data->beginDeleteNo + data->matchSkipYes[seqpos+1]);
	if (n_event < n_events)
	  fprintf (file, "io>m:%Lg, ", data->matchEventNo[seqpos] + data->beginDeleteNo + data->matchSkipNo[seqpos+1] + data->matchEventLogLike[inout_idx]);
      }
      if (n_event < n_events)
	fprintf (file, "out>m:%Lg", data->matchEventYes[seqpos] + data->matchEventLogLike[out_idx]);
      else
	fprintf (file, "null>e:%Lg", data->matchEventNo[seqpos]);
      fprintf (file, ")\n   skip %Lg", mxDelete[idx]);
      if (seqpos < seqlen) {
	fprintf (file, " (in>sk:%Lg", data->matchSkipYes[seqpos+1]);
	if (n_event < n_events)
	  fprintf (file, ", io>m:%Lg", data->matchSkipNo[seqpos+1] + data->matchEventLogLike[inout_idx]);
	fprintf (file, ")");
      }
      fprintf (file, "\n delete %Lg", mxDelete[idx]);
      if (seqpos < seqlen) {
	fprintf (file, " (in>d:%Lg", data->extendDeleteYes);
	if (n_event < n_events)
	  fprintf (file, ", io>m:%Lg", data->extendDeleteNo + data->matchEventLogLike[inout_idx]);
	fprintf (file, ")");
      }
      fprintf (file, "\n");
    }
  }
  fprintf (file, "Result %Lg\n\n", mxResult);

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
    xmlTextWriterBooleanCount (writer, XMLPREFIX(SKIP), counts->nMatchSkipYes[state], counts->nMatchSkipNo[state]);
    xmlTextWriterBooleanCount (writer, XMLPREFIX(EVENTS), counts->nMatchEventYes[state], counts->nMatchEventNo[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M0), "%Lg", counts->matchMoment0[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M1), "%Lg", counts->matchMoment1[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M2), "%Lg", counts->matchMoment2[state]);
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(START));
  xmlTextWriterBooleanCount (writer, XMLPREFIX(EVENTS), counts->nStartEventYes, counts->nStartEventNo);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(NULLMODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M0), "%Lg", counts->nullMoment0);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M1), "%Lg", counts->nullMoment1);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(M2), "%Lg", counts->nullMoment2);
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);

  SafeFree (id);

  return deleteXmlTextWriterLeavingText (writer);
}

xmlChar* make_squiggle_svg (Fast5_event_array *events, Seq_event_pair_model* model) {
  xmlTextWriterPtr writer;
  writer = newXmlTextWriter();

  if (events->n_events > 0) {
    // find viewport: width (total ticks), height (level range)
    double ticks = 0, min_lev = 0, max_lev = 0;
    int n_event;
    for (n_event = 0; n_event < events->n_events; ++n_event) {
      Fast5_event* event = events->event + n_event;
      double emin = event->mean - event->stdv, emax = event->mean + event->stdv;
      ticks += event->ticks;
      if (n_event == 0 || emin < min_lev)
	min_lev = emin;
      if (n_event == 0 || emax > max_lev)
	max_lev = emax;
    }
    long
      w = (long) (ticks * pixels_per_tick),
      h = (long) ((max_lev - min_lev) * pixels_per_lev);

    // <svg>
    xmlNode widthAttr, heightAttr, viewBoxAttr, versionAttr;
    char widthText[20], heightText[20], viewBoxText[100];

    widthAttr.name = (xmlChar*) "width";
    sprintf (widthText, "%ld", w);
    widthAttr.content = (xmlChar*) widthText;
    widthAttr.next = &heightAttr;

    heightAttr.name = (xmlChar*) "height";
    sprintf (heightText, "%ld", h);
    heightAttr.content = (xmlChar*) heightText;
    heightAttr.next = &viewBoxAttr;

    viewBoxAttr.name = (xmlChar*) "viewBox";
    sprintf (viewBoxText, "0 0 %ld %ld", w, h);
    viewBoxAttr.content = (xmlChar*) viewBoxText;
    viewBoxAttr.next = &versionAttr;

    versionAttr.name = (xmlChar*) "version";
    versionAttr.content = (xmlChar*) "1.1";
    versionAttr.next = NULL;

    xmlTextWriterStartElementWithAttrs (writer, (xmlChar*) "svg", &widthAttr);

    // loop over events, drawing model level in background
    int pass;
    for (pass = 0; pass < 3; ++pass) {
      ticks = 0;
      for (n_event = 0; n_event < events->n_events; ++n_event) {
	Fast5_event* event = events->event + n_event;

	long xmin = (long) (ticks * pixels_per_tick);
	ticks += event->ticks;
	long xmax = (long) (ticks * pixels_per_tick);

	int state = decode_state_identifier (model->order, event->model_state);
	double model_mean = model->matchMean[state];
	double model_stdv = 1 / sqrt (model->matchPrecision[state]);

	double mmin = model_mean - model_stdv, mmax = model_mean + model_stdv;
	long mymin = (long) ((mmin - min_lev) * pixels_per_lev);
	long mymax = (long) ((mmax - min_lev) * pixels_per_lev);

	double emin = event->mean - event->stdv, emax = event->mean + event->stdv;
	long eymin = (long) ((emin - min_lev) * pixels_per_lev);
	long eymax = (long) ((emax - min_lev) * pixels_per_lev);

	switch (pass) {
	case 0:
	  {
	    // rect for model level
	    xmlNode levXAttr, levYAttr, levWidthAttr, levHeightAttr, levStrokeWidthAttr, levFillAttr;
	    char levXText[20], levYText[20], levWidthText[20], levHeightText[20];

	    levXAttr.name = (xmlChar*) "x";
	    sprintf (levXText, "%ld", xmin);
	    levXAttr.content = (xmlChar*) levXText;
	    levXAttr.next = &levYAttr;

	    levYAttr.name = (xmlChar*) "y";
	    sprintf (levYText, "%ld", h - mymax);
	    levYAttr.content = (xmlChar*) levYText;
	    levYAttr.next = &levWidthAttr;

	    levWidthAttr.name = (xmlChar*) "width";
	    sprintf (levWidthText, "%ld", xmax - xmin + 1);
	    levWidthAttr.content = (xmlChar*) levWidthText;
	    levWidthAttr.next = &levHeightAttr;

	    levHeightAttr.name = (xmlChar*) "height";
	    sprintf (levHeightText, "%ld", mymax - mymin);
	    levHeightAttr.content = (xmlChar*) levHeightText;
	    levHeightAttr.next = &levStrokeWidthAttr;

	    levStrokeWidthAttr.name = (xmlChar*) "stroke-width";
	    levStrokeWidthAttr.content = (xmlChar*) "0";
	    levStrokeWidthAttr.next = &levFillAttr;

	    levFillAttr.name = (xmlChar*) "fill";
	    levFillAttr.content = (xmlChar*) "blue";
	    levFillAttr.next = NULL;

	    xmlTextWriterStartElementWithAttrs (writer, (xmlChar*) "rect", &levXAttr);
	    xmlTextWriterEndElement (writer);
	  }
	  break;

	case 1:
	  {
	    // rect for event
	    xmlNode evtXAttr, evtYAttr, evtWidthAttr, evtHeightAttr, evtStrokeAttr, evtStrokeWidthAttr, evtFillAttr;
	    char evtXText[20], evtYText[20], evtWidthText[20], evtHeightText[20];

	    evtXAttr.name = (xmlChar*) "x";
	    sprintf (evtXText, "%ld", xmin);
	    evtXAttr.content = (xmlChar*) evtXText;
	    evtXAttr.next = &evtYAttr;

	    evtYAttr.name = (xmlChar*) "y";
	    sprintf (evtYText, "%ld", h - eymax);
	    evtYAttr.content = (xmlChar*) evtYText;
	    evtYAttr.next = &evtWidthAttr;

	    evtWidthAttr.name = (xmlChar*) "width";
	    sprintf (evtWidthText, "%ld", xmax - xmin);
	    evtWidthAttr.content = (xmlChar*) evtWidthText;
	    evtWidthAttr.next = &evtHeightAttr;

	    evtHeightAttr.name = (xmlChar*) "height";
	    sprintf (evtHeightText, "%ld", eymax - eymin);
	    evtHeightAttr.content = (xmlChar*) evtHeightText;
	    evtHeightAttr.next = &evtStrokeAttr;

	    evtStrokeAttr.name = (xmlChar*) "stroke";
	    evtStrokeAttr.content = (xmlChar*) "black";
	    evtStrokeAttr.next = &evtStrokeWidthAttr;

	    evtStrokeWidthAttr.name = (xmlChar*) "stroke-width";
	    evtStrokeWidthAttr.content = (xmlChar*) "1";
	    evtStrokeWidthAttr.next = &evtFillAttr;

	    evtFillAttr.name = (xmlChar*) "fill";
	    evtFillAttr.content = (xmlChar*) "none";
	    evtFillAttr.next = NULL;

	    xmlTextWriterStartElementWithAttrs (writer, (xmlChar*) "rect", &evtXAttr);
	    xmlTextWriterEndElement (writer);
	  }
	  break;

	case 2:
	  {
	    // text for move
	    if (event->move > 0 || n_event == 0) {
	      xmlNode textXAttr, textYAttr, textFillAttr;
	      char textXText[20], textYText[20];

	      textXAttr.name = (xmlChar*) "x";
	      sprintf (textXText, "%ld", xmin);
	      textXAttr.content = (xmlChar*) textXText;
	      textXAttr.next = &textYAttr;

	      textYAttr.name = (xmlChar*) "y";
	      sprintf (textYText, "%ld", h - mymin);
	      textYAttr.content = (xmlChar*) textYText;
	      textYAttr.next = &textFillAttr;

	      textFillAttr.name = (xmlChar*) "fill";
	      textFillAttr.content = (xmlChar*) "lightblue";
	      textFillAttr.next = NULL;

	      xmlTextWriterStartElementWithAttrs (writer, (xmlChar*) "text", &textXAttr);
	      xmlTextWriterWriteFormatCDATA (writer, "%s", event->model_state + (n_event == 0 ? 0 : (strlen(event->model_state) - event->move)));
	      xmlTextWriterEndElement (writer);
	    }
	  }
	  break;

	  default:
	    break;
	}
      }
    }
  }

  xmlTextWriterEndElement (writer);
  return deleteXmlTextWriterLeavingText (writer);
}

void normalize_model (Seq_event_pair_model* model) {
  double sum = 0, sumsq = 0;
  for (int n = 0; n < model->states; ++n) {
    sum += model->matchMean[n];
    sumsq += model->matchMean[n] * model->matchMean[n];
  }
  double mean = sum / (double) model->states;
  double var = sumsq / (double) model->states - mean*mean;
  double sd = sqrt(var);
  for (int n = 0; n < model->states; ++n) {
    model->matchMean[n] = (model->matchMean[n] - mean) / sd;
    model->matchPrecision[n] *= var;
  }
  model->nullMean = (model->nullMean - mean) / sd;
  model->nullPrecision *= var;
}

Labeled_seq_event_pair_transition* new_labeled_seq_event_pair_transition (Seq_event_pair_model* model, Seq_event_pair_transition trans, int src_kmer, int dest_kmer, Fast5_event* emission) {
  Labeled_seq_event_pair_transition* lab = SafeMalloc (sizeof (Labeled_seq_event_pair_transition));
  lab->model = model;
  lab->trans = trans;
  lab->src_kmer = src_kmer;
  lab->dest_kmer = dest_kmer;
  lab->emission = emission;
  return lab;
}

void* copy_labeled_seq_event_pair_transition (void* vlab) {
  Labeled_seq_event_pair_transition* lab = (Labeled_seq_event_pair_transition*) vlab;
  return new_labeled_seq_event_pair_transition (lab->model, lab->trans, lab->src_kmer, lab->dest_kmer, lab->emission);
}

void delete_labeled_seq_event_pair_transition (void* vlab) {
  Labeled_seq_event_pair_transition* lab = (Labeled_seq_event_pair_transition*) vlab;
  SafeFree (lab);
}

static char* transition_name[] = { "Undefined", "Start->Start", "Start->Match", "Match->Match(out)", "Match->Skip", "Skip->Skip", "Match->Match(in/out)", "Skip->Match", "Match->Delete", "Delete->Delete", "Delete->Match", "Match->End" };

int seq_event_pair_transition_absorb_count (Seq_event_pair_transition trans) {
  return
    (trans == MatchSkipInTransition
     || trans == SkipSkipInTransition
     || trans == MatchMatchInOutTransition
     || trans == SkipMatchInOutTransition
     || trans == MatchDeleteInTransition
     || trans == DeleteDeleteInTransition
     || trans == DeleteMatchInOutTransition)
    ? 1
    : 0;
}

int seq_event_pair_transition_emit_count (Seq_event_pair_transition trans) {
  return
    (trans == StartStartOutTransition
     || trans == StartMatchOutTransition
     || trans == MatchMatchOutTransition
     || trans == MatchMatchInOutTransition
     || trans == SkipMatchInOutTransition
     || trans == DeleteMatchInOutTransition)
    ? 1
    : 0;
}

void print_labeled_seq_event_pair_transition (FILE *file, void* vlab) {
  Labeled_seq_event_pair_transition* lab = (Labeled_seq_event_pair_transition*) vlab;

  char* src_state_id = SafeMalloc ((lab->model->order + 1) * sizeof(char));
  char* dest_state_id = SafeMalloc ((lab->model->order + 1) * sizeof(char));
  encode_state_identifier (lab->src_kmer, lab->model->order, src_state_id);
  encode_state_identifier (lab->dest_kmer, lab->model->order, dest_state_id);

  int emit = seq_event_pair_transition_emit_count (lab->trans);
  int absorb = seq_event_pair_transition_absorb_count (lab->trans);
  
  fprintf (file, "%s %s %s", src_state_id, transition_name[lab->trans], dest_state_id);
  if (absorb)
    fprintf (file, " in:%c", token2base (lab->dest_kmer % AlphabetSize));
  if (emit)
    fprintf (file, " out:(%g %g %g)", lab->emission->ticks, lab->emission->mean, lab->emission->stdv);
  fprintf (file, " Score: transition %Lg", seq_event_pair_transition_trans_loglike (lab->model, lab->trans, lab->src_kmer, lab->dest_kmer));
  if (absorb)
    fprintf (file, ", in %Lg", seq_event_pair_transition_generator_loglike (lab->model, lab->trans, lab->dest_kmer));
  if (emit)
    fprintf (file, ", out %Lg", seq_event_pair_transition_emit_loglike (lab->model, lab->trans, lab->dest_kmer, lab->emission));
  fprintf (file, "\n");

  SafeFree(src_state_id);
  SafeFree(dest_state_id);
}

long double seq_event_pair_transition_trans_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int src_kmer_state, int dest_kmer_state) {
  long double t = -INFINITY;
  switch (trans) {

  case StartMatchOutTransition:
    t = log(1. - model->pStartEvent);
    break;

  case MatchMatchOutTransition:
    t = log(model->pMatchEvent[dest_kmer_state]);
    break;

  case MatchSkipInTransition:
    t = log(1. - model->pMatchEvent[src_kmer_state])
      + log(1. - model->pBeginDelete)
      + log(model->pMatchSkip[dest_kmer_state]);
    break;

  case SkipSkipInTransition:
    t = log(model->pMatchSkip[dest_kmer_state]);
    break;

  case MatchMatchInOutTransition:
    t = log(1. - model->pMatchEvent[src_kmer_state])
      + log(1. - model->pBeginDelete)
      + log(1. - model->pMatchSkip[dest_kmer_state]);
    break;

  case SkipMatchInOutTransition:
    t = log(1. - model->pMatchSkip[dest_kmer_state]);
    break;

  case MatchDeleteInTransition:
    t = log(1. - model->pMatchEvent[src_kmer_state])
      + log(model->pBeginDelete);
    break;

  case DeleteDeleteInTransition:
    t = log(model->pExtendDelete);
    break;

  case DeleteMatchInOutTransition:
    t = log(1. - model->pExtendDelete);
    break;

  case MatchEndTransition:
    t = log(1. - model->pMatchEvent[src_kmer_state]);
    break;
	
  case UndefinedTransition:
  default:
    Abort ("Unknown transition");
    break;
  }
  return t;
}

long double seq_event_pair_transition_emit_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int kmer_state, Fast5_event* emission) {
  if (seq_event_pair_transition_emit_count (trans)) {
    double mean, precision, pTick;
    if (trans == StartStartOutTransition) {
      mean = model->nullMean;
      precision = model->nullPrecision;
      pTick = model->pNullTick;
    } else {
      mean = model->matchMean[kmer_state];
      precision = model->matchPrecision[kmer_state];
    pTick = model->pMatchTick[kmer_state];
    }
    return log_event_density (emission,
			      mean,
			      precision,
			      log (precision),
			      log (pTick),
			      log (1. - pTick));
  }
  return 0;
}

long double seq_event_pair_transition_generator_loglike (Seq_event_pair_model* model, Seq_event_pair_transition trans, int dest_kmer_state) {
  if (seq_event_pair_transition_absorb_count (trans)) {
    int prefix = dest_kmer_state - (dest_kmer_state % AlphabetSize);
    double norm = 0;
    for (int offset = 0; offset < AlphabetSize; ++offset)
      norm += model->kmerProb[prefix + offset];
    return log (model->kmerProb[dest_kmer_state] / norm) + log (model->emitProb);
  }
  if (trans == MatchEndTransition)
    return log (1. - model->emitProb);
  return 0;
}


Labeled_seq_event_pair_path* new_labeled_seq_event_pair_path() {
  return newList (copy_labeled_seq_event_pair_transition,
		  delete_labeled_seq_event_pair_transition,
		  print_labeled_seq_event_pair_transition);
}

void delete_labeled_seq_event_pair_path (Labeled_seq_event_pair_path* path) {
  deleteList (path);
}

long double labeled_seq_event_pair_transition_conditional_loglike (Labeled_seq_event_pair_transition* lab) {
  return seq_event_pair_transition_trans_loglike (lab->model, lab->trans, lab->src_kmer, lab->dest_kmer) + seq_event_pair_transition_emit_loglike (lab->model, lab->trans, lab->dest_kmer, lab->emission);
}

long double labeled_seq_event_pair_transition_joint_loglike (Labeled_seq_event_pair_transition* lab) {
  return labeled_seq_event_pair_transition_conditional_loglike (lab) + seq_event_pair_transition_generator_loglike (lab->model, lab->trans, lab->dest_kmer);
}

long double labeled_seq_event_pair_path_conditional_loglike (Labeled_seq_event_pair_path* path) {
  double ll = 0;
  for (ListNode* node = path->head; node; node = node->next)
    ll += labeled_seq_event_pair_transition_conditional_loglike ((Labeled_seq_event_pair_transition*) node->value);
  return ll;
}

long double labeled_seq_event_pair_path_joint_loglike (Labeled_seq_event_pair_path* path) {
  double ll = 0;
  for (ListNode* node = path->head; node; node = node->next)
    ll += labeled_seq_event_pair_transition_joint_loglike ((Labeled_seq_event_pair_transition*) node->value);
  return ll;
}

void print_labeled_seq_event_pair_path (FILE* file, Labeled_seq_event_pair_path* path) {
  ListPrint (file, path);
  fprintf (file, "Log-likelihood: %Lg (conditional), %Lg (joint)\n", labeled_seq_event_pair_path_conditional_loglike (path),labeled_seq_event_pair_path_joint_loglike (path));
}
