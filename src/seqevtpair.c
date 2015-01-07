#include <math.h>
#include <string.h>

#include "seqevtpair.h"
#include "xmlutil.h"
#include "xmlkeywords.h"
#include "kseqcontainer.h"

int base2token (char base);
char token2base (int token);

void encode_state_identifier (int state, int order, char* state_id);
int decode_state_identifier (int order, char* state_id);

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
  return tokenize (base, dna_alphabet);
}

char token2base (int tok) {
  return tok < 0 || tok >= 4 ? 'N' : dna_alphabet[tok];
}

void encode_state_identifier (int state, int order, char* state_id) {
  int k;
  for (k = 0; k < order; ++k, state /= 4)
    state_id[order - k - 1] = token2base (state % 4);
  state_id[order] = '\0';
}

int decode_state_identifier (int order, char* state_id) {
  int k, token, p;
  for (token = 0, p = 1, k = 0; k < order; ++k, p *= 4)
    token += p * base2token (state_id[order - k - 1]);
  return token;
}

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml) {
  xmlNode *modelNode, *statesNode, *stateNode, *deleteNode, *startNode;
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
      state = decode_state_identifier (model->order, (char*) CHILDSTRING(stateNode,ID));
      model->pMatchEmit[state] = CHILDFLOAT(stateNode,EMIT);
      model->matchMean[state] = CHILDFLOAT(stateNode,MEAN);
      model->matchPrecision[state] = CHILDFLOAT(stateNode,PRECISION);
    }

  startNode = CHILD(modelNode,START);
  model->pStartEmit = CHILDFLOAT(startNode,EMIT);
  model->startMean = CHILDFLOAT(startNode,MEAN);
  model->startPrecision = CHILDFLOAT(startNode,PRECISION);

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
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(ID), "%s", id);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EMIT), "%g", model->pMatchEmit[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->matchMean[state]);
    xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(PRECISION), "%g", model->matchPrecision[state]);
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(START));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EMIT), "%g", model->pStartEmit);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->startMean);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(PRECISION), "%g", model->startPrecision);
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);
  return deleteXmlTextWriterLeavingText (writer);
}

/* Forward-backward matrix */

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_fb_matrix* mx;
  int n_events;

  mx = SafeMalloc (sizeof (Seq_event_pair_fb_matrix));
  mx->model = model;
  mx->seqlen = seqlen;
  mx->seq = seq;
  mx->events = events;

  n_events = events->n_events;

  mx->fwdStart = SafeMalloc ((n_events + 1) * sizeof(long double));
  mx->backStart = SafeMalloc ((n_events + 1) * sizeof(long double));

  mx->fwdMatch = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));
  mx->fwdDelete = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));
  mx->backMatch = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));
  mx->backDelete = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));

  mx->startEmitDensity = SafeMalloc ((n_events + 1) * sizeof(long double));
  mx->matchEmitDensity = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));
  mx->matchEmitProb = SafeMalloc ((n_events + 1) * (seqlen + 1) * sizeof(long double));

  mx->state = SafeMalloc ((seqlen + 1) * sizeof(int));

  return mx;
}

void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* mx) {
  SafeFree (mx->state);
  SafeFree (mx->startEmitDensity);
  SafeFree (mx->matchEmitDensity);
  SafeFree (mx->matchEmitProb);
  SafeFree (mx->fwdStart);
  SafeFree (mx->backStart);
  SafeFree (mx->fwdMatch);
  SafeFree (mx->fwdDelete);
  SafeFree (mx->backMatch);
  SafeFree (mx->backDelete);
  SafeFree (mx);
}

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model) {
  Seq_event_pair_counts* counts;
  counts = SafeMalloc (sizeof (Seq_event_pair_counts));
  counts->order = model->order;
  counts->states = model->states;
  counts->nEmit_yes = SafeMalloc (model->states * sizeof(long double));
  counts->nEmit_no = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment0 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment1 = SafeMalloc (model->states * sizeof(long double));
  counts->matchMoment2 = SafeMalloc (model->states * sizeof(long double));
  return counts;
}

void delete_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  SafeFree (counts->nEmit_yes);
  SafeFree (counts->nEmit_no);
  SafeFree (counts->matchMoment0);
  SafeFree (counts->matchMoment1);
  SafeFree (counts->matchMoment2);
  SafeFree (counts);
}

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts) {
  int state;
  for (state = 0; state < counts->states; ++state) {
    counts->nEmit_yes[state] = 0;
    counts->nEmit_no[state] = 0;
    counts->matchMoment0[state] = 0;
    counts->matchMoment1[state] = 0;
    counts->matchMoment2[state] = 0;
  }
  counts->nBeginDelete_yes = 0;
  counts->nBeginDelete_no = 0;
  counts->nExtendDelete_yes = 0;
  counts->nExtendDelete_no = 0;
  counts->startMoment0 = 0;
  counts->startMoment1 = 0;
  counts->startMoment2 = 0;
}

void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events) {
  /* MORE TO GO HERE ... pending inspection of basecalled FAST5 data */
}

static const double log_sqrt2pi = 1.83787706640935;
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
  int seqlen, n_events, order, seqpos, n_event, state;
  double logStartPrecision, *logMatchPrecision;
  Fast5_event* event;
  seqlen = matrix->seqlen;
  n_events = matrix->events->n_events;
  model = matrix->model;
  order = model->order;

  /* precalculate logs of precisions */
  logStartPrecision = log (model->startPrecision);
  logMatchPrecision = SafeMalloc (model->states * sizeof(double));
  for (state = 0; state < model->states; ++state)
    logMatchPrecision[state] = log (model->matchPrecision[state]);

  /* calculate states */
  for (seqpos = 0; seqpos < order; ++seqpos)
    matrix->state[seqpos] = -1;
  for (seqpos = order; seqpos <= seqlen; ++seqpos)
    matrix->state[seqpos] = decode_state_identifier (order, matrix->seq + seqpos - order);

  /* calculate emit densities */
  for (n_event = 0; n_event < n_events; ++n_event) {
    event = matrix->events->event[n_event];
    matrix->startEmitDensity[n_event] = log_event_density (event,
							   model->startMean,
							   model->startPrecision,
							   logStartPrecision);

    for (seqpos = 0; seqpos < order; ++seqpos)
      matrix->matchEmitDensity[Seq_event_pair_index(seqlen,seqpos,n_event)] = -INFINITY;

    for (seqpos = order; seqpos <= seqlen; ++seqpos) {
      state = matrix->state[seqpos];
      matrix->matchEmitDensity[Seq_event_pair_index(seqlen,seqpos,n_event)]
	= log_event_density (event,
			     model->matchMean[state],
			     model->matchPrecision[state],
			     logMatchPrecision[state]);
    }
  }

  /* fill forward */
  /* fill backward */
  /* update counts */

  /* free */
  SafeFree (logMatchPrecision);
}

void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* matrix, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Fast5_event_array* events, double minimum_fractional_log_likelihood_increase);
