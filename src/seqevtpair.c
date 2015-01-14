#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "seqevtpair.h"
#include "xmlutil.h"
#include "xmlkeywords.h"
#include "kseqcontainer.h"

static const double log_sqrt2pi = 1.83787706640935;

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

void update_max (long double *current_max, int* current_max_idx, long double candidate_max, int candidate_max_idx);

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
      state = decode_state_identifier (model->order, (char*) CHILDSTRING(stateNode,ID));
      model->pMatchEmit[state] = CHILDFLOAT(stateNode,EMIT);
      model->matchMean[state] = CHILDFLOAT(stateNode,MEAN);
      model->matchPrecision[state] = CHILDFLOAT(stateNode,PRECISION);
    }

  startNode = CHILD(modelNode,START);
  model->pStartEmit = CHILDFLOAT(startNode,EMIT);

  nullNode = CHILD(modelNode,NULLMODEL);
  model->pNullEmit = CHILDFLOAT(nullNode,EMIT);
  model->nullMean = CHILDFLOAT(nullNode,MEAN);
  model->nullPrecision = CHILDFLOAT(nullNode,PRECISION);

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
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, (xmlChar*) XMLPREFIX(NULLMODEL));
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(EMIT), "%g", model->pNullEmit);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(MEAN), "%g", model->nullMean);
  xmlTextWriterWriteFormatElement (writer, (xmlChar*) XMLPREFIX(PRECISION), "%g", model->nullPrecision);
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);
  return deleteXmlTextWriterLeavingText (writer);
}

/* Pairwise alignment data */

Seq_event_pair_data* new_seq_event_pair_data (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_data* data;
  int n_events, matrix_cells, seqpos;

  data = SafeMalloc (sizeof (Seq_event_pair_data));
  data->model = model;
  data->seqlen = seqlen;
  data->seq = seq;
  data->events = events;

  n_events = events->n_events;

  matrix_cells = (n_events + 1) * (seqlen - model->order + 1);
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
    loglike = log_event_density (&data->events->event[n_event],
				 model->nullMean,
				 model->nullPrecision,
				 logNullPrecision);
    data->nullEmitDensity[n_event + 1] = loglike;
    data->nullModel += loglike + data->nullEmitYes;
  }

  for (seqpos = 0; seqpos < order; ++seqpos) {
    data->matchEmitYes[seqpos] = -INFINITY;
    data->matchEmitNo[seqpos] = -INFINITY;
  }

  for (seqpos = order; seqpos <= seqlen; ++seqpos) {
    state = data->state[seqpos];
    data->matchEmitYes[seqpos] = log (model->pMatchEmit[state]);
    data->matchEmitNo[seqpos] = log (1. - model->pMatchEmit[state]);
    mean = model->matchMean[state];
    precision = model->matchPrecision[state];
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

/* Forward-backward matrix */

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_fb_matrix* mx;
  int n_events, matrix_cells;

  mx = SafeMalloc (sizeof (Seq_event_pair_fb_matrix));
  mx->data = new_seq_event_pair_data (model, seqlen, seq, events);

  n_events = mx->data->events->n_events;
  matrix_cells = mx->data->matrix_cells;

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
  counts->nNullEmitYes = 0;
  counts->nNullEmitNo = 0;
  counts->nullMoment0 = 0;
  counts->nullMoment1 = 0;
  counts->nullMoment2 = 0;
}

void inc_seq_event_pair_counts_from_fast5 (Seq_event_pair_counts* counts, Fast5_event_array* events) {
  int n_event, moves, new_moves, state;
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

double log_gaussian_density (double x, double mean, double precision, double log_precision) {
  double xz;
  xz = x - mean;
  return log_precision/2. - log_sqrt2pi - precision*xz*xz/2.;
}

double log_event_density (Fast5_event* event, double mean, double precision, double log_precision) {
  return event->ticks * (log_precision/2. - log_sqrt2pi - precision*mean*mean/2)
    - precision*(event->sumticks_cur_sq/2. - event->sumticks_cur*mean);
}

long double log_sum_exp (long double a, long double b) {
  double min, max, diff;
  if (a < b) { min = a; max = b; }
  else { min = b; max = a; }
  diff = max - min;
  return max + log (1. + exp(-diff));
}

void fill_seq_event_pair_fb_matrix_and_inc_counts (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event, state;
  long double mat, del;
  int idx, inputIdx, outputIdx;
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

  for (seqpos = order; seqpos <= seqlen; ++seqpos) {
    for (n_event = 0; n_event <= n_events; ++n_event) {
      idx = Seq_event_pair_index(seqpos,n_event);
      inputIdx = Seq_event_pair_index(seqpos-1,n_event);
      outputIdx = Seq_event_pair_index(seqpos,n_event-1);

      mat = matrix->fwdStart[n_event] + data->startEmitNo;   /* Start -> Match (input) */

      if (n_event > 0) {
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
    }

    matrix->fwdEnd = log_sum_exp
      (matrix->fwdEnd,
       matrix->fwdMatch[Seq_event_pair_index(seqpos,n_events)] + data->matchEmitNo[seqpos]);  /* Match -> End (input) */
  }

  /* fill backward & accumulate counts */
  for (n_event = n_events; n_event >= 0; --n_event) {
    event = n_event < n_events ? &data->events->event[n_event] : NULL;

    matrix->backStart[n_event] = -INFINITY;

    for (seqpos = seqlen; seqpos >= order; --seqpos) {
      state = data->state[seqpos];

      idx = Seq_event_pair_index(seqpos,n_event);
      inputIdx = Seq_event_pair_index(seqpos+1,n_event);
      outputIdx = Seq_event_pair_index(seqpos,n_event+1);

      if (n_event == n_events) {
	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEmitNo[seqpos],
			   0.,
			   matrix,
			   &counts->nMatchEmitNo[state], NULL,
			   NULL, NULL, NULL, NULL);   /* Match -> End (input) */
      } else {  /* n_event < n_events */
	mat = accum_count (-INFINITY,
			   matrix->fwdMatch[idx],
			   data->matchEmitYes[seqpos] * event->ticks
			   + data->matchEmitDensity[outputIdx],
			   matrix->backMatch[outputIdx],
			   matrix,
			   &counts->nMatchEmitYes[state], NULL,
			   event,
			   &counts->matchMoment0[state],
			   &counts->matchMoment1[state],
			   &counts->matchMoment2[state]);  /* Match -> Match (output) */
      }

      if (seqpos == seqlen) {
	del = -INFINITY;
      } else {
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
			   &counts->nMatchEmitNo[state],
			   &counts->nBeginDeleteYes,
			   NULL, NULL, NULL, NULL);  /* Match -> Delete (input) */
	mat = accum_count (mat,
			   matrix->fwdMatch[idx],
			   data->matchEmitNo[seqpos] + data->beginDeleteNo,
			   matrix->backMatch[inputIdx],
			   matrix,
			   &counts->nMatchEmitNo[state],
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

  /* update null model */
  for (n_event = 1; n_event <= n_events; ++n_event) {
    event = &data->events->event[n_event - 1];
    counts->nullMoment0 += event->ticks;
    counts->nullMoment1 += event->sumticks_cur;
    counts->nullMoment2 += event->sumticks_cur_sq;
  }
  counts->nNullEmitYes += n_events;
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

void optimize_seq_event_pair_model_for_counts (Seq_event_pair_model* model, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior) {
  int state;
  for (state = 0; state < model->states; ++state) {
    model->pMatchEmit[state] = (counts->nMatchEmitYes[state] + prior->nMatchEmitYes[state]) / (counts->nMatchEmitYes[state] + prior->nMatchEmitYes[state] + counts->nMatchEmitNo[state] + prior->nMatchEmitNo[state]);
    model->matchMean[state] = (counts->matchMoment1[state] + prior->matchMoment1[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]);
    model->matchPrecision[state] = 1. / (model->matchMean[state] * model->matchMean[state] - (counts->matchMoment2[state] + prior->matchMoment2[state]) / (counts->matchMoment0[state] + prior->matchMoment0[state]));
  }
  model->pBeginDelete = (counts->nBeginDeleteYes[state] + prior->nBeginDeleteYes[state]) / (counts->nBeginDeleteYes[state] + prior->nBeginDeleteYes[state] + counts->nBeginDeleteNo[state] + prior->nBeginDeleteNo[state]);
  model->pExtendDelete = (counts->nExtendDeleteYes[state] + prior->nExtendDeleteYes[state]) / (counts->nExtendDeleteYes[state] + prior->nExtendDeleteYes[state] + counts->nExtendDeleteNo[state] + prior->nExtendDeleteNo[state]);
  model->pStartEmit = (counts->nStartEmitYes + prior->nStartEmitYes) / (counts->nStartEmitYes + prior->nStartEmitYes + counts->nStartEmitNo + prior->nStartEmitNo);
  model->pNullEmit = (counts->nNullEmitYes + prior->nNullEmitYes) / (counts->nNullEmitYes + prior->nNullEmitYes + counts->nNullEmitNo + prior->nNullEmitNo);
  model->nullMean = (counts->nullMoment1 + prior->nullMoment1) / (counts->nullMoment0 + prior->nullMoment0);
  model->nullPrecision = 1. / (model->nullMean * model->nullMean - (counts->nullMoment2 + prior->nullMoment2) / (counts->nullMoment0 + prior->nullMoment0));
}

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Vector* event_arrays) {
  int iter;
  long double loglike, prev_loglike;
  for (iter = 0; iter < seq_evt_pair_EM_max_iterations; ++iter) {
    loglike = 0.;
    /* MORE TO GO HERE: loop through event_arrays */
    if (iter > 0 && prev_loglike != 0. && abs(prev_loglike) != INFINITY
	&& abs((loglike-prev_loglike)/prev_loglike) < seq_evt_pair_EM_min_fractional_loglike_increment)
      break;
    prev_loglike = loglike;
  }
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

void write_seq_event_pair_alignment_as_gff_cigar (Seq_event_pair_alignment* align, FILE* out) {
  int n, deleted;
  fprintf (out, "%s\t%s\t%s\t%d\t%d\t.\t.\t.\t%s=", align->seqname, gff3_source, gff3_feature, align->start_seqpos + 1, align->end_seqpos, gff3_gap_attribute);
  if (align->start_n_event > 0)
    fprintf (out, "I%d ", align->start_n_event);
  deleted = 0;
  for (n = 0; n < align->end_seqpos - align->start_seqpos; ++n) {
    ++deleted;
    if (align->events_at_pos[n] > 0) {
      fprintf (out, "D%d I%d ", deleted, align->events_at_pos[n]);
      deleted = 0;
    }
  }
  if (deleted)
    fprintf (out, "D%d ", deleted);
  fprintf (out, "\n");
}

Seq_event_pair_viterbi_matrix* new_seq_event_pair_viterbi_matrix (Seq_event_pair_model* model, int seqlen, char *seq, Fast5_event_array* events) {
  Seq_event_pair_viterbi_matrix* mx;
  int n_events, matrix_cells;

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
  int idx, inputIdx, outputIdx;
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
      inputIdx = Seq_event_pair_index(seqpos-1,n_event);
      outputIdx = Seq_event_pair_index(seqpos,n_event-1);

      mat = matrix->vitStart[n_event] + data->startEmitNo;   /* Start -> Match (input) */

      if (n_event > 0) {
	event = &data->events->event[n_event - 1];
	mat = max_func (mat,
			matrix->vitMatch[outputIdx]
			+ data->matchEmitYes[seqpos] * event->ticks
			+ data->matchEmitDensity[idx]);     /* Match -> Match (output) */
      }

      if (seqpos == order) {
	del = -INFINITY;
      } else {  /* seqpos > order */
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
}

Seq_event_pair_alignment* get_seq_event_pair_viterbi_matrix_traceback (Seq_event_pair_viterbi_matrix* matrix) {
  Seq_event_pair_model* model;
  Seq_event_pair_data* data;
  int seqlen, n_events, order, seqpos, n_event, end_seqpos, start_seqpos, start_n_event, zero, n, k;
  long double loglike;
  int idx, inputIdx, outputIdx;
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

  zero = 0;

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
  VectorPushBack (events_emitted, &zero);

  while (state != Start) {
    idx = Seq_event_pair_index(seqpos,n_event);
    inputIdx = Seq_event_pair_index(seqpos-1,n_event);
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
      state = Match;  /* which it already is */
      --n_event;
      ++*((int*) VectorBack (events_emitted));
      break;

    case MatchMatchIn:
      state = Match;  /* which it already is */
      --seqpos;
      VectorPushBack (events_emitted, &zero);
      break;

    case MatchDeleteIn:
      state = Match;
      --seqpos;
      VectorPushBack (events_emitted, &zero);
      break;

    case DeleteDeleteIn:
      state = Delete;  /* which it already is */
      --seqpos;
      VectorPushBack (events_emitted, &zero);
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
  align->start_seqpos = start_seqpos;
  align->end_seqpos = end_seqpos;
  align->start_n_event = start_n_event;
  align->events_at_pos = SafeMalloc (VectorSize(events_emitted) * sizeof(int));
  for (n = VectorSize(events_emitted) - 1, k = 0; n >= 0; --n, ++k)
    align->events_at_pos[k] = *((int*) VectorGet (events_emitted, n));

  deleteVector (events_emitted);

  return align;
}
