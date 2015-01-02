#include <math.h>
#include <string.h>

#include "seqevtpair.h"
#include "xmlutil.h"
#include "xmlkeywords.h"

const char* dna_alphabet = "ACGT";

int base2token (char base);
char token2base (int token);

void encode_state_identifier (int state, int order, char* state_id);
int decode_state_identifier (char* state_id);

Seq_event_pair_model* new_seq_event_pair_model (int order) {
  Seq_event_pair_model *model;
  model = SafeMalloc (sizeof (Seq_event_pair_model));
  model->order = order;
  model->states = pow(4,order);
  model->pEmitCurrent = SafeMalloc (model->states * sizeof(double));
  model->currentMean = SafeMalloc (model->states * sizeof(double));
  model->currentPrecision = SafeMalloc (model->states * sizeof(double));
  return model;
}

void delete_seq_event_pair_model (Seq_event_pair_model* model) {
  SafeFree (model->pEmitCurrent);
  SafeFree (model->currentMean);
  SafeFree (model->currentPrecision);
  SafeFree (model);
}

int base2token (char base) {
  int tok;
  tok = strchr (dna_alphabet, toupper(base)) - dna_alphabet;
  return tok >= 4 ? -1 : tok;
}

char token2base (int token) {
  return tok < 0 || tok >= 4 ? "N" : dna_alphabet[tok];
}

void encode_state_identifier (int state, int order, char* state_id) {
  int k;
  for (k = 0; k < order; ++k, state /= 4)
    state_id[order - k - 1] = token2base (state % 4);
  state_id[order] = '\0';
}

int decode_state_identifier (char* state_id) {
  int k, token, p;
  for (token = 0, p = 1, k = 0; k < order; ++k, p *= 4)
    token += p * base2token (state_id[order - k - 1]);
  return token;
}

Seq_event_pair_model* new_seq_event_pair_model_from_xml_string (const char* xml) {
  xmlNode *node, *modelNode, *statesNode, *stateNode, *deleteNode, *currentNode, *startNode;
  Seq_event_pair_model *model;
  int state;
  node = xmlTreeFromString (xml);
  modelNode = CHILD(node,MODEL);
  model = new_seq_event_pair_model (CHILDINT(modelNode,ORDER));

  deleteNode = CHILD(stateNode,DELETE);
  model->pBeginDelete = CHILDFLOAT(deleteNode,BEGIN);
  model->pExtendDelete = CHILDFLOAT(deleteNode,EXTEND);

  statesNode = CHILD(modelNode,STATES);
  for (stateNode = statesNode->children; stateNode; stateNode = stateNode->next)
    if (MATCHES(stateNode,STATE)) {
      state = decode_state_identifier (model->order, CHILDINT(stateNode,ID));
      model->pEmitCurrent[state] = CHILDFLOAT(stateNode,EMIT);
      model->currentMean[state] = CHILDFLOAT(stateNode,MEAN);
      model->currentPrecision[state] = CHILDFLOAT(stateNode,PRECISION);
    }
  startNode = CHILD(modelNode,START);
  model->pStartEmitCurrent = CHILDFLOAT(startNode,EMIT);
  model->startCurrentMean = CHILDFLOAT(startNode,MEAN);
  model->startCurrentPrecision = CHILDFLOAT(startNode,PRECISION);
}

char* convert_seq_event_pair_model_to_xml_string (Seq_event_pair_model* model) {
  xmlTextWriterPtr writer;
  char* id;
  int state;

  id = SafeMalloc ((model->order + 1) * sizeof(char));

  writer = newXmlTextWriter();
  xmlTextWriterStartElement (writer, XMLPREFIX(MODEL));
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(ORDER), "%d", model->order);

  xmlTextWriterStartElement (writer, XMLPREFIX(DELETE));
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(BEGIN), "%g", model->pBeginDelete);
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(EXTEND), "%g", model->pExtendDelete);
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, XMLPREFIX(STATES));
  for (state = 0; state < model->states; ++state) {
    xmlTextWriterStartElement (writer, XMLPREFIX(STATE));
    encode_state_identifier (state, model->order, id);
    xmlTextWriterWriteFormatElement (writer, XMLPREFIX(ID), "%s", id);
    xmlTextWriterWriteFormatElement (writer, XMLPREFIX(EMIT), "%g", model->pEmitCurrent[state]);
    xmlTextWriterWriteFormatElement (writer, XMLPREFIX(MEAN), "%g", model->currentMean[state]);
    xmlTextWriterWriteFormatElement (writer, XMLPREFIX(PRECISION), "%g", model->currentPrecision[state]);
    xmlTextWriterEndElement (writer);
  }
  xmlTextWriterEndElement (writer);

  xmlTextWriterStartElement (writer, XMLPREFIX(START));
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(EMIT), "%g", model->pStartEmitCurrent);
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(MEAN), "%g", model->startCurrentMean);
  xmlTextWriterWriteFormatElement (writer, XMLPREFIX(PRECISION), "%g", model->startCurrentPrecision);
  xmlTextWriterEndElement (writer);

  xmlTextWriterEndElement (writer);
  return deleteXmlTextWriterLeavingText (writer);
}

/* Forward-backward matrix */

Seq_event_pair_fb_matrix* new_seq_event_pair_fb_matrix (Seq_event_pair_model* model, Kseq_container* seq, Fast5_event_array* events);  /* allocates only, does not fill */
void delete_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix);

Seq_event_pair_counts* new_seq_event_pair_counts (Seq_event_pair_model* model);
void delete_seq_event_pair_counts (Seq_event_pair_counts* counts);

void reset_seq_event_pair_counts (Seq_event_pair_counts* counts);
void fill_seq_event_pair_fb_matrix (Seq_event_pair_fb_matrix* matrix, Seq_event_pair_counts* counts);
void update_seq_event_pair_model (Seq_event_pair_model* matrix, Seq_event_pair_counts* counts, Seq_event_pair_counts* prior);

void fit_seq_event_pair_model (Seq_event_pair_model* model, Kseq_container* seq, Fast5_event_array* events, double minimum_fractional_log_likelihood_increase);
