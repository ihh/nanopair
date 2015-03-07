#include <locale.h>
#include <stdio.h>

#include "../src/seqevtpair.h"

#define N_EVENTS 2
#define MODEL_ORDER 5

void delete_fast5_event_array_null (void* ptr) {
  delete_fast5_event_array ((Fast5_event_array*) ptr);
}

int main (int argc, char** argv) {
  Fast5_event_array *events;
  Fast5_event event[N_EVENTS] = { { .mean = 100, .stdv = 2, .length = .5, .model_state = "AAAAA", .mp_model_state = "AAAAA", .move = 0, .raw = 0 },
				  { .mean = 102, .stdv = 1, .length = .4, .model_state = "AAAAC", .mp_model_state = "AAAAC", .move = 1, .raw = 100 } };
  Vector* event_arrays;
  Kseq_container *seqs;
  Seq_event_pair_model *params;
  xmlChar *xml_params;

  setlocale(LC_ALL, "C");

  /* create Fast5_event_array */
  events = alloc_fast5_event_array (-1, N_EVENTS, DefaultFast5TickLength);

  SafeFree (events->event);
  events->event = event;

  event_arrays = newVector (NullCopyFunction, delete_fast5_event_array_null, NullPrintFunction);
  VectorPushBack (event_arrays, events);

  /* populate Kseq_container */
  seqs = (Kseq_container*) SafeMalloc (sizeof(Kseq_container));
  seqs->n = 1;
  seqs->name = (char**) SafeCalloc (1, sizeof(char*));
  seqs->len = (int*) SafeCalloc (1, sizeof(int));
  seqs->seq = (char**) SafeCalloc (1, sizeof(char*));

  seqs->name[0] = StringCopy ("tiny");
  seqs->seq[0] = StringCopy ("ACAGCTAGCTAGCTAG");
  seqs->len[0] = strlen (seqs->seq[0]);

  validate_kseq_container (seqs, dna_alphabet, stderr);

  /* main program goes here */

  /* do Baum-Welch */
  params = new_seq_event_pair_model (MODEL_ORDER);
  fit_seq_event_pair_model (params, seqs, event_arrays);

  /* output model */
  xml_params = convert_seq_event_pair_model_to_xml_string (params);
  fprintf (stdout, "%s", (char*) xml_params);

  /* free memory */
  SafeFree (xml_params);
  delete_seq_event_pair_model (params);

  /* clean up */
  free_kseq_container (seqs);

  events->event = NULL;  /* call this before deleteVector, to avoid deleting static Fast5_event array */
  deleteVector (event_arrays);

  return EXIT_SUCCESS;
}
