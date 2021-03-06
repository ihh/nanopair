#include <locale.h>
#include <stdio.h>

#include "../src/seqevtpair.h"
#include "../src/logsumexp.h"

#define N_EVENTS 2
#define MODEL_ORDER 5

void delete_fast5_event_array_null (void* ptr) {
  delete_fast5_event_array ((Fast5_event_array*) ptr);
}

int main (int argc, char** argv) {
  Fast5_event_array *events;
  Fast5_event event[N_EVENTS] = { { .mean = 100, .stdv = 2, .length = .5, .model_state = "AAAAA", .mp_model_state = "AAAAA", .move = 0 },
				  { .mean = 102, .stdv = 1, .length = .4, .model_state = "AAAAC", .mp_model_state = "AAAAC", .move = 1 } };
  Vector* event_arrays;
  Kseq_container *seqs;
  Seq_event_pair_model *params;
  xmlChar *xml_params;
  int n_event;

  setlocale(LC_ALL, "C");
  init_log_sum_exp_lookup();

  /* create Fast5_event_array */
  events = alloc_fast5_event_array (-1, N_EVENTS);

  SafeFree (events->event);

  for (n_event = 0; n_event < N_EVENTS; ++n_event)
    fast5_event_calc_moments (event + n_event, DefaultFast5TickLength);
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
  seqs->seq[0] = StringCopy ("ACAGCT");
  seqs->len[0] = strlen (seqs->seq[0]);

  validate_kseq_container (seqs, dna_alphabet, stderr);

  /* main program goes here */

  /* initialize model */
  params = new_seq_event_pair_model (MODEL_ORDER);
  optimize_seq_event_model_for_events (params, event_arrays);

  /* output initial model */
  FILE *initParamsFile = fopen ("initparams.xml", "w");
  xml_params = convert_seq_event_pair_model_to_xml_string (params);
  fprintf (initParamsFile, "%s", (char*) xml_params);
  SafeFree (xml_params);
  fclose (initParamsFile);

  /* do Baum-Welch */
  fit_seq_event_pair_model (params, seqs, event_arrays, 0);

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
