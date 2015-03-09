#include <locale.h>
#include <stdarg.h>

#include "../src/seqevtpair.h"
#include "../src/logsumexp.h"

const char* help_message = 
  "Usage:\n"
  " nanopair seed <read.fast5>  > params.xml\n"
  " nanopair count <read.fast5> [<read2.fast5> ...]  > params.xml\n"
  " nanopair train <params.xml> <refs.fasta> <read.fast5> [<read2.fast5> ...]  > newparams.xml\n"
  " nanopair align <params.xml> <refs.fasta> <read.fast5> [<read2.fast5> ...]  > hits.gff\n";

#define MODEL_ORDER 5

int help_failure (char* warning, ...) {
  va_list argptr;
  if (warning) {
    va_start (argptr, warning);
    fprintf (stderr, "\n");
    vfprintf (stderr, warning, argptr);
    va_end (argptr);
  }
  fprintf (stderr, "\n\n%s\n", help_message);
  return EXIT_FAILURE;
}

void delete_fast5_event_array_null (void* ptr) {
  delete_fast5_event_array ((Fast5_event_array*) ptr);
}

Vector* init_fast5_event_array_vector (int argc, char** argv) {
  Vector* vec;
  Fast5_event_array* events;
  vec = newVector (NullCopyFunction, delete_fast5_event_array_null, NullPrintFunction);
  for (; argc > 0; --argc, ++argv) {
    events = read_fast5_event_array (*argv, DefaultFast5TickLength);
    if (events == NULL) {
      deleteVector (vec);  /* automatically calls delete_fast5_event_array */
      (void) help_failure ("Couldn't open FAST5 file %s", *argv);
      return NULL;
    }
    VectorPushBack (vec, events);
  }
  if (VectorSize(vec) == 0) {
    deleteVector (vec);
    (void) help_failure ("Couldn't open any FAST5 files");
    return NULL;
  }
  return vec;
}

Seq_event_pair_model* init_params (char* filename) {
  xmlChar* xml_params = (xmlChar*) readFileAsString (filename);
  if (xml_params == NULL) {
    (void) help_failure ("Parameter file %s not found", filename);
    return NULL;
  }
  Seq_event_pair_model* params = new_seq_event_pair_model_from_xml_string ((char*) xml_params);
  SafeFree (xml_params);
  return params;
}

Kseq_container* init_seqs (char* filename) {
  Kseq_container* seqs = init_kseq_container (filename);
  if (seqs == NULL)
    (void) help_failure ("Couldn't open FASTA file %s", filename);
  validate_kseq_container (seqs, dna_alphabet, stderr);
  return seqs;
}

void write_params (Seq_event_pair_model *params) {
  xmlChar *xml_params = convert_seq_event_pair_model_to_xml_string (params);
  fprintf (stdout, "%s", (char*) xml_params);
  SafeFree (xml_params);
}

int main (int argc, char** argv) {
  int i, j;
  Vector *event_arrays;
  Kseq_container *seqs;
  Seq_event_pair_model *params;
  xmlChar *xml_params;

  init_log_sum_exp_lookup();

  if (argc < 2)
    return help_failure ("Please specify a command");

  else if (strcmp (argv[1], "seed") == 0) {
    /* SEED: initialize emit parameters from model in a FAST5 read file */
    if (argc != 3)
      return help_failure ("For count-seeding, please specify a FAST5 file");

    /* initialize model */
    params = new_seq_event_pair_model (MODEL_ORDER);
    copy_seq_event_model_params_from_fast5 (params, argv[2]);

    /* output model */
    xml_params = convert_seq_event_pair_model_to_xml_string (params);
    fprintf (stdout, "%s", (char*) xml_params);

    /* MORE SHOULD GO HERE: we need to set
       pBeginDelete, pExtendDelete,
       pStartEmit,
       pNullEmit, nullMean, nullPrecision
    */

    /* free memory */
    SafeFree (xml_params);
    delete_seq_event_pair_model (params);

  } else if (strcmp (argv[1], "count") == 0) {
    /* COUNT: initialize parameters from base-called reads */
    if (argc < 3)
      return help_failure ("For count-seeding, please specify at least one FAST5 file");

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 2, argv + 2);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* initialize model */
    params = new_seq_event_pair_model (MODEL_ORDER);
    optimize_seq_event_model_for_events (params, event_arrays);

    /* output model */
    write_params (params);

    /* free memory */
    delete_seq_event_pair_model (params);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "train") == 0) {
    /* TRAINING */
    if (argc < 5)
      return help_failure ("For training, please specify parameters, a FASTA reference sequence and at least one FAST5 file");

    /* read parameters */
    params = init_params (argv[2]);
    if (params == NULL)
      return EXIT_FAILURE;

    /* read sequences */
    seqs = init_seqs (argv[3]);
    if (seqs == NULL)
      return EXIT_FAILURE;

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 4, argv + 4);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* do Baum-Welch */
    fit_seq_event_pair_model (params, seqs, event_arrays);

    /* output model */
    write_params (params);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "align") == 0) {
    /* ALIGNMENT */
    if (argc < 5)
      return help_failure ("For alignment, please specify parameters, a FASTA reference sequence and at least one FAST5 read file");

    /* read parameters */
    params = init_params (argv[2]);
    if (params == NULL)
      return EXIT_FAILURE;

    /* read sequences */
    seqs = init_seqs (argv[3]);
    if (seqs == NULL)
      return EXIT_FAILURE;

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 4, argv + 4);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* loop through sequences, FAST5 files */
    for (i = 0; i < seqs->n; ++i)
      for (j = 0; j < (int) VectorSize(event_arrays); ++j)
	print_seq_evt_pair_alignments_as_gff_cigar (params, seqs->len[i], seqs->seq[i], seqs->name[i], (Fast5_event_array*) VectorGet(event_arrays,j), stdout, 0.);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "help") == 0) {
    return help_failure (NULL);

  } else
    return help_failure ("Unrecognized command: %s", argv[1]);

  return EXIT_SUCCESS;
}
