#include <locale.h>
#include <stdarg.h>

#include "../src/seqevtpair.h"
#include "../src/logsumexp.h"

const char* help_message = 
  "Usage: nanopair {seed,eventseed,count,train,align} <args>\n"
  "\n"
  " nanopair seed <read.fast5>  >params.xml\n"
  "  (to parameterize a model from the HMM in a FAST5 file)\n"
  "\n"
  " nanopair eventseed <read.fast5> [<read2.fast5> ...]  >params.xml\n"
  "  (to parameterize a model from the basecalled event data in a FAST5 file)\n"
  "\n"
  " nanopair count <params.xml> <refs.fasta> <read.fast5> [...]  >counts.xml\n"
  "  (to calculate expected counts under the posterior distribution,\n"
  "   as summary statistics or for distributed EM updates)\n"
  "\n"
  " nanopair train <params.xml> <refs.fasta> <read.fast5> [...]  >newparams.xml\n"
  "  (to re-parameterize a model via the Baum-Welch/EM algorithm,\n"
  "   aligning one or more FAST5 reads to one or more reference sequences)\n"
  "\n"
  " nanopair align <params.xml> <refs.fasta> <read.fast5> [...]  >hits.gff\n"
  "  (to align FAST5 reads to reference sequences via the Viterbi algorithm)\n"
  "\n"
  "For 'align', 'train' & 'count' commands, to bypass the appropriate seed step,\n"
  "use '-eventseed' or '-seed' in place of <params.xml>.\n";

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
    if (events != NULL)
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

Seq_event_pair_model* seed_params (char* filename) {
  Seq_event_pair_model* params = new_seq_event_pair_model (MODEL_ORDER);
  if (init_seq_event_model_from_fast5 (params, filename) != 0) {
    delete_seq_event_pair_model (params);
    (void) help_failure ("Could not read parameters from %s", filename);
    return NULL;
  }
  return params;
}

Seq_event_pair_model* eventseed_params (Vector* event_arrays) {
  Seq_event_pair_model* params = new_seq_event_pair_model (MODEL_ORDER);
  optimize_seq_event_model_for_events (params, event_arrays);
  return params;
}

Seq_event_pair_model* get_params (char** argv, Vector* event_arrays) {
  if (strcmp (argv[2], "-eventseed") == 0)
    return eventseed_params (event_arrays);
  if (strcmp (argv[2], "-seed") == 0)
    return seed_params (argv[4]);
  return init_params (argv[2]);
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

void write_counts (Seq_event_pair_counts *counts) {
  xmlChar *xml_counts = convert_seq_event_pair_counts_to_xml_string (counts);
  fprintf (stdout, "%s", (char*) xml_counts);
  SafeFree (xml_counts);
}

int main (int argc, char** argv) {
  int i, j;
  Vector *event_arrays;
  Kseq_container *seqs;
  Seq_event_pair_model *params;
  Seq_event_pair_counts *counts;
  xmlChar *xml_params;

  init_log_sum_exp_lookup();

  if (argc < 2)
    return help_failure ("Please specify a command.");

  else if (strcmp (argv[1], "seed") == 0) {
    /* SEED: initialize emit parameters from model in a FAST5 read file */
    if (argc != 3)
      return help_failure ("For seeding, please specify a FAST5 file");

    /* initialize model */
    params = seed_params (argv[2]);
    if (params == NULL)
      return EXIT_FAILURE;

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

  } else if (strcmp (argv[1], "eventseed") == 0) {
    /* EVENTSEED: initialize parameters from base-called reads */
    if (argc < 3)
      return help_failure ("For event-based seeding, please specify at least one FAST5 file");

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 2, argv + 2);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* initialize model */
    params = eventseed_params (event_arrays);
    if (params == NULL)
      return EXIT_FAILURE;

    /* output model */
    write_params (params);

    /* free memory */
    delete_seq_event_pair_model (params);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "count") == 0) {
    /* COUNT: single E-step of Baum-Welch/EM algorithm */
    if (argc < 5)
      return help_failure ("To calculate summary counts please specify parameters, a FASTA reference sequence and at least one FAST5 file");

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 4, argv + 4);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* read sequences */
    seqs = init_seqs (argv[3]);
    if (seqs == NULL)
      return EXIT_FAILURE;

    /* read parameters */
    params = get_params (argv, event_arrays);
    if (params == NULL)
      return EXIT_FAILURE;

    /* get counts */
    counts = get_seq_event_pair_counts (params, seqs, event_arrays);

    /* output counts */
    write_counts (counts);

    /* free memory */
    delete_seq_event_pair_counts (counts);
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "train") == 0) {
    /* TRAIN: Baum-Welch/EM algorithm */
    if (argc < 5)
      return help_failure ("For training, please specify parameters, a FASTA reference sequence and at least one FAST5 file");

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 4, argv + 4);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* read sequences */
    seqs = init_seqs (argv[3]);
    if (seqs == NULL)
      return EXIT_FAILURE;

    /* read parameters */
    params = get_params (argv, event_arrays);
    if (params == NULL)
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
    /* ALIGN: Viterbi algorithm */
    if (argc < 5)
      return help_failure ("For alignment, please specify parameters, a FASTA reference sequence and at least one FAST5 read file");

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 4, argv + 4);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* read sequences */
    seqs = init_seqs (argv[3]);
    if (seqs == NULL)
      return EXIT_FAILURE;

    /* read parameters */
    params = get_params (argv, event_arrays);
    if (params == NULL)
      return EXIT_FAILURE;

    /* loop through sequences, FAST5 files */
    for (i = 0; i < seqs->n; ++i)
      for (j = 0; j < (int) VectorSize(event_arrays); ++j)
	print_seq_evt_pair_alignments_as_stockholm (params, seqs->len[i], seqs->seq[i], seqs->name[i], (Fast5_event_array*) VectorGet(event_arrays,j), stdout, 0.);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "squiggle") == 0) {
    /* SQUIGGLE: draw squiggle plot as SVG wrapped in HTML */
    if (argc < 4)
      return help_failure ("For squiggle plots, please specify parameters and a FAST5 read file");

    /* read FAST5 file */
    event_arrays = init_fast5_event_array_vector (argc - 3, argv + 3);
    if (event_arrays == NULL)
      return EXIT_FAILURE;

    /* read parameters */
    params = (strcmp (argv[2], "-seed") == 0)
      ? seed_params (argv[3])
      : init_params (argv[2]);
      
    if (params == NULL)
      return EXIT_FAILURE;

    /* loop through event arrays, print squiggle SVGs */
    printf ("<html>\n<title>Squiggle plot</title>\n<body>\n");
    for (j = 0; j < (int) VectorSize(event_arrays); ++j) {
      Fast5_event_array* events = (Fast5_event_array*) VectorGet(event_arrays,j);
      xmlChar* svg = make_squiggle_svg (events, params);
      printf ("<p>\n%s\n<p>\n%s", events->name, svg);
      SafeFree (svg);
    }
    printf ("</body>\n</html>\n");

    /* free memory */
    delete_seq_event_pair_model (params);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "help") == 0) {
    return help_failure (NULL);

  } else
    return help_failure ("Unrecognized command: %s", argv[1]);

  return EXIT_SUCCESS;
}
