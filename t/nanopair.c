#include <locale.h>
#include <stdarg.h>

#include "../src/seqevtpair.h"
#include "../src/logsumexp.h"

const char* help_message = 
  "Usage: nanopair {seed,eventseed,normalize,count,train,align} <args>\n"
  "\n"
  " nanopair seed -fast5 <read.fast5>  >params.xml\n"
  "  (to parameterize a model from the HMM in a FAST5 file)\n"
  "\n"
  " nanopair eventseed -fast5 <read.fast5> [-fast5 <read2.fast5> ...]  >params.xml\n"
  "  (to parameterize a model from the basecalled event data in a FAST5 file)\n"
  "\n"
  " nanopair normalize -in <in.fast5> -out <out.fast5>\n"
  "  (to normalize events in a FAST5 file)\n"
  "\n"
  " nanopair count -params <params.xml> -fasta <refs.fasta> -fast5 <read.fast5> [more fast5...]  >counts.xml\n"
  "  (to calculate expected counts under the posterior distribution,\n"
  "   as summary statistics or for distributed EM updates)\n"
  "\n"
  " nanopair train -params <params.xml> -fasta <refs.fasta> -fast5 <read.fast5> [more fast5...]  >newparams.xml\n"
  "  (to re-parameterize a model via the Baum-Welch/EM algorithm,\n"
  "   aligning one or more FAST5 reads to one or more reference sequences)\n"
  "\n"
  " nanopair align -params <params.xml> -fasta <refs.fasta> -fast5 <read.fast5> [more fast5...]  >hits.gff\n"
  "  (to align FAST5 reads to reference sequences via the Viterbi algorithm)\n"
  "\n"
  "For 'align', 'train' & 'count' commands, to bypass the appropriate seed step,\n"
  "use '-eventseed' or '-seed' in place of '-params <params.xml>'.\n";

#define MODEL_ORDER 5

typedef enum SeedFlag { EventSeed, Seed, Params } SeedFlag;

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

Seq_event_pair_model* seed_params (const char* filename) {
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

void get_params (SeedFlag seedFlag, StringVector* fast5_filenames, Vector* event_arrays, Seq_event_pair_model** modelPtr) {
  if (seedFlag == EventSeed)
    *modelPtr = eventseed_params (event_arrays);
  else if (seedFlag == Seed)
    *modelPtr = seed_params (StringVectorGet (fast5_filenames, 0));
  else
    Assert (*modelPtr != NULL, "No model parameters specified");
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

void normalize_fast5 (const char* fnIn, const char* fnOut) {
  /* read, write & destroy event array */
  Fast5_event_array* event_array = read_fast5_event_array (fnIn);
  write_fast5_event_array (event_array, fnOut);
  delete_fast5_event_array (event_array);
}

void init_event_arrays (StringVector* filenames, Vector* event_arrays) {
  for (int n = 0; n < (int) VectorSize(filenames); ++n) {
    const char* filename = StringVectorGet (filenames, n);
    Fast5_event_array* events = read_fast5_event_array (filename);
    Assert (events != NULL, "Couldn't read fast5 file %s", filename);
    VectorPushBack (event_arrays, events);
  }
  Assert (VectorSize (event_arrays) > 0, "No fast5 files specified");
}

/* parsers */
int parse_params (int* argcPtr, char*** argvPtr, SeedFlag* seedFlag, Seq_event_pair_model** modelPtr) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-eventseed") == 0) {
      *seedFlag = EventSeed;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    } else if (strcmp (**argvPtr, "-seed") == 0) {
      *seedFlag = EventSeed;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    } else if (strcmp (**argvPtr, "-params") == 0) {
      Assert (*argcPtr > 1, "-params must have an argument");
      *modelPtr = init_params ((*argvPtr)[1]);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    }
  }
  return 0;
}

int parse_normalize (int* argcPtr, char*** argvPtr, const char** inFilenamePtr, const char** outFilenamePtr) {
  if (*argcPtr > 0) {
    const int isIn = strcmp (**argvPtr, "-in") == 0;
    const int isOut = strcmp (**argvPtr, "-out") == 0;
    if (isIn || isOut) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      const char** filenamePtr = isIn ? inFilenamePtr : outFilenamePtr;
      Assert (*filenamePtr==0, "Can't specify multiple fast5 files with %s", **argvPtr);
      *filenamePtr = (*argvPtr)[1];
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    }
  }
  return 0;
}

int parse_fast5_filename (int* argcPtr, char*** argvPtr, const char** filenamePtr) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-fast5") == 0) {
      Assert (*argcPtr > 1, "-fast5 must have an argument");
      Assert (*filenamePtr==0, "Can't specify multiple fast5 files with this command");
      *filenamePtr = (*argvPtr)[1];
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    }
  }
  return 0;
}

int parse_fast5 (int* argcPtr, char*** argvPtr, StringVector* filenames) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-fast5") == 0) {
      Assert (*argcPtr > 1, "-fast5 must have an argument");
      VectorPushBack (filenames, StringNew ((*argvPtr)[1]));
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    }
  }
  return 0;
}

int parse_seq (int* argcPtr, char*** argvPtr, Kseq_container **seqsPtr) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-fasta") == 0) {
      Assert (*argcPtr > 1, "-fasta must have an argument");
      Assert (*seqsPtr == NULL, "Can't specify multiple fasta files");
      *seqsPtr = init_seqs ((*argvPtr)[1]);
      Assert (*seqsPtr != NULL, "Couldn't read fasta file %s", (*argvPtr)[1]);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    }
  }
  return 0;
}

int parse_unknown (int argc, char** argv) {
  if (argc > 0) {
    (void) help_failure ("Unknown option: %s", *argv);
    Abort ("Error parsing command-line options");
  }
  return 0;
}

int parse_dp (int* argcPtr, char*** argvPtr, SeedFlag* seedFlag, Seq_event_pair_model** modelPtr, StringVector* fast5_filenames, Kseq_container **seqsPtr) {
  return parse_params (argcPtr, argvPtr, seedFlag, modelPtr)
    || parse_fast5 (argcPtr, argvPtr, fast5_filenames)
    || parse_seq (argcPtr, argvPtr, seqsPtr)
    || parse_unknown (*argcPtr, *argvPtr);
}

/* main */
int main (int argc, char** argv) {
  int i, j;
  Seq_event_pair_counts *counts;
  xmlChar *xml_params;

  Seq_event_pair_model *params = NULL;
  Kseq_container *seqs = NULL;
  StringVector *fast5_filenames = newStringVector();
  Vector *event_arrays = newVector (NullCopyFunction, delete_fast5_event_array_null, NullPrintFunction);
  SeedFlag seedFlag = Params;
  const char* fast5inFilename = NULL;
  const char* fast5outFilename = NULL;

  init_log_sum_exp_lookup();

  if (argc < 2)
    return help_failure ("Please specify a command.");

  if (argc < 3)
    return help_failure ("All commands need arguments.");

  const char* command = argv[1];
  argv += 2;
  argc -= 2;

  if (strcmp (command, "seed") == 0) {
    /* SEED: initialize emit parameters from model in a FAST5 read file */
    while (parse_fast5_filename (&argc, &argv, &fast5inFilename)
	   || parse_unknown (argc, argv))
      { }

    Assert (fast5inFilename != NULL, "No FAST5 file specified");

    /* initialize model */
    params = seed_params (fast5inFilename);
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

  } else if (strcmp (command, "eventseed") == 0) {
    /* EVENTSEED: initialize parameters from base-called reads */
    while (parse_fast5 (&argc, &argv, fast5_filenames)
	   || parse_unknown (argc, argv))
      { }

    init_event_arrays (fast5_filenames, event_arrays);

    /* initialize model */
    params = eventseed_params (event_arrays);
    if (params == NULL)
      return EXIT_FAILURE;

    /* output model */
    write_params (params);

    /* free memory */
    delete_seq_event_pair_model (params);

  } else if (strcmp (command, "normalize") == 0) {
    /* NORMALIZE: write out normalized fast5 file */
    while (parse_normalize (&argc, &argv, &fast5inFilename, &fast5outFilename)
	   || parse_unknown (argc, argv))
      { }

    Assert (fast5inFilename != NULL && fast5outFilename != NULL, "Both input & output FAST5 filenames must be specified");

    /* normalize */
    normalize_fast5 (fast5inFilename, fast5outFilename);

  } else if (strcmp (command, "count") == 0) {
    /* COUNT: single E-step of Baum-Welch/EM algorithm */
    while (parse_dp (&argc, &argv, &seedFlag, &params, fast5_filenames, &seqs))
      { }

    Assert (seqs != NULL, "Reference sequences not specified");
    init_event_arrays (fast5_filenames, event_arrays);
    get_params (seedFlag, fast5_filenames, event_arrays, &params);

    /* get counts */
    counts = get_seq_event_pair_counts (params, seqs, event_arrays);

    /* output counts */
    write_counts (counts);

    /* free memory */
    delete_seq_event_pair_counts (counts);
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);

  } else if (strcmp (command, "train") == 0) {
    /* TRAIN: Baum-Welch/EM algorithm */
    while (parse_dp (&argc, &argv, &seedFlag, &params, fast5_filenames, &seqs))
      { }

    Assert (seqs != NULL, "Reference sequences not specified");
    init_event_arrays (fast5_filenames, event_arrays);
    get_params (seedFlag, fast5_filenames, event_arrays, &params);

    /* do Baum-Welch */
    fit_seq_event_pair_model (params, seqs, event_arrays);

    /* output model */
    write_params (params);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);

  } else if (strcmp (command, "align") == 0) {
    /* ALIGN: Viterbi algorithm */
    while (parse_dp (&argc, &argv, &seedFlag, &params, fast5_filenames, &seqs))
      { }

    Assert (seqs != NULL, "Reference sequences not specified");
    init_event_arrays (fast5_filenames, event_arrays);
    get_params (seedFlag, fast5_filenames, event_arrays, &params);

    /* loop through sequences, FAST5 files */
    for (i = 0; i < seqs->n; ++i)
      for (j = 0; j < (int) VectorSize(event_arrays); ++j)
	print_seq_evt_pair_alignments_as_stockholm (params, seqs->len[i], seqs->seq[i], seqs->name[i], (Fast5_event_array*) VectorGet(event_arrays,j), stdout, 0.);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);

  } else if (strcmp (command, "squiggle") == 0) {
    /* SQUIGGLE: draw squiggle plot as SVG wrapped in HTML */
    while (parse_fast5_filename (&argc, &argv, &fast5inFilename)
	   || parse_params (&argc, &argv, &seedFlag, &params)
	   || parse_unknown (argc, argv))
      { }

    init_event_arrays (fast5_filenames, event_arrays);
    get_params (seedFlag, fast5_filenames, event_arrays, &params);

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

  } else if (strcmp (command, "help") == 0) {
    return help_failure (NULL);

  } else
    return help_failure ("Unrecognized command: %s", command);

  deleteVector (event_arrays);
  deleteStringVector (fast5_filenames);

  return EXIT_SUCCESS;
}
