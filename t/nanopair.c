#include <locale.h>
#include <stdarg.h>

#include "../src/seqevtpair.h"
#include "../src/basecaller.h"
#include "../src/logsumexp.h"
#include "../src/logger.h"

const char* help_message = 
  "Usage: nanopair {modelseed,eventseed,priorseed,normalize,count,train,align,squiggle} <args>\n"
  "\n"
  " nanopair modelseed -fast5 <read.fast5>  >params.xml\n"
  "  (to parameterize a model from the HMM in a FAST5 file)\n"
  "\n"
  " nanopair eventseed -fast5 <read.fast5> [-fast5 <read2.fast5> ...]  >params.xml\n"
  "  (to parameterize a model from the basecalled event data in a FAST5 file)\n"
  "\n"
  " nanopair priorseed  >params.xml\n"
  "  (to parameterize a model with the mode of the prior distribution)\n"
  "\n"
  " nanopair normalize -in <in.fast5> -out <out.fast5>\n"
  "  (to normalize events in a FAST5 file)\n"
  "\n"
  " nanopair squiggle -fast5 <read.fast5>  >plot.svg\n"
  "  (to make a crude squiggle plot)\n"
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
  " nanopair basecall -params <params.xml> -fast5 <read.fast5> [more fast5...]  >seqs.fasta\n"
  "  (to base-call FAST5 reads via the Viterbi algorithm)\n"
  "\n"
  "\n"
  "Options (some are only available in certain command modes):\n"
  " -verbose, -vv, -vvv, -v4, etc.\n"
  " -log <function_name>\n"
  "                 various levels of logging\n"
  " -bothstrands    use forward & reverse complement of FASTA (default is forward only)\n"
  " -allvsall       train all reads on all references, instead of trying to pair them up\n"
  " -mininc         minimum fractional log-likelihood increment for EM to proceed\n"
  " -maxiter        maximum number of iterations of EM\n"
  " -pseudo {[no_]event,tick,skip,delete,extend,emit} <count>\n"
  "                 override various pseudocounts from the command-line\n"
  "\n"
  "For 'align', 'train' & 'count' commands, in place of '-params <params.xml>',\n"
  "can use '-eventseed', '-priorseed' or '-modelseed'.\n";

#define MODEL_ORDER 5

typedef enum SeedFlag { EventSeed, ModelSeed, PriorSeed, ParamFile } SeedFlag;

typedef struct Nanopair_args_str {
  Seq_event_pair_model *params;
  StringDoubleMap *prior;
  Kseq_container *seqs;
  StringVector *fast5_filenames;
  const char *fast5inFilename, *fast5outFilename;
  Vector *event_arrays;
  SeedFlag seedFlag;
  int model_order;
  Logger *logger;
  Seq_event_pair_config config;
} Nanopair_args;

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

void copy_prior (Seq_event_pair_model *model, StringDoubleMap *prior) {
  for (StringDoubleMapNode *pseudoNode = RBTreeFirst(prior);
       !RBTreeIteratorFinished(prior,pseudoNode);
       pseudoNode = RBTreeSuccessor(prior,pseudoNode)) {
    StringDoubleMapSet (model->prior, (char*) pseudoNode->key, *(double*)pseudoNode->value);
  }
}

void get_params (SeedFlag seedFlag, StringVector* fast5_filenames, Vector* event_arrays, Seq_event_pair_model** modelPtr, int order, Logger* logger, Seq_event_pair_config* config, StringDoubleMap *prior) {
  if (seedFlag == EventSeed)
    *modelPtr = eventseed_params (event_arrays);
  else if (seedFlag == ModelSeed)
    *modelPtr = seed_params (StringVectorGet (fast5_filenames, 0));
  else if (seedFlag == PriorSeed)
    *modelPtr = new_seq_event_pair_model_from_prior (order, prior);
  else
    Assert (*modelPtr != NULL, "No model parameters specified");

  (*modelPtr)->logger = logger;
  (*modelPtr)->config = *config;

  copy_prior (*modelPtr, prior);
}

Kseq_container* init_seqs (char* filename) {
  Kseq_container* seqs = init_kseq_container (filename);
  if (seqs == NULL)
    (void) help_failure ("Couldn't open FASTA file %s", filename);
  else
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
int parse_pseudo (int* argcPtr, char*** argvPtr, StringDoubleMap *prior) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-pseudo") == 0) {
      Assert (*argcPtr > 2, "-pseudo must have two arguments");
      StringDoubleMapSet (prior, (*argvPtr)[1], atof((*argvPtr)[2]));
      *argvPtr += 3;
      *argcPtr -= 3;
      return 1;
    }
  }
  return 0;
}

int parse_params (int* argcPtr, char*** argvPtr, SeedFlag* seedFlag, Seq_event_pair_model** modelPtr, int* orderPtr, StringDoubleMap *prior) {
  if (*argcPtr > 0) {
    if (strcmp (**argvPtr, "-eventseed") == 0) {
      *seedFlag = EventSeed;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    } else if (strcmp (**argvPtr, "-modelseed") == 0) {
      *seedFlag = ModelSeed;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    } else if (strcmp (**argvPtr, "-priorseed") == 0) {
      *seedFlag = PriorSeed;
      ++*argvPtr;
      --*argcPtr;
      return 1;
    } else if (strcmp (**argvPtr, "-order") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      *orderPtr = atoi ((*argvPtr)[1]);
      *argvPtr += 2;
      *argcPtr -= 2;
    } else if (strcmp (**argvPtr, "-params") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      *modelPtr = init_params ((*argvPtr)[1]);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    } else
      return parse_pseudo (argcPtr, argvPtr, prior);
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
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
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
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
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
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
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

int parse_dp (int* argcPtr, char*** argvPtr, Nanopair_args* npargPtr) {
  return parse_params (argcPtr, argvPtr, &npargPtr->seedFlag, &npargPtr->params, &npargPtr->model_order, npargPtr->prior)
    || parse_fast5 (argcPtr, argvPtr, npargPtr->fast5_filenames)
    || parse_seq (argcPtr, argvPtr, &npargPtr->seqs)
    || parseLogArgs (argcPtr, argvPtr, npargPtr->logger)
    || parse_seq_event_pair_config_general (argcPtr, argvPtr, &npargPtr->config)
    || parse_unknown (*argcPtr, *argvPtr);
}

/* main */
int main (int argc, char** argv) {
  int i, j;

  Seq_event_pair_counts *counts = NULL;
  xmlChar *xml_params = NULL;

  Nanopair_args npargs;
  npargs.params = NULL;
  npargs.seqs = NULL;
  npargs.fast5_filenames = newStringVector();
  npargs.event_arrays = newVector (NullCopyFunction, delete_fast5_event_array_null, NullPrintFunction);
  npargs.seedFlag = ParamFile;
  npargs.fast5inFilename = NULL;
  npargs.fast5outFilename = NULL;
  npargs.model_order = 5;
  npargs.prior = new_seq_event_pair_model_default_prior();

  npargs.logger = newLogger();
  init_seq_event_pair_config (&npargs.config);
  
  init_log_sum_exp_lookup();

  if (argc < 2)
    return help_failure ("Please specify a command.");

  const char* command = argv[1];
  argv += 2;
  argc -= 2;

  if (strcmp (command, "modelseed") == 0) {
    /* MODELSEED: initialize emit parameters from model in a FAST5 read file */
    while (parse_fast5_filename (&argc, &argv, &npargs.fast5inFilename)
	   || parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_pseudo (&argc, &argv, npargs.prior)
	   || parse_unknown (argc, argv))
      { }

    Assert (npargs.fast5inFilename != NULL, "No FAST5 file specified");

    /* initialize model */
    npargs.params = seed_params (npargs.fast5inFilename);
    if (npargs.params == NULL)
      return EXIT_FAILURE;

    copy_prior (npargs.params, npargs.prior);

    /* output model */
    xml_params = convert_seq_event_pair_model_to_xml_string (npargs.params);
    fprintf (stdout, "%s", (char*) xml_params);

    /* MORE SHOULD GO HERE: we need to set
       pBeginDelete, pExtendDelete,
       pStartEmit,
       pNullEmit, nullMean, nullPrecision
    */

  } else if (strcmp (command, "eventseed") == 0) {
    /* EVENTSEED: initialize parameters from base-called reads */
    while (parse_fast5 (&argc, &argv, npargs.fast5_filenames)
	   || parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_pseudo (&argc, &argv, npargs.prior)
	   || parse_unknown (argc, argv))
      { }

    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);

    /* initialize model */
    npargs.params = eventseed_params (npargs.event_arrays);
    if (npargs.params == NULL)
      return EXIT_FAILURE;

    copy_prior (npargs.params, npargs.prior);

    /* output model */
    write_params (npargs.params);

  } else if (strcmp (command, "priorseed") == 0) {
    /* PRIORSEED: initialize parameters from prior */
    while (parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_pseudo (&argc, &argv, npargs.prior)
	   || parse_unknown (argc, argv))
      { }

    /* initialize model */
    npargs.params = new_seq_event_pair_model_from_prior (npargs.model_order, npargs.prior);
    if (npargs.params == NULL)
      return EXIT_FAILURE;

    /* output model */
    write_params (npargs.params);

  } else if (strcmp (command, "normalize") == 0) {
    /* NORMALIZE: write out normalized fast5 file */
    while (parse_normalize (&argc, &argv, &npargs.fast5inFilename, &npargs.fast5outFilename)
	   || parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_unknown (argc, argv))
      { }

    Assert (npargs.fast5inFilename != NULL && npargs.fast5outFilename != NULL, "Both input & output FAST5 filenames must be specified");

    /* normalize */
    normalize_fast5 (npargs.fast5inFilename, npargs.fast5outFilename);

  } else if (strcmp (command, "count") == 0) {
    /* COUNT: single E-step of Baum-Welch/EM algorithm */
    while (parse_dp (&argc, &argv, &npargs)
	   || parse_seq_event_pair_config_training (&argc, &argv, &npargs.config))
      { }

    Assert (npargs.seqs != NULL, "Reference sequences not specified");
    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);
    get_params (npargs.seedFlag, npargs.fast5_filenames, npargs.event_arrays, &npargs.params, npargs.model_order, npargs.logger, &npargs.config, npargs.prior);

    /* get counts */
    counts = get_seq_event_pair_counts (npargs.params, npargs.seqs, npargs.event_arrays);

    /* output counts */
    write_counts (counts);

  } else if (strcmp (command, "train") == 0) {
    /* TRAIN: Baum-Welch/EM algorithm */
    while (parse_dp (&argc, &argv, &npargs)
	   || parse_seq_event_pair_config_training (&argc, &argv, &npargs.config))
      { }

    Assert (npargs.seqs != NULL, "Reference sequences not specified");
    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);
    get_params (npargs.seedFlag, npargs.fast5_filenames, npargs.event_arrays, &npargs.params, npargs.model_order, npargs.logger, &npargs.config, npargs.prior);

    /* do Baum-Welch */
    fit_seq_event_kmer_model (npargs.params, npargs.seqs);
    fit_seq_event_null_model (npargs.params, npargs.event_arrays);
    fit_seq_event_pair_model (npargs.params, npargs.seqs, npargs.event_arrays);

    /* output model */
    write_params (npargs.params);

  } else if (strcmp (command, "align") == 0) {
    /* ALIGN: Viterbi algorithm */
    while (parse_dp (&argc, &argv, &npargs))
      { }

    Assert (npargs.seqs != NULL, "Reference sequences not specified");
    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);
    get_params (npargs.seedFlag, npargs.fast5_filenames, npargs.event_arrays, &npargs.params, npargs.model_order, npargs.logger, &npargs.config, npargs.prior);

    /* loop through sequences, FAST5 files */
    for (i = 0; i < npargs.seqs->n; ++i)
      for (j = 0; j < (int) VectorSize(npargs.event_arrays); ++j)
	print_seq_evt_pair_alignments_as_stockholm (npargs.params, npargs.seqs->len[i], npargs.seqs->seq[i], npargs.seqs->name[i], (Fast5_event_array*) VectorGet(npargs.event_arrays,j), stdout, 0.);

  } else if (strcmp (command, "basecall") == 0) {
    /* BASECALL: Viterbi algorithm on basecaller HMM */
    while (parse_params (&argc, &argv, &npargs.seedFlag, &npargs.params, &npargs.model_order, npargs.prior)
	   || parse_fast5 (&argc, &argv, npargs.fast5_filenames)
	   || parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_unknown (argc, argv))
      { }

    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);
    get_params (npargs.seedFlag, npargs.fast5_filenames, npargs.event_arrays, &npargs.params, npargs.model_order, npargs.logger, &npargs.config, npargs.prior);

    /* loop through FAST5 files */
    for (j = 0; j < (int) VectorSize(npargs.event_arrays); ++j) {
      char* seq = basecall_fast5_event_array (npargs.params, (Fast5_event_array*) VectorGet(npargs.event_arrays,j));
      printf (">%s\n%s\n", (char*) VectorGet(npargs.fast5_filenames,j), seq);
      SafeFree (seq);
    }

  } else if (strcmp (command, "squiggle") == 0) {
    /* SQUIGGLE: draw squiggle plot as SVG wrapped in HTML */
    while (parse_fast5_filename (&argc, &argv, &npargs.fast5inFilename)
	   || parse_params (&argc, &argv, &npargs.seedFlag, &npargs.params, &npargs.model_order, npargs.prior)
	   || parseLogArgs (&argc, &argv, npargs.logger)
	   || parse_unknown (argc, argv))
      { }

    init_event_arrays (npargs.fast5_filenames, npargs.event_arrays);
    get_params (npargs.seedFlag, npargs.fast5_filenames, npargs.event_arrays, &npargs.params, npargs.model_order, npargs.logger, &npargs.config, npargs.prior);

    /* loop through event arrays, print squiggle SVGs */
    printf ("<html>\n<title>Squiggle plot</title>\n<body>\n");
    for (j = 0; j < (int) VectorSize(npargs.event_arrays); ++j) {
      Fast5_event_array* events = (Fast5_event_array*) VectorGet(npargs.event_arrays,j);
      xmlChar* svg = make_squiggle_svg (events, npargs.params);
      printf ("<p>\n%s\n<p>\n%s", events->name, svg);
      SafeFree (svg);
    }
    printf ("</body>\n</html>\n");

  } else if (strcmp (command, "help") == 0) {
    return help_failure (NULL);

  } else
    return help_failure ("Unrecognized command: %s", command);

  SafeFreeOrNull (xml_params);
  if (counts)
    delete_seq_event_pair_counts (counts);

  if (npargs.params)
    delete_seq_event_pair_model (npargs.params);

  if (npargs.seqs)
    free_kseq_container (npargs.seqs);

  deleteStringMap (npargs.prior);
  deleteVector (npargs.event_arrays);
  deleteStringVector (npargs.fast5_filenames);
  deleteLogger (npargs.logger);
  
  return EXIT_SUCCESS;
}
