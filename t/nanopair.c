#include <locale.h>

#include "../src/seqevtpair.h"

const char* help_message = 
  "Usage:\n"
  " nanopair train ref.fasta run.fast5 [more.fast5 ...]   > paramfile.xml\n"
  " nanopair align -params paramfile.xml ref.fasta run.fast5 [more.fast5 ...]  > hits.gff\n";

#define MODEL_ORDER 5

int help_failure (char* warning) {
  if (warning)
    Warn ("%s\n", warning);
  Warn (help_message);
  return EXIT_FAILURE;
}

void delete_fast5_event_array_null (void* ptr) {
  delete_fast5_event_array ((Fast5_event_array*) ptr);
}

Vector* init_fast5_event_array_vector (int argc, char** argv) {
  Vector* vec;
  vec = newVector (NullCopyFunction, delete_fast5_event_array_null, NullPrintFunction);
  for (; argc > 0; --argc, ++argv)
    VectorPushBack (vec, read_fast5_event_array (*argv, DefaultFast5TickLength));
  return vec;
}

int main (int argc, char** argv) {
  int i, j;
  Vector *event_arrays;
  Kseq_container *seqs;
  Seq_event_pair_model *params;
  xmlChar *xml_params;

  if (argc < 2)
    return help_failure ("Please specify a command");

  else if (strcmp (argv[1], "train") == 0) {
    /* TRAINING */
    if (argc < 4)
      return help_failure ("For training, please specify a FASTA reference sequence and at least one FAST5 file");

    /* read sequences */
    seqs = init_kseq_container (argv[2]);
    validate_kseq_container (seqs, dna_alphabet, stderr);

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 3, argv + 3);

    /* do Baum-Welch */
    params = new_seq_event_pair_model (MODEL_ORDER);
    fit_seq_event_pair_model (params, seqs, event_arrays);

    /* output model */
    xml_params = convert_seq_event_pair_model_to_xml_string (params);
    fprintf (stdout, "%s", (char*) xml_params);

    /* free memory */
    SafeFree (xml_params);
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else if (strcmp (argv[1], "align") == 0) {
    /* ALIGNMENT */
    if (argc < 6)
      return help_failure ("For alignment, please specify parameters, a FASTA reference sequence and at least one FAST5 file");

    /* read parameters */
    if (strcmp (argv[2], "-params") != 0)
      return help_failure ("For alignment, please specify parameters using -params");
    xml_params = (xmlChar*) readFileAsString (argv[3]);
    if (xml_params == NULL)
      return help_failure ("Parameter file not found");
    params = new_seq_event_pair_model_from_xml_string ((char*) xml_params);
    SafeFree (xml_params);

    /* read sequences */
    seqs = init_kseq_container (argv[4]);
    validate_kseq_container (seqs, dna_alphabet, stderr);

    /* read FAST5 files */
    event_arrays = init_fast5_event_array_vector (argc - 5, argv + 5);

    /* loop through sequences, FAST5 files */
    for (i = 0; i < seqs->n; ++i)
      for (j = 0; j < (int) VectorSize(event_arrays); ++j)
	print_seq_evt_pair_alignments_as_gff_cigar (params, seqs->len[i], seqs->seq[i], seqs->name[i], (Fast5_event_array*) VectorGet(event_arrays,j), stdout, 0.);

    /* free memory */
    delete_seq_event_pair_model (params);
    free_kseq_container (seqs);
    deleteVector (event_arrays);  /* automatically calls delete_fast5_event_array */

  } else
    return help_failure("Unrecognized command");

  return EXIT_SUCCESS;
}
