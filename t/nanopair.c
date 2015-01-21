#include <locale.h>

#include "../src/seqevtpair.h"

void print_help() {
    printf ("Usage:\n"
	    "nanopair train <FASTA file> <Fast5 file> [<more Fast5 files...>]   > paramfile.xml\n"
	    "nanopair align <paramfile> <FASTA file> <Fast5file> [<more fast5 files...>]  > hits.gff\n");
}

Vector* init_fast5_event_array_vector (int argc, char** argv) {
  /* MORE TO GO HERE */
  return NULL;
}

int main (int argc, char** argv) {
  Vector* event_arrays;
  Kseq_container* seqs;
  Seq_event_pair_model* params;
  if (argc <= 1)
    print_help();
  else if (strcmp (argv[1], "train") == 0) {
    /* MORE TO GO HERE */
  } else if (strcmp (argv[1], "align") == 0) {
    /* MORE TO GO HERE */
  } else {
    printf ("Unrecognized command\n\n");
    print_help();
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
