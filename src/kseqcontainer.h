#ifndef KSEQCONTAINER_INCLUDED
#define KSEQCONTAINER_INCLUDED

#include <stdio.h>

#include "../kseq/kseq.h"

extern const char* dna_alphabet;

int tokenize (char c, const char* alphabet);

typedef struct Kseq_container {
  int n;
  char **name;
  int *len;
  char **seq;
} Kseq_container;

Kseq_container* init_kseq_container (const char* filename);
void free_kseq_container (Kseq_container* ksc);

int validate_kseq_container (Kseq_container* ksc, const char* alphabet, FILE* err);  /* returns number of bad chars */

char* new_revcomp_seq (const char* dna_seq, int len);

#endif /* KSEQCONTAINER_INCLUDED */
