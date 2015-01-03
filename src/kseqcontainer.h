#ifndef KSEQCONTAINER_INCLUDED
#define KSEQCONTAINER_INCLUDED

#include "../kseq/kseq.h"

extern const char* dna_alphabet;

int tokenize (char c, const char* alphabet);

typedef struct Kseq_container {
  int n;
  char **name;
  int *len;
  char **seq;
  int **dsq;
} Kseq_container;

Kseq_container* init_kseq_container (const char* filename, const char* alphabet);
void free_kseq_container (Kseq_container* ksc);


#endif /* KSEQCONTAINER_INCLUDED */
