#ifndef KSEQCONTAINER_INCLUDED
#define KSEQCONTAINER_INCLUDED

#include "../kseq/kseq.h"

typedef struct Kseq_container {
  int n;
  char **name;
  char **seq;
} Kseq_container;

Kseq_container* init_kseq_container (const char* filename);
void free_kseq_container (Kseq_container* ksc);


#endif /* KSEQCONTAINER_INCLUDED */
