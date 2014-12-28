#include <zlib.h>
#include <stdio.h>
#include "kseqcontainer.h"
#include "vector.h"

KSEQ_INIT(gzFile, gzread)

Kseq_container* init_kseq_container (const char* filename) {
  gzFile fp;
  kseq_t *seq;
  int l, i;

  StringVector *names, *seqs;
  Kseq_container* ksc;

  names = newStringVector();
  seqs = newStringVector();

  fp = gzopen(filename, "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    StringVectorPushBack(names,seq->name.s);
    StringVectorPushBack(seqs,seq->seq.s);
  }
  kseq_destroy(seq);
  gzclose(fp);

  ksc = SafeMalloc (sizeof (Kseq_container));
  ksc->n = StringVectorSize(names);
  ksc->name = SafeMalloc (ksc->n * sizeof(char*));
  ksc->seq = SafeMalloc (ksc->n * sizeof(char*));

  for (i = 0; i < ksc->n; ++i) {
    /* avoid an unnecessary copy by pilfering the StringVector elements & setting them to NULL */
    ksc->name[i] = names->begin[i];
    ksc->seq[i] = seqs->begin[i];
    names->begin[i] = seqs->begin[i] = NULL;
  }

  deleteStringVector(names);
  deleteStringVector(seqs);

  return ksc;
}

void free_kseq_container (Kseq_container* ksc) {
  for (i = 0; i < ksc->n; ++i) {
    SafeFree(ksc->name[i]);
    SafeFree(ksc->seq[i]);
  }
  SafeFree(ksc->name);
  SafeFree(ksc->seq);
  SafeFree(ksc);
}
