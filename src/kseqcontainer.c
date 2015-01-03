#include <zlib.h>
#include <stdio.h>
#include "kseqcontainer.h"
#include "vector.h"
#include "stringmap.h"

KSEQ_INIT(gzFile, gzread)

const char* dna_alphabet = "ACGT";

int tokenize (char c, const char* alphabet) {
  int tok;
  tok = strchr (alphabet, toupper(c)) - alphabet;
  return tok >= (int) strlen(alphabet) ? -1 : tok;
}

Kseq_container* init_kseq_container (const char* filename, const char* alphabet) {
  gzFile fp;
  kseq_t *seq;
  int l, i, pos, len;

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
  ksc->len = SafeMalloc (ksc->n * sizeof(int));
  ksc->seq = SafeMalloc (ksc->n * sizeof(char*));
  ksc->dsq = SafeMalloc (ksc->n * sizeof(int*));

  for (i = 0; i < ksc->n; ++i) {
    /* avoid an unnecessary copy by pilfering the StringVector elements & setting them to NULL */
    ksc->name[i] = names->begin[i];
    ksc->seq[i] = seqs->begin[i];
    names->begin[i] = seqs->begin[i] = NULL;

    len = strlen (ksc->seq[i]);
    ksc->len[i] = len;
    ksc->dsq[i] = SafeMalloc (len * sizeof(int));
    for (pos = 0; pos < len; ++pos)
      ksc->dsq[i][pos] = tokenize (ksc->seq[i][pos], alphabet);
  }

  deleteStringVector(names);
  deleteStringVector(seqs);

  return ksc;
}

void free_kseq_container (Kseq_container* ksc) {
  int i;
  for (i = 0; i < ksc->n; ++i) {
    SafeFree(ksc->name[i]);
    SafeFree(ksc->seq[i]);
    SafeFree(ksc->dsq[i]);
  }
  SafeFree(ksc->name);
  SafeFree(ksc->len);
  SafeFree(ksc->seq);
  SafeFree(ksc->dsq);
  SafeFree(ksc);
}
