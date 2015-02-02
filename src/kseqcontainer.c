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

Kseq_container* init_kseq_container (const char* filename) {
  gzFile fp;
  kseq_t *seq;
  int l, i, c, pos, len;

  StringVector *names, *seqs;
  Kseq_container* ksc;

  fp = gzopen(filename, "r");
  if (fp == Z_NULL)
    return NULL;

  names = newStringVector();
  seqs = newStringVector();

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

  for (c = 0; c < 256; ++c)
    ksc->freq[c] = 0;

  for (i = 0; i < ksc->n; ++i) {
    /* avoid an unnecessary copy by pilfering the StringVector elements & setting them to NULL */
    ksc->name[i] = names->begin[i];
    ksc->seq[i] = seqs->begin[i];
    names->begin[i] = seqs->begin[i] = NULL;

    len = strlen (ksc->seq[i]);
    ksc->len[i] = len;

    for (pos = 0; pos < len; ++pos)
      ++ksc->freq[(int) ksc->seq[i][pos]];
  }

  deleteStringVector(names);
  deleteStringVector(seqs);

  return ksc;
}

int validate_kseq_container (Kseq_container* ksc, const char *alphabet, FILE *err) {
  int c, count, badchars;

  for (badchars = c = 0; c < 256; ++c)
    if (ksc->freq[c] && tokenize (c, alphabet) < 0)
      badchars += ksc->freq[c];

  if (badchars && err != NULL) {
    fprintf (err, "Warning: sequence file contains bad characters");
    for (count = c = 0; c < 256; ++c)
      if (ksc->freq[c] && tokenize (c, alphabet) < 0)
	fprintf (err, "%s %c (%d)", count++ ? "" : ",", c, ksc->freq[c]);
    fprintf (err, "\n");
  }

  return badchars;
}

void free_kseq_container (Kseq_container* ksc) {
  int i;
  for (i = 0; i < ksc->n; ++i) {
    SafeFree(ksc->name[i]);
    SafeFree(ksc->seq[i]);
  }
  SafeFree(ksc->name);
  SafeFree(ksc->len);
  SafeFree(ksc->seq);
  SafeFree(ksc);
}

char* new_revcomp_seq (char* dna_seq, int len) {
  int i, tok;
  char *rev;
  rev = SafeMalloc ((len + 1) * sizeof(char));
  for (i = 0; i < len; ++i) {
    tok = tokenize (dna_seq[i], dna_alphabet);
    rev[len - 1 - i] = tok < 0 ? 'N' : dna_alphabet[4 - tok];
  }
  return rev;
}
