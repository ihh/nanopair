#include <string.h>
#include <stdio.h>
#include "stringmap.h"

/* String* functions */
void* StringNew(const char *a) {
  char *ptr;
  ptr = (char*) SafeMalloc ((strlen(a) + 1) * sizeof(char));
  (void) strcpy (ptr, a);
  return (void*) ptr;
}

void* StringCopy(void* a) {
  return StringNew ((char*) a);
}

void StringDelete(void* a) {
  SafeFree((int*)a);
}

int StringCompare(void* a, void* b) {
  int cmp = strcmp ((char*)a, (char*)b);
  return cmp > 0 ? +1 : (cmp < 0 ? -1 : 0);
}

void StringPrint(FILE* file, void* a) {
  fprintf(file,"%i",*(int*)a);
}

char* StringConcat (const char *a, const char *b) {
  char *s = SafeMalloc ((strlen(a) + strlen(b) + 1) * sizeof(char));
  sprintf (s, "%s%s", a, b);
  return s;
}

