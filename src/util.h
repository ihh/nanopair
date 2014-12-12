#ifndef FAST5UTIL_INCLUDED
#define FAST5UTIL_INCLUDED

#include <stdlib.h>
#include <stdarg.h>
#include <malloc/malloc.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <locale.h>

/* errors, warnings, assertions */
void Abort(char* error, ...);
void Assert(int assertion, char* error, ...);
void Warn(char* warning, ...);

/* generic alloc functions */
void *SafeMalloc(size_t size);
void *SafeCalloc(size_t count, size_t size);
#define SafeFree(PTR) free(PTR)
#define SafeFreeOrNull(PTR) if (PTR) SafeFree(PTR);


#endif /* FAST5UTIL_INCLUDED */
