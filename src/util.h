#ifndef UTIL_INCLUDED
#define UTIL_INCLUDED

#include <stdlib.h>

/* function pointer typedefs for generic containers */
typedef int (*CompareFunction) (void*, void*);
typedef void (*DestroyFunction) (void*);
typedef void* (*CopyFunction) (void*);
typedef void (*PrintFunction) (void*);

/* null functions for generic containers */
void NullDestroyFunction(void*);  /* does nothing */
void* NullCopyFunction(void*);  /* returns the supplied parameter without doing anything */
void NullPrintFunction(void*);  /* does nothing */

/* abort functions for when null functions should never be called */
void AbortDestroyFunction(void*);
void* AbortCopyFunction(void*);

/* typedef and container functions for 64-bit signed int's.
   (It's tempting to think that rather than allocating space,
   one could just use the (void*) pointer to store the int value;
   however, this risks platform-specific errors/warnings,
   due to differences in bytesize/signedness between void* and int.)
*/
typedef signed long long int Int64;

void* IntNew(Int64 a);
void* IntCopy(void* a);
void IntDelete(void* a);
int IntCompare(void* a, void* b);
void IntPrint(void* a);

/* Container functions for double's. */
void* DoubleNew(double a);
void* DoubleCopy(void* a);
void DoubleDelete(void* a);
int DoubleCompare(void* a, void* b);
void DoublePrint(void* a);

/* DUMP */
#undef DUMP
#define DUMP(x, fmt) printf("%s:%u: %s=" fmt, __FILE__, __LINE__, #x, x)

/* MIN, MAX, ABS */
#undef MIN
#undef MAX
#undef ABS
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define ABS(X)   ((X) >= 0 ? (X) : -(X))

#define MIN2(X,Y,Z) MIN(MIN(X,Y),Z)
#define MAX2(X,Y,Z) MAX(MAX(X,Y),Z)

#define MIN3(W,X,Y,Z) MIN(MIN(W,X),MIN(Y,Z))
#define MAX3(W,X,Y,Z) MAX(MAX(W,X),MAX(Y,Z))

long double max_func (long double x, long double y);
long double min_func (long double x, long double y);

/* errors, warnings, assertions */
void Abort(const char* error, ...);
void Assert(int assertion, const char* error, ...);
void Warn(const char* warning, ...);

/* alloc functions */
void *SafeMalloc(size_t size);
void *SafeCalloc(size_t count, size_t size);
#define SafeFree(PTR) free(PTR)
#define SafeFreeOrNull(PTR) if (PTR) SafeFree(PTR);

/* stupid quoting bollix */
#define QUOTEME_(x) #x
#define QUOTEME(x) QUOTEME_(x)

/* functions to convert decimal/hexadecimal strings to ints.
 */
int decToSignedInt (const char *);
unsigned int hexToUnsignedInt (const char *);

/* function to read file into a string */
char* readFileAsString (const char *);  /* must be freed by caller */

#endif /* UTIL_INCLUDED */
