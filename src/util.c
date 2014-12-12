#include "util.h"

/* function defs */
void Warn(char* warning, ...) {
  va_list argptr;
  va_start (argptr, warning);
  vfprintf(stderr,warning,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
}

void Abort(char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  printf("Abort: ");
  vprintf(error,argptr);
  printf("\n");
  va_end (argptr);
  exit(-1);
}

void Assert(int assertion, char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    printf("Assertion Failed: ");
    vprintf(error,argptr);
    printf("\n");
    va_end (argptr);
    exit(-1);
  }
}

void * SafeMalloc(size_t size) {
  void * result;

  if ( (result = malloc(size)) ) { /* assignment intentional */
    return(result);
  } else {
    printf("memory overflow: malloc failed in SafeMalloc.");
    printf("  Exiting Program.\n");
    exit(-1);
  }
  return(0);
}

void *SafeCalloc(size_t count, size_t size) {
  void * result;

  if ( (result = calloc(count,size)) ) { /* assignment intentional */
    return(result);
  } else {
    printf("memory overflow: calloc failed in SafeCalloc.");
    printf("  Exiting Program.\n");
    exit(-1);
  }
  return(0);
}
