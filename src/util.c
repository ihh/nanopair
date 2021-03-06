#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "util.h"

void Warn(const char* warning, ...) {
  va_list argptr;
  va_start (argptr, warning);
  vfprintf(stderr,warning,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
}

void Abort(const char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  fprintf(stderr,"Abort: ");
  vfprintf(stderr,error,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
  exit(-1);
}

void Assert(int assertion, const char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    fprintf(stderr,"Assertion Failed: ");
    vfprintf(stderr,error,argptr);
    fprintf(stderr,"\n");
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

/*  NullDestroyFunction & NullPrintFunction do nothing; they are included so that they can be passed
    as a function to newRBTree, etc, when no other suitable function has been defined */

void* NullCopyFunction(void * item) { return item; }
void NullDestroyFunction(void * junk) { ; }
void NullPrintFunction(FILE*file,void * junk) { ; }

/* Abort null functions */
void AbortDestroyFunction(void* junk) { Abort ("Unimplemented destroy function called"); }
void* AbortCopyFunction(void* junk) { Abort ("Unimplemented copy function called"); return NULL; }

/* Int functions */
void* IntNew(int a) {
  int *ptr;
  ptr = (int*) SafeMalloc (sizeof (int));
  *ptr = a;
  return (void*) ptr;
}

void* IntCopy(void* a) {
  return IntNew (*(int*)a);
}

void IntDelete(void* a) {
  SafeFree((int*)a);
}

int IntCompare(void* a, void* b) {
  if( *(int*)a > *(int*)b) return 1;
  if( *(int*)a < *(int*)b) return -1;
  return 0;
}

void IntPrint(FILE* file, void* a) {
  fprintf(file,"%d",*(int*)a);
}

/* Double functions */
void* DoubleNew(double a) {
  double *ptr;
  ptr = (double*) SafeMalloc (sizeof (double));
  *ptr = a;
  return (void*) ptr;
}

void* DoubleCopy(void* a) {
  return DoubleNew (*(double*)a);
}

void DoubleDelete(void* a) {
  SafeFree((double*)a);
}

int DoubleCompare(void* a, void* b) {
  if( *(double*)a > *(double*)b) return(1);
  if( *(double*)a < *(double*)b) return(-1);
  return(0);
}

void DoublePrint(FILE* file, void* a) {
  fprintf(file,"%g",*(double*)a);
}

/* char* to unsigned int conversions */
int decToSignedInt( const char *ca ) {
  int ig;
  char c;
  int sign;
  ig = 0;
  /* test for prefixing white space */
  while (*ca == ' ' || *ca == '\t' ) 
    ca++;
  /* Check sign entered or no */
  sign = 1;
  if ( *ca == '-' )
    sign = -1;
  /* convert string to int */
  while ((c = tolower(*ca++)) != '\0')
    if (c >= '0' && c <= '9')
      ig = ig * 10LL + (int) (c - '0');
  return ig * (int) sign;
}

unsigned int hexToUnsignedInt( const char *ca ) {
  unsigned int ig;
  char c;
  ig = 0;
  /* test for prefixing white space */
  while (*ca == ' ' || *ca == '\t' ) 
    ca++;
  /* convert string to int */
  while ((c = tolower(*ca++)) != '\0')
    if (c >= '0' && c <= '9')
      ig = ig * 16LL + (unsigned int) (c - '0');
    else if (c >= 'a' && c <= 'f')
      ig = ig * 16LL + 10LL + (unsigned int) (c - 'a');
  return ig;
}


long double max_func (long double x, long double y) {
  return MAX (x, y);
}

long double min_func (long double x, long double y) {
  return MIN (x, y);
}

char* readFileAsString (const char *filename) {
  char * buffer = 0;
  long length;
  FILE * f = fopen (filename, "rb");

  if (f) {
    fseek (f, 0, SEEK_END);
    length = ftell (f);
    fseek (f, 0, SEEK_SET);
    buffer = malloc (length);
    if (buffer)
      fread (buffer, 1, length, f);
    fclose (f);
  }

  return buffer;
}

clock_t progress_startTime;
double progress_lastElapsedSeconds, progress_reportInterval;
char* progress_desc = NULL;
void init_progress (const char* desc, ...) {
  progress_startTime = clock();
  progress_lastElapsedSeconds = 0;
  progress_reportInterval = 2;

  time_t rawtime;
  struct tm * timeinfo;

  time (&rawtime);
  timeinfo = localtime (&rawtime);
  
  SafeFreeOrNull (progress_desc);

  va_list argptr;
  va_start (argptr, desc);
  vasprintf (&progress_desc, desc, argptr);
  va_end (argptr);
  fprintf (stderr, "%s: started at %s", progress_desc, asctime(timeinfo));
}

void log_progress (double completedFraction, const char* desc, ...) {
  va_list argptr;
  const clock_t currentTime = clock();
  const double elapsedSeconds = ((double) (currentTime - progress_startTime)) / CLOCKS_PER_SEC;
  const double estimatedTotalSeconds = elapsedSeconds / completedFraction;
  if (elapsedSeconds > progress_lastElapsedSeconds + progress_reportInterval) {
    const double estimatedSecondsLeft = estimatedTotalSeconds - elapsedSeconds;
    const double estimatedMinutesLeft = estimatedSecondsLeft / 60;
    const double estimatedHoursLeft = estimatedMinutesLeft / 60;
    const double estimatedDaysLeft = estimatedHoursLeft / 24;
    fprintf (stderr, "%s: ", progress_desc);
    va_start (argptr, desc);
    vfprintf (stderr, desc, argptr);
    va_end (argptr);
    fprintf (stderr, ". Estimated time left: ");
    if (estimatedDaysLeft > 2)
      fprintf (stderr, "%g days", estimatedDaysLeft);
    else if (estimatedHoursLeft > 2)
      fprintf (stderr, "%g hrs", estimatedHoursLeft);
    else if (estimatedMinutesLeft > 2)
      fprintf (stderr, "%g mins", estimatedMinutesLeft);
    else
      fprintf (stderr, "%g secs", estimatedSecondsLeft);
    fprintf (stderr, " (%g%%)\n", 100*completedFraction);
    progress_lastElapsedSeconds = elapsedSeconds;
    progress_reportInterval = MIN (10, 2*progress_reportInterval);
  }
}
