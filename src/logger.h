#ifndef LOGGER_INCLUDED
#define LOGGER_INCLUDED

#include "stringmap.h"

typedef struct Logger {
  int verbosity;
  StringSet *logTags;
} Logger;

Logger* newLogger();
void deleteLogger (Logger*);

int parseLogArgs (int* argcPtr, char*** argvPtr, Logger*);

#define LogAt(V)     (logger != NULL && logger->verbosity >= (V))
#define LogWhen(TAG) (logger != NULL && logger->logTags != NULL && StringSetFind(logger->logTags,TAG) != NULL)
#define LogThis      LogWhen(__FUNCTION__)
#define LogThisAt(V) (LogAt(V) || LogThis)

#endif /* LOGGER_INCLUDED */

