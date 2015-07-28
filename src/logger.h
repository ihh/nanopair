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

#define LogVerb(LOGGER,V) ((LOGGER) != NULL && (LOGGER)->verbosity >= (V))
#define LogTag(LOGGER,TAG) ((LOGGER) != NULL && (LOGGER)->logTags != NULL && StringSetFind((LOGGER)->logTags,TAG) != NULL)
#define LogFuncTag(LOGGER) LogTag(LOGGER,__FUNCTION__)
#define LogFunc(LOGGER,V) (LogVerb(LOGGER,V) || LogFuncTag(LOGGER))

#endif /* LOGGER_INCLUDED */

