#include <string.h>
#include "logger.h"

Logger* newLogger() {
  Logger* logger = SafeMalloc (sizeof (Logger));
  logger->verbosity = 0;
  logger->logTags = NULL;
  return logger;
}

void deleteLogger (Logger* logger) {
  if (logger->logTags)
    deleteStringMap (logger->logTags);
  SafeFree (logger);
}

int parseLogArgs (int* argcPtr, char*** argvPtr, Logger* logger) {
  if (*argcPtr > 0) {
    const char* arg = **argvPtr;
    if (strcmp (arg, "-log") == 0) {
      Assert (*argcPtr > 1, "%s must have an argument", **argvPtr);
      const char* tag = (*argvPtr)[1];
      if (logger->logTags == NULL)
	logger->logTags = newStringSet();
      if (StringSetFind(logger->logTags,tag) == NULL)
	StringSetInsert (logger->logTags, tag);
      *argvPtr += 2;
      *argcPtr -= 2;
      return 1;
    } else if (strcmp (arg, "-verbose") == 0) {
      logger->verbosity = MAX (logger->verbosity, 1);
      *argvPtr += 1;
      *argcPtr -= 1;
      return 1;
    } else if (arg[0] == '-' && arg[1] == 'v') {
      int all_v = 1, v;
      for (v = 1; arg[v+1] != '\0'; ++v) {
	if (arg[v+1] != 'v') {
	  all_v = 0;
	  break;
	}
      }
      if (all_v) {
	logger->verbosity = MAX (logger->verbosity, v);
	*argvPtr += 1;
	*argcPtr -= 1;
	return 1;
      } else {
	v = atoi (arg + 2);
	if (v >= 1) {
	  logger->verbosity = MAX (logger->verbosity, v);
	  *argvPtr += 1;
	  *argcPtr -= 1;
	  return 1;
	}
      }
    }
  }
  return 0;
}
