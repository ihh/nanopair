#ifndef FAST5EVENTS_INCLUDED
#define FAST5EVENTS_INCLUDED

#include <hdf5.h>
#include <hdf5_hl.h>

/* fast5_event */
typedef struct fast5_event {
  double mean, stdv, length;
  char *model_state, *mp_model_state;
  long move, raw;
} fast5_event;

/* fast5_event_array */
typedef struct fast5_event_array {
  fast5_event* event;
  int n_events;
} fast5_event_array;

fast5_event_array* alloc_fast5_event_array (int model_order, int n_events);
void free_fast5_event_array (fast5_event_array* ev);

fast5_event_array* read_fast5_event_array (const char* filename);

#endif /* FAST5EVENTS_INCLUDED */
