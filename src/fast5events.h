#ifndef FAST5EVENTS_INCLUDED
#define FAST5EVENTS_INCLUDED

#include <hdf5.h>
#include <hdf5_hl.h>

/* fast5_event */
typedef struct Fast5_event {
  double mean, stdv, length;
  char *model_state, *mp_model_state;
  long move, raw;
} Fast5_event;

/* fast5_event_array */
typedef struct Fast5_event_array {
  Fast5_event* event;
  int n_events;
} Fast5_event_array;

Fast5_event_array* alloc_fast5_event_array (int model_order, int n_events);
void delete_fast5_event_array (Fast5_event_array* ev);

Fast5_event_array* read_fast5_event_array (const char* filename);

#endif /* FAST5EVENTS_INCLUDED */
