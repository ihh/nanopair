#ifndef FAST5EVENTS_INCLUDED
#define FAST5EVENTS_INCLUDED

#include <hdf5.h>
#include <hdf5_hl.h>

#include "vector.h"

/* Default tick length.
   This is a crude hack to allow us to treat variable-length segments as a sequence of discrete current samples ("ticks").
   In practise it should be the lowest common denominator of the event lengths.
   The value of 0.0002 is based on the usual sample rate of 5000.
   If working with raw data files, it should be possible to read the sample rate from the metadata.
 */
#define DefaultFast5TickLength 0.0002

/* fast5_event */
typedef struct Fast5_event {
  double mean, stdv, length;
  char *model_state, *mp_model_state;
  long move, raw;
  /* { mean, stdv, length } converted to zeroth, first & second moments of a series of "ticks" */
  double ticks,       /*           ticks = length / tick_length */
    sumticks_cur,     /*    sumticks_cur = ticks * mean */
    sumticks_cur_sq;  /* sumticks_cur_sq = ticks * (stdv^2 + mean^2) */
} Fast5_event;

/* fast5_event_array */
typedef struct Fast5_event_array {
  Fast5_event* event;
  int n_events;
  double tick_length;
} Fast5_event_array;

Fast5_event_array* alloc_fast5_event_array (int model_order, int n_events, double tick_length);
void delete_fast5_event_array (Fast5_event_array* ev);

Fast5_event_array* read_fast5_event_array (const char* filename, double tick_length);

#endif /* FAST5EVENTS_INCLUDED */
