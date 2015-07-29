#include <locale.h>
#include <stdio.h>
#include "../src/fast5events.h"

#define N_EVENTS 2

int main (int argc, char** argv) {
  char *fn;
  Fast5_event_array *events;
  Fast5_event event[N_EVENTS] = { { .mean = .8, .stdv = .5, .length = .5, .model_state = "AAAAA", .mp_model_state = "AAAAA", .move = 0 },
				  { .mean = .08, .stdv = .2, .length = .4, .model_state = "AAAAC", .mp_model_state = "AAAAC", .move = 1 } };

  /* check for output file name */
  if ( argc < 2 )
    {
      fprintf(stderr,"usage: %s <out.fast5>\n", argv[0]);
      return EXIT_FAILURE;
    }

  setlocale(LC_ALL, "C");

  /* set output file name */
  fn = argv[1];

  /* create & write Fast5_event_array */
  events = alloc_fast5_event_array (-1, N_EVENTS);

  SafeFree (events->event);
  events->event = event;

  write_fast5_event_array (events, fn);

  events->event = NULL;
  delete_fast5_event_array (events);

  return EXIT_SUCCESS;
}
