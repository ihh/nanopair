#include <stdio.h>
#include <locale.h>

#include "../src/fast5events.h"
#include "../src/util.h"

int main(int argc, char * argv[])
{
  /* file name */
  char *fnIn = NULL, *fnOut = NULL;

  /* check for input & output file names */
  if ( argc != 3 )
    {
      fprintf(stderr,"usage: %s <in.fast5> <out.fast5>\n", argv[0]);
      return EXIT_FAILURE;
    }

  setlocale(LC_ALL, "C");

  /* set input & output file names */
  fnIn = argv[1];
  fnOut = argv[2];

  /* read, write & destroy event array */
  Fast5_event_array* event_array = read_fast5_event_array (fnIn);
  write_fast5_event_array (event_array, fnOut);
  delete_fast5_event_array (event_array);
	
  return EXIT_SUCCESS;
}
