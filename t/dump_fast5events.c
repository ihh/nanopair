/*
    dump_fast5events
    Copyright (C) 2014 Ian Holmes
    Based on fast5tofastq by German Tischler.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <locale.h>

#include "../src/fast5events.h"
#include "../src/util.h"

int main(int argc, char * argv[])
{
	/* file name */
	char * fn = NULL;

	/* check for input file name */
	if ( !(1<argc) )
	{
		fprintf(stderr,"usage: %s <in.fast5>\n", argv[0]);
		return EXIT_FAILURE;
	}

	setlocale(LC_ALL, "C");

	/* set input file name */
	fn = argv[1];

	/* read & destroy event array */
	Fast5_event_array* event_array = read_fast5_event_array (fn, DefaultFast5TickLength);

	printf ("n mean stdv length model_state mp_model_state move raw ticks sumticks_cur sumticks_cur_sq\n");
	for (int n = 0; n < event_array->n_events; ++n) {
	  Fast5_event* ev = event_array->event + n;
	  printf("%d %g %g %g %s %s %ld %ld %g %g %g\n",n,ev->mean,ev->stdv,ev->length,ev->model_state,ev->mp_model_state,ev->move,ev->raw,ev->ticks,ev->sumticks_cur,ev->sumticks_cur_sq);
	}

	delete_fast5_event_array (event_array);
	
	return EXIT_SUCCESS;
}
