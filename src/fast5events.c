/*
    fast5events
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
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdlib.h>
#include <stdarg.h>
#include <malloc/malloc.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>

#define _XOPEN_SOURCE       /* See feature_test_macros(7) */
#define __USE_XOPEN
#include <time.h>
#include <locale.h>

/* errors, warnings, assertions */
void Abort(char* error, ...);
void Assert(int assertion, char* error, ...);
void Warn(char* warning, ...);

/* generic alloc functions */
void *SafeMalloc(size_t size);
void *SafeCalloc(size_t count, size_t size);
#define SafeFree(PTR) free(PTR)
#define SafeFreeOrNull(PTR) if (PTR) SafeFree(PTR);

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

/* fast5_event_array_iterator
   Used to populate a fast5_event_array */
typedef struct fast5_event_array_iterator {
  fast5_event_array* event_array;
  int event_array_index;
  size_t mean_offset, stdv_offset, length_offset, model_state_offset, move_offset, mp_model_state_offset, raw_offset;
  int model_order;
} fast5_event_array_iterator;

herr_t populate_event_array (void *elem, hid_t type_id, unsigned ndim, 
			     const hsize_t *point, void *operator_data)
{
  fast5_event_array_iterator *iter = (fast5_event_array_iterator*) operator_data;

  fast5_event* ev = iter->event_array->event + (iter->event_array_index++);
  ev->mean = *((double*) (elem + iter->mean_offset));
  ev->stdv = *((double*) (elem + iter->stdv_offset));
  ev->length = *((double*) (elem + iter->length_offset));
  for (int n = 0; n < iter->model_order; ++n) {
    ev->model_state[n] = *((char*) (elem + iter->model_state_offset + n));
    ev->mp_model_state[n] = *((char*) (elem + iter->mp_model_state_offset + n));
  }
  ev->model_state[iter->model_order] = ev->mp_model_state[iter->model_order] = '\0';
  ev->move = *((long*) (elem + iter->move_offset));
  ev->raw = *((long*) (elem + iter->raw_offset));
  fprintf(stderr,"%d %g %g %g %s %s %ld %ld\n",iter->event_array_index-1,ev->mean,ev->stdv,ev->length,ev->model_state,ev->mp_model_state,ev->move,ev->raw);
  return 0;
}



/* function defs */
void Warn(char* warning, ...) {
  va_list argptr;
  va_start (argptr, warning);
  vfprintf(stderr,warning,argptr);
  fprintf(stderr,"\n");
  va_end (argptr);
}

void Abort(char* error, ...) {
  va_list argptr;
  va_start (argptr, error);
  printf("Abort: ");
  vprintf(error,argptr);
  printf("\n");
  va_end (argptr);
  exit(-1);
}

void Assert(int assertion, char* error, ...) {
  va_list argptr;
  if(!assertion) {
    va_start (argptr, error);
    printf("Assertion Failed: ");
    vprintf(error,argptr);
    printf("\n");
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

fast5_event_array* alloc_fast5_event_array (int model_order, int n_events) {
  fast5_event_array* ev = SafeMalloc (sizeof (fast5_event_array));
  ev->n_events = n_events;
  ev->event = SafeMalloc (n_events * sizeof (fast5_event));
  for (int n = 0; n < n_events; ++n) {
    ev->event[n].model_state = SafeMalloc ((model_order + 1) * sizeof (char));
    ev->event[n].mp_model_state = SafeMalloc ((model_order + 1) * sizeof (char));
  }
  return ev;
};

void free_fast5_event_array (fast5_event_array* ev) {
  for (int n = 0; n < ev->n_events; ++n) {
    SafeFree (ev->event[n].model_state);
    SafeFree (ev->event[n].mp_model_state);
  }
  SafeFree (ev->event);
  SafeFree (ev);
};

int main(int argc, char * argv[])
{
	/* file name */
	char * fn = NULL;
	/* path we check for events data */
	const char* path = "/Analyses/Basecall_2D_000/BaseCalled_template/Events";
        /* file handle */
        hid_t file_id;
        /* op status */
        herr_t status;

	/* event_array */
	fast5_event_array* event_array = NULL;

	/* check for input file name */
	if ( !(1<argc) )
	{
		fprintf(stderr,"usage: %s <in.fast5>\n", argv[0]);
		return EXIT_FAILURE;
	}

	setlocale(LC_ALL, "C");

	/* set input file name */
	fn = argv[1];
	
	/* open file with default properties */
	file_id = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
        
        /* check if opening file was succeful */
        if ( file_id < 0 )
        {
        	fprintf(stderr,"Failed to open/init input file %s\n", fn);
        	return EXIT_FAILURE;
        }

	/* get root group */
	hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
	if (root_id < 0)
	  {
	    fprintf(stderr,"failed to open root group\n");
	    return EXIT_FAILURE;
	  }

	/* see if path exists */
	if (H5LTpath_valid ( file_id, path, 1))
	  {
	    fprintf(stderr,"path %s is valid\n",path);

	    /* get dataset */
	    hid_t events_id = H5Oopen(file_id, path, H5P_DEFAULT);
	    if ( events_id < 0 )
	      {
		fprintf(stderr,"failed to open dataset %s\n",path);
	      }
	    else
	      {

		/* get information about fields in an event */
		hid_t events_type_id = H5Dget_type( events_id );

		fast5_event_array_iterator iter;

		int mean_idx = H5Tget_member_index( events_type_id, "mean" );
		int stdv_idx = H5Tget_member_index( events_type_id, "stdv" );
		int length_idx = H5Tget_member_index( events_type_id, "length" );
		int model_state_idx = H5Tget_member_index( events_type_id, "model_state" );
		int move_idx = H5Tget_member_index( events_type_id, "move" );
		int mp_model_state_idx = H5Tget_member_index( events_type_id, "mp_state" );
		int raw_idx = H5Tget_member_index( events_type_id, "raw_index" );

		fprintf(stderr,"member indices: %d %d %d %d %d %d %d\n", mean_idx, stdv_idx, length_idx, model_state_idx, move_idx, mp_model_state_idx, raw_idx);

		iter.mean_offset = H5Tget_member_offset( events_type_id, mean_idx );
		iter.stdv_offset = H5Tget_member_offset( events_type_id, stdv_idx );
		iter.length_offset = H5Tget_member_offset( events_type_id, length_idx );
		iter.model_state_offset = H5Tget_member_offset( events_type_id, model_state_idx );
		iter.move_offset = H5Tget_member_offset( events_type_id, move_idx );
		iter.mp_model_state_offset = H5Tget_member_offset( events_type_id, mp_model_state_idx );
		iter.raw_offset = H5Tget_member_offset( events_type_id, raw_idx );

		fprintf(stderr,"member offsets: %lu %lu %lu %lu %lu %lu %lu\n", iter.mean_offset, iter.stdv_offset, iter.length_offset, iter.model_state_offset, iter.move_offset, iter.mp_model_state_offset, iter.raw_offset);

		iter.model_order =  (int) H5Tget_size (H5Tget_member_type( events_type_id, model_state_idx ));

		fprintf(stderr,"model_order is %d\n",iter.model_order);


		/* get dimensions */
		hid_t events_space_id = H5Dget_space( events_id );
		hssize_t events_npoints = H5Sget_simple_extent_npoints( events_space_id );
		size_t event_size = H5Tget_size( events_type_id );

		fprintf(stderr,"dataset has %lld points, event takes %lu bytes\n",events_npoints,event_size);


		/*
		 * Create variable-length string datatype.
		 */
		hid_t strtype = H5Tcopy (H5T_C_S1);
		H5Tset_size (strtype, H5T_VARIABLE);
		    
		/*
		 * Create the compound datatype for memory.
		hid_t fast5_event_type_id = H5Tcreate (H5T_COMPOUND, sizeof (fast5_event));
		H5Tinsert (fast5_event_type_id, "mean", HOFFSET (fast5_event, mean), H5T_NATIVE_DOUBLE);
		H5Tinsert (fast5_event_type_id, "stdv", HOFFSET (fast5_event, stdv), H5T_NATIVE_DOUBLE);
		H5Tinsert (fast5_event_type_id, "length", HOFFSET (fast5_event, length), H5T_NATIVE_DOUBLE);
		H5Tinsert (fast5_event_type_id, "model_state", HOFFSET (fast5_event, model_state), strtype);
		H5Tinsert (fast5_event_type_id, "mp_model_state", HOFFSET (fast5_event, mp_model_state), strtype);
		H5Tinsert (fast5_event_type_id, "move", HOFFSET (fast5_event, move), H5T_NATIVE_INT);
		H5Tinsert (fast5_event_type_id, "raw_index", HOFFSET (fast5_event, raw), H5T_NATIVE_INT);
		 */

		/* read into memory buffer */
		void *buf = SafeMalloc (events_npoints * event_size);
		H5Dread( events_id, events_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );

		/* convert */
		event_array = alloc_fast5_event_array (iter.model_order, events_npoints);
		iter.event_array = event_array;
		iter.event_array_index = 0;
		H5Diterate( buf, events_type_id, events_space_id, populate_event_array, &iter );
		      
		/* close */
		H5Sclose(events_space_id);
		H5Tclose(events_type_id);
		H5Dclose(events_id);
	      }
	  }
	else
	  fprintf(stderr,"path %s not valid\n",path);

	/* close root */
	H5Gclose(root_id);

	/* close file */
	status = H5Fclose(file_id);
	
	if ( status < 0 )
	{
		fprintf(stderr,"Failed to close HDF file %s\n", fn);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
