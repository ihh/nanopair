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
  int move, raw;
} fast5_event;

/* fast5_event_array */
typedef struct fast5_event_array {
  fast5_event* event;
  int n_events;
} fast5_event_array;

fast5_event_array* alloc_fast5_event_array (int model_order, int n_events);
void free_fast5_event_array (fast5_event_array* ev);


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
	/* set of hdf paths we check for events data */
        char const * events_paths[] =
        {
	        "/Analyses/Basecall_2D_000/BaseCalled_2D/Events",
	        "/Analyses/Basecall_2D_000/BaseCalled_template/Events",
	        "/Analyses/Basecall_2D_000/BaseCalled_complement/Events",
	        0
        };
	/* iterator for events_paths */
	char const ** evp = NULL;
	/* still searching for events */
        int running = 1;
        /* file handle */
        hid_t file_id;
        /* op status */
        herr_t status;

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

        /* try list of paths to events data */
        for ( evp = &events_paths[0]; running && *evp; ++evp )
        {
        	/* events path */
        	char const * path = *evp;
		
		/* get root group */
	        hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
		if (root_id < 0)
	        {
		  fprintf(stderr,"failed to open root group\n");
		  break;
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
			break;
		      }

		    H5O_info_t events_info;
		    herr_t status = H5Oget_info( events_id, &events_info );
		    fprintf(stderr,"events object has type %d\n",events_info.type);

		    hid_t events_type_id = H5Dget_type( events_id );

		    int strand_idx = H5Tget_member_index( events_type_id, "strand" );
		    int mean_idx = H5Tget_member_index( events_type_id, "mean" );
		    int stdv_idx = H5Tget_member_index( events_type_id, "stdv" );
		    int length_idx = H5Tget_member_index( events_type_id, "length" );
		    int model_state_idx = H5Tget_member_index( events_type_id, "model_state" );
		    int move_idx = H5Tget_member_index( events_type_id, "move" );
		    int mp_model_state_idx = H5Tget_member_index( events_type_id, "mp_model_state" );
		    int raw_idx = H5Tget_member_index( events_type_id, "raw" );

		    int model_order =  (int) H5Tget_size (H5Tget_member_type( events_type_id, model_state_idx ));

		    fprintf(stderr,"model_order is %d\n",model_order);

		    hid_t events_space_id = H5Dget_space( events_id );
		    hssize_t events_npoints = H5Sget_simple_extent_npoints( events_space_id );

		    fprintf(stderr,"dataset has %d points\n",events_npoints);

		      
		    /* close */
		    H5Tclose(events_type_id);
		    H5Dclose(events_id);
		    H5Gclose(root_id);

		    /* stop looking */
		    running = 0;
		  }
		else
		  fprintf(stderr,"path %s not valid\n",path);
	}

	/* close file */
	status = H5Fclose(file_id);
	
	if ( status < 0 )
	{
		fprintf(stderr,"Failed to close HDF file %s\n", fn);
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
