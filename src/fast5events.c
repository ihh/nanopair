#include "fast5events.h"
#include "util.h"

/* Fast5_event_array_iterator
   Used to populate a Fast5_event_array */
typedef struct Fast5_event_array_iterator {
  Fast5_event_array* event_array;
  int event_array_index;
  size_t mean_offset, stdv_offset, length_offset, model_state_offset, move_offset, mp_model_state_offset, raw_offset;
  int model_order;
} Fast5_event_array_iterator;

herr_t populate_event_array (void *elem, hid_t type_id, unsigned ndim, 
			     const hsize_t *point, void *operator_data)
{
  Fast5_event_array_iterator *iter = (Fast5_event_array_iterator*) operator_data;

  Fast5_event* ev = iter->event_array->event + (iter->event_array_index++);
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
  ev->ticks = ev->length / iter->event_array->tick_length;
  ev->sumticks_cur = ev->ticks * ev->mean;
  ev->sumticks_cur_sq = ev->ticks * (ev->stdv * ev->stdv + ev->mean * ev->mean);
  return 0;
}

Fast5_event_array* alloc_fast5_event_array (int model_order, int n_events, double tick_length) {
  Fast5_event_array* ev = SafeMalloc (sizeof (Fast5_event_array));
  ev->n_events = n_events;
  ev->event = SafeMalloc (n_events * sizeof (Fast5_event));
  for (int n = 0; n < n_events; ++n) {
    ev->event[n].model_state = SafeMalloc ((model_order + 1) * sizeof (char));
    ev->event[n].mp_model_state = SafeMalloc ((model_order + 1) * sizeof (char));
  }
  ev->tick_length = tick_length;
  return ev;
};

void delete_fast5_event_array (Fast5_event_array* ev) {
  for (int n = 0; n < ev->n_events; ++n) {
    SafeFree (ev->event[n].model_state);
    SafeFree (ev->event[n].mp_model_state);
  }
  SafeFree (ev->event);
  SafeFree (ev);
};

Fast5_event_array* read_fast5_event_array (const char* filename, double tick_length)
{
  /* path we check for events data */
  const char* path = "/Analyses/Basecall_2D_000/BaseCalled_template/Events";

  /* file handle */
  hid_t file_id;

  /* event_array */
  Fast5_event_array* event_array = NULL;

  /* open file with default properties */
  file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
        
  /* check if opening file was succeful */
  if ( file_id < 0 )
    {
      fprintf(stderr,"Failed to open/init input file %s\n", filename);
      return NULL;
    }

  /* get root group */
  hid_t root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
  if (root_id < 0)
    {
      fprintf(stderr,"failed to open root group\n");
      return NULL;
    }

  /* see if path exists */
  if (H5LTpath_valid ( file_id, path, 1))
    {
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

	  Fast5_event_array_iterator iter;

	  int mean_idx = H5Tget_member_index( events_type_id, "mean" );
	  int stdv_idx = H5Tget_member_index( events_type_id, "stdv" );
	  int length_idx = H5Tget_member_index( events_type_id, "length" );
	  int model_state_idx = H5Tget_member_index( events_type_id, "model_state" );
	  int move_idx = H5Tget_member_index( events_type_id, "move" );
	  int mp_model_state_idx = H5Tget_member_index( events_type_id, "mp_state" );
	  int raw_idx = H5Tget_member_index( events_type_id, "raw_index" );

	  iter.mean_offset = H5Tget_member_offset( events_type_id, mean_idx );
	  iter.stdv_offset = H5Tget_member_offset( events_type_id, stdv_idx );
	  iter.length_offset = H5Tget_member_offset( events_type_id, length_idx );
	  iter.model_state_offset = H5Tget_member_offset( events_type_id, model_state_idx );
	  iter.move_offset = H5Tget_member_offset( events_type_id, move_idx );
	  iter.mp_model_state_offset = H5Tget_member_offset( events_type_id, mp_model_state_idx );
	  iter.raw_offset = H5Tget_member_offset( events_type_id, raw_idx );

	  iter.model_order =  (int) H5Tget_size (H5Tget_member_type( events_type_id, model_state_idx ));

	  /* get dimensions */
	  hid_t events_space_id = H5Dget_space( events_id );
	  hssize_t events_npoints = H5Sget_simple_extent_npoints( events_space_id );
	  size_t event_size = H5Tget_size( events_type_id );

	  /* read into memory buffer */
	  void *buf = SafeMalloc (events_npoints * event_size);
	  H5Dread( events_id, events_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf );

	  /* convert */
	  event_array = alloc_fast5_event_array (iter.model_order, events_npoints, tick_length);
	  iter.event_array = event_array;
	  iter.event_array_index = 0;
	  H5Diterate( buf, events_type_id, events_space_id, populate_event_array, &iter );

	  /* free buffer */
	  SafeFree (buf);

	  /* close objects */
	  H5Sclose(events_space_id);
	  H5Tclose(events_type_id);
	  H5Dclose(events_id);
	}
    }
  else
    fprintf(stderr,"path %s not valid\n",path);

  /* close root group */
  H5Gclose(root_id);

  /* close file */
  H5Fclose(file_id);

  return event_array;
}
