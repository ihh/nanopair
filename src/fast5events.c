#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>

#include "fast5events.h"
#include "util.h"

/* path we check for events data */
const char* events_path = "/Analyses/Basecall_2D_000/BaseCalled_template/Events";

/* names of various member fields */
#define FAST5_EVENT_MEAN "mean"
#define FAST5_EVENT_STDV "stdv"
#define FAST5_EVENT_LENGTH "length"
#define FAST5_EVENT_MODEL_STATE "model_state"
#define FAST5_EVENT_MOVE "move"
#define FAST5_EVENT_MP_STATE "mp_state"
#define FAST5_EVENT_RAW_INDEX "raw_index"

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

  /* HACK: allow for variable-length state names */
  if (iter->model_order < 0) {
    char *model_state = *(char**)(elem + iter->model_state_offset);
    int model_state_len = strlen (model_state);
    ev->model_state = SafeMalloc ((model_state_len + 1) * sizeof (char));
    for (int n = 0; n <= model_state_len; ++n)
      ev->model_state[n] = model_state[n];

    char *mp_model_state = *(char**)(elem + iter->mp_model_state_offset);
    int mp_model_state_len = strlen (mp_model_state);
    ev->mp_model_state = SafeMalloc ((mp_model_state_len + 1) * sizeof (char));
    for (int n = 0; n <= mp_model_state_len; ++n)
      ev->mp_model_state[n] = mp_model_state[n];
  } else {  /* fixed-length state names */
    for (int n = 0; n < iter->model_order; ++n) {
      ev->model_state[n] = *((char*) (elem + iter->model_state_offset + n));
      ev->mp_model_state[n] = *((char*) (elem + iter->mp_model_state_offset + n));
    }
    ev->model_state[iter->model_order] = ev->mp_model_state[iter->model_order] = '\0';
  }

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

  /* HACK: if state name strings are variable length in HDF5 file, postpone allocation of state names */
  for (int n = 0; n < n_events; ++n)
    if (model_order < 0) {
      ev->event[n].model_state = NULL;
      ev->event[n].mp_model_state = NULL;
    } else {
      ev->event[n].model_state = SafeMalloc ((model_order + 1) * sizeof (char));
      ev->event[n].mp_model_state = SafeMalloc ((model_order + 1) * sizeof (char));
    }

  ev->tick_length = tick_length;
  return ev;
};

void delete_fast5_event_array (Fast5_event_array* ev) {
  for (int n = 0; n < ev->n_events; ++n) {
    SafeFreeOrNull (ev->event[n].model_state);
    SafeFreeOrNull (ev->event[n].mp_model_state);
  }
  SafeFree (ev->event);
  SafeFree (ev);
};

Fast5_event_array* read_fast5_event_array (const char* filename, double tick_length)
{
  hid_t file_id, strtype;

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
  if (H5LTpath_valid ( file_id, events_path, 1))
    {
      /* get dataset */
      hid_t events_id = H5Oopen(file_id, events_path, H5P_DEFAULT);
      if ( events_id < 0 )
	{
	  fprintf(stderr,"failed to open dataset %s\n",events_path);
	}
      else
	{

	  /* get information about fields in an event */
	  hid_t events_type_id = H5Dget_type( events_id );

	  Fast5_event_array_iterator iter;

	  int mean_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_MEAN );
	  int stdv_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_STDV );
	  int length_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_LENGTH );
	  int model_state_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_MODEL_STATE );
	  int move_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_MOVE );
	  int mp_model_state_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_MP_STATE );
	  int raw_idx = H5Tget_member_index( events_type_id, FAST5_EVENT_RAW_INDEX );

	  iter.mean_offset = H5Tget_member_offset( events_type_id, mean_idx );
	  iter.stdv_offset = H5Tget_member_offset( events_type_id, stdv_idx );
	  iter.length_offset = H5Tget_member_offset( events_type_id, length_idx );
	  iter.model_state_offset = H5Tget_member_offset( events_type_id, model_state_idx );
	  iter.move_offset = H5Tget_member_offset( events_type_id, move_idx );
	  iter.mp_model_state_offset = H5Tget_member_offset( events_type_id, mp_model_state_idx );
	  iter.raw_offset = H5Tget_member_offset( events_type_id, raw_idx );

	  /* HACK: detect variable-length strings, flag using a model_order of -1 */
	  strtype = H5Tget_member_type( events_type_id, model_state_idx );
	  iter.model_order = H5Tis_variable_str(strtype) ? -1 : (int) H5Tget_size (strtype);

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
    fprintf(stderr,"path %s not valid\n",events_path);

  /* close root group */
  H5Gclose(root_id);

  /* close file */
  H5Fclose(file_id);

  return event_array;
}

void write_fast5_event_array (Fast5_event_array* events, const char* filename) {
  hid_t   file, filetype, memtype, strtype, space, dset, loc;
  herr_t  status;
  hsize_t dims[1], strtype_size, dbl_size, int_size, filetype_size;
  char    *obj_name;
  int     i, j;

  if (events->n_events == 0) {

    Warn ("No events; cannot determine model order, not writing %s", filename);

  } else {

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate (filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create variable-length string datatype.
     */
    strtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (strtype, H5T_VARIABLE);

    /*
     * Create the compound datatype for memory.
     */
    memtype = H5Tcreate (H5T_COMPOUND, sizeof (Fast5_event));
    status = H5Tinsert (memtype, FAST5_EVENT_MEAN,
			HOFFSET (Fast5_event, mean),
			H5T_NATIVE_DOUBLE);
    status = H5Tinsert (memtype, FAST5_EVENT_STDV,
			HOFFSET (Fast5_event, stdv),
			H5T_NATIVE_DOUBLE);
    status = H5Tinsert (memtype, FAST5_EVENT_LENGTH,
			HOFFSET (Fast5_event, length),
			H5T_NATIVE_DOUBLE);
    status = H5Tinsert (memtype, FAST5_EVENT_MODEL_STATE,
			HOFFSET (Fast5_event, model_state),
			strtype);
    status = H5Tinsert (memtype, FAST5_EVENT_MOVE,
			HOFFSET (Fast5_event, move),
			H5T_NATIVE_INT);
    status = H5Tinsert (memtype, FAST5_EVENT_MP_STATE,
			HOFFSET (Fast5_event, mp_model_state),
			strtype);
    status = H5Tinsert (memtype, FAST5_EVENT_RAW_INDEX,
			HOFFSET (Fast5_event, raw),
			H5T_NATIVE_INT);

    /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */
    strtype_size = H5Tget_size (strtype);
    dbl_size = H5Tget_size (H5T_IEEE_F64LE);
    int_size = H5Tget_size (H5T_IEEE_F64LE);
    filetype_size = 3*dbl_size + 2*strtype_size + 2*int_size;

    filetype = H5Tcreate (H5T_COMPOUND, filetype_size);
    status = H5Tinsert (filetype, FAST5_EVENT_MEAN, 0, H5T_IEEE_F64LE);
    status = H5Tinsert (filetype, FAST5_EVENT_STDV, dbl_size, H5T_IEEE_F64LE);
    status = H5Tinsert (filetype, FAST5_EVENT_LENGTH, 2*dbl_size, H5T_IEEE_F64LE);
    status = H5Tinsert (filetype, FAST5_EVENT_MODEL_STATE, 3*dbl_size, strtype);
    status = H5Tinsert (filetype, FAST5_EVENT_MOVE, 3*dbl_size + strtype_size, H5T_STD_I64LE);
    status = H5Tinsert (filetype, FAST5_EVENT_MP_STATE, 3*dbl_size + strtype_size + int_size, strtype);
    status = H5Tinsert (filetype, FAST5_EVENT_RAW_INDEX, 3*dbl_size + 2*strtype_size + int_size, H5T_STD_I64LE);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    dims[0] = events->n_events;
    space = H5Screate_simple (1, dims, NULL);

    /* Create the necessary groups to ensure the path to the dataset exists */
    obj_name = SafeMalloc ((strlen(events_path) + 1) * sizeof(char));
    for (loc = file, i = j = 0; events_path[i] != '\0'; ++i) {
      if (events_path[i] == '/') {
	obj_name[j] = '\0';
	if (j > 0)
	  loc = H5Gcreate (loc, obj_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	j = 0;
      } else
	obj_name[j++] = events_path[i];
    }
    obj_name[j] = '\0';

    /*
     * Create the dataset and write the compound data to it.
     */
    dset = H5Dcreate (loc, obj_name, filetype, space, H5P_DEFAULT, H5P_DEFAULT,
		      H5P_DEFAULT);
    status = H5Dwrite (dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, events->event);

    /*
     * Close and release resources.
     */
    status = H5Dclose (dset);
    status = H5Sclose (space);
    status = H5Tclose (filetype);
    status = H5Fclose (file);

    SafeFree (obj_name);
  }
}
