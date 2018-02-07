/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * @file   grasp_settings_data_types.h
 * @author David Fuertes
 *
 * Created on 14 de noviembre de 2013, 14:59
 */

#ifndef GRASP_SETTINGS_DATA_TYPES_H
#define	GRASP_SETTINGS_DATA_TYPES_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "yamlsettings/yamlsettings_error.h"
#include "yamlsettings/yamlsettings_data_types.h"
    
#define GRASP_DATA_TYPE_FORTRAN_STRING(L) grasp_data_type_fortran_string_set,grasp_data_type_fortran_string_get,L
int grasp_data_type_fortran_string_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_fortran_string_get(void *mem_pos,  int position, int maxlength);

#define GRASP_DATA_TYPE_FORTRAN_FILE_PATH(L) grasp_data_type_fortran_file_path_set,grasp_data_type_fortran_file_path_get,L
int grasp_data_type_fortran_file_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_fortran_file_path_get(void *mem_pos,  int position, int maxlength);

#define GRASP_DATA_TYPE_FORTRAN_FOLDER_PATH(L) grasp_data_type_fortran_folder_path_set,grasp_data_type_fortran_folder_path_get,L
int grasp_data_type_fortran_folder_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_fortran_folder_path_get(void *mem_pos,  int position, int maxlength);

#define GRASP_DATA_TYPE_STREAM(L) grasp_data_type_stream_set,grasp_data_type_stream_get,L
int grasp_data_type_stream_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_stream_get(void *mem_pos,  int position, int maxlength);

#define GRASP_DATA_TYPE_MINIMIZATION grasp_data_type_minimization_set,grasp_data_type_minimization_get,-1
int grasp_data_type_minimization_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_minimization_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_IPPLANE grasp_data_type_ipplane_set,grasp_data_type_ipplane_get,-1
int grasp_data_type_ipplane_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_ipplane_get(void *mem_pos,int position, int maxlength);
    
#define GRASP_DATA_TYPE_PFITTING grasp_data_type_pfitting_set,grasp_data_type_pfitting_get,-1
int grasp_data_type_pfitting_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_pfitting_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_IFITTING grasp_data_type_ifitting_set,grasp_data_type_ifitting_get,-1
int grasp_data_type_ifitting_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_ifitting_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_BIN grasp_data_type_bin_set,grasp_data_type_bin_get,-1
int grasp_data_type_bin_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_bin_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_INVSING grasp_data_type_invsing_set,grasp_data_type_invsing_get,-1
int grasp_data_type_invsing_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_invsing_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_IMSC grasp_data_type_imsc_set,grasp_data_type_imsc_get,-1
int grasp_data_type_imsc_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_imsc_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_ERROR grasp_data_type_error_set,grasp_data_type_error_get,-1
int grasp_data_type_error_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_error_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_MEASURETYPES grasp_data_type_measuretypes_set,grasp_data_type_measuretypes_get,-1
int grasp_data_type_measuretypes_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_measuretypes_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_IMQ grasp_data_type_imq_set,grasp_data_type_imq_get,-1
int grasp_data_type_imq_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_imq_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_CHARACTERISTIC grasp_data_type_characteristic_set,grasp_data_type_characteristic_get,-1
extern yamlsettings_enumeration_definition grasp_data_type_characteristic;
int grasp_data_type_characteristic_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_characteristic_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_SURFACE grasp_data_type_surface_set,grasp_data_type_surface_get,-1
int grasp_data_type_surface_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_surface_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_MOL_PROF_VERT_TYPE grasp_data_type_molecular_profile_vertical_type_set,grasp_data_type_molecular_profile_vertical_type_get,-1
int grasp_data_type_molecular_profile_vertical_type_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_molecular_profile_vertical_type_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_LAND_OCEAN_FILTER grasp_data_type_land_ocean_filter_set,grasp_data_type_land_ocean_filter_get,-1
int grasp_data_type_land_ocean_filter_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_land_ocean_filter_get(void *mem_pos,int position, int maxlength);

#define GRASP_DATA_TYPE_AER_PROF_VERT_TYPE grasp_data_type_aerosol_profile_vertical_type_set,grasp_data_type_aerosol_profile_vertical_type_get,-1
int grasp_data_type_aerosol_profile_vertical_type_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *grasp_data_type_aerosol_profile_vertical_type_get(void *mem_pos,int position, int maxlength);

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_SETTINGS_DATA_TYPES_H */

