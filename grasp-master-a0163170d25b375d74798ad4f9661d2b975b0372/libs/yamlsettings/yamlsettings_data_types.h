/**
 * @file   yamlsettings_data_types.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Low-level data type definition for yamlsettings_dictionary
 *
 * This file contains basic (low-level) assignation functions
 */
#ifndef _YAMLSETTINGS_DATA_TYPES_H
#define _YAMLSETTINGS_DATA_TYPES_H

#ifdef	__cplusplus
extern "C" {
#endif


#include <stdbool.h>

#include "yamlsettings_error.h"

typedef int (*yamlsettings_data_type_set)(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);    
typedef char *(*yamlsettings_data_type_get)(void *mem_pos,int position, int nelements);

// Extract data in a memory position according to data type. Position let manipulate arrays. 
// 0 if not is an array, else a number of determined position. 
// Return 0 if the process is ok, otherwise number lower than 0. (-10 means that the number cant be allocated in the memory position, there is not memory position)
int yamlsettings_data_type_assign(yamlsettings_data_type_set function_type_set, yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlegnth, int file_index, const char *settings_file);

// Print data stored like a yamlsettings type. Transform data in screen text.
char *yamlsettings_data_type_get_string(yamlsettings_data_type_get function_type_get,void *mem_pos,int position, int nelements);

/* DATA TYPES */

#define YS_DATA_TYPE_BOOLEAN yamlsettings_data_type_boolean_set,yamlsettings_data_type_boolean_get,-1
int yamlsettings_data_type_boolean_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_boolean_get(void *mem_pos,  int position, int maxlength);

#define YS_DATA_TYPE_INT yamlsettings_data_type_int_set,yamlsettings_data_type_int_get,-1
int yamlsettings_data_type_int_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_int_get(void *mem_pos,  int position, int maxlength);

#define YS_DATA_TYPE_FLOAT yamlsettings_data_type_float_set,yamlsettings_data_type_float_get,-1
int yamlsettings_data_type_float_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_float_get(void *mem_pos,  int position, int maxlength);

#define YS_DATA_TYPE_STRING(L) yamlsettings_data_type_string_set,yamlsettings_data_type_string_get,L
int yamlsettings_data_type_string_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_string_get(void *mem_pos,  int position, int maxlength);

#define YS_DATA_TYPE_FILE_PATH(L) yamlsettings_data_type_file_path_set,yamlsettings_data_type_file_path_get,L
int yamlsettings_data_type_file_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_file_path_get(void *mem_pos,  int position, int maxlength);

// This data type it same than file path but add a slash (/) at the end of path if it doesnt exist
#define YS_DATA_TYPE_FOLDER_PATH(L) yamlsettings_data_type_folder_path_set,yamlsettings_data_type_folder_path_get,L
int yamlsettings_data_type_folder_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
char *yamlsettings_data_type_folder_path_get(void *mem_pos,  int position, int maxlength);

/* FOR ENUMERATION SYSTEM*/

// Maximum length when a type is transformed to string
#define YAMLSETTINGS_VAR_VALUE_MAX_LENGTH 60

// Definition of a greneric enumeration element
typedef struct yamlsettings_enumeration_element_{
    char name[YAMLSETTINGS_VAR_VALUE_MAX_LENGTH];
    char second_name[YAMLSETTINGS_VAR_VALUE_MAX_LENGTH];
    int value;
} yamlsettings_enumeration_element;

// Definition of a generic enumeration type that contains an array of enumeration elements
typedef struct yamlsettings_enumeration_definition_{
    char name[YAMLSETTINGS_VAR_VALUE_MAX_LENGTH];
    int nelements;
    yamlsettings_enumeration_element elements[45];
} yamlsettings_enumeration_definition;  

// Meta function that it will be use for enumeration

// funtion typedef for set a value from a enumeration. The value is the value in yamlsettings_enumeration_element structure. 
// With this way you can implement specific way to set your enumeration
typedef int (*yamlsettings_data_type_enumeration_save)(void *mem_pos,int position, int value);
// funtion typedef for retrieve a value from a enumeration. The value returned is the value in yamlsettings_enumeration_element structure. 
// With this way you can implement specific way to retrieved your enumeration
typedef int (*yamlsettings_data_type_enumeration_retrieve)(void *mem_pos,int position);
// Implementation of yamlsettings_data_type_enumeration_save for integer values
int yamlsettings_data_type_generic_enumeration_save_integer(void *mem_pos,int position, int value);
// Implementation of yamlsettings_data_type_enumeration_retrieve for integer values
int yamlsettings_data_type_generic_enumeration_retrieve_integer(void *mem_pos,int position);
// Generic function to set a enumeration. You can use for easily define enumeration. See documentation about how to define a enumeration
int yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_enumeration_save save_function, yamlsettings_enumeration_definition *ged, yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file);
// Generic function to get a value string from a enumeration. You can use this function for define quickly specific enumerations. See documentation about how to define enumerations.
char *yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_enumeration_retrieve function_retrieve, yamlsettings_enumeration_definition *ged, void *mem_pos,int position, int maxlength);


#ifdef	__cplusplus
}
#endif

#endif /* _YAMLSETTINGS_DATA_TYPES_H */ 
