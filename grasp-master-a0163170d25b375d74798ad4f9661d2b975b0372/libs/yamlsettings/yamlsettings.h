/**
 * @file   yamlssettings.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Main interface binded for input settings 
 *
 * Here is defined option parameters and main functions for retrieved it
 * from a configuration file
 * 
 */
#ifndef _YAMLSETTINGS_H
#define _YAMLSETTINGS_H


#ifdef	__cplusplus
extern "C" {
#endif
    
#include <stdbool.h>
#include <glib.h>  
#include "yamlsettings_dictionary.h"
#include "yamlsettings_parse_settings_file_mode.h"
    
// Function to be implemented by preprocess methods. This function will be call after read the file but before store it.
typedef void (*yamlsettings_preprocess)(GNode *root,yamlsettings_dictionary_t *dictionary); 

// Function to be implemented by postprocess methos. This functions will be call after store values in dictionary but before validate it.
typedef void (*yamlsettings_postprocess)(yamlsettings_dictionary_t *dictionary);
    

// Function that read a settings file in Yaml format and return a dictionary. 
// The dictionary must be initializated before. You can use yamlsettings_dictionary_default_set function. Arguments:
// yamlsettings_dictionary_get get_dictionary function to obtain a new dictionary
// number of  parameters that will be read from array of strings instead of from a file
// char const *parameters[]  first value is the main settings file, rest of values are string with parameters that will be added
// int npreprocess Number or preprocessor function that reader have to call
// yamlsettings_preprocess preprocess[] Array of functions
// int npostprocess Number of postprocessed function that reader have to call
// yamlsettings_postprocess postprocess[] Array of postprocessed functions
// Return number of errors in total (validators + parse) in negative
int yamlsettings_read_file(yamlsettings_dictionary_t *dictionary, int nparameters, char const *parameters[], yamlsettings_preprocess preprocess, yamlsettings_postprocess postprocess, yamlsettings_parser_settings_file_mode settings_file_mode);

// Print values that will be used for the application. The values read from file and stored in dictionary.
void yamlsettings_debug_dictionary(FILE *stream, yamlsettings_dictionary_t *dictionary, const char *filter);

// Print help information (parameter name and description) in stream.
void yamlsettings_print_help_information(FILE *stream, yamlsettings_dictionary_t *dictionary, const char *filter);

/**
 * Function which dump a dictionary in yaml format recursively.
 * @param dictionary A dictionary which will be dumped in yaml format
 * @param print_defaults If you want to print all nodes or only set nodes
 * @param reference It is a dictionary with all settings initializated with default values. It is used to know which parameter has different value than default
 * @return 
 */
int yamlsettings_dump(FILE *f, yamlsettings_dictionary_t *dictionary, bool print_defaults, yamlsettings_dictionary_t *reference);

// This function copy parameter in origin in dest
// This function is useful to use static dictionary definition and them move it to a dinamyc one
void yamlsettings_copy_parameter(yamlsettings_parameter *origin, yamlsettings_parameter *dest);

#ifdef	__cplusplus
}
#endif

#endif /*_YAMLSETTINGS_H */
