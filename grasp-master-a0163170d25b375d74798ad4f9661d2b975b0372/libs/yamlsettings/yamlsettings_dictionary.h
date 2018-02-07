/**
 * @file   yamlsettings_dictionary.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  List of valid configuration parameters and its description and validators
 *
 * This definitions helps to manipulate a dictionary. A dictionary is a list of parameters
 * (yamlsettings structure) iterable. Also it has a description information and validators 
 */
#ifndef _YAMLSETTINGS_DICTIONARY_H
#define _YAMLSETTINGS_DICTIONARY_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "yamlsettings_validators.h"
#include "yamlsettings_defaults.h"
#include "yamlsettings_data_types.h"
#include "yamlsettings_error.h"
#include "yamlsettings_input_yaml.h"
#include <glib.h>
    
#define YAMLSETTINGS_FILE_UNSET 0 
#define YAMLSETTINGS_FILE_UNKNOWN 1
#define YAMLSETTINGS_FILE_COMMAND_LINE 2
#define YAMLSETTINGS_FILE_MAIN_CONFIGURATION 3
    
// Maximum value of dimensions that a parameter can have. 
// It is use to dump settings in yaml format (reverse process)
#define YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS 10

typedef struct yamlsettings_parameter_dimensions_{
    int ndimensions;
    int dimension_size[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
}yamlsettings_parameter_dimensions;

// Represent if a paramter is a single scalar, an array or could be both
typedef enum { YS_PARAM_SCALAR, YS_PARAM_OPTIONAL_ARRAY, YS_PARAM_ARRAY } yamlsettings_parameter_allow_array;

/* SINGLE PARAMETER IN A DICTIONARY */
typedef struct yamlsettings_parameter_{
    const char *name;
    void *mem_pos;
    int  *counter_mem_pos; // Memory position of var that store the number of elements    
    yamlsettings_data_type_set data_type_set;//int (*data_type_set)(void *mem_pos, int nelements,  char *name, char *data, int position, int maxlength);   
    yamlsettings_data_type_get data_type_get;//char *(*data_type_get)(void *mem_pos,int position, int nelements);                           
    int maxlength;
    yamlsettings_default default_value;
    yamlsettings_parameter_dimensions  dimensions; // Number of values. 1 if a unique value. N if is an array. 0 if is dynamic array 
    yamlsettings_parameter_allow_array allow_array;
    const char *description;
    yamlsettings_validator validators[YAMLSETTINGS_MAX_VALIDATORS];
    //int  counter_nelements;
    void  *additional_information; // Private record for the developer. If you need extra information for your own custom input or debug function you can link it here to a void pointer
    int (*input_function)(yamlsettings_status *status, struct yamlsettings_parameter_ *param, GNode *node, void *settings);
    char *(*debug_function)(struct yamlsettings_parameter_ *param); 
    char **settings_file; // Developer has not to define it, can use it. The library is the responsible to set it
    bool name_deallocatable;  // if is true yamlsettings_dictionary_destroy will free memory of char name. You can use it if you don't use constant char. The function yamlsettings_copy_parameter can change this value automatically
} yamlsettings_parameter;


/* DICTIONARY */

// YAMLSETTINGS_END_VAR is used in dictionary definition to symplify it
#define YAMLSETTINGS_END_VAR NULL, yamlsettings_parameter_default_set, yamlsettings_parameter_default_debug, 0, false

// Dictionary definition as an array of parameter with nparameters length. Each parameter
// has an associated yamlsetting var
struct yamlsettings_dictionary_t_{
    char *files[YAMLSETTINGS_MAXIMPORTEDFILES]; // 0: unset (YAMLSETTINGS_FILE_UNSET), 1: unknown (YAMLSETTINGS_FILE_UNKNOWN), 2: command line, working directory(YAMLSETTINGS_FILE_COMMAND_LINE), 3: main configuration file(YAMLSETTINGS_FILE_MAIN_CONFIGURATION), 4...255: files
    yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES];
    int nparameters;
    void *settings;
    yamlsettings_status status;
    yamlsettings_parameter *parameters;
};  
    
// Set initial values for parameters that it will be used in read process (for example, files array)
void yamlsettings_initialize_dictionary(yamlsettings_dictionary_t *dictionary, const char *inputFile);

// Deallocate dictionary memory
void yamlsettings_dictionary_complete_destroy(yamlsettings_dictionary_t *dictionary);

// Deallocate dictionary keeping settings structure
void yamlsettings_dictionary_destroy(yamlsettings_dictionary_t *dictionary);

// Set default values in a dictionary structure. Each parameter of yamlsettings structure
// will have default value especified in dictionary.
void yamlsettings_dictionary_set_defaults(yamlsettings_dictionary_t *dictionary);

// Find an option in yamlsettings_parameters array and return its index (-1 if can't find)
int yamlsettings_dictionary_find_parameter_by_name(char *name, yamlsettings_dictionary_t *dictionary);

// Return index of a file in dictionary given filename. Return -1 if file does not exist
int yamlsettings_dictionary_index_of_file(yamlsettings_dictionary_t *dictionary, const char *filename);

// Function that validate entirely a dictionary
int yamlsettings_dictionary_validate(yamlsettings_dictionary_t *dictionary);

// This function calculate the total number of elements of a parameter, using yamlsettings_parameter_dimensions
int yamlsettings_parameter_number_of_elements(yamlsettings_parameter *param);

// This function calculate the total number of elements of the counter of a parameter, using yamlsettings_parameter_dimensions
int yamlsettings_parameter_counter_number_of_elements(yamlsettings_parameter *param);

// Return number of dimensions of the counter. This function take care of don't return negative numbers
int yamlsettings_parameter_counter_number_of_dimensions(yamlsettings_parameter *param);

// Return number of elements that are indexes in the name of the parameter. Example if parameters is x[1].y[2].v=[1,2,3] the result will be 2
int yamlsettings_parameter_number_of_dimensions_in_string(yamlsettings_parameter *param);

// Return an allocated string with the parameter name adding indexes between brackets []
char *yamlsettings_parameter_get_with_indexes(const char *parameter_name, int nindexes, int *indexes);

// Return position indexes and maximums from a parameter given a specific node
void yamlsettings_parameter_get_indexes(int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], int maximums[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], yamlsettings_parameter *param,GNode *node);

// Return the position of the element in a linear array (for set in memory) for a specific node that represents a parameter
int yamlsettings_parameter_get_position_in_linear_array(yamlsettings_parameter *param,GNode *node);

// Return the index position of the counter in linear array for a specific node that represents a parameter
int yamlsettings_parameter_get_counter_position_in_linear_array(yamlsettings_parameter *param,GNode *node);

// Return true if the parameter in completely unset (in all its dimensions)
bool yamlsettings_parameter_is_unset(int param_index, yamlsettings_dictionary_t *dictionary);

// Return the index of file which this parameter came from. This the parameter come from different settings files it will return YAMLSETTINGS_FILE_UNKNOWN
int yamlsettings_parameter_file_index(int param_index, yamlsettings_dictionary_t *dictionary);

/**
 * This function return the category name from a parameter of a specific level
 * @param parameter Parameter name that is going to be analyzed
 * @param level Number of de level which will be returned
 * @return allocated string with category name of the parameter at the specific level
 */
char *yamlsettings_parameter_name_at_level(const char *parameter, int level);

/**
 * This function split a parameter name getting the first (as much as level) categories 
 * @param parameter The name of the parameter
 * @param level Number of category blocks to get before split the parameter name
 * @return Allocate string with the parameter name until level blocks
 */
char *yamlsettings_parameter_name_until_level(const char *parameter, int level);

/**
 * This function returns the number of categories that a parameter name contains. For example, in case of the parameter x.y.z it will return 3
 * @param param_name string will be analyzed
 * @return number of blocks in the parameter
 */
int yamlsettings_parameter_number_of_blocks(const char *param_name);

/**
 * This function counts how many arrays are in the parameter name (only analyzing the name). For example x.y[] will return 1.
 * @param param_name the name which will be studied
 * @return number of array elements in the string
 */
int yamlsettings_parameter_name_number_of_dimensions(const char *param_name);


/**
 * This function return number of parameters which are inside certain block
 * For example, this yaml:
 * xxx:
 *     yyy:
 *     zzz:
 * it is defined like this:
 * xxx.yyy
 * xxx.zzz
 * Given this dictionary, this function for parameter 1 and level 1 will return two elements (yyy and zzz). 
 * It is mandatory to ask for first element with the title. If not an assertion will launched. 
 * @param dictionary
 * @param parameter
 * @param level
 * @return 
 */
int yamlsettings_parameter_number_of_parameters_in_a_block(yamlsettings_dictionary_t *dictionary,int parameter, int level);

/**
 * This functions return true if in a parameter in a specific level has not information different from default value and it is not mandatory (checked by validator)
 * for example, if we have parameters "x.y.z" and "x.y.w" when level is 2 this function will check if z and w has values different from default value
 * It is mandatory to call this function with first parameter which contain this title. In our example "x.y.z" because it is first definition of y
 * @param dictionary Current dictionary
 * @param parameter Parameter number in dictionary
 * @param level Level which will be studied.
 * @return true if it has no default values (the parameter contains useful information), otherwise false
 */
bool yamlsetting_parameter_is_not_omissible(yamlsettings_dictionary_t *dictionary,int parameter, int level);

/**
 * Private function which dump a dictionary in yaml format recursively.
 * @param dictionary A dictionary which will be dumped in yaml format
 * @param level Current indentation level you are studying
 * @param ilevelstart from which parameter you start to dump
 * @param indexes Current indexes visited (if it is an array)
 * @param indexstop to which parameter you are going to dump
 * @param print_defaults If you want to print all nodes or only set nodes
 * @return 0 if everything was ok
 */
int yamlsettings_dictionary_dump_recursive(FILE *f, yamlsettings_dictionary_t *dictionary,int level, int ilevelstart, int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], int indexstop, bool print_defaults);

#ifdef	__cplusplus
}
#endif

#endif /* _YAMLSETTINGS_DICTIONARY_H */
