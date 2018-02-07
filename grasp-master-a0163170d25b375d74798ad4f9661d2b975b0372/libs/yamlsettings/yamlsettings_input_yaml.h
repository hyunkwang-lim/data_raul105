/**
 * @file   yamlsettings_input_yaml.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  High-level API to read yaml files
 *
 * This functions read yamls files in a file and return a GTree
 */

#ifndef YAMLSETTINGS_INPUT_YAML_H
#define	YAMLSETTINGS_INPUT_YAML_H

#ifdef	__cplusplus
extern "C" {
#endif

    
#include <yaml.h>
#include <glib.h>
#include "yamlsettings_error.h"
#include "yamlsettings_parse_settings_file_mode.h"
    
// Define maximum length of parameter (e.g. retrieval_input_structure.knsing) name.
#define YAMLSETTINGS_MAXSIZEPARAMETERNAME 215
// Maximum of files that can be imported in total, adding imports recursive
//#define YAMLSETTINGS_MAXIMPORTEDFILES 252
#define YAMLSETTINGS_MAXIMPORTEDFILES 64
// Maximum of nodes returned by yamlsettings_find_node
#define YAMLSETTINGS_PROCESS_FIND_MAX_NODES 50

// Types of nodes: variable, value, sequence
// VAR = Node
// VAL = Value, leaf
// SEQ = Array in leaf
enum storage_flags { VAR, VAL, SEQ }; // "Store as" switch

// This enumeration defines where insert the node in the tree, before or after. import=after and template=before. 
// If the value is set to YS_PARSE_MODE_PREPEND the library will allow to overwrite it
typedef enum  { YS_PARSE_MODE_PREPEND , YS_PARSE_MODE_APPEND } yamlsettings_parser_mode;

// Definition of information can be stored in a GNode
typedef struct yamlsettings_node_element_{
    int type; //{ VAR, VAL, SEQ }
    int nelements; // If value is SEQ, this element contains the size of array
    int element_position; // IF is a SEQ, this is the number of the element in this array
    unsigned char file_index;
    gchar *data;
}yamlsettings_node_element;

/** 
 * @fn int yamlsettings_input_yaml_read(char *inputFile, GNode **cfg)
 * @brief opens a YAML configuration file and return a GTree
 * 
 * @param[out] *status List of errors during the parse and validation process
 * @param[out] **cfg Pointer to config in a GTree
 * @param[in] nparamters Number of input paramters
 * @param[in] parameters Input parameters. First one have to be input file
 * @param[out] files List of processed files (including import and template files)
 * @param[out] file_modes If file is read before or after main settings file 
 * @param[in] settings_file_mode specified if the settings file is mandatory or it can work only with arguments (parameters)
 * @return Status of process. 0 if it is ok. -1 if there are errors.
 */
int yamlsettings_input_yaml_read(yamlsettings_status *status, GNode **cfg, int nparameters, char const *parameters[], char *files[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_settings_file_mode settings_file_mode);

/** 
 * @fn char *yamlsettings_node_name(GNode *node)
 * @brief return a name of a parameter in a node in a GTree 
 * 
 * @param[in] *node Pointer to specific node
 * @return String containing real name
 */
char *yamlsettings_parameter_name(GNode *node);

/** 
 * @fn char *yamlsettings_node_name(GNode *node)
 * @brief return a name of a parameter in a node in a GTree removing number of element in arrays ( [x] by [] )
 * 
 * @param[in] *node Pointer to specific node
 * @return String containing name normalized
 */
char *yamlsettings_parameter_name_normalized(GNode *node);

// Return name of node without array elements (array[2] => array)
char *yamlsettings_name_of_array_node (const char *name);

// Return index of array node. If node is not an array, return -1. If index is invalid return 0. 
int yamlsettings_index_of_array_node(const char *name);

// Return true if node has array name, for example: param[3] is an array.
gboolean yamlsettings_is_node_array(const char *name);

// Return true if node has array in complete name, taking into account parent nodes.
gboolean yamlsettings_is_parameter_array(GNode *node);


// Find a node by name (using: yamlsettings_find_node_traverse that use yamlsettings_node_name_normalized so you can use for example ".node[]."). Return number of nodes 
int yamlsettings_find_node(GNode *root, const char *name, GNode **dest);

// Return first parent of node that is a array eg. element[3]
GNode *yamlsettings_first_array_parent(GNode *node);


#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_INPUT_YAML_H */
