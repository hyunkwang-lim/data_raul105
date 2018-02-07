/**
 * @file   yamlsettings_process_yaml.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Low-level functions to read yaml files
 *
 * This functions read yamls files in a file and return a GTree
 */
#ifndef YAMLSETTINGS_PROCESS_YAML_H
#define	YAMLSETTINGS_PROCESS_YAML_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include <yaml.h>
#include <glib.h>
#include "yamlsettings_input_yaml.h"
#include <stdbool.h>

    
// Recursive function that read yaml file using libyaml
int yamlsettings_yaml_file_process(yamlsettings_status *status, yaml_parser_t *parser, GNode *data, char file_index, yamlsettings_parser_mode mode);

// Extract imported files from a tree and insert it in the list of files for process. If ignore is true the import statament will be removed from true but no added to list of processing files.
int yamlsettings_process_imported_files(bool ignore, yamlsettings_status *status, GNode *cfg, char *files[YAMLSETTINGS_MAXIMPORTEDFILES],yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path);

// Extract template files from a tree and insert it in the list of files for process. If ignore is true the import statament will be removed from true but no added to list of processing files.
int yamlsettings_process_template_files(bool ignore, yamlsettings_status *status, GNode *cfg, char *files[YAMLSETTINGS_MAXIMPORTEDFILES],yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path);

// Add a file to array of files
int yamlsettings_yaml_file_add(yamlsettings_status *status, const char *file, char *files[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path, yamlsettings_parser_mode mode);

// Main function to check if all arrays are completed ( i.e: if there is node 2, must exist node 1)
// This function assumes that the tree is normalized
void yamlsettings_check_arrays(yamlsettings_status *status, GNode *root);

// Function call for each array node in order to check if the array is completed (used by g_node_traverse)
gboolean yamlsettings_check_array_elements(GNode *node, gpointer data);

// Function call for each node find by yamlsettings_check_array_elements for check if there are other similar
gboolean yamlsettings_check_array_elements_find_name(GNode *node, gpointer data);

// This function return a GNode with input arguments like data
GNode *yamlsettings_create_node(int type, int nelements, int element_position, char file_index, char *data);

// This function remove the data assigned to the nodes. Call it before destroy a tree for each node. The function has the interface for being calling by
// g_node_traverse function but don't need any extra information, data can be NULL
gboolean yamlsettings_destroy_node_data(GNode *node, gpointer data);

// This function deallocate a tree taking into account the memory taken by the nodes.
void yamlsettings_destroy_tree(GNode *node);


// Internal structure used in yamlsettings_find_node and yamlsettings_find_node_traverse functions
typedef struct yamlsettings_find_node_traverse_t_{
    const char *name;
    int nnodes;
    GNode **nodes;
}yamlsettings_find_node_traverse_t;

// Function that is call for each node by yamlsettings_find_node and store in yamlsettings_find_node_traverse_t (data)
// each node with the same name (skipping leaves)
gboolean yamlsettings_find_node_traverse(GNode *node, gpointer data);

// This function will extract the nodes in string_parameter and add them to root.
void yamlsettings_add_parameter_to_a_tree(GNode *root, const char string_parameter[]);

// This function receive an array of extra parameters defined in arrays (Parameter example: retrieval.number_wavelengths=3) an return a tree without normalize with all parameters
GNode *yamlsettings_tree_from_extra_parameters(int nparameters, char const *parameters[]);

// This function return a root of tree transforming parameter string. Parameter example: retrieval.number_wavelengths=3
GNode *yamlsettings_parameter_to_tree(char const *parameter);

// Debug function to print tree information
gboolean yamlsettings_tree_debug(GNode *node, gpointer data);


#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_PROCESS_YAML_H */
