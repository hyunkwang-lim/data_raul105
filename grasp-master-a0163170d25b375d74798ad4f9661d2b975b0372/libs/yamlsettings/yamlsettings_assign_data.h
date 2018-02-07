/**
 * @file   yamlsettings_assign_data.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  High-level function to store parameter information
 *
 * This file contains assignation functions to store parameters
 */
#ifndef YAMLSETTINGS_ASSIGN_DATA_H
#define	YAMLSETTINGS_ASSIGN_DATA_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include <glib.h>
#include "yamlsettings_dictionary.h"


// This function can be called for each node (g_node_traverse) in order to save in yamlsetting structure
// the information stored in the tree. 
gboolean yamlsettings_assign_node(GNode *node, gpointer data);

// Default function to store a parameter in a mem_position. Simply, find a position a save the value. If is an array also store counter
int yamlsettings_parameter_default_set(yamlsettings_status *status, yamlsettings_parameter *param,GNode *node,void *settings);

// Default function to print values of a parameter. Iterate over each value and call yamlsettings_data_type_print
char *yamlsettings_parameter_default_debug(yamlsettings_parameter *param);



#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_ASSIGN_DATA_H */
