/**
 * @file   yamlsettings_defaults.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2016
 * @brief  Functions to define default values of the tree
 *
 * To define default value of a dictonary a set of functions are defined by
 * yamlsettings library. The user can use this functions or define its own one
 * 
 ***** IMPORTANT NOTE
 * first argument of all default function has to be the "fill value"
 * because it is used in yamlsettings_assign_data:yamlsettings_parameter_default_set
 * in case a nose is an empty array [] (filling with param->default_value.arguments[0]) 
 */

#ifndef YAMLSETTINGS_DEFAULTS_H
#define	YAMLSETTINGS_DEFAULTS_H

#ifdef	__cplusplus
extern "C" {
#endif

    
#define YAMLSETTINGS_DEFAULTS_MAX_ARGUMENTS 5
#define YAMLSETTINGS_DEFAULTS_MAX_ARGUMENT_SIZE 500
    
typedef struct yamlsettings_parameter_ yamlsettings_parameter ;
typedef struct yamlsettings_status_ yamlsettings_status;

typedef struct{
    int (*defaults)(int param_index, yamlsettings_parameter *parameter, yamlsettings_status *status);
    char arguments[YAMLSETTINGS_DEFAULTS_MAX_ARGUMENTS][YAMLSETTINGS_DEFAULTS_MAX_ARGUMENT_SIZE];
}yamlsettings_default;


/**
 * This set default function set all position with same values which is received as first argument
 * @param param_index
 * @param parameter
 * @param status
 * @return 
 */
int ys_d_all (int param_index,yamlsettings_parameter *parameter, yamlsettings_status *status);

/**
 * To set an array. Set all values with same value (first argument) and then set specifically rest of values in order.
 * After first arguemnt with the fill value, then has the be the number of elements to be set and then
 * pairs of position, value. Example: "0","2","0:true","1:false"
 * @param param_index
 * @param parameter
 * @param status
 * @return 
 */
int ys_d_array (int param_index,yamlsettings_parameter *parameter, yamlsettings_status *status);



#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_DEFAULTS_H */
