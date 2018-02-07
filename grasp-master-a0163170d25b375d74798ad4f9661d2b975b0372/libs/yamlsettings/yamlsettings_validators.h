/**
 * @file   yamlsettings_validators.h
 * @Author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Validator functions used in dictionary
 *
 */

#ifndef YAMLSETTINGS_VALIDATORS_H
#define	YAMLSETTINGS_VALIDATORS_H

#ifdef	__cplusplus
extern "C" {
#endif

    
#define YAMLSETTINGS_MAX_VALIDATORS 5
#define YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS 2
#define YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE 500
    
typedef struct{
    int (*validator)(); 
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE];
}yamlsettings_validator;


typedef struct yamlsettings_dictionary_t_ yamlsettings_dictionary_t ;

// The parameter is mandatory. It don't need arguments
int yamlsettings_validator_mandatory         (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is mandatory if other parameter is true. Name of other parameter first argument name
int yamlsettings_validator_mandatory_if_parameter_true (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is mandatory if other parameter is parameter enum has value. Argument 1=argument name. Argument 2=enum value 
int yamlsettings_validator_mandatory_if_parameter_enum (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is mandatory if other parameter is parameter enum has not value. Argument 1=argument name. Argument 2=enum value 
int yamlsettings_validator_mandatory_if_parameter_enum_not (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is an integer in a strict range. Argument 1=min. Argument 2=max
int yamlsettings_validator_int               (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is an double in a strict range. Argument 1=min. Argument 2=max
int yamlsettings_validator_double            (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is an double in a strict range. Argument 1=min. Argument 2=max
int yamlsettings_validator_float             (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter has a number of elements in a strict range. Argument 1=min. Argument 2=max
int yamlsettings_validator_nelements         (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter has the samennumber of elements that other element. Argument 1=name of second element 
int yamlsettings_validator_same_nelements    (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// Check if arrays always grown (or equal)
int yamlsettings_validator_ascendant_int_array(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);
// The parameter is a directory and it must be exist
int yamlsettings_validator_directory         (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]);

#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_VALIDATORS_H */
