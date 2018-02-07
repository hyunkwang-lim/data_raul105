/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file   grasp_settings_assign_data.h
 * @author David Fuertes (david.fuertes@univ-lille1.fr)
 * @date   August, 2013
 * @brief  Specific function to store parameters
 *
 * Function to store special parameters highly dependent on the problem
 */

#ifndef GRASP_SETTINGS_ASSIGN_DATA_CUSTOM_H
#define	GRASP_SETTINGS_ASSIGN_DATA_CUSTOM_H

#ifdef	__cplusplus
extern "C" {
#endif
    


// Pre-processed

// Saving dimension of ndim part
//void grasp_settings_save_ndim_part(GNode *root,yamlsettings_dictionary_t *dictionary);

// Post-processed
// Calculate and set values of NDIM dimensions and transform contraints to arrays of KPARS
void grasp_settings_calculate_ndim_part(yamlsettings_dictionary_t *dictionary);
// Calculate iwl and key values based on keyBIN value
void grasp_settings_calculate_iwl_and_key(yamlsettings_dictionary_t *dictionary);
// Set input method: coordinates_refernece and cordinates_type from values in dictionary.
void grasp_settings_input_method(yamlsettings_dictionary_t *dictionary);
// Set to true ISTOP if simulated data has been specified 
void grasp_settings_simulated_sdata(yamlsettings_dictionary_t *dictionary);
// Deactivate output if controller is not going to perform retrieval
void grasp_settings_controller_perform_retrieval(yamlsettings_dictionary_t *dictionary);
#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_SETTINGS_ASSIGN_DATA_CUSTOM_H */