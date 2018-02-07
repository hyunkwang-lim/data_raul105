/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <string.h>
#include <stdlib.h>
#include "yamlsettings/yamlsettings_input_yaml.h"
#include "yamlsettings/yamlsettings_data_types.h"
#include <grasp/utils.h>
#include "math.h"
#include "yamlsettings/yamlsettings_dictionary.h"
#include "mod_par_inv.inc"
#include "grasp_settings.h"
#include "../global/grasp_retrieval_characteristic_type.h"


void grasp_settings_calculate_ndim_part(yamlsettings_dictionary_t *dictionary){
    int idim1,idim2,idim3;
    int tmp, tmpf, nsd, index, sd_index;
    grasp_settings *settings;
    char error[512];
    
    settings=(grasp_settings *)dictionary->settings;    
    
    // Initialize NDIM part
    settings->retrieval.NDIM.n1=0;    
    for(idim1=0;idim1<_KIDIM1; idim1++){
        settings->retrieval.NDIM.n2[idim1]=0;
        for(idim2=0;idim2<_KIDIM2; idim2++){
            settings->retrieval.NDIM.n3[idim1][idim2]=0;
            settings->retrieval.NDIM.ISTARSING[idim1][idim2]=0;
        }
    }        
    
    // Finding maximums for NDIM
    for (idim1 = 0; idim1 < _KIDIM1; idim1++) {
        for (idim2 = 0; idim2 < _KIDIM2; idim2++) {
            for (idim3 = 0; idim3 < _KIDIM3; idim3++) {
                if(settings->tmp.TAPSING[idim1][idim2][idim3]>FLT_MIN){
                    settings->retrieval.NDIM.n3[idim1][idim2]++;
                }else{
                    break;
                }
            }
            if(settings->tmp.TAPSING[idim1][idim2][0]>FLT_MIN){
                settings->retrieval.NDIM.n2[idim1]++;
            }else{
                break;
            }
        }
        if(settings->tmp.TAPSING[idim1][0][0]>FLT_MIN){
            settings->retrieval.NDIM.n1++;
        }else{
            break;
        }
    }
    
    // Set other parameters
    tmp=1; 
    tmpf=1;
    nsd=0;   
    settings->retrieval.KNSING=0;
    index=0;
    for(idim1=0;idim1<settings->retrieval.NDIM.n1; idim1++){
        if(settings->retrieval.NDIM.n2[idim1]>_KIDIM2){
            sprintf(error,"Error allocating NDIM parameters. They are more than constants allows (KIDIM2[%d]=%d)",idim1,settings->retrieval.NDIM.n2[idim1]);
            yamlsettings_error_add_parse_error(&(dictionary->status), error, YAMLSETTINGS_FILE_UNKNOWN, YS_ERROR); 
            break;
        }   
        
        for(idim2=0; idim2<settings->retrieval.NDIM.n2[idim1]; idim2++){
            if(settings->retrieval.NDIM.n3[idim1][idim2]>_KIDIM3){
                sprintf(error,"Error allocating NDIM parameters. There are more values than constants allows (KIDIM3[%d][%d]=%d)",idim1,idim2,settings->retrieval.NDIM.n3[idim1][idim2]);
                yamlsettings_error_add_parse_error(&(dictionary->status), error, YAMLSETTINGS_FILE_UNKNOWN, YS_ERROR); 
                break;
            }    
            settings->retrieval.NDIM.ISTARSING[idim1][idim2]=tmp;
            tmp=tmp+settings->retrieval.NDIM.n3[idim1][idim2];
            if(settings->retrieval.NDIM.par_retr[idim1]==true){
                tmpf=tmpf+settings->retrieval.NDIM.n3[idim1][idim2];
            }
                
            if(idim1==0) nsd++;
            for (idim3 = 0; idim3 < settings->retrieval.NDIM.n3[idim1][idim2]; idim3++) {
                // Settings kndim parameters
                assert(index>=0 && index<_KPARS);
                settings->retrieval.APSING[index]=settings->tmp.TAPSING[idim1][idim2][idim3];
                settings->retrieval.APSMIN[index]=settings->tmp.TAPSMIN[idim1][idim2][idim3];
                settings->retrieval.APSMAX[index]=settings->tmp.TAPSMAX[idim1][idim2][idim3];
                settings->retrieval.IWW_SINGL[index]=settings->tmp.TIWW_SINGL[idim1][idim2][idim3];
                index++;
            }            
        }
    }
    settings->retrieval.KNSING=tmp-1;
    settings->retrieval.KNSINGF=tmpf-1;
    sd_index=grasp_parameters_index_of_parameter_type_by_kind_of_parameter(&settings->retrieval.NDIM, par_type_SD_beg, par_type_SD_end);
    if (sd_index>=0){
        settings->retrieval.NSD=grasp_parameters_number_of_modes_of_parameter(&settings->retrieval.NDIM, settings->retrieval.NDIM.par_type[sd_index]);
    }else{
        // It means there is an error in settings definition but here
        // we are not going to trigger any action. Other validator will
        // show the error
        settings->retrieval.NSD=1; 
    }
}


void grasp_settings_calculate_iwl_and_key(yamlsettings_dictionary_t *dictionary){
    int idim1;
    grasp_settings *settings;
    
    settings=(grasp_settings *)dictionary->settings;
    
      for(idim1=0;idim1<settings->retrieval.NDIM.n1;idim1++){
        if(settings->retrieval.NDIM.par_type[idim1] > par_type_SD_beg && settings->retrieval.NDIM.par_type[idim1] < par_type_SD_end){      
          if(settings->retrieval.NDIM.par_type[idim1] == par_type_SD_TB) { // triangle  bins
            settings->retrieval.DLSF.IWL = 0;
            settings->retrieval.DLSF.key = 0;       
          }else if(settings->retrieval.NDIM.par_type[idim1] == par_type_SD_LB) { // precalculated lognormal bins      
            settings->retrieval.DLSF.IWL = 1;
            settings->retrieval.DLSF.key = 2;               
          }else if(settings->retrieval.NDIM.par_type[idim1] == par_type_SD_LN){ // parameters of bi-modal lognormal SD
            settings->retrieval.DLSF.IWL = 0;
            settings->retrieval.DLSF.key = 3;
          }
        }
    }    
}

void grasp_settings_input_method(yamlsettings_dictionary_t *dictionary){
    int center_latitude_index, center_longitude_index, corner_row_index, corner_col_index;
    
    center_latitude_index=yamlsettings_dictionary_find_parameter_by_name("input.center.latitude", dictionary);
    center_longitude_index=yamlsettings_dictionary_find_parameter_by_name("input.center.latitude", dictionary);
    corner_row_index=yamlsettings_dictionary_find_parameter_by_name("input.corner.row", dictionary);
    corner_col_index=yamlsettings_dictionary_find_parameter_by_name("input.corner.column", dictionary);
    
    assert(center_latitude_index>=0);
    assert(center_longitude_index>=0);
    assert(corner_row_index>=0);
    assert(corner_col_index>=0);
    
    if(yamlsettings_parameter_is_unset(center_latitude_index, dictionary)==false && yamlsettings_parameter_is_unset(center_longitude_index, dictionary)==false){
        strcpy(((grasp_settings *)dictionary->settings)->input.coordinates_refence,"center");
        strcpy(((grasp_settings *)dictionary->settings)->input.coordinates_type,"latlon");
    }
    if(yamlsettings_parameter_is_unset(corner_row_index, dictionary)==false && yamlsettings_parameter_is_unset(corner_col_index, dictionary)==false){
        strcpy(((grasp_settings *)dictionary->settings)->input.coordinates_refence,"corner");
        strcpy(((grasp_settings *)dictionary->settings)->input.coordinates_type,"rowcol");
    }        
}


void grasp_settings_simulated_sdata(yamlsettings_dictionary_t *dictionary){
    int sim,i;

    sim=yamlsettings_dictionary_find_parameter_by_name("retrieval.debug.simulated_sdata_file", dictionary);        
    assert(sim>=0);

    if(yamlsettings_parameter_is_unset(sim, dictionary)==false){ // If sim sdata is activated (set, no default value), istop is set to true 
        ((grasp_settings *)(dictionary->settings))->retrieval.ISTOP=true;
    }else{ // If the parameter is not set, it is set to empty fortran string
        for (i = 0; i < _GBL_FILE_PATH_LEN; i++) {
            ((grasp_settings *)(dictionary->settings))->retrieval.sdata_sim_file[i]=' ';
        }
    }
}

void grasp_settings_controller_perform_retrieval(yamlsettings_dictionary_t *dictionary){
    grasp_settings *settings;
    
    settings=(grasp_settings *)dictionary->settings; 
    if(settings->controller.perform_retrieval==false){
        settings->output.ncurrent_output_function=0;
        settings->output.nsegment_output_function=0;
        settings->output.ntile_output_function=0;
    }
}
