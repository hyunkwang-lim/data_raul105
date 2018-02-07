/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */
#include <stdio.h>
#include "grasp_input_transformer_segment_imagedat.h"
#include "../../../input/grasp_input.h"
#include "../../../input/grasp_input_functions.h"

#include "yamlsettings/yamlsettings.h"
#include "yamlsettings/yamlsettings_dictionary.h"
#include "yamlsettings/yamlsettings_assign_data.h"
#include "yamlsettings/yamlsettings_validators.h"

grasp_input_transformer_t grasp_input_transformer_segment_imagedat(){
    grasp_input_transformer_t x;
    
    x.init=grasp_input_transformer_segment_imagedat_init;
    x.function=grasp_input_transformer_segment_imagedat_function;
    x.close=grasp_input_transformer_segment_imagedat_close;

    return x;
}

int grasp_input_transformer_segment_imagedat_init(grasp_settings *settings, grasp_tile_description_t *input_information){
    if(settings->retrieval.INPUT==true) {
        fprintf(stderr,"SEGMENT_IMAGEDAT TRANSFORMER ERROR: The transformer is not compatible with retrieval.debug.use_internal_initial_guess equal to TRUE. Please, set that setting to FALSE if you want to use the SEGMENT_IMAGEDAT transformer");
        abort();
    }

    return 0;
}

int grasp_input_transformer_segment_imagedat_close(void){
    return 0;
}

int grasp_input_transformer_segment_imagedat_function(grasp_settings *settings,grasp_segment_t *segment){
    FILE *f ;
    int iparam;
    float read_value;
    int read_nparam;
    int ipixel;
    
    f = fopen( settings->input.transformer.segment_imagedat.input_file,"r");
    
    if (f==NULL){
        fprintf(stderr,"SEGMENT_IMAGEDAT TRANSFORMER ERROR: %s does not exist or it is not a valid file", settings->input.transformer.segment_imagedat.input_file);
        abort();
    }
    
    for(iparam=0;iparam<settings->retrieval.KNSING;iparam++){
        if(fscanf(f, "%d", &read_nparam )!=1){
            fprintf(stderr,"SEGMENT_IMAGEDAT TRANSFORMER ERROR: Problem parsing imagedat file, please, review the format.");
            abort();
        }
        if(iparam+1!=read_nparam){
            fprintf(stderr,"SEGMENT_IMAGEDAT TRANSFORMER ERROR: Problem parsing imagedat file, please, review the format. Expected to read parameter %d but %d was found",iparam+1, read_nparam );
            abort();
        }

        for(ipixel=0;ipixel<segment->sdata.npixels;ipixel++){
            if(fscanf(f, "%f", &read_value ) != 1 )
            {
                fprintf(stderr,"SEGMENT_IMAGEDAT TRANSFORMER ERROR: Problem parsing imagedat file, please, review the format. File finished before expected");
                abort();
            }
            if(read_value==0){
                read_value=0.1e-9;
            }
            segment->iguess[ipixel][iparam]=read_value;
        }
    }
    
    fclose(f);

    return 0;
}

grasp_settings_parameter_array *grasp_input_transformer_settings_segment_imagedat(grasp_settings *settings){
    int i;
    
    // Static definition of a dictionary
    yamlsettings_parameter parameters[]= {
         /*                 name                                     , memory direction                               , counter memory direction                     , func to set variable type (length of value)             , initial value                               ,  number of elements          , allow_array             , parameter description                                                                                                                                                                                                                                                                                                                ,                 validator 1                 ,                   validator 2                         ,                    validator 3                            , {input function, output function, assigned}   */
        {"file"            , &settings->input.transformer.segment_imagedat.input_file      , NULL                                  , YS_DATA_TYPE_FILE_PATH(_GBL_FILE_PATH_LEN)              , {ys_d_all,{""}}, {0,{}}       , YS_PARAM_SCALAR         , "File which contains initial guess information in classic input.dat format"                                                                                                                                                                                                                                                          , {                                                                                                                                                              },YAMLSETTINGS_END_VAR  },         
    };
    grasp_settings_parameter_array *result;
  
    result=grasp_settings_parameter_array_allocate(sizeof(parameters)/sizeof(yamlsettings_parameter));
    
    for(i=0;i<sizeof(parameters)/sizeof(yamlsettings_parameter);i++){
       yamlsettings_copy_parameter(&(parameters[i]),&(result->parameters[i]));
    }
    
    return result;
}
