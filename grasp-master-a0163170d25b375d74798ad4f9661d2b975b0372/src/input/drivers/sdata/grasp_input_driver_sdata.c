/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "sdata.h"
#include <string.h>
#include <time.h>
#include "../../grasp_input.h"
#include <grasp/utils.h>
#include <math.h>
#include "../../util/grasp_box.h"
#include "sdata-impl.h"
#include "grasp_input_driver_sdata.h"
#include "yamlsettings/yamlsettings.h"
#include "yamlsettings/yamlsettings_dictionary.h"
#include "yamlsettings/yamlsettings_assign_data.h"
#include "yamlsettings/yamlsettings_validators.h"

static SDATA_HANDLE *handle;


grasp_input_driver_t grasp_input_driver_sdata(){
    grasp_input_driver_t x;
    
    x.init=grasp_input_driver_sdata_init;
    x.get_segment=grasp_input_driver_sdata_get_segment;
    x.close=grasp_input_driver_sdata_close;
    
    return x;
}

int grasp_input_driver_sdata_get_segment(grasp_settings *settings, grasp_segment_t *segment, int col, int row, int itime){    
    grasp_box_t *box;
    bool success;
    SDATA_PIXEL pixel_data;
    SDATA_HEADER header;
    int iwl,ivalidmeas,isurf,ipixel=0;
    int ip; /* index for looping on polarized components (up to 1 for a non polarized channel, up to 3 for polarized channels) */
    int ivm; /* index of valid measurement */
    sensor_data_t *sdata;
    
    sdata=&segment->sdata;    
    box = sdata_load_box(handle, NULL /* default user settings */);
    
    if(box==NULL){
        printf("ERROR: Box could not be loaded\n");
        exit(-1);
    }
    
    sdata_get_header(handle,&header);
    
    grasp_segment_set_nt(sdata, header.nrecords);
    grasp_segment_set_nx(sdata, header.nx);
    grasp_segment_set_ny(sdata, header.ny);

    do{
        success=grasp_box_get_next_pixel(box,&pixel_data);
        if(success==true){
            grasp_segment_set_npixels(sdata,ipixel+1);
            grasp_segment_set_pixel_ix(sdata, ipixel, pixel_data.ix);
            grasp_segment_set_pixel_iy(sdata, ipixel, pixel_data.iy);
            grasp_segment_set_pixel_it(sdata, ipixel, pixel_data.it);
            grasp_segment_set_pixel_cloudy(sdata, ipixel, pixel_data.cloud_flag);
            grasp_segment_set_pixel_irow(sdata, ipixel, pixel_data.irow);
            grasp_segment_set_pixel_icol(sdata, ipixel, pixel_data.icol);
            grasp_segment_set_pixel_file_index(sdata, ipixel, GRASP_INPUT_UNKNOWN_FILEINDEX);
              
            grasp_segment_set_pixel_x(sdata, ipixel, pixel_data.longitude);
            grasp_segment_set_pixel_y(sdata, ipixel, pixel_data.latitude);
            grasp_segment_set_pixel_t(sdata, ipixel, pixel_data.unix_time);
            grasp_segment_set_pixel_out_t(sdata, ipixel, pixel_data.it-1);
            grasp_segment_set_pixel_out_y(sdata, ipixel, pixel_data.iy-1);
            grasp_segment_set_pixel_out_x(sdata, ipixel, pixel_data.ix-1);
            grasp_segment_set_pixel_masl(sdata, ipixel, pixel_data.hgr);
            grasp_segment_set_pixel_hobs(sdata, ipixel, pixel_data.hobs);
            grasp_segment_set_pixel_land_percent(sdata, ipixel, pixel_data.land_percent);
            grasp_segment_set_pixel_nwl(sdata, ipixel, pixel_data.num_wavelengths);
            grasp_segment_set_pixel_ifgas(sdata, ipixel, pixel_data.ifgas);
	    
	    for (ivm = 0 ; ivm < pixel_data.num_hvp ; ivm++) {
                grasp_segment_set_pixel_hvp(sdata, ipixel, ivm, pixel_data.hvp[ivm]);
	    }
	    
            for(iwl=0;iwl<sdata->pixel[ipixel].nwl;iwl++){
                //grasp_segment_set_pixel_meas_ind_wl(sdata, ipixel, iwl, iwl);
                grasp_segment_set_pixel_meas_wl(sdata, ipixel, iwl,pixel_data.wavelengths[iwl]);
                
                grasp_segment_set_pixel_meas_nsurf(sdata, ipixel, iwl, pixel_data.nsurf);
                grasp_segment_set_pixel_meas_gaspar(sdata, ipixel, iwl, pixel_data.gaspar[iwl]);
                grasp_segment_set_pixel_meas_sza(sdata ,ipixel, iwl, pixel_data.thetas[iwl]);
                
                for(isurf=0;isurf<sdata->pixel[ipixel].meas[iwl].nsurf;isurf++){
                    grasp_segment_set_pixel_meas_groundpar(sdata, ipixel, iwl, isurf, pixel_data.groundpar[iwl][isurf]);
                }                
                      
                grasp_segment_set_pixel_meas_nip(sdata, ipixel, iwl, pixel_data.num_meas_types[iwl]);
                for(ip=0;ip<pixel_data.num_meas_types[iwl];ip++){
		  int meas_type;
		  grasp_segment_set_pixel_meas_meas_type(sdata, ipixel, iwl, ip, pixel_data.meas_types[iwl][ip]);
                    
		  meas_type =  pixel_data.meas_types[iwl][ip];

		  if (meas_type == MEAS_TYPE_LS || meas_type == MEAS_TYPE_DP || meas_type == MEAS_TYPE_RL) {
                    grasp_segment_set_pixel_meas_nbvm(sdata, ipixel, iwl, ip, pixel_data.num_hvp);
		  }
		  else {
                    grasp_segment_set_pixel_meas_nbvm(sdata, ipixel, iwl, ip, pixel_data.num_valid_meas[iwl][ip]);
		  }

                grasp_segment_set_pixel_meas_ifcov(sdata, ipixel, iwl, ip, pixel_data.ifcov[iwl][ip]);
                grasp_segment_set_pixel_meas_ifmp(sdata, ipixel, iwl, ip, pixel_data.ifmp[iwl][ip]);
                for(ivalidmeas=0;ivalidmeas<sdata->pixel[ipixel].meas[iwl].nbvm[ip] ;ivalidmeas++){
                    grasp_segment_set_pixel_meas_thetav(sdata, ipixel, iwl, ip, ivalidmeas, pixel_data.thetav[iwl][ip][ivalidmeas]);

                    grasp_segment_set_pixel_meas_phi(sdata, ipixel, iwl, ip, ivalidmeas, pixel_data.dphi[iwl][ip][ivalidmeas]);

                    if(sdata->pixel[ipixel].meas[iwl].ifcov[ip]==1){
                        grasp_segment_set_pixel_meas_cmtrx(sdata, ipixel, iwl, ip, ivalidmeas, pixel_data.cmtrx[iwl][ip][ivalidmeas]);
                    }

                    if (sdata->pixel[ipixel].meas[iwl].ifmp[ip]==1) {
                        grasp_segment_set_pixel_meas_mprof(sdata, ipixel, iwl, ip, ivalidmeas, pixel_data.mprof[iwl][ip][ivalidmeas]);
                    }
                }
                     
                    

                }
                
                for(ivalidmeas=0;ivalidmeas<_NBVM ;ivalidmeas++){ 
                    grasp_segment_set_pixel_meas_tau(sdata, ipixel, iwl, ivalidmeas, pixel_data.tau[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p11(sdata, ipixel, iwl, ivalidmeas, pixel_data.p11[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p12(sdata, ipixel, iwl, ivalidmeas, pixel_data.p12[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p22(sdata, ipixel, iwl, ivalidmeas, pixel_data.p22[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p33(sdata, ipixel, iwl, ivalidmeas, pixel_data.p33[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p34(sdata, ipixel, iwl, ivalidmeas, pixel_data.p34[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p44(sdata, ipixel, iwl, ivalidmeas, pixel_data.p44[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_i(sdata, ipixel, iwl, ivalidmeas, pixel_data.I[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_q(sdata, ipixel, iwl, ivalidmeas, pixel_data.Q[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_u(sdata, ipixel, iwl, ivalidmeas, pixel_data.U[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_p(sdata, ipixel, iwl, ivalidmeas, pixel_data.P[iwl][ivalidmeas]);
                }
                for(ivalidmeas=0;ivalidmeas<_KVERTM ;ivalidmeas++){
                    grasp_segment_set_pixel_meas_ls(sdata, ipixel, iwl, ivalidmeas, pixel_data.LS[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_dp(sdata, ipixel, iwl, ivalidmeas, pixel_data.DP[iwl][ivalidmeas]);
                    grasp_segment_set_pixel_meas_rl(sdata, ipixel, iwl, ivalidmeas, pixel_data.RL[iwl][ivalidmeas]);
                }                
                
            }
            
            ipixel++;
        }
    }while(success==true);
            
    grasp_box_delete(box);
    
    return ipixel;
}

int grasp_input_driver_sdata_close(void){
    int err;
    
    err=sdata_close(handle);
    
    if(err==SDATA_ERROR){
        printf("ERROR: Could not close file: ");
        perror(sdata_get_file_name(handle));
        exit(-1);
    }
    
    return 0; 
}

int grasp_input_driver_sdata_init(grasp_settings *settings, grasp_tile_description_t *input_information){   
   
    // Some validations...
    if(input_information->ninput_files!=1){
        printf("ERROR: SData driver only support 1 input file yet\n");
        abort();
    }
    
    // Read data
    handle=sdata_open(input_information->input_files[0]);
    
    if(handle==NULL){
        printf("ERROR: Could not open file: ");
        perror(input_information->input_files[0]);
        exit(-1);
    }
        
    // Number of segment that will be processed
    input_information->dimensions.segment_ncols=1;
    input_information->dimensions.segment_nrows=1; 
    input_information->dimensions.segment_ntimes=1;
    
    // Number of pixels that will be processed
    input_information->dimensions.npixel_estimated=handle->header.npixels;
    input_information->dimensions.tile_nt=handle->header.nrecords;
    input_information->dimensions.tile_nx=handle->header.nx;
    input_information->dimensions.tile_ny=handle->header.ny;
    
    // Used files
    input_information->nused_files=1;
    input_information->used_files=(char **)trackmem_malloc(sizeof(char *));
    assert(input_information->used_files!=NULL);
    input_information->used_files[0]=(char *)trackmem_malloc(sizeof(char) * (strlen(input_information->input_files[0]) +1));
    assert(input_information->used_files[0]!=NULL);
    strcpy(input_information->used_files[0],input_information->input_files[0]);
    
     
    // Return iterator
    return 0;    
}

grasp_settings_parameter_array *grasp_input_driver_settings_sdata(grasp_settings *settings){
    int i;
    
    // Static definition of a dictionary
    yamlsettings_parameter parameters[]= {
         /*                 name                                                                                                                                                                                   , memory direction                               , counter memory direction              , func to set variable type (length of value)             , initial value                               ,  number of elements          , allow_array             , parameter description                                                                                                                                                                                                                                                                                                                ,                 validator 1                 ,                   validator 2                         ,                    validator 3                            , {input function, output function, assigned}   */
        {"debug"                                                                                                                                                                       , &settings->input.driver.sdata.debug             , NULL                                  , YS_DATA_TYPE_BOOLEAN                                    , {ys_d_all,{"false"}}                                     , {0,{}}                       , YS_PARAM_SCALAR         , "Print debug information from sdata reader subsystem"                                                                                                                                                                                                                                                                                , {                                                                                                                                                              },YAMLSETTINGS_END_VAR  },                  
    };
    grasp_settings_parameter_array *result;
    
    result = (grasp_settings_parameter_array *) trackmem_malloc(sizeof (grasp_settings_parameter_array));
    assert(result!=NULL);
    
    result->nparameters=sizeof(parameters)/sizeof(yamlsettings_parameter);
    
    result->parameters = (yamlsettings_parameter *) trackmem_malloc(sizeof (yamlsettings_parameter)*result->nparameters);
    assert(result->parameters!=NULL);
    
    for(i=0;i<sizeof(parameters)/sizeof(yamlsettings_parameter);i++){
       yamlsettings_copy_parameter(&(parameters[i]),&(result->parameters[i]));
    }
    
    return result;   
}

