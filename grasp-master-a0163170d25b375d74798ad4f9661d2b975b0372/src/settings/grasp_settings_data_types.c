/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>
#include "yamlsettings/yamlsettings_data_types.h"
#include "grasp_settings_data_types.h"
#include <grasp/utils.h>
#include "../input/grasp_input_settings.h"
#include "../global/grasp_retrieval_meas_type.h"
#include "../global/grasp_retrieval_characteristic_type.h"

int grasp_data_type_fortran_string_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max,i;
    char error[512];
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0];
    max = strlen(data);
    if (max <= maxlength) {
        // Copy the string in fortran format
        for (i = 0; i < max; i++) {
            c[i] = data[i];
        }
        for (i = max; i < maxlength; i++) {
            c[i] = ' ';
        }
    } else {
        sprintf(error,"String too long in parameter %s\n", name);
        yamlsettings_error_add_parse_error(status, error, file_index, YS_ERROR); 
        return -1;
    }

    return 0;
}
char *grasp_data_type_fortran_string_get(void *mem_pos,  int position, int maxlength){
    int nvalues;
    char *c,*result;
    
    result = (char *) trackmem_malloc(sizeof (char)*maxlength+1);
    assert( result!= NULL);
    strcpy(result,"");
    
    if(mem_pos!=NULL){
        nvalues = strfortranlen(mem_pos, maxlength);
        c = (char *) trackmem_malloc(sizeof (char)*nvalues + 1);
        assert( c!= NULL);
        strncpy(c, mem_pos, nvalues);
        c[nvalues] = '\0';
        strcpy(result, c);
        trackmem_free(c);  
    }
    
    return result;
}

int grasp_data_type_fortran_file_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max_data,i,max_base=0;
    char *base_path;
    char error[512];
    
    base_path=pathoffile(settings_file);
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0];
    max_data = strlen(data);
    
    if(isabsolute(data)==false){
        max_base=strlen(base_path);
    }

    if (max_base+max_data <= maxlength) {
        // Copy the string in fortran format
        for (i = 0; i < max_base; i++) { // Copy base if its necessary (if not max_base is equal to 0)
            c[i] = base_path[i];
        }
        for (i = 0; i < max_data; i++) { // Copy data
            c[max_base+i] = data[i];
        }
        for (i = 0; i < maxlength-(max_base+max_data); i++) { // Filling with zeros
            c[max_base+max_data+i] = ' ';
        }
        
    } else {
        printf("String too long in parameter %s\n", name);
        yamlsettings_error_add_parse_error(status, error, file_index, YS_ERROR);
        trackmem_free(base_path);
        return -1;
    }

    trackmem_free(base_path);
    return 0;
}
char *grasp_data_type_fortran_file_path_get(void *mem_pos,  int position, int maxlength){
    int nvalues;
    char *c,*result;
    
    result = (char *) trackmem_malloc(sizeof (char)*maxlength+1);
    assert( result!= NULL);
    strcpy(result,"");
    
    if(mem_pos!=NULL){
        nvalues = strfortranlen(mem_pos, maxlength);
        c = (char *) trackmem_malloc(sizeof (char)*nvalues + 1);
        assert( c!= NULL);
        strncpy(c, mem_pos, nvalues);
        c[nvalues] = '\0';
        strcpy(result, c);
        trackmem_free(c);  
    }
    
    return result;
}

int grasp_data_type_fortran_folder_path_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max_data,i,max_base=0;
    char *base_path;
    bool add_slash=false;
    char error[512];
    
    base_path=pathoffile(settings_file);
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0];
    max_data = strlen(data);
    
    if(data[max_data-1]!='/'){
        max_data+=1;
        add_slash=true;
    }
    
    if(isabsolute(data)==false){
        max_base=strlen(base_path);
    }

    if (max_base+max_data <= maxlength) {
        // Copy the string in fortran format
        for (i = 0; i < max_base; i++) { // Copy base if its necessary (if not max_base is equal to 0)
            c[i] = base_path[i];
        }
        for (i = 0; i < max_data; i++) { // Copy data
            c[max_base+i] = data[i];
        }
        for (i = 0; i < maxlength-(max_base+max_data); i++) { // Filling with zeros
            c[max_base+max_data+i] = ' ';
        }
        if(add_slash==true){
            c[max_base+max_data-1]='/';
        }        
    } else {
        printf("String too long in parameter %s\n", name);
        yamlsettings_error_add_parse_error(status, error, file_index, YS_ERROR);
        trackmem_free(base_path);
        return -1;
    }

    trackmem_free(base_path);
    return 0;
}
char *grasp_data_type_fortran_folder_path_get(void *mem_pos,  int position, int maxlength){
    int nvalues;
    char *c,*result;
    
    result = (char *) trackmem_malloc(sizeof (char)*maxlength+1);
    assert( result!= NULL);
    strcpy(result,"");
    
    if(mem_pos!=NULL){
        nvalues = strfortranlen(mem_pos, maxlength);
        c = (char *) trackmem_malloc(sizeof (char)*nvalues + 1);
        assert( c!= NULL);
        strncpy(c, mem_pos, nvalues);
        c[nvalues] = '\0';
        strcpy(result, c);
        trackmem_free(c);  
    }
    
    return result;
}

int grasp_data_type_stream_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c, *d;
    
    c = strtolower(data);
    
    // calculating position
    d = (char *) mem_pos;
    d = &d[(position*maxlength)+0];
    
    if (strcmp(c, "screen") == 0 || strcmp(c, "stdout") == 0 || strcmp(c, "true") == 0|| strcmp(c, "t") == 0 || strcmp(c, "1") == 0 || strcmp(c, "") == 0) {
        assert(maxlength>7);
        strcpy(d,"screen");
        trackmem_free(c);
        return 0;
    } else if (strcmp(c, "none") == 0 || strcmp(c, "null") == 0|| strcmp(c, "false") == 0|| strcmp(c, "f") == 0 || strcmp(c, "0") == 0) {
        assert(maxlength>5);
        strcpy(d,"none");
        trackmem_free(c);
        return 0;
    }else{
        if((strlen(c)>5 && c[0]=='{' && c[1]=='p' && c[2]=='w' && c[3]=='d' && c[4]=='}') ||  (strlen(c)>5 && c[0]=='{' && c[1]=='y' && c[2]=='m' && c[3]=='l' && c[4]=='}')){
            trackmem_free(c);
            return yamlsettings_data_type_string_set(status, mem_pos, nelements,  name, data, position, maxlength, file_index, settings_file);
        }else{
            trackmem_free(c);
            return yamlsettings_data_type_file_path_set(status, mem_pos, nelements,  name, data, position, maxlength, file_index, settings_file);
        }
    }
}
char *grasp_data_type_stream_get(void *mem_pos,  int position, int maxlength){
    return yamlsettings_data_type_file_path_get(mem_pos, position, maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_minimization={"minimization", 2,
        {
            {"absolute", "abs", 0},
            {"logarithm", "log", 1}
        }};
int grasp_data_type_minimization_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_minimization, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_minimization_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_minimization,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_ipplane={"ipplane", 2,
        {
            {"principal_plane", "pp", 0},
            {"meridian", "mer", 1}
        }};
int grasp_data_type_ipplane_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_ipplane, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_ipplane_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_ipplane,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_pfitting={"polarization_fitting", 5,
        {
            {"absolute_polarization_components", "abs_pol", 1},
            {"relative_polarization_components", "rel_pol", 2},
            {"polarized_reflectance", "pol_ref", 3},
            {"linear_polarization", "lin_pol", 4},
            {"relative_linear_polarization", "rel_lin_pol", 5}
        }};
int grasp_data_type_pfitting_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_pfitting, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_pfitting_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_pfitting,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_ifitting={"irradiance_fitting", 2,
        {
            {"radiances", "rad", 1},
            {"relative_radiances", "rel_rad", 2}
        }};
int grasp_data_type_ifitting_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_ifitting, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_ifitting_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_ifitting,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_bin={"bin", 2,
        {
            {"absolute", "abs", 1},
            {"logarithm", "log", -1}
        }};
int grasp_data_type_bin_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_bin, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_bin_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_bin,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_invsing={"inversion_regime", 3,
        {
            {"single_pixel", "single", 0},
            {"multi_pixel_followed_by_single_pixel", "mp_then_sp", 1},
            {"multi_pixel", "multi", 2}
        }};
int grasp_data_type_invsing_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_invsing, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_invsing_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_invsing,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_imsc={"imsc", 3,
        {
            {"multiple_scattering", "multiple", 0},
            {"single_scattering", "single", 1},
            {"derivatives","derivatives",2}
        }};
int grasp_data_type_imsc_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_imsc, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_imsc_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_imsc,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_error={"error", 2,
        {
            {"absolute", "abs", 1},
            {"relative", "rel", 0}
        }};
int grasp_data_type_error_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_error, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_error_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_error,mem_pos,position,maxlength);
}

#ifdef WARN_DRY
#warning "__MEAS_TYPE__ duplicated"
#endif   

yamlsettings_enumeration_definition grasp_data_type_measuretypes={"measuretypes", 15,
        {
            {"tod", "total_optical_depth",   MEAS_TYPE_TOD},
            {"aod", "aerosol_optical_depth", MEAS_TYPE_AOD},
            {"p11", "p11",                   MEAS_TYPE_P11},
            {"p12", "p12",                   MEAS_TYPE_P12},
            {"p22", "p22",                   MEAS_TYPE_P22},
            {"p33", "p33",                   MEAS_TYPE_P33},
            {"p34", "p34",                   MEAS_TYPE_P34},
            {"p44", "p44",                   MEAS_TYPE_P44},
            {"ls", "ls",                     MEAS_TYPE_LS},
            {"dp", "dp",                     MEAS_TYPE_DP},
            {"rl", "r",                      MEAS_TYPE_RL},
            {"i", "i",                       MEAS_TYPE_I},
            {"q", "q",                       MEAS_TYPE_Q},
            {"u", "u",                       MEAS_TYPE_U},
            {"p", "p",                       MEAS_TYPE_P}
            
        }};
int grasp_data_type_measuretypes_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_measuretypes, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_measuretypes_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_measuretypes,mem_pos,position,maxlength);
}

yamlsettings_enumeration_definition grasp_data_type_imq={"normal_system_solver", 3,
        {
            {"simple_linear_iterations", "sli", 1},
            {"singular_value_decomposition", "svd", 2},
            {"sparse_matrix_solver", "sms", 3}
        }};
int grasp_data_type_imq_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_imq, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_imq_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_imq,mem_pos,position,maxlength);
}

#ifdef WARN_DRY
#warning "__CHARACTERISTIC_TYPE__ duplicated"
#endif  

yamlsettings_enumeration_definition grasp_data_type_characteristic={"characteristic_type", 25,
        {
            {"size_distribution_triangle_bins",                       "SizeDistrTriangBin",    par_type_SD_TB},
            {"size_distribution_precalculated_lognormal",             "SizeDistrLogNormBin",   par_type_SD_LB},
            {"size_distribution_lognormal",                           "SizeDistrLogNorm",      par_type_SD_LN},
            {"real_part_of_refractive_index_spectral_dependent",      "RealRefIndSpect",       par_type_RERI_spect},
            {"real_part_of_refractive_index_constant",                "RealRefIndConst",       par_type_RERI_const},
            {"particle_component_volume_fractions_linear_mixture",    "VolFractLinearMix",     par_type_CXRI_nmix},
            {"particle_component_fractions_chemical_mixture",         "FractChemMix",          par_type_CXRI_chem},
            {"imaginary_part_of_refractive_index_spectral_dependent", "ImagRefIndSpect",       par_type_IMRI_spect},
            {"imaginary_part_of_refractive_index_constant",           "ImagRefIndConst",       par_type_IMRI_const},
            {"sphere_fraction",                                       "SphereFraction",        par_type_SHD_fsph},
            {"aspect_ratio_distribution",                             "AspectRatioDistrib",    par_type_SHD_distr},
            {"vertical_profile_parameter_height",                     "VertProfileHeight",     par_type_AVP_par_height},
            {"vertical_profile_normalized",                           "VertProfileNormalized", par_type_AVP_prof},
            {"aerosol_concentration",                                 "AerosolConcentration",  par_type_Cv},                    
            {"lidar_calibration_coefficient",                         "LidarCalibrCoeff",      par_type_CL},  
            {"vertical_profile_parameter_standard_deviation",         "VertProfileStdDev",     par_type_AVP_par_std},   
            {"surface_land_brdf_ross_li",                             "LandBRDFRossLi",        par_type_SURF1_land_Ross_Li_BRDF},
            {"surface_land_brdf_rpv",                                 "LandBRDFRPV",           par_type_SURF1_land_RPV_BRDF},
            {"surface_land_litvinov",                                 "LandBRMLitvinov",       par_type_SURF1_land_Litvinov},
            {"surface_land_litvinov_fast",                            "LandBRDFfastLitvinov",  par_type_SURF1_land_Litvinov_fast},
            {"surface_land_polarized_maignan_breon",                  "LandBPDFMaignanBreon",  par_type_SURF2_land_Maignan_Breon},
            {"surface_land_polarized_litvinov",                       "LandBPDFLitvinov",      par_type_SURF2_land_Litvinov},                      
            {"surface_water_cox_munk_iso",                            "WaterBRMCoxMunkIso",    par_type_SURF_water_Cox_Munk_iso},
            {"surface_water_cox_munk_ani",                            "WaterBRMCoxMunkAni",    par_type_SURF_water_Cox_Munk_ani},   
            {"surface_water_litvinov",                                "WaterBRMLitvinov",      par_type_SURF_water_Litvinov},
        }};

int grasp_data_type_characteristic_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {    
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_characteristic, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_characteristic_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_characteristic,mem_pos,position,maxlength);
}


yamlsettings_enumeration_definition grasp_data_type_surface={"surface", 2,
        {
            {"ocean", "o", 1},
            {"land", "l", 2}
        }};
int grasp_data_type_surface_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_surface, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_surface_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_surface,mem_pos,position,maxlength);
}


yamlsettings_enumeration_definition grasp_data_type_molecular_profile_vertical_type={"molecular_profile_vertical_type", 2,
        {
            {"exponential", "exp", 0},
            {"standard_atmosphere", "std_atm", 1},
        }};
int grasp_data_type_molecular_profile_vertical_type_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_molecular_profile_vertical_type, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_molecular_profile_vertical_type_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_molecular_profile_vertical_type,mem_pos,position,maxlength);
}


yamlsettings_enumeration_definition grasp_data_type_aerosol_profile_vertical_type={"aerosol_profile_vertical_type", 2,
        {
            {"exponential", "exp", 0},
            {"gaussian", "gauss", 1},
        }};
int grasp_data_type_aerosol_profile_vertical_type_set(yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file) {
    return yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_generic_enumeration_save_integer, &grasp_data_type_aerosol_profile_vertical_type, status, mem_pos,nelements,name,data,position,maxlength, file_index, settings_file);
}
char *grasp_data_type_aerosol_profile_vertical_type_get(void *mem_pos,int position, int maxlength){
    return yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_generic_enumeration_retrieve_integer, &grasp_data_type_aerosol_profile_vertical_type,mem_pos,position,maxlength);
}

