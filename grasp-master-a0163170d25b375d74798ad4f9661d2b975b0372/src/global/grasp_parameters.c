/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "grasp_parameters.h"
#include "../settings/grasp_settings_data_types.h"

void grasp_parameters_initialize(float iguess[_KIMAGE][_KPARS]){
    int iimage, ipar;

    for (iimage = 0; iimage < _KIMAGE; iimage++) {
        for (ipar = 0; ipar < _KPARS; ipar++) {
            iguess[iimage][ipar] = -999.0;
        }
    }
}

int grasp_parameters_get_characteristic_index_by_parameter_number(par_number_NDIM *dimensions, int parameter_number){
    int i;
    
    // Look for characteristic index of the parameter
    for (i = dimensions->n1-1; i >= 0; i--) { 
        if(parameter_number>=dimensions->ISTARSING[i][0]-1){
            break;
        }
    }
    
    if(i<0 || i>=_KIDIM1){
        return -1;
    }
    
    return i; // index found
}

int grasp_parameters_get_characteristic_code_by_parameter_number(par_number_NDIM *dimensions, int parameter_number){
    int i;
    
    i=grasp_parameters_get_characteristic_index_by_parameter_number(dimensions, parameter_number);
    
    if(i<0 || i>=_KIDIM1){
        return -1;
    }
    
    return dimensions->par_type[i]; // Code of the characteristic found  
}


int grasp_parameters_get_mode_by_parameter_number(par_number_NDIM *dimensions, int parameter_number){
    int i, code;
    int from, to;
    
    code=grasp_parameters_get_characteristic_index_by_parameter_number(dimensions,parameter_number);
    
    // Look for characteristic index of the parameter
    for (i = 0; i < _KIDIM2-1; i++) { 
        from=dimensions->ISTARSING[code][i];
        to=dimensions->ISTARSING[code][i+1];
        if(to==0) to=9999;
        if(parameter_number+1 >= from &&
           parameter_number+1 < to ){
            break;
        }
    }
    
    if(i<0 || i>=_KIDIM2){
        return -1;
    }
    
    return i+1; // Mode of the characteristic found  
}

int grasp_parameters_get_position_by_parameter_number(par_number_NDIM *dimensions, int parameter_number){
    int code, mode;

    code=grasp_parameters_get_characteristic_index_by_parameter_number(dimensions,parameter_number);
    mode=grasp_parameters_get_mode_by_parameter_number(dimensions, parameter_number);
    // Look for characteristic index of the parameter
    
    return parameter_number + 1 - dimensions->ISTARSING[code][mode-1];
}

void grasp_parameters_get_characteristic_type_longname_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, char *characteristic_name, int size_characteristic_name){
    int i,code;

    code=grasp_parameters_get_characteristic_code_by_parameter_number(dimensions,parameter_number);                       
    
    assert( code > -1 );
    
    // Looking for the name of the characteristic
    for (i = 0; i < grasp_data_type_characteristic.nelements; i++) {
        if(grasp_data_type_characteristic.elements[i].value==code){
            assert(strlen(grasp_data_type_characteristic.elements[i].name)<size_characteristic_name);
            sprintf(characteristic_name,"%s",grasp_data_type_characteristic.elements[i].name); 
            break;
        }
    }         
}

void grasp_parameters_get_characteristic_type_shortname_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, char *characteristic_name, int size_characteristic_name){
    int i,code;
                           
    code=grasp_parameters_get_characteristic_code_by_parameter_number(dimensions,parameter_number);                       
    
    assert( code > -1 );
    
    // Looking for the name of the characteristic
    for (i = 0; i < grasp_data_type_characteristic.nelements; i++) {
        if(grasp_data_type_characteristic.elements[i].value==code){
            assert(strlen(grasp_data_type_characteristic.elements[i].name)<size_characteristic_name);
            sprintf(characteristic_name,"%s",grasp_data_type_characteristic.elements[i].second_name); 
            break;
        }
    }         
}


void grasp_parameters_get_characteric_type_pretty_name_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, bool longname, float wavelenghts[_KW], int wavelenghts_involved[_KPARS], char *characteristic_name, int size_characteristic_name){
    char *tmp;
    int code;
    int mode,modes;
    int position;

    // We get the name and contact it
    if (longname==true){
        grasp_parameters_get_characteristic_type_longname_by_parameter_number(dimensions,parameter_number,characteristic_name,size_characteristic_name);
    }else{
        grasp_parameters_get_characteristic_type_shortname_by_parameter_number(dimensions,parameter_number,characteristic_name,size_characteristic_name);
    }

    // We check if it is wavelength dependent
    code=grasp_parameters_get_characteristic_code_by_parameter_number(dimensions,parameter_number);
    if(wavelenghts_involved[parameter_number]>0){
        tmp = (char *) malloc(sizeof (char)*size_characteristic_name);
        if(longname==true){
            strcat(characteristic_name,"_");
        }
        sprintf(tmp,"%d", (int)(wavelenghts[wavelenghts_involved[parameter_number]-1]*1000)); // WARNING: Here we truncate the wavelength to have a integer part. I could be problematic
        strcat(characteristic_name,tmp);
        free(tmp);
    }else{
        mode=grasp_parameters_get_mode_by_parameter_number(dimensions,parameter_number);
        // Number inside the mode.
        position=grasp_parameters_get_position_by_parameter_number(dimensions,parameter_number);
        if(grasp_parameters_number_of_elements_of_parameter(dimensions, code, mode)>1){ // More than one elements, otherwise we don't need to specify
            tmp = (char *) malloc(sizeof (char)*size_characteristic_name);
            if(longname==true){
                strcat(characteristic_name,"_");
            }
            sprintf(tmp,"%d", position+1);
            strcat(characteristic_name,tmp);
            free(tmp);
        }
    }
    
    // If the characteristic has more than one mode I concat it at the begining of the value
    // We check how many modes has the parameter
    modes=grasp_parameters_number_of_modes_of_parameter(dimensions, code);
    assert(modes>-1);
    if(modes==1){ 
        // We do nothing
    }else{
        // We extract mode number of current parameter and concate it at the begining of the string
        tmp = (char *) malloc(sizeof (char)*size_characteristic_name);
        sprintf(tmp,"_%d", grasp_parameters_get_mode_by_parameter_number(dimensions,parameter_number));
        strcat(characteristic_name,tmp);
        free(tmp);
    }    
}

int grasp_parameters_index_of_parameter_type(par_number_NDIM *dimensions, int characteristic_type){
    int i;
    
    for (i = 0; i < dimensions->n1; i++) {
        if(dimensions->par_type[i]==characteristic_type){
            return i;
        }
    }
    
    return -1;    
}

int grasp_parameters_index_of_parameter_type_by_kind_of_parameter(par_number_NDIM *dimensions, int begin_characteristic_type, int end_characteristic_type){
    int i;
    
    for (i = 0; i < dimensions->n1; i++) {
        if(dimensions->par_type[i]>begin_characteristic_type && dimensions->par_type[i]<end_characteristic_type){
            return i;
        }
    }
    
    return -1;    
}

int grasp_parameters_characteristic_code_present_of_kind_of_parameter(par_number_NDIM *dimensions, int begin_characteristic_type, int end_characteristic_type){
    int i;
    
    for (i = 0; i < dimensions->n1; i++) {
        if(dimensions->par_type[i]>begin_characteristic_type && dimensions->par_type[i]<end_characteristic_type){
            return dimensions->par_type[i];
        }
    }
    
    return -1;    
}


bool grasp_parameters_has_parameter_type(par_number_NDIM *dimensions, int characteristic_type){
    if(grasp_parameters_index_of_parameter_type(dimensions, characteristic_type)<0){
        return false;
    }else{
        return true;
    }
}

int grasp_parameters_set_value(par_number_NDIM *dimensions, float iguess[_KPARS], int characteristic_type, int mode , int pos, float value){
    int n1;
    int index;
    
    n1=grasp_parameters_index_of_parameter_type(dimensions,characteristic_type);
    
    assert(n1<_KIDIM1 && mode<=_KIDIM2 && pos<dimensions->n3[n1][mode-1]);

    index=dimensions->ISTARSING[n1][mode-1]+pos-1;
            
    if(value <= 0){
        return -1;
    }
    if(index < 0 || index >= _KPARS){
        return -2;
    }

    iguess[index]=value;
    
    return index; 
}

int grasp_parameters_number_of_modes_of_parameter(par_number_NDIM *dimensions, int characteristic_type){
    int n1;
    
    n1=grasp_parameters_index_of_parameter_type(dimensions,characteristic_type);
    
    assert(n1<_KIDIM1); 
    
    if(n1<0){
        return -1;
    }

    return dimensions->n2[n1];

}

int grasp_parameters_number_of_elements_of_parameter(par_number_NDIM *dimensions, int characteristic_type, int mode){
    int n1;
    
    n1=grasp_parameters_index_of_parameter_type(dimensions,characteristic_type);
    
    assert(n1<_KIDIM1); 
    
    if(n1<0){
        return -1;
    }

    if(mode>dimensions->n2[n1]){
        return -2;
    }
    
    return dimensions->n3[n1][mode-1];
}

float grasp_parameters_output_get_value(par_number_NDIM *dimensions, const float parameters[_KPARS], int characteristic_type, int mode , int pos){
    int index;
    
    index=grasp_parameters_get_position(dimensions, characteristic_type, mode, pos);
    
    if(index<0){
        return index;
    }
    
    assert(index<_KPARS);
    
    return parameters[index];
}

int grasp_parameters_get_position(par_number_NDIM *dimensions, int characteristic_type, int mode , int pos){
    int n1;
    int index;
    
    n1=grasp_parameters_index_of_parameter_type(dimensions,characteristic_type);
    
    if(n1<0){
        return -1;
    }
    
    if(grasp_parameters_number_of_modes_of_parameter(dimensions, characteristic_type)<mode) {
        return -2;
    }
    
    assert(mode<=_KIDIM2 && pos<dimensions->n3[n1][mode-1]);

    index=dimensions->ISTARSING[n1][mode-1]+pos-1;
      
    return index;
}

float grasp_parameters_iguess_get_value(par_number_NDIM *dimensions, float APSING[_KPARS], float iguess[_KPARS], float APSMIN[_KPARS], float APSMAX[_KPARS],  int characteristic_type, int mode , int pos){
    int n1;
    int index;
    
    n1=grasp_parameters_index_of_parameter_type(dimensions,characteristic_type);
    
    if(n1<0){
        return -1.0;
    }
    
    if(grasp_parameters_number_of_modes_of_parameter(dimensions, characteristic_type)>mode){
        return -2.0;
    }
    
    assert(mode<=_KIDIM2 && pos<dimensions->n3[n1][mode-1]);

    index=dimensions->ISTARSING[n1][mode-1]+pos-1;
            
    assert(index<_KPARS);
    
    if(iguess[index]!=-999){
        if(iguess[index]<APSMIN[index]){
            return APSMAX[index];
        }else{
            if(iguess[index]>APSMAX[index]){
                return APSMAX[index];
            }else{
                return iguess[index];
            }
        }
    }else{
        return APSING[index];
    }
}
