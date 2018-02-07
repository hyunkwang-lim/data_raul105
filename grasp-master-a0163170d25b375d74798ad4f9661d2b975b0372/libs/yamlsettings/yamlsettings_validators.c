#include <dirent.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <grasp/utils.h>
#include "yamlsettings_validators.h"
#include "yamlsettings_dictionary.h"



int yamlsettings_validator_int(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int *value,min,max,i;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    param=&(dictionary->parameters[param_index]);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        value=(int *)param->mem_pos;
        min=atoi(arguments[0]);
        max=atoi(arguments[1]);
        for(i=0;i<yamlsettings_parameter_number_of_elements(param);i++){
            if(!(value[i]>=min && value[i]<=max)){
                sprintf(error,"%s must be a number between %d and %d (value: %d)",param->name, min, max,value[i]);
                yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
                nerror--;
            }
        }
    }
    
    return nerror;
}

int yamlsettings_validator_double(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int i;
    double *value,min,max;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    param=&(dictionary->parameters[param_index]);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        value=(double *)param->mem_pos;
        min=atof(arguments[0]);
        max=atof(arguments[1]);
        for(i=0;i<yamlsettings_parameter_number_of_elements(param);i++){
            if(!(value[i]>=min && value[i]<=max)){
                sprintf(error,"%s must be a number between %f and %f (read: %f)",param->name, min, max, value[i]);
                yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
                nerror--;
            }
        }
    }
    
    return nerror;
}

int yamlsettings_validator_float(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int i;
    float *value,min,max;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        value=(float *)param->mem_pos;
        min=atof(arguments[0]);
        max=atof(arguments[1]);
        for(i=0;i<yamlsettings_parameter_number_of_elements(param);i++){
            if(!(value[i]>=min && value[i]<=max)){
                sprintf(error,"%s must be a number between %e and %e (read: %e)",param->name, min, max, value[i]);
                yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
                nerror--;
            }
        }
    }
    
    return nerror;
}

int yamlsettings_validator_mandatory (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    param=&(dictionary->parameters[param_index]);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==true){
        sprintf(error,"%s is required",param->name);
        yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
        nerror--;
    }
    
    return nerror;
}

int yamlsettings_validator_mandatory_if_parameter_true (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    yamlsettings_parameter *param;
    int index;
    bool reference_boolean;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    index=yamlsettings_dictionary_find_parameter_by_name(arguments[0],dictionary);
    assert(index>=0);
    reference_boolean= *((bool *)dictionary->parameters[index].mem_pos);
        
    if(reference_boolean==true){
        if(yamlsettings_parameter_is_unset(param_index, dictionary)==true){
            sprintf(error,"%s is required if %s is true",param->name, dictionary->parameters[index].name);
            yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
            nerror--;
        }
    }

    return nerror;
    
}

int yamlsettings_validator_mandatory_if_parameter_enum (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    yamlsettings_parameter *param;
    int index;
    char *value;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    index=yamlsettings_dictionary_find_parameter_by_name(arguments[0],dictionary);
    assert(index>=0);
    value= dictionary->parameters[index].data_type_get(dictionary->parameters[index].mem_pos,0,dictionary->parameters[index].maxlength);

    if(strcmp(value,arguments[1])==0){
        if(yamlsettings_parameter_is_unset(param_index, dictionary)==true){
            sprintf(error,"%s is required if %s is %s",param->name, dictionary->parameters[index].name,value);
            yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
            nerror--;
        }
    }   
    
    free(value);
    return nerror;
}

int yamlsettings_validator_mandatory_if_parameter_enum_not (int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    yamlsettings_parameter *param;
    int index;
    char *value;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    index=yamlsettings_dictionary_find_parameter_by_name(arguments[0],dictionary);
    assert(index>=0);
    value= dictionary->parameters[index].data_type_get(dictionary->parameters[index].mem_pos,0,dictionary->parameters[index].maxlength);

    if(strcmp(value,arguments[1])!=0){
        if(yamlsettings_parameter_is_unset(param_index, dictionary)==true){
            sprintf(error,"%s is required if %s is %s",param->name, dictionary->parameters[index].name,value);
            yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
            nerror--;
        }
    }
    
    free(value);
    return nerror;
       
}

int yamlsettings_validator_directory(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    char *path;
    DIR *dir;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    
    // Retrieve parameter
    param=&(dictionary->parameters[param_index]);    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        // Get pointer to first char
        path=(char *)param->mem_pos;

        dir= opendir(path);

        if (dir) { //Directory exists.
            closedir(dir);
        } else {
            sprintf(error,"%s must be an existing directory", param->name);
            yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
            nerror--;
        } 

        //free(path);
    }
    
    return nerror;
}

int yamlsettings_validator_nelements(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int *value,min,max;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        value=(int *)param->counter_mem_pos;
        min=atoi(arguments[0]);
        max=atoi(arguments[1]);
        if(!(*value>=min && *value<=max)){
            sprintf(error,"%s must be an array having from %d to %d elements.",param->name, min, max);
            yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
            nerror--;
        }
        
    }
    
    return nerror;    
}

int yamlsettings_validator_same_nelements(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int *value1,*value2,index2,i,indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    yamlsettings_parameter *param;
    char error[2048], *param_name1, *param_name2;  
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);    
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        value1=(int *)param->counter_mem_pos;
        index2=yamlsettings_dictionary_find_parameter_by_name(arguments[0],dictionary);
        assert(index2>=0);
        value2=dictionary->parameters[index2].counter_mem_pos;
        for (i = 0; i < yamlsettings_parameter_counter_number_of_elements(param); i++) {
            if(value1[i]!=value2[i]){
                indexesofND(i,indexes,yamlsettings_parameter_number_of_dimensions_in_string(param),param->dimensions.dimension_size);
                param_name1=yamlsettings_parameter_get_with_indexes(param->name, yamlsettings_parameter_number_of_dimensions_in_string(param), indexes);
                param_name2=yamlsettings_parameter_get_with_indexes(dictionary->parameters[index2].name, yamlsettings_parameter_number_of_dimensions_in_string(param), indexes);
                sprintf(error,"%s (%d elements) must be an array the same number of elements that %s (%d elements)."
                        ,param_name1,value1[i]/**(param->counter_mem_pos)*/, 
                        param_name2, value2[i] /**(dictionary->parameters[index2].counter_mem_pos)*/);
                yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
                nerror--;
                free(param_name1);
                free(param_name2);
            }
        }        
    }
    
    return nerror;    
}

int yamlsettings_validator_ascendant_int_array(int param_index,yamlsettings_dictionary_t *dictionary, char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE]){
    int *value,*counter,i;
    yamlsettings_parameter *param;
    char error[2048];
    int nerror=0;
    
    param=&(dictionary->parameters[param_index]);
    
    // This validator is only valid (for the moment) for validating one dimensional arrays
    assert(param->dimensions.ndimensions==1);
    
    if(yamlsettings_parameter_is_unset(param_index, dictionary)==false){
        counter=(int *)param->counter_mem_pos;
        value=(int *)param->mem_pos;
        for (i = 1; i < *counter; i++) {
            if(value[i-1]>value[i]){
                sprintf(error,"%s must be an ascendant array",param->name);
                yamlsettings_error_add_validation_error(&(dictionary->status), yamlsettings_parameter_file_index(param_index, dictionary), error); 
                nerror--;
            }
        }
    }
    
    return nerror;    
}

