#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include "yamlsettings_data_types.h"
#include <grasp/utils.h> // for strtolower, safe_atoi and safe_atofloat
#include <glib.h>
#include "yamlsettings_error.h"

int yamlsettings_data_type_boolean_set(yamlsettings_status *status,void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    bool *d;
    GString *error;
    
    d = (bool *) mem_pos;

    c = strtolower(data);

    if (strcmp(c, "true") == 0 || strcmp(c, "yes") == 0 || strcmp(c, "1") == 0 || strcmp(c, "t") == 0 || strcmp(c, "") == 0) {
        d[position] = true;
    } else if (strcmp(c, "false") == 0 || strcmp(c, "no") == 0 || strcmp(c, "0") == 0 || strcmp(c, "f") == 0) {
        d[position] = false;
    } else {  
        error=g_string_new(NULL);
        g_string_printf (error,
                 "Data can not be assigned. Boolean value unrecognised in parameter %s (read: %s)",
                 name, c);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);        
        free(c);
        g_string_free(error,true);
        return -1;
    }
    free(c);
    
    return 0;
}
char *yamlsettings_data_type_boolean_get(void *mem_pos,int position, int maxlength){
    char *result;
    
    result = (char *) malloc(sizeof (char)*YAMLSETTINGS_VAR_VALUE_MAX_LENGTH);
    assert(result!=NULL);
    strcpy(result,"");
    if(mem_pos!=NULL){
        if (((bool *) mem_pos)[position] == true) {
            strcpy(result, "true");
        } else {
            strcpy(result, "false");
        }    
    }
    
    return result;
}

int yamlsettings_data_type_int_set(yamlsettings_status *status, void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    int *a;
    int result, v;
    GString *error;
    
    v=safe_atoi(data, &result);
    if(v==0){
        a = (int *) mem_pos;
        a[position] = result;
        return 0;
    }else{
        error=g_string_new(NULL);
        g_string_printf (error,"%s has a invalid integer (read: %s)", name, data);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error,true);
        return -1;        
    }
}
char *yamlsettings_data_type_int_get(void *mem_pos,  int position, int maxlength){
    char *result;

    result = (char *) malloc(sizeof (char)*YAMLSETTINGS_VAR_VALUE_MAX_LENGTH);
    assert(result!=NULL);
    strcpy(result,"");
    if(mem_pos!=NULL){
        sprintf(result, "%d", ((int *) mem_pos)[position]);
    }
    
    return result;
}

int yamlsettings_data_type_float_set(yamlsettings_status *status,void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    float *a;
    float result;
    int v;
    GString *error;
    
    v=safe_atofloat(data, &result);
    if(v==0){
        a = (float *) mem_pos;
        a[position] = result;
        return 0;
    }else{
        error=g_string_new(NULL);
        g_string_printf (error,"%s has a invalid float (read: %s)", name, data);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error,true);
        return -1;        
    }    
}
char *yamlsettings_data_type_float_get(void *mem_pos,  int position, int maxlength){    
    char *result;

    result = (char *) malloc(sizeof (char)*YAMLSETTINGS_VAR_VALUE_MAX_LENGTH);
    assert(result!=NULL);
    strcpy(result,"");
    if(mem_pos!=NULL){
        sprintf(result, "%e", ((float *) mem_pos)[position]);
    }
    
    return result;
}

int yamlsettings_data_type_string_set(yamlsettings_status *status,void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max;
    GString *error;
    
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0];
    max = strlen(data);
    if (max < maxlength) {
        // Copy the string in fortran format
        strcpy(c,data);
    } else {
        error=g_string_new(NULL);
        g_string_printf (error,"String too long in parameter %s", name);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error, true);
        return -1;
    }

    return 0;
}
char *yamlsettings_data_type_string_get(void *mem_pos,  int position, int maxlength){
    char *result;
    
    result = (char *) malloc(sizeof (char)*maxlength+1);
    assert(result!=NULL);
    strcpy(result,"");
        
    if(mem_pos!=NULL){
        strcpy(result, mem_pos);
    }
    
    return result;
}

int yamlsettings_data_type_file_path_set(yamlsettings_status *status,void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max;
    char *base_path;
    GString *error;
    
    base_path=pathoffile(settings_file);
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0]; //index2D(position,0,maxlength)
    max = strlen(data);
    
    if(isabsolute(data)==false){
        max=max+strlen(base_path);
    }

    if (max < maxlength) {
        // Copy the string in fortran format     
        if(isabsolute(data)==true){
            strcpy(c,"");
        }else{
            sprintf(c,"%s",base_path);
        }
        strcat(c,data);
    } else {
        error=g_string_new(NULL);
        g_string_printf (error,"String too long in parameter %s", name);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error, true);
        free(base_path);
        return -1;
    }
    
    free(base_path);
    return 0;
}
char *yamlsettings_data_type_file_path_get(void *mem_pos,  int position, int maxlength){
    char *result;
    
    result = (char *) malloc(sizeof (char)*maxlength+1);
    assert(result!=NULL);
    strcpy(result,"");
        
    if(mem_pos!=NULL){
        strcpy(result, mem_pos);
    }
    
    return result;
}

int yamlsettings_data_type_folder_path_set(yamlsettings_status *status,void *mem_pos, int nelements,  const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int max;
    char *base_path;
    bool add_slash=false;
    GString *error;
    
    base_path=pathoffile(settings_file);
    c = (char *) mem_pos;
    c = &c[(position*maxlength)+0];
    max = strlen(data);
    
    if(data[max-1]!='/'){
        max+=1;
        add_slash=true;
    }
    
    if(isabsolute(data)==false){
        max=max+strlen(base_path);
    }

    if (max < maxlength) {
        // Copy the string in fortran format     
        if(isabsolute(data)==true){
            strcpy(c,"");
        }else{
            strcpy(c,base_path);
        }
        strcat(c,data);
        if(add_slash==true){
            strcat(c,"/");
        }
    } else {
        error=g_string_new(NULL);
        g_string_printf (error, "String too long in parameter %s", name);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error,true);
        free(base_path);
        return -1;
    }
    
    free(base_path);
    return 0;
}
char *yamlsettings_data_type_folder_path_get(void *mem_pos,  int position, int maxlength){
    char *result;
    
    result = (char *) malloc(sizeof (char)*maxlength+1);
    assert(result!=NULL);
    strcpy(result,"");
        
    if(mem_pos!=NULL){
        strcpy(result, mem_pos);
    }
    
    return result;
}


int yamlsettings_data_type_generic_enumeration_save_integer(void *mem_pos,int position, int value){
    int *a;
    
    a=(int *)mem_pos;
    a[position]=value;
    
    return 0;
}

int yamlsettings_data_type_generic_enumeration_retrieve_integer(void *mem_pos,int position){
    return ((int *) mem_pos)[position];
}

int yamlsettings_data_type_generic_enumeration_set(yamlsettings_data_type_enumeration_save save_function, yamlsettings_enumeration_definition *ged,yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    char *c;
    int i;
    GString *error;
    
    c = strtolower(data);
    
    for (i = 0; i < ged->nelements; i++) {
        if (strcmp(ged->elements[i].name, c) == 0 || strcmp(ged->elements[i].second_name, c) == 0) {
            save_function(mem_pos,position,ged->elements[i].value);
            free(c);
            return 0;
        }
    }
    error=g_string_new(NULL);
    g_string_printf (error, "Data can not be assigned. %s is unrecognised value in parameter %s (values accepted are:", c, name);
    for (i = 0; i < ged->nelements; i++) {
        if (i != 0) g_string_append(error, ",");
        g_string_append(error, " ");
        g_string_append(error, ged->elements[i].name);
    }
    g_string_append(error, ")");
    yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
    g_string_free(error,true);
    free(c);
    return -1;
}

char *yamlsettings_data_type_generic_enumeration_get(yamlsettings_data_type_enumeration_retrieve function_retrieve, yamlsettings_enumeration_definition *ged, void *mem_pos,int position, int maxlength){
    int i, value;
    char *result;
    if(mem_pos!=NULL){
        result = (char *) malloc(sizeof (char)*YAMLSETTINGS_VAR_VALUE_MAX_LENGTH);
        assert(result!=NULL);
        strcpy(result,"");
        for (i = 0; i < ged->nelements; i++) {
            value=function_retrieve(mem_pos,position);
            if (ged->elements[i].value == value) {
                strcpy(result, ged->elements[i].name);
            }
        }
        if(strcmp(result,"")==0){
            strcpy(result,"Unknown data value");
        }
    }else{
        result = (char *) malloc(sizeof (char)*YAMLSETTINGS_VAR_VALUE_MAX_LENGTH*ged->nelements);
        assert(result!=NULL);
        strcpy(result,"");
        for (i = 0; i < ged->nelements; i++) {
            if(i!=0){
                strcat(result,", ");
            }
            strcat(result,ged->elements[i].name);
            
        }
    }
    
    return result; 
}

int yamlsettings_data_type_assign(yamlsettings_data_type_set function_type_set,yamlsettings_status *status, void *mem_pos, int nelements, const char *name, const char *data, int position, int maxlength, int file_index, const char *settings_file){
    GString *error;
    
    // This condition always should be true because it is already checked in the previous function
    if(position<nelements){ // If the position is bigger than the max number of elements the value must not be assigned
        return function_type_set(status, mem_pos, nelements,  name, data, position, maxlength, file_index, settings_file);
    }else{
        error=g_string_new(NULL);
        g_string_printf (error, "%s: Array contains more elements than size (maximum: %d, element: %d)", name,nelements,position);
        yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
        g_string_free(error,true);
        return -10;
    }
}
char *yamlsettings_data_type_get_string(yamlsettings_data_type_get function_type_get, void *mem_pos, int position, int maxlength) {
    return function_type_get( mem_pos, position, maxlength);
}

