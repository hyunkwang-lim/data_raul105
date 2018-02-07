#include <stdio.h>
#include <stdlib.h>
#include "yamlsettings_error.h"
#include <string.h>
#include <assert.h>



void yamlsettings_error_add_parse_error(yamlsettings_status *status, char *message, int file_index, yamlsettings_error_type type){
    yamlsettings_parse_error *error;
    char *error_message;
    
    error = (yamlsettings_parse_error *) malloc(sizeof (yamlsettings_parse_error));
    assert(error!=NULL);
    error_message = (char *) malloc(sizeof (char)*(strlen(message)+1));
    assert(error_message!=NULL);
    strcpy(error_message, message);
    error->message=error_message;
    error->file_index=file_index;
    error->type=type;
    g_ptr_array_add(status->parse_errors, error);       
}

void yamlsettings_error_add_validation_error(yamlsettings_status *status, int file_index, char *message){
    yamlsettings_validation_error *error;
    char *error_message;
    
    error = (yamlsettings_validation_error *) malloc(sizeof (yamlsettings_validation_error));
    assert(error!=NULL);
    error_message = (char *) malloc(sizeof (char)*(strlen(message)+1));
    assert(error_message!=NULL);
    strcpy(error_message, message);
    error->message=error_message;
    error->file_index=file_index;
    g_ptr_array_add(status->validation_errors, error);       
}

int yamlsettings_error_number_parse_error(yamlsettings_status *status, yamlsettings_error_type searching_error){
    int result = 0;
    int i;
    
    for (i = 0; i < status->parse_errors->len; i++) {
        if( ((yamlsettings_parse_error *)g_ptr_array_index(status->parse_errors,i))->type==searching_error ){
            result++;
        }
    }
    
    return result;
}

int yamlsettings_error_number_validation_error(yamlsettings_status *status){
    int result = 0;
    int i;
    
    for (i = 0; i < status->validation_errors->len; i++) {
        result++;
    }
    
    return result;
}

void yaml_settings_status_parse_error_destroy(yamlsettings_parse_error *x){
    free(x->message);
    free(x);    
}

void yaml_settings_status_validation_error_destroy(yamlsettings_validation_error *x){
    free(x->message);
    free(x);        
}

void yamlsettings_error_status_destroy(yamlsettings_status *status){
    int i;
    
    // For parse errors
    for (i = 0; i < status->parse_errors->len; i++) {
        yaml_settings_status_parse_error_destroy(((yamlsettings_parse_error *)g_ptr_array_index(status->parse_errors,i)));
    }
    g_ptr_array_free (status->parse_errors, FALSE);
    
    // For validation errors
    for (i = 0; i < status->validation_errors->len; i++) {
        yaml_settings_status_validation_error_destroy(((yamlsettings_validation_error *)g_ptr_array_index(status->validation_errors,i)));
    }   
    g_ptr_array_free (status->validation_errors, FALSE);
    
    
}
