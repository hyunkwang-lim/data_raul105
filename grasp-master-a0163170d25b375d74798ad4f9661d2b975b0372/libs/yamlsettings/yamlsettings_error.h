/* 
 * File:   yamlsettings_error.h
 * Author: david
 *
 * Created on 24 de marzo de 2014, 10:11
 */

#ifndef YAMLSETTINGS_ERROR_H
#define	YAMLSETTINGS_ERROR_H

#ifdef	__cplusplus
extern "C" {
#endif
    
#include <glib.h>
    
    typedef enum  { YS_FATAL_ERROR , YS_ERROR, YS_WARNING, YS_INFO } yamlsettings_error_type;
    
    typedef struct yamlsettings_parse_error_{
        yamlsettings_error_type type;
        char *message;
        int file_index;
    }yamlsettings_parse_error;
    
    typedef struct yamlsettings_validation_error_{
        char *message;
        int file_index;
    }yamlsettings_validation_error;            
    
    typedef struct yamlsettings_status_{
        GPtrArray *parse_errors;
        GPtrArray *validation_errors;        
    } yamlsettings_status;

    // Add a parse error to yamlsettings status element
    void yamlsettings_error_add_parse_error(yamlsettings_status *status, char *message, int file_index, yamlsettings_error_type type);
    
    // Add a validator error to yamlsettings status element
    void yamlsettings_error_add_validation_error(yamlsettings_status *status, int file_index, char *message);    
    
    // Return number of error type in a status structure 
    int yamlsettings_error_number_parse_error(yamlsettings_status *status, yamlsettings_error_type searching_error);
    
    // Return number of validation errors
    int yamlsettings_error_number_validation_error(yamlsettings_status *status);
    
    // Free memory from a yamlsettings_parse_error structure
    void yaml_settings_status_parse_error_destroy(yamlsettings_parse_error *x);
    
    // Free memory from a yamlsettings_validation_error structure
    void yaml_settings_status_validation_error_destroy(yamlsettings_validation_error *x);
    
    // Free memory form yamlsettings_status structure
    void yamlsettings_error_status_destroy(yamlsettings_status *status);
    

#ifdef	__cplusplus
}
#endif

#endif	/* YAMLSETTINGS_ERROR_H */

