#include "yamlsettings_defaults.h"
#include "yamlsettings_dictionary.h"


int ys_d_all (int param_index, yamlsettings_parameter *parameter, yamlsettings_status *status){
    int i;
    
    for(i=0; i<yamlsettings_parameter_number_of_elements(parameter); i++){ // Initilize values. Loop for initilize arrays
        if(yamlsettings_data_type_assign(parameter->data_type_set, status, parameter->mem_pos,yamlsettings_parameter_number_of_elements(parameter),parameter->name,  parameter->default_value.arguments[0], i, parameter->maxlength, YAMLSETTINGS_FILE_UNSET, parameter->settings_file[0])<0){    
            printf("ERROR: Invalid initial value of parameter %s\n", parameter->name);
            abort();
        }
    }    
    
    return 0;
}

int ys_d_array (int param_index, yamlsettings_parameter *parameter, yamlsettings_status *status){
    int ipars, npars, position;
    char arg[YAMLSETTINGS_DEFAULTS_MAX_ARGUMENT_SIZE],*colon, *value;
    
    ys_d_all(param_index,parameter,status);
    
    npars=atoi(parameter->default_value.arguments[1]);
    
    for (ipars = 0; ipars < npars; ipars++) {
        // Get a copy of the argument
        strcpy(arg,parameter->default_value.arguments[2+ipars]);
        
        // Position of the colon symbol
        colon=strchr(arg,':');
        // The string "value" starts in next element, just after colon
        value=colon+1;
        // Replacing colon by EOL  we transform arg into position
        *colon='\0';
        // Just casting it to INT
        position=atoi(arg);
        // We save the element
        if(yamlsettings_data_type_assign(parameter->data_type_set, status, parameter->mem_pos,yamlsettings_parameter_number_of_elements(parameter),parameter->name,  value, position, parameter->maxlength, YAMLSETTINGS_FILE_UNSET, parameter->settings_file[0])<0){    
            printf("ERROR: Invalid initial value of parameter %s\n", parameter->name);
            abort();
        }
    }
    
    parameter->counter_mem_pos[0]=npars;
    
    return 0;
}
