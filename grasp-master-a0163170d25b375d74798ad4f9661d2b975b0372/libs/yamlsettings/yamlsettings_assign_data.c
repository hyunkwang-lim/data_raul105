#include <stdlib.h>
#include <stdio.h>
#include "yamlsettings_assign_data.h"
#include "yamlsettings_input_yaml.h"
#include "yamlsettings_data_types.h"
#include "yamlsettings_dictionary.h"
#include "yamlsettings_error.h"
#include <grasp/utils.h>
#include <string.h>
#include <glib/gstring.h>

gboolean yamlsettings_assign_node(GNode *node, gpointer data) {
    yamlsettings_node_element *ne;
    char *complete_name_normalized;
    char *complete_name;
    int option;
    yamlsettings_dictionary_t *dictionary;
    char *unset; // pointer to YAMLSETTINGS_FILE_UNSET file name
    GString *error;
    int index_settings_file;
    int status;
    int previous_index_file;
    
    dictionary = (yamlsettings_dictionary_t *) data;
    unset = dictionary->files[YAMLSETTINGS_FILE_UNSET];

    // Complete real name
    ne = (yamlsettings_node_element *) node->data;

    //if (ne->type > 0) { // This node need to be saved
        complete_name_normalized = yamlsettings_parameter_name_normalized(node);
        option = yamlsettings_dictionary_find_parameter_by_name(complete_name_normalized, dictionary);
        free(complete_name_normalized);
        if (option >= 0) {
            // Set from which file come data to this parameter
            index_settings_file = yamlsettings_parameter_get_position_in_linear_array(&dictionary->parameters[option], node);
            
            if (index_settings_file < yamlsettings_parameter_number_of_elements(&dictionary->parameters[option])) { // If the position is bigger than the max number of elements the value must not be assigned
                // If the value is already assigned we will inform  the user
                if (dictionary->parameters[option].settings_file[index_settings_file] != unset) {
                    complete_name = yamlsettings_parameter_name(node);
                    //fprintf(stderr, "%s:%d: debug: %s : %d\n", __FILE__, __LINE__, complete_name, index_settings_file);
                    if (ne->file_index == YAMLSETTINGS_FILE_COMMAND_LINE && (dictionary->parameters[option].allow_array == YS_PARAM_SCALAR || ne->element_position == 0)) {
                        if (dictionary->parameters[option].settings_file[index_settings_file] == dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE]) {
                            error=g_string_new(NULL);
                            g_string_printf(error, "Parameter %s overwritten by command line repeatedly", complete_name);
                            yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_ERROR);
                            g_string_free(error,true);
                        } else {
                            dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE];
                            error=g_string_new(NULL);
                            g_string_printf(error, "Parameter %s overwritten by command line", complete_name);
                            yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_INFO);
                            g_string_free(error, true);
                        }
                    }
                    if (ne->file_index != YAMLSETTINGS_FILE_COMMAND_LINE && (dictionary->parameters[option].allow_array == YS_PARAM_SCALAR || ne->element_position == 0)) {
                        if (strcmp(dictionary->parameters[option].settings_file[index_settings_file], dictionary->files[ne->file_index]) == 0) {
                            error=g_string_new(NULL);
                            g_string_printf(error, "Parameter %s repeatedly in same file", complete_name);
                            yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_ERROR);
                            g_string_free(error,true);
                        } else {
                            previous_index_file=yamlsettings_dictionary_index_of_file(dictionary,dictionary->parameters[option].settings_file[index_settings_file]);
                            assert(previous_index_file>=0);
                            dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[ne->file_index];
                            // Checking if previuos set was from a overwritable parameter (it comes from template) or not
                            if(dictionary->file_modes[previous_index_file]==YS_PARSE_MODE_APPEND){
                                error=g_string_new(NULL);
                                g_string_printf(error, "Parameter %s repeatedly in different files (%s and %s)", 
                                        complete_name,
                                        dictionary->files[previous_index_file],
                                        dictionary->files[ne->file_index]);
                                yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_WARNING); 
                                g_string_free(error, true);
                            }else if(dictionary->file_modes[previous_index_file]==YS_PARSE_MODE_PREPEND){
                                error=g_string_new(NULL);
                                g_string_printf(error, "Parameter %s was defined in %s and overwritten in %s", 
                                        complete_name,
                                        dictionary->files[previous_index_file],
                                        dictionary->files[ne->file_index]);
                                yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_INFO);      
                                g_string_free(error, true);
                            }else{
                              abort();  
                            }
                        }
                    }
                    free(complete_name);

                } else {
                    dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[ne->file_index];
                }

                if (ne->type == SEQ && dictionary->parameters[option].allow_array == YS_PARAM_SCALAR) { // Add an error if is an array when the option is not defined to be that
                    complete_name = yamlsettings_parameter_name(node);
                    dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[ne->file_index];
                    error=g_string_new(NULL);
                    g_string_printf(error, "Parameter %s cannot be an array", complete_name);
                    yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_WARNING);
                    g_string_free(error,true);
                    free(complete_name);
                }
                if (ne->type == VAL && dictionary->parameters[option].allow_array == YS_PARAM_ARRAY) { // Add an error if is an scalar when the option is not defined to be that
                    complete_name = yamlsettings_parameter_name(node);
                    dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[ne->file_index];
                    error=g_string_new(NULL);
                    g_string_printf(error, "Parameter %s cannot be a scalar. If it is a single value surround it in brackets [ ]", complete_name);
                    yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_WARNING);
                    g_string_free(error,true);
                    free(complete_name);
                }

                // Set parameter data calling assigning function
                status=dictionary->parameters[option].input_function(&(dictionary->status), &(dictionary->parameters[option]), node, dictionary->settings);
              
                // If it was an error I mark the parameter like unset
                if(status<0){
                    dictionary->parameters[option].settings_file[index_settings_file] = dictionary->files[YAMLSETTINGS_FILE_UNSET];
                }

            } else {
                complete_name = yamlsettings_parameter_name(node);
                error=g_string_new(NULL);
                g_string_printf(error, "%s: Array contains more elements than size (maximum: %d, element: %d)", complete_name, yamlsettings_parameter_number_of_elements(&dictionary->parameters[option]), index_settings_file+1);
                yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_ERROR);
                g_string_free(error,true);
                free(complete_name);
            }

        } else {
            complete_name = yamlsettings_parameter_name(node);
            error=g_string_new(NULL);
            g_string_printf(error, "Option %s has not been recognized. Ignored", complete_name);
            yamlsettings_error_add_parse_error(&dictionary->status, (char *)error->str, ne->file_index, YS_WARNING);
            g_string_free(error,true);
            free(complete_name);
        }

    return (FALSE);
}

int yamlsettings_parameter_default_set(yamlsettings_status *status, yamlsettings_parameter *param,GNode *node,void *settings){
    yamlsettings_node_element *ne;    
    char *error, *complete_name,*tmp_string;
    int v=0,i;
    int position;
    int counter;
    int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    int maximums[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    
    assert(param->counter_mem_pos!=NULL || param->dimensions.ndimensions==0);
    
    ne=(yamlsettings_node_element *)node->data;
    
    // Obtain memory position
    yamlsettings_parameter_get_indexes(indexes, maximums, param, node);
    position=indexND(param->dimensions.ndimensions,indexes,maximums);
    if(indexNDisvalid(param->dimensions.ndimensions,indexes,maximums)!=0){ 
            error = (char *) malloc(sizeof (char)*512);
            assert(error!=NULL);
            complete_name = yamlsettings_parameter_name(node);
            tmp_string = (char *) malloc(sizeof (char)*50);
            assert(tmp_string!=NULL);
            sprintf(error, "%s: Invalid indexes are out of bounds (maximums: [", complete_name);
            for (i = 0; i < param->dimensions.ndimensions; i++) {
                if(i!=0) strcat(error,",");
                sprintf(tmp_string,"%d",maximums[i]);
                strcat(error,tmp_string);
            }
            strcat(error,"], indexes: [");
            for (i = 0; i < param->dimensions.ndimensions; i++) {
                if(i!=0) strcat(error,",");
                sprintf(tmp_string,"%d",indexes[i]+1);
                strcat(error,tmp_string);
            }            
            strcat(error,"])");
            yamlsettings_error_add_parse_error(status, error, ne->file_index, YS_ERROR);
            free(error);
            free(tmp_string);
            free(complete_name);
            return -11;          
    }
    if(ne->type > 0){ // if this node is savable 
        v=yamlsettings_data_type_assign(param->data_type_set, 
                status, 
                param->mem_pos, 
                yamlsettings_parameter_number_of_elements(param) ,
                param->name,  
                (gpointer) ne->data,
                position, 
                param->maxlength,
                ne->file_index,
                param->settings_file[yamlsettings_parameter_get_counter_position_in_linear_array(param, node)]); 
        
        ne=(yamlsettings_node_element *)node->parent->data;
        counter=ne->nelements;

    }else{ // if the node is not a leaf of the tree means that it is an empty array
        // we will set it as default value
        for (i = position; i < param->dimensions.dimension_size[param->dimensions.ndimensions-1]; i++) {
                v=yamlsettings_data_type_assign(param->data_type_set, 
                status, 
                param->mem_pos, 
                yamlsettings_parameter_number_of_elements(param) ,
                param->name,  
                param->default_value.arguments[0], // We fill with default value without taking into account the position. It is like if it is erased.
                i, 
                param->maxlength,
                ne->file_index,
                param->settings_file[yamlsettings_parameter_get_counter_position_in_linear_array(param, node)]);
        }

        // and we will set the counter to 0
        counter=0;
    }
    
    // Update counter 
    if(v==0){
        if(param->counter_mem_pos!=NULL){ // If it is the first element and is and array it is necessary to set counter value
            position=indexND(yamlsettings_parameter_counter_number_of_dimensions(param),indexes,maximums);
            if(position<yamlsettings_parameter_counter_number_of_elements(param)){
                // If the parameter is a scalar, arrays should be completed and maximum index never can be removed. So counter value is the maximum.
                if(param->allow_array==YS_PARAM_SCALAR){
                    if(param->counter_mem_pos[position]<indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]+1){          
                        param->counter_mem_pos[position]=indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]+1;
                    }
                }else{ // If the parameter is or can be an array, the maximum can be overwritten if by the user. Example: this array [1, 2, 3] can be overwritten by [1, 2]
                    param->counter_mem_pos[position]=counter;
                }
            }else{
                error = (char *) malloc(sizeof (char)*512);
                assert(error!=NULL);
                complete_name = yamlsettings_parameter_name(node);
                sprintf(error, "%s: Counter can not be updated (maximum: %d, element: %d)", complete_name,yamlsettings_parameter_counter_number_of_elements(param),position+1);
                yamlsettings_error_add_parse_error(status, error, ne->file_index, YS_ERROR);
                free(error);
                free(complete_name);
                return -10;                
            }
        } 
    }    

    return v;
}

char *yamlsettings_parameter_default_debug(yamlsettings_parameter *param) {
    int i,j,tmp=0,count;
    char *tmpchar;
    GString *gresult;
    int trash=1;
    int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    
    if(param->counter_mem_pos==NULL){
       param->counter_mem_pos=&trash; 
    }
    
    gresult = g_string_new(NULL);
    count=0;
    for(i=0;i<yamlsettings_parameter_counter_number_of_elements(param);i++){
        if(param->counter_mem_pos[i]>0){
            if(count!=0){
                g_string_append(gresult, "\n     ");
            }
            
            switch(param->allow_array){
                case YS_PARAM_SCALAR:
                    for(j=0;j<param->counter_mem_pos[i];j++){
                        if(j!=0){
                            g_string_append(gresult, "\n     ");
                        }
                        // Obtaining parameter name
                        indexesofND(i, indexes, yamlsettings_parameter_counter_number_of_dimensions(param),  param->dimensions.dimension_size);
                        indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]=j;//param->counter_mem_pos[i]-1;
                        tmpchar=yamlsettings_parameter_get_with_indexes(param->name, param->dimensions.ndimensions, indexes);
                        g_string_append(gresult, tmpchar);
                        free(tmpchar);
                        g_string_append(gresult, ": ");    
                        tmpchar=yamlsettings_data_type_get_string(param->data_type_get, param->mem_pos, indexND(param->dimensions.ndimensions, indexes ,param->dimensions.dimension_size)/*0*/, param->maxlength);
                        g_string_append(gresult, tmpchar);
                        free(tmpchar);
                    }
                    break;
                case YS_PARAM_OPTIONAL_ARRAY:
                    // Obtaining parameter name
                    indexesofND(i, indexes, yamlsettings_parameter_counter_number_of_dimensions(param),  param->dimensions.dimension_size);
                    indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]=param->counter_mem_pos[i]-1;

                    tmpchar=yamlsettings_parameter_get_with_indexes(param->name, yamlsettings_parameter_counter_number_of_dimensions(param), indexes);
                    g_string_append(gresult, tmpchar);
                    free(tmpchar);
                    g_string_append(gresult, ": ");                    
                    
                    tmp=indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]+1;
                    if(tmp>1){
                        g_string_append(gresult, "[");
                        for(j=0;j<tmp;j++){
                            if(j!=0){
                                g_string_append(gresult, ", ");
                            }
                            indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]=j;
                            tmpchar=yamlsettings_data_type_get_string(param->data_type_get, param->mem_pos, indexND(param->dimensions.ndimensions, indexes ,param->dimensions.dimension_size), param->maxlength);
                            g_string_append_printf(gresult," %s", tmpchar);
                            free(tmpchar);
                        }
                        g_string_append_printf(gresult," ] #%d", tmp);
                    }else{
                        tmpchar=yamlsettings_data_type_get_string(param->data_type_get, param->mem_pos, 0, param->maxlength);
                        g_string_append(gresult, tmpchar);
                        free(tmpchar);
                    }
                    break;
                case YS_PARAM_ARRAY:
                    // Obtaining parameter name
                    indexesofND(i, indexes, yamlsettings_parameter_counter_number_of_dimensions(param),  param->dimensions.dimension_size);
                    indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]=param->counter_mem_pos[i]-1;

                    tmpchar=yamlsettings_parameter_get_with_indexes(param->name, yamlsettings_parameter_counter_number_of_dimensions(param), indexes);
                    g_string_append(gresult, tmpchar);
                    free(tmpchar);
                    g_string_append(gresult, ": ");  
                    
                    tmp=indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]+1;
                        g_string_append(gresult, "[");
                        for(j=0;j<tmp;j++){
                            if(j!=0){
                                 g_string_append(gresult, ", ");
                            }
                            indexes[yamlsettings_parameter_counter_number_of_dimensions(param)]=j;
                            tmpchar=yamlsettings_data_type_get_string(param->data_type_get, param->mem_pos, indexND(param->dimensions.ndimensions, indexes ,param->dimensions.dimension_size), param->maxlength);
                            g_string_append_printf(gresult," %s", tmpchar);
                            free(tmpchar);
                        }
                        g_string_append_printf(gresult," ] #%d", tmp);                
                    break;
                default:
                    fprintf(stderr,"%s: %d: FATAL ERROR: Cannot print parameters", __FILE__, __LINE__);
                    break;
            }
            
            count++;
        }
    }
    if(count==0){
        // Initilize indexes
        for (i = 0; i < YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS; i++) {
            indexes[i]=0;
        }
        tmp=1;
        if(param->allow_array==YS_PARAM_SCALAR){
            tmp=0;
        }

        tmpchar=yamlsettings_parameter_get_with_indexes(param->name, param->dimensions.ndimensions-tmp, indexes);
        g_string_append(gresult, tmpchar);
        free(tmpchar);
        g_string_append(gresult, ": [] #0");  
    }
    
     // Unset trash if it was used.
    if(param->counter_mem_pos==&trash){
        param->counter_mem_pos=NULL;
    }    
    
    tmpchar=(char *) g_string_free(gresult, FALSE);
    return tmpchar;    
}
