#include <assert.h>
#include <string.h>
#include "yamlsettings_process_yaml.h"
#include "yamlsettings_input_yaml.h"
#include "yamlsettings_dictionary.h"
#include <grasp/utils.h>

int yamlsettings_input_yaml_read(yamlsettings_status *status, GNode **cfg, int nparameters, char const *parameters[], char *files[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_settings_file_mode settings_file_mode) {
    int value_returned = 0; // Defined status returned
    GNode *extraparams=NULL;
    GNode *extraparams_root=NULL;
    // Define parser variable
    yaml_parser_t parser;
    // Define a file
    FILE *source;
    int i;
    char *base_path;
    GString *error;
    int nerror=0;
    *cfg = yamlsettings_create_node(VAR, 0, 0, YAMLSETTINGS_FILE_MAIN_CONFIGURATION, "");
    bool ignore_import_files=false, ignore_template_files=false;
    bool nosettingsfile=false;
    yamlsettings_parse_error *parse_error;
    
    assert(nparameters>=1);
    
    // If there are extra parameters... we extract them
    if (nparameters > 1) { 
        // Get a tree with extra parameters
        extraparams_root = yamlsettings_tree_from_extra_parameters(nparameters-1, &parameters[1]);
        base_path = pathoffile(files[YAMLSETTINGS_FILE_COMMAND_LINE]);
        value_returned = yamlsettings_process_imported_files(false, status, extraparams_root, files, file_modes, base_path);              
        if (value_returned < 0) { // I return a value of cut a cicle adding files
            return value_returned;
        }   
        if (value_returned>0){ // Import parameter is replaced by command line
            ignore_import_files=true;
        }
        value_returned = yamlsettings_process_template_files(false, status, extraparams_root, files, file_modes, base_path);              
        if (value_returned < 0) { // I return a value of cut a cicle adding files
            return value_returned;
        }   
        if (value_returned>0){ // Template parameter is replaced by command line
            ignore_template_files=true;
        }  
        free(base_path); 
    }
    
    if(settings_file_mode!=YS_PARSE_SETTINGS_FILE_FORBIDDEN){
        i = YAMLSETTINGS_FILE_MAIN_CONFIGURATION; 
        while (files[i] != NULL) {
            source = fopen(files[i], "rb");
            if (!source) {
                error=g_string_new(NULL);
                
                if(i==YAMLSETTINGS_FILE_MAIN_CONFIGURATION){ // Error if file can not be read is mai file
                    g_string_printf(error,"Can't open file %s", files[i]);
                    nosettingsfile=true;
                }else{
                    if(file_modes[i]==YS_PARSE_MODE_PREPEND){
                        g_string_printf(error,"Can't open template file %s", files[i]);
                    }else if(file_modes[i]==YS_PARSE_MODE_APPEND){
                        g_string_printf(error,"Can't open import file %s", files[i]);
                    }
                }
                
                yamlsettings_error_add_parse_error(status, (char *)error->str, i, YS_FATAL_ERROR);  
                g_string_free(error,true);
                
                nerror--;
            }else{
                yaml_parser_initialize(&parser);
                yaml_parser_set_input_file(&parser, source);

                value_returned = yamlsettings_yaml_file_process(status, &parser, *cfg, i, file_modes[i]); // Recursive parser

                // Closing parser
                yaml_parser_delete(&parser);

                // Closing input file
                fclose(source);

                if(i==YAMLSETTINGS_FILE_MAIN_CONFIGURATION && ignore_import_files==true){ 
                    ignore_import_files=true;
                }else{
                    ignore_import_files=false;
                }
                base_path = pathoffile(files[i]);
                value_returned = yamlsettings_process_imported_files(ignore_import_files, status, *cfg, files, file_modes, base_path);
                if(ignore_import_files==true && value_returned>0){ // If imported files is overwritten by command line
                    yamlsettings_error_add_parse_error(status, "Import statement overwritten by command line", i, YS_INFO);                
                }
                if (value_returned < 0) { // I return a value of cut a cicle adding files
                    free(base_path);
                    return value_returned;
                }
                value_returned = yamlsettings_process_template_files(ignore_template_files, status, *cfg, files, file_modes, base_path);
                if(ignore_template_files==true && value_returned>0){ // If imported files is overwritten by command line
                    yamlsettings_error_add_parse_error(status, "Template statement overwritten by command line", i, YS_INFO);                
                }
                if (value_returned < 0) { // I return a value of cut a cicle adding files
                    free(base_path);
                    return value_returned;
                }            
                free(base_path);
            }
            i++; // Prepare for read next file
        }
    }
    
    // If the settings file is optional we will remove the error or try to open it
    // and can not.
    if(settings_file_mode==YS_PARSE_SETTINGS_FILE_OPTIONAL){
        for (i = 0; i < status->parse_errors->len; i++) {
            if(((yamlsettings_parse_error *)(status->parse_errors->pdata[i]))->type==YS_FATAL_ERROR &&
                ((yamlsettings_parse_error *)(status->parse_errors->pdata[i]))->file_index==YAMLSETTINGS_FILE_MAIN_CONFIGURATION){
                parse_error=g_ptr_array_remove_index (status->parse_errors,i);
                yaml_settings_status_parse_error_destroy(parse_error);
                break;
            }
        }      
    }
    
    // If there was problems opening the main settings file we will add first argument like if it is a settings
    // If someone does not provide any settings file it could be a settings like "help"
    if(nosettingsfile==true){
        yamlsettings_add_parameter_to_a_tree(*cfg,parameters[0]);
    }

    // If there are extra parameters we add them to the tree
    if (nparameters > 1) { 
        // ... and join with cfg
        for (i = 0; i < ((yamlsettings_node_element *) extraparams_root->data)->nelements; i++) {
            extraparams = extraparams_root->children;
            g_node_unlink(extraparams);
            // Add to tree
            g_node_append(*cfg, extraparams);
            // Update information
            ((yamlsettings_node_element *) extraparams->data)->element_position = ((yamlsettings_node_element *) (*cfg)->data)->nelements;
            ((yamlsettings_node_element *) (*cfg)->data)->nelements++;
        }
        // Remove root
        yamlsettings_destroy_tree(extraparams_root);
    }
    
    // Check arrays
    yamlsettings_check_arrays(status, *cfg);
    
    //g_node_traverse(*cfg, G_LEVEL_ORDER, G_TRAVERSE_ALL, -1, yamlsettings_tree_debug, NULL);
    
    return nerror;

}

char *yamlsettings_parameter_name(GNode *node) {
    int i, j;
    char *complete_name;
    GNode *current_node;
    int node_deph = g_node_depth(node);
    int steps = 1;
    yamlsettings_node_element *ne;
    
    ne = (yamlsettings_node_element *) node->data;

    if (ne->type > 0) {
        steps = 2;
    }
    // Initialize vars
    complete_name = (char *) malloc(sizeof (char) * YAMLSETTINGS_MAXSIZEPARAMETERNAME);
    assert(complete_name!=NULL);
    strcpy(complete_name, "");

    if(((yamlsettings_node_element *)(node->data))->data==NULL){
        return complete_name;
    }
    
    for (i = 0; i < node_deph - steps; i++) {
        current_node = node;
        for (j = 0; j < node_deph - 2 - i; j++) {
            current_node = current_node->parent;
        }

        if (g_node_depth(current_node) > 2) strcat(complete_name, ".");
        strcat(complete_name, ((yamlsettings_node_element *) current_node->data)->data);

    }

    return complete_name;
}

char *yamlsettings_parameter_name_normalized(GNode *node) {
    int i, j;
    char *complete_name;
    char normalized_name[YAMLSETTINGS_MAXSIZEPARAMETERNAME];
    char *tmp;
    GNode *current_node;
    int node_deph = g_node_depth(node);
    int steps = 1;
    yamlsettings_node_element *ne;

    ne = (yamlsettings_node_element *) node->data;

    if (ne->type > 0) {
        steps = 2;
    }
    // Initialize vars
    complete_name = (char *) malloc(sizeof (char) * YAMLSETTINGS_MAXSIZEPARAMETERNAME);
    assert(complete_name!=NULL);
    strcpy(complete_name, "");

    for (i = 0; i < node_deph - steps; i++) {
        strcpy(normalized_name, "");
        current_node = node;
        for (j = 0; j < node_deph - 2 - i; j++) {
            current_node = current_node->parent;
        }

        if (g_node_depth(current_node) > 2) strcat(complete_name, ".");
        // Normalizing name
        strcpy(normalized_name, ((yamlsettings_node_element *) current_node->data)->data);
        tmp = strchr(normalized_name, '[');
        if (tmp != NULL) {
            tmp[1] = ']';
            tmp[2] = '\0';
        }
        strcat(complete_name, normalized_name);

    }

    return complete_name;
}

char *yamlsettings_name_of_array_node(const char *name) {
    char *result, *end;

    result = (char *) malloc(sizeof (char)*(strlen(name) + 1));

    strcpy(result, name);

    if(yamlsettings_is_node_array(name)==true){
        end = strchr(result, '[');

        if (end != NULL) {
            end[0] = '\0';
        }
    }
    
    return result;
}

gboolean yamlsettings_is_node_array(const char *name){
    char *tmp;

    tmp=strchr(name,'[');
    if(tmp!=NULL){
        tmp=strchr(tmp,']');
        if(tmp!=NULL){
            return true;
        }
    }
    
    return false;
}

gboolean yamlsettings_is_parameter_array(GNode *node){
    GNode *parent;
    
    parent=node;
    while(!G_NODE_IS_ROOT(parent)){
       parent=parent->parent;
       if(yamlsettings_is_node_array((char *) ((yamlsettings_node_element *) (parent->data))->data)==TRUE){           
           return TRUE; //FOUND!
       }  
    }
        
    return FALSE;

}

int yamlsettings_index_of_array_node(const char *name) {
    char *start, *end, strindex[5];
    int i, index;

    start = strchr(name, '[');

    if (start != NULL) {

        end = strchr(start, ']');
        if(end==NULL){
            return -1; // not is an array
        }
        // Retrieve index number
        start++;
        i = 0;
        while (start != end) {
            strindex[i] = *start;
            i++;
            start++;
        }
        strindex[i] = '\0';

        index = atoi(strindex);
        if (index > 0) {
            return index;
        } else {
            return 0;
        }
    } else {
        return -1;
    }
}

int yamlsettings_find_node(GNode *root, const char *name, GNode **dest) {
    yamlsettings_find_node_traverse_t data;

    data.nnodes = 0;
    data.name = name;
    data.nodes = dest;
    g_node_traverse(root, G_LEVEL_ORDER, G_TRAVERSE_ALL, -1, yamlsettings_find_node_traverse, (gpointer) (&data));

    return data.nnodes;
}

GNode *yamlsettings_first_array_parent(GNode *node){
    GNode *parent;
    
    parent=node;
    while(!G_NODE_IS_ROOT(parent)){
       parent=parent->parent;
       if(yamlsettings_is_node_array((char *) ((yamlsettings_node_element *) (parent->data))->data)==TRUE){
           
           break; //FOUND!
       }  
    }
        
    return parent;
}
