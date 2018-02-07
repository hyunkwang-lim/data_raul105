#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <assert.h>
#include <grasp/utils.h> // for isabsolute
#include "yamlsettings_process_yaml.h"
#include "yamlsettings.h"

int yamlsettings_yaml_file_add(yamlsettings_status *status, const char *file, char *files[YAMLSETTINGS_MAXIMPORTEDFILES], yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path, yamlsettings_parser_mode mode){
    int i=0;
    char *copystring;
    int copylength;
    char error[256];
    
    // Finding last position in array
    while(files[i]!=NULL && i<YAMLSETTINGS_MAXIMPORTEDFILES+1){
        i++;
    }
    
    if(i==YAMLSETTINGS_MAXIMPORTEDFILES){ // Error, array is full
        sprintf(error, "Reached maximum of imported files. Try to reduce the number of input files or perhaps exist some loop");
        yamlsettings_error_add_parse_error(status, error, YAMLSETTINGS_FILE_MAIN_CONFIGURATION, YS_ERROR); 
        //return -1;
    }else{
        // Calculating size of result string
        copylength=strlen(file)+1;
        if(isabsolute(file)==false){
            copylength=copylength+strlen(base_path);
        }
        // Allocating string
        copystring=(char *)malloc(sizeof(char)*copylength);
        assert(copystring!=NULL);    
        
        // Add a copy of string
        if(isabsolute(file)==true){
            strcpy(copystring, "");
        }else{
            strcpy(copystring, base_path );
        }
        strcat(copystring, file);
        files[i]=copystring;
        file_modes[i]=mode;
    }
    
    return 0;
}

int yamlsettings_process_imported_and_template_files(bool ignore, yamlsettings_status *status, GNode *cfg, char *files[YAMLSETTINGS_MAXIMPORTEDFILES],yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path, yamlsettings_parser_mode mode){
    GNode *current; 
    GNode *node_file;
    yamlsettings_node_element *ne;
    char *tmp;
    int found=0;
    
    current=cfg->children;
    while(current!=NULL){
        tmp=yamlsettings_parameter_name(current);
         
        if(mode==YS_PARSE_MODE_APPEND && strcmp(tmp,"import")==0){ // If I found import node
            found=1;
            g_node_unlink(current); 
             
            node_file = current->children;
            while(node_file!=NULL){
                ne=(yamlsettings_node_element *)node_file->data;
                
                if(ne->type==SEQ || ne->type==VAL){
                    if(ignore==false){
                        yamlsettings_yaml_file_add(status, (char *)ne->data, files,file_modes, base_path, YS_PARSE_MODE_APPEND);  
                    }
                }
                node_file=node_file->next;
            }
            ((yamlsettings_node_element *) cfg->data)->nelements--;
            yamlsettings_destroy_tree(current); // Freeing memory
            current=NULL;
        }else if(mode==YS_PARSE_MODE_PREPEND && strcmp(tmp,"template")==0){ // If I found import node
            found=1;
            g_node_unlink(current); 
             
            node_file = current->children;
            while(node_file!=NULL){
                ne=(yamlsettings_node_element *)node_file->data;
                
                if(ne->type==SEQ || ne->type==VAL){
                    if(ignore==false){
                        yamlsettings_yaml_file_add(status, (char *)ne->data, files, file_modes, base_path, YS_PARSE_MODE_PREPEND);  
                    }
                }
                node_file=node_file->next;
            }
            ((yamlsettings_node_element *) cfg->data)->nelements--;
            yamlsettings_destroy_tree(current); // Freeing memory
            current=NULL;            
        }else{
            current=current->next;
        }
        free(tmp);
    }
    
    return found;
}

int yamlsettings_process_imported_files(bool ignore, yamlsettings_status *status, GNode *cfg, char *files[YAMLSETTINGS_MAXIMPORTEDFILES],yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path){
    return yamlsettings_process_imported_and_template_files(ignore,status,cfg,files,file_modes,base_path, YS_PARSE_MODE_APPEND);
}

int yamlsettings_process_template_files(bool ignore, yamlsettings_status *status, GNode *cfg, char *files[YAMLSETTINGS_MAXIMPORTEDFILES],yamlsettings_parser_mode file_modes[YAMLSETTINGS_MAXIMPORTEDFILES], const char *base_path){
    return yamlsettings_process_imported_and_template_files(ignore,status,cfg,files,file_modes,base_path, YS_PARSE_MODE_PREPEND);
}

int yamlsettings_yaml_file_process(yamlsettings_status *status, yaml_parser_t *parser, GNode *data, char file_index, yamlsettings_parser_mode mode) {
    GNode *last_leaf = data;
    yaml_event_t event;
    int storage = VAR; // mapping cannot start with VAL definition without VAR key
    yamlsettings_node_element *real_data;
    int returned_value=0;
    char *token_name;
    GString *error;
    
    while (1) {
    	yaml_parser_parse(parser, &event);

    	if (event.type == YAML_NO_EVENT) {
            token_name=yamlsettings_parameter_name(last_leaf);
            error=g_string_new(NULL);
            g_string_printf(error,"Yaml file can't be parsed. Last token read: %s", token_name);
            yamlsettings_error_add_parse_error(status, (char *)error->str, file_index, YS_ERROR);
            g_string_free(error,true);
            free(token_name);
            return -3;
    	}                
        
    	// Parse value either as a new leaf in the mapping
    	//  or as a leaf value (one of them, in case it's a sequence)
    	if (event.type == YAML_SCALAR_EVENT) {
                real_data=(yamlsettings_node_element *)malloc(sizeof(yamlsettings_node_element));
                assert(real_data!=NULL);
                real_data->type=storage;
                real_data->nelements=0;
                real_data->element_position=0;
                real_data->file_index=file_index;
                real_data->data=g_strdup((gchar*) event.data.scalar.value);
                
    		if (storage){ // If its a value
                    //if(storage==SEQ){ // If it is uncommented we can do differences between arrays and simple values.
                        real_data->element_position=((yamlsettings_node_element *)last_leaf->data)->nelements;
                        ((yamlsettings_node_element *)last_leaf->data)->nelements++;
                    //}    

                    g_node_append_data(last_leaf, (gpointer *)real_data);
    		}else{
                    if(mode==YS_PARSE_MODE_APPEND){
                        last_leaf = g_node_append(data,  g_node_new((gpointer *)real_data));
                    }else if(mode==YS_PARSE_MODE_PREPEND){
                        last_leaf = g_node_prepend(data, g_node_new((gpointer *)real_data));
                    }
    		}
                storage ^= VAL; // Flip VAR/VAL switch for the next event
                if(storage>2){
                    storage=2;
                }
    	
    	// Sequence - all the following scalars will be appended to the last_leaf
        }else if (event.type == YAML_SEQUENCE_START_EVENT){
            storage = SEQ;
    	}else if (event.type == YAML_SEQUENCE_END_EVENT){
            storage = VAR;
        }else if (event.type == YAML_MAPPING_START_EVENT) { // depth += 1
    		returned_value=yamlsettings_yaml_file_process(status, parser, last_leaf, file_index, mode);
                if (returned_value<0){
                    return returned_value;
                }
    		storage ^= VAL; // Flip VAR/VAL, w/o touching SEQ
    	}else if (event.type == YAML_MAPPING_END_EVENT || event.type == YAML_STREAM_END_EVENT){ // depth -= 1
            break; // END OF LOOP    
        } 

    	yaml_event_delete(&event);
    }
    
    return 0;
}         
 
typedef struct yamlsettings_check_arrays_struct_{
    yamlsettings_status *status;
    GNode *root;
}yamlsettings_check_arrays_struct;

void yamlsettings_check_arrays(yamlsettings_status *status, GNode *root){
    yamlsettings_check_arrays_struct data;
    data.status=status;
    data.root=root;
    g_node_traverse(root, G_LEVEL_ORDER, G_TRAVERSE_NON_LEAVES, -1, yamlsettings_check_array_elements, (gpointer)(&data));
}

gboolean yamlsettings_check_array_elements(GNode *node, gpointer data){
    char *nameorigin, *parameter_origin, paramnamedest[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENT_SIZE], *tmp;
    int indexorigin;
    GNode *root;
    gpointer searching_data[2];
    yamlsettings_check_arrays_struct *indata;
    yamlsettings_status *status;
    char error [512];
    
    indata=data;
    root=(GNode *)(indata->root);
    status=indata->status;
    
    indexorigin=yamlsettings_index_of_array_node((char *) ((yamlsettings_node_element *) (node->data))->data);    
        
    if(indexorigin==0){
        sprintf(error,"Index of %s must be greater than 0", yamlsettings_parameter_name(node));
        yamlsettings_error_add_parse_error(status, error, ((yamlsettings_node_element *)node->data)->file_index, YS_ERROR);
    }

    if(indexorigin>1){ // It is an array
        // Prepare name for searching
        parameter_origin=yamlsettings_parameter_name(node);
        nameorigin=yamlsettings_name_of_array_node( (char *) ((yamlsettings_node_element *) (node->data))->data);
        tmp=strrchr(parameter_origin,'.');
        *tmp='\0';
        sprintf(paramnamedest,"%s.%s[%d]",parameter_origin,nameorigin, indexorigin-1);
        free(parameter_origin);
        free(nameorigin);

        // Prepare data for traverse function
        searching_data[0] = (gpointer)paramnamedest;
        searching_data[1] = NULL;
        
        // Look for
        g_node_traverse(root, G_LEVEL_ORDER, G_TRAVERSE_NON_LEAVES, -1, yamlsettings_check_array_elements_find_name, (gpointer)searching_data);

        // Check results
        if(searching_data[1]==NULL){
            sprintf(error,"Invalid settings: %s is an incomplete array because element %d doesn't exist",yamlsettings_parameter_name(node),indexorigin-1);
            yamlsettings_error_add_parse_error(status, error, ((yamlsettings_node_element *)node->data)->file_index, YS_ERROR);
        }
    }
    
    return FALSE;
}

gboolean yamlsettings_check_array_elements_find_name(GNode *node, gpointer data){
    char *searching_name;
    gpointer *input_data;
    char *current_parameter_name;
    input_data=data;
    
    // Prepare data for comparison
    searching_name=(char *)(*input_data);
    current_parameter_name=yamlsettings_parameter_name(node);
    
    // Compare
    if(strcmp(searching_name,current_parameter_name)==0){   
        // If it is found, stop and return value
        *(input_data+1)=node;     
        free(current_parameter_name);
        return TRUE;
    }
    
    // It is not found
    free(current_parameter_name);
    return FALSE;
}

gboolean yamlsettings_find_node_traverse(GNode *node, gpointer data){
    char *tmp;
    yamlsettings_find_node_traverse_t *gsfntt;
    
    gsfntt=(yamlsettings_find_node_traverse_t *)data;
    
    tmp=yamlsettings_parameter_name_normalized(node);
    
    if(!G_NODE_IS_LEAF(node) && strcmp(tmp,gsfntt->name)==0){
        gsfntt->nodes[gsfntt->nnodes]=node;
        gsfntt->nnodes++;
        assert(gsfntt->nnodes<YAMLSETTINGS_PROCESS_FIND_MAX_NODES);
    }
    
    free(tmp);
    
    return FALSE;
}

void yamlsettings_add_parameter_to_a_tree(GNode *root, const char string_parameter[]){
    GNode *parameter;
    GNode *parameter_root;
    
    // Get tree
    parameter_root=yamlsettings_parameter_to_tree(string_parameter);
    // Remove root
    parameter=parameter_root->children;
    g_node_unlink(parameter);
    yamlsettings_destroy_tree(parameter_root); 
    // Add to tree
    g_node_append(root,parameter);
    // Update information
    ((yamlsettings_node_element *)parameter->data)->element_position=((yamlsettings_node_element *)root->data)->nelements; 
    ((yamlsettings_node_element *)root->data)->nelements++;  
}

GNode *yamlsettings_tree_from_extra_parameters(int nparameters, char const *parameters[]){
    int i;
    GNode *root;

    root=yamlsettings_create_node(VAR, 0, 0, YAMLSETTINGS_FILE_COMMAND_LINE, "");

    for(i=0;i<nparameters;i++){
        yamlsettings_add_parameter_to_a_tree(root, parameters[i]);
    }
    
    return root;
}

GNode *yamlsettings_parameter_to_tree(char const *parameter){
    GNode *root;
    GNode *current;
    GNode *last_node;
    char *equal;
    char *parameter_tokenizer;
    char *ptr;
    
    // Duplicate parameter
    parameter_tokenizer=(char *)malloc(sizeof(char)*(strlen(parameter)+8));    
    assert(parameter_tokenizer!=NULL);
    strcpy(parameter_tokenizer,parameter);
    
    root=yamlsettings_create_node(VAR, 1, 0, YAMLSETTINGS_FILE_COMMAND_LINE, "");   

    // Simulate that the string end in the = symbol
    equal=strchr(parameter_tokenizer,'=');
    if(equal==NULL){
        // Adding end to string because can be that string don't have it
        strcat(parameter_tokenizer, "="); //empty string ("")
        equal=strchr(parameter_tokenizer,'=');
    }
    *equal='\0';
    
    last_node=root;
    ptr = strtok (parameter_tokenizer,".");
    while (ptr != NULL){
        current=yamlsettings_create_node(VAR, 1, 0, YAMLSETTINGS_FILE_COMMAND_LINE, ptr);
        current=g_node_append(last_node, current); 
        
        last_node=current;
        ptr = strtok (NULL, ".");
    }    
    
    // Add last node with the information after equal symbol
    equal++;
    current=yamlsettings_create_node(VAR, 1, 0, YAMLSETTINGS_FILE_COMMAND_LINE, equal);
    current=g_node_append(last_node, current); 
    last_node=current;    
   
    if(((yamlsettings_node_element *)(last_node->data))->data[0]=='['){ // It is an array
        ((yamlsettings_node_element *)(last_node->parent->data))->nelements=0;
        ptr = strtok (((yamlsettings_node_element *)(last_node->data))->data,"[, ]");
        while (ptr != NULL){ 
            current=yamlsettings_create_node(SEQ, 0, ((yamlsettings_node_element *)(last_node->parent->data))->nelements, YAMLSETTINGS_FILE_COMMAND_LINE, ptr);
            current=g_node_append(last_node->parent, current);
            ((yamlsettings_node_element *)(last_node->parent->data))->nelements++;

            ptr = strtok (NULL, "[, ]");
        }
        
        yamlsettings_destroy_tree(last_node);
    }else{
        ((yamlsettings_node_element *)(last_node->data))->type=VAL;
    }
    free(parameter_tokenizer);
    return root;
}

GNode *yamlsettings_create_node(int type, int nelements, int element_position, char file_index, char *data){
    GNode  *result;
    yamlsettings_node_element *resultdata;
            
    // Make a root
    resultdata=(yamlsettings_node_element *)malloc(sizeof(yamlsettings_node_element));
    assert(resultdata!=NULL);
    resultdata->type=type;
    resultdata->nelements=nelements;
    resultdata->element_position=element_position;
    resultdata->file_index=file_index;
    resultdata->data=(char *)malloc(sizeof(char)*(strlen(data)+1));
    assert(resultdata->data!=NULL);
    strcpy(resultdata->data,data);
    
    result = g_node_new((gpointer *)resultdata);
    
    return result;
}

void yamlsettings_destroy_tree(GNode *node){
    g_node_traverse(node, G_PRE_ORDER, G_TRAVERSE_ALL, -1, yamlsettings_destroy_node_data, NULL);
    g_node_destroy(node);    
}

gboolean yamlsettings_destroy_node_data(GNode *node, gpointer in){
    yamlsettings_node_element *data;
    
    data=(yamlsettings_node_element *)node->data;
    free(data->data);
    free(node->data);
    return (FALSE);
}

gboolean yamlsettings_tree_debug(GNode *node, gpointer data) {
    yamlsettings_node_element *ne;
    int i = g_node_depth(node);
    
    ne=(yamlsettings_node_element *)node->data;
    
    
    while (--i) printf("  ");
    printf("[depth=%d][type=%d][nelements=%d][element_position=%d][file=%d] %s\n", g_node_depth(node), ne->type, ne->nelements, ne->element_position,ne->file_index, ne->data);
    return(FALSE);
}
