#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>
#include "yamlsettings_dictionary.h"
#include "yamlsettings_validators.h"
#include "yamlsettings_data_types.h"
#include "yamlsettings_input_yaml.h"
#include <grasp/utils.h>

void yamlsettings_initialize_dictionary(yamlsettings_dictionary_t *dictionary, const char *inputFile){
    int i,j;
    char *ptr;
    
    if(YAMLSETTINGS_MAXIMPORTEDFILES<5 || YAMLSETTINGS_MAXIMPORTEDFILES>252){
        fprintf(stderr, "FATAL ERROR: Bad yamlsettings library compilation. YAMLSETTINGS_MAXIMPORTEDFILES constant have to be between 5 and 252\n");
        abort();
    }
    
    // initialize files;
    for(i=0;i<YAMLSETTINGS_MAXIMPORTEDFILES;i++){
        dictionary->files[i]=NULL;
        dictionary->file_modes[i]=YS_PARSE_MODE_APPEND;
    }

    dictionary->files[YAMLSETTINGS_FILE_UNSET] = (char *) malloc(sizeof (char)*12);
    assert(dictionary->files[YAMLSETTINGS_FILE_UNSET]!=NULL);
    strcpy(dictionary->files[YAMLSETTINGS_FILE_UNSET],"__unset__");
    dictionary->files[YAMLSETTINGS_FILE_UNKNOWN] = (char *) malloc(sizeof (char)*13);
    assert(dictionary->files[YAMLSETTINGS_FILE_UNKNOWN]!=NULL);
    strcpy(dictionary->files[YAMLSETTINGS_FILE_UNKNOWN],"__unknown__");
    dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE] = (char *) malloc(sizeof (char)*1000);
    assert(dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE]!=NULL);
    ptr = getcwd(dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE], 1000);
    
    if(ptr==NULL){
        printf("ERROR: libyamlsettings can not retrieve working directory\n");
        exit(-1);
    }
    strcat(dictionary->files[YAMLSETTINGS_FILE_COMMAND_LINE],"/application_name");

    dictionary->files[YAMLSETTINGS_FILE_MAIN_CONFIGURATION] = (char *) malloc(sizeof (char)*(strlen(inputFile)+1));
    assert(dictionary->files[YAMLSETTINGS_FILE_MAIN_CONFIGURATION]!=NULL);
    strcpy(dictionary->files[YAMLSETTINGS_FILE_MAIN_CONFIGURATION],inputFile);
    
    for(i=0;i<dictionary->nparameters;i++){
        dictionary->parameters[i].settings_file = (char **) malloc(sizeof (char *)* yamlsettings_parameter_number_of_elements(&dictionary->parameters[i]));
        assert(dictionary->parameters[i].settings_file!=NULL);
    }
    
    for(i=0;i<dictionary->nparameters;i++){
        for (j = 0; j < yamlsettings_parameter_number_of_elements(&dictionary->parameters[i]); j++) {
            dictionary->parameters[i].settings_file[j]=dictionary->files[YAMLSETTINGS_FILE_UNSET];
        }
    }     
    
    // For status of parse
    dictionary->status.parse_errors=g_ptr_array_new();
    dictionary->status.validation_errors=g_ptr_array_new();
}

void yamlsettings_dictionary_complete_destroy(yamlsettings_dictionary_t *dictionary){
    if(dictionary->settings!=NULL){
        free(dictionary->settings);
    }    
    yamlsettings_dictionary_destroy(dictionary);
}

void yamlsettings_dictionary_destroy(yamlsettings_dictionary_t *dictionary){
    int i;

    if(dictionary->parameters!=NULL){
        for(i=0;i<dictionary->nparameters;i++){
            if(dictionary->parameters[i].name_deallocatable){
                free((char *)dictionary->parameters[i].name);
            }
            free(dictionary->parameters[i].settings_file);
        }
        free(dictionary->parameters);
    }
    
    for(i=0;i<YAMLSETTINGS_MAXIMPORTEDFILES;i++){
        if(dictionary->files[i]!=NULL){
            free(dictionary->files[i]);
        }
    }
        
    yamlsettings_error_status_destroy(&(dictionary->status));

    free(dictionary);
}

void yamlsettings_dictionary_set_defaults(yamlsettings_dictionary_t *dictionary){
    int i,j,max, status;
  
    max=dictionary->nparameters;
    for (i=0; i<max; i++){ 
        //Set counters to 0
        if(dictionary->parameters[i].counter_mem_pos!=NULL){ 
            for(j=0;j<yamlsettings_parameter_counter_number_of_elements(&(dictionary->parameters[i]));j++){ // Inicialize counter if the var is an array
                dictionary->parameters[i].counter_mem_pos[j]=0;
            }
        }
        
        // Set defaults
        status=dictionary->parameters[i].default_value.defaults(i,&(dictionary->parameters[i]),&dictionary->status);
        assert(status==0);
    }

}

int yamlsettings_dictionary_find_parameter_by_name( char *name, yamlsettings_dictionary_t *dictionary){
    int i,max;
    
    max=dictionary->nparameters;
    for (i=0; i<max; i++){
        if(strcmp(name,dictionary->parameters[i].name)==0){
            return i;
        }
    }
    
    return -1;    
}


int yamlsettings_dictionary_validate(yamlsettings_dictionary_t *dictionary){
    int iparam,ivalidator,max;
    int status, errors=0;
    
    max=dictionary->nparameters;
    for (iparam=0; iparam<max; iparam++){
        for(ivalidator=0;ivalidator<YAMLSETTINGS_MAX_VALIDATORS;ivalidator++){
            if(dictionary->parameters[iparam].validators[ivalidator].validator!=NULL){
                status=dictionary->parameters[iparam].validators[ivalidator].validator(iparam,dictionary,dictionary->parameters[iparam].validators[ivalidator].arguments);
                if(status!=0){
                    errors++;
                }
            }
        }
    }
    
    return errors;
}

int yamlsettings_dictionary_index_of_file(yamlsettings_dictionary_t *dictionary, const char *filename){
    int i=0;
        
    while(dictionary->files[i]!=NULL){
        if(strcmp(dictionary->files[i],filename)==0){
            return i;
        }
        i++;
    }
    
    return -1;
}

int yamlsettings_parameter_number_of_elements(yamlsettings_parameter *param){
    int i,result=1;
    
    for(i=0;i<param->dimensions.ndimensions;i++){
        result*=param->dimensions.dimension_size[i];
    }
    
    return result;
}

int yamlsettings_parameter_counter_number_of_elements(yamlsettings_parameter *param){
    int i,result=1;
    
    for(i=0;i<param->dimensions.ndimensions-1;i++){
        result*=param->dimensions.dimension_size[i];
    }
    
    return result;
}

int yamlsettings_parameter_counter_number_of_dimensions(yamlsettings_parameter *param){
    int result;
    
    result=param->dimensions.ndimensions-1;
    if(result<0){
        result=0;
    }
    
    return result;
}

int yamlsettings_parameter_number_of_dimensions_in_string(yamlsettings_parameter *param){
    int result;
    
    result=param->dimensions.ndimensions-1;
    
    if(param->allow_array==YS_PARAM_SCALAR){
        result++;
    }
    
    if(result<0){
        result=0;
    }
    
    return result;    
}

char *yamlsettings_parameter_get_with_indexes(const char *parameter_name, int nindexes, int *indexes){
    char *param;
    char *rest, *token, *ptr; // params for main string tokenizer
    char **tokens;
    int i=0;
    int length=0;
    char *result;
    char tmp[64];

    //assert(nindexes>0);
    // Doing a copy of parameter_name because strtok function will destroy the string over it works
    param = (char *) malloc(sizeof (char)*(strlen(parameter_name)+5));
    assert(param!=NULL);
    strcpy(param,parameter_name);
    strcat(param,";");
    ptr=param;
    // Taking place for point to all tokens
    tokens = (char **) malloc(sizeof (char *)*(nindexes+1));
    assert(tokens!=NULL);
    // Retrieving all parts of the parameter
    while((token = strtok_r(ptr, "[]", &rest))) {
        tokens[i]=token;
        length+=strlen(token)+10; // +10 for take into account brackets and index
        i++;
        if(i>nindexes){
            break;
        }
        ptr = rest;
    }

    //printf("%d == %d en %s\n",i-1,nindexes,parameter_name );
    assert(i-1==nindexes);

    // Preparing result
    result = (char *) malloc(sizeof (char)*(length+3));
    assert(result!=NULL);
    strcpy(result,tokens[0]);

    for(i=0;i<nindexes;i++){
        sprintf(tmp,"[%d]",indexes[i]+1);
        strcat(result,tmp);
        strcat(result,tokens[i+1]);
    }
    
    result[strlen(result)-1]='\0'; // We remove the symbol added at the begin of the routine
    
    // Cleaning memory
    free(param);
    free(tokens);

    return result;
}


void yamlsettings_parameter_get_indexes(int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], int maximums[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], yamlsettings_parameter *param,GNode *node){
    yamlsettings_node_element *ne;  
    int position,tmp,i;
    GNode *parent;
    int index;
    
    ne=(yamlsettings_node_element *)node->data;
    
    indexes[0]=0;
    maximums[0]=1;
    tmp=0;
    if(param->allow_array!=YS_PARAM_SCALAR){ // It could be an array. Last dimension is an array of elements
        tmp=1;
        position=yamlsettings_parameter_counter_number_of_dimensions(param);
        assert(position>=0);
        indexes[position]=ne->element_position; // Set array position
        maximums[position]=param->dimensions.dimension_size[position];
    }
    
    position=param->dimensions.ndimensions-tmp-1; // dimension that will be set
    parent=node;
    
    for(i=0;i<param->dimensions.ndimensions-tmp;i++){ // All indexes that we should obtain from var name
        parent=yamlsettings_first_array_parent(parent);
        index=yamlsettings_index_of_array_node((char *) ((yamlsettings_node_element *) (parent->data))->data); 
        indexes[position]=index-1; // Index in settings start in 1
        maximums[position]=param->dimensions.dimension_size[position];
        position--;
    }       
}

int yamlsettings_parameter_get_position_in_linear_array(yamlsettings_parameter *param,GNode *node){
    int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], maximums[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    int position;
    
    yamlsettings_parameter_get_indexes(indexes, maximums, param, node);
    
    position=indexND(param->dimensions.ndimensions,indexes,maximums);
    
    return position; 
}

int yamlsettings_parameter_get_counter_position_in_linear_array(yamlsettings_parameter *param,GNode *node){
    int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], maximums[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS];
    int position;
    //int tmp=0;
    
    yamlsettings_parameter_get_indexes(indexes, maximums, param, node);

    //if(param->allow_array!=YS_PARAM_SCALAR){ 
    //    tmp=1;
    //}
    
    position=indexND(yamlsettings_parameter_counter_number_of_dimensions(param),indexes,maximums);
    
    return position;    
}

bool yamlsettings_parameter_is_unset(int param_index, yamlsettings_dictionary_t *dictionary){
    int i;
    
    for (i = 0; i < yamlsettings_parameter_number_of_elements(&dictionary->parameters[param_index]); i++) {
        if(dictionary->parameters[param_index].settings_file[i] != dictionary->files[YAMLSETTINGS_FILE_UNSET]){
            return false;
        }
    }
                
    return true;
}


int yamlsettings_parameter_file_index(int param_index, yamlsettings_dictionary_t *dictionary){
    int i;
    int found_file_index=YAMLSETTINGS_FILE_UNSET;
    int parameter_file_index;
    
    
    for (i = 0; i < yamlsettings_parameter_number_of_elements(&dictionary->parameters[param_index]); i++) {
        if(dictionary->parameters[param_index].settings_file[i] != dictionary->files[YAMLSETTINGS_FILE_UNSET]){
            parameter_file_index=yamlsettings_dictionary_index_of_file(dictionary, dictionary->parameters[param_index].settings_file[i]);
            if(found_file_index!=YAMLSETTINGS_FILE_UNSET && found_file_index!=parameter_file_index){
                // If the file is set from different files we will return unknown
                return YAMLSETTINGS_FILE_UNKNOWN;
            }
            found_file_index=parameter_file_index;
        }
    }
                
    return found_file_index;    
}

char *yamlsettings_parameter_name_at_level(const char *parameter, int level){
    char *result=NULL;
    char *parameter_tokenizer;
    char *ptr;
    int l;
    
    // Duplicate parameter
    parameter_tokenizer=(char *)malloc(sizeof(char)*(strlen(parameter)+1));    
    strcpy(parameter_tokenizer,parameter);
    
    l=0;
    ptr = strtok (parameter_tokenizer,".");
    while (ptr != NULL){
        if(l==level){
            result = (char *) malloc(sizeof (char)*(strlen(ptr)+1));
            strcpy(result,ptr);
        }
        ptr = strtok (NULL, ".");
        l++;
    }    
    
    free(parameter_tokenizer);
    
    return result;
}

char *yamlsettings_parameter_name_until_level(const char *parameter, int level){
    char *result=NULL;
    char *parameter_tokenizer;
    char *ptr;
    int l;
        
    // Duplicate parameter
    parameter_tokenizer=(char *)malloc(sizeof(char)*(strlen(parameter)+5));    
    strcpy(parameter_tokenizer,parameter);
    
    result = (char *) malloc(sizeof(char)*(strlen(parameter)+5));
    strcpy(result,"");
    
    l=0;
    ptr = strtok (parameter_tokenizer,".");
    while (ptr != NULL){
        if(l<=level){
            strcat(result,ptr);
            strcat(result,".");
        }
        ptr = strtok (NULL, ".");
        l++;
    }    
    
    free(parameter_tokenizer);
    
    // We remove the last "." character
    result[strlen(result)-1]='\0';
    
    return result;
}

int yamlsettings_parameter_number_of_blocks(const char *param_name){
    char *parameter_tokenizer;
    char *ptr;
    int l;
    
    // Duplicate parameter
    parameter_tokenizer=(char *)malloc(sizeof(char)*(strlen(param_name)+1));    
    strcpy(parameter_tokenizer,param_name);
    
    l=0;
    ptr = strtok (parameter_tokenizer,".");
    while (ptr != NULL){
        ptr = strtok (NULL, ".");
        l++;
    }    
    
    free(parameter_tokenizer);
    
    return l;    
}

int yamlsettings_parameter_name_number_of_dimensions(const char *param_name){
    char *parameter_tokenizer;
    char *ptr;
    int l;
    
    // Duplicate parameter
    parameter_tokenizer=(char *)malloc(sizeof(char)*(strlen(param_name)+4));    
    strcpy(parameter_tokenizer,param_name);
    strcat(parameter_tokenizer,"x"); // We force to not finish in []
    
    l=0;
    ptr = strtok (parameter_tokenizer,"[]");
    while (ptr != NULL){
        ptr = strtok (NULL, "[]");
        l++;
    }    
    
    free(parameter_tokenizer);
    
    return l-1;    
}

int yamlsettings_parameter_number_of_parameters_in_a_block(yamlsettings_dictionary_t *dictionary,int parameter, int level){
    int result=0,i;
    char *l_before=NULL,*l;
    char *tmp;
    
    // ASSERT: it has to be first parameter which has that block
    if(parameter>0){
        l_before=yamlsettings_parameter_name_at_level(dictionary->parameters[parameter-1].name, level);
    }
    l=yamlsettings_parameter_name_at_level(dictionary->parameters[parameter].name, level);
    assert(parameter>0 || (l_before==NULL || (l!=NULL && strcmp(l_before,l)!=0)));
    free(l_before);
    for (i = parameter+1; i < dictionary->nparameters; i++) {
        tmp=yamlsettings_parameter_name_at_level(dictionary->parameters[i].name, level);
        if(tmp!=NULL && strcmp(tmp,l)==0){
            result++;
        }else{
            break;
        }
        free(tmp);
        tmp=NULL;
    }
    
    free(l);
    free(tmp);
    return result+1;
}

bool yamlsetting_parameter_is_not_omissible(yamlsettings_dictionary_t *dictionary,int parameter, int level){
    int npars;
    int ipar,ivalue;
    int (*validator_mandatory)();
    
    validator_mandatory=yamlsettings_validator_mandatory;
    
    npars=yamlsettings_parameter_number_of_parameters_in_a_block(dictionary,parameter,level);
    for (ipar = parameter; ipar < npars+parameter; ipar++) {
        for (ivalue = 0; ivalue < YAMLSETTINGS_MAX_VALIDATORS; ivalue++) {
            if(dictionary->parameters[ipar].validators[ivalue].validator==validator_mandatory ){
                return true;
            }
        }
        for(ivalue=0;ivalue<yamlsettings_parameter_number_of_elements(&dictionary->parameters[ipar]);ivalue++){
            if(dictionary->parameters[ipar].settings_file[ivalue]!=dictionary->files[YAMLSETTINGS_FILE_UNSET]){
                return true;
            }
        }
    }
    return false;
}

/*
 * Private function used by yamlsettings_dictionary_dump_recursive to print the indentation of each parameter
 */
void yamlsettings_dictionary_dump_print_indentation(FILE *f, int level) {
    int i;
    
    for (i = 0; i < level; i++) {
        fprintf(f, "    ");
    }
}

int yamlsettings_dictionary_dump_recursive(FILE *f, yamlsettings_dictionary_t *dictionary,int level, int iparameterstart, int indexes[YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS], int iparameterstop, bool print_defaults){
    int i,j,next_indexstop,iparameter, next_level, nindexes;
    int nelements, nelements_in_block;
    char *title, *title_until, *level_name, *tmp, *tmp1;    
  
    if(iparameterstart>=iparameterstop){
        return 0;
    }

    level_name=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameterstart].name,level); 
    
    // We iterate to print all parameters between
    for (iparameter = iparameterstart; iparameter < iparameterstop; iparameter++) {
        // We analyze the current category of the element
        title=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameter].name,level);    
        title_until=yamlsettings_parameter_name_until_level(dictionary->parameters[iparameter].name,level); // We have to free memory here
        nindexes=yamlsettings_parameter_name_number_of_dimensions(title_until);
        free(title_until);
        assert(nindexes<=YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS);
        
        if(title==NULL){ // This means that we have finish of exploring this level of blocks
            if(iparameter+1<iparameterstop){ 
                // if it is last element of block
                // We have to calculate the level of next parameter to jump
                next_level=0;
                tmp = (char *) malloc(sizeof (char)*1);
                tmp1 = (char *) malloc(sizeof (char)*1);
                do{
                    free(tmp);
                    free(tmp1);
                    tmp=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameter].name,next_level);    
                    tmp1=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameter+1].name,next_level);
                    next_level++; 
               }while(strcmp(tmp,tmp1)==0);
               free(tmp);
               free(tmp1);
               next_level--;
               yamlsettings_dictionary_dump_recursive(f, dictionary,next_level,iparameter+1, indexes,iparameterstop, print_defaults);  
            }
            
            free(title);
            return 0;
        }
        
        if(yamlsettings_is_node_array(title)){ // If this indentation block is an array  
            tmp=yamlsettings_name_of_array_node(title);
            next_indexstop=yamlsettings_parameter_number_of_parameters_in_a_block(dictionary, iparameter, level);
            
            // Initialize indexes
            for (i = nindexes-1; i < YAMLSETTINGS_MAX_PARAMETER_DIMENSIONS; i++) {
               indexes[i]=0;
            }   
            
            // Now we have to know how many elements there are inside this block
            // It will be the maximum in all elements it has
            nelements=0;
            for (j = iparameter; j < iparameter+next_indexstop; j++) {
                indexes[nindexes-1]=0;
                if(dictionary->parameters[j].allow_array==YS_PARAM_SCALAR){
                    nelements_in_block=dictionary->parameters[j].counter_mem_pos[
                        indexND(
                            yamlsettings_parameter_counter_number_of_dimensions(&dictionary->parameters[j]),
                            indexes,
                            dictionary->parameters[j].dimensions.dimension_size)
                        ];
  

                }else{
                    nelements_in_block=0;
                    while(dictionary->parameters[j].counter_mem_pos[  
                            indexND(
                                yamlsettings_parameter_counter_number_of_dimensions(&dictionary->parameters[j]),
                                indexes,
                                dictionary->parameters[j].dimensions.dimension_size)
                            ]!=0){ 
                        nelements_in_block++;
                        indexes[nindexes-1]++;
                    }  
                }   
                if(nelements_in_block>nelements){
                    nelements=nelements_in_block;
                }                      
            }
             
            indexes[nindexes-1]=0;
            if(nelements>0){
                for (j = 0; j < nelements; j++){ 
                    if(print_defaults==true || yamlsetting_parameter_is_not_omissible(dictionary, iparameter, level)==true){
                        // We prepare the indentation
                        yamlsettings_dictionary_dump_print_indentation(f, level);                        
                        fprintf(f, "%s[%d]:\n", tmp, j+1); 
                    }

                    yamlsettings_dictionary_dump_recursive(f, dictionary,level+1,iparameterstart,indexes, iparameter+next_indexstop, print_defaults);
                    indexes[nindexes-1]++;
                }  
            }else{
                if(print_defaults==true){
                    yamlsettings_dictionary_dump_print_indentation(f, level);
                    fprintf(f, "%s[1]:\n", tmp);                     
                    yamlsettings_dictionary_dump_recursive(f, dictionary,level+1,iparameterstart,indexes, iparameter+next_indexstop, print_defaults);
                }
            }
            
            free(tmp);
            iparameterstart=iparameterstart+next_indexstop;
            
            
            // Obtaining next level
            next_level=0;
            tmp = (char *) malloc(sizeof (char)*1);
            tmp1 = (char *) malloc(sizeof (char)*1);
            do{
                free(tmp);
                free(tmp1);
                tmp=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameterstart-1].name,next_level);    
                tmp1=yamlsettings_parameter_name_at_level(dictionary->parameters[iparameterstart].name,next_level);
                next_level++; 
           }while(strcmp(tmp,tmp1)==0);
           level=next_level-1;
        }else{ // if this indentation block is a regular block 
            // We prepare the indentation
            if(print_defaults==true || yamlsetting_parameter_is_not_omissible(dictionary, iparameter, level)==true){  //here we can filter if we print only set values
                yamlsettings_dictionary_dump_print_indentation(f, level);           
                fprintf(f, "%s: ", title);
                
                // If it is the last block we have to print the values
                if(level+1==yamlsettings_parameter_number_of_blocks(dictionary->parameters[iparameter].name)){
                    if(dictionary->parameters[iparameter].allow_array==YS_PARAM_ARRAY || (dictionary->parameters[iparameter].allow_array==YS_PARAM_OPTIONAL_ARRAY && dictionary->parameters[iparameter].counter_mem_pos[
                        indexND(
                            yamlsettings_parameter_counter_number_of_dimensions(&dictionary->parameters[iparameter]),
                            indexes,
                            dictionary->parameters[iparameter].dimensions.dimension_size)
                        ]>1)){
                        fprintf(f, "[");
                        for (i = 0; i < dictionary->parameters[iparameter].counter_mem_pos[
                        indexND(
                            yamlsettings_parameter_counter_number_of_dimensions(&dictionary->parameters[iparameter]),
                            indexes,
                            dictionary->parameters[iparameter].dimensions.dimension_size)
                        ]; i++) {
                            if(i!=0){
                                fprintf(f, ", ");
                            }
                            tmp=yamlsettings_data_type_get_string(dictionary->parameters[iparameter].data_type_get, dictionary->parameters[iparameter].mem_pos, 
                                    indexND(
                                        dictionary->parameters[iparameter].dimensions.ndimensions,
                                        indexes,
                                        dictionary->parameters[iparameter].dimensions.dimension_size)+i                    
                                    , dictionary->parameters[iparameter].maxlength);
                            fprintf(f, "%s", tmp); 
                            free(tmp);                        
                        }
                        fprintf(f, "]");
                    }else{
                        tmp=yamlsettings_data_type_get_string(dictionary->parameters[iparameter].data_type_get, dictionary->parameters[iparameter].mem_pos, 
                                indexND(
                                        dictionary->parameters[iparameter].dimensions.ndimensions,
                                        indexes,
                                        dictionary->parameters[iparameter].dimensions.dimension_size)                   
                                    , dictionary->parameters[iparameter].maxlength);                        
                        fprintf(f, "%s", tmp); 
                        free(tmp);
                    }
                }

                fprintf(f, "\n"); 
            }
            level++;
        }
        yamlsettings_dictionary_dump_recursive(f, dictionary,level,iparameterstart,indexes,iparameterstop, print_defaults);
        if(strcmp(level_name,title)==0){
            return 0;
        }
        free(title);
    }    
    
    return 0;
}




