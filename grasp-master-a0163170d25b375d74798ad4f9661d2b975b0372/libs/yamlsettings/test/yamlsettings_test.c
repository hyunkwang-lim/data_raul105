/*
 * File:   newcunittest.c
 * Author: david
 *
 * Created on 21-jul-2012, 21:37:59
 */

#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/Console.h>
#include <CUnit/Automated.h>
#include "../yamlsettings_process_yaml.h"
#include "../yamlsettings.h"
#include "../yamlsettings_dictionary.h"
#include "../yamlsettings_validators.h"
#include "../yamlsettings_data_types.h"
#include "../yamlsettings_assign_data.h"
#include <grasp/utils.h>

/*
 * CUnit Test Suite
 */

int init_suite(void) {
    return 0;
}

int clean_suite(void) {
    return 0;
}


void libyamelinput_testNodeFunctions(){
    GNode  *root;
    yamlsettings_node_element *rootdata;
    char badarray[]="array[a]", invalidindexarray1[]="array1[-1]", invalidindexarray2[]="array2[0]", noarray[]="noarray", goodarray[]="v[3]", badformedarray[]="v[4", badformedarray2[]="v][4";
    
    // Make a node
    rootdata=(yamlsettings_node_element *)malloc(sizeof(yamlsettings_node_element));
    assert(rootdata!=NULL);
    rootdata->type=0;
    rootdata->nelements=0;
    rootdata->data=badarray;
    root = g_node_new((gpointer *)rootdata); 
    
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"array")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==0);
    CU_ASSERT(yamlsettings_is_node_array(root)==true);

    rootdata->data=invalidindexarray1;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"array1")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==0);  
    CU_ASSERT(yamlsettings_is_node_array(root)==true);
    
    rootdata->data=invalidindexarray2;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"array2")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==0);      
    CU_ASSERT(yamlsettings_is_node_array(root)==true);
    
    rootdata->data=noarray;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"noarray")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==-1);   
    CU_ASSERT(yamlsettings_is_node_array(root)==false);
    
    rootdata->data=goodarray;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"v")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==3);   
    CU_ASSERT(yamlsettings_is_node_array(root)==true);
    
    rootdata->data=badformedarray;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"v[4")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==-1);      
    CU_ASSERT(yamlsettings_is_node_array(root)==false);
    
    rootdata->data=badformedarray2;
    CU_ASSERT(strcmp(yamlsettings_name_of_array_node(root),"v][4")==0);
    CU_ASSERT(yamlsettings_index_of_array_node(root)==-1);      
    CU_ASSERT(yamlsettings_is_node_array(root)==false);    

}


void libyamlinput_testCreateNode() {
    GNode *x=NULL;
    
    CU_ASSERT(x==NULL);
    
    x=yamlsettings_create_node(VAR, 0, 0, 1, "test");
    
    CU_ASSERT( ((yamlsettings_node_element *)x->data)->type == VAR );
    CU_ASSERT( ((yamlsettings_node_element *)x->data)->nelements == 0 );
    CU_ASSERT( ((yamlsettings_node_element *)x->data)->element_position == 0 );
    CU_ASSERT( strcmp(((yamlsettings_node_element *)x->data)->data,"test")==0 );
}

void libyamlinput_testGetname() {
    GNode  *root, *a, *b, *c;
    /*
     * Preparing appropiate enviroment
     *               root (root)
     *                |
     *   retrieval_input_structure (a)
     *                |
     *              knsing (b)
     *                |
     *                2 (c)
     */
    
    // Make a root
    root=yamlsettings_create_node(VAR, 1, 0, 1,"");
    
    // Add a node
    a=yamlsettings_create_node(VAR, 1, 0,1, "retrieval_input_structure[1]");
    a=g_node_append(root,  a);
    
    // Add a node
    b=yamlsettings_create_node(VAR, 1, 0, 1,"knsing");
    b=g_node_append(a,  b);    
   
    
    // Add a leaf
    c=yamlsettings_create_node(VAL, 1, 0,1, "2");
    c=g_node_append(b,  c);           
    
    
    // Making assertions
    CU_ASSERT(strcmp(yamlsettings_node_name(c),"retrieval_input_structure[1].knsing") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name(b),"retrieval_input_structure[1].knsing") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name(a),"retrieval_input_structure[1]") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name(root),"") == 0);
    
    // Checking normalized version
    CU_ASSERT(strcmp(yamlsettings_node_name_normalized(c),"retrieval_input_structure[].knsing") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name_normalized(b),"retrieval_input_structure[].knsing") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name_normalized(a),"retrieval_input_structure[]") == 0);
    CU_ASSERT(strcmp(yamlsettings_node_name_normalized(root),"") == 0);    
}

void libyamelinput_testParameterTree(){
    GNode *cfg;
    char parameter[]="retrieval.number_of_wavelengths=3";
    char parameter2[]="retrieval.wavelengths=[8 , 4, 5]";
    char parameter3[]="retrieval.constraints.surface_reflectance.BRDF.wavelength[3].initial_guess.wavelength_involved=[8,9,10,11,12,13]";
    
    cfg = yamlsettings_parameter_to_tree (parameter);
    
    
    CU_ASSERT(g_node_n_nodes (cfg,G_TRAVERSE_ALL)==4);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->data))->data , "" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->data))->data , "retrieval" )==0);
    CU_ASSERT(g_node_n_children (cfg->children)==1);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->data))->type == VAR);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->data))->data , "number_of_wavelengths" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->data))->data , "3" )==0);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->children->children->data))->type == VAL);
    
    yamlsettings_destroy_tree(cfg);
    
    cfg = yamlsettings_parameter_to_tree (parameter2);
    
    CU_ASSERT(g_node_n_nodes (cfg,G_TRAVERSE_ALL)==6);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->data))->data , "" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->data))->data , "retrieval" )==0);
    CU_ASSERT(g_node_n_children (cfg->children)==1);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->data))->type == VAR);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->data))->data , "wavelengths" )==0);
    CU_ASSERT(g_node_n_children (cfg->children->children)==3);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->data))->data , "8" )==0);    
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->next->data))->data , "4" )==0);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->children->children->data))->type == SEQ);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->children->children->next->data))->type == SEQ);

    yamlsettings_destroy_tree(cfg);
    
    cfg = yamlsettings_parameter_to_tree (parameter3);
    
    CU_ASSERT(g_node_n_nodes (cfg,G_TRAVERSE_ALL)==14);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->data))->data , "" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->data))->data , "retrieval" )==0);
    CU_ASSERT(g_node_n_children (cfg->children)==1);
    CU_ASSERT( ((yamlsettings_node_element *)(cfg->children->data))->type == VAR);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->children->children->data))->data , "wavelength[3]" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->children->children->children->children->data))->data , "wavelength_involved" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->children->children->children->children->children->data))->data , "8" )==0);
    CU_ASSERT(strcmp( ((yamlsettings_node_element *)(cfg->children->children->children->children->children->children->children->children->next->data))->data , "9" )==0);

    yamlsettings_destroy_tree(cfg);    
}

void libyamlinput_testNormalize() {
    GNode  *root, *a, *b, *c, *d, *e, *f, *g, *h;
    GNode *node;
    char *name;
    /*
     * Preparing appropiate enviroment
     *                root (root)
     *                /         \ 
     *          general (a)     general(f)
     *          /       \           |
     *      knsing (b)  iloop (d)  knsing (g)
     *          |         |         |
     *         2 (c)     3 (e)     5 (h)
     */ 
    
    // Make a root
    root=yamlsettings_create_node(VAR, 2, 0, 1,"");
    
    // Add a node
    a=yamlsettings_create_node(VAR, 2, 0, 1,"general");
    a=g_node_append(root,  a);
    
    
    // Add a node
    b=yamlsettings_create_node(VAR, 1, 0, 1,"knsing");
    b=g_node_append(a,  b);
     
    
    // Add a leaf
    c=yamlsettings_create_node(VAL, 1, 0,1, "2");
    c=g_node_append(b,  c); 

    // Add a node
    d=yamlsettings_create_node(VAR, 1, 1,1, "iloop");
    d=g_node_append(a,  d);

    // Add a node
    e=yamlsettings_create_node(VAL, 0, 0,1, "3");
    e=g_node_append(d,  e);    
    
    // Add a node
    f=yamlsettings_create_node(VAR, 1, 1,1, "general");
    f=g_node_append(root,  f);

    // Add a node
    g=yamlsettings_create_node(VAR, 1, 0,1, "knsing");
    g=g_node_append(f,  g);

    // Add a leaf
    h=yamlsettings_create_node(VAL, 0, 0,1, "5");
    h=g_node_append(g,  h);    
    
    // Making assertions
    CU_ASSERT(g_node_n_nodes(root,G_TRAVERSE_ALL) == 9);
       
    // Remove duplicate and check
    yamlsettings_remove_duplicate_paramters(root);
    
    CU_ASSERT(g_node_n_nodes(root,G_TRAVERSE_ALL) == 7);
     
    //node=root->children->next->children->children;
    node=root->children->next->children->children;
    
    name=yamlsettings_node_name(node);
   
    CU_ASSERT(strcmp(name,"general.knsing")==0);
    
    free(name);
    
    CU_ASSERT( strcmp(((yamlsettings_node_element *)(node->data))->data,"5" )==0);
    
    yamlsettings_normalize_tree(root);
    
    CU_ASSERT(g_node_n_nodes(root,G_TRAVERSE_ALL) == 6);    
    
}
/*
void grasp_settings_testDictionaryInitialize(){
    yamlsettings_dictionary_t *dictionary;
    grasp_settings *settings;
    int index;
    
    // Allocate settings
    settings=(grasp_settings *)malloc(sizeof(grasp_settings));
    
    // Get dictionary for allocated settings
    dictionary=grasp_settings_dictionary_get(settings);

    // Manipulate dictionary
    index=yamlsettings_dictionary_find_parameter_by_name("retrieval.general.threshold_for_stopping", dictionary);
    dictionary->parameters[index].initial_value="3.6";
    index=yamlsettings_dictionary_find_parameter_by_name("retrieval.kernel.phase_matrix_package.path_of_new_kernels", dictionary);
    dictionary->parameters[index].initial_value="test string";
    
    // Initialize
    yamlsettings_dictionary_set_defaults(dictionary);
   
    // Make test
    CU_ASSERT_DOUBLE_EQUAL(dictionary->settings->RIN.EPSP, 3.6, 0.001);
    CU_ASSERT(dictionary->settings->RIN.KL == 0);
    CU_ASSERT(settings->RIN.NOISE.IWLP[0][0][0]==0);
    CU_ASSERT(settings->RIN.NOISE.IWLP[0][0][1]==0);
    CU_ASSERT(settings->RIN.NOISE.NWLP[0][0]==0);
    CU_ASSERT( strcmp( fstr2c(dictionary->settings->RIN.DLSF.distname_N, 255),"test string") ==0);
    CU_ASSERT( strcmp( fstr2c(dictionary->settings->RIN.DLSF.distname_O, 255),"") ==0);
    CU_ASSERT_DOUBLE_EQUAL(dictionary->settings->RIN.SHIFT, 1.0, 0.001);
    
    // Free dictionary memory but keeping grasp_settings structure
    yamlsettings_dictionary_complete_destroy(dictionary);
}

void grasp_settings_testFindRinInputOption() {
    grasp_settings_dictionary_t *d;    
    grasp_settings s;
    d=grasp_settings_dictionary_get(&s);
    CU_ASSERT(grasp_settings_dictionary_find_parameter_by_name("retrieval.general.minimization_convention",d)>=0);
    CU_ASSERT(grasp_settings_dictionary_find_parameter_by_name("retrieval_input_structure.NOEXIST",d)<0);
}
*/

void grasp_settings_validator_test_int(){
    yamlsettings_dictionary_t d;

    yamlsettings_parameter param[]={{"test_param" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   1 ,  "description" ,  { }  ,YAMLSETTINGS_END_VAR}};
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE]={"0","3"};
    int *value;
    int status;
    d.nparameters=1;
    d.parameters=param;
    
    value=(int *)malloc(sizeof(int));
    assert(value!=NULL);
    
    *value=3;
    
    param[0].settings_file="test_file";
    param[0].mem_pos=value;
    
    // Checking one good value
    status=yamlsettings_validator_int(0,&d, arguments);
    CU_ASSERT(status==0);
    
    // Checking one bad value
    *value=4;
    status=yamlsettings_validator_int(0,&d, arguments);
    CU_ASSERT(status==-1);
    
    // Checking multiple values, all of them are good.
    param[0].nelements=3;
    free(value);
    value=(int *)malloc(sizeof(int)*3);
    assert(value!=NULL);
    value[0]=1;
    value[1]=1;
    value[2]=1;
    param[0].mem_pos=value;
    status=yamlsettings_validator_int(0,&d, arguments);
    CU_ASSERT(status==0);    
    
    // Checking multple values, one is bad value
    value[1]=5;
    status=yamlsettings_validator_int(0,&d, arguments);
    CU_ASSERT(status==-1);      
}

void grasp_settings_validator_test_double(){
    yamlsettings_dictionary_t d;
    
    yamlsettings_parameter param[]={{"test_param" , NULL,NULL,NULL , DATA_TYPE_FLOAT  , "0"  ,   1 ,  "description" ,  { }  ,YAMLSETTINGS_END_VAR}};
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE]={"0.0","3.5"};
    double *value;
    int status;
    
    d.nparameters=1;
    d.parameters=param;
     
    value=(double *)malloc(sizeof(double));
    assert(value!=NULL);
    
    *value=3.0;
    
    param[0].settings_file="test_file";
    param[0].mem_pos=value;
    
    // Checking one good value
    status=yamlsettings_validator_double(0,&d, arguments);
    CU_ASSERT(status==0);
    
    // Checking one bad value
    *value=4.0;
    status=yamlsettings_validator_double(0,&d, arguments);
    CU_ASSERT(status==-1);
    
    // Checking in the limits
    *value=3.5;
    status=yamlsettings_validator_double(0,&d, arguments);
    CU_ASSERT(status==0);
    
    // Checking multiple values, all of them are good.
    param[0].nelements=3;
    free(value);
    value=(double *)malloc(sizeof(double)*3);
    assert(value!=NULL);
    value[0]=1.0;
    value[1]=1.0;
    value[2]=1.0;
    param[0].mem_pos=value;
    status=yamlsettings_validator_double(0,&d, arguments);
    CU_ASSERT(status==0);    
    
    // Checking multple values, one is bad value
    value[1]=5.0;
    status=yamlsettings_validator_double(0,&d, arguments);
    CU_ASSERT(status==-1);          
}


void grasp_settings_validator_test_nelements(){
    yamlsettings_dictionary_t d;
    yamlsettings_parameter param[]={{"test_param" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   5 ,  "description" ,  { } ,YAMLSETTINGS_END_VAR }};
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE]={"0","3"};
    int counter;
    int status;
    
    d.nparameters=1;
    d.parameters=param;    
    
    counter=2;
    param[0].settings_file="test_file";
    param[0].counter_mem_pos=&counter;
    
    // Checking one good value
    status=yamlsettings_validator_nelements(0,&d, arguments);
    CU_ASSERT(status==0);
    
    // Checking one bad value
    counter=6;
    status=yamlsettings_validator_nelements(0,&d, arguments);
    CU_ASSERT(status==-1);
    
    // Checking one limit
    counter=0;
    status=yamlsettings_validator_nelements(0,&d, arguments);
    CU_ASSERT(status==0);    
}


void grasp_settings_validator_test_mandatory(){
    yamlsettings_dictionary_t d;
    yamlsettings_parameter param[]={{"test_param" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   1 ,  "description" ,  { } ,YAMLSETTINGS_END_VAR }};
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE];
    int status;
      
    d.nparameters=1;
    d.parameters=param;  
    
    // without assign
    status=yamlsettings_validator_mandatory(0,&d, arguments);
    CU_ASSERT(status==-1);     
    
    // Assigned one time
    param[0].settings_file="test1";
    status=yamlsettings_validator_mandatory(0,&d, arguments);
    CU_ASSERT(status==0);      
    
    // Assigned two times
    param[0].settings_file="test2";
    status=yamlsettings_validator_mandatory(0,&d, arguments);
    CU_ASSERT(status==0);          
}

void grasp_settings_validator_test_directory(){
    yamlsettings_dictionary_t d;
    yamlsettings_parameter param[]={{"test_param" ,  NULL,NULL,NULL , DATA_TYPE_STRING(4)  , "0"  ,   1 ,  "description" ,  { }  ,YAMLSETTINGS_END_VAR} };
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE];
    char path[4]="./";
    int status;
   
    d.nparameters=1;
    d.parameters=param;  
    
    // Transform the path in fortran string
    param[0].mem_pos=(void *)path;
    param[0].settings_file="test_file";
    status=yamlsettings_validator_directory(0,&d, arguments);
    CU_ASSERT(status==0); 
    
    path[2]='x';
    path[3]='\0';
    status=yamlsettings_validator_directory(0,&d, arguments);
    CU_ASSERT(status==-1); 
}

void grasp_settings_validator_test_same_nelements(){
    yamlsettings_dictionary_t d;
    yamlsettings_parameter param[]={
                                        {"test_param1" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   3 ,  "description" ,  { } ,YAMLSETTINGS_END_VAR },
                                        {"test_param2" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   3 ,  "description" ,  { } ,YAMLSETTINGS_END_VAR },
                                  };
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE]={"test_param1"};
    int counter1, counter2;
    int status;
    
    d.nparameters=2;
    d.parameters=param;
    

    counter1=3;
    counter2=3;
    
    param[0].settings_file="test_file";
    param[0].counter_mem_pos=&counter1;
    param[1].settings_file="test_file";
    param[1].counter_mem_pos=&counter2;    
    
    // Checking one good value
    status=yamlsettings_validator_same_nelements(1,&d, arguments);
    CU_ASSERT(status==0);
    
    counter2=2;
    // Checking one bad value
    status=yamlsettings_validator_same_nelements(1,&d, arguments);
    CU_ASSERT(status==-1);
    
    
}


void grasp_settings_validator_test_forbidden(){
    yamlsettings_dictionary_t d;
    yamlsettings_parameter param[]={{"test_param" ,  NULL,NULL,NULL , DATA_TYPE_INT  , "0"  ,   1 ,  "description" ,  { } ,YAMLSETTINGS_END_VAR }};
    char arguments[YAMLSETTINGS_VALIDATOR_MAX_ARGUMENTS][YAMLSETTINGS_MAX_ARGUMENT_SIZE];
    int status;
      
    d.nparameters=1;
    d.parameters=param;  
    
    // without assign
    status=yamlsettings_validator_forbidden(0,&d, arguments);
    CU_ASSERT(status==0);     
    
    // Assigned one time
    param[0].settings_file="test_file1";
    status=yamlsettings_validator_forbidden(0,&d, arguments);
    CU_ASSERT(status==-1);      
    
    // Assigned two times
    param[0].settings_file="test_file2";
    status=yamlsettings_validator_forbidden(0,&d, arguments);
    CU_ASSERT(status==-1);          
}

int main() {
    CU_pSuite libyamlinput = NULL;
    CU_pSuite grasp_settingsSuite = NULL;
    CU_pSuite grasp_settingsValidatorsSuite = NULL;
    

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    libyamlinput = CU_add_suite("libyamlinput_tests", init_suite, clean_suite);
    if (NULL == libyamlinput) {
        CU_cleanup_registry();
        return CU_get_error();
    }
    
    grasp_settingsSuite = CU_add_suite("grasp_settings_test", init_suite, clean_suite);
    if (NULL == grasp_settingsSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }   
    
    grasp_settingsValidatorsSuite = CU_add_suite("grasp_settings_validators_test", init_suite, clean_suite);
    if (NULL == grasp_settingsValidatorsSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }       

    /* Add the tests to the suite */
    if ((NULL == CU_add_test(libyamlinput                 , "libyamlinput_testGetname"                   , libyamlinput_testGetname                    )) ||
        (NULL == CU_add_test(libyamlinput                 , "libyamlinput_testNormalize"                     , libyamlinput_testNormalize                      )) ||
        (NULL == CU_add_test(libyamlinput                 , "libyamelinput_testNodeFunctions"                , libyamelinput_testNodeFunctions                 )) ||
        (NULL == CU_add_test(libyamlinput                 , "libyamelinput_testParameterTree"                , libyamelinput_testParameterTree                 )) ||
        (NULL == CU_add_test(libyamlinput                 , "libyamlinput_testCreateNode"                    , libyamlinput_testCreateNode                     )) ||
            
//        (NULL == CU_add_test(grasp_settingsSuite          , "grasp_settings_testFindRinInputOption"          , grasp_settings_testFindRinInputOption           )) ||
//        (NULL == CU_add_test(grasp_settingsSuite          , "grasp_settings_testDictionaryInitialize"        , grasp_settings_testDictionaryInitialize         )) ||
            
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_int"              , grasp_settings_validator_test_int               )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_double"           , grasp_settings_validator_test_double            )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_nelements"        , grasp_settings_validator_test_nelements         )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_same_nelements"   , grasp_settings_validator_test_same_nelements    )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_directory"        , grasp_settings_validator_test_directory         )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_mandatory"        , grasp_settings_validator_test_mandatory         )) ||
        (NULL == CU_add_test(grasp_settingsValidatorsSuite, "grasp_settings_validator_test_forbidden"        , grasp_settings_validator_test_forbidden         )) 
       ) {
        CU_cleanup_registry();
        return CU_get_error();
    }
  
    /* Run all tests using the CUnit Basic interface */
    //CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
}
