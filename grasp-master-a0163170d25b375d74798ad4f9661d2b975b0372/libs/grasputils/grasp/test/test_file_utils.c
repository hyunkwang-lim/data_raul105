/*
 * File:   newcunittest.c
 * Author: david
 *
 * Created on 21-jul-2012, 21:37:59
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/Console.h>
#include <CUnit/Automated.h>
#include <string.h>
#include "../utils/file_utils.h"

/*
 * CUnit Test Suite
 */

int init_suite(void) {
    return 0;
}

int clean_suite(void) {
    return 0;
}

void test_pathoffile(){
    char a[]="a/b/c.d";
    char b[]="/a/b/c.d";
    char c[]="c.d";
    char d[]="";
    char e[]="/c.d";
    char *path;
    
    path=pathoffile(a);
    CU_ASSERT(strcmp("a/b/",path)==0);
    free(path);
    
    path=pathoffile(b);
    CU_ASSERT(strcmp("/a/b/",path)==0);
    free(path);
    
    path=pathoffile(c);
    CU_ASSERT(strcmp("",path)==0);
    free(path);
    
    path=pathoffile(d);
    CU_ASSERT(strcmp("",path)==0);
    free(path);
    
    path=pathoffile(e);
    CU_ASSERT(strcmp("/",path)==0);
    free(path);    
    
}

void test_isabsolute(){
    char a[]="a/b/c.d";
    char b[]="/a/b/c.d";
    char c[]="c.d";
    char d[]="";
    char e[]="/c.d";    

    CU_ASSERT(isabsolute(a)==false);
    CU_ASSERT(isabsolute(b)==true);
    CU_ASSERT(isabsolute(c)==false);
    CU_ASSERT(isabsolute(d)==false);
    CU_ASSERT(isabsolute(e)==true);
}




int main() {
    CU_pSuite frameworkSuite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    
    frameworkSuite = CU_add_suite("file_utils_tests", init_suite, clean_suite);
    if (NULL == frameworkSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }   
    
  
    /* Add the tests to the suite */
    if ((NULL == CU_add_test(frameworkSuite   , "test_pathoffile"       , test_pathoffile                )) ||
        (NULL == CU_add_test(frameworkSuite   , "test_isabsolute"       , test_isabsolute                )) 
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

