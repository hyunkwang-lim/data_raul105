#include <stdio.h>
#include <stdlib.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <CUnit/Console.h>
#include <CUnit/Automated.h>
#include <time.h>
#include "../utils/time_utils.h"

/*
 * CUnit Test Suite
 */

int init_suite(void) {
    return 0;
}

int clean_suite(void) {
    return 0;
}

void test_strfortranlen(){
    char s[10]="hello";
    int length;
    
    // Check length of "usual" string
    s[5]=' ';
    s[6]=' ';
    s[7]=' ';
    s[8]=' ';
    s[9]=' ';    
    length=strfortranlen(s,10);   
    CU_ASSERT(length==5);
    
    // Check length of empty string
    s[0]=' ';
    s[1]=' ';
    s[2]=' ';
    s[3]=' ';
    s[4]=' ';
    length=strfortranlen(s,10);    
    CU_ASSERT(length==0);    
    
    // Check length of full string
    s[9]='x';
    length=strfortranlen(s,10);    
    CU_ASSERT(length==10);    
    
}



int main() {
    CU_pSuite frameworkSuite = NULL;

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* Add a suite to the registry */
    
    frameworkSuite = CU_add_suite("time_utils_tests", init_suite, clean_suite);
    if (NULL == frameworkSuite) {
        CU_cleanup_registry();
        return CU_get_error();
    }   
    
  
    /* Add the tests to the suite */
    if (
        (NULL == CU_add_test(frameworkSuite   , "test_strfortranlen"       , test_strfortranlen                )) 
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
    
