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

void test_datesFormat(){
    char a[]="asd";
    char b[]="01-10-2013";
    char *c;
    
    c=dateDDMMAAAA2AAAAMMDD(a);
    
    CU_ASSERT(strcmp(c,"")==0);
        
    free(c);
    
    c=dateDDMMAAAA2AAAAMMDD(b);
    
    CU_ASSERT(strcmp(c,"2013-10-01")==0);
        
    free(c);  
}

void test_julianDate(){
    time_t time;
    

    CU_ASSERT(convert_string_to_time("2013-01-01T00:00:00Z", &time, 0)==0);
    
    CU_ASSERT_DOUBLE_EQUAL(getJulianDate(time),1.0,0.001);
    
    CU_ASSERT(convert_string_to_time("20130101120000", &time, 1)==0);
    
    CU_ASSERT_DOUBLE_EQUAL(getJulianDate(time),1.5,0.001);    
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
        (NULL == CU_add_test(frameworkSuite   , "test_datesFormat"      , test_datesFormat               )) ||  
        (NULL == CU_add_test(frameworkSuite   , "test_julianDate"       , test_julianDate                ))
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
    
