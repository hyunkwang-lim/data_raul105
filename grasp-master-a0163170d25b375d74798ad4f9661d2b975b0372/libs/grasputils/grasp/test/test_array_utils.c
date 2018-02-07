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

void test_indexArray(){
    int ndimension;
    int indexes[4];
    int dimensions[4];
    
    //v[3][5]
    CU_ASSERT(index2D(0,0,0)==0);    
    CU_ASSERT(index2D(1,3,5)==8);  
    CU_ASSERT(index2D(2,4,5)==14);  
    
    //v[3][5][2]
    CU_ASSERT(index3D(0,0,0,5,2)==0);
    CU_ASSERT(index3D(1,0,0,5,2)==10);
    CU_ASSERT(index3D(1,1,1,5,2)==13);
    CU_ASSERT(index3D(2,4,1,5,2)==29);
    
    //same test using arrayND function
    // first definition
    ndimension=2;
    dimensions[0]=3;
    dimensions[1]=5;
    
    indexes[0]=0;
    indexes[1]=0;
    CU_ASSERT(index2D(0,0,0)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index2D(0,0,0),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==0);
    CU_ASSERT(indexes[1]==0);
    
    indexes[0]=1;
    indexes[1]=3;
    CU_ASSERT(index2D(1,3,5)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index2D(1,3,5),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==1);
    CU_ASSERT(indexes[1]==3);
    
    indexes[0]=2;
    indexes[1]=4;
    CU_ASSERT(index2D(2,4,5)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index2D(2,4,5),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==2);
    CU_ASSERT(indexes[1]==4);
    
    // second definition
    ndimension=3;
    dimensions[0]=3;
    dimensions[1]=5;
    dimensions[2]=2;
    
    indexes[0]=0;
    indexes[1]=0;
    indexes[2]=0;
    CU_ASSERT(index3D(0,0,0,5,2)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index3D(0,0,0,5,2),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==0);
    CU_ASSERT(indexes[1]==0);
    CU_ASSERT(indexes[2]==0);
    
    
    indexes[0]=1;
    indexes[1]=0;
    indexes[2]=0;
    CU_ASSERT(index3D(1,0,0,5,2)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index3D(1,0,0,5,2),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==1);
    CU_ASSERT(indexes[1]==0);
    CU_ASSERT(indexes[2]==0);
    
    indexes[0]=1;
    indexes[1]=1;
    indexes[2]=1;
    CU_ASSERT(index3D(1,1,1,5,2)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index3D(1,1,1,5,2),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==1);
    CU_ASSERT(indexes[1]==1);
    CU_ASSERT(indexes[2]==1);
    
    indexes[0]=2;
    indexes[1]=4;
    indexes[2]=1;
    CU_ASSERT(index3D(2,4,1,5,2)==indexND(ndimension,indexes,dimensions));
    //testing reverse function
    indexesofND(index3D(2,4,1,5,2),indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==2);
    CU_ASSERT(indexes[1]==4);
    CU_ASSERT(indexes[2]==1);
    
    // third definition
    ndimension=4;
    dimensions[0]=1;
    dimensions[1]=2;
    dimensions[2]=2;
    dimensions[3]=48;
    
    indexes[0]=0;
    indexes[1]=0;
    indexes[2]=0;
    indexes[3]=0;
    CU_ASSERT(indexND(ndimension,indexes,dimensions)==0);
    //testing reverse function
    indexesofND(0,indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==0);
    CU_ASSERT(indexes[1]==0);
    CU_ASSERT(indexes[2]==0);
    CU_ASSERT(indexes[3]==0);
    
    indexes[0]=0;
    indexes[1]=0;
    indexes[2]=0;
    indexes[3]=47;
    CU_ASSERT(indexND(ndimension,indexes,dimensions)==47);
    //testing reverse function
    indexesofND(47,indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==0);
    CU_ASSERT(indexes[1]==0);
    CU_ASSERT(indexes[2]==0);
    CU_ASSERT(indexes[3]==47);
    
    indexes[0]=0;
    indexes[1]=1;
    indexes[2]=0;
    indexes[3]=0;
    CU_ASSERT(indexND(ndimension,indexes,dimensions)==96);
    //testing reverse function
    indexesofND(96,indexes,ndimension,dimensions);
    CU_ASSERT(indexes[0]==0);
    CU_ASSERT(indexes[1]==1);
    CU_ASSERT(indexes[2]==0);
    CU_ASSERT(indexes[3]==0);
    
    
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
        (NULL == CU_add_test(frameworkSuite   , "test_indexArray"       , test_indexArray                )) 
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
    
