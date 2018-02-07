/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/**
 * @file grasp_parameters.h
 * @author David Fuertes
 * @date 9 Aug 2014
 * @brief This library help you to work with the parameters array that the inversion code use for retrieving and like initial guess
 *
 * The inversion code works over an array of parameters which describe different atmosphere characteristics.
 * This information is an array for each pixel, so it will be a matrix. You can reproduce all output calling
 * to the forward model with this array. The initial guess of the retrieval code is an array of parameter which 
 * each value represent the initial guess for each caracteristic. This library helps you to work with this matrix
 * and you can use it for working with initial guess or output parameters because the structure is the same.
 */

#ifndef GRASP_PARAMETERS_H
#define	GRASP_PARAMETERS_H

#ifdef	__cplusplus
extern "C" {
#endif

#include "mod_par_OS.inc"
#include "mod_par_inv.inc"    
#include <stdbool.h>

    typedef struct{
        int  n1                          ;
        int  n2        [_KIDIM1]         ;
        int  n3        [_KIDIM1][_KIDIM2];
        int  ISTARSING [_KIDIM1][_KIDIM2];
        int  par_type  [_KIDIM1];
        bool par_retr  [_KIDIM1];
    }par_number_NDIM;    
    

    /**
     * @brief Initialize a matrix structure of parameters
     * @param iguess Matrix of parameters
     * 
     * This function will set -999 for all values of iguess matrix. This function is
     * useful for working with input parameters (initial guess) because if you don't want to define
     * initial guess for each pixel you have to set everything to -999. Then, you can replace (overwrite)
     * the positions that you want to change if you want to specify a specific initial guess for each pixel.
     */
    void grasp_parameters_initialize(float iguess[_KIMAGE][_KPARS]);
    
    /**
     * @brief Retrieve the name (long) of a characteristic knowing its position
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @param characteristic_name Returned name of the characteristic
     * @param size_characteristic_name Size of allocated characteristic_name argument in order to check that no memory leak is produced
     * 
     * Given a parameter number and the dimensions of the parameter array this funtion 
     * will get the characteristic name in characteristic_name string where 
     * size_characteristic_name is the maximum size of the string and should be 
     * bigger than the characteristic name returned
     */
    void grasp_parameters_get_characteristic_type_longname_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, char *characteristic_name, int size_characteristic_name);
    
    /**
     * @brief Retrieve the name (short) of a characteristic knowing its position
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @param characteristic_name Returned name of the characteristic
     * @param size_characteristic_name Size of allocated characteristic_name argument in order to check that no memory leak is produced
     * 
     * Given a parameter number and the dimensions of the parameter array this funtion 
     * will get the characteristic name in characteristic_name string where 
     * size_characteristic_name is the maximum size of the string and should be 
     * bigger than the characteristic name returned
     */
    void grasp_parameters_get_characteristic_type_shortname_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, char *characteristic_name, int size_characteristic_name);    
    
    /**
     * @brief Retrieve the name of the characteristic of a parameter like a unique string take into account the wavelength (or the position of the parameter) and the mode
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @param longname If you want to get a long name or not. If you don't want a long name you'll obtain a short one
     * @param wavelenghts Array of wavelengths from settings (settings->retrieval.WAVE)
     * @param wavelenghts_involved Array of wavelengths involved from settings (settings->retrieval.IWW_SINGL)
     * @param characteristic_name Returned value. This string will be set with a pretty name which describe the parameter
     * @param size_characteristic_name Size of characteristic_name to check the string in in the limits.
     */
    void grasp_parameters_get_characteric_type_pretty_name_by_parameter_number(par_number_NDIM *dimensions, int parameter_number, bool longname, float wavelenghts[_KW], int wavelenghts_involved[_KPARS], char *characteristic_name, int size_characteristic_name);
    
    /**
     * @brief Get the index of the characteristic in NDIM.N1 array
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @return Number of the characteristic or -1 if there was an error
     */
    int grasp_parameters_get_characteristic_index_by_parameter_number(par_number_NDIM *dimensions, int parameter_number);    
    
    /**
     * Get the code of characteristic (example, if characteristic is SizeDistribBin it returns par_type_SD_TB=10101) of a specific parameter knowing its position
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @return The code of characteristic (a valid constant defined in grasp_retrieval_characteristic_type.h) or -1 if there was an error
     */
    int grasp_parameters_get_characteristic_code_by_parameter_number(par_number_NDIM *dimensions, int parameter_number);    
    
    /**
     * Get the mode in which a specific parameter is (starting in 1)
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @return The mode in which a specific parameter is (starting in 1)
     */
    int grasp_parameters_get_mode_by_parameter_number(par_number_NDIM *dimensions, int parameter_number);
    
    /**
     * Get the position inside the block of a specific parameter. This is the position of the parameter in its mode and characteristic tipe
     * @param dimensions Definition of parameter structure
     * @param parameter_number Number of the parameter which the user want to retrieve the name
     * @return The position of the parameter starting to count in its characteristic type and mode
     */
    int grasp_parameters_get_position_by_parameter_number(par_number_NDIM *dimensions, int parameter_number);
    
    /**
     * @brief Look for the index of a characteristic type
     * @param dimensions NDIM block from settings
     * @param characteristic_type Number of the characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @return Index of position of characteristic type in the array of characteristics. If the characteristic is not present it will return -1
     */
    int grasp_parameters_index_of_parameter_type(par_number_NDIM *dimensions, int characteristic_type);    
    
    /**
     * Look for the index of a kind of characteristic. For example, if you want to know where is defined size distribution
     * but it does not matter with kind of size distribution you can call this function with par_type_SD_beg and par_type_SD_end as arguments
     * and you'll get the index of the sd characteristic
     * @param dimensions NDIM block from settings
     * @param begin_characteristic_type Begin of block of characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @param end_characteristic_type End of block of characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @return Index of position of characteristic type in the array of characteristics. If the characteristic is not present it will return -1
     */
    int grasp_parameters_index_of_parameter_type_by_kind_of_parameter(par_number_NDIM *dimensions, int begin_characteristic_type, int end_characteristic_type);   
    
    /**
     * Look for the type of characteristics of the group. For example, if you want to know what kind of size distribution is defined
     * you'll do a call with begin=par_type_SD_beg and end=par_type_SD_end and it will return the code of present size distribution
     * otherwise, if it is not present, -1
     * @param dimensions NDIM block from settings
     * @param begin_characteristic_type Begin of block of characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @param end_characteristic_type End of block of characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @return Index of position of characteristic type in the array of characteristics. If the characteristic is not present it will return -1
     */
    int grasp_parameters_characteristic_code_present_of_kind_of_parameter(par_number_NDIM *dimensions, int begin_characteristic_type, int end_characteristic_type);   
    
    /**
     * @brief Checks if a characteristic was set in settings
     * @param dimensions NDIM block from settings
     * @param characteristic_type Number of the characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @return true if the characteristic is present in current settings definition
     */
    bool grasp_parameters_has_parameter_type(par_number_NDIM *dimensions, int characteristic_type);
    
    /**
     * @brief Set in initial guess array a value
     * @param dimensions Definition of initial guess (NDIM Block)
     * @param iguess Initial guess values array of the pixel
     * @param characteristic_type Characteristic index (base on grasp_retrieval_characteristic_type.h) where the value will be set
     * @param mode Mode which will be set (starting in 1 like in yml settings file)
     * @param pos Index of the value inside the characteristic and mode which will be set (starting in 0 like C array)
     * @param value Value to be set
     * @return index set if the value is well set, -1 if value can not be set because is not strictly positive or -2 if index obtained is outside the valid range [0,_KPARS)
     */
    int grasp_parameters_set_value(par_number_NDIM *dimensions, float iguess[_KPARS], int characteristic_type, int mode , int pos, float value);

    /**
     * This function return number of modes of specific parameter
     * @param dimensions NDIM block from settings
     * @param characteristic_type Number of the characteristic (grasp_retrieval_characteristic_type.h) that you want to look for
     * @return It has to return a number bigger than 0 because when a characteristic type is defined at least has to be one mode. If characteristic type is not defined it will return -1
     */
    int grasp_parameters_number_of_modes_of_parameter(par_number_NDIM *dimensions, int characteristic_type);

    /**
     * Return number of elements of specific mode (mode starts in 1 like in setting file) of a characteristic.
     * @param dimensions NDIM block from settings
     * @param characteristic_type Characteristic index (base on grasp_retrieval_characteristic_type.h) where the value will be set
     * @param mode Mode which will be set (starting in 1 like in yml settings file)
     * @return 
     */
    int grasp_parameters_number_of_elements_of_parameter(par_number_NDIM *dimensions, int characteristic_type, int mode);

    /**
     * This function return a output value of a specific characteristic from an array of output parameters
     * @param dimensions NDIM block from settings
     * @param parameters Array of output parameters
     * @param characteristic_type Characteristic index (base on grasp_retrieval_characteristic_type.h) where the value will be set
     * @param mode Mode which will be set (starting in 1 like in yml settings file)
     * @param pos Index of the value inside the characteristic and mode which will be set (starting in 0 like C array)
     * @return Value from output requested. -1 if Characteristic type does not exist and -2 if characteristic exist but it does not have specific mode
     */
    float grasp_parameters_output_get_value(par_number_NDIM *dimensions, const float parameters[_KPARS], int characteristic_type, int mode , int pos);
    
    /**
     * This function return a the position of a specific characteristic from an array of output parameters
     * @param dimensions NDIM block from settings
     * @param parameters Array of output parameters
     * @param characteristic_type Characteristic index (base on grasp_retrieval_characteristic_type.h) where the value will be set
     * @param mode Mode which will be set (starting in 1 like in yml settings file)
     * @param pos Index of the value inside the characteristic and mode which will be set (starting in 0 like C array)
     * @return Value from output requested. -1 if Characteristic type does not exist and -2 if characteristic exist but it does not have specific mode
     */
    int grasp_parameters_get_position(par_number_NDIM *dimensions, int characteristic_type, int mode , int pos);
    
    /**
     * This function return a value which will be used like initial guess. It is the value of apsing expect if iguess if different from -999.0.
     * If iguess is different of -999 it will be returned except if it is outside of range, in this case min or max will be returned (depending which limit is overflowed)
     * @param dimensions NDIM block from settings
     * @param APSING Initial guess array
     * @param iguess Initial guess array for a specific pixel
     * @param APSMIN Minimum value of initial guess
     * @param APSMAX Maximum value of initial guess
     * @param characteristic_type Characteristic index (base on grasp_retrieval_characteristic_type.h) where the value will be set
     * @param mode Mode which will be set (starting in 1 like in yml settings file)
     * @param pos Index of the value inside the characteristic and mode which will be set (starting in 0 like C array)
     * @return Value from output requested. -1 if Characteristic type does not exist and -2 if characteristic exist but it does not have specific mode
     */
    float grasp_parameters_iguess_get_value(par_number_NDIM *dimensions, float APSING[_KPARS], float iguess[_KPARS], float APSMIN[_KPARS], float APSMAX[_KPARS],  int characteristic_type, int mode , int pos);
    


    

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_PARAMETERS_H */

