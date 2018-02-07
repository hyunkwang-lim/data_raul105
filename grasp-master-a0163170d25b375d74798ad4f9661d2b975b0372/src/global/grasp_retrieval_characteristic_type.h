/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_retrieval_characteristic_type.h
 * Author: fuertes
 *
 * Created on 28 de octubre de 2013, 10:34
 */

#ifndef GRASP_RETRIEVAL_CHARACTERISTIC_TYPE_H
#define	GRASP_RETRIEVAL_CHARACTERISTIC_TYPE_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef WARN_DRY
#warning "__CHARACTERISTIC_TYPE__ duplicated"
#endif  
    
#define par_type_aerosol_beg 10000    
    
    // Size Distribution 
#define par_type_SD_beg  10100
#define par_type_SD_TB  10101   //Normalized Size Distribution dV / dlnr at "triangle" bins
#define par_type_SD_LB  10102   //Normalized Size Distribution dV / dlnr for precalculated lognormal bins
#define par_type_SD_LN  10103   //Parameters of bi - modal Lognormal Size Distribution dV / dlnr
#define par_type_SD_end  10200

//Real part of complex Refractive Index or Chemistry
#define par_type_RERI_beg  10200
#define par_type_RERI_spect  10201 //Spectral dependent Real and Imaginary parts of complex refractive index
#define par_type_RERI_const  10202 //Complex Refractive Index is constant
#define par_type_CXRI_nmix  10203  //Real part of complex refractive index is mixture
#define par_type_CXRI_chem  10204  //Chemistry : fwtr, fslbl, finslbl, fsoot, firon
#define par_type_RERI_end  10300

//Imaginary part of complex Refractive Index
#define par_type_IMRI_beg  10300
#define par_type_IMRI_spect  10301 //Spectral dependent Real and Imaginary parts of complex refractive index
#define par_type_IMRI_const  10302 //Complex Refractive Index is constant
//#define par_type_IMRI_nmix  10303  //Real part of complex refractive index is mixture
#define par_type_IMRI_end  10400

//Particles shape(nonsphericity)
#define par_type_SHD_beg  10400
#define par_type_SHD_fsph  10401  //Fraction of spherical particles
#define par_type_SHD_distr  10402 //Axis Ratio Distribution
#define par_type_SHD_end  10500

//Aerosol profile
#define par_type_AVP_beg  10500
#define par_type_AVP_par_height  10501  //Parameter of Aerosol Vertical Profile(...)
#define par_type_AVP_prof  10502 //Normalized Aerosol Vertical Profile
#define par_type_AVP_end  10600

//Aerosol concentration
#define par_type_Cv_beg  10600
#define par_type_Cv  10601    //Aerosol concentration
#define par_type_Cv_end  10700

//Additional parameters
#define par_type_ADD_beg  10700
#define par_type_CL  10701     //Calibration coefficient for lidar
#define	par_type_AVP_par_std 10702  //Standard deviation for vertical profile 
#define par_type_ADD_end 10800
    
#define par_type_aerosol_end 20000

#define par_type_surface_beg 20000
    
//Surface models
#define par_type_SURF1_land_beg  20100
#define par_type_SURF1_land_Ross_Li_BRDF  20101
#define par_type_SURF1_land_RPV_BRDF  20102
#define par_type_SURF1_land_Litvinov  20103
#define par_type_SURF1_land_Litvinov_fast 20104    
#define par_type_SURF1_land_end  20200

#define par_type_SURF2_land_beg  20200
#define par_type_SURF2_land_Maignan_Breon  20201
#define par_type_SURF2_land_Litvinov  20202
#define par_type_SURF2_land_end  20300
    
#define par_type_SURF_water_beg  20300    
#define par_type_SURF_water_Cox_Munk_iso  20301
#define par_type_SURF_water_Cox_Munk_ani  20302
#define par_type_SURF_water_Litvinov  20303
#define par_type_SURF_water_end  20400

#define par_type_surface_end 30000

#define par_type_cloud_beg 30000

//Size Distribution
#define par_type_SD_cloud_beg  30100
#define par_type_SD_cloud_par  30101   // Parameters of Cloud Size Distributi
#define par_type_SD_cloud_end  30200

// Cloud profile
#define par_type_CVP_beg  30400
#define par_type_CVP_par  30401  // Parameters of Cloud Vertical Profile (
#define par_type_CVP_end  30500

// Cloud fraction
#define par_type_fcloud_beg 30500
#define par_type_fcloud_par 30501  // Pixel Cloud Fraction
#define par_type_fcloud_end 30600

#define par_type_cloud_end  40000

#ifdef	__cplusplus
        }
#endif

#endif	/* GRASP_RETRIEVAL_CHARACTERISTIC_TYPE_H */

