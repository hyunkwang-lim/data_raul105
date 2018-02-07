/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_segment_result.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 14:56
 */

#ifndef GRASP_OUTPUT_SEGMENT_RESULT_H
#define	GRASP_OUTPUT_SEGMENT_RESULT_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdbool.h>
#include "mod_par_inv.inc"
#include "mod_par_OS.inc"
#include "mod_par_DLS.inc"
#include "mod_par_DLS_bin.inc"
#include "../settings/grasp_products.h"
#include "../input/grasp_input_segment.h"
    
#ifdef WARN_DRY
#warning "__RETRIEVAL_OUTPUT_DEFINITION__ binded"
#endif          
    
#define index_ext 1
#define index_ssa 2
#define index_lr  3

// Retrieval output
// Structure contains retrieval output data :
// res       - single pixel total residual (meas + smothness constrain )
// resa,resr - detailed absolut and relative measurement residuals for single pixel 
// niter       - number of iterations
// rest        - total residual for multi-pixel retrieval
// resat,resrt - detailed absolut and relative measurement residuals for segment
// par         - retrieved aerosol and surface reflectance parameters
// ngrid       - number of grid radii for SD (saved to print output)
// radius      - grid radii (saved to print output)
// SDL         - if applicable, precalculated lognormal bins SD (saved to print output)
   
    typedef struct output_pixel_residual_{
         int   niter;
         float res;
         float resa[_KKNOISE];
         float resr[_KKNOISE];
    } output_pixel_residual;
    
    typedef struct output_segment_residual_{
         int   niter;
         float rest;  
         float resat[_KKNOISE];
         float resrt[_KKNOISE];
         output_pixel_residual pixel[_KIMAGE];  
    } output_segment_residual;

    typedef struct output_pixel_retr_par_{
         float par[_KPARS];
    } output_pixel_retr_par;
    
    typedef struct output_segment_retr_par_{
         int   ngrid;           // ???
         float radius[_NRR];     // ???
         float SDL[_NRC][_NRR];   // ??? 
         output_pixel_retr_par pixel[_KIMAGE];  
    } output_segment_retr_par;

    typedef struct otput_segment_fitting_{
        sensor_data_t segment_fit;
    }output_segment_fitting;
    
    typedef struct output_segment_retrieval_ {
         output_segment_residual res;  
         output_segment_retr_par par;  
         output_segment_fitting fit;
    } output_segment_retrieval;

    // Optical characteristics      

// ext,ssa,aext   - spectral ext, sca and aext for each aerosol component
// Aexp      - Angstrom exponent (wl(4)/wl(5))
      typedef struct output_pixel_opt_wl_ {
         float extt;
         float ssat;
         float aextt;
         float ext[_KSD];
         float ssa[_KSD];
         float aext[_KSD];
      } output_pixel_opt_wl;
      
      typedef struct output_pixel_opt_ {
         float Aexp;
         output_pixel_opt_wl wl[_KW];
      } output_pixel_opt;
      
      typedef struct output_segment_opt_ {
         output_pixel_opt pixel[_KIMAGE];
      } output_segment_opt;
      
// Refractive index
      typedef struct output_pixel_rindex_wl_ {
         float mreal[_KSD];
         float mimag[_KSD];
      } output_pixel_rindex_wl;
      
      typedef struct output_pixel_rindex_ {
         output_pixel_rindex_wl wl[_KW];
      } output_pixel_rindex;
      
      typedef struct output_segment_rindex_ {
         output_pixel_rindex pixel[_KIMAGE];
      } output_segment_rindex;

// Phase matrix :
      typedef struct output_pixel_ph_matrix_wl_ {
         float ph11[_KSD][_KMpar];
         float ph12[_KSD][_KMpar];
         float ph22[_KSD][_KMpar];
         float ph33[_KSD][_KMpar];
         float ph34[_KSD][_KMpar];
         float ph44[_KSD][_KMpar];
         
         float pht11[_KMpar];
         float pht12[_KMpar];
         float pht22[_KMpar];
         float pht33[_KMpar];
         float pht34[_KMpar];
         float pht44[_KMpar];         
         
      } output_pixel_ph_matrix_wl;
      
      typedef struct output_pixel_ph_matrix_ {   
         output_pixel_ph_matrix_wl wl[_KW];
      }output_pixel_ph_matrix;
      
      typedef struct output_segment_ph_matrix_ {
         int  nangle;
         float angle[_KMpar];
         output_pixel_ph_matrix pixel[_KIMAGE];
      } output_segment_ph_matrix;

//Lidar and depolarization ratios :
      
      typedef struct output_pixel_lidar_ratio_wl_ {
         float lr[_KSD];
         float ldpr[_KSD];
         float lrt;
         float ldprt;
      } output_pixel_lidar_ratio_wl;
      
      typedef struct output_pixel_lidar_ratio_ {
         output_pixel_lidar_ratio_wl wl[_KW];
      } output_pixel_lidar_ratio;
      
      typedef struct output_segment_lidar_ratio_ {
         output_pixel_lidar_ratio pixel[_KIMAGE];
      } output_segment_lidar_ratio;

// Chemistry parameters
      typedef struct output_pixel_chemistry_ {
         float rh[_KSD];
         float fwtr[_KSD];
         float fslbl[_KSD];
         float finslbl[_KSD];
         float fsoot[_KSD];
         float firon[_KSD];
         float fbrc[_KSD];   /* add by lei on 15/11/2016 */
      } output_pixel_chemistry;
              
      typedef struct output_segment_chemistry_ { 
         output_pixel_chemistry pixel[_KIMAGE];
      } output_segment_chemistry;

// Two mode aerosol characteristics
// Structure contains output data (microphysical parameters):
// 0 - total, 1 - fine mode, 2 - coarse mode 
// reff        - volume median radius
// std         - standard deviation
// cv          - concentration 
// rm          - median radius 
// ext         - ext each aerosol component
      typedef struct output_pixel_sd2m_mph_ {
         float cv[3];
         float std[3];
         float rm[3];
         float reff[3];
      } output_pixel_sd2m_mph; 
                 
      typedef struct output_segment_sd2m_mph_ {
         output_pixel_sd2m_mph pixel[_KIMAGE];
      }output_segment_sd2m_mph;

      typedef struct output_pixel_sd2m_opt_ {	
         float ext[2][_KW];
      } output_pixel_sd2m_opt;
      
      typedef struct output_segment_sd2m_opt_ {
         output_pixel_sd2m_opt pixel[_KIMAGE];
      } output_segment_sd2m_opt;

      typedef struct output_segment_sd2m_ {
         output_segment_sd2m_mph mph;
         output_segment_sd2m_opt opt;
      } output_segment_sd2m;
//-----------------------------------------------------------------------------------------
// AL Particulate Matter
      typedef struct output_pixel_PM_ {
         float  pm[2];
      } output_pixel_PM;
      typedef struct output_segment_PM_ {
         output_pixel_PM   pixel[_KIMAGE];
      } output_segment_PM;
// AL aerosol or cloud type
// index value — Aerosol type
// 0 – Complex mixture
// 1 – Background Aerosol
// 2 – Water/Maritime
// 3 — Urban Polluted
// 4 – Mixed aerosol
// 5 – Urban Clean
// 6 – Smoke Smoldering
// 7 – Smoke flaming
// 8 – Mineral dust
      typedef struct output_pixel_types_ {
         int  index;
      } output_pixel_types;
      typedef struct output_segment_types_ {
         output_pixel_types   pixel[_KIMAGE];
      } output_segment_types;
// Surface      
      
// Surface albedo 
      
      typedef struct output_pixel_surface_wl_{
         float salbedo;
      } output_pixel_surface_wl;
      
      typedef struct output_pixel_surface_{
         float ndvi;
         output_pixel_surface_wl wl[_KW];
      }output_pixel_surface;
      
      typedef struct output_segment_surface_ {
         output_pixel_surface pixel[_KIMAGE];
      } output_segment_surface;

// Radiative forcing

// Radiative broad band flux and forcing
      typedef struct output_pixel_bbflux_ {
         int nhlv;
         float bbufx0[_KNT];
         float bbdfx0[_KNT];
         float bbufxa[_KNT];
         float bbdfxa[_KNT];
         float hlv[_KNT];
      } output_pixel_bbflux;
  
      typedef struct output_segment_bbflux_ {
         output_pixel_bbflux pixel[_KIMAGE];
      } output_segment_bbflux;

      typedef struct output_pixel_forcing_ {
         int   nhlv;
         float netforc[_KNT];
         float forceff[_KNT];
         float hlv[_KNT];
      } output_pixel_forcing;
      
      typedef struct output_segment_forcing_ {
         output_pixel_forcing pixel[_KIMAGE];
      } output_segment_forcing;

      typedef struct output_segment_rad_forcing_ {
         output_segment_bbflux  bbflux;
         output_segment_forcing forcing;
      } output_segment_rad_forcing;

// Error estimations

// ERRP - Standard deviations of retrieved parameter logarithms (~relative errors)
// BIASP - Standard deviation of systematic errors of retrieved parameter logarithms
// ERR_ - Standard deviations of retrieved optical characteristic logarithms (~relative errors) 
// BIAS_ - Standard deviations of systematic errors of retrieved optical characteristic logarithms
// structure par      contains BIAS & ERR for all retrieved (of aerosol, clouds and surface)
// structure aerosol1 contains BIAS & ERR for ext & ssa - optical thickness and single scattering albedo of aerosol
// structure aerosol2 contains BIAS & ERR for lr        - lidar ratio of aerosol
// structure cloud1   contains BIAS & ERR for ext & ssa - optical thickness and single scattering albedo of clouds
// structure cloud2   contains BIAS & ERR for lr        - lidar ratio of clouds
      typedef struct output_pixel_err_estim_par_ {
         float ERRP[_KPARS];
         float BIASP[_KPARS];
      } output_pixel_err_estim_par;
  
      typedef struct output_segment_err_estim_par_ {
         output_pixel_err_estim_par pixel[_KIMAGE];
      } output_segment_err_estim_par;

      typedef struct output_pixel_err_estim_particles_opt_ {
         float ERR_ext[_KSD][_KW];
         float BIAS_ext[_KSD][_KW];
         float ERR_extt[_KW];
         float BIAS_extt[_KW];
         float ERR_ssa[_KSD][_KW];
         float BIAS_ssa[_KSD][_KW];
         float ERR_ssat[_KW];
         float BIAS_ssat[_KW];
      } output_pixel_err_estim_particles_opt;
      
      typedef struct output_segment_err_estim_particles_opt_ {
         output_pixel_err_estim_particles_opt pixel[_KIMAGE];
      } output_segment_err_estim_particles_opt;
      
      typedef struct output_pixel_err_estim_particles_lidar_ {
         float ERR_lr[_KSD][_KW];
         float BIAS_lr[_KSD][_KW];     
         float ERR_lrt[_KW];
         float BIAS_lrt[_KW];
      } output_pixel_err_estim_particles_lidar;
      
      typedef struct output_segment_err_estim_particles_lidar_ {
         output_pixel_err_estim_particles_lidar pixel[_KIMAGE];
      } output_segment_err_estim_particles_lidar;

      typedef struct output_segment_err_estim_particles_ { 
         output_segment_err_estim_particles_opt    opt;
         output_segment_err_estim_particles_lidar  lidar;
      } output_segment_err_estim_particles;

      typedef struct output_segment_err_estim_ { 
         output_segment_err_estim_par      par;
         output_segment_err_estim_particles aerosol ;  
         output_segment_err_estim_particles clouds  ;   
      } output_segment_err_estim;

      
// Aerosol
      
      typedef struct output_segment_particles_ {
         output_segment_opt        opt;
         output_segment_rindex     rind;
         output_segment_ph_matrix  phmx;
         output_segment_lidar_ratio lidar;
         output_segment_sd2m       sd2m;
         output_segment_chemistry  chem;
         output_segment_PM         pm;
         output_segment_types      types;
      } output_segment_particles;


// General output
      typedef struct output_segment_general_ {
          output_segment_products     products;
          output_segment_retrieval    retrieval;
          output_segment_particles    aerosol;
          output_segment_particles    clouds;
          output_segment_surface      surface;
          output_segment_err_estim    errest;
          output_segment_rad_forcing  forcing;
      } output_segment_general;
    
      
      
/**
 * Return array of retrieved parameters. It can be used within grasp_parameters 
 * library to extract easily information in retrieved parameters 
 * @param output Output of retrieval from a segment
 * @param ipix Index of pixel
 * @return array of retrieved parameters
 */      
const float *grasp_output_segment_parameters(const output_segment_general *output, int ipix);

/**
 * Return true if output retrieval residuals is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */        
bool grasp_output_segment_products_retrieval_res (const output_segment_general *output);

/**
 * Return true if output retrieval parameters is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_retrieval_par (const output_segment_general *output);

/**
 * Return true if output retrieval fitting is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_retrieval_fit (const output_segment_general *output);

/**
 * Return true if output for aerosol optical products is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_opt (const output_segment_general *output);

/**
 * Return true if output for aerosol refractive index is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_rind (const output_segment_general *output);

/**
 * Return true if output of aerosol chemistry is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_chem (const output_segment_general *output);

/**
 * Return true if output of aerosol phase matrix is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_phmx (const output_segment_general *output);

/**
 * Return true if output for aerosol lidar products is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_lidar (const output_segment_general *output);

/**
 * Return true if output products for simulated two modes is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_sd2m_mph (const output_segment_general *output);

/**
 * Return true if output extinction for aerosol two modes simulated is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_sd2m_ext (const output_segment_general *output);

/**
 * Return true if output aerosol particular matter is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_pm (const output_segment_general *output);

/**
 * Return true if output for aerosol types is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_aerosol_types (const output_segment_general *output);

/**
 * Return true if output for cloud optical properties is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_opt (const output_segment_general *output);

/**
 * Return true if output of cloud refractive index is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_rind (const output_segment_general *output);

/**
 * Return true if output cloud chemistry is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_chem (const output_segment_general *output);

/**
 * Return true if output cloud phase matrix is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_phmx (const output_segment_general *output);

/**
 * Return true if output cloud lidar is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_lidar (const output_segment_general *output);

/**
 * Return true if output clouds in two simulated modes is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_sd2m_mph (const output_segment_general *output);

/**
 * Return true if output extinction for clouds in two simulated modes is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_sd2m_ext (const output_segment_general *output);

/**
 * Return true if output for cloud particular matter is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_pm (const output_segment_general *output);

/**
 * Return true if cloud type is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_clouds_types (const output_segment_general *output);

/**
 * If surface products is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_surface_surf (const output_segment_general *output);

/**
 * Return true if output error estimation for parameters is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_errest_par (const output_segment_general *output);

/**
 * Return true if output error estimation for aerosol optical products is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_errest_aerosol_opt (const output_segment_general *output);

/**
 * Return true if output error estimation for aerosol lidar signal is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_errest_aerosol_lidar (const output_segment_general *output);

/**
 * Return true if output for error estimation of clouds optical properties is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_errest_clouds_opt (const output_segment_general *output);

/**
 * Return true if output for cloud lidar signal is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_errest_clouds_lidar (const output_segment_general *output);

/**
 * Return true if output for forcing flux is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_forcing_bbflux (const output_segment_general *output);

/**
 * Return true if output for forcing is available
 * @param output Output of retrieval from a segment
 * @return if it is available
 */
bool grasp_output_segment_products_forcing_forcing (const output_segment_general *output);

/**
 * Return number of iterations in multipixel scenario
 * @param output Output of retrieval from a segment
 * @return number of iterations in multipixel scenario
 */
int grasp_output_segment_retrieval_res_niter (const output_segment_general *output);

/**
 * Return total residual for multi-pixel retrieval
 * @param output Output of retrieval from a segment
 * @return total residual for multi-pixel retrieval
 */
float grasp_output_segment_retrieval_res_rest (const output_segment_general *output);

/**
 * Return detailed absolute measurement residuals for segment
 * @param output Output of retrieval from a segment 
 * @param inoise Number of noise
 * @return detailed absolute measurement residuals for segment
 */
float grasp_output_segment_retrieval_res_resat (const output_segment_general *output, int inoise);

/**
 * Return detailed relative measurement residuals for segment
 * @param output Output of retrieval from a segment
 * @param inoise Number of noise
 * @return detailed relative and relative measurement residuals for segment
 */
float grasp_output_segment_retrieval_res_resrt (const output_segment_general *output, int inoise);

/**
 * Return number of iteration in single pixel scenario
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return number of iteration in single pixel scenario
 */
int grasp_output_segment_retrieval_res_pixel_niter (const output_segment_general *output, int ipix);

/**
 * Return single pixel total residual (meas + smothness constrain )
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return single pixel total residual (meas + smothness constrain )
 */
float grasp_output_segment_retrieval_res_pixel_res (const output_segment_general *output, int ipix);

/**
 * Return detailed absolute measurement residuals for single pixel 
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param inoise Number of noise
 * @return detailed absolute measurement residuals for single pixel 
 */
float grasp_output_segment_retrieval_res_pixel_resa (const output_segment_general *output, int ipix, int inoise);

/**
 * Return detailed relative measurement residuals for single pixel 
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param inoise Number of noise
 * @return detailed relative measurement residuals for single pixel 
 */
float grasp_output_segment_retrieval_res_pixel_resr (const output_segment_general *output, int ipix, int inoise);

/**
 * Return number of grid radii for SD (saved to print output)
 * @param output Output of retrieval from a segment
 * @return number of grid radii for SD (saved to print output)
 */
int grasp_output_segment_retrieval_par_ngrid (const output_segment_general *output);

/**
 * Return grid radii (saved to print output)
 * @param output Output of retrieval from a segment
 * @param irr
 * @return grid radii (saved to print output)
 */
float grasp_output_segment_retrieval_par_radius (const output_segment_general *output, int irr);

/**
 * Return if applicable, precalculated lognormal bins SD (saved to print output)
 * @param output Output of retrieval from a segment
 * @param irr
 * @param irc
 * @return if applicable, precalculated lognormal bins SD (saved to print output)
 */
float grasp_output_segment_retrieval_par_sdl (const output_segment_general *output, int irr, int irc);

/**
 * Return retrieved aerosol and surface reflectance parameters
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param ipar
 * @return retrieved aerosol and surface reflectance parameters
 */
float grasp_output_segment_retrieval_par_parameters (const output_segment_general *output, int ipix, int ipar);

/**
 * Return segment fit
 * @param output Output of retrieval from a segment
 * @return Fit segment. It is like sdata structure but filled with fitting information
 */
const sensor_data_t *grasp_output_segment_retrieval_fit_segment_fit (const output_segment_general *output);

///////////////

/**
 * Return angstrom exponent for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return angstrom exponent for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_aexp (const output_segment_general *output, int ipix);

/**
 * Return total extinction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total extinction for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total single scattering albedo for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total single scattering albedo for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total absorption extinction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total absorption extinction for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_aextt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return extinction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return extinction for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return single scattering albedo for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return single scattering albedo for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return absorption extinction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return absorption extinction for aerosol optical properties
 */
float grasp_output_segment_aerosol_opt_aext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return real part refractive index for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return real part refractive index for aerosol optical properties
 */
float grasp_output_segment_aerosol_rind_mreal (const output_segment_general *output, int ipix, int iwl, int isd);


/**
 * Return imaginary part refractive index for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return imaginary part refractive index for aerosol optical properties
 */
float grasp_output_segment_aerosol_rind_mimag (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return number of angles of phase matrix for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @return number of angles of phase matrix for aerosol optical properties
 */
int grasp_output_segment_aerosol_phmx_nangle (const output_segment_general *output);

/**
 * Return angles of the phase matrix for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param impar
 * @return angles of the phase matrix for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_angle (const output_segment_general *output, int impar);

/**
 * Return p11 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p11 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph11 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p12 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p12 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph12 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p22 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p22 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph22 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p33 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p33 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph33 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p34 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p34 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph34 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p44 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p44 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_ph44 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return total p11 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p11 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht11 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p12 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p12 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht12 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p22 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p22 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht22 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p33 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p33 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht33 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p34 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p34 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht34 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p44 for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p44 for aerosol optical properties
 */
float grasp_output_segment_aerosol_phmx_pht44 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return lidar ratio for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return lidar ratio for aerosol optical properties
 */
float grasp_output_segment_aerosol_lidar_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return lidar depolarization profile for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return lidar depolarization profile for aerosol optical properties
 */
float grasp_output_segment_aerosol_lidar_ldpr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return total lidar ratio for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total lidar ratio for aerosol optical properties
 */
float grasp_output_segment_aerosol_lidar_lrt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total lidar depolarization for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total lidar depolarization for aerosol optical properties
 */
float grasp_output_segment_aerosol_lidar_ldprt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return concentration for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_cv (const output_segment_general *output, int ipix, int i);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_std (const output_segment_general *output, int ipix, int i);

/**
 * Return XXX for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_rm (const output_segment_general *output, int ipix, int i);

/**
 * Return refractive index for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_reff (const output_segment_general *output, int ipix, int i);

/**
 * Return extinction for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param i mode. 0 for fine and 1 for coarse
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_opt_ext (const output_segment_general *output, int ipix, int iwl, int i);

/**
 * Return concentration for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_cv with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_cv_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_std with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_std_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return XXX for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_rm with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_rm_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return refractive index for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_reff with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_reff_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return extinction for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_opt_ext with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_opt_ext_fine_mode (const output_segment_general *output, int ipix, int iwl);

/**
 * Return concentration for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_cv with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_cv_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_std with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_std_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return XXX for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_rm with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_rm_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return refractive index for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_mph_reff with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_mph_reff_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return extinction for two simulated modes for aerosol optical properties. Alias to grasp_output_segment_aerosol_sd2m_opt_ext with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_segment_aerosol_sd2m_opt_ext_coarse_mode (const output_segment_general *output, int ipix, int iwl);

/**
 * Return relative humidity for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return relative humidity for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_rh (const output_segment_general *output, int ipix, int isd);

/**
 * Return water fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return water fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_fwtr (const output_segment_general *output, int ipix, int isd);

/**
 * Return soluble fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return soluble fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_fslbl (const output_segment_general *output, int ipix, int isd);

/**
 * Return insoluble fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return  insoluble fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_finslbl (const output_segment_general *output, int ipix, int isd);

/**
 * Return soot fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return soot fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_fsoot (const output_segment_general *output, int ipix, int isd);

/**
 * Return iron fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return iron fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_firon (const output_segment_general *output, int ipix, int isd);

    
/**  add by lei on 15/11/2016
 * Return brc fraction for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of aerosol component
 * @return brc fraction for aerosol optical properties
 */
float grasp_output_segment_aerosol_chem_fbrc (const output_segment_general *output, int ipix, int isd);




/**
 * Return particular matter for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i
 * @return particular matter for aerosol optical properties
 */
float grasp_output_segment_aerosol_pm_pm (const output_segment_general *output, int ipix, int i);

/**
 * Return aerosol type for aerosol optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return aerosol type: 0 – Complex mixture; 1 – Background Aerosol; 2 – Water/Maritime; 3 — Urban Polluted; 4 – Mixed aerosol;  5 – Urban Clean;  6 – Smoke Smoldering;  7 – Smoke flaming;  8 – Mineral dust
 */
int grasp_output_segment_aerosol_types_index (const output_segment_general *output, int ipix);


///////////////

/**
 * Return angstrom exponent for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return angstrom exponent for clouds optical properties
 */
float grasp_output_segment_clouds_opt_aexp (const output_segment_general *output, int ipix);

/**
 * Return total extinction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total extinction for clouds optical properties
 */
float grasp_output_segment_clouds_opt_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total single scattering albedo for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total single scattering albedo for clouds optical properties
 */
float grasp_output_segment_clouds_opt_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total absorption extinction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total absorption extinction for clouds optical properties
 */
float grasp_output_segment_clouds_opt_aextt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return extinction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return extinction for clouds optical properties
 */
float grasp_output_segment_clouds_opt_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return single scattering albedo for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return single scattering albedo for clouds optical properties
 */
float grasp_output_segment_clouds_opt_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return absorption extinction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return absorption extinction for clouds optical properties
 */
float grasp_output_segment_clouds_opt_aext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return real part refractive index for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return real part refractive index for clouds optical properties
 */
float grasp_output_segment_clouds_rind_mreal (const output_segment_general *output, int ipix, int iwl, int isd);


/**
 * Return imaginary part refractive index for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return imaginary part refractive index for clouds optical properties
 */
float grasp_output_segment_clouds_rind_mimag (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return number of angles of phase matrix for clouds optical properties
 * @param output Output of retrieval from a segment
 * @return number of angles of phase matrix for clouds optical properties
 */
int grasp_output_segment_clouds_phmx_nangle (const output_segment_general *output);

/**
 * Return angles of the phase matrix for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param impar
 * @return angles of the phase matrix for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_angle (const output_segment_general *output, int impar);

/**
 * Return p11 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p11 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph11 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p12 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p12 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph12 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p22 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p22 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph22 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p33 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p33 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph33 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p34 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p34 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph34 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return p44 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p44 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_ph44 (const output_segment_general *output, int ipix, int iwl, int impar, int isd);

/**
 * Return total p11 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p11 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht11 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p12 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p12 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht12 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p22 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p22 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht22 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p33 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p33 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht33 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p34 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p34 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht34 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return total p44 for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p44 for clouds optical properties
 */
float grasp_output_segment_clouds_phmx_pht44 (const output_segment_general *output, int ipix, int iwl, int impar);

/**
 * Return lidar ratio for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return lidar ratio for clouds optical properties
 */
float grasp_output_segment_clouds_lidar_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return lidar depolarization profile for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return lidar depolarization profile for clouds optical properties
 */
float grasp_output_segment_clouds_lidar_ldpr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return total lidar ratio for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total lidar ratio for clouds optical properties
 */
float grasp_output_segment_clouds_lidar_lrt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return total lidar depolarization for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return total lidar depolarization for clouds optical properties
 */
float grasp_output_segment_clouds_lidar_ldprt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return concentration for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_cv (const output_segment_general *output, int ipix, int i);

/**
 * Return standard deviation for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_std (const output_segment_general *output, int ipix, int i);

/**
 * Return XXX for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_rm (const output_segment_general *output, int ipix, int i);

/**
 * Return refractive index for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i mode. 0 for fine and 1 for coarse
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_reff (const output_segment_general *output, int ipix, int i);

/**
 * Return extinction for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param i mode. 0 for fine and 1 for coarse
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_opt_ext (const output_segment_general *output, int ipix, int iwl, int i);

/**
 * Return concentration for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_cvfine_mode with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_cv_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return standard deviation for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_std with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_std_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return XXX for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_std with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_rm_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return refractive index for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_reff with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_reff_fine_mode (const output_segment_general *output, int ipix);

/**
 * Return extinction for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_opt_ext with i = 0
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_opt_ext_fine_mode (const output_segment_general *output, int ipix, int iwl);

/**
 * Return concentration for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_cv with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_cv_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return standard deviation for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_std with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_std_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return XXX for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_rm with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_rm_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return refractive index for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_mph_reff with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_mph_reff_coarse_mode (const output_segment_general *output, int ipix);

/**
 * Return extinction for two simulated modes for clouds optical properties. Alias to grasp_output_segment_clouds_sd2m_opt_ext with i = 1
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_segment_clouds_sd2m_opt_ext_coarse_mode (const output_segment_general *output, int ipix, int iwl);

/**
 * Return relative humidity for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return relative humidity for clouds optical properties
 */
float grasp_output_segment_clouds_chem_rh (const output_segment_general *output, int ipix, int isd);

/**
 * Return water fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return water fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_fwtr (const output_segment_general *output, int ipix, int isd);

/**
 * Return soluble fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return soluble fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_fslbl (const output_segment_general *output, int ipix, int isd);

/**
 * Return insoluble fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return  insoluble fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_finslbl (const output_segment_general *output, int ipix, int isd);

/**
 * Return soot fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return soot fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_fsoot (const output_segment_general *output, int ipix, int isd);

/**
 * Return iron fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return iron fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_firon (const output_segment_general *output, int ipix, int isd);

/**add by lei on 15/11/2016 
 * Return brc fraction for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param isd Index of clouds component
 * @return brc fraction for clouds optical properties
 */
float grasp_output_segment_clouds_chem_fbrc (const output_segment_general *output, int ipix, int isd);

/**
 * Return particular matter for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param i
 * @return particular matter for clouds optical properties
 */

float grasp_output_segment_clouds_pm_pm (const output_segment_general *output, int ipix, int i);

/**
 * Return clouds type for clouds optical properties
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return clouds type: 0 – Complex mixture; 1 – Background Aerosol; 2 – Water/Maritime; 3 — Urban Polluted; 4 – Mixed clouds;  5 – Urban Clean;  6 – Smoke Smoldering;  7 – Smoke flaming;  8 – Mineral dust
 */
int grasp_output_segment_clouds_types_index (const output_segment_general *output, int ipix);


///////////////

/**
 * Return surface ndvi
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return surface ndvi
 */
float grasp_output_segment_surface_ndvi (const output_segment_general *output, int ipix);

/**
 * Return surface albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return surface albedo
 */
float grasp_output_segment_surface_salbedo (const output_segment_general *output, int ipix, int iwl);

/////////////////

/**
 * Return error estimation for parameter
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param ipar
 * @return error estimation for parameter
 */
float grasp_output_segment_errest_par_errp (const output_segment_general *output, int ipix, int ipar);

/**
 * Return bias of error estimation of parameters
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param ipar
 * @return bias of error estimation of parameters
 */
float grasp_output_segment_errest_par_biasp (const output_segment_general *output, int ipix, int ipar);

/**
 * Return error estimation of extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of extinction
 */
float grasp_output_segment_errest_aerosol_opt_err_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of extinction
 */
float grasp_output_segment_errest_aerosol_opt_bias_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total extinction
 */
float grasp_output_segment_errest_aerosol_opt_err_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total extinction
 */
float grasp_output_segment_errest_aerosol_opt_bias_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return error estimation of single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of single scattering albedo
 */
float grasp_output_segment_errest_aerosol_opt_err_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of single scattering albedo
 */
float grasp_output_segment_errest_aerosol_opt_bias_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total single scattering albedo
 */
float grasp_output_segment_errest_aerosol_opt_err_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total single scattering albedo 
 */
float grasp_output_segment_errest_aerosol_opt_bias_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return error estimation of lidar ration
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of lidar ration
 */
float grasp_output_segment_errest_aerosol_lidar_err_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of lidar ratio
 */
float grasp_output_segment_errest_aerosol_lidar_bias_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total lidar ratio
 */
float grasp_output_segment_errest_aerosol_lidar_err_lrt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total lidar ratio
 */
float grasp_output_segment_errest_aerosol_lidar_bias_lrt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return error estimation of extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of extinction
 */
float grasp_output_segment_errest_clouds_opt_err_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of extinction
 */
float grasp_output_segment_errest_clouds_opt_bias_ext (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total extinction
 */
float grasp_output_segment_errest_clouds_opt_err_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total extinction
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total extinction
 */
float grasp_output_segment_errest_clouds_opt_bias_extt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return error estimation of single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of single scattering albedo
 */
float grasp_output_segment_errest_clouds_opt_err_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of single scattering albedo
 */
float grasp_output_segment_errest_clouds_opt_bias_ssa (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total single scattering albedo
 */
float grasp_output_segment_errest_clouds_opt_err_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total single scattering albedo
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total single scattering albedo 
 */
float grasp_output_segment_errest_clouds_opt_bias_ssat (const output_segment_general *output, int ipix, int iwl);

/**
 * Return error estimation of lidar ration
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of lidar ration
 */
float grasp_output_segment_errest_clouds_lidar_err_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return bias of error estimation of lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of lidar ratio
 */
float grasp_output_segment_errest_clouds_lidar_bias_lr (const output_segment_general *output, int ipix, int iwl, int isd);

/**
 * Return error estimation of total lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return error estimation of total lidar ratio
 */
float grasp_output_segment_errest_clouds_lidar_err_lrt (const output_segment_general *output, int ipix, int iwl);

/**
 * Return bias of error estimation of total lidar ratio
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total lidar ratio
 */
float grasp_output_segment_errest_clouds_lidar_bias_lrt (const output_segment_general *output, int ipix, int iwl);

///////////////

/**
 * Return number of heights of forcing fluxes
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return number of heights of forcing fluxes
 */
int grasp_output_segment_forcing_bbflux_nhlv (const output_segment_general *output, int ipix);

/**
 * Return broad band up-ward flux without aerosol at each height
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return broad band up-ward flux without aerosol at each height
 */
float grasp_output_segment_forcing_bbflux_bbufx0 (const output_segment_general *output, int ipix, int iknt);

/**
 * Return broad band down-ward flux without aerosol at each height
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return broad band down-ward flux without aerosol at each height
 */
float grasp_output_segment_forcing_bbflux_bbdfx0 (const output_segment_general *output, int ipix, int iknt);

/**
 * Return broad band up-ward flux with aerosol at each height
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return broad band up-ward flux with aerosol at each height
 */
float grasp_output_segment_forcing_bbflux_bbufxa (const output_segment_general *output, int ipix, int iknt);

/**
 * Return broad band down-ward flux with aerosol at each height
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return broad band down-ward flux with aerosol at each height
 */
float grasp_output_segment_forcing_bbflux_bbdfxa (const output_segment_general *output, int ipix, int iknt);

/**
 * Return heights of forcing fluxes
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return heights of forcing fluxes
 */
float grasp_output_segment_forcing_bbflux_hlv (const output_segment_general *output, int ipix, int iknt);

/**
 * Return number of heights of forcing calculations
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @return number of heights of forcing calculations
 */
int grasp_output_segment_forcing_forcing_nhlv (const output_segment_general *output, int ipix);

/**
 * Return net forcing
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return net forcing
 */
float grasp_output_segment_forcing_forcing_netforc (const output_segment_general *output, int ipix, int iknt);

/**
 * Return forcing efficiency
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return forcing efficiency
 */
float grasp_output_segment_forcing_forcing_forceff (const output_segment_general *output, int ipix, int iknt);

/**
 * Return heights of forcing calculations
 * @param output Output of retrieval from a segment
 * @param ipix Number of pixel inside the segment
 * @param iknt
 * @return heights of forcing calculations
 */
float grasp_output_segment_forcing_forcing_hlv (const output_segment_general *output, int ipix, int iknt);

      
      
#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_SEGMENT_RESULT_H */

