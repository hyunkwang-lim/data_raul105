/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_output_tile_result.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 15:06
 */

#ifndef GRASP_OUTPUT_TILE_RESULT_H
#define	GRASP_OUTPUT_TILE_RESULT_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <inttypes.h>
//#include "mod_par_inv.inc"
//#include "mod_par_OS.inc"
#include "mod_par_DLS.inc"
#include "mod_par_DLS_bin.inc"
#include "../settings/grasp_products.h"
#include "../input/grasp_input_segment.h"
    
  typedef struct grasp_output_tile_pixel_information_ {
        // identification
        int segment_time; // segment position
        int segment_col; 
        int segment_row; 
        int it;           // index inside segment
        int ix;
        int iy;
        int out_x;        // index in output grid
        int out_y;
        int out_t; 
        float latitude;
        float longitude;
        int grid_col;    // index in input grid.
        int grid_row;
        int64_t time;    
        // content information
        int nwl;   
        int cloud_flag;
        float land_percent;
        int file_index;
        float masl;
        float *sza; /**< @brief [wln] */
  }grasp_output_tile_pixel_information;  
  
  typedef struct grasp_output_tile_retrieval_res_{
      int niter;
      float rest;
      float *resa; //[nnoises]
      float *resr; //[nnoises]
  }grasp_output_tile_retrieval_res;
  
  typedef struct grasp_output_tile_retrieval_par_{
      float *parameters; //[npars]
  }grasp_output_tile_retrieval_par;
  
  typedef struct grasp_output_tile_retrieval_fit_{
      pixel_t original_pixel;
      pixel_t fitting_pixel;
  }grasp_output_tile_retrieval_fit;
  
  typedef struct grasp_output_tile_aerosol_opt_{
      float Aexp;
      float *extt;  /**< @brief [wln] */   
      float *ssat;  /**< @brief [wln] */  
      float *aextt; /**< @brief [wln] */  
      float *ext;   /**< @brief [wln][nsd] */  
      float *ssa;   /**< @brief [wln][nsd] */  
      float *aext;  /**< @brief [wln][nsd] */  
  }grasp_output_tile_aerosol_opt;
  
  typedef struct grasp_output_tile_aerosol_rind_{
      float *mreal; /**< @brief [wln][nsd] */ 
      float *mimag; /**< @brief [wln][nsd] */ 
  }grasp_output_tile_aerosol_rind;
  
  typedef struct grasp_output_tile_aerosol_phmx_{
      float *ph11; /**< @brief [wln][nsd][nmpar] */ 
      float *ph12; /**< @brief [wln][nsd][nmpar] */ 
      float *ph22; /**< @brief [wln][nsd][nmpar] */ 
      float *ph33; /**< @brief [wln][nsd][nmpar] */ 
      float *ph34; /**< @brief [wln][nsd][nmpar] */ 
      float *ph44; /**< @brief [wln][nsd][nmpar] */ 
      float *pht11; /**< @brief [wln][nmpar] */ 
      float *pht12; /**< @brief [wln][nmpar] */ 
      float *pht22; /**< @brief [wln][nmpar] */ 
      float *pht33; /**< @brief [wln][nmpar] */ 
      float *pht34; /**< @brief [wln][nmpar] */ 
      float *pht44; /**< @brief [wln][nmpar] */ 
  }grasp_output_tile_aerosol_phmx;
  
  typedef struct grasp_output_tile_aerosol_lidar_{
      float *lr;    /**< @brief [wln][nsd] */
      float *ldpr;  /**< @brief [wln][nsd] */
      float *lrt;   /**< @brief [wln] */
      float *ldprt; /**< @brief [wln] */
  }grasp_output_tile_aerosol_lidar;
  
  typedef struct grasp_output_tile_aerosol_sd2m_mph_{
      float cv[3];   /**< @brief [3] */
      float std[3];  /**< @brief [3] */
      float rm[3];   /**< @brief [3] */
      float reff[3]; /**< @brief [3] */
  }grasp_output_tile_aerosol_sd2m_mph;  
  
  typedef struct grasp_output_tile_aerosol_sd2m_ext_{
      float *ext;   /**< @brief [wln][2] */
  }grasp_output_tile_aerosol_sd2m_ext;
  
  typedef struct grasp_output_tile_aerosol_chem_{
      float *rh;      /**< @brief [nsd] */
      float *fwrt;    /**< @brief [nsd] */
      float *fslbl;   /**< @brief [nsd] */
      float *finslbl; /**< @brief [nsd] */
      float *fsoot;   /**< @brief [nsd] */
      float *firon;   /**< @brief [nsd] */
      float *fbrc;   /**< @brief [nsd] */ /* add by lei on 15/11/2016 */
  }grasp_output_tile_aerosol_chem;
  
  typedef struct grasp_output_tile_aerosol_pm_{
      float *pm;  /**< @brief [nPM_diam] */
  }grasp_output_tile_aerosol_pm;
  
  typedef struct grasp_output_tile_aerosol_types_{
      int index;
  }grasp_output_tile_aerosol_types;
  //
  
  typedef struct grasp_output_tile_clouds_opt_{
      float Aexp;
      float *extt;  /**< @brief [wln] */   
      float *ssat;  /**< @brief [wln] */  
      float *aextt; /**< @brief [wln] */  
      float *ext;   /**< @brief [wln][nsd] */  
      float *ssa;   /**< @brief [wln][nsd] */  
      float *aext;  /**< @brief [wln][nsd] */  
  }grasp_output_tile_clouds_opt;
  
  typedef struct grasp_output_tile_clouds_rind_{
      float *mreal; /**< @brief [wln][nsd] */ 
      float *mimag; /**< @brief [wln][nsd] */ 
  }grasp_output_tile_clouds_rind;
  
  typedef struct grasp_output_tile_clouds_phmx_{
      float *ph11; /**< @brief [wln][nsd][nmpar] */ 
      float *ph12; /**< @brief [wln][nsd][nmpar] */ 
      float *ph22; /**< @brief [wln][nsd][nmpar] */ 
      float *ph33; /**< @brief [wln][nsd][nmpar] */ 
      float *ph34; /**< @brief [wln][nsd][nmpar] */ 
      float *ph44; /**< @brief [wln][nsd][nmpar] */ 
      float *pht11; /**< @brief [wln][nmpar] */ 
      float *pht12; /**< @brief [wln][nmpar] */ 
      float *pht22; /**< @brief [wln][nmpar] */ 
      float *pht33; /**< @brief [wln][nmpar] */ 
      float *pht34; /**< @brief [wln][nmpar] */ 
      float *pht44; /**< @brief [wln][nmpar] */ 
  }grasp_output_tile_clouds_phmx;
  
  typedef struct grasp_output_tile_clouds_lidar_{
      float *lr;    /**< @brief [wln][nsd] */
      float *ldpr;  /**< @brief [wln][nsd] */
      float *lrt;   /**< @brief [wln] */
      float *ldprt; /**< @brief [wln] */
  }grasp_output_tile_clouds_lidar;
  
  typedef struct grasp_output_tile_clouds_sd2m_mph_{
      float cv[3];   /**< @brief [3] */
      float std[3];  /**< @brief [3] */
      float rm[3];   /**< @brief [3] */
      float reff[3]; /**< @brief [3] */
  }grasp_output_tile_clouds_sd2m_mph;  
  
  typedef struct grasp_output_tile_clouds_sd2m_ext_{
      float *ext;   /**< @brief [wln][2] */
  }grasp_output_tile_clouds_sd2m_ext;
  
  typedef struct grasp_output_tile_clouds_chem_{
      float *rh;      /**< @brief [nsd] */
      float *fwrt;    /**< @brief [nsd] */
      float *fslbl;   /**< @brief [nsd] */
      float *finslbl; /**< @brief [nsd] */
      float *fsoot;   /**< @brief [nsd] */
      float *firon;   /**< @brief [nsd] */
      float *fbrc;   /**< @brief [nsd] */ /* add by lei on 15/11/2016 */
  }grasp_output_tile_clouds_chem;
  
  typedef struct grasp_output_tile_clouds_pm_{
      float *pm;  /**< @brief [nPM_diam] */
  }grasp_output_tile_clouds_pm;
  
  typedef struct grasp_output_tile_clouds_types_{
      int index;
  }grasp_output_tile_clouds_types;
  
  //
  
  typedef struct grasp_output_tile_surface_surf_{
      float ndvi;
      float *salbedo; /**< @brief [nwl] */
  }grasp_output_tile_surface_surf;
  
  typedef struct grasp_output_tile_errest_par_{
      float *ERRP; /**< @brief [npars] */
      float *BIASP; /**< @brief [npars] */
  }grasp_output_tile_errest_par;
  
  typedef struct grasp_output_tile_errest_aerosol_aerosol_opt_{
      float *ERR_ext; /**< @brief [nsd][nwl] */
      float *BIAS_ext; /**< @brief [nsd][nwl] */
      float *ERR_extt; /**< @brief [nwl] */
      float *BIAS_extt; /**< @brief [nwl] */
      float *ERR_ssa; /**< @brief [nsd][nwl] */
      float *BIAS_ssa; /**< @brief [nsd][nwl] */
      float *ERR_ssat; /**< @brief [nwl] */
      float *BIAS_ssat;  /**< @brief [nwl] */     
  }grasp_output_tile_errest_aerosol_opt;  
  
  typedef struct grasp_output_tile_errest_aerosol_aerosol_lidar_{
      float *ERR_lr; /**< @brief [nsd][nwl] */
      float *BIAS_lr; /**< @brief [nsd][nwl] */
      float *ERR_lrt; /**< @brief [nwl] */
      float *BIAS_lrt; /**< @brief [nwl] */
  }grasp_output_tile_errest_aerosol_lidar;
  
  typedef struct grasp_output_tile_errest_clouds_opt_{
      float *ERR_ext; /**< @brief [nsd][nwl] */
      float *BIAS_ext; /**< @brief [nsd][nwl] */
      float *ERR_extt; /**< @brief [nwl] */
      float *BIAS_extt; /**< @brief [nwl] */
      float *ERR_ssa; /**< @brief [nsd][nwl] */
      float *BIAS_ssa; /**< @brief [nsd][nwl] */
      float *ERR_ssat; /**< @brief [nwl] */
      float *BIAS_ssat; /**< @brief [nwl] */
  }grasp_output_tile_errest_clouds_opt;
  
  typedef struct grasp_output_tile_errest_clouds_lidar_{
      float *ERR_lr; /**< @brief [nsd][nwl] */
      float *BIAS_lr; /**< @brief [nsd][nwl] */
      float *ERR_lrt; /**< @brief [nwl] */
      float *BIAS_lrt; /**< @brief [nwl] */
  }grasp_output_tile_errest_clouds_lidar;
 
  typedef struct grasp_output_tile_forcing_bbflux_{
      int nhlv;
      float *bbufx0;  /**< @brief [nhlv] */
      float *bbdfx0; /**< @brief [nhlv] */
      float *bbufxa; /**< @brief [nhlv] */
      float *bbdfxa; /**< @brief [nhlv] */
      float *hlv; /**< @brief [nhlv] */
  }grasp_output_tile_forcing_bbflux;
 
  
  typedef struct grasp_output_tile_forcing_forcing_{
        int   nhlv;
        float *netforc;   /**< @brief [nhlv] */
        float *forceff;   /**< @brief [nhlv] */
        float *hlv; /**< @brief [nhlv] */
  }grasp_output_tile_forcing_forcing;
  
  typedef struct pixel_result_t_{
      grasp_output_tile_pixel_information    information;
      grasp_output_tile_retrieval_res        *retrieval_res;
      grasp_output_tile_retrieval_par        *retrieval_par;
      grasp_output_tile_retrieval_fit        *retrieval_fit;
      grasp_output_tile_aerosol_opt          *aerosol_opt;
      grasp_output_tile_aerosol_rind         *aerosol_rind;
      grasp_output_tile_aerosol_phmx         *aerosol_phmx;
      grasp_output_tile_aerosol_lidar        *aerosol_lidar;
      grasp_output_tile_aerosol_sd2m_mph     *aerosol_sd2m_mph;
      grasp_output_tile_aerosol_sd2m_ext     *aerosol_sd2m_ext;
      grasp_output_tile_aerosol_chem         *aerosol_chem;
      grasp_output_tile_aerosol_pm           *aerosol_pm;
      grasp_output_tile_aerosol_types        *aerosol_types;
      grasp_output_tile_clouds_opt           *clouds_opt;
      grasp_output_tile_clouds_rind          *clouds_rind;
      grasp_output_tile_clouds_phmx          *clouds_phmx;
      grasp_output_tile_clouds_lidar         *clouds_lidar;
      grasp_output_tile_clouds_sd2m_mph      *clouds_sd2m_mph;
      grasp_output_tile_clouds_sd2m_ext      *clouds_sd2m_ext;
      grasp_output_tile_clouds_chem          *clouds_chem;
      grasp_output_tile_clouds_pm            *clouds_pm;
      grasp_output_tile_clouds_types         *clouds_types;
      grasp_output_tile_surface_surf         *surface_surf;
      grasp_output_tile_errest_par           *errest_par;
      grasp_output_tile_errest_aerosol_opt   *errest_aerosol_opt;
      grasp_output_tile_errest_aerosol_lidar *errest_aerosol_lidar;
      grasp_output_tile_errest_clouds_opt    *errest_clouds_opt;
      grasp_output_tile_errest_clouds_lidar  *errest_clouds_lidar;
      grasp_output_tile_forcing_bbflux       *forcing_bbflux;
      grasp_output_tile_forcing_forcing      *forcing_forcing;
    }pixel_result_t;

    typedef struct grasp_output_tile_segment_result_t_{
        int npixel;
        pixel_result_t *pixel_result;
    }grasp_output_tile_segment_result_t;      

    typedef struct grasp_output_tile_retrieval_information_{
        int tile_npixels; // Total number of pixels
        int tile_npixels_t;
        int tile_npixels_x;
        int tile_npixels_y;
        int npars;
        int nrr;
        int nrc;
        int nmpar;
        int nsd;
        int retrieval_par_ngrid;
        float retrieval_par_radius[_NRR];
        float retrieval_par_SDL[_NRC][_NRR];
        int phmx_angle[_KMpar];
        int nPM_diam;
        int nnoises;
    } grasp_output_tile_retrieval_information;
    
    typedef struct grasp_results_t_{
        output_segment_products products; // Available products in result
        grasp_output_tile_segment_result_t *segment_result;   // dimensions [t, x, y ] [itime,icol,irow]
        grasp_output_tile_retrieval_information information;
        pixel_result_t **tile_result_map;  // dimensions [t, x, y ] [itime,icol,irow]. results reindexed    
    }grasp_results_t;     

    
/**
 * 
 * @param output Output of retrieval from a tile
 * @param it index t
 * @param ix index x
 * @param iy index y
 * @return output Output of retrieval from a tile
 */    
bool grasp_output_tile_is_pixel(const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return true if output retrieval residuals is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */        
bool grasp_output_tile_products_retrieval_res (const grasp_results_t *output);

/**
 * Return true if output retrieval parameters is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_retrieval_par (const grasp_results_t *output);

/**
 * Return true if output retrieval fitting is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_retrieval_fit (const grasp_results_t *output);

/**
 * Return true if output for aerosol optical products is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_opt (const grasp_results_t *output);

/**
 * Return true if output for aerosol refractive index is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_rind (const grasp_results_t *output);

/**
 * Return true if output of aerosol chemistry is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_chem (const grasp_results_t *output);

/**
 * Return true if output of aerosol phase matrix is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_phmx (const grasp_results_t *output);

/**
 * Return true if output for aerosol lidar products is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_lidar (const grasp_results_t *output);

/**
 * Return true if output products for simulated two modes is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_sd2m_mph (const grasp_results_t *output);

/**
 * Return true if output extinction for aerosol two modes simulated is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_sd2m_ext (const grasp_results_t *output);

/**
 * Return true if output aerosol particular matter is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_pm (const grasp_results_t *output);

/**
 * Return true if output for aerosol types is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_aerosol_types (const grasp_results_t *output);

/**
 * Return true if output for cloud optical properties is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_opt (const grasp_results_t *output);

/**
 * Return true if output of cloud refractive index is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_rind (const grasp_results_t *output);

/**
 * Return true if output cloud chemistry is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_chem (const grasp_results_t *output);

/**
 * Return true if output cloud phase matrix is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_phmx (const grasp_results_t *output);

/**
 * Return true if output cloud lidar is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_lidar (const grasp_results_t *output);

/**
 * Return true if output clouds in two simulated modes is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_sd2m_mph (const grasp_results_t *output);

/**
 * Return true if output extinction for clouds in two simulated modes is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_sd2m_ext (const grasp_results_t *output);

/**
 * Return true if output for cloud particular matter is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_pm (const grasp_results_t *output);

/**
 * Return true if cloud type is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_clouds_types (const grasp_results_t *output);

/**
 * If surface products is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_surface_surf (const grasp_results_t *output);

/**
 * Return true if output error estimation for parameters is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_errest_par (const grasp_results_t *output);

/**
 * Return true if output error estimation for aerosol optical products is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_errest_aerosol_opt (const grasp_results_t *output);

/**
 * Return true if output error estimation for aerosol lidar signal is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_errest_aerosol_lidar (const grasp_results_t *output);

/**
 * Return true if output for error estimation of clouds optical properties is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_errest_clouds_opt (const grasp_results_t *output);

/**
 * Return true if output for cloud lidar signal is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_errest_clouds_lidar (const grasp_results_t *output);

/**
 * Return true if output for forcing flux is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_forcing_bbflux (const grasp_results_t *output);

/**
 * Return true if output for forcing is available
 * @param output Output of retrieval from a tile
 * @return if it is available
 */
bool grasp_output_tile_products_forcing_forcing (const grasp_results_t *output);

/**
 * Return total number of pixels availables in output
 * @param output Output of retrieval from a tile
 * @return total number of pixels availables in output
 */
int grasp_output_tile_information_tile_npixels(const grasp_results_t *output);

/**
 * Return total number of times in the tile
 * @param output Output of retrieval from a tile
 * @return total number of times in the tile
 */
int grasp_output_tile_information_tile_npixels_t(const grasp_results_t *output);

/**
 * Return width of the tile
 * @param output Output of retrieval from a tile
 * @return width of the tile
 */
int grasp_output_tile_information_tile_npixels_x(const grasp_results_t *output);

/**
 * Return height of the tile
 * @param output Output of retrieval from a tile
 * @return height of the tile
 */
int grasp_output_tile_information_tile_npixels_y(const grasp_results_t *output);

/**
 * Return number of total retrieved parameter for each pixel
 * @param output Output of retrieval from a tile
 * @return umber of total retrieved parameter for each pixel
 */
int grasp_output_tile_information_npars(const grasp_results_t *output);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @return 
 */
int grasp_output_tile_information_nrr(const grasp_results_t *output);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @return 
 */
int grasp_output_tile_information_nrc(const grasp_results_t *output);

/**
 * Return number of valid angles in output phase matrix
 * @param output Output of retrieval from a tile
 * @return number of valid angles in output phase matrix
 */
int grasp_output_tile_information_nmpar(const grasp_results_t *output);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @return 
 */
int grasp_output_tile_information_nsd(const grasp_results_t *output);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @return 
 */
float grasp_output_tile_information_retrieval_par_ridius(const grasp_results_t *output);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @param irr Index of rr of SDL
 * @param irc Inxed of rc of SDL
 * @return 
 */
float grasp_output_tile_information_retrieval_par_SDL(const grasp_results_t *output, int irr, int irc);

/**
 * 
 * @param output Output of retrieval from a tile
 * @return 
 */
int grasp_output_tile_information_npm_diam(const grasp_results_t *output);

/**
 * Return a specific angle (indexed by iangle) of output phase function
 * @param iangle Index of phase matrix angle
 * @param output Output of retrieval from a tile
 * @return a specific angle (indexed by iangle) of output phase function
 */
int grasp_output_tile_information_phmx_angle(const grasp_results_t *output, int iangle);

/**
 * Return 
 * @param output Output of retrieval from a tile
 * @return 
 */
int grasp_output_tile_information_nnoises(const grasp_results_t *output);


//////////

  
/**
 * Return position of the segment, dimension T
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return position of the segment, dimension T
 */
int grasp_output_tile_pixel_information_segment_time(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return position of the segment, dimension X
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return position of the segment, dimension X
 */
int grasp_output_tile_pixel_information_segment_col(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return position of the segment, dimension Y
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return position of the segment, dimension Y
 */
int grasp_output_tile_pixel_information_segment_row(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Internal T index of pixel inside the segment starting in 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Internal T index of pixel inside the segment starting in 1
 */
int grasp_output_tile_pixel_information_it(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Internal X index of pixel inside the segment starting in 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Internal X index of pixel inside the segment starting in 1
 */
int grasp_output_tile_pixel_information_ix(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Internal Y index of pixel inside the segment starting in 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Internal Y index of pixel inside the segment starting in 1
 */
int grasp_output_tile_pixel_information_iy(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Index of X index in output grid.
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Index of X index in output grid.
 */
int grasp_output_tile_pixel_information_out_x(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Index of Y index in output grid.
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Index of Y index in output grid.
 */
int grasp_output_tile_pixel_information_out_y(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Index of T index in output grid.
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Index of T index in output grid.
 */
int grasp_output_tile_pixel_information_out_t(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return latitude of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return latitude of the pixel
 */
float grasp_output_tile_pixel_information_latitude(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return longitude of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return longitude of the pixel
 */
float grasp_output_tile_pixel_information_longitude(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return index in input grid for cols
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return index in input grid for cols
 */
int grasp_output_tile_pixel_information_grid_col(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return index in input grid for rows
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return index in input grid for rows
 */
int grasp_output_tile_pixel_information_grid_row(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return index in input grid for time
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return index in input grid for time
 */
int64_t grasp_output_tile_pixel_information_time(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return number of wavelengths for the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return number of wavelengths for the pixel
 */
int grasp_output_tile_pixel_information_nwl(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return cloud flag of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return cloud flag of the pixel
 */
int grasp_output_tile_pixel_information_cloud_flag(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return land percent information of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return land percent information of the pixel
 */
float grasp_output_tile_pixel_information_land_percent(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return file index (origin) of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return file index (origin) of the pixel
 */
int grasp_output_tile_pixel_information_file_index(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return meters above sea-level of the pixel
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return meters above sea-leve of the pixel
 */
float grasp_output_tile_pixel_information_masl(const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return Solar zenit angle of current wavelength in degrees
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Solar zenit angle of current wavelength in degrees
 */
float grasp_output_tile_pixel_information_sza(const grasp_results_t *output,int it, int ix, int iy, int iwl);

/**
 * Return the number of valid pixels in this pixel's segment
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param nrCol number of columns of segments in tile
 * @param nrRow number if rows of segments in tile
 * @return the number of valid pixels in this pixel's segment
 */
int grasp_output_tile_pixel_segment_npixels(const grasp_results_t *output, int it, int ix, int iy, int nrCol, int nrRow);

//////////////

/**
 * Return number of iterations in multipixel scenario
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return number of iterations in multipixel scenario
 */
int grasp_output_tile_retrieval_res_niter (const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return total residual for multi-pixel retrieval
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return total residual for multi-pixel retrieval
 */
float grasp_output_tile_retrieval_res_rest (const grasp_results_t *output,int it, int ix, int iy);

/**
 * Return detailed absolute measurement residuals for tile
 * @param output Output of retrieval from a tile 
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param inoise Number of noise
 * @return detailed absolute measurement residuals for tile
 */
float grasp_output_tile_retrieval_res_resa (const grasp_results_t *output,int it, int ix, int iy, int inoise);

/**
 * Return detailed relative measurement residuals for tile
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param inoise Number of noise
 * @return detailed relative and relative measurement residuals for tile
 */
float grasp_output_tile_retrieval_res_resr (const grasp_results_t *output,int it, int ix, int iy, int inoise);


/**
 * Return retrieved aerosol and surface reflectance parameters
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param ipar
 * @return retrieved aerosol and surface reflectance parameters
 */
float grasp_output_tile_retrieval_par_parameters (const grasp_results_t *output, int it, int ix, int iy, int ipar);

/**
 * Return original pixel data
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Original pixel data. It is like sdata pixel with input information
 */
const pixel_t *grasp_output_tile_retrieval_fit_pixel_original (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return fit pixel data
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return Fit pixel. It is like sdata pixel but filled with fitting information
 */
const pixel_t *grasp_output_tile_retrieval_fit_pixel_fit (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return angstrom exponent for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return angstrom exponent for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_aexp (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return total extinction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total extinction for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total single scattering albedo for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total single scattering albedo for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total absorption extinction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total absorption extinction for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_aextt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return extinction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return extinction for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return single scattering albedo for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return single scattering albedo for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return absorption extinction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return absorption extinction for aerosol optical properties
 */
float grasp_output_tile_aerosol_opt_aext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return real part refractive index 
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return real part refractive index 
 */
float grasp_output_tile_aerosol_rind_mreal (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return real part refractive index 
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return real part refractive index 
 */
float grasp_output_tile_aerosol_rind_mimag (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

//////

/**
 * Return p11 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p11 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p12 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p12 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p22 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p22 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p33 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p33 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p34 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p34 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p44 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of aerosol component
 * @return p44 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_ph44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return total p11 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p11 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p12 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p12 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p22 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p22 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p33 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p33 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p34 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p34 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p44 for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p44 for aerosol optical properties
 */
float grasp_output_tile_aerosol_phmx_pht44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/////

/**
 * Return lidar ratio for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return lidar ratio for aerosol optical properties
 */
float grasp_output_tile_aerosol_lidar_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return lidar depolarization profile for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return lidar depolarization profile for aerosol optical properties
 */
float grasp_output_tile_aerosol_lidar_ldpr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return total lidar ratio for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total lidar ratio for aerosol optical properties
 */
float grasp_output_tile_aerosol_lidar_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total lidar depolarization for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total lidar depolarization for aerosol optical properties
 */
float grasp_output_tile_aerosol_lidar_ldprt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return concentration for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_cv (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_std (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return XXX for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_rm (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return refractive index for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_reff (const grasp_results_t *output, int it, int ix, int iy, int i);


/**
 * Return extinction for two simulated modes for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param i mode. 0 for fine and 1 for coarse
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int i);

/**
 * Return concentration for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_cv with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_cv_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_std with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_std_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return XXX for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_rm with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_rm_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return refractive index for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_reff with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_reff_fine_mode (const grasp_results_t *output, int it, int ix, int iy);


/**
 * Return extinction for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_opt_ext with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_opt_ext_fine_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl);


/**
 * Return concentration for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_cv with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return concentration for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_cv_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return standard deviation for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_std with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return standard deviation for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_std_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return XXX for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_rm with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return XXX for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_rm_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return refractive index for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_mph_reff with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return refractive index for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_mph_reff_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);


/**
 * Return extinction for two simulated modes for aerosol optical properties. Alias to grasp_output_tile_aerosol_sd2m_opt_ext with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for aerosol optical properties
 */
float grasp_output_tile_aerosol_sd2m_opt_ext_coarse_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl);


/**
 * Return relative humidity for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return relative humidity for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_rh (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return water fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return water fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_fwtr (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return soluble fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return soluble fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_fslbl (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return insoluble fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return  insoluble fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_finslbl (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return soot fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return soot fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_fsoot (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return iron fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return iron fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_firon (const grasp_results_t *output, int it, int ix, int iy, int isd);


/** add by lei on 15/11/2016
 * Return brc fraction for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of aerosol component
 * @return iron fraction for aerosol optical properties
 */
float grasp_output_tile_aerosol_chem_fbrc (const grasp_results_t *output, int it, int ix, int iy, int isd);



/**
 * Return particular matter for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i
 * @return particular matter for aerosol optical properties
 */

float grasp_output_tile_aerosol_pm_pm (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return aerosol type for aerosol optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return aerosol type: 0 – Complex mixture; 1 – Background Aerosol; 2 – Water/Maritime; 3 — Urban Polluted; 4 – Mixed aerosol;  5 – Urban Clean;  6 – Smoke Smoldering;  7 – Smoke flaming;  8 – Mineral dust
 */
int grasp_output_tile_aerosol_types_index (const grasp_results_t *output, int it, int ix, int iy);

/////////

//////


/**
 * Return angstrom exponent for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return angstrom exponent for clouds optical properties
 */
float grasp_output_tile_clouds_opt_aexp (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return total extinction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total extinction for clouds optical properties
 */
float grasp_output_tile_clouds_opt_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total single scattering albedo for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total single scattering albedo for clouds optical properties
 */
float grasp_output_tile_clouds_opt_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total absorption extinction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total absorption extinction for clouds optical properties
 */
float grasp_output_tile_clouds_opt_aextt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return extinction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return extinction for clouds optical properties
 */
float grasp_output_tile_clouds_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return single scattering albedo for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return single scattering albedo for clouds optical properties
 */
float grasp_output_tile_clouds_opt_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return absorption extinction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return absorption extinction for clouds optical properties
 */
float grasp_output_tile_clouds_opt_aext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return real part refractive index 
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return real part refractive index 
 */
float grasp_output_tile_clouds_rind_mreal (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return real part refractive index 
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return real part refractive index 
 */
float grasp_output_tile_clouds_rind_mimag (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

//////

/**
 * Return p11 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p11 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p12 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p12 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p22 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p22 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p33 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p33 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p34 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p34 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return p44 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @param isd Index of clouds component
 * @return p44 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_ph44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar, int isd);

/**
 * Return total p11 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p11 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht11 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p12 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p12 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht12 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p22 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p22 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht22 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p33 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p33 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht33 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p34 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p34 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht34 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/**
 * Return total p44 for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param impar
 * @return total p44 for clouds optical properties
 */
float grasp_output_tile_clouds_phmx_pht44 (const grasp_results_t *output, int it, int ix, int iy, int iwl, int impar);

/////

/**
 * Return lidar ratio for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return lidar ratio for clouds optical properties
 */
float grasp_output_tile_clouds_lidar_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return lidar depolarization profile for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return lidar depolarization profile for clouds optical properties
 */
float grasp_output_tile_clouds_lidar_ldpr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return total lidar ratio for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total lidar ratio for clouds optical properties
 */
float grasp_output_tile_clouds_lidar_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return total lidar depolarization for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return total lidar depolarization for clouds optical properties
 */
float grasp_output_tile_clouds_lidar_ldprt (const grasp_results_t *output, int it, int ix, int iy, int iwl);


/**
 * Return concentration for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_cv (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return standard deviation for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_std (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return XXX for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_rm (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return refractive index for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i mode. 0 for fine and 1 for coarse
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_reff (const grasp_results_t *output, int it, int ix, int iy, int i);


/**
 * Return extinction for two simulated modes for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param i mode. 0 for fine and 1 for coarse
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_opt_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int i);


/**
 * Return concentration for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_cv with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_cv_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return standard deviation for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_std with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_std_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return XXX for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_rm with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_rm_fine_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return refractive index for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_reff with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_reff_fine_mode (const grasp_results_t *output, int it, int ix, int iy);


/**
 * Return extinction for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_opt_ext with i = 0
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_opt_ext_fine_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl);


/**
 * Return concentration for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_cv with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return concentration for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_cv_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return standard deviation for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_std with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return standard deviation for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_std_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return XXX for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_rm with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return XXX for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_rm_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return refractive index for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_mph_reff with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return refractive index for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_mph_reff_coarse_mode (const grasp_results_t *output, int it, int ix, int iy);


/**
 * Return extinction for two simulated modes for clouds optical properties. Alias to grasp_output_tile_clouds_sd2m_opt_ext with i = 1
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return extinction for two simulated modes for clouds optical properties
 */
float grasp_output_tile_clouds_sd2m_opt_ext_coarse_mode (const grasp_results_t *output, int it, int ix, int iy, int iwl);


/**
 * Return relative humidity for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return relative humidity for clouds optical properties
 */
float grasp_output_tile_clouds_chem_rh (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return water fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return water fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_fwtr (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return soluble fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return soluble fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_fslbl (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return insoluble fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return  insoluble fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_finslbl (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return soot fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return soot fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_fsoot (const grasp_results_t *output, int it, int ix, int iy, int isd);

/**
 * Return iron fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return iron fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_firon (const grasp_results_t *output, int it, int ix, int iy, int isd);


/**
 * Return brc fraction for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param isd Index of clouds component
 * @return iron fraction for clouds optical properties
 */
float grasp_output_tile_clouds_chem_fbrc (const grasp_results_t *output, int it, int ix, int iy, int isd);



/**
 * Return particular matter for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param i
 * @return particular matter for clouds optical properties
 */
float grasp_output_tile_clouds_pm_pm (const grasp_results_t *output, int it, int ix, int iy, int i);

/**
 * Return clouds type for clouds optical properties
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return clouds type: 0 – Complex mixture; 1 – Background Aerosol; 2 – Water/Maritime; 3 — Urban Polluted; 4 – Mixed clouds;  5 – Urban Clean;  6 – Smoke Smoldering;  7 – Smoke flaming;  8 – Mineral dust
 */
int grasp_output_tile_clouds_types_index (const grasp_results_t *output, int it, int ix, int iy);


////////

/**
 * Return surface ndvi
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return surface ndvi
 */
float grasp_output_tile_surface_ndvi (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return surface albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return surface albedo
 */
float grasp_output_tile_surface_salbedo (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/////////////////

/**
 * Return error estimation for parameter
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param ipar
 * @return error estimation for parameter
 */
float grasp_output_tile_errest_par_errp (const grasp_results_t *output, int it, int ix, int iy, int ipar);

/**
 * Return bias of error estimation of parameters
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param ipar
 * @return bias of error estimation of parameters
 */
float grasp_output_tile_errest_par_biasp (const grasp_results_t *output, int it, int ix, int iy, int ipar);

/**
 * Return error estimation of extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of extinction
 */
float grasp_output_tile_errest_aerosol_opt_err_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of extinction
 */
float grasp_output_tile_errest_aerosol_opt_bias_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total extinction
 */
float grasp_output_tile_errest_aerosol_opt_err_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total extinction
 */
float grasp_output_tile_errest_aerosol_opt_bias_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return error estimation of single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of single scattering albedo
 */
float grasp_output_tile_errest_aerosol_opt_err_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of single scattering albedo
 */
float grasp_output_tile_errest_aerosol_opt_bias_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total single scattering albedo
 */
float grasp_output_tile_errest_aerosol_opt_err_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total single scattering albedo 
 */
float grasp_output_tile_errest_aerosol_opt_bias_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return error estimation of lidar ration
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return error estimation of lidar ration
 */
float grasp_output_tile_errest_aerosol_lidar_err_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of aerosol component
 * @return bias of error estimation of lidar ratio
 */
float grasp_output_tile_errest_aerosol_lidar_bias_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total lidar ratio
 */
float grasp_output_tile_errest_aerosol_lidar_err_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total lidar ratio
 */
float grasp_output_tile_errest_aerosol_lidar_bias_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return error estimation of extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of extinction
 */
float grasp_output_tile_errest_clouds_opt_err_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of extinction
 */
float grasp_output_tile_errest_clouds_opt_bias_ext (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total extinction
 */
float grasp_output_tile_errest_clouds_opt_err_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total extinction
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total extinction
 */
float grasp_output_tile_errest_clouds_opt_bias_extt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return error estimation of single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of single scattering albedo
 */
float grasp_output_tile_errest_clouds_opt_err_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of single scattering albedo
 */
float grasp_output_tile_errest_clouds_opt_bias_ssa (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total single scattering albedo
 */
float grasp_output_tile_errest_clouds_opt_err_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total single scattering albedo
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total single scattering albedo 
 */
float grasp_output_tile_errest_clouds_opt_bias_ssat (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return error estimation of lidar ration
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return error estimation of lidar ration
 */
float grasp_output_tile_errest_clouds_lidar_err_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return bias of error estimation of lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @param isd Index of clouds component
 * @return bias of error estimation of lidar ratio
 */
float grasp_output_tile_errest_clouds_lidar_bias_lr (const grasp_results_t *output, int it, int ix, int iy, int iwl, int isd);

/**
 * Return error estimation of total lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return error estimation of total lidar ratio
 */
float grasp_output_tile_errest_clouds_lidar_err_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

/**
 * Return bias of error estimation of total lidar ratio
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iwl Index of the wavelength
 * @return bias of error estimation of total lidar ratio
 */
float grasp_output_tile_errest_clouds_lidar_bias_lrt (const grasp_results_t *output, int it, int ix, int iy, int iwl);

///////////////

/**
 * Return number of heights of forcing fluxes
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return number of heights of forcing fluxes
 */
int grasp_output_tile_forcing_bbflux_nhlv (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return broad band up-ward flux without aerosol at each height
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return broad band up-ward flux without aerosol at each height
 */
float grasp_output_tile_forcing_bbflux_bbufx0 (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return broad band down-ward flux without aerosol at each height
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return broad band down-ward flux without aerosol at each height
 */
float grasp_output_tile_forcing_bbflux_bbdfx0 (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return broad band up-ward flux with aerosol at each height
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return broad band up-ward flux with aerosol at each height
 */
float grasp_output_tile_forcing_bbflux_bbufxa (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return broad band down-ward flux with aerosol at each height
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return broad band down-ward flux with aerosol at each height
 */
float grasp_output_tile_forcing_bbflux_bbdfxa (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return heights of forcing fluxes
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return heights of forcing fluxes
 */
float grasp_output_tile_forcing_bbflux_hlv (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return number of heights of forcing calculations
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @return number of heights of forcing calculations
 */
int grasp_output_tile_forcing_forcing_nhlv (const grasp_results_t *output, int it, int ix, int iy);

/**
 * Return net forcing
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return net forcing
 */
float grasp_output_tile_forcing_forcing_netforc (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return forcing efficiency
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return forcing efficiency
 */
float grasp_output_tile_forcing_forcing_forceff (const grasp_results_t *output, int it, int ix, int iy, int iknt);

/**
 * Return heights of forcing calculations
 * @param output Output of retrieval from a tile
 * @param it Time index of the pixel in the tile
 * @param ix X index of the pixel in the tile
 * @param iy Y index of the pixel in the tile
 * @param iknt
 * @return heights of forcing calculations
 */
float grasp_output_tile_forcing_forcing_hlv (const grasp_results_t *output, int it, int ix, int iy, int iknt);


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_OUTPUT_TILE_RESULT_H */

