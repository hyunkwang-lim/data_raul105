/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_settings_t.h
 * Author: fuertes
 *
 * Created on 3 de octubre de 2014, 14:46
 */

#ifndef GRASP_SETTINGS_T_H
#define	GRASP_SETTINGS_T_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef WARN_DRY
#warning "__RETRIEVAL_SETTINGS_DEFINITION__ binded"
#endif    
    
typedef struct dls_{
        int IWL                                    ;
        int key                                    ;
        int keyEL                                  ;
        char distname_O[_GBL_FILE_PATH_LEN]        ;   
        char distname_N[_GBL_FILE_PATH_LEN]        ;       
        char internal_file_path[_GBL_FILE_PATH_LEN];   
        char external_file_path[_GBL_FILE_PATH_LEN];       
} dls;


typedef struct osh_{
    int IMSC;
    int NG;
    int NN;
    int NF;
} osh;

typedef struct single_pix_contraints_apriori_{
    int   IO  [_KIDIM1][_KIDIM2][_KIDIM3];
    float GSM [_KIDIM1][_KIDIM2][_KIDIM3];
} single_pix_contraints_apriori;

typedef struct single_pix_contraints_smoothness_{
    int   IO  [_KIDIM1][_KIDIM2];
    float GSM [_KIDIM1][_KIDIM2];
} single_pix_contraints_smoothness;

typedef struct multi_pix_constraints_{
    int   IOT  [_KIDIM1][_KIDIM2];
    int   IOX  [_KIDIM1][_KIDIM2];
    int   IOY  [_KIDIM1][_KIDIM2];
    float GSMT [_KIDIM1][_KIDIM2];
    float GSMX [_KIDIM1][_KIDIM2];
    float GSMY [_KIDIM1][_KIDIM2];
}multi_pix_constraints;
           
typedef struct NOISE_param_{
    int INOISE;
    float SGMS [_KKNOISE];
    int INN    [_KKNOISE];
    float DNN  [_KKNOISE];
    int NMT    [_KKNOISE];
    int MT     [_KKNOISE][_KIP];
    int NWLP   [_KKNOISE][_KIP];
    int IWLP   [_KKNOISE][_KIP][_KWM];
} NOISE_param;

typedef struct{
    int INVSING;
    bool TXY_group;
    float TTEST;
    float XTEST;
    float YTEST;
}inter_pixel_fit;

//  Edge sizes
typedef struct edges_size_{
         int nx;
         int ny;        
         int nt;         
}edges_size;

typedef struct retr_input_{
        int   KNSING      ;
        int   KNSINGF     ;           
        int   KL          ;      
        bool  ISTOP       ;          
        bool  IPRI_additional_info; 
        bool  IPRI_iter_result;
        bool  IPRI_verbose;
        int   NSD         ;       
        int   NLYRS       ;         
        int   ipplane     ;            
        int   iPOBS       ; 
        int   iIOBS       ;
        int   isurf_land[2];
        int   isurf_water ;
        int   Aexp_iwl[2] ;
        int   ndvi_iwl[2] ;
        float SHIFT       ;
        int   NW          ;
        float WAVE [_KW]  ;
        int   IBIN        ;
        int   IMQ         ;
        int   IPSTOP      ;
        bool  INPUT       ; 
        bool  ITRONC      ;
        int   MAXP        ;
        float EPSP        ;
        float EPSQ        ;
        float DL          ;
        int  mol_prof_type;
        int  aer_prof_type; // type of the aerosol vertical profile used
        float  PM_diam[2] ; // diameters at wich PM is calculated
        int   nPM_diam    ; // number of PM diamteres
        bool  use_models  ; // use modeled phase matrices instead of LB Kernels
  
        dls DLSF          ;
        par_number_NDIM NDIM; 
        osh OSHF          ;
        osh OSHD          ;
        single_pix_contraints_apriori SPCA;
        single_pix_contraints_smoothness SPCS;
        multi_pix_constraints MPCS;        
        NOISE_param NOISE ;   
        inter_pixel_fit IPFP;
        
        float APSING[_KPARS];
        float APSMIN[_KPARS];
        float APSMAX[_KPARS];

        float RMIN[_KSD];
        float RMAX[_KSD];
	      float RATIO1   [_KSD][_KIDIM3];
	      float RADIUS1  [_KSD][_KIDIM3];
        int   IWW_SINGL[_KPARS];
        int   NBIN       [_KSD];
        
        edges_size  edges;
        
	      char plotting_output_file [_GBL_FILE_PATH_LEN];
	      char main_output_file     [_GBL_FILE_PATH_LEN];
	      char sdata_sim_file       [_GBL_FILE_PATH_LEN]; 
                
        float eps_err;
                        
        output_segment_products products;            
} retr_input;

typedef struct temporal_data_{
    int NOISE_ISGMS;
    int NOISE_IDNN;   
    int IAPSING;
    int IAPSMIN;
    int IAPSMAX;
    int SPCA_IIO[_KIDIM1][_KIDIM2];
    int SPCA_IGSM[_KIDIM1][_KIDIM2];
    int SPCS_IIO[_KIDIM1];
    int SPCS_IGSM[_KIDIM1]; 
    
    int MPCS_IIOT[_KIDIM1];
    int MPCS_IGSMT[_KIDIM1];
    int MPCS_IIOX[_KIDIM1];
    int MPCS_IGSMX[_KIDIM1];
    int MPCS_IIOY[_KIDIM1];
    int MPCS_IGSMY[_KIDIM1];    
    
    int IIWW_SINGL;    
    
    int NDIM_npar_type;
    int NDIM_npar_retr;
    
    float TAPSING[_KIDIM1][_KIDIM2][_KPARS];
    int NTAPSING[_KIDIM1][_KIDIM2];
    
    float TAPSMIN[_KIDIM1][_KIDIM2][_KPARS];
    int NTAPSMIN[_KIDIM1][_KIDIM2];
    
    float TAPSMAX[_KIDIM1][_KIDIM2][_KPARS];
    int NTAPSMAX[_KIDIM1][_KIDIM2];
    
    int TIWW_SINGL[_KIDIM1][_KIDIM2][_KPARS];
    int NTIWW_SINGL[_KIDIM1][_KIDIM2];
    
    int NRMIN;
    int NRMAX;
    int NRATIO1[_KIDIM3];
    
    int nisurf_land;
    
    int nAexp_iwl;
    
    int nndvi_iwl;

} temporal_data;

/*
 * Settings of grasp settings module
 */
typedef struct grasp_settings_settings_{
    // True if debug information have to be dumped
    char debug[_GBL_FILE_PATH_LEN];
    // True if help information will be printed
    char help[_GBL_FILE_PATH_LEN];
    // if there was warning during the reading process of settings file and this
    // variable is false we will continue the workflow without stop (but informing the user)
    bool strict;   
    // Output stream patter to dump setting loaded in short format
    char short_dump[_GBL_FILE_PATH_LEN];
    // Output stream patter to dump setting loaded in long format
    char long_dump[_GBL_FILE_PATH_LEN];
    // If it is true version code version will be printed in 
    bool version;
}grasp_settings_settings;

typedef struct grasp_global_t_{
    char resources_path[_GBL_FILE_PATH_LEN];     
}grasp_global_t;

typedef struct grasp_settings_{
    retr_input              retrieval;
    input_settings_t        input;
    controller_settings_t   controller;
    output_settings_t       output;
    grasp_settings_settings settings;
    grasp_global_t          global;
    temporal_data           tmp;
}grasp_settings;


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_SETTINGS_T_H */

