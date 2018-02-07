/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_products.h
 * Author: fuertes
 *
 * Created on May 30, 2014, 10:35 AM
 */

#ifndef GRASP_PRODUCTS_H
#define	GRASP_PRODUCTS_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef WARN_DRY
#warning "__RETRIEVAL_PRODUCTS_DEFINITION__ binded"
#endif      
    
//  Output products                 
      typedef struct output_segment_products_retrieval_ {
         bool res;
         bool par;
         bool fit;
      } output_segment_products_retrieval;

      typedef struct output_segment_products_particles_ {
         bool opt ; //ext ssa Aexp
         bool rind;
         bool chem;
         bool phmx;
         bool lidar; //lr ldpr
         bool sd2m_mph;
         bool sd2m_ext;
         bool pm; // particulate matter
         bool types;
      } output_segment_products_particles;

      typedef struct output_segment_products_surface_ {
         bool surf;
      } output_segment_products_surface;

      typedef struct output_segment_products_forcing_ {
         bool bbflux;
         bool forcing;
      } output_segment_products_forcing;

      typedef struct output_segment_products_errest_particles_ {
         bool opt;        // ext ssa
         bool lidar;      // lr     
      } output_segment_products_errest_particles;

      typedef struct output_segment_products_errest_ {
         bool par;
         output_segment_products_errest_particles aerosol;
         output_segment_products_errest_particles clouds; 
      } output_segment_products_errest;

      typedef struct output_segment_products_ {
         output_segment_products_retrieval retrieval;
         output_segment_products_particles aerosol;
         output_segment_products_particles clouds;
         output_segment_products_surface   surface;
         output_segment_products_forcing   forcing;
         output_segment_products_errest    errest;
      } output_segment_products;    
    

#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_PRODUCTS_H */

