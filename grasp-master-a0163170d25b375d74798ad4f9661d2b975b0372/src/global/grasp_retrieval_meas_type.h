/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/* 
 * File:   grasp_retrieval_meas_type.h
 * Author: fuertes
 *
 * Created on 28 de octubre de 2013, 10:34
 */

#ifndef GRASP_RETRIEVAL_MEAS_TYPE_H
#define	GRASP_RETRIEVAL_MEAS_TYPE_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef WARN_DRY
#warning "__MEAS_TYPE__ duplicated"
#endif    
    
enum {
  MEAS_TYPE_UNKNOWN = 0,
  // tau beginning
  MEAS_TYPE_TOD = 11,
  MEAS_TYPE_AOD = 12,
  // tau end
  // phase matrix beginning
  MEAS_TYPE_P11 = 21,
  MEAS_TYPE_P12 = 22,
  MEAS_TYPE_P22 = 23,
  MEAS_TYPE_P33 = 24,
  MEAS_TYPE_P34 = 25,
  MEAS_TYPE_P44 = 26,
  // phase matrix end
  // lidar beginning
  MEAS_TYPE_LS  = 31,
  MEAS_TYPE_DP  = 32,
  MEAS_TYPE_RL  = 33,
  // lidar end
  // SvR beginning
  MEAS_TYPE_I   = 41,
  MEAS_TYPE_Q   = 42,
  MEAS_TYPE_U   = 43,
  MEAS_TYPE_P   = 44
  // SvR end
};
#define GRASP_NMEAS_TYPES 14


#ifdef	__cplusplus
}
#endif

#endif	/* GRASP_RETRIEVAL_MEAS_TYPE_H */

