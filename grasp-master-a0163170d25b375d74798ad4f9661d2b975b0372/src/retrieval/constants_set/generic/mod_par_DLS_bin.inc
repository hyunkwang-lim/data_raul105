/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */

/*      MODULE mod_par_DLS_bin
!c ********************************************************** c
!c **  Parameters for _bin version of DLS code             ** c
!c ********************************************************** c
!c **  KWLpar - number of wavelengths in fixed kernels     ** c
!c **  KCpar  - number of concentrations for recalculated  ** c
!c **           bins                                       ** c
!c ********************************************************** c
*/


/***  DLS_bin part of code ***/

/* parameter for number of wavelengths */
#define _KWLpar 12
/* parameter for number of lognormal bins */
#define _KCpar 9  

/*** Inversion code ***/

/* max number of lognormal bins or number of grid radii */
/* #define _NRC 22 */
#define _NRC 9

/* number of radii for printing SD  */
#define _NRR 45
	  

