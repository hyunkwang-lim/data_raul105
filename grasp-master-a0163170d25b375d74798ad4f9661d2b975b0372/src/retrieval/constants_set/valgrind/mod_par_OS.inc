/*
 *  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.
 *  
 *  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
 */


/* NGO: number of terms in phase matrix expansion */
#define _NG0 41

/* NNO: number of terms in directional integration */
#define _NN0 20
 
/*
!PLTL NG0T = NG0*4 if ground based observations are present
!PLTL NG0T = NG0   for satellite observations
      INTEGER, PARAMETER ::   NG0T = _NG0T  ! ground based observations are present
!      INTEGER, PARAMETER ::   NG0T = NG0   ! for satellite observations  
 */
#define _NG0T (_NG0*4)
 
/*
!PL NF0 is maximum number of Fourier expansion.
!PL In general NF0<=2*NG0-2 but for aerosol we don't need such accuracy.
!PL Therefore, we take:		
 */
#define _NF0 _NN0 

/* KSD: number of aerosol components */
#define _KSD 2

/* KNT: number of atmospheric layers */
#define _KNT 51 

/* NMM: number of atmospheric component in OS */
#define _NMM (KSD+1)

/* NBVM: number of observation angles for each measure */
#define _NBVM 300 

/* NMG: number of gas component with profile */
#define _NMG 0 
      
/* HMAX_atm (m) : atmosphere max height */
#define _HMAX_atm 40000.0 

/* KVERT_WD: number of profile grid heights */
#define _KVERT_WD 100 
