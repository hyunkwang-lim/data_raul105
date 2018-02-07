! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "mod_par_OS.inc"  
    MODULE mod_par_OS
!c ********************************************************** c
!c **   MAXIMUM DIMENSIONS as follows:                     ** c
!c ********************************************************** c
!C **                                                      ** c
!C ** : related to measurements specifications :           ** c        
!c **                                                      ** c
!c **  NG0    - number of terms in phase matrix expansion  ** c
!c **  NN0    - number of terms in directional integration ** c
!c **  KSD    - number of aerosol components               ** c
!c **  KNT    - number of atmospheric layers               ** c
!c **  NMM    - number of atmospheric component in OS      ** c
!c **  NBVM   - number of observation angles for each      ** c
!c **                                          wavelength  ** c
!c **  NF     - number of                                  ** c
!c **  IIX    - parameter for "gauss"                      ** c
!c **                                                      ** c
!c ********************************************************** c

      IMPLICIT NONE

      INTEGER, PARAMETER ::   NG0  = _NG0
      INTEGER, PARAMETER ::   NN0  = _NN0

!PLTL NG0T = NG0*4 if ground based observations are present
!PLTL NG0T = NG0   for satellite observations
      INTEGER, PARAMETER ::   NG0T = _NG0T  ! ground based observations are present
!      INTEGER, PARAMETER ::   NG0T = NG0   ! for satellite observations  
      
!PL NF0 is maximum number of Fourier expansion.
!PL In general NF0<=2*NG0-2 but for aerosol we don't need such accuracy.
!PL Therefore, we take:		
      INTEGER, PARAMETER ::   NF0  = _NF0 
      INTEGER, PARAMETER ::   KSD  = _KSD
      INTEGER, PARAMETER ::   KNT  = _KNT
      INTEGER, PARAMETER ::   NMM  = _NMM
      INTEGER, PARAMETER ::   NBVM = _NBVM
      INTEGER, PARAMETER ::   NMG  = _NMG ! number of gas component with profile
      !INTEGER, PARAMETER ::   IIX  = NG0 !  =600 
      REAL,    PARAMETER ::   HMAX_atm = _HMAX_atm ! atmosphere max height 
      INTEGER, PARAMETER ::   KVERT_WD = _KVERT_WD ! number of profile grid heights
      
      contains
      
        subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
   	end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
        
      END MODULE mod_par_OS


