! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "mod_par_DLS.inc"
    MODULE mod_par_DLS
! ********************************************************** c
! **  Parameters for matrix_...f                          ** c
! ********************************************************** c
! **  KN1par - number of grid radii in original kernels   ** c
! **           and fixed kernels                          ** c
! **  KM1par - number of scattering angles in org.kernels ** c
! **  KR1par - number of grid aspect ratios for original  ** c
! **                                            kernels   ** c
! **  KREpar - number of refractive index real parts      ** c
! **  KIMpar - number of refractive index imaginary parts ** c
! **                                                      ** c
! **  KNpar  - number of grid radii for optical charact.  ** c
! **  KRpar  - number of grid aspect ratios for axial     ** c
! **           ratio distribution                         ** c
! **  KMpar  - number of scattering angles in fixed       ** c
! **           kernels and for opt.characteristics        ** c
! **  KMD    - number of modes in size distribution       ** c
! **           (up to 2)                                  ** c
! **  rootdir - kernel directory name                     ** c
! ********************************************************** c
                      
!        IMPLICIT NONE

! parameters for Original or Fixed Kernels

      INTEGER, PARAMETER ::  KN1par = _KN1par
      INTEGER, PARAMETER ::  KR1par = _KR1par  
      INTEGER, PARAMETER ::  KM1par = _KM1par 
      INTEGER, PARAMETER ::  KREpar = _KREpar
      INTEGER, PARAMETER ::  KIMpar = _KIMpar 

! parameters for Optical Characteristics

      INTEGER, PARAMETER ::  KNpar  = _KNpar
      INTEGER, PARAMETER ::  KRpar  = _KRpar
      INTEGER, PARAMETER ::  KMpar  = _KMpar 
! parameters for Size Distribution

      INTEGER, PARAMETER ::  KMD  = _KMD 
      CHARACTER(len=60), PARAMETER :: rootdir=_rootdir
      
      contains
      
        subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
   	end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
        
      END MODULE mod_par_DLS
