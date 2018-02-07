! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "mod_par_DLS_bin.inc"
    MODULE mod_par_DLS_bin
!c ********************************************************** c
!c **  Parameters for _bin version of DLS code             ** c
!c ********************************************************** c
!c **  KWLpar - number of wavelengths in fixed kernels     ** c
!c **  KCpar  - number of concentrations for recalculated  ** c
!c **           bins                                       ** c
!c ********************************************************** c

!      IMPLICIT NONE

! *** DLS_bin part of code

! parameter for number of wavelengths 
      INTEGER, PARAMETER ::  KWLpar = _KWLpar
! parameter for number of lognormal bins 
      INTEGER, PARAMETER ::  KCpar  = _KCpar

! *** Inversion code

! max number of lognormal bins or number of grid radii 
!      INTEGER, PARAMETER ::	 NRC    = 22 
      INTEGER, PARAMETER ::	 NRC    = _NRC
! number of radii for printing SD 
      INTEGER, PARAMETER ::	 NRR    = _NRR
	 
      contains
      
        subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
   	end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
        
	END MODULE mod_par_DLS_bin
