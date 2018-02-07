! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
module mod_sdata_meas_type

      implicit none 	
!	------------------------------------------------------------------------------------------------------
#ifdef WARN_DRY
#warning "__MEAS_TYPE__ duplicated"
#endif   
        
      integer,parameter ::  meas_type_tau_beg = 10
      integer,parameter ::	meas_type_tod = 11   ! tod(wl) = aer+mol+gas - total optical depth
      integer,parameter ::	meas_type_aod = 12   ! aod(wl) - aerosol optical depth
      integer,parameter ::  meas_type_tau_end = 20

      integer,parameter ::  meas_type_phm_beg = 20
      integer,parameter ::	meas_type_p11 = 21   ! p11(angle,wl)  - phase matrix element P11
      integer,parameter ::	meas_type_p12 = 22   ! p12(angle,wl)  - phase matrix element P12
      integer,parameter ::	meas_type_p22 = 23   ! p22(angle,wl)  - phase matrix element P22
      integer,parameter ::	meas_type_p33 = 24   ! p33(angle,wl)  - phase matrix element P33
      integer,parameter ::	meas_type_p34 = 25   ! p34(angle,wl)  - phase matrix element P34
      integer,parameter ::	meas_type_p44 = 26   ! p44(angle,wl)  - phase matrix element P44
      integer,parameter ::  meas_type_phm_end = 30
      
      integer,parameter ::  meas_type_lid_beg = 30
      integer,parameter ::	meas_type_LS  = 31   ! LS(height,wl)  - Lidar Signal
      integer,parameter ::	meas_type_DP  = 32   ! DP(height,wl)  - DePolarization ratio
      integer,parameter ::	meas_type_RL  = 33   ! RL(height,wl)  - Raman Lidar signal
      integer,parameter ::  meas_type_lid_end = 40
      
      integer,parameter ::  meas_type_SvR_beg = 40
      integer,parameter ::	meas_type_I   = 41  ! I(angle,wl)    - Stokes parameter I
      integer,parameter ::	meas_type_Q   = 42  ! Q(angle,wl)    - Stokes parameter Q
      integer,parameter ::	meas_type_U   = 43  ! U(angle,wl)    - Stokes parameter U
      integer,parameter ::	meas_type_P   = 44  ! P(angle,wl)    - linear Polarization sqrt(Q*Q+U*U) or 
                                                !                  degree of linear Polarization sqrt(Q*Q+U*U)/I 
      integer,parameter ::  meas_type_SvR_end = 50

!	------------------------------------------------------------------------------------------------------

      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module                

end module mod_sdata_meas_type

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
