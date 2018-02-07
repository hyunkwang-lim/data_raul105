! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "mod_par_inv.inc"
    module mod_par_inv
!c ********************************************************** c
!c **   MAXIMUM DIMENSIONS as follows:                     ** c
!c ********************************************************** c
!C **                                                      ** c
!C ** : related to measurements specifications :           ** c      
!c **                                                      ** c
!c **  KW     - number of wavelengths                      ** c
!c **  KSHAPE - number of aerosol components  with         ** c
!c **           different shape distributions              ** c
!c **  KBF    - number of parameters in BRF model          ** c
!c **  KANG   - number of scattering angles used in base   ** c
!c **           phase matrices                             ** c
!C **                                                      ** c
!C ** : related to INVERSION organization :                ** c      
!C **                                                      ** c
!c **  KITIME - number of pixels corresponding to          ** c
!c **           different times                            ** c
!c **  KIX    - number of pixels corresponding with        ** c
!c **           different X - coordinates                  ** c
!c **  KIY    - number of pixels corresponding with        ** c
!c **           different Y - coordinates                  ** c
!c **  KIDIM1 - -number of retrieved characteristics       ** c
!c **           (SD,REAL,IMAG,SHAPE, height )              ** c
!c **  KIDIM2 - number of sub-component in KIDIM1          ** c
!c **  KIDIM3 - number of sub/sub-component in KIDIM1      ** c
!c **  KIDIM4 - number of sub/sub/sub-component in KIDIM1  ** c
!c **  KPARS  - number of retrieved parameters             ** c
!c **           for each pixel                             ** c
!c **  KKNOISE- number of measured optical characteristics ** c
!c **             with different noise distribution        ** c
!c ********************************************************** c
      use mod_par_OS, only : NBVM   ! NBVM - number of valid directions for  
                                   ! angle dependend measurements 

      implicit none

!-----------------------------------------------------------------------------------------------
! ask Oleg kwm = 6 refr index in forw_model
      integer, parameter  ::  KWM     = _KWM ! max number of wave lengths for measurement
      integer, parameter  ::  KW      = _KW  ! max number of wave lengths for wave dependent 
                                             ! retrieved parameters (total number of differtent 
                                             ! wave lengths in inversion)
 
      integer, parameter  ::  KSHAPE  = _KSHAPE 
      integer, parameter  ::  KBF     = _KBF 
!-----------------------------------------------------------------------------------------------
      integer, parameter  ::  KITIME  = _KITIME ! 45x2x2; 45x3x3; 24x4x4; 18x5x5
      integer, parameter  ::  KIX     = _KIX
      integer, parameter  ::  KIY     = _KIY
      integer, parameter  ::  KIMAGE  = _KIMAGE ! number of pixels
!-----------------------------------------------------------------------------------------------       
      integer, parameter  ::  KIEDGE  = _KIEDGE ! max number of pixels for different 
                                                ! X,Y,T - coordinates of inverted segment edges
!-----------------------------------------------------------------------------------------------      
      integer, parameter  ::  KIDIM1  = _KIDIM1 
      integer, parameter  ::  KIDIM2  = _KIDIM2
      integer, parameter  ::  KIDIM3  = _KIDIM3   ! max of radius bins and number of retrieved grid 
                                                  ! heights for vertical profile (if present)

      integer, parameter  ::  KPARS   = _KPARS    ! number of parameters driving forward model
      integer, parameter  ::  KKNOISE = _KKNOISE  
!-----------------------------------------------------------------------------------------------
      INTEGER, PARAMETER  ::  KVERTM  = _KVERTM ! maximum of number of measured profile vert. heights and 
                                                ! retrieved profile vert. heights
                                                
      integer, parameter  ::  KNBVM   = _KNBVM  ! parameter for size of measurement 
                                                ! arrays (max of number of angles 
                                                ! or heights)
!-----------------------------------------------------------------------------------------------
      integer, parameter  ::  KIP     = _KIP    ! max number of meas types for single pixel

      integer, parameter  ::  KMESS   = _KMESS 
      integer, parameter  ::  KPAR    = _KPAR
      integer, parameter  ::  nnz_par = _nnz_par ! parameter for sparse matrix 
!-----------------------------------------------------------------------------------------------
      integer, parameter  ::  KMPSM   = _KMPSM  ! parameter for multi pixel smoothness
      integer, parameter  ::  KDIF    = _KDIF   ! parameter for differences 
                                                ! in smoothness KDIF=MAX of 
                                                ! (KPARS,KTIME,KIX,KIY) 
!-----------------------------------------------------------------------------------------------
      integer, parameter  ::  KSURF  = _KSURF  
      !integer, parameter  ::  KGAS   = _KGAS       
!-----------------------------------------------------------------------------------------------

      contains
      
        subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
   	end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
          
      end module mod_par_inv
