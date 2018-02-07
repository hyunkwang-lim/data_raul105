! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module mod_index_cloud

    use mod_time_utils, only : KIND_TIME
    use mod_par_inv, only : KITIME,KIX,KIY,KIMAGE !,KW,KPARS            
	  
    implicit none
! -----------------------------------------------------------------------
!  Indices for recurcive function CHECK

    type ind_clouds

        integer :: NX
        integer	:: NY
        integer :: NT
        integer,dimension(KIX,KIY,KITIME) :: INIMAGE	
        integer,dimension(KIX,KIY,KITIME) :: ICLOUD
        integer,dimension(KIMAGE) :: ITIMAGE
        integer,dimension(KIMAGE) :: IYIMAGE
        integer,dimension(KIMAGE) :: IXIMAGE
        real,dimension(KIX,KIY,KITIME) :: X	
        real,dimension(KIX,KIY,KITIME) :: Y	
        integer(KIND_TIME),dimension(KIX,KIY,KITIME) :: T	

    end type ind_clouds
!	---------------------------------------------------------------------
!  Multi pixel smoothness for each pixel

		!type mp_smoothness

			!integer					        :: ksm
			!integer,dimension(KIMAGE,KPARS)	::	imsm
			!real,   dimension(KIMAGE,KPARS)	::	smim	

		!end type mp_smoothness
!	---------------------------------------------------------------------
      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module
        
end module mod_index_cloud

