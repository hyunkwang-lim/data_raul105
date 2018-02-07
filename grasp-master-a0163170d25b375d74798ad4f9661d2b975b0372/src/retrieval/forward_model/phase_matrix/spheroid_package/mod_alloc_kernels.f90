! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

 
MODULE mod_alloc_kernels
	  
      use mod_par_DLS,     only: KN1par,KR1par,KM1par,KRpar,KMpar,KREpar,KIMpar
      use mod_par_DLS_bin, only: KWLpar,KCpar

      implicit none

! -----------------------------------------------------------------------
!  Kernels 1

      type kernels_triangle_bin
     
         real,dimension(:,:,:,:,:),allocatable   ::  U11, U12, &
                                                     U22, U33, &
                                                     U34, U44, & 
                                                     UEA         
      end type kernels_triangle_bin


! -----------------------------------------------------------------------
!  Kernels 2

      type kernels_lognormal_bin
     
         real,dimension(:,:,:,:,:,:),allocatable  :: UO11, UO12, &
                                                     UO22, UO33, &
                                                     UO34, UO44, &
                                                     UOEA
         end type kernels_lognormal_bin

contains

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

      subroutine alloc_kernel_arrays(key,keyEL,KERNELS1,KERNELS2)
	  
      !use mod_par_DLS,     only: KN1par,KR1par,KM1par,KRpar,KMpar,KREpar,KIMpar
      !use mod_par_DLS_bin, only: KWLpar,KCpar
      
      implicit none	  	  
         
      integer,                     intent(in)     ::  key,keyEL
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2

      integer                                     ::  ierr


!c*********************************************************
!c ** FOR single scattering aerosol OPTICAL properties !!!*

      IF(key .eq. 0 .or. key .eq. 3) THEN 
!c
!c *** ALLOCATE ARRAYS (alloc1.mod)
!c
         ALLOCATE(KERNELS1%UEA(KR1par,2,KN1par,KIMpar,KREpar),    stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%UEA array'
         KERNELS1%UEA=0.

         if(keyEL .gt. 0) then
         ALLOCATE(KERNELS1%U11(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U11 array'
         KERNELS1%U11=0. 
         endif
         
         if(keyEL .gt. 1) then
         ALLOCATE(KERNELS1%U12(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U12 array'
         KERNELS1%U12=0.
         endif
         
         if(keyEL .gt. 2) then
         ALLOCATE(KERNELS1%U22(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U22 array'
         KERNELS1%U22=0.
         endif
         
         if(keyEL .gt. 3) then
         ALLOCATE(KERNELS1%U33(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U33 array'
         KERNELS1%U33=0.
         endif
         
         if(keyEL .gt. 4) then
         ALLOCATE(KERNELS1%U34(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U34 array'
         KERNELS1%U34=0.
         endif
         
         if(keyEL .gt. 5) then
         ALLOCATE(KERNELS1%U44(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS1%U44 array'
         KERNELS1%U44=0.
         endif 
         
      ELSE !
       
!c *** ALLOCATE ARRAYS UFWL to be used in subroutine USMATRIX 
!c *** (in matrix_intrpl.wl1.f)
!c
         ALLOCATE(KERNELS2%UOEA(KRpar,2,KCpar,KIMpar,KREpar,KWLpar),    stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UOEA array'
         KERNELS2%UOEA=0.
         
         if(keyEL .gt. 0) then        
         ALLOCATE(KERNELS2%UO11(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO11 array'
         KERNELS2%UO11=0.
         endif
        
         if(keyEL .gt. 1) then
         ALLOCATE(KERNELS2%UO12(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO12 array'
         KERNELS2%UO12=0.
         endif
         
         if(keyEL .gt. 2) then
         ALLOCATE(KERNELS2%UO22(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO22 array'
         KERNELS2%UO22=0.
         endif
         
         if(keyEL .gt. 3) then
         ALLOCATE(KERNELS2%UO33(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO33 array'
         KERNELS2%UO33=0.
         endif
         
         if(keyEL .gt. 4) then 
         ALLOCATE(KERNELS2%UO34(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO34 array'
         KERNELS2%UO34=0.
         endif
         
         if(keyEL .gt. 5) then
         ALLOCATE(KERNELS2%UO44(KRpar,KMpar,KCpar,KIMpar,KREpar,KWLpar),stat=ierr)
         if(ierr/=0) stop 'Can not allocate KERNELS2%UO44 array'
         KERNELS2%UO44=0.
         endif 

      ENDIF ! key    
	        
      return
      end subroutine alloc_kernel_arrays

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

      subroutine dealloc_kernel_arrays(key,keyEL,KERNELS1,KERNELS2)
	  
      !use mod_par_DLS,     only: KN1par,KR1par,KM1par,KRpar,KMpar,KREpar,KIMpar
      !use mod_par_DLS_bin, only: KWLpar,KCpar
      
      implicit none	  	  
         
      integer,                     intent(in)     ::  key,keyEL
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2

      integer                                     ::  ierr

!c*********************************************************
!c ** FOR single scattering aerosol OPTICAL properties !!!*
!c
!c *** DEALLOCATE ARRAYS (alloc.mod) in subroutine USMATRIX 
!c *** (in matrix_intrpl.f)
!c
      IF(key .eq. 0 .or. key .eq. 3) THEN 
!c
!c *** DEALLOCATE ARRAYS (alloc1.mod)
!c
         DEALLOCATE(KERNELS1%UEA,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%UEA array'
         
         if(keyEL .gt. 0) then
         DEALLOCATE(KERNELS1%U11,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U11 array'
         endif
         
         if(keyEL .gt. 1) then
         DEALLOCATE(KERNELS1%U12,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U12 array'
         endif
         
         if(keyEL .gt. 2) then
         DEALLOCATE(KERNELS1%U22,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U22 array'
         endif
         
         if(keyEL .gt. 3) then
         DEALLOCATE(KERNELS1%U33,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U33 array'
         endif
         
         if(keyEL .gt. 4) then
         DEALLOCATE(KERNELS1%U34,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U34 array'
         endif
         
         if(keyEL .gt. 5) then
         DEALLOCATE(KERNELS1%U44,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS1%U44 array'
         endif
          
      ELSE 
!c
!c *** DEALLOCATE ARRAYS (alloc2.mod) in subroutine USMATRIX 
!c *** (in matrix_intrpl.wl.f)
!c
         DEALLOCATE(KERNELS2%UOEA,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UOEA array'
         if(keyEL .gt. 0) then        
         DEALLOCATE(KERNELS2%UO11,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO11 array'
         endif
        
         if(keyEL .gt. 1) then
         DEALLOCATE(KERNELS2%UO12,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO12 array'
         endif
         
         if(keyEL .gt. 2) then        
         DEALLOCATE(KERNELS2%UO22,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO22 array'
         endif
         
         if(keyEL .gt. 3) then
         DEALLOCATE(KERNELS2%UO33,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO33 array'
         endif
         
         if(keyEL .gt. 4) then
         DEALLOCATE(KERNELS2%UO34,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO34 array'
         endif
         
         if(keyEL .gt. 5) then
         DEALLOCATE(KERNELS2%UO44,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate KERNELS2%UO44 array'
         endif 
         
        ENDIF ! key   

      return
      end subroutine dealloc_kernel_arrays

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

END MODULE mod_alloc_kernels
