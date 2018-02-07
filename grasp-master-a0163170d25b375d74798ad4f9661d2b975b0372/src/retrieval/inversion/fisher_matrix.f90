! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains :

! SUBROUTINE FISHMX
! subroutine modified_normal_system
! subroutine UF1_matrix
! subroutine UF_matrix
!

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE FISHMX (                   &
                           KM,KN,           &
                           FS,FPS,CS,US,    & 
                           UFS,FFS          &
                        )
!C**************************************************
!C  THIS SUBROUTINE CALCULATES "FISHER MATRIX":
!C            (U)**T (C) (U)
!C                and
!C  "gradient" (U)**T (C) (F(p)-F(*))
!C                and
!C   "residual" (F(p)-F(*))**T (C) (F(p)-F(*))
!C**************************************************
!C  INPUT:
!C        KN  I        - number of lines
!C        KM  I        - number of columns
!C        US   R(KM,KN) - (derivatives) matrix
!C        CS   R(KM)    - matrix inverse to covariance
!C        FPS  R(KM)    - P-th approximation of vector F
!C        FS   R(KM)    - "measurements" vector 
!C  OUTPUT:
!C         UFS R(KN,KN) -"Fisher matrix" normalized by
!C                      maximum diagonal ellement
!C         FFS R(KN)    - "gradient" vector normalized by
!C                      maximum diagonal ellement of UF  
!C***************************************************
      USE mod_par_inv, only : KMESS,KPARS
	  
      IMPLICIT NONE	  

! ----------------------------------------------------------
! IN :
      INTEGER,                    INTENT(IN)  :: KM,KN	  
      REAL,DIMENSION(KMESS,KPARS),INTENT(IN)  :: US
      REAL,DIMENSION(KMESS),      INTENT(IN)  :: CS
      REAL,DIMENSION(KMESS),      INTENT(IN)  :: FS,FPS
! ----------------------------------------------------------
! OUT :
      REAL,DIMENSION(KPARS,KPARS),INTENT(OUT) :: UFS
      REAL,DIMENSION(KPARS),      INTENT(OUT) :: FFS
! ----------------------------------------------------------
! LOCAL :	  
      INTEGER :: I,I1
! ----------------------------------------------------------

!C*** calculating "Fisher matrix" 
!write(*,*) 'in FISHMX'
!      write(*,*) 'US(1:13,14)', US(1:13,14)
!write(*,*) 'CS', CS(1:13)
!write(*,*) 'KM, KN', KM, KN
      DO I=1,KN
      DO I1=1,KN
         UFS(I,I1)=SUM(CS(1:KM)*US(1:KM,I)*US(1:KM,I1))
      ENDDO ! I1
      ENDDO ! I

      DO I=1,KN
         FFS(I)=SUM(US(1:KM,I)*CS(1:KM)*(FPS(1:KM)-FS(1:KM)))
      ENDDO ! I
!write(*,*) 'UFS(14,1:32)', UFS(14,1:32)
      RETURN
      END SUBROUTINE FISHMX 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
	 
      subroutine UF_matrix (                              & ! IN
                              RIN,                        &
                              KN,npixels,UFS,             &
                              TCCORG,SMIM1,NISM,KSM,IMSM, &	
                              UF,nnz                      & ! OUT 
                           )
	 
      use mod_par_inv, only : KIMAGE,KMPSM,KPARS,KPAR
      use mod_retr_settings_derived_type

      implicit none

! ---------------------------------------------------	  
! IN :
      type(retr_input_settings),          intent(in)  :: RIN
      integer,                            intent(in)  :: KN,npixels,NISM      
      integer,dimension(KMPSM),           intent(in)  :: KSM
      integer,dimension(KMPSM,KIMAGE),    intent(in)  :: IMSM 
      real,dimension(KPARS,KPARS,KIMAGE), intent(in)  :: UFS
      real,                               intent(in)  :: TCCORG      
      real,dimension(KIMAGE,KPARS,KIMAGE),intent(in)  :: SMIM1
! ---------------------------------------------------
! OUT :
      real, dimension(KPAR,KPAR),         intent(out) :: UF
      integer,                            intent(out) :: nnz      
! ---------------------------------------------------
! LOCAL :
      integer :: ipix,I,I1,IS,IS1,IS2
! ---------------------------------------------------	  
      UF(:,:) = 0.0
      nnz = 0
      do ipix=1,npixels
        do I1=1,RIN%KNSING
          do I=1,RIN%KNSING
               UF((ipix-1)*RIN%KNSING+I,(ipix-1)*RIN%KNSING+I1) =  &
               UFS(I,I1,ipix)
               if(UFS(I,I1,ipix) .ne. 0) &
               nnz = nnz+1
          enddo ! I
        enddo ! I1
      enddo ! ipix      
      do IS=1,NISM
        do I=1,RIN%KNSING
          do IS1=1,KSM(IS)
            do IS2=1,KSM(IS)
              if(IS1 .ne. IS2) then
                 UF((IMSM(IS,IS1)-1)*RIN%KNSING+I,(IMSM(IS,IS2)-1)*RIN%KNSING+I) =  & 
                 TCCORG*SMIM1(IMSM(IS,IS1),I,IMSM(IS,IS2))
              endif
            enddo ! IS2
          enddo ! IS1      
        enddo ! I
      enddo ! IS

    return
    end subroutine UF_matrix 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
	 
      subroutine UF1_matrix ( KN,UF,          & ! IN
                              IP2,NIPPN,IPPN, &      
                              UF1             & ! OUT 
                            )
	 
      use mod_par_inv, only : KPAR

      implicit none

! ---------------------------------------------------	  
! IN :
      integer,                      intent(in)  :: KN,IP2      
      integer,dimension(KPAR),      intent(in)  :: NIPPN  
      integer,dimension(KPAR,KPAR), intent(in)  :: IPPN           
      real, dimension(KPAR,KPAR),   intent(in)  :: UF
! ---------------------------------------------------
! OUT :
      real, dimension(KPAR,KPAR),   intent(out) :: UF1
! ---------------------------------------------------
! LOCAL :
      real, dimension(:),   allocatable :: AATEM
      integer :: IIPP,I,I1,IP3,IP4,J,alloc_stat
      real    :: AA
! ---------------------------------------------------	  
      allocate(AATEM(KPAR),stat=alloc_stat)
      if (alloc_stat /= 0) stop 'error while trying to allocate AATEM'

      UF1(:,:) = 0.0
      IIPP = 0
      
      if(NIPPN(IP2+1) .gt. 0) then
        IIPP=NIPPN(IP2+1)
        do J=1,NIPPN(IP2+1)		   
          do I=1,NIPPN(IP2+1)
            UF1(I,J) = UF(IPPN(IP2+1,I),IPPN(IP2+1,J))
          enddo ! I
        enddo ! J
      endif
      do IP3=1,IP2     
        AATEM(I) = 0.0 
        do I=1,KN
          AA = 0.0
          do I1=1,NIPPN(IP3)
            AA = AA+UF(I,IPPN(IP3,I1))
          enddo
          AATEM(I) = AA
        enddo
        if(NIPPN(IP2+1) .gt. 0) then
          do I=1,NIPPN(IP2+1)
            UF1(I,NIPPN(IP2+1)+IP3) = AATEM(IPPN(IP2+1,I))
            UF1(NIPPN(IP2+1)+IP3,I) = AATEM(IPPN(IP2+1,I))
          enddo
        endif
        do IP4=1,IP2
          AA = 0.0
          do I1=1,NIPPN(IP4)
            AA = AA+AATEM(IPPN(IP4,I1))
          enddo
          UF1(IIPP+IP4,IIPP+IP3) = AA
        enddo
      enddo ! DO IP3=1,IP2

      deallocate(AATEM,stat=alloc_stat)
      if (alloc_stat /= 0) stop 'error while trying to deallocate AATEM'

      return
      end subroutine UF1_matrix 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! (Appendix C in Dubovik et al. 2011) - modifying the normal system 
!  for the case when some retrieval parameters are assumed the same
! Start

      subroutine modified_normal_system ( iu_main_output,RIN,KN, & ! IN
                                          npixels,IP2,NIPPN,IPPN,   &                                                
                                          test_UFNZ,UFNZ,UFNZ1,UF,FFMI0,QIM,  &
                                          solver_timer,                 &
                                          sparse_solver_name            & ! OUT 
                                        )
                                        
	 
      use mod_retr_settings_derived_type
      use mod_par_inv, only : KPAR
      use mod_fisher_matrix_ccs
      
      implicit none

! ---------------------------------------------------	  
! IN :
      integer,                      intent(in)  :: iu_main_output      
      integer,                      intent(in)  :: KN,IP2,npixels
      type(retr_input_settings),    intent(in)  :: RIN	

      integer,dimension(KPAR),      intent(in)  :: NIPPN  
      integer,dimension(KPAR,KPAR), intent(in)  :: IPPN           
      real, dimension(KPAR,KPAR),   intent(in)  :: UF
      real, dimension(KPAR),        intent(in)  :: FFMI0
      type(nonzero),                intent(inout)  :: UFNZ,UFNZ1
      logical,                      intent(in)  :: test_UFNZ
      character(*),                 intent(in)  :: sparse_solver_name

! ---------------------------------------------------
! OUT :
      real, dimension(KPAR),      intent(out)   :: QIM 
      real,                       intent(inout) :: solver_timer

! ---------------------------------------------------
! LOCAL :
      real, dimension(:),   allocatable :: AATEM
      integer :: IIPP,I,I1,IP3,IP4,J,alloc_stat
      real    :: AA
      real, dimension(:,:), allocatable :: UF1

! ---------------------------------------------------	  
      real, dimension(KPAR)      :: FFMI01,QIM1 
      !type(nonzero)              :: UFNZ1
      real*8, dimension(KPAR)    :: b	
      integer                    :: KNN1,nnz
! ---------------------------------------------------	  
      QIM(:) = 0.0

      allocate(UF1(KPAR,KPAR),stat=alloc_stat)
      if (alloc_stat /= 0) stop 'error while trying to allocate UF1 in modified_normal_system'

      IIPP=0

      IF(NIPPN(IP2+1) .GT. 0) THEN
         IIPP = NIPPN(IP2+1)
         DO I=1,NIPPN(IP2+1)
         FFMI01(I) = FFMI0(IPPN(IP2+1,I))
		   ENDDO ! I
      ENDIF ! NIPPN(IP2+1) .GT. 0

      DO IP3=1,IP2
         AA = 0.0
         DO I=1,NIPPN(IP3)
         AA = AA + FFMI0(IPPN(IP3,I))
         ENDDO ! I
		   FFMI01(IIPP+IP3) = AA        
      ENDDO ! IP3

! ******************
!tl  do not delete
   IF(RIN%IMQ .eq. 2 .or. test_UFNZ) THEN

      call UF1_matrix ( KN,UF,          & ! IN
                        IP2,NIPPN,IPPN, &      
                        UF1             & ! OUT 
                      )

      !if(test_UFNZ) then
        !write(*,*) 'IP2=',IP2,'  NIPPN(IP2+1)+IP2=',NIPPN(IP2+1)+IP2
        !do J=1,NIPPN(IP2+1)+IP2
          !if(j .eq. -365 .or. j .eq. -366) then
          !do I=1,NIPPN(IP2+1)+IP2
            !write(*,*) 'i=',i,'  j=',j,'  UF1(i,j)=',UF1(i,j)
          !enddo ! I
        !endif
        !i=count(UF1(1:NIPPN(IP2+1)+ip2,j) .ne. 0 )
        !write(*,*) 'US1 test: j=',j,'  #nz=',i
        !enddo !j
      !endif ! test_UFNZ
 
   ENDIF ! RIN%IMQ .eq. 2 .or. 
!tl  do not delete to

      KNN1 = IIPP+IP2
      if(RIN%IPRI_verbose)  write(iu_main_output,*) 'KNN1=',KNN1
 
! ******************

!tl do not delete
      IF(RIN%IMQ .eq. 2) THEN
      if(RIN%IPRI_verbose)  write(iu_main_output,*) 'UF1 inversion by by SVD'  
         CALL ITERQ (                              &
                     2,KNN1,                       &
                     UF1(1:KNN1,1:KNN1),           &
                     FFMI01(1:KNN1),               &
                     QIM1(1:KNN1)                  & ! OUT
                    )						 
      ENDIF !RIN%IMQ .eq. 2
!tl  do not delete to
              write(*,*) '1:   IP2=',IP2

      IF(RIN%IMQ .EQ. 3) THEN
         CALL UF1_nonzero (                        &
                           RIN%KNSING,npixels,     & ! IN
                           IP2,NIPPN,IPPN,         &
                           UFNZ,                   &
                           UFNZ1,nnz               & ! OUT
                          )

         b(1:KNN1) = FFMI01(1:KNN1)

if(test_UFNZ) then
UF1(:,:) = 0.0
do i=1, nnz
UF1(UFNZ1%row(i),UFNZ1%col(i))=UFNZ1%val(i)
enddo ! i

!write(*,*) 'IP2=',IP2,'  NIPPN(IP2+1)+IP2=',NIPPN(IP2+1)+IP2
do J=1,NIPPN(IP2+1)+IP2
if(j .eq. -365 .or. j .eq. -366) then
do I=1,NIPPN(IP2+1)+IP2
   !write(*,*) 'i=',i,'  j=',j,'  UF1(i,j)=',UF1(i,j)
enddo ! I
endif
i=count(UF1(1:NIPPN(IP2+1)+ip2,j) .ne. 0 )
   !write(*,*) 'US1 test2: j=',j,'  #nz=',i
enddo !j
endif ! test_UFNZ

         if(RIN%IPRI_verbose)  write(iu_main_output,*)   & 
         'UF1 inversion, before sparse_solver ',trim(sparse_solver_name),'  nnz=',nnz  

         call sparse_solver ( iu_main_output,  & ! IN 
	                            KNN1,nnz,UFNZ1,     & 
                              b,                  & ! INOUT
                              solver_timer        & 
                            )     
         if(RIN%IPRI_verbose)  write(iu_main_output,*)   & 
         'UF1 inversion, after  sparse_solver ',trim(sparse_solver_name)
	   					
         QIM1(1:KNN1) = b(1:KNN1)

      ENDIF ! RIN%IMQ .EQ. 3

! ***************************************

      IF(NIPPN(IP2+1) .GT. 0) THEN
         DO IP3=1,IP2+1
         DO I=1,NIPPN(IP3)
         IF(IP3 .LT. IP2+1) THEN
         QIM(IPPN(IP3,I)) = QIM1(NIPPN(IP2+1)+IP3)
         ELSE
         QIM(IPPN(IP3,I)) = QIM1(I)
         ENDIF ! IP3 .LT. IP2+1
         ENDDO ! I
         ENDDO ! IP3
      ELSE
         DO IP3=1,IP2
         DO I=1,NIPPN(IP3)
         QIM(IPPN(IP3,I)) = QIM1(IP3)
         ENDDO ! I
         ENDDO ! IP3
      ENDIF ! NIPPN(IP2+1) .GT. 0

! End
!  (Appendix C in Dubovik et al. 2011) - modifying the normal system 
!  for the case when some retrieval parameters are assumed the same

      deallocate(UF1,stat=alloc_stat)
      if (alloc_stat /= 0) stop 'error while trying to deallocate UF1 in modified_normal_system'

      return
      end subroutine modified_normal_system 

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

