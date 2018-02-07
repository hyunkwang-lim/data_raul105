! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains:


! SUBROUTINE RESIDUAL
! subroutine residual_meas_term
! subroutine gradient_apriori_smooth_term
! subroutine gradient_apriori_estim_term
! subroutine residual_apriori_smooth_term
! subroutine residual_apriori_estim_term
! subroutine residual_inter_pixel_term
! subroutine residual_details

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual(IPRI,iu,KN,F,FP,AA)
!C******************************************************
!C*** Subroutine for calculating measurement fit term residual ***
!C***
!C
!C   INPUT:
!C      KM     I       - number of "measurements" 
!C      F     R(KM)    - vector of "measurements"
!C      FP    R(KM)    - vector of fitted "measurements"
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - residual F^T FP
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)             ::  KN,iu
      logical,intent(in)             ::  IPRI
      real,dimension(KN),intent(in)  ::  F,FP
!------------------------------------------------------------
! OUT :
      real,intent(out)                 :: AA
!------------------------------------------------------------
      AA = 0.0      
      AA = DOT_PRODUCT(F(1:KN),FP(1:KN))
      if(AA .lt. 0.) then
        if(IPRI) write(iu,*) AA,'-AA, negative, if(AA .lt. 0.) AA=0.; in residual !!!'
        AA = 0. ! OD
      endif ! AA

      return 
      end subroutine residual
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual_meas_term(IPRI,iu,KM,F,FP,CS,AA)
!C******************************************************
!C*** Subroutine for calculating measurement fit term residual ***
!C***
!C
!C   INPUT:
!C      KM     I       - number of "measurements" 
!C      F     R(KM)    - vector of "measurements"
!C      FP    R(KM)    - vector of fitted "measurements"
!C      CS    R(KM)    - vector of weights
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - residual (F-FP)^T C (F-FP)
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)             ::  KM,iu
      logical,intent(in)             ::  IPRI
      real,dimension(KM),intent(in)  ::  F,FP,CS
!------------------------------------------------------------
! OUT :
      real,intent(out)                 :: AA
!------------------------------------------------------------
! LOCAL ::
      real,dimension(KM)               :: DIF,temp
!------------------------------------------------------------	
      AA = 0.0      
      DIF(1:KM) = FP(1:KM)-F(1:KM)
      temp(1:KM) = CS(1:KM)*DIF(1:KM)      
      AA = DOT_PRODUCT(DIF(1:KM),temp(1:KM))
      if(AA .lt. 0.) then
        if(IPRI) write(iu,*) AA,'-AA, negative, if(AA .lt. 0.) AA = 0.; in residual_meas_term !!!'
        AA = 0.
      endif ! AA
      
      return 
      end subroutine residual_meas_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine gradient_apriori_smooth_term(KN,APS,SMSING,FFS)
!C******************************************************
!C*** Subroutine for calculating smoothness term of residual ***
!C***
!C
!C   INPUT:
!C      KN      I        - number of parameters 
!C      APS     R(KN)    - vector of parameters
!C      SMSING  R(KN,KN) - single pixel matrix of smoothness 
!C
!C***
!C  
!C   OUTPUT:
!C       FFS     R      - smoothness term in gradient of minimized residual
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)               :: KN
      real,dimension(KN),intent(in)    :: APS
      real,dimension(KN,KN),intent(in) :: SMSING
!------------------------------------------------------------
! OUT :
      real,dimension(KN)               :: FFS
!------------------------------------------------------------
      FFS(1:KN) = 0.0
      FFS(1:KN) = MATMUL(SMSING(1:KN,1:KN),APS(1:KN))

      return 
      end subroutine gradient_apriori_smooth_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine gradient_apriori_estim_term(KN,APS,APS0,SMSING0,FFS)
!C******************************************************
!C*** Subroutine for calculating smoothness term of residual ***
!C***
!C
!C   INPUT:
!C      KN       I        - number of parameters 
!C      APS      R(KN)    - vector of parameters
!C      APS0     R(KN)    - vector of a priori estimate parameters
!C      SMSING0  R(KN,KN) - single pixel matrix of smoothness 
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - residual (APS-APS0)^T SMSING0 (APS-APS0)
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)               :: KN
      real,dimension(KN),intent(in)    :: APS,APS0
      real,dimension(KN,KN),intent(in) :: SMSING0
!------------------------------------------------------------
! OUT :
      real,dimension(KN)               :: FFS
!------------------------------------------------------------
! LOCAL ::
      real,dimension(KN)               :: DIF
!------------------------------------------------------------	
      DIF(1:KN) = 0.0
      FFS(1:KN) = 0.0
      
      DIF(1:KN) = APS(1:KN)-APS0(1:KN)
      FFS(1:KN) = MATMUL(SMSING0(1:KN,1:KN),DIF(1:KN))
      
      return 
      end subroutine gradient_apriori_estim_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual_apriori_smooth_term(IPRI,iu,KN,APS,SMSING,AA)
!C******************************************************
!C*** Subroutine for calculating smoothness term of residual ***
!C***
!C
!C   INPUT:
!C      KN      I        - number of parameters 
!C      APS     R(KN)    - vector of parameters
!C      SMSING  R(KN,KN) - single pixel matrix of smoothness 
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - residual APS^T SMSING APS
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)               :: KN,iu
      logical,intent(in)               :: IPRI
      real,dimension(KN),intent(in)    :: APS
      real,dimension(KN,KN),intent(in) :: SMSING
!------------------------------------------------------------
! OUT :
      real,intent(out)                 :: AA
!------------------------------------------------------------
! LOCAL ::
      real,dimension(KN)               :: FFS
!------------------------------------------------------------	
      AA = 0.0      
      FFS(1:KN) = MATMUL(SMSING(1:KN,1:KN),APS(1:KN))
      AA = DOT_PRODUCT(FFS(1:KN),APS(1:KN))
      if(AA .lt. 0.) then
        if(IPRI) write(iu,*) AA,'-AA, negative, if(AA .lt. 0.) AA=0.; in residual_apriori_smooth_term !!!'
        AA = 0. ! OD
      endif ! AA

      return 
      end subroutine residual_apriori_smooth_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual_apriori_estim_term(IPRI,iu,KN,APS,APS0,SMSING0,AA)
!C******************************************************
!C*** Subroutine for calculating smoothness term of residual ***
!C***
!C
!C   INPUT:
!C      KN       I        - number of parameters 
!C      APS      R(KN)    - vector of parameters
!C      APS0     R(KN)    - vector of a priori estimate parameters
!C      SMSING0  R(KN,KN) - single pixel matrix of smoothness 
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - residual (APS-APS0)^T SMSING0 (APS-APS0)
!C*******************************************************

      implicit none
!------------------------------------------------------------
! IN :
      integer,intent(in)               :: KN, iu
      logical,intent(in)               :: IPRI
      real,dimension(KN),intent(in)    :: APS, APS0
      real,dimension(KN,KN),intent(in) :: SMSING0
!------------------------------------------------------------
! OUT :
      real,intent(out)                 :: AA
!------------------------------------------------------------
! LOCAL ::
      real,dimension(KN)               :: FFS,DIF
!------------------------------------------------------------	
      AA = 0.0      
      DIF(1:KN) = APS(1:KN)-APS0(1:KN)
      FFS(1:KN) = MATMUL(SMSING0(1:KN,1:KN),DIF(1:KN))
      AA = DOT_PRODUCT(FFS(1:KN),DIF(1:KN))
      if(AA .lt. 0.) then
        if(IPRI) write(iu,*) AA,'-AA, negative, if(AA .lt. 0.) AA=0.; in residual_apriori_estim_term !!!'
        AA = 0. ! OD
      endif ! AA

      return 
      end subroutine residual_apriori_estim_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual_inter_pixel_term(IPRI,iu,KN,npixels,APS,SMIM1,ledges,FF_edge,AB_edge,AA)
!C******************************************************
!C*** Subroutine for calculating smoothness term of residual ***
!C***
!C
!C   INPUT:
!C      KN      I        - number of parameters 
!C      APS     R(KN)    - vector of parameters
!C      SMIM1   R(KIMAGE,KPARS,KIMAGE) - matrix of smoothness 
!C
!C***
!C  
!C   OUTPUT:
!C       AA     R      - 
!C      It may include "inter-pixel smoothness term" +
!C                     "inter-pixel edge_term"
!C      inter-pixel smoothness term: (AP)^T (D_inter)^T(D_inter) AP)
!C                                     See  Dubovik et al. [2011]
!C      inter-pixel edge_term:
!C                   ((D1)* AP - (D2) * AB)^T ((D1)* AP - (D2) * AB)
!C                    AP - values outside of segment near edges
!C      this is equivalent to:
!C            (AP)^T (D1)^T (D1) AP - 2 (AP)^T (D1)^T (D2) AB + (AB)^T (D2)^T (D2) AB:
!C
!C           - (AP)^T (D1)^T (D1) AP is included in (AP)^T (D_inter)^T(D_inter) AP)
!C           - FF_edges = (D1)^T (D2) AB
!C           - AB_edge = (AB)^T (D2)^T (D2) AB:
!C***************************************************************************************
      use mod_par_inv, only : KIMAGE,KPARS
            
      implicit none
!------------------------------------------------------------
! IN :
      integer,                     intent(in)  ::  KN,iu,npixels
      logical,                     intent(in)  ::  IPRI
      real,dimension(KPARS,KIMAGE),intent(in)  ::  APS,FF_edge
      real,dimension(KIMAGE,KPARS,KIMAGE),intent(in) ::  SMIM1
      real,                        intent(in)   :: AB_edge
      logical,                     intent(in)   :: ledges

!------------------------------------------------------------
! OUT :
      real,intent(out)       :: AA  !,AA_edge
!------------------------------------------------------------
! LOCAL ::
      integer                :: IS1,IS2,i
      real,dimension(KN)     :: temp
      real                   :: AA_edge  ! ask Oleg
!------------------------------------------------------------	
      AA = 0.0
      temp(:) = 0.0
!C      inter-pixel smoothness term: (AP)^T (D_inter)^T(D_inter) AP)
!C  *NOTE(:       this term includes edge term (AP)^T (D1)^T (D1) AP
      do IS1=1,npixels
         do IS2=1,npixels
            temp(1:KN) = SMIM1(IS1,1:KN,IS2)*APS(1:KN,IS2)
            AA = AA+DOT_PRODUCT(APS(1:KN,IS1),temp(1:KN))
            !AA = AA+SUM(APS(1:KN,IS1)*SMIM1(IS1,1:KN,IS2)*APS(1:KN,IS2))
!         do i=1,KN
!         if(SMIM1(IS1,I,IS2) .ne. 0.) then
!         write(*,*) '** IS1,IS2,i,SMIM1(IS1,i,IS2),APS(I,IS1),APS(i,IS2),AA: ', &
!                        IS1,IS2,i,SMIM1(IS1,I,IS2),APS(I,IS1),APS(I,IS2),AA,    &
!                        '  in residual_inter_pixel_term'
!         write(*,*) 'temp:'
!         write(*,'(10e14.4)') temp(1:KN)
!         endif
!         enddo ! i
         enddo ! IS2
         !write(*,*) '** IS1,AA: ',IS1,AA,'  in residual_inter_pixel_term'
      enddo ! IS1  
!write(*,*) 'test: AA=',AA
!C************************************************************************************
      if ( ledges ) then
        AA_edge = 0.0
!C***    Calculaton of part of edge term: (AP)^T (D1)^T (D2) AB
        do IS1=1,npixels
          AA_edge = AA_edge+DOT_PRODUCT(FF_edge(1:KN,IS1),APS(1:KN,IS1))
!write(*,*) 'IS1=',IS1,'  FF_edge(1:KN,IS1): '
!write(*,'(10e14.4)')  FF_edge(1:KN,IS1)
!write(*,*) 'IS1=',IS1,'  APS(1:KN,IS1): '
!write(*,'(10e14.4)')  EXP(APS(1:KN,IS1))

        enddo ! IS1
!C adding:  [- 2 (AP)^T (D1)^T (D2) AB + (AB)^T (D2)^T (D2) AB]  
!C         to (AP)^T (D_inter)^T(D_inter) AP)
!write(*,*) 'residual_inter_pixel_term test: AA=',AA,' AA_edge=',AA_edge,' AB_edge=',AB_edge
!od and tl        AA = AA-2.*AA_edge+AB_edge
        AA = AA+2.*AA_edge+AB_edge
!od&tl        AA = AA+AA_edge+AB_edge
      endif ! ledges
!write(*,*) 'residual_inter_pixel_term test: AA=',AA
if(AA .ne. AA) stop '!!!!!!! in residual_inter_pixel_term: NaN was detected'
      if(AA .lt. 0.) then
        if(IPRI) write(iu,*) AA,'-AA, if(AA .lt. 0.) AA = 0.; residual_inter_pixel_term !!!'
        AA = 0. ! OD
      endif ! AA

      return
      end subroutine residual_inter_pixel_term

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine residual_details (                         & 
                                    RIN,                    &
                                    NIMAGE,IIMAGE,          & 
                                    IKI,IKI_shift,KNOISEI,  & 						 
                                    segment_vec_meas,       &
                                    segment_vec_fit,        &
                                    RESA,RESR,RESAT,RESRT   & ! OUT 		
                                  )

! *****************************************************************
! ***  Accounting for different accuracy levels in the data     ***
! ***       AND   modeling   RANDOM NOISE                       ***
! *****************************************************************
! *** INOISE  - the number of different noise sources           ***
! *** SGMS(I) - std of noise in i -th source                    ***
! *** INN(I)  - EQ.1.THEN error is absolute with                ***
! ***         - EQ.0 THEN error assumed relative                ***
! *** DNN(I)  - variation of the noise of the I-th source       ***
! *** IK(I)   - total number of measurements of i-th source     ***
! *** KNOISE(1,K)-specific numbers affected by i-th source      ***
! *** All the measurements which where not listed in  the INOISE-1 
! ***  first INOISE-1 sources, they belong to last source of noise 
! *****************************************************************

      use mod_par_inv, only : KIMAGE,KKNOISE,KMESS
      use mod_retr_settings_derived_type
      use mod_sdata
	  
      implicit none
! -----------------------------------------------------------------
! IN :
      type(retr_input_settings),               intent(in) :: RIN
      integer,                                 intent(in) :: NIMAGE, IIMAGE
      integer,dimension(KKNOISE,KIMAGE),       intent(in) :: IKI
      logical,dimension(KKNOISE,KIMAGE),       intent(in) :: IKI_shift
      integer,dimension(KMESS,KKNOISE,KIMAGE), intent(in) :: KNOISEI
      type(pixel_vector),dimension(KIMAGE),    intent(in) :: segment_vec_meas
      type(pixel_vector),dimension(KIMAGE),    intent(in) :: segment_vec_fit

! -----------------------------------------------------------------
! OUT :
      real,dimension(KKNOISE),                intent(out) :: RESA, RESR
      real,dimension(KKNOISE),                intent(out) :: RESAT,RESRT
! -----------------------------------------------------------------
! LOCAL : 
      integer                      :: I,IIN
      real                         :: DIFA,DIFR,tiny,TEMP
      real,dimension(KKNOISE),save :: TEMPA = 0.0,TEMPR = 0.0
! -----------------------------------------------------------------

      tiny    = 1e-6
      RESA(:) = 0.0
      RESR(:) = 0.0

      IF(IIMAGE .EQ. 1) THEN 
        TEMPA(:) = 0.0
        TEMPR(:) = 0.0
      ENDIF ! IIMAGE.EQ.1

      SELECT CASE(RIN%KL)
      CASE(0) ! absol.	

         DO IIN=1,RIN%NOISE%INOISE
         SELECT CASE(IKI_shift(IIN,IIMAGE))
         CASE(.false.) ! no SHIFT		 
            DO I=1,IKI(IIN,IIMAGE)
            DIFA = segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)) -      &
			          segment_vec_fit (IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)) 
            RESA(IIN)=RESA(IIN)+DIFA*DIFA
            DIFR = DIFA /  &
				       segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))
            RESR(IIN)=RESR(IIN)+DIFR*DIFR
            ENDDO ! I		 
            TEMPA(IIN)=TEMPA(IIN)+RESA(IIN)
            TEMPR(IIN)=TEMPR(IIN)+RESR(IIN)
         CASE(.true.) ! SHIFT added
            DO I=1,IKI(IIN,IIMAGE)
            DIFA = segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)) -      &
                   segment_vec_fit (IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)) 					
            RESA(IIN)=RESA(IIN)+DIFA*DIFA
            !TEMP=segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))-RIN%SHIFT            
            TEMP=( (segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))-RIN%SHIFT) +  &
                   (segment_vec_fit(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))-RIN%SHIFT) )*0.5            
            IF(ABS(TEMP) .LT. tiny)  TEMP = tiny  
            DIFR = DIFA/TEMP
            RESR(IIN)=RESR(IIN)+DIFR*DIFR
            ENDDO ! I		 
            TEMPA(IIN)=TEMPA(IIN)+RESA(IIN)
            TEMPR(IIN)=TEMPR(IIN)+RESR(IIN)
         END SELECT ! IIN
         ENDDO ! IIN 

      CASE(1) ! ln

         DO IIN=1,RIN%NOISE%INOISE

         SELECT CASE(IKI_shift(IIN,IIMAGE))
         CASE(.false.) ! no SHIFT		 
            DO I=1,IKI(IIN,IIMAGE)
            DIFA = EXP(segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) - &
                   EXP(segment_vec_fit(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) 
            RESA(IIN)=RESA(IIN)+DIFA*DIFA
            TEMP=( EXP(segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)))+    &  ! OD
                   EXP(segment_vec_fit(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) )*0.5	                                            
            DIFR = DIFA/TEMP  
            RESR(IIN)=RESR(IIN)+DIFR*DIFR
            ENDDO ! I		 
            TEMPA(IIN)=TEMPA(IIN)+RESA(IIN)
            TEMPR(IIN)=TEMPR(IIN)+RESR(IIN)
         CASE(.true.) ! SHIFT added
            DO I=1,IKI(IIN,IIMAGE)
            DIFA = EXP(segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) - &
                   EXP(segment_vec_fit(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) 					
            RESA(IIN)=RESA(IIN)+DIFA*DIFA
            !TEMP=EXP(segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)))-RIN%SHIFT 
            TEMP=( (EXP(segment_vec_meas(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE)))-RIN%SHIFT) +  & ! OD
                   (EXP(segment_vec_fit(IIMAGE)%FS(KNOISEI(I,IIN,IIMAGE))) -RIN%SHIFT) )*0.5
            IF(ABS(TEMP) .LT. tiny)  TEMP = tiny  
            DIFR = DIFA/TEMP
            RESR(IIN)=RESR(IIN)+DIFR*DIFR
            ENDDO ! I		 
            TEMPA(IIN)=TEMPA(IIN)+RESA(IIN)
            TEMPR(IIN)=TEMPR(IIN)+RESR(IIN)
         END SELECT ! IIN
         ENDDO ! IIN 

      END SELECT ! RIN%KL

      DO IIN=1,RIN%NOISE%INOISE
         IF(IKI(IIN,IIMAGE) .GT. 0) THEN
         RESA (IIN) = SQRT(RESA(IIN)/IKI(IIN,IIMAGE))
         RESR (IIN) = SQRT(RESR(IIN)/IKI(IIN,IIMAGE))
!tl         ELSE
!tl         WRITE(*,*) 'IIMAGE=',IIMAGE,'  IIN=',IIN,'  IKI=',IKI(IIN,IIMAGE),' IKI can not be 0' 
!tl         STOP 'STOP in residual_details' 
         ENDIF ! IKI(IIN,IIMAGE) .gt. 0
      ENDDO ! IIN

      IF(IIMAGE .EQ. NIMAGE) THEN
      DO IIN=1,RIN%NOISE%INOISE
            RESAT(IIN) = SQRT(TEMPA(IIN)/SUM(IKI(IIN,1:NIMAGE)))
            RESRT(IIN) = SQRT(TEMPR(IIN)/SUM(IKI(IIN,1:NIMAGE)))
      ENDDO ! IIN 
      ENDIF ! IIMAGE .EQ. NIMAGE
	   
!	   write(*,*) RESA,IIMAGE	   
!	   write(*,*) RESR,IIMAGE	   

!	   write(*,*)
!	   write(*,*) RESAT	   
!	   write(*,*) RESRT	   

!	   write(*,*)
!stop 'stop in residual_details'	   
	   
      return
      end subroutine residual_details
 
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
