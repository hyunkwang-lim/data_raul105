! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!> @file smoothness.f90
        !> File contains routines related to smoothness constrains.
        !>
        !> @authors Oleg Dubovik and Tatsiana Lapionak
        !> 

! file contains:
! subroutine inter_pix_smooth_constr_incl
! subroutine ordsing_for_smooth_single_pixel
! SUBROUTINE smoothterm_single_pixel_apriori
! SUBROUTINE smoothterm_single_pixel_smoothness
! SUBROUTINE smoothterm_multi_pixel
! SUBROUTINE SMOOM
! SUBROUTINE DIFERW
! RECURSIVE  SUBROUTINE CHECK
! SUBROUTINE time_space_group
! SUBROUTINE smoothterm_mutli_pixel_edges
! SUBROUTINE MAT_T_MAT

#include "../constants_set/mod_globals.inc"
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
        !> @brief Routine calculates the change of quadratic form gradient at p - interation
        !> @brief for the entire segment of pixels with inclusion of inter-pixel constraints.
        !> 
        !!> @param[in]      KNSING  - number of parameters driving forwar model
        !!> @param[in]      npixels - number of inverted pixels
        !!> @param[in]      UFS     - Fisher matrix 
        !!> @param[in]      QS      - 
        !!> @param[in]      FFMS    - gradient: contribution of measurements
        !!> @param[in]      FFMS0   - contribution of a priroi estimates to the gradient
        !!> @param[in]      SMMULTI - multipixel smoothness matrix
        !!> @param[in]      NISM    - number of groups in which smoothness can be implemented
        !!> @param[in]      KSM     - number of elements in each of NISM groups
        !!> @param[in]      IMSM    - pixel numbers which belong to the group
        !!> @param[in]      ccor    -
        !!> @param[out]     FFMSQ   - 
        !>
               
      subroutine inter_pix_smooth_constr_incl ( KNSING,npixels,           &
                                                UFS,QS,FFMS,FFMS0,        &
                                                SMMULTI,NISM,KSM,IMSM,ccor, &
                                                FFMSQ                     &
                                              )

      use mod_par_inv, only : KPARS,KIMAGE,KMPSM
            
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!***********************************************************************************
!> This subroutine calculates the change of quadratic form gradient at p - interation
!> for the entire segment of pixels with inclusion of inter-pixel constraints :
!
! the normal system: gradient = (UF) QS - Fp
!
!
! where:  UF = (Up)^T W^(-1) (Up) + ((D_singl)^T(D_singl)) + Wa^(-1) + (D_inter)^T (D_inter)
!         FP =  (UF) Ap - [(Up)^T W^(-1)(fp - f) + Wa^(-1)(Ap-A0) + (D_inter)^T (D_inter)Ap + AB_edge]
!
!     *NOTE* we represent: UF=UF_diag + UF_non_diag, i.e.:
!         (UF) QS = (UF_diag + UF_non_diag) QS
!         FP = FP1 + FP2; 
!            - FP1- all contribution used already in sinngle pixel retrieval;
!              FP1 =(Up)^T W^(-1)(fp - f)+Wa^(-1)(Ap-A0) +((D_singl)^T(D_singl))Ap;
!            - FP2- the contribution using multi-pixel constraints and "edge" information
!              FP2 = (D_inter)^T (D_inter)Ap + FF_edge
!***********************************************************************************
      implicit none
! -------------------------------------------------------------------------------        
      integer,                            intent(in)  :: KNSING,npixels,NISM
      real,                               intent(in)  :: ccor
      real,dimension(KPARS,KPARS,KIMAGE), intent(in)  :: UFS
      integer,dimension(KMPSM),           intent(in)  :: KSM
      integer,dimension(KMPSM,KIMAGE),    intent(in)  :: IMSM 
      real,dimension(KIMAGE,KPARS,KIMAGE),intent(in)  :: SMMULTI
      real,dimension(KPARS,KIMAGE),       intent(in)  :: QS
      real,dimension(KPARS,KIMAGE),       intent(in)  :: FFMS,FFMS0
      real,dimension(KPARS,KIMAGE),       intent(out) :: FFMSQ
! -------------------------------------------------------------------------------  
      integer                      :: ipix,I,IS,IS1,IS2
      real                         :: TAAQ1,TAAQ2,TAAQ3,TAAQ4
      real,dimension(KPARS,KIMAGE) :: FMSU0

! -------------------------------------------------------------------------------
      FMSU0(:,:) = 0.0
      do IS=1,NISM
        do I=1,KNSING
          do IS1=1,KSM(IS)
            do IS2=1,KSM(IS)
              if(IS1 .ne. IS2) FMSU0(I,IMSM(IS,IS1))=FMSU0(I,IMSM(IS,IS1))+  &          ! (UF_non_diag) QS
                               ccor*SMMULTI(IMSM(IS,IS1),I,IMSM(IS,IS2))*QS(I,IMSM(IS,IS2))
            enddo ! IS2
          enddo ! IS1      
        ENDDO ! I
      enddo ! IS

      FFMSQ(:,:) = 0.0
      do ipix=1,npixels
        do I=1,KNSING
          TAAQ1 = DOT_PRODUCT(UFS(I,1:KNSING,ipix),QS(1:KNSING,ipix)) !(UF_diag) QS
          TAAQ2 = FMSU0(I,ipix) ! (UF_non_diag) QS
          TAAQ3 = FFMS0(I,ipix) ! FP2
          TAAQ4 = FFMS(I,ipix)  ! FP1
          FFMSQ(I,ipix) = TAAQ1+TAAQ2-TAAQ3-TAAQ4
!write(100,*) 'ipix= ',ipix,' I= ',I
!write(100,*) TAAQ1,TAAQ2,TAAQ3,TAAQ4,' - TAAQ1,TAAQ2,TAAQ3,TAAQ4'
        enddo ! I 
      enddo  ! ipix

      return
      end subroutine inter_pix_smooth_constr_incl

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss      
        !!> @brief Routine calculates 
        !> 
        !!> @param[in]      KNSING  - number of parameters driving forwar model
        !!> @param[in]      npixels - number of inverted pixels
        !!> @param[in]      UFS     - Fisher matrix
        !!> @param[in]      SMMULTI - multipixel smoothness matrix
        !!> @param[in]      NISM    - number of groups in which smoothness can be implemented
        !!> @param[in]      KSM     - number of elements in each of NISM groups
        !!> @param[in]      IMSM    - pixel numbers which belong to the group
        !!> @param[in]      ccor    -
        !!> @param[out]     AAI     - 
        !>
         
      subroutine matrix_Q_iter (KNSING,npixels,UFS,SMMULTI,NISM,KSM,IMSM,ccor,AAI)

      use mod_par_inv, only : KPARS,KIMAGE,KMPSM
            
      implicit none
! -------------------------------------------------------------------------------        
      integer,                            intent(in)  :: KNSING,npixels,NISM
      real,                               intent(in)  :: ccor
      real,dimension(KPARS,KPARS,KIMAGE), intent(in)  :: UFS
      integer,dimension(KMPSM),           intent(in)  :: KSM
      integer,dimension(KMPSM,KIMAGE),    intent(in)  :: IMSM 
      real,dimension(KIMAGE,KPARS,KIMAGE),intent(in)  :: SMMULTI
      real,dimension(KPARS,KIMAGE),       intent(out) :: AAI
! -------------------------------------------------------------------------------  
      integer :: ipix,I,I1,IS,IS1,IS2
! -------------------------------------------------------------------------------  
      
      AAI(:,:) = 0.0
      do ipix=1,npixels
        do I=1,KNSING
          do I1=1,KNSING
            AAI(I,ipix) = AAI(I,ipix)+ABS(UFS(I,I1,ipix))
          enddo
        enddo
        do IS=1,NISM
          do I=1,KNSING
            do IS1=1,KSM(IS)
              do IS2=1,KSM(IS)
                if(IS1 .NE. IS2) AAI(I,IMSM(IS,IS1))=AAI(I,IMSM(IS,IS1))+  &
                                 ccor*ABS(SMMULTI(IMSM(IS,IS1),I,IMSM(IS,IS2)))
              enddo
            enddo      
          enddo
        enddo
      enddo ! ipix 

      return
      end subroutine matrix_Q_iter

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates ordinates to be used for non-equal binning 
        !!> @brief (aerosol particle radii and wave lengths) in smoothness.
        !> 
        !!> @param[in]      RIN     - invertion setting structure
        !!> @param[inout]   ORDSING (KIDIM3,KIDIM2,KIDIM1) is used for non-equal binning (aerosol particle radii and wave lengths) in smoothness.
        !
        
      subroutine ordsing_for_smooth_single_pixel(RIN,ORDSING)
      
! ORDSING will be used for non-equal binning (aerosol particle radii and wave lengths) in smoothness.

      use mod_retr_settings_derived_type 
      use mod_par_inv, only : KIDIM1,KIDIM2,KIDIM3
      
      implicit none 
!	-------------------------------------------------------------------------
      type(retr_input_settings),            intent(in)     :: RIN
      real,dimension(KIDIM3,KIDIM2,KIDIM1), intent(inout)  :: ORDSING

!	-------------------------------------------------------------------------
      integer   ::  IDIM1,IDIM2,IDIM3,par_type
!	-------------------------------------------------------------------------
!	-------------------------------------------------------------------------
      ORDSING(:,:,:) = 0.0 

      do IDIM1=1,RIN%NDIM%n1
      par_type = RIN%NDIM%par_type(IDIM1)
      do IDIM2=1,RIN%NDIM%n2(IDIM1) 
      do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
        if(par_type .gt. par_type_SD_beg .and. par_type .lt. par_type_SD_end) then
          if(RIN%IBIN .eq. 1)  ORDSING(IDIM3,IDIM2,IDIM1) = RIN%RADIUS1(IDIM3,IDIM2)
          if(RIN%IBIN .eq. -1) ORDSING(IDIM3,IDIM2,IDIM1) = LOG(RIN%RADIUS1(IDIM3,IDIM2)) 
        elseif(par_type .gt. par_type_RERI_beg .and. par_type .lt. par_type_IMRI_end) then
!OD       ORDSING(IDIM3,IDIM2,IDIM1)=RIN%WAVE(IDIM3)
          ORDSING(IDIM3,IDIM2,IDIM1)=LOG(RIN%WAVE(IDIM3))
        elseif(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then
          ORDSING(IDIM3,IDIM2,IDIM1)=RIN%RATIO1(IDIM3,IDIM2)
        !elseif(par_type .gt. par_type_AVP_beg  .and. par_type .lt. par_type_AVP_end) then
!OD       ORDSING(IDIM3,IDIM2,IDIM1)=RIN%WAVE(IDIM3)
        elseif(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF_water_end) then
!OD       ORDSING(IDIM3,IDIM2,IDIM1)=RIN%WAVE(IDIM3)
          ORDSING(IDIM3,IDIM2,IDIM1)=LOG(RIN%WAVE(IDIM3))
        endif ! par_type .gt. par_type_SD_beg .and.
      enddo  ! IDIM3
      enddo  ! IDIM2   
      enddo  ! IDIM1

      return
      end subroutine ordsing_for_smooth_single_pixel

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine determines smoothness term for single pixel (a priori estimates)
        !>
        !> @param[in]    RIN          - invertion setting structure
        !> @param[in]    SPCA         -
        !> @param[out]   KSM          - 
        !> @param[out]   SMSING       - a priori estimates term of smoothness matrix for single pixel
        !> @param[out]   apriori_present - presence of a priori estimates
        !>

      SUBROUTINE smoothterm_single_pixel_apriori ( RIN, SPCA,         & ! IN
                                                   KSM, SMSING, apriori_present  & ! OUT
                                                 )
!************************************************************
!*** Determining smoothness term for single pixel         ***
!************************************************************
      USE mod_par_inv, only : KIDIM1,KIDIM2,KIDIM3,KPARS,KDIF 
      use mod_retr_settings_derived_type

      IMPLICIT NONE
! ---------------------------------------------------------------------
! IN
      type(retr_input_settings), intent(in)  :: RIN
      type(single_pixel_constraints_apriori), intent(in) :: SPCA
! ---------------------------------------------------------------------
! OUT
      INTEGER,                          INTENT(OUT) :: KSM
      REAL,DIMENSION(KPARS,KPARS),      INTENT(OUT) :: SMSING
      logical,                          intent(out) :: apriori_present
! ---------------------------------------------------------------------
! LOCAL 
      REAL,DIMENSION(KDIF,KDIF)   :: DIF,SM1
      INTEGER :: I,I1,J,J1,KM1,IDIM1,IDIM2,IDIM3
      LOGICAL :: IPRINT2
      
      INTEGER :: NDIM3,IO,ISTARSING
      REAL    :: GSM
      REAL    :: G, xnorm
! ---------------------------------------------------------------------

      IPRINT2 = .false.
      IF(IPRINT2) WRITE(*,*) "in smoothterm_single_pixel_apriori"
      KSM         = 0
      SMSING(:,:) = 0.0
      DIF(:,:) = 0.0
      xnorm = 0.0  ! in SMOOM smoothness matrix (SM1) is normalized

      DO IDIM1=1,RIN%NDIM%n1
      IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
      DO IDIM2=1,RIN%NDIM%n2(IDIM1)
         NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
         ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
         DO IDIM3=1,NDIM3
            IO  = SPCA%IO(IDIM3,IDIM2,IDIM1)
            GSM = SPCA%GSM(IDIM3,IDIM2,IDIM1)
            IF(GSM .GT. 0.) THEN
              G = 1.0
              DO I=1,NDIM3
                DIF(I,I) = 1.0
              ENDDO
              KM1 = NDIM3-IO
              KSM = KSM+KM1
              IF(IPRINT2) WRITE(*,'(2(a,i0))') "before SMOOM KM1 = ",KM1,"  NDIM3 = ",NDIM3
              CALL SMOOM( KM1,NDIM3,      & ! IN
                          DIF(:,:),       &
                          SM1(:,:), xnorm & ! OUT
                        )
              IF(IPRINT2) WRITE(*,*) "after SMOOM"
              DO J1=1,NDIM3
              DO I1=1,NDIM3
                SMSING(ISTARSING+I1-1,ISTARSING+J1-1) = &
                SMSING(ISTARSING+I1-1,ISTARSING+J1-1)+GSM*SM1(I1,J1)
!                WRITE(*,*) ISTARSING+I1-1,ISTARSING+J1-1,SM1(I1,J1),' I1,J1,SM1'
              ENDDO ! I1
              ENDDO ! J1
            ENDIF ! GSM .GT. 0.
         ENDDO ! IDIM3
      ENDDO ! IDIM2
      ENDDO ! IDIM1
	   
      apriori_present = .false.
      DO I=1,RIN%KNSINGF
         IF(SMSING(I,I) .LE. 0.0) CYCLE
         apriori_present = .true.
         EXIT
      ENDDO

      IF(IPRINT2) WRITE(*,*) 'before return from smoothterm_single_pixel_apriori'

      RETURN
      END SUBROUTINE smoothterm_single_pixel_apriori

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine determines smoothness term for single pixel (a priori smoothness)
        !>
        !> @param[in]    RIN          - invertion setting structure
        !> @param[in]    SPCS         -
        !> @param[out]   KSM          - 
        !> @param[out]   SMSING       - smoothness term of smoothness matrix for single pixel
        !> @param[out]   apriori_present - presence of a priori estimates
        !>

      SUBROUTINE smoothterm_single_pixel_smoothness ( RIN, SPCS,         & ! IN
                                                      KSM,SMSING,apriori_present  & ! OUT
                                                    )
!************************************************************
!*** Determining smoothness term for single pixel         ***
!************************************************************
      USE mod_par_inv, only : KIDIM1,KIDIM2,KIDIM3,KPARS,KDIF 
      use mod_retr_settings_derived_type
      use mod_stop_report

      IMPLICIT NONE
! ---------------------------------------------------------------------
! IN
      type(retr_input_settings), intent(in)  :: RIN
      type(single_pixel_constraints_smoothness), intent(in)  :: SPCS
! ---------------------------------------------------------------------
! OUT
      INTEGER,                          INTENT(OUT) :: KSM
      REAL,DIMENSION(KPARS,KPARS),      INTENT(OUT) :: SMSING
      logical,                          intent(out) :: apriori_present
! ---------------------------------------------------------------------
! LOCAL 
      real,dimension(:,:,:),allocatable :: ORDSING
      REAL,DIMENSION(KDIF,KDIF)   :: DIF,SM1
      REAL,DIMENSION(KDIF)        :: XNEW
      INTEGER :: I,I1,J,J1,KM1,IDIM1,IDIM2,IDIM3
      LOGICAL :: IPRINT2
      
      INTEGER :: NDIM3,IO,ISTARSING
      REAL    :: GSM
      REAL    :: G, xnorm
      integer :: alloc_stat
      logical :: IPRI_verbose
! ---------------------------------------------------------------------

      IPRI_verbose = RIN%IPRI_verbose

! ---------------------------------------------------------------------

!C*****************************************************************************************
!C*** ODRSING will be used for non-equal binning in smoothness:
!C*****************************************************************************************
      allocate(ORDSING(KIDIM3,KIDIM2,KIDIM1),stat=alloc_stat)
      if (alloc_stat /= 0) then
        write(tmp_message,'(a)') &
        ' error while trying to allocate ORDSING in smoothterm_single_pixel_smoothness'
        G_ERROR(trim(tmp_message))
      endif

      call ordsing_for_smooth_single_pixel(RIN,ORDSING)

      IPRINT2 = .false. ! .true. !
      IF(IPRINT2) WRITE(*,*) "in smoothterm_single_pixel_smoothness"
      KSM         = 0
      SMSING(:,:) = 0.0
      XNEW(:)     = 0.0
      xnorm = 0.0  ! in SMOOM smoothness matrix is normalized

      DO IDIM1=1,RIN%NDIM%n1
      IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
      DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
         NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
         ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
         IO        = SPCS%IO(IDIM2,IDIM1)
         GSM       = SPCS%GSM(IDIM2,IDIM1)
         if(IO .eq. 0 .and. GSM .gt. 0)  then
           write(*,*) IDIM1,IDIM2,IDIM3,IO,GSM,'  IDIM1,IDIM2,IDIM3,IO,GSM'
           write(*,*) 'For single pixel smoothness constrains: IO=0 and GSM>0 combination is not valid.'
           stop 'stop in smoothterm_single_pixel_smoothness'
         endif

         IF(GSM .GT. 0.) THEN
            DO IDIM3=1,NDIM3
            XNEW(IDIM3)=ORDSING(IDIM3,IDIM2,IDIM1)
            IF(IPRINT2) WRITE(*,*) IDIM3,XNEW(IDIM3) ,' IDIM3,XNEW(IDIM3)'
            ENDDO ! IDIM3
            IF(IPRINT2) THEN
            WRITE(*,*) "before DIFERW" 
            WRITE(*,*) IDIM3,IDIM2,IDIM1,IO,  &
	          '  IDIM3,IDIM2,IDIM1, NDIM3(..), IO(..)'
            ENDIF ! IPRINT2
			
            CALL DIFERW ( IPRI_verbose, &
                          NDIM3,IO,     & ! IN
                          XNEW(:),      &
                          G,DIF(:,:)    & ! OUT
                        )

            IF(IPRINT2) WRITE(*,*) "after DIFERW"
!OD?        IF (G.NE.1) AA=SQRT(0.05*0.05/(G*SPCS%GSM(I)))

            KM1=NDIM3-IO
            KSM=KSM+KM1
            IF(IPRINT2) WRITE(*,*) "before SMOOM"
			
            CALL SMOOM( KM1,NDIM3,      & ! IN
                        DIF(:,:),       &
                        SM1(:,:), xnorm & ! OUT
                      )
				 
            IF(IPRINT2) WRITE(*,*) "after SMOOM"
            DO I1=1,NDIM3
            DO J1=1,NDIM3
               SMSING(ISTARSING+I1-1,ISTARSING+J1-1)=               &
               SMSING(ISTARSING+I1-1,ISTARSING+J1-1)+GSM*SM1(I1,J1)
!CD            WRITE(*,*) NDIM%ISTARSING(IDIM2,IDIM1)+I1-1,   &
!CD            NDIM%ISTARSING(IDIM2,IDIM1)+J1-1,SM1(I1,J1),'   I,J'
            ENDDO ! J1
            ENDDO ! I1
		  ENDIF ! GSM .GT. 0.
      ENDDO ! IDIM2  
      ENDDO ! IDIM1
	   
      apriori_present = .false.
      DO I=1,RIN%KNSINGF
         IF(SMSING(I,I) .LE. 0.0) CYCLE
         apriori_present = .true.
         EXIT
      ENDDO

      deallocate(ORDSING,stat=alloc_stat)
      if (alloc_stat /= 0) stop &
      'error while trying to deallocate ORDSING in smoothterm_single_pixel_smoothness'

      IF(IPRINT2) WRITE(*,*) 'before return from smoothterm_single_pixel_smoothness'

      RETURN
      END SUBROUTINE smoothterm_single_pixel_smoothness

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine determines smoothness term for milti pixel inversion
        !> 
        !> @param[in]    iu_main_output - type of constrains 'smoothness' or 'apriori' estimates
        !> @param[in]    RIN            - invertion setting structure
        !> @param[in]    index_clouds   - contains indices for recurcive function CHECK  
        !> @param[out]   NISM           - number of groups in which smoothness can be implemented
        !> @param[out]   KSM            - number of elements in each of NISM groups 
        !> @param[out]   IMSM           - pixel numbers which belong to the group
        !> @param[out]   IKSIM          - number of equivalent measurements defining smoothness (used for residual)
        !> @param[out]   SMMULTI        - smoothness term of smoothness matrix for multi pixel
        !>

SUBROUTINE smoothterm_multi_pixel ( iu_main_output,RIN,     & ! IN
                                    index_clouds,           & 
                                    NISM,KSM,IMSM,          & ! OUT
                                    IKSIM,SMMULTI             &
                                  )

      use mod_par_inv, only : KIMAGE,KITIME,KIX,KIY,KPARS,KIMAGE,KMPSM,KDIF
      use mod_index_cloud
      use mod_retr_settings_derived_type
      
      IMPLICIT NONE

! ------------------------------------------------------------
! IN :
      integer,                             intent(in)  ::  iu_main_output
      type(ind_clouds),                    intent(in)  ::  index_clouds
      type(retr_input_settings),           intent(in)  ::  RIN	
! ------------------------------------------------------------
! OUT :
      INTEGER,                             INTENT(OUT) :: NISM,IKSIM
      INTEGER,DIMENSION(KMPSM),            INTENT(OUT) :: KSM   
      INTEGER,DIMENSION(KMPSM,KIMAGE),     INTENT(OUT) :: IMSM  
      REAL,DIMENSION(KIMAGE,KPARS,KIMAGE), INTENT(OUT) :: SMMULTI 

      !type(mp_smoothness),dimension(KIMAGE)                    :: mps
! ------------------------------------------------------------
! LOCAL :	  	  
      INTEGER                   :: NYSM,NXSM,NTSM 
      INTEGER,DIMENSION(KIX)    :: IXX
      INTEGER,DIMENSION(KIY)    :: IYY
      INTEGER,DIMENSION(KITIME) :: ITT

      REAL,DIMENSION(KDIF,KDIF) :: DIF,SM1
      REAL,DIMENSION(KDIF)      :: XNEW

      INTEGER :: ITIME,IX,IY,INX,INY,IDIM1,IDIM2,IDIM3,   &
	               I1,I2,J1,KM1,INT,IS,IS1,IS2
      LOGICAL :: IPRINT2, IPRI_verbose
      REAL    :: G, xnorm
      REAL    :: GSMX,GSMY,GSMT
      INTEGER :: NDIM3,ISTARSING,IOX,IOY,IOT
! ------------------------------------------------------------

      IPRI_verbose = RIN%IPRI_verbose

! ------------------------------------------------------------
! KMPSM - parameter, max number of groups for multi-pixel smoothness
! NISM  - number of groups in which smoothness can be implemented
! IKSIM - number of equivalent measurements defining smoothness 
!         (used for residual)   
! KSM (1:KMPSM)           - number of elements in each of NISM groups
! IMSM(1:KMPSM,1:KIMAGE)  - pixel numbers which belong to the group
! SMMULTI(1:KIMAGE,1:KPARS,1:KIMAGE)  - smoothness matrix
! ------------------------------------------------------------
!moved into read inversion setting input
      !if(KITIME .gt. KDIF) then 
	    !write(*,*) 'KITIME=',KITIME,' .gt. KDIF=',KDIF
	    !write(*,*) 'Check inversion parameters: KDIF should be max(KPARS,KITIME,KIX,KIY).'
	    !stop 'stop in smoothterm_multi_pixel'
      !endif
      !if(KIX .gt. KDIF) then 
	    !write(*,*) 'KIX=',KIX,' .gt. KDIF=',KDIF
	    !write(*,*) 'Check inversion parameters: KDIF should be max(KPARS,KITIME,KIX,KIY).'
	    !stop 'stop in smoothterm_multi_pixel'
      !endif
      !if(KIY .gt. KDIF) then 
	    !write(*,*) 'KIY=',KIY,' .gt. KDIF=',KDIF
	    !write(*,*) 'Check inversion parameters: KDIF should be max(KPARS,KITIME,KIX,KIY).'
	    !stop 'stop in smoothterm_multi_pixel'
      !endif	  

      xnorm = 0.0  ! smoothnes matrix is normalized
      XNEW(:)   = 0.0

      IPRINT2   = .false.
      NISM      = 0
      IKSIM     = 0
      KSM(:)    = 0
      IMSM(:,:) = 0
      SMMULTI(:,:,:)  = 0.0
      !SMIM2(:,:,:)  = 0.0
	  
      IF(IPRINT2) WRITE(*,*) "in smoothterm_multi_pixel" 

! ** ITIME, IX

      DO ITIME=1,index_clouds%NT
      DO IX=1,index_clouds%NX
        NYSM = 0

        DO IY=1,index_clouds%NY
        IF(index_clouds%ICLOUD(IX,IY,ITIME) .GT. 0) THEN
        NYSM = NYSM + 1
        XNEW(NYSM) = index_clouds%Y(IX,IY,ITIME)
        IYY(NYSM)  = index_clouds%INIMAGE(IX,IY,ITIME)
        IF(IPRINT2) WRITE(*,*) ITIME,IX,IY,NYSM,XNEW(NYSM),'  1: ITIME,IX,IY,NYSM,XNEW(NYSM)'
        ENDIF
        ENDDO ! IY

        IF(NYSM .GT. 1) THEN
          NISM=NISM+1
          KSM(NISM)=NYSM
          DO INY=1,NYSM
          IMSM(NISM,INY)=IYY(INY)
          ENDDO ! INY

          DO IDIM1=1,RIN%NDIM%n1
          IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
          DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
            !WRITE(*,'(7i5,/a)') ITIME,IX,IDIM2,IDIM1,'  ITIME,IX,IDIM2,IDIM1'
            NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
            ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
            IOY       = RIN%MPCS%IOY(IDIM2,IDIM1)
            GSMY      = RIN%MPCS%GSMY(IDIM2,IDIM1)
            IF(GSMY .GT. 0. .AND. IOY .GT. 0) THEN
            CALL DIFERW ( IPRI_verbose, &
                          NYSM,IOY,     & ! IN
                          XNEW(:),      &
                          G,DIF(:,:)    & ! OUT
                        )
            DO I2=1,NDIM3
              KM1 = NYSM-IOY
              IF(KM1 .LE. 0) KM1 = 1
              CALL SMOOM( KM1,NYSM,       & ! IN
                          DIF(:,:),       &
                          SM1(:,:), xnorm & ! OUT
                         )
              IKSIM=IKSIM+1
              DO I1=1,NYSM
              DO J1=1,NYSM             
              !SMIM(NISM,ISTARSING+I2-1,I1,J1)=                &
              !SMIM(NISM,ISTARSING+I2-1,I1,J1)+GSMY*SM1(I1,J1)
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) =  &
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) + GSMY*SM1(I1,J1)
              ENDDO ! J1
              ENDDO ! I1
            ENDDO ! I2
            ENDIF ! GSMY .GT. 0. .AND. IOY .GT. 0
          ENDDO ! IDIM2
          ENDDO ! IDIM1
        ENDIF ! NYSM .GT. 1

      ENDDO ! IX
      ENDDO ! ITIME				
      IF(IPRINT2) write(*,'(a,i5)') '1: nism=',nism 

! ** ITIME, IY
	   
      DO ITIME=1,index_clouds%NT
      DO IY=1,index_clouds%NY

        NXSM = 0
        DO IX=1,index_clouds%NX
        IF(index_clouds%ICLOUD(IX,IY,ITIME) .GT. 0) THEN
        NXSM=NXSM+1
        XNEW(NXSM)=index_clouds%X(IX,IY,ITIME)
        IXX(NXSM)=index_clouds%INIMAGE(IX,IY,ITIME)
        IF(IPRINT2) WRITE(*,*) ITIME,IY,IX,NXSM,XNEW(NXSM),'  2: ITIME,IY,IX,NXSM,XNEW(NXSM)'
        ENDIF
        ENDDO ! IX
         
        IF(NXSM .GT. 1) THEN
          NISM=NISM+1
          KSM(NISM)=NXSM
          DO INX=1,NXSM
          IMSM(NISM,INX)=IXX(INX)
          ENDDO ! INX
		 
          DO IDIM1=1,RIN%NDIM%n1
          IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
          DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
            !WRITE(*,'(7i5,/a)') ITIME,IY,IDIM2,IDIM1,'  ITIME,IY,IDIM2,IDIM1'
            NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
            ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
            IOX       = RIN%MPCS%IOX(IDIM2,IDIM1)
            GSMX      = RIN%MPCS%GSMX(IDIM2,IDIM1)
            IF(GSMX .GT. 0. .AND. IOX .GT. 0) THEN
            CALL DIFERW ( IPRI_verbose, &
                          NXSM,IOX,     & ! IN
                          XNEW(:),      &
                          G,DIF(:,:)    & ! OUT
                        )
            DO I2=1,NDIM3
              KM1 = NXSM-IOX
              IF(KM1 .LE. 0) KM1 = 1
              CALL SMOOM( KM1,NXSM,       & ! IN
                          DIF(:,:),       &
                          SM1(:,:), xnorm & ! OUT
                        )
              IKSIM=IKSIM+1
              DO I1=1,NXSM
              DO J1=1,NXSM
              !SMIM(NISM,ISTARSING+I2-1,I1,J1)=    &
              !SMIM(NISM,ISTARSING+I2-1,I1,J1)+GSMX*SM1(I1,J1)
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) =  &
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) + GSMX*SM1(I1,J1)		 
              ENDDO ! J1
              ENDDO ! I1
            ENDDO ! I2
            ENDIF ! GSMX .GT. 0. .AND. IOX .GT. 0
          ENDDO ! IDIM2
          ENDDO ! IDIM1
        ENDIF ! NXSM .GT. 1

      ENDDO ! IY
      ENDDO ! ITIME
      IF(IPRINT2) write(*,'(a,i5)') '2: nism=',nism

! ** IX, IY

      DO IX=1,index_clouds%NX
      DO IY=1,index_clouds%NY

        NTSM=0
        DO ITIME=1,index_clouds%NT
        IF(index_clouds%ICLOUD(IX,IY,ITIME) .GT. 0) THEN
        NTSM=NTSM+1
        !write(*,*) 'real(T)=',real(index_clouds%T(IX,IY,ITIME))
        XNEW(NTSM)=real(index_clouds%T(IX,IY,ITIME))
        ITT(NTSM)=index_clouds%INIMAGE(IX,IY,ITIME)
        IF(IPRINT2) WRITE(*,*) IX,IY,ITIME,NTSM,XNEW(NTSM),'  3: IX,IY,ITIME,NTSM,XNEW(NTSM)'
        ENDIF
        ENDDO ! ITIME

        IF(NTSM .GT. 1) THEN
          NISM=NISM+1
          KSM(NISM)=NTSM
          DO INT=1,NTSM
          IMSM(NISM,INT)=ITT(INT)
          ENDDO ! INT

          DO IDIM1=1,RIN%NDIM%n1
          IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
          DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
            !WRITE(*,'(7i5,/a)') IX,IY,IDIM2,IDIM1,'  IX,IY,IDIM2,IDIM1'
            NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
            ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
            IOT       = RIN%MPCS%IOT(IDIM2,IDIM1)
            GSMT      = RIN%MPCS%GSMT(IDIM2,IDIM1)
            IF(GSMT .GT. 0. .AND. IOT .GT. 0) THEN
            CALL DIFERW ( IPRI_verbose, &
                          NTSM,IOT,     & ! IN
                          XNEW(:),      &
                          G,DIF(:,:)    & ! OUT
                        )
            DO I2=1,NDIM3
              KM1 = NTSM-IOT
              IF(KM1 .LE. 0) KM1 = 1
              CALL SMOOM( KM1,NTSM,       & ! IN
                          DIF(:,:),       &
                          SM1(:,:), xnorm & ! OUT
                        )
              IKSIM=IKSIM+1
              DO I1=1,NTSM
              DO J1=1,NTSM        
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) =  &
              SMMULTI(IMSM(NISM,I1),ISTARSING+I2-1,IMSM(NISM,J1)) + GSMT*SM1(I1,J1) 
              ENDDO ! J1
              ENDDO ! I1
            ENDDO ! I2
            ENDIF ! GSMT .GT. 0. .AND. IOT .GT. 0
          ENDDO ! IDIM2
          ENDDO ! IDIM1
        ENDIF ! NTSM .GT. 1

      ENDDO ! IY
      ENDDO ! IX
		IF(IPRINT2) write(*,'(a,i5)') '3: nism=',nism		

      IF(NISM .GT. KMPSM) THEN
         WRITE(*,'(2(a,i5))') 'NISM =',NISM,' .GT. KMPSM =',KMPSM
         STOP 'STOP in SUBROUTINE smoothterm_multi_pixel'
      ENDIF
	  
      IF(IPRINT2) THEN
         WRITE(*,*) 'END of smoothterm_multi_pixel!!!'
         WRITE(*,'(a12,2i5)') 'NISM,IKSIM: ',NISM,IKSIM
         DO IS=1,NISM
         write(*,'(a12,2i5)') 'KSM(IS),IS: ',KSM(IS),IS
         DO IS1=1,KSM(IS)
         write(*,*) IMSM(IS,IS1),IS, IS1,' IMSM(IS,IS1),IS, IS1' 
         ENDDO ! IS1
         ENDDO ! IS
      ENDIF ! IPRINT2

!**************
! TL
		!type mp_smoothness

			!integer					             ::  ksm
			!integer,dimension(KIMAGE,KPARS)	 ::  imsm	
			!real,   dimension(KIMAGE,KPARS)  ::  smim	

		!end type mp_smoothness

      !mps(1:KIMAGE)%ksm =       
      !mps(1:KIMAGE)%imsm(1:KIMAGE,1:KPARS) = 0
      !do i1=1,KIMAGE
      !mps(i1)%smim(1:KIMAGE,1:KPARS) = 0.0
      !enddo ! i1

!goto 66
      !SMMULTI(:,:,:)  = 0.0
      !DO IS=1,NISM
      !DO I1=1,KNSINGF
      !DO IS1=1,KSM(IS)
      !DO IS2=1,KSM(IS)
         !IF(IS1 .NE. IS2) THEN
            !SMMULTI(IMSM(IS,IS1),I1,IMSM(IS,IS2)) = SMIM(IS,I1,IS1,IS2)
            !mps(IMSM(IS,IS2))%smim(IMSM(IS,IS1),I1) = SMMULTI(IMSM(IS,IS1),I1,IMSM(IS,IS2))
            !IF(IMSM(IS,IS1) .GT. KIMAGE .OR. IMSM(IS,IS2) .GT. KIMAGE)  &
            !STOP 'STOP in smoothterm_multi_pixel: SMMULTI index problem'	
            !write(*,*) IMSM(IS,IS1),I1,IMSM(IS,IS2),	' - IMSM(IS,IS1),I1,IMSM(IS,IS2)'	   
         !ENDIF ! IS1 .NE. IS2
      !ENDDO ! IS2
      !ENDDO ! IS1      
		!ENDDO ! I1
      !ENDDO ! IS
!66 continue
!**************
!stop 'stop: test in smoothterm_multi_pixel'

! NISM  - number of groups in which smoothness can be implemented
! IKSIM - number of equivalent measurements defining smoothness 
      if(RIN%IPRI_verbose) then
        write(iu_main_output,'(a)') 'in smoothterm_multi_pixel:'
        write(iu_main_output,'(4x,2(a,i0,2x))') 'NISM = ',NISM,'IKSIM = ',IKSIM
      endif

RETURN
END SUBROUTINE smoothterm_multi_pixel

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine calculates smoothness matrix
        !> 
        !> @param[in]    KM  - number of rows
        !> @param[in]    KN  - number of columns
        !> @param[in]    DIF - matrix of differences 
        !> @param[out]   SM1 - smoothness matrix
        !> @param[inout] xnorm - smoothness normalization coefficient
        !>

      SUBROUTINE SMOOM( KM, KN, DIF, SM1, xnorm )

!C**************************************************
!C  THIS SUBROUTINE CALCULATES "SMOOTHNESS MATRIX":
!C           SM1=(DIF)^T (DIF)
!C**************************************************
!C  INPUT:
!C        KM  I        - number of lines
!C        KN  I        - number of columns
!C        DIF R(KM,KN) - matrix of differences
!C        C   R(KM)    - vector for putting weights in
!C                       smoothness matrix
!C  OUTPUT:
!C        SM R(KN,KN)  - smoothness matrix
!C                       SM=(DIF)**T(C)(DIF)
!C        SM is normalized if input xnorm=0.0
!C           is not normalized if input xnorm=-999.0
!C        GF           - the parameter coming from
!C              DELTA(ordinate), this parameter can be
!C              can be used for indicating norm of
!C              of derivatives A assumed in Lagrange 
!C              parameter: 
!C               Gamma=Sigma(mes)**2/(GF*A**2)
!C                see Dubovik & King [2000]
!C***************************************************
      USE mod_par_inv, only : KDIF
	  
      IMPLICIT NONE
! IN :
      INTEGER, INTENT(IN) ::  KM, KN
      REAL, DIMENSION(KDIF,KDIF), INTENT(IN) :: DIF
! OUT :
      REAL, DIMENSION(KDIF,KDIF), INTENT(OUT) :: SM1
      REAL, INTENT(INOUT) :: xnorm
! LOCAL :
      REAL,DIMENSION(KN) :: C
      INTEGER :: I,I1
      LOGICAL :: lstop = .false.
!C***************************************************

      IF(KDIF .LT. KM) THEN
        WRITE(*,'(3(a,i0))') 'KDIF = ',KDIF,' .lt. KM = ',KM, &
        ' not valid combination; KDIF must be greater than KM'
         lstop = .true.
      ENDIF
      IF(KDIF .LT. KN) THEN
        WRITE(*,'(3(a,i0))') 'KDIF = ',KDIF,' .lt. KN = ',KN, &
        ' not valid combination; KDIF must be greater than KN'
         lstop = .true.
      ENDIF
      if(lstop) &
      stop 'stop in SMOOM'

      C(1:KN)        = 0.
      SM1(1:KN,1:KN) = 0.
!do i=1,2
!write(*,*) 'in SMOOM: DIF',DIF(1:2,I)
!enddo

      IF(KN .GE. 1) THEN  ! CD 06/09/2013
        IF(C(1) .EQ. 0) THEN
          DO I =1,KN
          DO I1=1,KN
            SM1(I,I1) = SUM(DIF(1:KM,I)*DIF(1:KM,I1)) 
            !write(*,*) i,i1,SM1(I,I1),C(1),'  - i,i1,SM1(I,I1),C(1)'
          ENDDO ! I1
          ENDDO ! I
        ELSE
          DO I =1,KN
          DO I1=1,KN
            SM1(I,I1) = SUM(DIF(1:KM,I)*DIF(1:KM,I1)*C(1:KM))
            !write(*,*) i,i1,SM1(I,I1),C(1),'  - i,i1,SM1(I,I1),C(1)'
          ENDDO ! I1
          ENDDO ! I
        ENDIF ! C(1).EQ.0

!OD 2015-09-15 Smoothness matrix normalization
! For edges normalization is also used in subroutines DIF_before_calcul
! and DIF_after_calcul
        IF( xnorm .NE. -999.0 ) THEN
          xnorm = 0.
          DO I=1,KN
            IF(xnorm .LT. SM1(I,I)) xnorm = SM1(I,I)
          ENDDO ! I

          IF(xnorm .NE. 0.) THEN
            DO I =1,KN
            DO I1=1,KN
              SM1(I,I1) = SM1(I,I1)/xnorm
            ENDDO ! I1
!write(*,*) 'in SMOOM: I=',I,'  SM1=',SM1(I,1:KN)
            ENDDO ! I
!write(*,*) 'in SMOOM: xnorm=',xnorm
          ENDIF ! xnorm .NE. 0.
!write(*,*) 'in SMOOM: xnorm=',xnorm
        ENDIF ! xnorm .NE. -999.0
      ENDIF ! KN.GT.1

      RETURN
      END SUBROUTINE SMOOM  

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine defines matrix of n-th differences
        !> 
        !> @param[in]    KN  - number of parameters
        !> @param[in]    X  - values of ordinates Xi of the grid point
        !> @param[inout] IO  - order of differences
        !> @param[out]   G  - coefficient relating differences with derivatives
        !> @param[out]   DIF - matrix of differences 
        !>

      SUBROUTINE DIFERW(IPRI_verbose,KN,IO,X,G,DIF)
!C**********************************************************
!C*** THIS SUBROUTINE DEFINES MATRIX of N-th DIFFERENCES ***
!C*** for the function Y(X) defined in the points Xi     ***
!C    INPUT:
!C        KN     I   - number of parameters
!C        IO     I   - order of differences
!C        X(KN)  R   - the values of ordinates Xi of the ***
!C                     grid point
!C
!C***
!C    OUTPUT:
!C      G    R   - coefficient relating differences with
!C                 derivatives
!C            - G.NE.1 if Xi+1-Xi is non constant 
!C            - G=1/(A**IO) if Xi+1-Xi=const=A
!C              this parameter can be used for indicating norm
!C              of derivatives A assumed in Lagrange parameter: 
!C               Gamma=Sigma(mes)**2/(GF*A**2)
!C                see Dubovik & King [2000]
!C      DIF  R(KN-IO,KN) -matrix of IO-th differences     
!C****                                                   ***
!C**********************************************************
      USE mod_par_inv, only : KDIF
	  
      IMPLICIT NONE
! ----------------------------------------------------------
! IN :
      LOGICAL,             INTENT(IN)       ::  IPRI_verbose
      INTEGER,             INTENT(IN)       ::  KN
      INTEGER,             INTENT(INOUT)    ::  IO ! do not change INOUT
      REAL,DIMENSION(KDIF),INTENT(IN)       ::  X	
! ----------------------------------------------------------
! OUT :
      REAL,                     INTENT(OUT)  ::  G
      REAL,DIMENSION(KDIF,KDIF),INTENT(OUT)  ::  DIF	  
    
! ----------------------------------------------------------
! LOCAL :
      REAL,DIMENSION(KN)         ::  DX,DX1
      REAL,DIMENSION(KN,KN)      ::  DIFT      
      INTEGER                    ::  I,J,II,IIO
      LOGICAL                    ::  diff_const	  
! ----------------------------------------------------------
!      write(*,*) "AH KDIF = ",KDIF,"  size(X) = ",size(X)
      G = 0.0
      DIFT(:,:) = 0.  ! CD 06/09/2013
      DIF (:,:) = 0.  ! CD 06/09/2013 
      
      IF(IO .GE. KN) THEN 
         IF(IPRI_verbose .EQV. .TRUE.) &
         WRITE(*,*) IO,KN,' IO,KN, IO.GE.KN - PROBLEM with SMOOTHNESS!!!'
         IO=KN-1
         IF(IPRI_verbose .EQV. .TRUE.) &
         WRITE(*,*) IO,KN,' IO,KN  - IO has been changed    !!!'
         IF(IO .eq. 0) RETURN 
      ENDIF ! IO
    
      IF(KN .GE. 1) THEN ! CD 06/09/2013
         DO J=1,KN
            DIFT(J,J) = 1.0
            DX(J)     = X(J)
         ENDDO ! J

!write(*,*) 'in DIFERW IO=',IO,'  KN=',KN
         IF(IO .GT. 0) THEN
!C**********************************************************
!C*** Checking if Xi+1-Xi=const or not                   ***
			   diff_const = .true.
         IF(KN .GT. 2) THEN
            DO J=2,KN-1
            IF(ABS(((X(J-1)-X(J))-(X(J)-X(J+1)))/(X(J-1)-X(J))) .GT. 1.0e-3) THEN
            diff_const = .false.
            EXIT 
            ENDIF
            ENDDO ! J
         ELSEIF(KN .GT. 1) THEN
!write(*,*) 'in DIFERW X:',X(1:2)
            IF(ABS((X(2)-X(1))/X(1)) .GT. 1.0e-3) THEN
            diff_const = .false. 
            ENDIF    
         ENDIF ! KN.GT.2
         IF(diff_const) THEN
!C**********************************************************
!C*** Calculating matrix of differences and G            ***
!C***   for  Xi+1-Xi=const                               ***
               G = 1.0/(ABS((X(1)-X(2)))**IO)
               IF(G .NE. G) THEN
               WRITE(*,*) 'G=',G,' !!!  (X(1)-X(2))=',(X(1)-X(2)),' in DIFERW'
               STOP 'STOP IN DIFERW'
               ENDIF
!CD            WRITE(6,*) KN,IO,' KN,IO'
               DO IIO=1,IO 
               DO I=1,KN-IIO
               DO J=1,KN
               DIF(I,J)  = DIFT(I,J)-DIFT(I+1,J)
               ENDDO ! J
               ENDDO ! I
               
               DO J=1,KN
               DO I=1,KN-IIO
               DIFT(I,J) = DIF(I,J)
               ENDDO ! I
               ENDDO ! J
               ENDDO ! IIO

	        ELSE ! diff_const
!C**********************************************************
!C*** Calculating matrix of differences and G            ***
!C***   for  Xi+1-Xi not const                           ***
               G = 1.0
!cs            write(6,*)'check5'
               DO I=1,KN
               DX(I)  = X(I)
               DX1(I) = X(I)
               ENDDO ! I
               DO IIO=1,IO
               DO I=1,KN-IIO
               DO J=1,KN
               IF((DX(I)-DX(I+1)) .EQ. 0.0) THEN
               WRITE(*,*) DX(I),DX(I+1),' DX(I),DX(I+1),DX=0 !!!'
               STOP 'STOP IN DIFERW'
               ENDIF
               !DIF(I,J) = 1.0/ABS((DX(I)-DX(I+1)))*(DIFT(I,J)-DIFT(I+1,J))
               DIF(I,J) = (DIFT(I,J)-DIFT(I+1,J))/ABS(DX(I)-DX(I+1))
               ENDDO ! J
               ENDDO ! I
                
               DO I=1,KN-IIO
               DO J=1,KN
               DIFT(I,J)=DIF(I,J)
               ENDDO ! J
               ENDDO ! I
               DO I=1,KN-IIO
               DX(I)  = (DX1(I)+DX1(I+1))*0.5
               ENDDO ! I
               DO I=1,KN-IIO
               DX1(I) = DX(I)
               ENDDO ! I
               ENDDO ! IIO
	        ENDIF ! diff_const

! OD 2015-07-20
!          DIF(:,:) = 0.0
!          DO I=1,KN-IO
!          DO J=1,KN
!            DIF(I,J) = DIFT(I,J)
!          ENDDO ! J
!          ENDDO ! I

         ELSE ! IO 
!C*********For a priori estimates *************************
            DO J=1,KN
            DIF(J,J) = 1.0  !CD DIF(J,J)=0.0 !18-01-2013
            ENDDO ! J
            G = 1.0
         ENDIF ! IO .GT. 0
      ENDIF ! KN .GT. 1

      RETURN
      END SUBROUTINE DIFERW
      
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Recursive procedure searches pixels to form time-space groups
      RECURSIVE SUBROUTINE CHECK(IT,IX,IY,ITS,IXS,IYS,ISEARCH,  &
      NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,               &
      T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS)  

      USE mod_par_inv, only : KITIME,KIX,KIY,KIMAGE

      IMPLICIT NONE 
! --------------------------------------------------------------	  
! IN:
      INTEGER,DIMENSION(KIX,KIY,KITIME), INTENT(IN) :: ICLOUD
      REAL,DIMENSION(KIX,KIY,KITIME),    INTENT(IN) :: T,X,Y                    
      REAL,                              INTENT(IN) :: TTEST,XTEST,YTEST

! --------------------------------------------------------------
! INOUT :
      INTEGER,DIMENSION(KIX,KIY,KITIME) :: INDS
      INTEGER,DIMENSION(KIMAGE,KIMAGE)  :: ITP,IXP,IYP
      INTEGER,DIMENSION(KIMAGE)         :: NIPN

      INTEGER,                        INTENT(IN)    :: NT,NX,NY
      INTEGER,                        INTENT(INOUT) :: IT,IX,IY,     &
						                               ITS,IXS,IYS,  &
						                               ISEARCH,IP1
! --------------------------------------------------------------                                                      
! LOCAL :
      INTEGER ::  IT1,IX1,IY1
      LOGICAL ::  IPRINT2
      REAL    ::  GT,GX,GY,    &
	              DT,DX,DY  
! --------------------------------------------------------------
	       
       IPRINT2 = .false.
       IF(IPRINT2) THEN
        WRITE(*,*) IP1, "IP1 - IN CHECK"
        write(*,*) IT, IX, IY,'IT, IX, IY'
        WRITE(*,*) INDS(IX,IY,IT),'INDS(IX,IY,IT)'        
        WRITE(*,*) ISEARCH,'ISEARCH BEGINNING'
       ENDIF
       IF(INDS(IX,IY,IT).EQ.0) THEN
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
        INDS(IX,IY,IT)=IP1
       ENDIF !IF(INDS(IX,IY,IT).EQ.0)
       IT1=IT+1
    1  CONTINUE
        IF(IT1 .GE. 1 .AND. IT1 .LE. NT) THEN
         IF(INDS(IX,IY,IT1) .EQ. 0) THEN
         IF(ICLOUD(IX,IY,IT1) .EQ. 0) THEN
         IT1=IT1+1
         GOTO 1
         ENDIF
       DT=T(IX,IY,IT1)-T(IX,IY,IT)
       DT=ABS(DT)
       IF(IPRINT2) &
       WRITE(*,*) IT,TTEST,GT/DT,GT,DT,'  - IT,TTEST,GT/DT,GT,DT'
        IF(DT .EQ. 0 .OR. TTEST .LE. GT/DT) THEN
        IF(IPRINT2) WRITE(*,*) 'AFTER TEST 1'
        NIPN(IP1)=NIPN(IP1)+1
        IT=IT1
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH .EQ. 1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH = 0
         ENDIF
!CD        INDS(IX,IY,IT)=IP1     
!CD        write(*,*) NIPN(IP1),'NIPN(IP1), before 1'
        CALL CHECK (                                             &
                    IT,IX,IY,ITS,IXS,IYS,ISEARCH,                &
                    NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,  &
                    T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS              &
                   )
        ELSE 
         IF(ISEARCH .EQ. 0) THEN
         ISEARCH=1
         ITS=IT1
         IXS=IX
         IYS=IY
         IF(IPRINT2) WRITE(*,*) ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF !IF(INDS(IX,IY,IT1).EQ.0)  
          IT1=IT-1

   2   CONTINUE
        IF(IT1 .GE. 1 .AND. IT1 .LE. NT) THEN
         IF(INDS(IX,IY,IT1) .EQ. 0) THEN
         IF(ICLOUD(IX,IY,IT1) .EQ. 0) THEN
         IT1=IT1-1
         GOTO 2
         ENDIF
       DT=T(IX,IY,IT1)-T(IX,IY,IT)
       DT=ABS(DT)
       IF(IPRINT2) &
       WRITE(*,*) IT,TTEST,GT/DT,GT,DT,'  - IT,TTEST,GT/DT,GT,DT'
        IF(DT .EQ. 0 .OR. TTEST .LE. GT/DT) THEN
        IF(IPRINT2)        WRITE(*,*) 'AFTER TEST 2'
        NIPN(IP1)=NIPN(IP1)+1
        IT=IT1
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH.EQ.1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH=0
         ENDIF
!CD        INDS(IX,IY,IT)=IP1
!CD        write(*,*) IT, IX, IY,'IT, IX, IY'
!CD        write(*,*) NIPN(IP1), 'NIPN(IP1), before 2'
         CALL CHECK (                                             &
		             IT,IX,IY,ITS,IXS,IYS,ISEARCH,                &
                     NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,  &
	                 T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS              &
					)  
        ELSE 
         IF(ISEARCH .EQ. 0) THEN
         ISEARCH=1
         ITS=IT1
         IXS=IX
         IYS=IY
       IF(IPRINT2) WRITE(*,*)ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF ! IF(INDS(IX,IY,IT1).EQ.0) 
       IX1=IX+1

    3  CONTINUE
         IF(IX1 .GE. 1 .AND. IX1 .LE. NX) THEN
         IF(INDS(IX1,IY,IT) .EQ. 0) THEN
         IF(ICLOUD(IX1,IY,IT) .EQ. 0) THEN
         IX1=IX1+1
         GOTO 3
         ENDIF
       DX=X(IX1,IY,IT)-X(IX,IY,IT)
       DX=ABS(DX)
       IF(IPRINT2) WRITE(*,*) IX,XTEST,GX/DX, GX, DX,  &
      'IX,XTEST,GX/DX, GX, DX,'
        IF(DX .EQ. 0 .OR. XTEST .LE. GX/DX) THEN
        IF(IPRINT2) WRITE(*,*) 'AFTER TEST 3'
        NIPN(IP1)=NIPN(IP1)+1
        IX=IX1
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH.EQ.1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH = 0
         ENDIF
!CD        INDS(IX,IY,IT)=IP1
!CD       WRITE(*,*) NIPN(IP1),IP1,'NIPN(IP1),IP1'
!CD        write(*,*) 'before 3'
!CD        WRITE(*,*) NIPN(IP1),IP1,'NIPN(IP1),IP1'
         CALL CHECK (                                              &
                     IT,IX,IY,ITS,IXS,IYS,ISEARCH,                 &
                     NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,   &
                     T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS               &
                    )
        ELSE 
         IF(ISEARCH .EQ. 0) THEN
         ISEARCH=1
         ITS=IT
         IXS=IX1
         IYS=IY
      IF(IPRINT2) WRITE(*,*)ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF ! IF(INDS(IX1,IY,IT).EQ.0)
       IX1=IX-1

    4  CONTINUE
         IF(IX1 .GE. 1 .AND. IX1 .LE. NX) THEN
         IF(INDS(IX1,IY,IT) .EQ. 0) THEN
         IF(ICLOUD(IX1,IY,IT) .EQ. 0) THEN
         IX1=IX1-1
         GOTO 4
         ENDIF
       DX=X(IX1,IY,IT)-X(IX,IY,IT)
       DX=ABS(DX)
      IF(IPRINT2) WRITE(*,*) IX,XTEST,GX/DX, GX, DX, &
      'IX,XTEST,GX/DX, GX, DX,'
        IF(DX .EQ. 0 .OR. XTEST .LE. GX/DX) THEN
      IF(IPRINT2)   WRITE(*,*) 'AFTER TEST 4'
        NIPN(IP1)=NIPN(IP1)+1
        IX=IX1
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH.EQ.1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH = 0
         ENDIF
!CD        INDS(IX,IY,IT)=IP1
!CD        write(*,*) 'before 4'
!CD        WRITE(*,*) NIPN(IP1),IP1,'NIPN(IP1),IP1'
         CALL CHECK (                                               &
                     IT,IX,IY,ITS,IXS,IYS,ISEARCH,                  &
                     NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,    &
                     T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS                &
					          )
        ELSE 
         IF(ISEARCH .EQ. 0) THEN
         ISEARCH=1
         ITS=IT
         IXS=IX1
         IYS=IY
      IF(IPRINT2) WRITE(*,*)ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF !IF(INDS(IX1,IY,IT).EQ.0)
       IY1=IY+1

    5  CONTINUE
         IF(IY1 .GE. 1 .AND. IY1 .LE. NY) THEN
!CD         WRITE(*,*) IT,IX,IY1,'IT,IX,IY1'
!CD         WRITE(*,*) INDS(IX,IY1,IT),'INDS(IX,IY1,IT)'
         IF(INDS(IX,IY1,IT) .EQ. 0) THEN
         IF(ICLOUD(IX,IY1,IT) .EQ. 0) THEN
         IY1=IY1+1
         GOTO 5
         ENDIF
       DY=ABS(Y(IX,IY1,IT)-Y(IX,IY,IT))
       DY=ABS(DY)
       IF(IPRINT2) WRITE(*,*) IY,YTEST,GY/DY,GY,DY,  &
       'IY,YTEST,GY/DY,GY,DY'
        IF(DY .EQ. 0 .OR. YTEST .LE. GY/DY) THEN
       IF(IPRINT2) WRITE(*,*) 'after test 5'
        NIPN(IP1)=NIPN(IP1)+1
        IY=IY1
       IXP(IP1,NIPN(IP1))=IX
       IYP(IP1,NIPN(IP1))=IY
       ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH.EQ.1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH = 0
         ENDIF
!CD       INDS(IX,IY,IT)=IP1
!CD       write(*,*) 'before 5'
!CD       WRITE(*,*) NIPN(IP1),IP1,'NIPN(IP1),IP1'
         CALL CHECK (                                               &
		                 IT,IX,IY,ITS,IXS,IYS,ISEARCH,                  &
                     NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,    &
                     T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS                &
					          )
        ELSE 
         IF(ISEARCH .EQ. 0) THEN
         ISEARCH=1
         ITS=IT
         IXS=IX
         IYS=IY1
      IF(IPRINT2) WRITE(*,*) ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF ! IF(INDS(IX,IY1,IT).EQ.0) THEN
       IY1=IY-1

    6  CONTINUE
         IF(IY1 .GE. 1 .AND. IY1 .LE. NY) THEN
         IF(INDS(IX,IY1,IT) .EQ. 0) THEN
         IF(ICLOUD(IX,IY1,IT) .EQ. 0) THEN
         IY1=IY1-1
         GOTO 6
         ENDIF
       DY=ABS(Y(IX,IY1,IT)-Y(IX,IY,IT))
       DY=ABS(DY)
       IF(IPRINT2) WRITE(*,*) IY,YTEST,GY/DY,GY,DY,  &
       'IY,YTEST,GY/DY,GY,DY'
        IF(DY .EQ. 0 .OR. YTEST .LE. GY/DY) THEN
        IF(IPRINT2) WRITE(*,*) 'after test 6'
        NIPN(IP1)=NIPN(IP1)+1
        IY=IY1
        IXP(IP1,NIPN(IP1))=IX
        IYP(IP1,NIPN(IP1))=IY
        ITP(IP1,NIPN(IP1))=IT
         IF(ISEARCH .EQ. 1) THEN
       IF(ITS .EQ. IT .AND. IXS .EQ. IX .AND. IYS .EQ. IY) ISEARCH = 0
         ENDIF
!CD        INDS(IX,IY,IT)=IP1
!CD        write(*,*) 'before 6'
!CD        WRITE(*,*) NIPN(IP1),IP1,'NIPN(IP1),IP1'
         CALL CHECK (                                               &
                     IT,IX,IY,ITS,IXS,IYS,ISEARCH,                  &
                     NT,NX,NY,GT,GX,GY,TTEST,XTEST,YTEST,ICLOUD,    &
                     T,X,Y,IP1,NIPN,ITP,IXP,IYP,INDS                &
					) 
        ELSE 
!CD        WRITE(*,*) ISEARCH,'ISEARCH after 6'
         IF(ISEARCH .EQ.  0) THEN
         ISEARCH=1
         ITS=IT
         IXS=IX
         IYS=IY1
         IF(IPRINT2) WRITE(*,*) ITS,IXS,IYS,'  ITS,IXS,IYS IN'
         ENDIF
        ENDIF 
       ENDIF
       ENDIF ! IF(INDS(IX,IY1,IT).EQ.0) THEN
!CD       WRITE(*,*) ISEARCH,ITS,IXS,IYS, 'ISEARCH END, ITS,IXS,IYS'
       RETURN
       END SUBROUTINE CHECK
	   
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine selects Time and Space qroups of the same parameters in order to reduce Jacobian size. 
        !> 
        !> @param[in]      iu_main_output - unit number for main output file
        !> @param[in]      RIN            - invertion setting structure
        !> @param[in]      index_clouds   - contains indices for recurcive function CHECK  
        !> @param[in]      npixels        - number of inverted pixels
        !!> @param[out]    IP2            - 
        !!> @param[inout]  NIPPN          - 
        !!> @param[inout]  IPPN           -                
        !>

      subroutine time_space_group ( iu_main_output,     & ! IN
                                    RIN,                & 
                                    index_clouds,       &
                                    npixels,            &
                                    IP2,NIPPN,IPPN      & ! OUT           
                                  )

      use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_par_inv, only : KIMAGE,KITIME,KIX,KIY,KPAR

      implicit none
!-------------------------------------------------------------
! IN :
      type(retr_input_settings), intent(in) :: RIN
      type(ind_clouds),          intent(in) :: index_clouds
      integer,                   intent(in) :: npixels,iu_main_output
!-------------------------------------------------------------
! OUT :
      integer,                       intent(out)   :: IP2
      integer,dimension(KPAR),       intent(inout) :: NIPPN                         
      integer,dimension(KPAR,KPAR),  intent(inout) :: IPPN

!-------------------------------------------------------------
! LOCAL :
      integer  :: IPS,IT,IX,IY,IDIM1,IDIM2,IDIM3,           &
                  IX0,IY0,IP1,ICHECK,IPCHECK,ISEARCH,       &
                  IX1,IY1,IP,I1,IT1,I,IIA,IIIA,I3,IP3
      integer  :: ITS,IXS,IYS,IB,IA
      real     :: GT,GX,GY,TTEST,XTEST,YTEST
      integer, dimension(KIMAGE) :: NIPN
      integer, dimension(KIMAGE,KIMAGE) :: IXP,IYP,ITP,IXPA,IYPA,ITPA 
      integer :: IT0=0
      INTEGER,DIMENSION(KPAR) :: IPPS
      INTEGER,DIMENSION(KIX,KIY,KITIME) :: INDS
!-------------------------------------------------------------
!**************
! OD ! do not delete
		!integer                                         :: ipix,IMAGE_IP
		!integer, dimension(KIMAGE)                      :: NIP_IM
		!integer, dimension(KIDIM1*KIDIM2,KIMAGE)        :: ISTART,ISTOP,NIMAGE_IP
		!integer, dimension(KIMAGE,KIDIM1*KIDIM2,KIMAGE) :: IMAGE_I

!-------------------------------------------------------------

		!NIP_IM(:)=0
		!ISTART(:,:)=0
		!ISTOP(:,:)=0
		!NIMAGE_IP(:,:)=0
		!IMAGE_I(:,:,:)=0
! OD
!**************		

!C*** Selection of the group of the same parameters

      TTEST=RIN%IPFP%TTEST
      XTEST=RIN%IPFP%XTEST
      YTEST=RIN%IPFP%YTEST

      IF(XTEST .GT. 0 .OR. YTEST .GT. 0 .OR. TTEST .GT. 0) THEN
      YTEST=YTEST/0.05
      XTEST=XTEST/0.05
      IP2=0
      IPS=0
      DO IDIM1=1,RIN%NDIM%n1
        IF(.not. RIN%NDIM%par_retr(IDIM1)) EXIT
        DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
          GT=RIN%MPCS%GSMT(IDIM2,IDIM1)
          GX=RIN%MPCS%GSMX(IDIM2,IDIM1)
          GY=RIN%MPCS%GSMY(IDIM2,IDIM1)
          DO IX=1,index_clouds%NX
          DO IY=1,index_clouds%NY 
          DO IT=1,index_clouds%NT
          INDS(IX,IY,IT)=0
          IF(IT0 .EQ. 0) THEN
          IF(index_clouds%ICLOUD(IX,IY,IT).EQ.1) THEN
          IT0=IT
          IX0=IX
          IY0=IY
          ENDIF ! ICLOUD
          ENDIF ! IT0
          ENDDO ! IT
          ENDDO ! IY
          ENDDO ! IX
          IP1=1
          NIPN(1)=1
          IT=IT0
          IX=IX0
          IY=IY0
          IPCHECK=0
          ISEARCH=0

! Recursive procedure 
		 
17        CALL CHECK(IT,IX,IY,ITS,IXS,IYS,ISEARCH,            &
          index_clouds%NT,index_clouds%NX,index_clouds%NY,    &
					GT,GX,GY,TTEST,XTEST,YTEST,index_clouds%ICLOUD,     &
          real(index_clouds%T),index_clouds%X,index_clouds%Y, & 
					IP1,NIPN,ITP,IXP,IYP,INDS)
					 
!CD         WRITE(iu_main_output,*) 'after search'
!CD         WRITE(iu_main_output,*) IP1,NIPN(IP1),'IP1,NIPN(IP1)'
!CD         IF(NIPN(IP1).EQ.1) THEN
!CD           WRITE(iu_main_output,*) IXP(IP1,NIPN(IP1)),IYP(IP1,NIPN(IP1)), &
!CD           ITP(IP1,NIPN(IP1)),'IXP,IYP,ITP'
!CD         ENDIF
          IF(ISEARCH .EQ. 1) THEN
          IPCHECK=0
          ISEARCH=0
          IT=ITS
          IX=IXS
          IY=IYS
          IP1=IP1+1
          NIPN(IP1)=1
          GO TO 17
          ENDIF ! ISEARCH

          DO IY1=1,index_clouds%NY
          DO IX1=1,index_clouds%NX
          DO IT1=1,index_clouds%NT
               IF(index_clouds%INIMAGE(IX1,IY1,IT1).GE.1) THEN
               IF(INDS(IX1,IY1,IT1).EQ.0) THEN
               IPCHECK=IPCHECK+1
               IF((NIPN(IP1)-IPCHECK).GE.1) THEN
               IT=ITP(IP1,NIPN(IP1)-IPCHECK)
               IX=IXP(IP1,NIPN(IP1)-IPCHECK)
               IY=IYP(IP1,NIPN(IP1)-IPCHECK)
               ELSE
               IT=IT1
               IX=IX1
               IY=IY1
               IPCHECK=0
               IP1=IP1+1
               NIPN(IP1)=1
               ENDIF ! (NIPN(IP1)-IPCHECK).GE.1
               GO TO 17
               ENDIF ! INDS(IX1,IY1,IT1).EQ.0
               ENDIF ! INIMAGE(IX1,IY1,IT1).GE.1
          ENDDO ! IT1
          ENDDO ! IX1
          ENDDO ! IY1

!        K_GROUP(IDIM2,IDIM2)=NIPN(IP1)
          DO IP=1, IP1
          IF(NIPN(IP) .GT. 1) THEN
          DO I1=1,NIPN(IP)
          IA=npixels
          IF (I1.EQ.1) IB=0
          DO I=1,NIPN(IP)
          IIIA=index_clouds%INIMAGE(IXP(IP,I),IYP(IP,I),ITP(IP,I))
          IF(IIIA .LE. IA .AND. IIIA .GT. IB) THEN
          IA=IIIA
          IIA=I
          ENDIF ! IIIA
          ENDDO ! I
          ITPA(IP,I1) = ITP(IP,IIA)
          IXPA(IP,I1) = IXP(IP,IIA)
          IYPA(IP,I1) = IYP(IP,IIA)
          IB=IA
          ENDDO ! I1
          DO I=1,NIPN(IP)
          ITP(IP,I)=ITPA(IP,I)
          IXP(IP,I)=IXPA(IP,I)
          IYP(IP,I)=IYPA(IP,I)
          ENDDO ! 
          ENDIF ! NIPN(IP).GT.1
!           N_GROUP(IP,IDIM2,IDIM1)=NIPN(IP)
!           DO I=1,NIPN(IP)
!             IM_GROUP(IP,I,IDIM2,IDIM1)=           &
!             INIMAGE(ITP(IP,I),IXP(IP,I),IYP(IP,I))
!           ENDDO ! I

           WRITE(iu_main_output,'(a3,i4,3x,a8,2x,1000i5)')      &
           'IP=',IP,' INIMAG:',                                    &
           (index_clouds%INIMAGE(IXP(IP,I),IYP(IP,I),ITP(IP,I)),   &
		                                           I=1,NIPN(IP))
!**************
! OD ! do not delete
	    !DO ipix=1,npixels
         !IF(NIPN(IP).GT.1.AND.   &
		 !index_clouds%INIMAGE(IXP(IP,1),IYP(IP,1),ITP(IP,1)).EQ.IIMAGE) THEN
		   !NIP_IM(ipix)=NIP_IM(ipix)+1
		   !ISTART(NIP_IM(ipix),ipix)=RIN%NDIM%ISTARSING(IDIM2,IDIM1)
		   !ISTOP(NIP_IM(ipix),ipix)=RIN%NDIM%ISTARSING(IDIM2,IDIM1)+  &
		   !                             RIN%NDIM%n3(IDIM2,IDIM1)-1
		   !NIMAGE_IP(NIP_IM(ipix),ipix)=NIPN(IP)
		   !DO I=1,NIPN(IP)
		   !  IMAGE_I(ipix,NIP_IM(ipix),I)=  &
		   !          index_clouds%INIMAGE(IXP(IP,I),IYP(IP,I),ITP(IP,I))
		   !ENDDO
		 !ENDIF
		!ENDDO ! ipix
! OD
!**************
		   
      ENDDO ! IP

      DO IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
		  DO IP=1,IP1
			IF(NIPN(IP) .GT. 1) THEN
			IP2=IP2+1
			I1=0
			DO I=1,NIPN(IP)
			  I1=I1+1
        IPPN(IP2,I1)=RIN%KNSING *                                    &
        (index_clouds%INIMAGE(IXP(IP,I),IYP(IP,I),ITP(IP,I))-1) +    &
        RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1
			ENDDO ! I
			NIPPN(IP2)=I1
			ELSE
			IPS=IPS+1
			IPPS(IPS)=RIN%KNSING *                                          &
            (index_clouds%INIMAGE(IXP(IP,1),IYP(IP,1),ITP(IP,1))-1) + &
            RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1
			ENDIF ! NIPN(IP).GT.1
		  ENDDO ! IP
      ENDDO ! IDIM3

      ENDDO ! IDIM2
      ENDDO ! IDIM1

      IF(IP2 .GT. 0) THEN

      IF(IPS .GT. 0) THEN
        NIPPN(IP2+1)=IPS
        DO I=1,IPS
        IPPN(IP2+1,I)=IPPS(I)
        ENDDO ! I
      ENDIF ! IPS
      ICHECK=0
      DO IP=1,IP2+1
        ICHECK=ICHECK+NIPPN(IP)
      ENDDO ! IP
      IF(RIN%IPRI_additional_info) THEN
      WRITE(iu_main_output,*) IP2,'  IP2 - FINAL'
      DO IP3=1,IP2 
        WRITE(iu_main_output,*)                &
        'IP3=',IP3,'  NIPPN(IP3)=',NIPPN(IP3)     &
        ,'  IPPN(IP3,I3):',                       &
        (IPPN(IP3,I3),I3=1,NIPPN(IP3)),'  IPPN(IP2,I3)'
      ENDDO ! IP3
      ENDIF ! RIN%IPRI_additional_info
	   
      IF(NIPPN(IP2+1) .GT. 0) THEN
        WRITE(iu_main_output,*) '  NIPPN(IP2+1).GT.0'
        IP3=IP2+1
        IF(RIN%IPRI_additional_info )                                &
        WRITE(iu_main_output,*)                    &
        'IP3=',IP3,'NIPPN(IP3)=',NIPPN(IP3)           &
       ,'IPPN(IP3,I3):',                              &
        (IPPN(IP3,I3),I3=1,NIPPN(IP3)),'  IPPN(IP2,I3)'
      ENDIF ! NIPPN(IP2+1)

      IF(RIN%IPRI_additional_info)                   &
        WRITE(iu_main_output,*)    &
        ICHECK,npixels*RIN%KNSING,'  ICHECK,npixels*KNSING'
      ELSE !  IF(IP2 .GT. 0)
          WRITE(iu_main_output,*) IP2, " = IP2 - ZERO !!!"
      ENDIF !  IF(IP2 .GT. 0)

      ENDIF ! IF(XTEST .GT. 0 .OR. YTEST .GT. 0 .OR. TTEST .GT. 0)

!**************
! OD ! do not delete
	   !DO ipix=1,npixels
		 !WRITE(*,*) ipix,' ipix'
		 !DO IP=1,NIP_IM(ipix)
		  !WRITE(*,*) ' ', IP,' IP'
		  !DO I=ISTART(IP,ipix),ISTOP(IP,ipix)
		  !WRITE(*,*) '  ', I,(IMAGE_I(ipix,IP,IMAGE_IP),IMAGE_IP=1,NIMAGE_IP(IP,ipix)),' I,IMAGE_I)'
		  !ENDDO
		 !ENDDO
		!ENDDO
! OD
!**************

return
end subroutine time_space_group

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine calculates matrices for applying constraints forcing
        !> @brief smooth connection between the parameters at the edge of the block of 
        !> @brief inverted pixels and those retrieved perviouly outside of this block 
        !> @brief next to the edge 
        !>
        !> @author Oleg Dubovik and Tatsiana Lapionak
        !> @date 20 APR 2015
        !>
        !> @param[in]    iu_main_output - unit number of main output file of retrieval
        !> @param[in]    RIN - structure of settings for retrieval
        !> @param[in]    index_clouds   TO_BE_DOCUMENTED
        !> @param[in]    edges - structure of egdes data
        !> @param[inout] NISM  - number of groups in which smoothness can be implemented
        !> @param[inout] KSM   - number of elements in each of NISM groups 
        !> @param[inout] IKSIM - number of equivalent measurements defining smoothness (used for residual)
        !> @param[inout] IMSM  - pixel numbers which belong to the group
        !> @param[inout] SMMULTI - smoothness matrix
        !> @param[out]   FF_edge        TO_BE_DOCUMENTED
        !> @param[out]   AB_edge        TO_BE_DOCUMENTED 
               
      subroutine smoothterm_mutli_pixel_edges(  iu_main_output,RIN,index_clouds,edges,   & ! IN
                                                NISM,KSM,IKSIM,IMSM,SMMULTI,               & ! INOUT
                                                FF_edge,AB_edge )                ! OUT
!C*********************************************************
!C the subroutine prepares matrices for applying constraints forcing
!C smooth connection between the parameters at the edge of the block of 
!C inverted pixels and those retrieved perviouly outside of this block 
!C next to the edge
!C********************************************************* 

      use mod_par_inv, only : KIX,KIY,KITIME,KIMAGE,   &
                              KIEDGE,KPARS,KMPSM,KDIF
      use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_time_utils
      use mod_edges
      
      implicit none
! ------------------------------------------------------------
! IN :
      integer,                             intent(in)  ::  iu_main_output
      type(ind_clouds),                    intent(in)  ::  index_clouds
      type(retr_input_settings),           intent(in)  ::  RIN
      type(segment_edges),                 intent(in)  ::  edges
! ------------------------------------------------------------
! IN / OUT:
      integer,                             intent(inout) :: NISM,IKSIM
      real,dimension(KIMAGE,KPARS,KIMAGE), intent(inout) :: SMMULTI
      integer,dimension(KMPSM),            intent(inout) :: KSM
      integer,dimension(KMPSM,KIMAGE),     intent(inout) :: IMSM  
      real,dimension(KPARS,KIMAGE),        intent(out)   :: FF_edge
      real,                                intent(out)   :: AB_edge
! ------------------------------------------------------------
! LOCAL :
!      INTEGER,DIMENSION(KIX_E)                :: IXE
!	     INTEGER,DIMENSION(KIY_E)                :: IYE
!      INTEGER,PARAMETER :: KMPSM_E = KIX_E*KIMAGE_E+KIY_E*KIMAGE_E+KIX_E*KIY_E
!      INTEGER,DIMENSION(KMPSM_E,KIMAGE) :: IMSM_edge
!      INTEGER,DIMENSION(KMPSM_E) :: NIO,KSM_AP,KSM_A
!      INTEGER,PARAMETER :: KMPSM_E = KIX_E*KIMAGE_E+KIY_E*KIMAGE_E+KIX_E*KIY_E
      integer,parameter :: KIX_E    = max(KIEDGE,KIX)
      integer,parameter :: KIY_E    = max(KIEDGE,KIY)
      integer,parameter :: KITIME_E = max(KIEDGE,KITIME)
      integer,parameter :: KIMAGE_E = 2*(3*KIX*KITIME+3*KIY*KITIME+3*KIX*KIY)

      integer :: NISM_edge
      integer,dimension(KMPSM,KIMAGE) :: IMSM_edge
      integer,dimension(KMPSM)        :: NIO,KSM_AP,KSM_A
!      integer,dimension(:,:), allocatable :: IMSM_edge
!      integer,dimension(:),   allocatable :: NIO,KSM_AP,KSM_A

      integer,dimension(max(KIMAGE,KIMAGE_E))  :: IX_A,IX_AP,IY_A,IY_AP,IT_A,IT_AP
!      integer,dimension(max(KIX_E,KIY_E,KIMAGE_E),max(KIX_max,KIY_max,KIMAGE_E))  ::    &
      integer,dimension(max(KIX_E,KIY_E,KITIME_E),max(KIX_E,KIY_E,KITIME_E))  ::    &  ! ???
                                                          NX_EDGE_A,NX_EDGE_AP,     &
                                                          NY_EDGE_A,NY_EDGE_AP,     & 
                                                          NT_EDGE_A,NT_EDGE_AP
      real,dimension(KDIF,KDIF)               :: DIF,DIF_A,DIF_AP,SM_A,SM_AP
      real,dimension(KIMAGE,KPARS,KIMAGE)     :: SMIM_edge !???

      real,dimension(KDIF)                    :: XX_edge 
      real,dimension(KIMAGE_E)                :: XX_edge_A !???
      real,dimension(KIMAGE)                  :: XX_edge_AP !???       
      real,dimension(KPARS,KIMAGE_E)          :: APX !???

      real,dimension(:,:,:),  allocatable  :: X_edge,Y_edge
      real,dimension(:),      allocatable  :: T_edge
      integer,dimension(:),   allocatable  :: IT
      real,dimension(:,:,:,:),allocatable  :: AP_edge
      integer                              :: alloc_stat

      integer :: IX,IY,IDIM1,IDIM2,IDIM3,NDIM3,I1,I2,J1,I,J,ii1,KNF
      real    :: G, FF, AB, xnorm
      integer :: NSUM 
      real    :: GSM
      integer :: ISTARSING,IOX,IOY,IOT,NX_edge,NY_edge,    &
                 NT_edge,NIM,NNN_E,IOI,N_I_EDGE,I_EDGE,II_E,IT_edge,IT_main, &
                 IS1, IS2, IKSIM_E,IX_edge,IY_edge,IS_edge,NPIXT,ipix   
      real    :: AB_edge_XYT 
      logical :: IPRINT2, lnism, IPRI_verbose
! ------------------------------------------------------------

      IPRI_verbose = RIN%IPRI_verbose

! ------------------------------------------------------------
!  The below text needs to be updated for "edges"
! ------------------------------------------------------------
! KMPSM - parameter, max number of groups for multi-pixel smoothness
! NISM  - number of groups in which smoothness can be implemented
! IKSIM - number of equivalent measurements defining smoothness
!         (used for residual)
! KSM (1:KMPSM)           - number of elements in each of NISM groups
! IMSM(1:KMPSM,1:KIMAGE)  - pixel numbers which belong to the group
! SMMULTI(1:KIMAGE,1:KPARS,1:KIMAGE)  - smoothness matrix
! ------------------------------------------------------------
!moved into read inversion setting input
!if(KITIME .gt. KDIF) then
!write(*,*) 'KITIME=',KITIME,' .gt. KDIF=',KDIF
!write(*,*) 'Check inversion parameters: KDIF has to be max(KPARS,KITIME,KIX,KIY).'
!stop 'stop in smoothterm_multi_pixel'
!endif
!if(KIX .gt. KDIF) then
!write(*,*) 'KIX=',KIX,' .gt. KDIF=',KDIF
!write(*,*) 'Check inversion parameters: KDIF has to be max(KPARS,KITIME,KIX,KIY).'
!stop 'stop in smoothterm_multi_pixel'
!endif
!if(KIY .gt. KDIF) then
!write(*,*) 'KIY=',KIY,' .gt. KDIF=',KDIF
!write(*,*) 'Check inversion parameters: KDIF has to be max(KPARS,KITIME,KIX,KIY).'
!stop 'stop in smoothterm_multi_pixel'
!endif
!      XNEW(:)   = 0.0

!Allocate X_edge array
      if(.not. allocated(X_edge)) then
        allocate(X_edge(KIX_E,KIY_E,KITIME_E),stat=alloc_stat)       
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate X_edge(:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Allocate Y_edge array
      if(.not. allocated(Y_edge)) then
        allocate(Y_edge(KIX_E,KIY_E,KITIME_E),stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate Y_edge(:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Allocate T_edge array
      if(.not. allocated(T_edge)) then
        allocate(T_edge(KITIME_E),stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate T_edge(:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Allocate IT array
      if(.not. allocated(IT)) then
        allocate(IT(KITIME_E),stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate IT(:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Allocate AP_edge array
      if(.not. allocated(AP_edge)) then
        allocate(AP_edge(KPARS,KIX_E,KIY_E,KITIME_E),stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to allocate AP_edge(:,:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif

      IPRINT2 = .false.
      IKSIM_E = 0
      NISM_edge = 0
      IMSM_edge(:,:) = 0
      FF_edge  (:,:) = 0.0
      AB_edge = 0.0
      SMIM_edge(:,:,:) = 0.0
      xnorm = -999.0 ! smoothness matrix is normalized

      KNF = RIN%KNSINGF
      N_I_EDGE = edges%N_I_EDGE
!write(*,*)  'in smoothterm_mutli_pixel_edges: ',KNF,N_I_EDGE,'  - KNF,N_I_EDGE'
!write(*,*)  'in smoothterm_mutli_pixel_edges: I_EDGE(1:N_I_EDGE):',edges%I_EDGE(1:N_I_EDGE)
!*      NX_before NX_after NY_before NY_after NT_before NT_after

LOOP_I_EDGE: DO ii1=1,N_I_EDGE ! MAX of N_I_EDGE = 6
!*  I_EDGE=1,6 - corresponds to X_before X_after Y_before Y_after T_before T_after
! X - longitude, Y - latitude
          X_edge(:,:,:)    = 0.0
          Y_edge(:,:,:)    = 0.0
          T_edge(:)        = 0.0
          AP_edge(:,:,:,:) = 0.0
          I_EDGE  = edges%I_EDGE(ii1)
!write(*,*) 'in smoothterm_mutli_pixel_edges: ',I_EDGE          
          select case(I_EDGE)
          case(1,2)
! X edges
            II_E = I_EDGE
            NX_edge = edges%group_X(II_E)%nx
            NY_edge = edges%group_X(II_E)%ny
            NT_edge = edges%group_X(II_E)%nt
!write(*,*) 'in smoothterm_mutli_pixel_edges: ',II_E, NX_edge, NY_edge, NT_edge, &
!                                         '   - II_E, NX_edge, NY_edge, NT_edge'
            DO IT_edge=1,NT_edge
              IT(IT_edge) = edges%group_X(II_E)%it(IT_edge)
!write(*,*) 'IT_edge=',IT_edge,'  IT(IT_edge)=',IT(IT_edge)
              X_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_X(II_E)%x(1:NX_edge,1:NY_edge,IT(IT_edge))
              Y_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_X(II_E)%y(1:NX_edge,1:NY_edge,IT(IT_edge)) 
              do IY=1,NY_edge
              do IX=1,NX_edge
                if(edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge)) .eq. 1)  &
                AP_edge(1:KNF,IX,IY,IT(IT_edge)) = LOG(edges%group_X(II_E)%AP(1:KNF,IX,IY,IT(IT_edge)))
!write(*,*) IX,IY,edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge)),  & 
!      '  - IX,IY,ICLOUD_edge(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge))'
              enddo ! IX
              enddo ! IY
            enddo ! IT_edge  
          case(3,4)
! Y_edge
            II_E = I_EDGE-2
            NX_edge = edges%group_Y(II_E)%nx
            NY_edge = edges%group_Y(II_E)%ny
            NT_edge = edges%group_Y(II_E)%nt
!write(*,*) 'in smoothterm_mutli_pixel_edges: ',II_E, NX_edge, NY_edge, NT_edge, &
!                                         '   - II_E, NX_edge, NY_edge, NT_edge'
            DO IT_edge=1,NT_edge
              IT(IT_edge) = edges%group_Y(II_E)%it(IT_edge)
!write(*,*) 'IT_edge=',IT_edge,'  IT(IT_edge)=',IT(IT_edge)
              X_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_Y(II_E)%x(1:NX_edge,1:NY_edge,IT(IT_edge))
              Y_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_Y(II_E)%y(1:NX_edge,1:NY_edge,IT(IT_edge)) 
              do IY=1,NY_edge
              do IX=1,NX_edge
                if(edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge)) .eq. 1)  &
                AP_edge(1:KNF,IX,IY,IT(IT_edge)) = LOG(edges%group_Y(II_E)%AP(1:KNF,IX,IY,IT(IT_edge)))
!write(*,*) IX,IY,edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge)),  & 
!      '  - IX,IY,ICLOUD_edge(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge))'
              enddo ! IX
              enddo ! IY
            enddo ! IT_edge            
          case(5,6)
! T_edge          
            II_E = I_EDGE-4
            NX_edge = edges%group_T(II_E)%nx
            NY_edge = edges%group_T(II_E)%ny
            NT_edge = edges%group_T(II_E)%nt
!write(*,*) 'in smoothterm_mutli_pixel_edges: ',II_E, NX_edge, NY_edge, NT_edge, &
!                                         '   - II_E, NX_edge, NY_edge, NT_edge'
            DO IT_edge=1,NT_edge
              IT(IT_edge) = edges%group_T(II_E)%it(IT_edge)
!write(*,*) 'IT_edge=',IT_edge,'  IT(IT_edge)=',IT(IT_edge)
              X_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_T(II_E)%x(1:NX_edge,1:NY_edge,IT(IT_edge))
              Y_edge(1:NX_edge,1:NY_edge,IT(IT_edge)) = edges%group_T(II_E)%y(1:NX_edge,1:NY_edge,IT(IT_edge)) 
              do IY=1,NY_edge
              do IX=1,NX_edge
                if(edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge)) .eq. 1)  &
                AP_edge(1:KNF,IX,IY,IT(IT_edge)) = LOG(edges%group_T(II_E)%AP(1:KNF,IX,IY,IT(IT_edge)))
!write(*,*) IX,IY,edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge)),  & 
!      '  - IX,IY,ICLOUD_edge(IX,IY,IT(IT_edge)),X_edge(IX,IY,IT(IT_edge)),Y_edge(IX,IY,IT(IT_edge))'
              enddo ! IX
              enddo ! IY
            enddo ! IT_edge            
          end select
          
          select case(I_EDGE)
          case(1) ! I_EDGE=1, X_before **********************************************
!            II_E = ii1
            II_E = I_EDGE
            NX_EDGE_A(:,:)  = 0
            NX_EDGE_AP(:,:) = 0 
            DO IT_edge=1,NT_edge
            IT(IT_edge) = edges%group_X(II_E)%it(IT_edge)
            DO IY=1,NY_edge
              DO IX=NX_edge,1,-1
!              IF(NX_EDGE_A(IY,IT(IT_edge)) .LT. 3) THEN
              IF(NX_EDGE_A(IY,IT(IT_edge)) .LE. KIEDGE) THEN
                IF(edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NX_EDGE_A(IY,IT(IT_edge)) = NX_EDGE_A(IY,IT(IT_edge))+edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge))
                  IX_A(NX_EDGE_A(IY,IT(IT_edge))) = IX
                  XX_edge_A(NX_EDGE_A(IY,IT(IT_edge))) = X_edge(IX,IY,IT(IT_edge))
!write(*,*) IT_edge,IY,IX,IT(IT_edge),XX_edge_A(NX_EDGE_A(IY,IT(IT_edge))),  & 
!        '  IT_edge,IY,IX,IT(IT_edge),XX_edge_A(NX_EDGE_A(IY,IT(IT_edge)))'
                ENDIF
              ENDIF
              ENDDO ! IX
              DO IX=1,index_clouds%NX
!              IF(NX_EDGE_AP(IY,IT(IT_edge)) .LT. 3) THEN
              IF(NX_EDGE_AP(IY,IT(IT_edge)) .LE. KIEDGE) THEN
!write(*,*) IX,IY,IT(IT_edge),index_clouds%ICLOUD(IX,IY,IT(IT_edge)),  & 
!        '  IX,IY,IT(IT_edge),index_clouds%ICLOUD(IX,IY,IT(IT_edge))'
		            IF(index_clouds%ICLOUD(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NX_EDGE_AP(IY,IT(IT_edge)) = NX_EDGE_AP(IY,IT(IT_edge))+index_clouds%ICLOUD(IX,IY,IT(IT_edge))
                  IX_AP(NX_EDGE_AP(IY,IT(IT_edge))) = IX
                  XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge))) = index_clouds%X(IX,IY,IT(IT_edge))
!write(*,*) IT_edge,IY,IX,IT(IT_edge),XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge))),  & 
!        '  IT_edge,IY,IX,IT(IT_edge),XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge)))'
                ENDIF
              ENDIF
              ENDDO ! IX
!write(*,*) 'X_before   **************** IY,IT(IT_edge): ',IY,IT(IT_edge)
!write(*,*) index_clouds%NX,NX_EDGE_A(IY,IT(IT_edge)),'  index_clouds%NX,NX_EDGE_A'
!write(*,*) NX_EDGE_AP(IY,IT(IT_edge)),NX_EDGE_A(IY,IT(IT_edge)),  &
!        '  NX_EDGE_AP(IY,IT(IT_edge)),NX_EDGE_A(IY,IT(IT_edge))'
!write(*,*) 'X_before   ****************'
!stop 'edges test stop'
              IF(NX_EDGE_AP(IY,IT(IT_edge)) .GE. 1 .AND. NX_EDGE_A(IY,IT(IT_edge)) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOX(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMX(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                    lnism = .true.
                    NSUM = NX_EDGE_AP(IY,IT(IT_edge))+NX_EDGE_A(IY,IT(IT_edge))
                    call KSM_calcul(NSUM,NX_EDGE_A(IY,IT(IT_edge)),NX_EDGE_AP(IY,IT(IT_edge)),  &
                                    IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                    XX_edge(:) = 0.0
                    APX(:,:) = 0.0
                    DO IX=KSM_A(NISM_edge),1,-1
                      XX_edge(KSM_A(NISM_edge)-IX+1) = XX_edge_A(IX)
                      APX(1:KNF,KSM_A(NISM_edge)-IX+1) = AP_edge(1:KNF,IX_A(IX),IY,IT(IT_edge))
                    ENDDO ! IX
                    DO IX=1,KSM_AP(NISM_edge)
                      XX_edge(KSM_A(NISM_edge)+IX) = XX_edge_AP(IX)
                      IMSM_edge(NISM_edge,IX) = index_clouds%INIMAGE(IX_AP(IX),IY,IT(IT_edge))
!write(*,*) 'NISM_edge,IX,IMSM_edge(NISM_edge,IX),IT(IT_edge):  ',NISM_edge,IX,IMSM_edge(NISM_edge,IX),IT(IT_edge)
                    ENDDO ! IX
                    DIF(:,:) = 0.0
                    NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                    NIO(NISM_edge) = NNN_E-IOI
!write(*,*) 'X_before, before DIFERW XX_edge: '
!write(*,*) 'NNN_E,IOI:  ',NNN_E,IOI
!write(*,'(10e14.4)') XX_edge(1:2)
!stop
                    CALL DIFERW ( IPRI_verbose, &
                                  NNN_E,IOI,    &
                                  XX_edge(:),G, &
                                  DIF(:,:) )
!write(*,*) 'G=',G
!do i=1,IOI
!write(*,'(a,i4,10e14.4)') 'i,DIF - ',i,DIF(i,1:NNN_E)
!enddo
                    CALL DIF_before_calcul(NIO(NISM_edge),KSM_A(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!***********************************************************************************************
!**** The following quadratic form is calculated for the residal :
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (D1)T (D1) ap - 2 (ap)T (D1)T (D2) ab + (ab)T (D2)T (D2)ab
!****
!**** BELOW we alculate the followibg vlaues:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                    SM_AP(:,:) = 0.0
                    CALL SMOOM( NIO(NISM_edge),KSM_AP(NISM_edge), &
                                DIF_AP(:,:), &
                                SM_AP(:,:), xnorm )
!write(*,*) 'NIO(NISM_edge),KSM_AP(NISM_edge): ',NIO(NISM_edge),KSM_AP(NISM_edge)
!write(*,*) 'SM_AP:'
!do i=1,KSM_AP(NISM_edge)
!write(*,*) 'i=',i,SM_AP(i,1:KSM_AP(NISM_edge))
!enddo
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                    SM_A(:,:) = 0.0
                    CALL MAT_T_MAT ( KSM_AP(NISM_edge),KSM_A(NISM_edge),   &
                                     NIO(NISM_edge),                       &
                                     DIF_AP(:,:),DIF_A(:,:),SM_A(:,:) )
!write(*,*) 'KSM_AP(NISM_edge),KSM_A(NISM_edge): ',KSM_AP(NISM_edge),KSM_A(NISM_edge)
!write(*,*) 'SM_A:'
!do i=1,KSM_AP(NISM_edge)
!write(*,*) 'i=',i,SM_A(i,1:KSM_A(NISM_edge))
!enddo

                    IKSIM_E = IKSIM_E+NIO(NISM_edge)
                    DO I2=1,NDIM3
                      AB_edge = AB_edge +  &
                      AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                      CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                          KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Here we define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be added to general a priori term         ****
                      CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                                KSM_AP(NISM_edge),SM_AP,SMIM_edge)
!write(*,*) 'case(1) I2=',I2,'  after CALL smoothness_edge_term'
                    ENDDO ! I2=1,NDIM3
                  ENDIF ! IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if(.not. lnism) NISM_edge = NISM_edge-1
              ENDIF ! NX_EDGE_AP(IY,IT(IT_edge)) .GE. 1 .AND. 
            ENDDO ! DO IY=1,NY_edge
            ENDDO ! DO IT_edge=1,NT_edge
!write(*,*) 'IS1=',IS1,'  FF_edge(:,1): '
!write(*,'(10e14.4)')  FF_edge(:,1)
!stop 'case(1): test stop in smoothterm_mutli_pixel_edges'
          case(2) ! I_EDGE=2, X_after **********************************************
!            II_E = ii1
            II_E = I_EDGE
            NX_EDGE_A(:,:)  = 0
            NX_EDGE_AP(:,:) = 0
            DO IT_edge=1,NT_edge
            IT(IT_edge) = edges%group_X(II_E)%it(IT_edge)
            DO IY=1,NY_edge
              DO IX=1,NX_edge
!              IF(NX_EDGE_A(IY,IT(IT_edge)) .LT. 3) THEN
              IF(NX_EDGE_A(IY,IT(IT_edge)) .LE. KIEDGE) THEN
                IF(edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NX_EDGE_A(IY,IT(IT_edge)) = NX_EDGE_A(IY,IT(IT_edge))+edges%group_X(II_E)%icloud(IX,IY,IT(IT_edge))
                  IX_A(NX_EDGE_A(IY,IT(IT_edge))) = IX
                  XX_edge_A(NX_EDGE_A(IY,IT(IT_edge))) = X_edge(IX,IY,IT(IT_edge))
!write(*,*) IT_edge,IY,IX,IT(IT_edge),XX_edge_A(NX_EDGE_A(IY,IT(IT_edge))),  & 
!        '  IT_edge,IY,IX,IT(IT_edge),XX_edge_A(NX_EDGE_A(IY,IT(IT_edge)))'
                ENDIF
              ENDIF
              ENDDO ! IX
              DO IX=index_clouds%NX,1,-1
!              IF(NX_EDGE_AP(IY,IT(IT_edge)) .LT. 3) THEN
              IF(NX_EDGE_AP(IY,IT(IT_edge)) .LE. KIEDGE) THEN
                IF(index_clouds%ICLOUD(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NX_EDGE_AP(IY,IT(IT_edge)) = NX_EDGE_AP(IY,IT(IT_edge))+index_clouds%ICLOUD(IX,IY,IT(IT_edge))
                  IX_AP(NX_EDGE_AP(IY,IT(IT_edge))) = IX
                  XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge))) = index_clouds%X(IX,IY,IT(IT_edge))
!write(*,*) IT_edge,IY,IX,IT(IT_edge),XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge))),  & 
!        '  IT_edge,IY,IX,IT(IT_edge),XX_edge_AP(NX_EDGE_AP(IY,IT(IT_edge)))'
                ENDIF
              ENDIF
              ENDDO ! IX
!write(*,*) 'X_after    **************** IY,IT(IT_edge): ',IY,IT(IT_edge)
!write(*,*) index_clouds%NX,NX_EDGE_A(IY,IT(IT_edge)),'  NX_EDGE_AP,NX_EDGE_A'
!write(*,*) NX_EDGE_AP(IY,IT(IT_edge)),NX_EDGE_A(IY,IT(IT_edge)),  & 
!        '  NX_EDGE_AP(IY,IT(IT_edge)),NX_EDGE_A(IY,IT(IT_edge))'
!write(*,*) 'X_after    ****************'
!stop 'X_after edges test stop'
              IF(NX_EDGE_AP(IY,IT(IT_edge)) .GE. 1 .AND. NX_EDGE_A(IY,IT(IT_edge)) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOX(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMX(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                    lnism = .true.
                    NSUM = NX_EDGE_AP(IY,IT(IT_edge))+NX_EDGE_A(IY,IT(IT_edge))
                    call KSM_calcul(NSUM,NX_EDGE_A(IY,IT(IT_edge)),NX_EDGE_AP(IY,IT(IT_edge)),  &
                                    IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                    XX_edge(:) = 0.0
                    APX(:,:) = 0.0
                    DO IX=1,KSM_A(NISM_edge)
                      XX_edge(KSM_AP(NISM_edge)+IX) = XX_edge_A(IX)
                      APX(1:KNF,IX) = AP_edge(1:KNF,IX_A(IX),IY,IT(IT_edge))
                    ENDDO ! IX

                    DO IX=KSM_AP(NISM_edge),1,-1
!tlod                    DO IX=1,KSM_AP(NISM_edge)
!tlod                      XX_edge(KSM_AP(NISM_edge)-IX+1) = XX_edge_AP(index_clouds%NX-IX+1)
                      XX_edge(IX) = XX_edge_AP(KSM_AP(NISM_edge)-IX+1)
                      IMSM_edge(NISM_edge,IX) = index_clouds%INIMAGE(IX_AP(IX),IY,IT(IT_edge))
!write(*,*) 'NISM_edge,KSM_AP(NISM_edge)-IX+1,IMSM_edge(NISM_edge,KSM_AP(NISM_edge)-IX+1),IT(IT_edge):  ',  & 
!            NISM_edge,KSM_AP(NISM_edge)-IX+1,IMSM_edge(NISM_edge,KSM_AP(NISM_edge)-IX+1),IT(IT_edge)
!write(*,*)  IX,NISM_edge,IMSM_edge(NISM_edge,KSM_AP(NISM_edge)-IX+1),index_clouds%INIMAGE(IX_AP(IX),IY,IT(IT_edge)), &
!         '  IX,NISM_edge,IMSM_edge(NISM_edge,KSM_AP(NISM_edge)-IX+1),index_clouds%INIMAGE(IX_AP(IX),IY,IT(IT_edge))'
                    ENDDO

                    DIF(:,:) = 0.0
                    !NNN_E = NISM_edge+KSM_A(NISM_edge)
                    NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                    NIO(NISM_edge) = NNN_E-IOI
!write(*,*) 'X_after : ',NISM_edge,KSM_AP(NISM_edge),KSM_A(NISM_edge),'  NISM_edge,KSM_AP(NISM_edge),KSM_A(NISM_edge)'
!write(*,*) 'NNN_E,IOI  ',NNN_E,IOI
!write(*,*) 'XX_edge:'
!write(*,'(10e14.4)') XX_edge(1:NNN_E) 
                    CALL DIFERW ( IPRI_verbose, &
                                  NNN_E,IOI,    &
                                  XX_edge(:),G, &
                                  DIF(:,:) )
!write(*,*) 'G=',G,'  DIF: '
!do i=1,IOI
!write(*,'(10e14.4)') DIF(i,1:NNN_E)
!enddo
                    CALL DIF_after_calcul(NIO(NISM_edge),KSM_AP(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!***********************************************************************************************
!**** For the residal we need calculate the following quadratic form:
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (S1)T (S1) ap - 2 (ap)T (S1)T (S2) ab + (ab)T (S2)T (S2)ab
!****
!**** BELOW we calculate the following values:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                    SM_AP(:,:) = 0.0
                    CALL SMOOM( NIO(NISM_edge),NNN_E,  &
                                DIF_AP(:,:),           &
                                SM_AP(:,:), xnorm )
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                    SM_A(:,:) = 0.0
                    CALL MAT_T_MAT( KSM_AP(NISM_edge),KSM_A(NISM_edge),   &
                                    NIO(NISM_edge),                       &
                                    DIF_AP(:,:),DIF_A(:,:),SM_A(:,:) )
                    IKSIM_E = IKSIM_E+NIO(NISM_edge)
                    DO I2=1,NDIM3
                      AB_edge = AB_edge +  &
                      AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                      CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                          KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Here we define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be addred to general a priori term         ****
                      CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                                KSM_AP(NISM_edge),SM_AP,SMIM_edge)
                    ENDDO ! I2=1,NDIM3
!                   ENDDO
                  ENDIF ! IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if (.not. lnism) NISM_edge =NISM_edge-1
              ENDIF ! NX_EDGE_AP(IT,IY) .GE. 1 .AND. 
            ENDDO ! DO IY=1,NY_edge
            ENDDO ! DO IT_edge=1,NT_edge
!stop 'case(2): test stop in smoothterm_mutli_pixel_edges'
          case(3) ! I_EDGE=3, Y_before *********************************************
!            II_E = ii1-2
            II_E = I_EDGE-2
            NY_EDGE_A(:,:)  = 0
            NY_EDGE_AP(:,:) = 0
            DO IT_edge=1,NT_edge
!            write(*,*) 'II_E,IT_edge,NT_edge: ',II_E,IT_edge,NT_edge
            IT(IT_edge) = edges%group_Y(II_E)%it(IT_edge)
            DO IX=1,NX_edge
              DO IY=NY_edge,1,-1
!              IF(NY_EDGE_A(IX,IT(IT_edge)) .LT. 3) THEN
              IF(NY_EDGE_A(IX,IT(IT_edge)) .LE. KIEDGE) THEN
                IF(edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NY_EDGE_A(IX,IT(IT_edge)) = NY_EDGE_A(IX,IT(IT_edge))+edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge))
                  IY_A(NY_EDGE_A(IX,IT(IT_edge))) = IY
                  XX_edge_A(NY_EDGE_A(IX,IT(IT_edge))) = Y_edge(IX,IY,IT(IT_edge))
!write(*,*) NX_EDGE_A(IX,IT(IT_edge)),IY_A(NX_EDGE_A(IX,IT(IT_edge))),       &
!        '  NX_EDGE_A(IX,IT(IT_edge)),IY_A(NX_EDGE_A(IX,IT(IT_edge)))'
!write(*,*) IT_edge,IX,IY,IT(IT_edge),XX_edge_A(NY_EDGE_A(IX,IT(IT_edge))),  & 
!        '  IT_edge,IX,IY,IT(IT_edge),,XX_edge_A(NY_EDGE_A(IX,IT(IT_edge)))'

                ENDIF
              ENDIF
              ENDDO ! IY
              DO IY=1,index_clouds%NY
!              IF(NY_EDGE_AP(IX,IT(IT_edge)) .LT. 3) THEN
              IF(NY_EDGE_AP(IX,IT(IT_edge)) .LE. KIEDGE) THEN
                IF(index_clouds%ICLOUD(IX,IY,IT(IT_edge)) .GT. 0) THEN
                  NY_EDGE_AP(IX,IT(IT_edge)) = NY_EDGE_AP(IX,IT(IT_edge))+index_clouds%ICLOUD(IX,IY,IT(IT_edge))
                  IY_AP(NY_EDGE_AP(IX,IT(IT_edge))) = IY
                  XX_edge_AP(NY_EDGE_AP(IX,IT(IT_edge))) = index_clouds%Y(IX,IY,IT(IT_edge))
!write(*,*) NY_EDGE_AP(IX,IT(IT_edge)),IY_AP(NY_EDGE_AP(IX,IT(IT_edge))),      & 
!        '  NY_EDGE_AP(IX,IT(IT_edge)),IY_AP(NY_EDGE_AP(IX,IT(IT_edge)))'
!write(*,*) IT_edge,IX,IY,IT(IT_edge),XX_edge_AP(NY_EDGE_AP(IX,IT(IT_edge))),  & 
!        '  IT_edge,IX,IY,IT(IT_edge),XX_edge_AP(NY_EDGE_AP(IX,IT(IT_edge)))'
                ENDIF
              ENDIF
              ENDDO ! IY
!write(*,*) 'Y_before   **************** IX,IT(IT_edge): ',IX,IT(IT_edge)
!write(*,*) index_clouds%NY,NY_EDGE_A(IX,IT(IT_edge)),'  NY_EDGE_AP,NY_EDGE_A'
!write(*,*) NY_EDGE_AP(IX,IT(IT_edge)),NY_EDGE_A(IX,IT(IT_edge)),  & 
!        '  NY_EDGE_AP(IX,IT(IT_edge)),NY_EDGE_A(IX,IT(IT_edge))'
!write(*,*) 'Y_before   ****************'
!stop 'Y_before edges test stop'
              IF(NY_EDGE_AP(IX,IT(IT_edge)) .GE. 1 .AND. NY_EDGE_A(IX,IT(IT_edge)) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOY(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMY(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                  lnism = .true.
                  NSUM = NY_EDGE_AP(IX,IT(IT_edge))+NY_EDGE_A(IX,IT(IT_edge))
                  call KSM_calcul(NSUM,NY_EDGE_A(IX,IT(IT_edge)),NY_EDGE_AP(IX,IT(IT_edge)),  &
                                  IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                  XX_edge(:) = 0.0
                  APX(:,:) = 0.0
                  DO IY=KSM_A(NISM_edge),1,-1
                    XX_edge(KSM_A(NISM_edge)-IY+1) = XX_edge_A(IY)
                    APX(1:KNF,KSM_A(NISM_edge)-IY+1) = AP_edge(1:KNF,IX,IY_A(IY),IT(IT_edge))
                  ENDDO
                  DO IY=1,KSM_AP(NISM_edge)
                    XX_edge(KSM_A(NISM_edge)+IY)=XX_edge_AP(IY)
                    IMSM_edge(NISM_edge,IY) = index_clouds%INIMAGE(IX,IY_AP(IY),IT(IT_edge))
                  ENDDO

                  DIF(:,:) = 0.0
                  NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                  NIO(NISM_edge) = NNN_E-IOI
                  CALL DIFERW ( IPRI_verbose, &
                                NNN_E,IOI,    &
                                XX_edge(:),G, &
                                DIF(:,:) )
                  CALL DIF_before_calcul(NIO(NISM_edge),KSM_A(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!write(*,*) 'IKSIM_E, NIO(NISM_edge), NNN_E, IOI - ',IKSIM_E,NIO(NISM_edge),NNN_E,IOI
!do i=1,IOI
!write(*,'(a,i4,10e14.4)') 'i,DIF_A - ',i,DIF_A(i,1:NNN_E)
!enddo
!do i=1,IOI
!write(*,'(a,i4,10e14.4)') 'i,DIF_AP - ',i,DIF_AP(i,1:NNN_E)
!enddo
!***********************************************************************************************
!**** For the residal we need calculate the following quadratic form:
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (S1)T (S1) ap - 2 (ap)T (S1)T (S2) ab + (ab)T (S2)T (S2)ab
!****
!**** BELOW we alculate the followibg vlaues:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                  SM_AP(:,:) = 0.0
!od&tl                  CALL SMOOM (  NIO(NISM_edge),NNN_E,  &
                  CALL SMOOM( NIO(NISM_edge),KSM_AP(NISM_edge), &
                              DIF_AP(:,:), &
                              SM_AP(:,:), xnorm )
!write(*,*) 'NIO(NISM_edge),KSM_AP(NISM_edge): ',NIO(NISM_edge),KSM_AP(NISM_edge)
!write(*,*) 'SM_AP:'
!do i=1,KSM_AP(NISM_edge)
!write(*,*) 'i=',i,SM_AP(i,1:KSM_AP(NISM_edge))
!enddo
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                  SM_A(:,:) = 0.0
                  CALL MAT_T_MAT( KSM_AP(NISM_edge),KSM_A(NISM_edge),  &
                                  NIO(NISM_edge),                      &
                                  DIF_AP(:,:),DIF_A(:,:),SM_A(:,:) )
!write(*,*) 'KSM_AP(NISM_edge),KSM_A(NISM_edge): ',KSM_AP(NISM_edge),KSM_A(NISM_edge)
!write(*,*) 'SM_A:'
!do i=1,KSM_AP(NISM_edge)
!write(*,*) 'i=',i,SM_A(i,1:KSM_A(NISM_edge))
!enddo ! i
                  IKSIM_E = IKSIM_E+NIO(NISM_edge)
                  DO I2=1,NDIM3
                    AB_edge = AB_edge +  &
                    AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                    CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                        KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Herewe define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be addred to general a priori term         ****
                    CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                              KSM_AP(NISM_edge),SM_AP,SMIM_edge)
!write(*,*) 'case(3) I2=',I2,'  after CALL smoothness_edge_term'
                  ENDDO ! I2=1,NDIM3!                   ENDDO
                ENDIF ! NY_EDGE_AP(IT(IT_edge),IX) .GE. 1 .AND.     ?????
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if(.not. lnism) NISM_edge = NISM_edge-1
              ENDIF ! NY_EDGE_AP(IT(IT_edge),IX) .GE. 1 .AND.
            ENDDO ! DO IX=1,NX_edge
            ENDDO ! DO IT_edge=1,NT_edge
!stop 'Y_before edges test stop'
          case(4) ! I_EDGE=4, YB after main retrieved block ************************
!            II_E = ii1-2
            II_E = I_EDGE-2
            NY_EDGE_A(:,:)  = 0.0
            NY_EDGE_AP(:,:) = 0.0
            DO IT_edge=1,NT_edge
            IT(IT_edge) = edges%group_Y(II_E)%it(IT_edge)
            IT(IT_edge) = IT(IT_edge)
            DO IX=1,NX_edge
              DO IY=1,NY_edge
!                IF(NY_EDGE_A(IX,IT(IT_edge)) .LT. 3) THEN
                IF(NY_EDGE_A(IX,IT(IT_edge)) .LE. KIEDGE) THEN
                  IF(edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
                    NY_EDGE_A(IX,IT(IT_edge)) = NY_EDGE_A(IX,IT(IT_edge))+edges%group_Y(II_E)%icloud(IX,IY,IT(IT_edge))
                    IY_A(NY_EDGE_A(IX,IT(IT_edge))) = IY
                    XX_edge_A(NY_EDGE_A(IX,IT(IT_edge))) = X_edge(IX,IY,IT(IT_edge))
                  ENDIF
                ENDIF
              ENDDO
              DO IY=index_clouds%NY,1,-1
!                IF(NY_EDGE_AP(IX,IT(IT_edge)) .LT. 3) THEN
                IF(NY_EDGE_AP(IX,IT(IT_edge)) .LE. KIEDGE) THEN
                  IF(index_clouds%ICLOUD(IX,IY,IT(IT_edge)) .GT. 0) THEN
                    NY_EDGE_AP(IX,IT(IT_edge)) = NY_EDGE_AP(IX,IT(IT_edge))+index_clouds%ICLOUD(IX,IY,IT(IT_edge))
                    IY_AP(NY_EDGE_AP(IX,IT(IT_edge))) = IY
                    XX_edge_AP(NY_EDGE_AP(IX,IT(IT_edge))) = index_clouds%Y(IX,IY,IT(IT_edge))
                  ENDIF
                ENDIF
              ENDDO
!write(*,*) 'Y_after   **************** IX,IT(IT_edge): ',IX,IT(IT_edge)
!write(*,*) index_clouds%NY,NY_EDGE_A(IX,IT(IT_edge)),'  NY_EDGE_AP,NY_EDGE_A'
!write(*,*) NY_EDGE_AP(IX,IT(IT_edge)),NY_EDGE_A(IX,IT(IT_edge)),  & 
!        '  NY_EDGE_AP(IX,IT(IT_edge)),NY_EDGE_A(IX,IT(IT_edge))'
!write(*,*) 'Y_after   ****************'
!stop 'Y_before edges test stop'
              IF(NY_EDGE_AP(IX,IT(IT_edge)) .GE. 1 .AND. NY_EDGE_A(IX,IT(IT_edge)) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOY(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMY(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                    lnism = .true.
                    NSUM = NY_EDGE_AP(IX,IT(IT_edge))+NY_EDGE_A(IX,IT(IT_edge))
                    call KSM_calcul(NSUM,NY_EDGE_A(IX,IT(IT_edge)),NY_EDGE_AP(IX,IT(IT_edge)),  &
                                    IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                    XX_edge(:) = 0.0
                    APX(:,:) = 0.0
                    DO IY=1,KSM_A(NISM_edge)
                      XX_edge(KSM_AP(NISM_edge)+IY) = XX_edge_A(IY)
                      APX(1:KNF,IX) = AP_edge(1:KNF,IX,IY_A(IY),IT(IT_edge))
                    ENDDO
!                    DO IY=KSM_AP(NISM_edge),1,-1
                    DO IY=1,KSM_AP(NISM_edge)
                      XX_edge(IY) = XX_edge_AP(KSM_AP(NISM_edge)-IY+1)
                      IMSM_edge(NISM_edge,IY) = index_clouds%INIMAGE(IY,IY_AP(IY),IT(IT_edge))
                    ENDDO
                    DIF(:,:) = 0.0
                    NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                    NIO(NISM_edge) = NNN_E-IOI
                    CALL DIFERW ( IPRI_verbose, &
                                  NNN_E,IOI,    &
                                  XX_edge(:),G, &
                                  DIF(:,:) )
                    CALL DIF_after_calcul(NIO(NISM_edge),KSM_AP(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!***********************************************************************************************
!**** For the residal we need calculate the following quadratic form:
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (S1)T (S1) ap - 2 (ap)T (S1)T (S2) ab + (ab)T (S2)T (S2)ab
!****
!**** BELOW we alculate the followibg vlaues:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                    SM_AP(:,:) = 0.0
                    CALL SMOOM( NIO(NISM_edge),NNN_E, &
                                DIF_AP(:,:), &
                                SM_AP(:,:), xnorm )
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                    SM_A(:,:) = 0.0
                    CALL MAT_T_MAT( KSM_AP(NISM_edge),KSM_A(NISM_edge),   &
                                    NIO(NISM_edge),                       &
                                    DIF_AP(:,:),DIF_A(:,:),SM_A(:,:) )
                    IKSIM_E = IKSIM_E+NIO(NISM_edge)
                    DO I2=1,NDIM3
                      AB_edge = AB_edge +  &
                      AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                      CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                          KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Here we define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be addred to general a priori term         ****
                      CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                                KSM_AP(NISM_edge),SM_AP,SMIM_edge)
                    ENDDO ! I2=1,NDIM3

                  ENDIF ! IF((GSM.GT.0.).AND.IOY.GT.0) THEN
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if(.not. lnism) NISM_edge = NISM_edge-1 
              ENDIF ! NX_EDGE_AP(IT,IY) .GE. 1 .AND. 
            ENDDO ! DO IX=1,NX_edge
            ENDDO ! DO IT=1,NT_edge
!stop 'Y_after edges test stop'
          case(5) ! I_EDGE=5, TB before main retrieved block !****************************
!            II_E = ii1-4
            II_E = I_EDGE-4
            NT_EDGE_A(:,:)  = 0
            NT_EDGE_AP(:,:) = 0
!tl            DO IX=1,NX_edge
            DO IY=1,NY_edge
            DO IX=1,NX_edge
              DO IT_edge=NT_edge,1,-1
!                IF(NT_EDGE_A(IX,IY) .LT. 3) THEN
                IF(NT_EDGE_A(IX,IY) .LE. KIEDGE) THEN
                  IF(edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
!write(*,*) 'IX=',IX,'  IY=',IY,'  NX_EDGE_A(IX,IY)=',NX_EDGE_A(IX,IY)
!write(*,*) 'IT_edge=',IT_edge,'  IT(IT_edge)=',IT(IT_edge)
                    NT_EDGE_A(IX,IY) = NT_EDGE_A(IX,IY)+edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge))
                    IT_A(NT_EDGE_A(IX,IY)) = IT(IT_edge)
                    XX_edge_A(NT_EDGE_A(IX,IY)) = T_edge(IT_edge)
                  ENDIF
                ENDIF
              ENDDO
              DO IT_main=1,index_clouds%NT
!                IF(NT_EDGE_AP(IX,IY) .LT. 3) THEN
                IF(NT_EDGE_AP(IX,IY) .LE. KIEDGE) THEN
                  IF(index_clouds%ICLOUD(IX,IY,IT_main) .GT. 0) THEN
                    NT_EDGE_AP(IX,IY) = NT_EDGE_AP(IX,IY)+index_clouds%ICLOUD(IX,IY,IT_main)
                    IT_AP(NT_EDGE_AP(IX,IY)) = IT_main
                    XX_edge_AP(NT_EDGE_AP(IX,IY)) = index_clouds%T(IX,IY,IT_main)
                  ENDIF
                ENDIF
              ENDDO
!write(*,*) 'T_before   **************** IX,IY: ',IX,IY
!write(*,*) index_clouds%NT,NT_EDGE_A(IX,IY),'  NT_EDGE_AP,NT_EDGE_A'
!write(*,*) NT_EDGE_AP(IX,IY),NT_EDGE_A(IX,IY),  & 
!        '  NT_EDGE_AP(IX,IY),NT_EDGE_A(IX,IY)'
!write(*,*) 'T_before   ****************'
!stop 'T_before edges test stop'
              IF(NT_EDGE_AP(IX,IY) .GE. 1 .AND. NT_EDGE_A(IX,IY) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOT(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMT(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                    lnism = .true.
                    NSUM = NT_EDGE_AP(IX,IY)+NT_EDGE_A(IX,IY)
                    call KSM_calcul(NSUM,NT_EDGE_A(IX,IY),NT_EDGE_AP(IX,IY),  &
                                    IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                    XX_edge(:) = 0.0
                    APX(:,:) = 0.0
                    DO IT_edge=KSM_A(NISM_edge),1,-1
                      XX_edge(KSM_A(NISM_edge)-IT_edge+1) = XX_edge_A(IT_edge)
                      APX(1:KNF,KSM_A(NISM_edge)-IT_edge+1) = AP_edge(1:KNF,IX,IY,IT_edge)
                    ENDDO
!write(*,*) 'NISM_edge,KSM_AP(NISM_edge): ',NISM_edge,KSM_AP(NISM_edge)
                    DO IT_main=1,KSM_AP(NISM_edge)
                      XX_edge(KSM_A(NISM_edge)+(RIN%edges%nt+IT_main)) = XX_edge_AP(IT_main)
                      IMSM_edge(NISM_edge,IT_main) = index_clouds%INIMAGE(IX,IY,IT_AP(IT_main))
!write(*,*) 'NISM_edge,IMSM_edge(NISM_edge,IT_main),IT_main:  ',  &
!            NISM_edge,IMSM_edge(NISM_edge,IT_main),IT_main
                    ENDDO
                    DIF(:,:) = 0.0
                    NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                    NIO(NISM_edge) = NNN_E-IOI
                    CALL DIFERW ( IPRI_verbose, &
                                  NNN_E,IOI,    &
                                  XX_edge(:),G, &
                                  DIF(:,:) )
                    CALL DIF_before_calcul(NIO(NISM_edge),KSM_A(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!***********************************************************************************************
!**** For the residal we need calculate the following quadratic form:
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (S1)T (S1) ap - 2 (ap)T (S1)T (S2) ab + (ab)T (S2)T (S2)ab
!****
!**** BELOW we alculate the followibg vlaues:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                    SM_AP(:,:) = 0.0
!                    CALL SMOOM (  NIO(NISM_edge),NNN_E,  &
                    CALL SMOOM( NIO(NISM_edge),KSM_AP(NISM_edge), &
                                DIF_AP(:,:), &
                                SM_AP(:,:), xnorm )
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                    SM_A(:,:) = 0.0
                    CALL MAT_T_MAT (  KSM_AP(NISM_edge),KSM_A(NISM_edge),  &
                                      NIO(NISM_edge),                      &
                                      DIF_AP(:,:),DIF_A(:,:),SM_A(:,:) )
                    IKSIM_E = IKSIM_E+NIO(NISM_edge)
                    DO I2=1,NDIM3
                      AB_edge = AB_edge +  &
                      AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                      CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                          KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Here we define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be addred to general a priori term         ****
                      CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                                KSM_AP(NISM_edge),SM_AP,SMIM_edge)
                    ENDDO ! I2=1,NDIM3
                  ENDIF ! GSM .GT. 0. .AND.
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if(.not. lnism) NISM_edge = NISM_edge
              ENDIF ! NX_EDGE_AP(IT,IY) .GE. 1 .AND. 
            ENDDO ! IY=1,NY_edge
            ENDDO ! DO IX=1,NX_edge
!stop 'T_before edges test stop'
          case(6) ! I_EDGE=6, T_after ******************************************************
!            II_E = ii1-4
            II_E = I_EDGE-4
            NT_EDGE_A(:,:)  = 0
            NT_EDGE_AP(:,:) = 0
            DO IY=1,NY_edge
            DO IX=1,NX_edge
              DO IT_edge=1,NT_edge
!                IF(NT_EDGE_A(IX,IY) .LT. 3) THEN
                IF(NT_EDGE_A(IX,IY) .LE. KIEDGE) THEN
                  IF(edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge)) .GT. 0) THEN
                    NT_EDGE_A(IX,IY) = NT_EDGE_A(IX,IY)+edges%group_T(II_E)%icloud(IX,IY,IT(IT_edge))
                    IT_A(NT_EDGE_A(IX,IY)) = IT(IT_edge)
                    XX_edge_A(NT_EDGE_A(IX,IY)) = T_edge(IT_edge)
                  ENDIF
                ENDIF
              ENDDO
              DO IT_main=index_clouds%NT,1,-1
!                IF(NT_EDGE_AP(IX,IY) .LT. 3) THEN
                IF(NT_EDGE_AP(IX,IY) .LE. KIEDGE) THEN
                  IF(index_clouds%ICLOUD(IX,IY,IT_main) .GT. 0) THEN
                    NT_EDGE_AP(IX,IY) = NT_EDGE_AP(IX,IY)+index_clouds%ICLOUD(IX,IY,IT_main)
                    IT_AP(NT_EDGE_AP(IX,IY)) = IT_main
                    XX_edge_AP(NT_EDGE_AP(IX,IY)) = index_clouds%T(IX,IY,IT_main)
                  ENDIF
                ENDIF
              ENDDO
!write(*,*) 'T_after   **************** IX,IY: ',IX,IY
!write(*,*) index_clouds%NT,NT_EDGE_A(IX,IY),'  NT_EDGE_AP,NT_EDGE_A'
!write(*,*) NT_EDGE_AP(IX,IY),NT_EDGE_A(IX,IY),  & 
!        '  NT_EDGE_AP(IX,IY),NT_EDGE_A(IX,IY)'
!write(*,*) 'T_after   ****************'
!stop 'T_after edges test stop'
              IF(NT_EDGE_AP(IX,IY) .GE. 1 .AND. NT_EDGE_A(IX,IY) .GE. 1) THEN
                lnism = .false.
                NISM_edge = NISM_edge+1
                DO IDIM1=1,RIN%NDIM%n1
                DO IDIM2=1,RIN%NDIM%n2(IDIM1)
                  NDIM3     = RIN%NDIM%n3(IDIM2,IDIM1)
                  ISTARSING = RIN%NDIM%ISTARSING(IDIM2,IDIM1)
                  IOI       = RIN%MPCS%IOT(IDIM2,IDIM1)
                  GSM       = RIN%MPCS%GSMT(IDIM2,IDIM1)
                  IF(GSM .GT. 0. .AND. IOI .GT. 0) THEN
                    lnism = .true.
                    NSUM = NT_EDGE_AP(IX,IY)+NT_EDGE_A(IX,IY)
                    call KSM_calcul(NSUM,NT_EDGE_A(IX,IY),NT_EDGE_AP(IX,IY),  &
                                    IOI,KSM_A(NISM_edge),KSM_AP(NISM_edge))
                    XX_edge(:) = 0.0
                    APX(:,:) = 0.0
                    DO IT_edge=1,KSM_A(NISM_edge)
                      XX_edge(KSM_AP(NISM_edge)+IT_edge) = XX_edge_A(IT_edge)
                      APX(1:KNF,IT_edge) = AP_edge(1:KNF,IX,IY,IT_edge)
                    ENDDO
                    DO IT_main=1,KSM_AP(NISM_edge)
!tl                      XX_edge(KSM_AP(NISM_edge)-IT_m(IT_main)+1) = XX_edge_AP(index_clouds%NT-IT_main+1)
                      XX_edge(IT_main) = XX_edge_AP(KSM_AP(NISM_edge)-IT_main+1)
!write(*,*) 'IT_main=',IT_main,'  IT_AP(IT_main)=',IT_AP(IT_main),  & 
!   '  index_clouds%INIMAGE(IX,IY,IT_AP(IT_main))=',index_clouds%INIMAGE(IX,IY,IT_AP(IT_main))
!write(*,*) 'NISM_edge=',NISM_edge,'  index_clouds%NT-IT_main+1=',index_clouds%NT-IT_main+1
                      !IMSM_edge(NISM_edge,index_clouds%NT-IT_main+1) = index_clouds%INIMAGE(IX,IY,IT_AP(IT_main))
                      IMSM_edge(NISM_edge,IT_main) = index_clouds%INIMAGE(IX,IY,IT_AP(IT_main))
                    ENDDO

                    DIF(:,:) = 0.0
                    NNN_E = KSM_AP(NISM_edge)+KSM_A(NISM_edge)
                    NIO(NISM_edge) = NNN_E-IOI
                    CALL DIFERW ( IPRI_verbose, &
                                  NNN_E,IOI,    &
                                  XX_edge(:),G, &
                                  DIF(:,:) )
                    CALL DIF_after_calcul(NIO(NISM_edge),KSM_AP(NISM_edge),NNN_E,DIF,DIF_A,DIF_AP)
!***********************************************************************************************
!**** For the residal we need calculate the following quadratic form:
!**** ((D1)* ap - (D2) * ab)T ((D1)* ap - (D2) * ab)
!**** The equivalent is:  (ap)T (S1)T (S1) ap - 2 (ap)T (S1)T (S2) ab + (ab)T (S2)T (S2)ab
!****
!**** BELOW we alculate the followibg vlaues:
!**** - SM_AP = (D1)T (D1), SMIM_edge - represent the total term from all I_EDGE;
!**** - FF_B = (ab)T (D2)T (D2)ab = ((D2)ab)**2- the total term from all I_EDGE;
!**** - FF   - (D1)T (D2) ab      - the total term from all I_EDGE;
!****   (ap)t FF - needs to be calcated for contribtuion of edge constraints into total residual;
!**** - Aslo: (-1) (D1)T (D2) ab  is the  a priori edge  term contribtuion to the right side  ****
!****                                                                     of Normal System:
!*************************************************************************************************
!**** Here we define SM_AP= (DIF_AP)**T (DIF_AP)   *********
                    SM_AP(:,:) = 0.0
                    CALL SMOOM( NIO(NISM_edge),NNN_E, &
                                DIF_AP(:,:), &
                                SM_AP(:,:), xnorm )
!**** Here we define SM_A= (DIF_AP)**T DIF_A        *********
                    SM_A(:,:) = 0.0
                    CALL MAT_T_MAT (  KSM_AP(NISM_edge),KSM_A(NISM_edge),  &
                                      NIO(NISM_edge),                      &
                                      DIF_AP(:,:),DIF_A(:,:),SM_A(:,:)  )
                    IKSIM_E = IKSIM_E+NIO(NISM_edge)
                    DO I2=1,NDIM3
                      AB_edge = AB_edge +  &
                      AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO(NISM_edge),KSM_A(NISM_edge),DIF_A,APX)
                      CALL FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                          KSM_A(NISM_edge),KSM_AP(NISM_edge),SM_A,APX,FF_edge)
!**** Here we define the endge a priori term contribtuion for the right side of Noramal System: ****
!**** SMIM_edge - smoothness  edge term that should be added to general a priori term         ****
                      CALL smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge(NISM_edge,:),  &
                                                KSM_AP(NISM_edge),SM_AP,SMIM_edge)
                    ENDDO ! I2=1,NDIM3

                  ENDIF ! GSM .GT. 0. .AND.
                ENDDO ! IDIM2
                ENDDO ! IDIM1
                if(.not. lnism) NISM_edge = NISM_edge-1
              ENDIF ! NX_EDGE_AP(IT,IY) .GE. 1 .AND. 
            ENDDO ! DO IX=1,NY_edge
            ENDDO ! DO IX=1,NY_edge
!stop 'T_after edges test stop'
          end select ! I_EDGE  !***********************************************
          if(RIN%IPRI_verbose)  &
          write(iu_main_output,'(a,i2,a,i5,a)') 'I_EDGE=',I_EDGE,'  NISM_edge=',NISM_edge, &
          '     in smoothterm_mutli_pixel_edges '
ENDDO LOOP_I_EDGE    ! MAX of N_I_EDGE = 6

! Deallocate X_edge array
      if(allocated(X_edge)) then
        deallocate(X_edge,stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate X_edge(:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Deallocate Y_edge array
      if(allocated(Y_edge)) then
        deallocate(Y_edge,stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate Y_edge(:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Deallocate T_edge array
      if(allocated(T_edge)) then
        deallocate(T_edge,stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate T_edge(:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Deallocate IT array
      if(allocated(IT)) then
        deallocate(IT,stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate IT(:) in smoothterm_mutli_pixel_edges'
        endif              
      endif
! Deallocate AP_edge array
      if(allocated(AP_edge)) then
        deallocate(AP_edge,stat=alloc_stat)
        if (alloc_stat /= 0) then
        stop 'error while trying to deallocate AP_edge(:,:,:,:) in smoothterm_mutli_pixel_edges'
        endif              
      endif

!write(*,*) 'NISM_edge=',NISM_edge
goto 55
      if(NISM_edge .GT. 0) THEN
        IKSIM = IKSIM+IKSIM_E
        do IS_edge=1,NISM_edge
          NISM = NISM+1
!          IS = NISM
          KSM(NISM) = KSM_AP(NISM_edge)
          IMSM(NISM,1:KSM(NISM)) = IMSM_edge(IS_edge,1:KSM(NISM))
          do I=1,RIN%KNSINGF
            do IS1=1,KSM_AP(IS_edge)
!write(*,*) I,IS_edge,IS1,IMSM(IS_edge,IS1),'  - I,IS_edge,IS1,IMSM(IS_edge,IS1) in smoothterm_mutli_pixel_edges'
            do IS2=1,KSM_AP(IS_edge)
!write(*,*) I,IS_edge,IS2,IMSM(IS_edge,IS2),'  - I,IS_edge,IS2,IMSM(IS_edge,IS2) in smoothterm_mutli_pixel_edges'
              SMMULTI(IMSM(NISM,IS1),I,IMSM(NISM,IS2)) = SMMULTI(IMSM(NISM,IS1),I,IMSM(NISM,IS2))+ &
              SMIM_edge(IMSM_edge(IS_edge,IS1),I,IMSM_edge(IS_edge,IS2))
!write(*,*) 'IS_edge,IS1,IS2,IMSM(NISM,IS1),I,IMSM(NISM,IS2),SMMULTI(IMSM(NISM,IS1),I,IMSM(NISM,IS2)),', &
!'IMSM_edge(IS_edge,IS1),IMSM_edge(IS_edge,IS2) : ', &
!IS_edge,IS1,IS2,IMSM(NISM,IS1),I,IMSM(NISM,IS2),SMMULTI(IMSM(NISM,IS1),I,IMSM(NISM,IS2)), &
!IMSM_edge(IS_edge,IS1),IMSM_edge(IS_edge,IS2)
            enddo ! IS2
            enddo ! IS1      
          enddo ! I
        enddo ! IS_edge
      endif ! NISM_edge .GT. 0
55 continue

! OD 2015-09-16 Testing needed for IO=2 multi pixel constrains

      if(NISM_edge .GT. 0) THEN
        IKSIM = IKSIM+IKSIM_E
        do IS_edge=1,NISM_edge
          NISM = NISM+1
! bug ???          KSM(NISM) = KSM_AP(NISM_edge)
          KSM(NISM) = KSM_AP(IS_edge)
          IMSM(NISM,1:KSM(NISM)) = IMSM_edge(IS_edge,1:KSM(NISM))
        enddo ! IS_edge
        SMMULTI(:,:,:) = SMMULTI(:,:,:) + SMIM_edge(:,:,:)
      endif ! NISM_edge .GT. 0

      if(IPRINT2) then
      endif
!stop 'test stop in smoothterm_mutli_pixel_edges'

      END SUBROUTINE smoothterm_mutli_pixel_edges
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates Fisher matrix
        !>
        !> @param[in]    KM1 - number of columns of U1
        !> @param[in]    KM2 - number of columns of U2
        !> @param[in]    KN  - number of rows
        !> @param[in]    U1  - matrix
        !> @param[in]    U2  - matrix
        !> @param[out]   UFS - Fisher matrix 
        !>
      SUBROUTINE MAT_T_MAT ( KM1,KM2,KN,U1,U2,UFS )
!C**************************************************
!C  THIS SUBROUTINE CALCULATES "FISHER MATRIX":
!C            (U1)**T ((U2)
!C                and
!C!C**************************************************
!C  INPUT:
!C        KN  I         - number of lines
!C        KM1 I         - number of columns of U1
!C        KM2 I         - number of columns of U2
!C        U1   R(KM1,KN) -  matrix
!C        U2   R(KM2,KN) -  matrix 
!C  OUTPUT:
!C         UF R(KN,KN) - matrix(U1)**T ((U2)
!C***************************************************
      USE mod_par_inv, only : KDIF

      IMPLICIT NONE

! ----------------------------------------------------------
! IN :
      INTEGER,                       INTENT(IN)  :: KM1,KM2,KN
      REAL,DIMENSION(KDIF,KDIF),     INTENT(IN)  :: U1,U2
!      REAL,DIMENSION(KMESS),        INTENT(IN)  :: CS
!      REAL,DIMENSION(KMESS),        INTENT(IN)  :: FS,FPS
! ----------------------------------------------------------
! OUT :
      REAL,DIMENSION(KDIF,KDIF),     INTENT(OUT) :: UFS
! ----------------------------------------------------------
! LOCAL :	  
      INTEGER :: I,I1
      REAL :: AM
! ----------------------------------------------------------

!C*** calculating "Fisher matrix" 
!write(*,*) 'in MAT_T_MAT: KM1=',KM1,'  KM2=',KM2,'  KN=',KN
!write(*,*) 'in MAT_T_MAT: U1:'
!do i=1,KM1
!write(*,'(10e14.4)') U1(1:KN,i)
!enddo
!write(*,*) 'in MAT_T_MAT: U2:'
!do i=1,KM2
!write(*,'(10e14.4)') U2(1:KN,i)
!enddo
      DO I=1,KM1
       DO I1=1,KM2
          UFS(I,I1)=SUM(U1(1:KN,I)*U2(1:KN,I1)) 
       ENDDO ! I1
!write(*,*) 'in MAT_T_MAT 1: I=',I,'  UFS(I,1:KM2):'
!write(*,'(10e14.4)') UFS(I,1:KM2)
      ENDDO ! I

!      AM = 0.
!      DO I=1,KDIF
!        IF(AM .LT. abs(UFS(I,I))) AM = abs(UFS(I,I))
!      ENDDO ! I
!      IF(AM .NE. 0.) THEN
!        DO I =1,KDIF
!          DO I1=1,KDIF
!            UFS(I,I1) = UFS(I,I1)/AM
!          ENDDO ! I1
!        ENDDO ! I
!write(*,*) 'in MAT_T_MAT 12: I=',1,'  UFS=',UFS(1,1:1)
!      ENDIF ! AM.NE.0.

      RETURN
      END SUBROUTINE MAT_T_MAT

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Function returns 
        !>
        !> @param[in]    KIMAGE_E   - 
        !> @param[in]    GSM        - 
        !> @param[in]    I2         - 
        !> @param[in]    ISTARSING  - 
        !> @param[in]    NIO        - 
        !> @param[in]    KSM_A      -  
        !> @param[in]    DIF_A      -
        !> @param[in]    APX        -
        !>

      real function AB_edge_XYT(KIMAGE_E,GSM,I2,ISTARSING,NIO,KSM_A,DIF_A,APX)

      use mod_par_inv, only : KDIF,KPARS

      implicit none
! ----------------------------------------------------------
! IN :
      integer,                       intent(in)  :: KIMAGE_E,I2,ISTARSING
      integer,                       intent(in)  :: NIO,KSM_A
      real,                          intent(in)  :: GSM
      real,dimension(KDIF,KDIF),     intent(in)  :: DIF_A
      real,dimension(KPARS,KIMAGE_E),intent(in)  :: APX
! ----------------------------------------------------------
! LOCAL :
      integer :: I
      real    :: AB
! ----------------------------------------------------------
!OD&TL ???
      AB_edge_XYT = 0.0
      do I=1,NIO
        AB = DOT_PRODUCT(DIF_A(I,1:KSM_A),APX(ISTARSING+I2-1,1:KSM_A))
        AB_edge_XYT = AB_edge_XYT+GSM*AB*AB
      enddo ! I

      end function AB_edge_XYT

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates 
        !>
        !> @param[in]    KIMAGE_E   - 
        !> @param[in]    GSM        - 
        !> @param[in]    I2         - 
        !> @param[in]    ISTARSING  - 
        !> @param[in]    IMSM_edge  - 
        !> @param[in]    KSM_A      -  
        !> @param[in]    KSM_AP     -
        !> @param[in]    SM_A       -
        !> @param[in]    APX        -
        !> @param[inout] FF_edge    -
        !>
      subroutine FF_edge_calcul(KIMAGE_E,GSM,I2,ISTARSING,IMSM_edge,KSM_A,KSM_AP,SM_A,APX,FF_edge)

      use mod_par_inv, only : KDIF,KPARS,KMPSM,KIMAGE

      implicit none
! ----------------------------------------------------------
! IN :
      integer,                        intent(in)  :: I2,ISTARSING,KIMAGE_E
      integer,                        intent(in)  :: KSM_A,KSM_AP
      real,                           intent(in)  :: GSM
      real,dimension(KDIF,KDIF),      intent(in)  :: SM_A
      real,dimension(KPARS,KIMAGE_E), intent(in)  :: APX
      integer,dimension(KIMAGE),      intent(in)  :: IMSM_edge
      real,dimension(KPARS,KIMAGE),intent(inout)  :: FF_edge
! ----------------------------------------------------------
! LOCAL :
      integer :: I, ii
      real    :: FF
! ----------------------------------------------------------
      FF = 0.0
      do I=1,KSM_AP
        FF = GSM*DOT_PRODUCT(SM_A(I,1:KSM_A),APX(ISTARSING+I2-1,1:KSM_A))
!write(*,*) 'I=',I,'  IMSM_edge(I)=',IMSM_edge(I)
!do ii=1, KSM_A
!write(*,*) I,I2,ISTARSING+I2-1,GSM,SM_A(I,ii),exp(APX(ISTARSING+I2-1,ii)), &
!'  - I,I2,ISTARSING+I2-1,GSM,SM_A(I,ii),exp(APX(ISTARSING+I2-1,ii))'
!enddo ! ii
        FF_edge(ISTARSING+I2-1,IMSM_edge(I)) =  &
        FF_edge(ISTARSING+I2-1,IMSM_edge(I)) + FF
      enddo

      end subroutine FF_edge_calcul

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine calculates smoothness term for edges of inverted segment.
        !>
        !> @param[in]    GSM        - 
        !> @param[in]    I2         - 
        !> @param[in]    ISTARSING  - 
        !> @param[in]    IMSM_edge  - 
        !> @param[in]    KSM_AP     -
        !> @param[in]    SM_AP       -
        !> @param[inout] SMIM_edge    -
        !>
      subroutine smoothness_edge_term(GSM,I2,ISTARSING,IMSM_edge,KSM_AP,SM_AP,SMIM_edge)

      use mod_par_inv, only : KDIF,KPARS,KIMAGE

      implicit none
! ----------------------------------------------------------
! IN :
      integer,                        intent(in)  :: I2,ISTARSING
      integer,                        intent(in)  :: KSM_AP
      real,                           intent(in)  :: GSM
      real,dimension(KDIF,KDIF),      intent(in)  :: SM_AP
      integer,dimension(KIMAGE),      intent(in)  :: IMSM_edge
      real,dimension(KIMAGE,KPARS,KIMAGE),intent(inout)  :: SMIM_edge
! ----------------------------------------------------------
! LOCAL :
      integer :: I1,J1
! ----------------------------------------------------------
      do I1=1,KSM_AP
      do J1=1,KSM_AP
        SMIM_edge(IMSM_edge(I1),ISTARSING+I2-1,IMSM_edge(J1)) =  &
        SMIM_edge(IMSM_edge(I1),ISTARSING+I2-1,IMSM_edge(J1)) +  & 
        GSM*SM_AP(I1,J1)
!write(*,*) 'I1=',I1,'  J1=',J1,'  GSM=',GSM,'  SM_AP(I1,J1)=',SM_AP(I1,J1),  &
!           '  SMIM_edge(IMSM_edge(I1),ISTARSING+I2-1,IMSM_edge(J1))=',       &
!              SMIM_edge(IMSM_edge(I1),ISTARSING+I2-1,IMSM_edge(J1))

      enddo ! J1
      enddo ! I1

      return
      end subroutine smoothness_edge_term

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates 
        !>
        !> @param[in]    NSUM   - 
        !> @param[in]    N_EDGE_A        - 
        !> @param[in]    N_EDGE_AP         - 
        !> @param[inout] IOI  - 
        !> @param[out]   KSM_A      -  
        !> @param[out]   KSM_AP     -
        !>
      subroutine KSM_calcul(NSUM,N_EDGE_A,N_EDGE_AP,IOI,KSM_A,KSM_AP)

      implicit none
! ----------------------------------------------------------
! IN :
      integer,intent(in)     :: NSUM
      integer,intent(in)     :: N_EDGE_A,N_EDGE_AP
      integer,intent(inout)  :: IOI 
      integer,intent(out)    :: KSM_A,KSM_AP
! ----------------------------------------------------------
! ----------------------------------------------------------

      if(NSUM .lt. IOI+1) IOI = NSUM-1

      if(N_EDGE_A .gt. IOI) then
        KSM_A = IOI
      else
        KSM_A = N_EDGE_A
      endif

      if(N_EDGE_AP .gt. IOI) then
        KSM_AP = IOI
      else
        KSM_AP = N_EDGE_AP
      endif

      return
      end subroutine KSM_calcul 
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates 
        !>
        !> @param[in]    NIO    - 
        !> @param[in]    KSM_A  - 
        !> @param[in]    NNN_E  - 
        !> @param[in]    DIF    - 
        !> @param[in]    DIF_A  - 
        !> @param[in]    DIF_AP -  
        !>
      subroutine DIF_before_calcul( NIO, KSM_A, NNN_E, DIF, DIF_A, DIF_AP )

      use mod_par_inv, only : KDIF
      
      implicit none
! ----------------------------------------------------------
! IN :
      integer,intent(in)                    :: NIO,KSM_A,NNN_E
      real,dimension(KDIF,KDIF),intent(in)  :: DIF
      real,dimension(KDIF,KDIF),intent(out) :: DIF_A,DIF_AP
! ----------------------------------------------------------
! LOCAL:
      real,dimension(KDIF,KDIF) :: SM1
      real,dimension(KDIF,KDIF) :: DIF_temp
      integer  :: I,I1,J
      real :: xnorm
! ----------------------------------------------------------
        DIF_temp(:,:) = DIF(:,:)
        xnorm = 0.0 ! smoothnes matrix is normalized

        CALL SMOOM ( KDIF,KDIF, & ! IN
                     DIF(:,:), &
                     SM1(:,:), xnorm  & ! OUT
                  )
        IF(xnorm .NE. 0.0) THEN
          DO I =1,KDIF
          DO I1=1,KDIF
            DIF_temp(I,I1) = DIF_temp(I,I1)/sqrt(xnorm)
          ENDDO ! I1
!write(*,*) 'in DIF_before_calcul: I=',I,'  DIF_temp=',DIF_temp(I,1:1)
          ENDDO ! I
        ENDIF ! xnorm .NE. 0.0

      DIF_A(:,:)  = 0.0
      do I=1,NIO
      do J=1,KSM_A
        DIF_A(I,J) = DIF_temp(I,J)
!write(*,*) 'in DIF_before_calcul: I=',I,'  J=',J,'  DIF_A(I,J)=',DIF_A(I,J)
      enddo
      enddo

      DIF_AP(:,:) = 0.0
      do I=1,NIO
      do J=KSM_A+1,NNN_E
        DIF_AP(I,J-KSM_A) = DIF_temp(I,J)
!write(*,*) 'in DIF_before_calcul: I=',I,'  J=',J,'  J-KSM_A(NISM_edge)=',J-KSM_A,'  KSM_A=',KSM_A,  &
!'  DIF_AP(I,J-KSM_A)=',DIF_AP(I,J-KSM_A)
      enddo
      enddo

      return
      end subroutine DIF_before_calcul 
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !!> @brief Routine calculates 
        !>
        !> @param[in]    NIO    - 
        !> @param[in]    KSM_AP  - 
        !> @param[in]    NNN_E  - 
        !> @param[in]    DIF    - 
        !> @param[in]    DIF_A  - 
        !> @param[in]    DIF_AP -  
        !>
      subroutine DIF_after_calcul( NIO, KSM_AP, NNN_E, DIF, DIF_A, DIF_AP )

      use mod_par_inv, only : KDIF
      
      implicit none
! ----------------------------------------------------------
! IN :
      integer,intent(in)                    :: NIO,KSM_AP,NNN_E
      real,dimension(KDIF,KDIF),intent(in)  :: DIF
      real,dimension(KDIF,KDIF),intent(out) :: DIF_A,DIF_AP
! ----------------------------------------------------------
! LOCAL:
      real,dimension(KDIF,KDIF) :: SM1
      real,dimension(KDIF,KDIF) :: DIF_temp
      integer  :: I,I1,J
      real :: xnorm
! ----------------------------------------------------------
        DIF_temp(:,:) = DIF(:,:)
        xnorm = 0.0  ! smoothnes matrix is normalized

        CALL SMOOM ( KDIF,KDIF, & ! IN
                     DIF(:,:), &
                     SM1(:,:), xnorm & ! OUT
                    )
        IF(xnorm .NE. 0.0) THEN
          DO I =1,KDIF
          DO I1=1,KDIF
            DIF_temp(I,I1) = DIF_temp(I,I1)/sqrt(xnorm)
          ENDDO ! I1
!write(*,*) 'in DIF_after_calcul: I=',I,'  DIF_temp=',DIF_temp(I,1:KN)
          ENDDO ! I
        ENDIF ! xnorm .NE. 0.0

      DIF_AP(:,:) = 0.0
      do I=1,NIO
      do J=1,KSM_AP
        DIF_AP(I,J) = DIF_temp(I,J)
!write(*,*) 'in DIF_after_calcul: I=',I,'  J=',J,'  DIF_AP(I,J)=',DIF_AP(I,J)
      enddo
      enddo

      DIF_A(:,:)  = 0.0
      do I=1,NIO
      do J=KSM_AP+1,NNN_E
        DIF_A(I,J-KSM_AP) = DIF_temp(I,J)
!write(*,*) 'in DIF_after_calcul: I=',I,'  J=',J,'  J-KSM_AP=',J-KSM_AP,'  KSM_AP=',KSM_AP,  & 
!'  DIF_A(I,J-KSM_AP)=',DIF_A(I,J-KSM_AP)
      enddo
      enddo

      return
      end subroutine DIF_after_calcul 
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
