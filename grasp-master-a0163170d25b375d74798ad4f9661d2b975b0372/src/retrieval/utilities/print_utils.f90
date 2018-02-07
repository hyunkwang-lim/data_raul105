! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains :

! subroutine print_input_settings
! subroutine write_rad_transf_input
! subroutine iP_write_output_retrieval
! subroutine print_classic_plot
! subroutine print_phmx

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine print_input_settings ( RIN )

      use mod_par_inv
      use mod_retr_settings_derived_type
	      
      implicit none
! --------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) ::  RIN
! --------------------------------------------------------------------------------------
! LOCAL :
      integer                         ::  IDIM1,IDIM2,IDIM3
      character(len=255)              ::  file_name
      integer                         ::  iu_iP_input
!      real,dimension(KIDIM3,KIDIM2)   ::  RMIN,RMAX
!      real                            ::  AD
      integer                         ::  I,II,par_type !,IW,keyBIN
      character(len=20)               ::  mpl
! ---------------------------------------------------------------------------------------
	  	  
         file_name = trim(RIN%DLSF%external_file_path)//'INPUT_Settings_test.txt'
         OPEN (newunit=iu_iP_input, FILE=trim(file_name), STATUS='replace',recl=20000)
	        
!C*****************************************************************************************

      WRITE (iu_iP_input,'(i4,a)') RIN%DLSF%keyEL,'  keyEL' 
      
      WRITE (iu_iP_input,'(i4,3L3,a)') RIN%KL,RIN%ISTOP,RIN%IPRI_additional_info,RIN%IPRI_iter_result, &
                      '     - KL,ISTOP,IPRI_additional_info,IPRI_iter_result'
      WRITE (iu_iP_input,'(3i4,f7.2,2i4,l6,e12.4,i4,a)') &
      RIN%NSD,RIN%NLYRS,RIN%ipplane,RIN%SHIFT,   &
      RIN%iPOBS,RIN%iIOBS,                       &
      RIN%ITRONC,RIN%eps_err,RIN%mol_prof_type,  &
      '     - NSD,NLYRS,ipplane,SHIFT,iPOBS,iIOBS,ITRONC,eps_err,mol_prof_type'

!C*****************************************************************************************

      WRITE (iu_iP_input,'(8i4,a)') RIN%OSHF%IMSC,RIN%OSHF%NG,RIN%OSHF%NN,RIN%OSHF%NF,   &
                            RIN%OSHD%IMSC,RIN%OSHD%NG,RIN%OSHD%NN,RIN%OSHD%NF,   &
                            '     - IMSC_F,NG_F,NN_F,NF_F, IMSC_D,NG_D,NN_D,NF_D' 

!C*****************************************************************************************
      WRITE(mpl,*) RIN%NDIM%n1
      WRITE (iu_iP_input,'(i4,'//trim(mpl)//'(i8,l3),a)') RIN%NDIM%n1,  & 
      (RIN%NDIM%par_type(I),RIN%NDIM%par_retr(I),I=1,RIN%NDIM%n1),'  - NDIM1,par_type,par_retr'
	  
      DO IDIM1=1,RIN%NDIM%n1
      WRITE (iu_iP_input,'(i4,a)') RIN%NDIM%n2(IDIM1),'     - NDIM2'
      WRITE(mpl,*) RIN%NDIM%n2(IDIM1)
      WRITE (iu_iP_input,'('//trim(mpl)//'i4,a)') (RIN%NDIM%n3(IDIM2,IDIM1),IDIM2=1,RIN%NDIM%n2(IDIM1)),'     - NDIM3'
      ENDDO ! IDIM1

!C*****************************************************************************************
 
      WRITE (iu_iP_input,'(i4,a)')   RIN%NOISE%INOISE,'     - INOISE'		
      DO I=1,RIN%NOISE%INOISE
        WRITE (iu_iP_input,*)  &
        RIN%NOISE%SGMS(I),RIN%NOISE%INN(I),RIN%NOISE%DNN(I) ,               &
        RIN%NOISE%NMT(I),(RIN%NOISE%MT(II,I),RIN%NOISE%NWLP(II,I),          & 
        RIN%NOISE%IWLP(1:RIN%NOISE%NWLP(II,I),II,I),II=1,RIN%NOISE%NMT(I)), &
        '     - SGMS,INN,DNN, NMT,(MT(IP,I),NWLP(IP,I),IWLP(1:NWLP(IP,I),IP,1),IP=1,NMT(I))'
      ENDDO ! I

!C*****************************************************************************************

      WRITE (iu_iP_input,*)   RIN%IBIN,'     - IBIN'
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .gt. par_type_SD_beg  .and. par_type .lt. par_type_SD_end) then      
          if(par_type .eq. par_type_SD_TB) then
            do IDIM2=1,RIN%NDIM%n2(IDIM1)
              write(iu_iP_input,'(2e13.5,a)') RIN%RMIN(IDIM2),RIN%RMAX(IDIM2),'     - RMIN,RMAX'
            enddo ! IDIM2
          else 
            do IDIM2=1,RIN%NDIM%n2(IDIM1)
              write(iu_iP_input,'(2e13.5,3x,22e13.5)')                          &
              RIN%RMIN(IDIM2),RIN%RMAX(IDIM2), &
              (RIN%RADIUS1(IDIM3,IDIM2),IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1))  
            enddo ! IDIM2
          endif ! par_type .gt. par_type_SD_TB
        exit
        endif ! par_type .gt. par_type_SD_beg  .and. 
      enddo ! IDIM1

!C*****************************************************************************************
      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then
          if(RIN%NDIM%par_type(IDIM1) .eq. par_type_SHD_distr) then
            write(mpl,*) RIN%NDIM%n3(IDIM2,IDIM1)
            do IDIM2=1,RIN%NDIM%n2(IDIM1) 
              write(iu_iP_input,'('//trim(mpl)//'e12.4,a)')  & 
              (RIN%RATIO1(IDIM3,IDIM2),IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)),'    - RATIO1(:,:)'
            enddo  ! IDIM2   
          endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SHD_SHD_distr
        exit
        endif ! par_type .gt. par_type_SHD_beg  .and. 
      enddo ! IDIM1      
!C*****************************************************************************************
      WRITE (iu_iP_input,'(a)') 'PARAMETERS OF MATRIX INVERSION (or Q-iterations):'

      WRITE (iu_iP_input,'(2i4,L4,3e12.4,a)')              &
          RIN%IMQ,RIN%IPFP%INVSING,RIN%IPFP%TXY_group,     &
          RIN%IPFP%TTEST,RIN%IPFP%XTEST,RIN%IPFP%YTEST,    &   
          '     - IMQ,INVSING,TXY_group,TTEST,XTEST,YTEST'
!C*****************************************************************************************      
      WRITE(iu_iP_input,'(a)') 'PARAMETER FOR Levenb.-Marquardt correction:'
      WRITE(iu_iP_input,'(i4,a)')  RIN%IPSTOP,'     - IPSTOP'
!C*****************************************************************************************
      WRITE(iu_iP_input,'(a)') 'SINGLE PIXEL A PRIORI ESTIMATES:'
      DO IDIM1=1,RIN%NDIM%n1
       DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
         DO IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
           WRITE(iu_iP_input,'(i4,e12.4,3i4,a)') &
           RIN%SPCA%IO(IDIM3,IDIM2,IDIM1),RIN%SPCA%GSM(IDIM3,IDIM2,IDIM1), &
           IDIM3,IDIM2,IDIM1,'     - IO,GSM, IDIM3,IDIM2,IDIM1'

         ENDDO
        ENDDO
      ENDDO
!C*****************************************************************************************
      WRITE(iu_iP_input,'(a)') 'SINGLE PIXEL SMOOTHNESS ESTIMATES:'
      DO IDIM1=1,RIN%NDIM%n1
       DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
         WRITE(iu_iP_input,'(i4,e12.4,2i4,a)') &
         RIN%SPCS%IO(IDIM2,IDIM1),RIN%SPCS%GSM(IDIM2,IDIM1), &
         IDIM2,IDIM1,'     - IO,GSM, IDIM2,IDIM1'
       ENDDO ! IDIM2   
      ENDDO ! IDIM1

      WRITE(iu_iP_input,'(a)') 'MULTI-pixes smoothness constraints:'
      DO IDIM1=1,RIN%NDIM%n1
       DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
         WRITE(iu_iP_input,'(3(i4,e12.4),2i4,a)') &
         RIN%MPCS%IOT(IDIM2,IDIM1),RIN%MPCS%GSMT(IDIM2,IDIM1), &
		     RIN%MPCS%IOX(IDIM2,IDIM1),RIN%MPCS%GSMX(IDIM2,IDIM1), &
         RIN%MPCS%IOY(IDIM2,IDIM1),RIN%MPCS%GSMY(IDIM2,IDIM1), &
         IDIM2,IDIM1,'     - IOT,GSMT,IOX,GSMX,IOY,GSMY, IDIM2,IDIM1'
        
       ENDDO   
      ENDDO
!C*****************************************************************************************

      WRITE (iu_iP_input,'(i5,3e12.4,a)') RIN%MAXP,RIN%EPSP,RIN%EPSQ,RIN%DL,'     - MAXP,EPSP,EPSQ,DL'

!C*****************************************************************************************
      WRITE (iu_iP_input,'(a)') 'INITIAL GUESS  and "MIN" and "MAX" values FOR THE SOLUTION:'	  
      DO IDIM1=1,RIN%NDIM%n1
       DO IDIM2=1,RIN%NDIM%n2(IDIM1) 
        DO IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
         II= RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1
         WRITE(iu_iP_input,'(3e12.4,i4,a)')                               &
         RIN%APSING(II),RIN%APSMIN(II),RIN%APSMAX(II),RIN%IWW_SINGL(II),  &
         '     - APSING,APSMIN,APSMAX,IWW_SINGL'
        ENDDO
       ENDDO   
      ENDDO
!C*****************************************************************************************
      write (iu_iP_input,'(a)') 'SIZES OF SEGMENT EDGES:'
      write(iu_iP_input,'(3i5,a)') RIN%edges%nx,RIN%edges%ny,RIN%edges%nt, &
                                   '  - edges%nx,edges%ny,edges%nt'     
!C*****************************************************************************************
      write(iu_iP_input,'(a)') 'GENERAL OUTPUT products'
      ! Retrieval output
      write(iu_iP_input,'(3L4,a)') RIN%products%retrieval,'  - retrieval'  
      ! Aerosol output
      write(iu_iP_input,'(7L4,a)') RIN%products%aerosol,'  - aerosol'  
      ! Clouds output
      write(iu_iP_input,'(7L4,a)') RIN%products%clouds,'  - clouds'
      ! Surface output
      write(iu_iP_input,'(1L4,a)') RIN%products%surface,'  - surface'  
      ! Forcing output
      write(iu_iP_input,'(2L4,a)') RIN%products%forcing,'  - forcing'  
      ! Error estimates output
      write(iu_iP_input,'(5L4,a)') RIN%products%errest,'  - errest'  
!C*****************************************************************************************
      WRITE(iu_iP_input,'(1L4,4i4,a)')  RIN%INPUT,RIN%Aexp_iwl(1:2),RIN%ndvi_iwl(1:2),  &
                                       '     - INPUT, Aexp_iwl(1:2), ndvi_iwl(1:2)'
!C*****************************************************************************************

      WRITE(iu_iP_input,'(2a)') RIN%DLSF%internal_file_path,'     - internal_file_path'
      WRITE(iu_iP_input,'(2a)') RIN%DLSF%external_file_path,'     - external_file_path'

      WRITE(iu_iP_input,'(2a)') RIN%DLSF%distname_O,'     - distname_O'
      WRITE(iu_iP_input,'(2a)') RIN%DLSF%distname_N,'     - distname_N'

!C*****************************************************************************************

      close(iu_iP_input)
      
      return
      end subroutine print_input_settings

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE write_rad_transf_input ( external_file_path, IPRI_additional_info, &
                                          IMSC, NG1, NG2, NF, eps_err,              &
                                          iBRDF_land, iBPDF_land, iBRM_water, land_percent, &
                                          NQDR, ITRONC, iop, NT1,                   &
                                          sza, NBV, vis, fis,                       &
                                          IW, WAVE,                                 &
                                          surf_land_par_num, surf_land_par_vect,    &
                                          surf_water_par_num, surf_water_par_vect,  &
                                          !NM, ngas,					                        &
                                          EXT, SSA,                                 &
                                          NANG, ANGL,                               &
                                          HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,    &
                                          !KVERTM, NH0, H01, PROF0,                  &
                                          laerosol, lsurface,           	          &
                                          !mol_prof_type, aer_prof_type,             &
                                          PF11, PF12, PF22, PF33                    &
                                        )


!C********************************************************************
!C*** This subroutine writes INPUT for OS_HERMAN subroutine into file
!C*** if IPRI=1 (INPUT check)
!C********************************************************************
!C
!C   INPUT:
!C     IMSC  I    -defines the order of scattering:
!C             = 0 -multiple scattering
!C             = 1 -single scattering
!C     HOBS
!C     IW    I - order number of the spectral channel
!C     WAVE  R - value of wavelength
!C     HRG   R - hight of the the surface over the see level
!C     SZA   R - solar zenith angle
!C     NT1-1 I - number of layers in radiance transfer calculations
!C     NM    I - number of aerosol fractions
!C     NANG  i - number of angles for S-matrix
!C     ANGL  R(NANG) - values of scattering angles for phase function
!C     EXT   R(1:KM+1) - EXT for each aerosol component + molecular 
!C     SSA   R(1:KM+1) - SSA for each aerosol component + molecular
!C     JP    I(1:KM+1)
!C     ECH   R(1:KM+1)
!C     sigma R(1:KM+1)
!C     PFII  R(1:NANG,1:1:KM) - scattering matrix for each aerosol component
!C
!C     NG1
!C     NG2
!C     NF
!C     iBRDF
!C     iBPDF

!C     NBF  I - number of parameters used for modeling BRF
!C     NBRP I - number of parameters used for modeling polarized BRF
!C     BRF   R(KBRF) - parameters  for modeling  BRF
!C     BRP  R(KBRF) - parameters  for modeling polarized BRF
!C     H0    R(KSD)  - median height in gaussian profile of vertical distribution
!C                     of each aerosol component
!C     W     R(KSD)  - width of gaussian profile of vertical distribution
!C                     of each aerosol component
!C     iop I - determines how we use polarization
!C             = 0 - in scattering plane
!C             = 1 - in meridian (vertical) plane

!C     NBV
!C     vis   R(NBVM) - zenith observation angles
!C     fis   R(NBVM) - azimuth observation angles
!C     ITRONC
!C     NQDR
!C     isaut,wind,iwat,rsurf,igrd
!C     vk,gteta,rho,rhoH,anb

!C********************************************************************

      use mod_par_OS,  only : NBVM, NMM, NMG, KSD
      use mod_par_inv, only : KBF
      use mod_vertical_distr_derived_type

	    implicit none
!	-----------------------------------------------------------------------------------
! IN :
      character(*),               intent(in)  ::  external_file_path
      logical,                    intent(in)  ::  IPRI_additional_info
      real,dimension(2*NBVM),     intent(in)  ::  vis,fis
      real,dimension(NMM+NMG),    intent(in)  ::  EXT, SSA ! ,ECH, sigma
      !integer,dimension(NMM+NMG), intent(in)  ::  JP
      !integer,                    intent(in)  ::  KVERTM
      !integer,dimension(NMM+NMG), intent(in)  ::  NH0
      !real,dimension(KVERTM,NMM+NMG),intent(in)  ::  H01
      !real,dimension(KVERTM,NMM+NMG),intent(in)  ::  PROF0
      type(discret_vertical_distribution), intent(in) :: DISCRVD
      integer,                    INTENT(IN)  ::  NBV, ITRONC, iop,       &
	                                                NQDR, NT1, NANG, IMSC,  &
                                                  NG1, NG2, NF
      real,                       intent(in)  ::  eps_err
      real,                       intent(in)  ::  sza, WAVE, HOBS_km, HGR_km, HMAX_atm_km

      real,dimension(NANG),       intent(in)  ::  ANGL
      real,dimension(NANG,KSD),   intent(in)  ::  PF11, PF12, PF22, PF33
      integer,                    intent(in)  ::  IW, iBRDF_land, iBPDF_land, iBRM_water
      logical,                    intent(in)  ::  laerosol, lsurface
      !integer,                    intent(in)  ::  mol_prof_type, aer_prof_type
      real,                       intent(in)  ::  land_percent
      integer,dimension(2),       intent(in)  ::  surf_land_par_num
      real, dimension(2*KBF),     intent(in)  ::  surf_land_par_vect
      integer,                    intent(in)  ::  surf_water_par_num
      real, dimension(KBF),       intent(in)  ::  surf_water_par_vect
!	-----------------------------------------------------------------------------------
      integer :: natm, naer, ncld, nmol, ngas
      integer        ::  IA, IV, ISD, iu_check_OS
      integer, save  ::  counter=1
!	-----------------------------------------------------------------------------------
!	-----------------------------------------------------------------------------------

!c***************************************************************
!      DIMENSION  FWAVE(122), WFBEAM(122)
!OD:  DO NOT DELETE DATA FWAVE and WFBEAM ( data can be used)
!c***************************************************************
!      DATA FWAVE /0.30, 0.31, 0.31, 0.32, 0.32, 0.33,
!     & 0.33, 0.34, 0.34, 0.35, 0.35, 0.36, 0.37, 0.38, 0.39,
!     & 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48,
!     & 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.57, 0.59,
!     & 0.61, 0.63, 0.66, 0.67, 0.69, 0.71, 0.72, 0.72, 0.74,
!     & 0.75, 0.76, 0.76, 0.77, 0.78, 0.80, 0.82, 0.82, 0.83,
!     & 0.84, 0.86, 0.88, 0.91, 0.92, 0.93, 0.93, 0.94, 0.95,
!     & 0.97, 0.98, 0.99, 1.04, 1.07, 1.10, 1.12, 1.13, 1.15,
!     & 1.16, 1.17, 1.20, 1.24, 1.27, 1.29, 1.32, 1.35, 1.40,
!     & 1.44, 1.46, 1.48, 1.50, 1.52, 1.54, 1.56, 1.58, 1.59,
!     & 1.61, 1.63, 1.65, 1.68, 1.74, 1.80, 1.86, 1.92, 1.96,
!     & 1.99, 2.01, 2.04, 2.07, 2.10, 2.15, 2.20, 2.27, 2.36,
!     & 2.45, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20,
!     & 3.30, 3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00/ 

!      DATA WFBEAM /535.90,	558.30,	622.00,	692.70,
!     & 715.10, 832.90, 961.90, 931.90, 900.60, 911.30, 
!     & 975.50, 975.90, 1119.90,	1103.80, 1033.80, 1479.10,	
!     & 1701.30,	1740.40, 1587.20, 1837.00, 2005.00,	2043.00,
!     & 1987.00,	2027.00, 1896.00, 1909.00, 1927.00, 1831.00,
!     & 1891.00,	1898.00, 1892.00, 1840.00, 1768.00,	1728.00,
!     & 1658.00,	1524.00, 1531.00, 1420.00, 1399.00,	1374.00,
!     & 1373.00,	1298.00, 1269.00, 1245.00, 1223.00,	1205.00,
!     & 1183.00,	1148.00, 1091.00, 1062.00, 1038.00,	1022.00,
!     & 998.70, 947.20, 893.20, 868.20, 829.70, 830.30, 
!     & 814.00, 786.90, 768.30, 767.00, 757.60, 688.10,
!     & 640.70, 606.20, 585.90, 570.20, 564.10, 544.20,
!     & 533.40, 501.60, 477.50, 442.70, 440.00, 416.80,
!     & 391.40, 358.90, 327.50, 317.50, 307.30, 300.40, 
!     & 292.80, 275.50, 272.10, 259.30, 246.90, 244.00,	
!     & 243.50, 234.80, 220.00, 190.80, 171.10, 144.50,
!     & 135.70, 123.00, 123.80, 113.00, 108.50, 97.50,
!     & 92.40, 82.40, 74.60,	68.30, 63.80, 49.50, 48.50,
!     & 38.60, 36.60, 32.00, 28.10, 24.80, 22.10, 19.60,
!     & 17.50, 15.70, 14.10, 12.70, 11.50, 10.40, 9.50, 8.60/
!c***************************************************************
      natm = DISCRVD%natm
      naer = DISCRVD%naer
      ncld = DISCRVD%ncld
      nmol = DISCRVD%nmol
      ngas = DISCRVD%ngas

!       IF(counter .le. 10) THEN
!        IF(IW .eq. 1 .and. counter .eq. 1) THEN
! Print RT Inputs for 10 last RT calls
      iu_check_OS = 22
      IF(counter .eq. 1 .or. (counter .ge. 10 .and. MOD(counter,10) .eq. 0)) THEN
      !IF(counter .eq. 1) THEN
        open(iu_check_OS,file=trim(external_file_path)//'check_OS_input.dat',status='replace',recl=20000)
      ELSE
        open(iu_check_OS,file=trim(external_file_path)//'check_OS_input.dat',status='old',      &
                                          position='append',recl=20000)	   
      ENDIF
      counter=counter+1
      WRITE(iu_check_OS,*) '**********************  counter=',counter-1

      WRITE(iu_check_OS,*) WAVE, "   WAVE"
      WRITE(iu_check_OS,*) IMSC, "   IMSC"
      WRITE(iu_check_OS,*) NBV,  "   NBV"

      DO IV=1,NBV
        WRITE(iu_check_OS,'(2f11.5,i5,a23)') vis(IV),fis(IV),IV,     &
        '   - vis(IV),fis(IV),IV'
      ENDDO
      WRITE(iu_check_OS,*) sza,' solar_zenith_angl_THETAS'
      WRITE(iu_check_OS,*)  surf_land_par_num(1:2),'   surf_par_num(1:2)'
      WRITE(iu_check_OS,*)  surf_land_par_vect(1:surf_land_par_num(1)), &
                             '   surf_land_par_vect(1:npar)'
      WRITE(iu_check_OS,*)  surf_water_par_num,"   surf_water_par_num"
      WRITE(iu_check_OS,*)  surf_water_par_vect(1:surf_water_par_num), &
                             '   surf_water_par_vect'

      WRITE(iu_check_OS,*)  iop,   '   iop'
      WRITE(iu_check_OS,*)  NG1,   '   NG1'
      WRITE(iu_check_OS,*)  NG2,   '   NG2'
      WRITE(iu_check_OS,*)  NF,    '   NF1'
      WRITE(iu_check_OS,*)  iBRDF_land,  '   iBRDF_land'
      WRITE(iu_check_OS,*)  iBPDF_land,  '   iBPDF_land'
      WRITE(iu_check_OS,*)  iBRM_water,  '   iBRM_water'
      WRITE(iu_check_OS,*)  land_percent,'   land_percent'
              
      WRITE(iu_check_OS,*)  naer,     '   naer'
      WRITE(iu_check_OS,*)  NANG,   '   NANG'
      WRITE(iu_check_OS,*)  NQDR,   '   NQDR'
      WRITE(iu_check_OS,*)  ITRONC, '   ITRONC'
      WRITE(iu_check_OS,*)  eps_err, '   eps_err'

      WRITE(iu_check_OS,'(5(i0,2x),a)') natm, naer, ncld, nmol, ngas, &
                                       '  natm, naer, ncld, nmol, ngas'
      WRITE(iu_check_OS,'(a)') 'EXT(1:natm)'
      WRITE(iu_check_OS,'(10e16.7)') EXT(1:natm)
      WRITE(iu_check_OS,'(a)') 'SSA(1:natm)'
      WRITE(iu_check_OS,'(10e16.7)') SSA(1:natm)

      WRITE(iu_check_OS,'(a)') 'iv, h_km(iv), VD(iv,1:natm)'
      DO iv=1,DISCRVD%nh
        WRITE(iu_check_OS,'(x,i02x,f16.7,10e16.7,a)') &
                              iv,DISCRVD%h_km(iv),DISCRVD%val(iv,1:natm)
      ENDDO ! iv
      DO ISD=1,naer
        DO IA=1,NANG
          WRITE(iu_check_OS,'(f9.3,4e18.8,a)')                           &
          ANGL(IA),PF11(IA,ISD),PF12(IA,ISD),PF22(IA,ISD),PF33(IA,ISD),  &
          '   - ANGL(IA),PF11(IA,ISD),PF12(IA,ISD),PF22(IA,ISD),PF33(IA,ISD)'
        ENDDO ! IA
       ENDDO ! ISD

       WRITE(iu_check_OS,*) NT1,'   NT1'
       WRITE(iu_check_OS,*) HOBS_km,HGR_km,HMAX_atm_km,  &
                         '  HOBS_km,HGR_km,HMAX_atm_km,'
       WRITE(iu_check_OS,*) laerosol,lsurface,'  laerosol,lsurface'

       close(iu_check_OS)
	  
!      ENDIF ! counter .le. 10
	   
      END SUBROUTINE write_rad_transf_input
	  
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! This routine finds and returns the indices i1, i2, i3 of a retrieved
! parameter (identified by the index ipar in the GOUT_retrieval_par structure)
! after unpacking into the usual multidimensional array.
subroutine get_i1_i2_i3_from_ipar(RIN, ipar, i1, i2, i3)
  use mod_retr_settings_derived_type
  use mod_par_inv, only : KPARS
  implicit none

  ! input and output
  type(retr_input_settings),   intent(in)  :: RIN
  integer,                     intent(in)  :: ipar
  integer,                     intent(out) :: i1, i2, i3
  ! local variables
  integer                         :: n1, n2, n3, i3offs
  integer, dimension(KPARS), save :: i1_cache = 0, i2_cache = 0, i3_cache = 0

  ! Take indices from cache, if present, and return immediately.
  if (i1_cache(ipar) .ne. 0) then
    i1 = i1_cache(ipar)
    i2 = i2_cache(ipar)
    i3 = i3_cache(ipar)
  end if

  ! Otherwise we find the corresponding indices,
  ! save them in the cache, and return them.

  n1 = RIN%NDIM%n1
  do i1 = 1, n1
    n2 = RIN%NDIM%n2(i1)
    do i2=1,n2
      n3 = RIN%NDIM%n3(i2,i1)
      i3offs = RIN%NDIM%ISTARSING(i2,i1) - 1
      do i3 = 1, n3
        !par(i3,i2,:) = GOUT_retrieval_par%pixel(:)%par(i3offs + i3)
        if (ipar .eq. i3offs + i3) then
          i1_cache(ipar) = i1
          i2_cache(ipar) = i2
          i3_cache(ipar) = i3
          return
        end if
      enddo
    enddo
  enddo

  ! We should not arrive here, since we expect to find the indices for ipar.
  stop 'stop in get_i1_i2_i3_from_ipar: invalid index ipar'
end subroutine

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine par_type_to_string(RIN, par_type, par_type_str)
  use mod_retr_settings_derived_type
  implicit none

  ! input and output
  type(retr_input_settings),   intent(in)  :: RIN
  integer,                     intent(in)  :: par_type
  character(LEN=50),           intent(out) :: par_type_str

  select case (par_type)
    case(par_type_SD_TB)
      par_type_str = 'SD_TB'
    case(par_type_SD_LB)
      par_type_str = 'SD_LB'
    case(par_type_SD_LN)
      par_type_str = 'SD_LN'
    case(par_type_RERI_spect)
      par_type_str = 'RERI_spect'
    case(par_type_RERI_const)
      par_type_str = 'RERI_const'
    case(par_type_CXRI_nmix)
      par_type_str = 'CXRI_nmix'
    case(par_type_CXRI_chem)
      par_type_str = 'CXRI_chem'
    case(par_type_IMRI_spect)
      par_type_str = 'IMRI_spect'
    case(par_type_IMRI_const)
      par_type_str = 'IMRI_const'
    case(par_type_SHD_fsph)
      par_type_str = 'SHD_fsph'
    case(par_type_SHD_distr)
      par_type_str = 'SHD_distr'
    case(par_type_AVP_par_height)
      par_type_str = 'AVP_par_height'
    case(par_type_AVP_prof)
      par_type_str = 'AVP_prof'
    case(par_type_Cv)
      par_type_str = 'Cv'
    case(par_type_CL)
      par_type_str = 'CL'
    case(par_type_AVP_par_std)
      par_type_str = 'AVP_par_std'
    case(par_type_SURF1_land_Ross_Li_BRDF)
      par_type_str = 'SURF1_land_Ross_Li_BRDF'
    case(par_type_SURF1_land_RPV_BRDF)
      par_type_str = 'SURF1_land_RPV_BRDF'
    case(par_type_SURF1_land_Litvinov)
      par_type_str = 'SURF1_land_Litvinov'
    case(par_type_SURF1_land_Litvinov_fast)
      par_type_str = 'SURF1_land_Litvinov_fast'
    case(par_type_SURF2_land_Maignan_Breon)
      par_type_str = 'SURF2_land_Maignan_Breon'
    case(par_type_SURF2_land_Litvinov)
      par_type_str = 'SURF2_land_Litvinov'
    case(par_type_SURF_water_Cox_Munk_iso)
      par_type_str = 'SURF_water_Cox_Munk_iso'
    case(par_type_SURF_water_Cox_Munk_ani)
      par_type_str = 'SURF_water_Cox_Munk_ani'
    case(par_type_SURF_water_Litvinov)
      par_type_str = 'SURF_water_Litvinov'
    case(par_type_SD_cloud_par)
      par_type_str = 'SD_cloud_par'
    case(par_type_CVP_par)
      par_type_str = 'CVP_par'
    case(par_type_fcloud_par)
      par_type_str = 'fcloud_par'
    case default
      write(par_type_str,'(i15)') par_type
  end select

end subroutine

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine i3_to_label(RIN, par_type, i3, i3label)
  use mod_retr_settings_derived_type
  implicit none

  ! input and output
  type(retr_input_settings), intent(in)  :: RIN
  integer,                   intent(in)  :: par_type, i3
  character(LEN=50),         intent(out) :: i3label

  select case(par_type)
    case(par_type_RERI_beg:par_type_RERI_end-1,             &
         par_type_IMRI_beg:par_type_IMRI_end-1,             &
         par_type_SURF1_land_beg:par_type_SURF1_land_end-1, &
         par_type_SURF2_land_beg:par_type_SURF2_land_end-1, &
         par_type_SURF_water_beg:par_type_SURF_water_end-1)
      ! Use the wavelength in nm as the label.
      write(i3label, '(i15)') nint(1000*RIN%wave(i3))
    case default
      ! By default, just use the index as the label.
      write(i3label, '(i15)') i3
  end select
end subroutine

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine ipar_to_label(RIN, ipar, strpar)
  use mod_retr_settings_derived_type

  implicit none

  ! input and output
  type(retr_input_settings), intent(in)  :: RIN
  integer,                   intent(in)  :: ipar
  CHARACTER(LEN=50),         intent(out) :: strpar
  ! local variables
  integer                                :: i1, i2, i3
  integer                                :: par_type
  character(LEN=50)                      :: i1label, i2label, i3label

  call get_i1_i2_i3_from_ipar(RIN, ipar, i1, i2, i3)

  ! The first component of the label is the parameter type
  par_type = RIN%NDIM%par_type(i1)
  call par_type_to_string(RIN, par_type, i1label)

  ! The second comoent of the label is a number used for grouping.
  write(i2label, '(i15)') i2

  ! The third component of the label refers to some attribute,
  ! usually (but not always) the wavelength.
  ! In the case of the aerosol size distributions this is the radius.
  call i3_to_label(RIN, par_type, i3, i3label)

  write(strpar, '(a)') &
    trim(adjustl(i1label)) // '_' // &
    trim(adjustl(i2label)) // '_' // &
    trim(adjustl(i3label))

end subroutine

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! development needed for sunphotometer+lidar retrievals

      subroutine print_classic_plot(RIN,segment_meas,GOUT)

      use iso_fortran_env, only : output_unit, error_unit
      use mod_par_inv,      only : KIDIM3, KIDIM2, KIMAGE, KSHAPE
      USE mod_par_OS,       only : KSD
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type

      implicit none  
!	------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------
! IN :
      type(retr_input_settings),   intent(in) :: RIN
      type(output_segment_general),intent(in) :: GOUT
      type(segment_data),          intent(in) :: segment_meas
!	------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------
! LOCAL :	  	  
      integer                             :: IDIM1, IDIM2, IDIM3, KN,  &
                                             IW, ipix, npixels, NDIM3, &
                                             I, i1, i2
      character(LEN=20)                   :: KNC, INOISEC
      character(LEN=150)                  :: CFMT 
      character(LEN=255)                  :: file_name
      character(LEN=12),dimension(KIMAGE) :: pdate
      character(LEN=12),dimension(KIMAGE) :: ptime
      real,dimension(KIMAGE)              :: pjday
!	    real                                :: XX,ZZ
!      character(LEN=2)                    :: HH,MM,SS
!      character(LEN=3)                    :: DD
      real,dimension(KSD,KIMAGE)          :: h01
      real,dimension(KSHAPE,KSD,KIMAGE)   :: sph
      integer                             :: id_wr_TL_output	
      real,dimension(NRR,KIMAGE)          :: SDout 

      integer                             :: err_code
      character(LEN=255)                  :: err_msg
      integer                             :: par_type, isurf1, isurf2
      character(LEN=50)                   :: strpar

!	------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------
      isurf1 = 0
      isurf2 = 0
      
      !SDout(:,:) = 0.0
      KN=GOUT%retrieval%par%ngrid
      npixels=segment_meas%npixels

      file_name = RIN%plotting_output_file
      open( newunit=id_wr_TL_output,file=trim(file_name),status='replace',recl=20000)

      do IDIM1=1,RIN%NDIM%n1
      par_type = RIN%NDIM%par_type(IDIM1)
      if(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then      
        do ipix=1,npixels
          do IDIM2=1,RIN%NDIM%n2(IDIM1) 
              NDIM3=RIN%NDIM%n3(IDIM2,IDIM1)
              i1=RIN%NDIM%ISTARSING(IDIM2,IDIM1)
              i2=i1+NDIM3-1
              sph(1:NDIM3,IDIM2,ipix) =   & 
			                    GOUT%retrieval%par%pixel(ipix)%par(i1:i2)
          enddo ! IDIM2
        enddo ! ipix
      elseif(par_type .gt. par_type_AVP_beg  .and. par_type .lt. par_type_AVP_end) then
        do ipix=1,npixels
          do IDIM2=1,RIN%NDIM%n2(IDIM1) 
              NDIM3=RIN%NDIM%n3(IDIM2,IDIM1)
              i1=RIN%NDIM%ISTARSING(IDIM2,IDIM1)
              i2=i1+NDIM3-1
              h01(1:NDIM3,ipix) = GOUT%retrieval%par%pixel(ipix)%par(i1:i2)	
          enddo ! IDIM2
        enddo ! ipix
      elseif(par_type .gt. par_type_SURF1_land_beg  .and. par_type .lt. par_type_SURF1_land_end) then
              isurf1=RIN%NDIM%ISTARSING(1,IDIM1)
      elseif(par_type .gt. par_type_SURF2_land_beg  .and. par_type .lt. par_type_SURF2_land_end) then
              isurf2=RIN%NDIM%ISTARSING(1,IDIM1)
      endif ! par_type .gt. par_type_SHD_beg  .and.
      enddo ! IDIM1
      
      write( id_wr_TL_output,'(500i5)') RIN%DLSF%IWL,RIN%NDIM%n1,  &
               (RIN%NDIM%n2(IDIM1),(RIN%NDIM%n3(IDIM2,IDIM1),      & 
				       IDIM2=1,RIN%NDIM%n2(IDIM1)),IDIM1=1,RIN%NDIM%n1)

! ***  Print Retrieval Results for Visualization     **************	  
      write(KNC,*) RIN%nw
      CFMT = '(8i5,'//trim(adjustl(KNC))//'f15.6,a)'
      write( id_wr_TL_output,trim(CFMT))      &
      npixels,RIN%KNSING,RIN%KNSINGF,RIN%NSD,KN,isurf1,isurf2,RIN%nw,RIN%wave(1:RIN%nw), &
      ' - npixels,KNSING,KNSINGF,NSD,NBIN,isurf1,isurf2,NW,wl(1:NW)'

      write(KNC,*) KN
      write( id_wr_TL_output, '(3A)', advance='no')  &
      ' pixel yyyy-mm-dd  hh:mm:ss     Day_number         ', &
      '  lon            lat       Residual      '

      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0,A)', advance='no') ' AOT_', nint(1000*RIN%wave(iw)), '       '
      enddo
      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0,A)', advance='no') ' SSA_', nint(1000*RIN%wave(iw)), '       '
      enddo
      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0,A)', advance='no') ' SALB_', nint(1000*RIN%wave(iw)), '      '
      enddo
      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0,A)', advance='no') ' RRE_', nint(1000*RIN%wave(iw)), '       '
      enddo
      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0,A)', advance='no') ' RIM_', nint(1000*RIN%wave(iw)), '       '
      enddo

      write( id_wr_TL_output, '(A)', advance='no') ' %sphere        height      '

      do I=1, RIN%KNSING
        call ipar_to_label(RIN, I, strpar)
        write( id_wr_TL_output, '(A35)', advance='no') adjustr(strpar(1:35))
      enddo

!      (ROUT%radius(I),I=1,KN),(I,I=1,RIN%KNSING)

      do iw = 1,RIN%nw
        write(id_wr_TL_output, '(A,i0)', advance='no') '        sza_', nint(1000*RIN%wave(iw))
      enddo
      write(id_wr_TL_output, '(A)', advance='no') '   MASL'
      write(id_wr_TL_output, '(A)') '      x         y'

      call yyyymmdd_hhmmss ( segment_meas,pdate,ptime )
      call day_number_curr_year ( segment_meas,pjday )
  LOOP_PIXEL: DO ipix=1,npixels
  
!      WRITE( id_wr_TL_output,'(i4,2x,a10,2x,a8,f15.7,2i10,1000e15.6)', advance='no') &
      WRITE( id_wr_TL_output,'(i4,2x,a10,2x,a8,f15.7,1000e15.6)', advance='no') &
         ipix,                                             & !pixel
         trim(pdate(ipix)),                                & !yyyy-mm-dd
         trim(ptime(ipix)),                                & !hh:mm:ss
         pjday(ipix),                                      & !Day_number
!tl         segment_meas%pixels(ipix)%irow,                   & !x
!tl         segment_meas%pixels(ipix)%icol,                   & !y
         segment_meas%pixels(ipix)%x,                      & !lon
         segment_meas%pixels(ipix)%y,                      & !lat
         GOUT%retrieval%res%pixel(ipix)%resr(1),           & !Residual
         (GOUT%aerosol%opt%pixel(ipix)%wl(iw)%extt, iw=1,RIN%nw), & !AOT
         (GOUT%aerosol%opt%pixel(ipix)%wl(iw)%ssat, iw=1,RIN%nw), & !SSA
         (GOUT%surface%pixel(ipix)%wl(iw)%salbedo,  iw=1,RIN%nw), & !SALB
         ((GOUT%aerosol%rind%pixel(ipix)%wl(iw)%mreal(I),IW=1,RIN%nw),I=1,RIN%NSD),   & !RRE
         ((GOUT%aerosol%rind%pixel(ipix)%wl(iw)%mimag(I),IW=1,RIN%nw),I=1,RIN%NSD),   & !RIM
         sph(1,1,ipix)*100., & !%sphere
         h01(1,ipix)           !height

!        (SDout(I,ipix),I=1,KN),                            &

      WRITE( id_wr_TL_output,'(1000e35.6)', advance='no') &
         (GOUT%retrieval%par%pixel(ipix)%par(I),I=1,RIN%KNSING)  !parameters i

      WRITE( id_wr_TL_output,'(1000e15.6)', advance='no') &
         (segment_meas%pixels(ipix)%meas(iw)%sza,IW=1,RIN%nw),        & !sza
         segment_meas%pixels(ipix)%MASL                                 !masl

      WRITE( id_wr_TL_output,'(2i10)' ) &
         segment_meas%pixels(ipix)%irow, & !x
         segment_meas%pixels(ipix)%icol    !y

   ENDDO LOOP_PIXEL  ! ipix
! ***  Write Size Distribution     *********************************	  
goto 33
      WRITE( id_wr_TL_output,*) & 
       'Size Distribution :'
      WRITE( id_wr_TL_output,'(a11,4x,1000f15.7)')      & 
       'Day_number :',                                  &
      (pjday(ipix),ipix=1,npixels) 
      WRITE( id_wr_TL_output,'(a5,10x,20000(5x,A10))')  &
       'Date:',                                         &
       (pdate(ipix),ipix=1,npixels)
      WRITE( id_wr_TL_output,'(a5,10x,20000(7x,A8))')   & 
       'Time:',                                         &
      (ptime(ipix),ipix=1,npixels) 

      DO I=1,KN
      WRITE(id_wr_TL_output,'(20000e15.6)')             &
		  GOUT%retrieval%par%radius(I),                     &
		 (SDout(i,ipix),ipix=1,npixels)
      ENDDO ! I  
33 continue
! ***  Print Residual              *********************************
 
      WRITE(id_wr_TL_output,*)
	   
      IF(RIN%IPFP%INVSING .ge. 2) THEN
  LOOP_PIXEL1: DO ipix=1,npixels
      WRITE(INOISEC,*) RIN%NOISE%INOISE
      !CFMT = '(f10.5,'//trim(adjustl(INOISEC))//  &
	    !                              '(i5,a,e14.5,f15.5,a),i5,4x,a,i5,a,2x,i5,a)'
      !WRITE(id_wr_TL_output,TRIM(CFMT)) GOUT%retrieval%res%pixel(ipix)%res, &
      !(I,':',GOUT%retrieval%res%pixel(ipix)%resa(I), &
      !       GOUT%retrieval%res%pixel(ipix)%resr(I), &
      !' %',I=1,RIN%NOISE%INOISE),ipix,'  - Residual after LP=',  &
      !GOUT%retrieval%res%niter,' - iteration', ipix,' - PIXEL'
      CFMT = '(f10.5,'//trim(adjustl(INOISEC))//  &
	                                  '(i5,a,e14.5,f15.5),i5,4x,a,i5,a,2x,i5,a)'
      WRITE(id_wr_TL_output,TRIM(CFMT)) GOUT%retrieval%res%pixel(ipix)%res, &
      (I,':',GOUT%retrieval%res%pixel(ipix)%resa(I), &
             GOUT%retrieval%res%pixel(ipix)%resr(I), &
      I=1,RIN%NOISE%INOISE),ipix,'  - Residual after LP=',  &
      GOUT%retrieval%res%niter,' - iteration', ipix,' - PIXEL'

  ENDDO LOOP_PIXEL1  ! IMAGE
      WRITE(id_wr_TL_output,*)
!      CFMT = '(f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5,a),9x,a,i5,a)'
!      WRITE(id_wr_TL_output,TRIM(CFMT))  GOUT%retrieval%res%rest,     &
!      (I,':',GOUT%retrieval%res%resat(I),GOUT%retrieval%res%resrt(I), &
!      ' %',I=1,RIN%NOISE%INOISE),  &
!      '  - Residual after LP=',GOUT%retrieval%res%niter,' - iteration TOTAL IMAGE'
      CFMT = '(f10.5,'//trim(adjustl(INOISEC))//'(i5,a,e14.5,f15.5),9x,a,i5,a)'
      WRITE(id_wr_TL_output,TRIM(CFMT))  GOUT%retrieval%res%rest,     &
      (I,':',GOUT%retrieval%res%resat(I),GOUT%retrieval%res%resrt(I), &
      I=1,RIN%NOISE%INOISE),  &
      '  - Residual after LP=',GOUT%retrieval%res%niter,' - iteration TOTAL IMAGE'

      ELSEIF(RIN%IPFP%INVSING .ge. 0 .and. RIN%IPFP%INVSING .le. 1) THEN
  LOOP_PIXEL2: DO ipix=1,npixels
      WRITE(INOISEC,*) RIN%NOISE%INOISE
!      CFMT = '(f10.5,'//trim(adjustl(INOISEC))//  &
!	                                  '(i5,a,e14.5,f15.5,a),i5,a)'
!      WRITE(id_wr_TL_output,TRIM(CFMT)) GOUT%retrieval%res%pixel(ipix)%res,              &
!      (I,':',GOUT%retrieval%res%pixel(ipix)%resa(I),GOUT%retrieval%res%pixel(ipix)%resr(I),' %',  &
!      I=1,RIN%NOISE%INOISE),ipix,' - PIXEL'
      CFMT = '(f10.5,'//trim(adjustl(INOISEC))//  &
	                                  '(i5,a,e14.5,f15.5),i5,a)'
      WRITE(id_wr_TL_output,TRIM(CFMT)) GOUT%retrieval%res%pixel(ipix)%res,              &
      (I,':',GOUT%retrieval%res%pixel(ipix)%resa(I),GOUT%retrieval%res%pixel(ipix)%resr(I),  &
      I=1,RIN%NOISE%INOISE),ipix,' - PIXEL'
  ENDDO LOOP_PIXEL2  ! IMAGE
      ENDIF ! RIN%IPFP%INVSING .ge. 2
      ! Size Distribution   
      call print_size_distribution(id_wr_TL_output,RIN,segment_meas,GOUT)

      
      close( id_wr_TL_output, iostat=err_code, iomsg=err_msg)
	  
      if (err_code /= 0) then
         write(error_unit, '("iP_write_output.f90: attempt to close ",A," failed: ",A)') &
         trim(RIN%plotting_output_file), trim(err_msg)
      else
         write(output_unit, '("iP_write_output.f90: ",A," done ")') &
         trim(RIN%plotting_output_file)
      end if
	  
      flush(output_unit)
      flush(error_unit)
      return
      end subroutine print_classic_plot

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine print_phmx_file ( RIN, segment_meas, GOUT )

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_time_utils
      use mod_sdata_derived_type
      	  
      implicit none
! -----------------------------------------------------------------------------------
      type(retr_input_settings),    intent(in) :: RIN
      type(output_segment_general), intent(in) :: GOUT
      type(segment_data),           intent(in) :: segment_meas
! -----------------------------------------------------------------------------------
      integer :: II, ipix, IW, ISD, iu_phmx, IDIM1, NDIM3, npixels, par_type
      character(LEN=255) :: file_name
      character(LEN=12), dimension(KIMAGE) :: pdate, ptime
      real,dimension(KIDIM3,KIDIM2,KIMAGE) :: sph
! -----------------------------------------------------------------------------------
      npixels = segment_meas%npixels

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)        
        if(par_type .gt. par_type_SHD_beg  .and. par_type .lt. par_type_SHD_end) then
          call unpack_retr_param_to_print ( RIN,GOUT%retrieval%par,IDIM1,npixels,sph )
          if(RIN%NDIM%n2(IDIM1) .eq. 1 .and. RIN%NDIM%n3(1,IDIM1) .eq. 2) & 
          sph(2,:,:) = sph(1,:,:)     	  
        exit
        endif ! par_type .gt. par_type_SHD_beg  .and. 
      enddo ! IDIM1


      file_name = trim(RIN%DLSF%external_file_path)//"phase_mtrx"//".txt"
      open( newunit=iu_phmx, file=trim(file_name), status='replace')

      write(iu_phmx,*)
      write(iu_phmx,*)	"								Phase Matrix"
      write(iu_phmx,*)

      call yyyymmdd_hhmmss ( segment_meas,pdate,ptime )

      select case (RIN%DLSF%keyEL)
      case(1)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(2)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca        ph12/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(3)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca        ph12/sca        ph22/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD)

              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(4)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(5)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca        ph34/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph34(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      case(6)
        do ipix=1,segment_meas%npixels
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          write(iu_phmx,'(a,i0,a,a10,a,a8)') 'ipix=',ipix,'  yymmdd = ',pdate(ipix),'  hhmmss = ',ptime(ipix)
          write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
          do IW=1,RIN%nw
            do ISD=1,RIN%NSD
              NDIM3 = RIN%NDIM%n3(ISD,IDIM1)
              write(iu_phmx,10) &
              'wl=',RIN%wave(IW),'  ISD=',ISD, &
              '  sca=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD)* &
                       GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ext(ISD), &
              '  ssa=',GOUT%aerosol%opt%pixel(ipix)%wl(IW)%ssa(ISD), &
              '  n=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mreal(ISD), &
              '  k=',GOUT%aerosol%rind%pixel(ipix)%wl(IW)%mimag(ISD), &
              '  fsphere=',sph(1:NDIM3,ISD,ipix)
              write(iu_phmx,'(2a)')	"         # 	        angle        ", &
              "ph11/sca        ph12/sca        ph22/sca        ph33/sca        ph34/sca        ph44/sca"
              do II = 1,GOUT%aerosol%phmx%nangle
                write(iu_phmx,11) II,GOUT%aerosol%phmx%angle(II),                    &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph11(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph12(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph22(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph33(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph34(II,ISD), &
                                  GOUT%aerosol%phmx%pixel(ipix)%wl(IW)%ph44(II,ISD)
              enddo ! II
            enddo ! ISD
          enddo ! IW
        enddo ! ipix
      end select

      write(iu_phmx,'(a)') 	"--------------------------------------------------------------------------"
10    format(a,f9.6,a,i0,4(a,e12.5),a,25e12.5)
11    format(6x,I4,3x,f12.3,2x,6(e14.5,2x))

      close(iu_phmx)

      return
      end subroutine print_phmx_file

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine yyyymmdd_hhmmss ( segment_meas,pdate,ptime )

      use mod_par_inv, only : KIMAGE
      use mod_sdata_derived_type
      
      implicit none
! ----------------------------------------------------------------	   
      type(segment_data),                 intent(in)   ::  segment_meas
      character(len=12),dimension(KIMAGE),intent(out)  ::  pdate,ptime
! ----------------------------------------------------------------	  
      integer  ::  ipix
! ----------------------------------------------------------------	  
! http://www.cplusplus.com/reference/ctime/strftime/

      pdate(:) = ' '
      ptime(:) = ' '
      do ipix=1,segment_meas%npixels
         call convert_time_to_string(segment_meas%pixels(ipix)%t, "%F", pdate(ipix))
         call convert_time_to_string(segment_meas%pixels(ipix)%t, "%X", ptime(ipix))
      enddo ! ipix

      return
      end subroutine yyyymmdd_hhmmss
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine day_number_curr_year ( segment_meas,pjday )

      use mod_par_inv, only : KIMAGE
      use mod_sdata_derived_type
      
      implicit none
! ----------------------------------------------------------------	   
      type(segment_data),      intent(in)   ::  segment_meas
      real,dimension(KIMAGE),  intent(out)  ::  pjday
! ----------------------------------------------------------------	  
      integer           ::  ipix
	    real              ::  XX,ZZ
      character(LEN=2)  ::  HH,MM,SS
      character(LEN=3)  ::  DD
! ----------------------------------------------------------------	  
! http://www.cplusplus.com/reference/ctime/strftime/

! calculate day number in current year (like proveded in Aeronet retrieval product)
        pjday(:) = 0.0
        do ipix=1,segment_meas%npixels
          call convert_time_to_string(segment_meas%pixels(ipix)%t, "%j", DD)
          read(DD,*) XX
          pjday(ipix) = XX
          HH = ' '
          MM = ' '
          SS = ' '      
          call convert_time_to_string(segment_meas%pixels(ipix)%t, "%H", HH)
          call convert_time_to_string(segment_meas%pixels(ipix)%t, "%M", MM)
          call convert_time_to_string(segment_meas%pixels(ipix)%t, "%S", SS)
          !write(*,*) 'HH=',HH,'  MM=',MM,'  SS=',SS
          ZZ = 0.
          read(HH,*) XX
          ZZ=ZZ+XX*3600.
          read(MM,*) XX
          ZZ=ZZ+XX*60.		 		 
          read(SS,*) XX
          ZZ=ZZ+XX
          ZZ=ZZ/(24.*3600.)		
          pjday(ipix)=pjday(ipix)+ZZ
           !write(*,*) 'ipix=','  pjday(ipix)=',pjday(ipix)
        enddo ! ipix

      return
      end subroutine day_number_curr_year

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

