! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../../../constants_set/mod_globals.inc"
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine spheroid_package (                          &
                                    NSD,NBIN,RADIUS,SD,      & ! IN
                                    SHD,NSHAPE1,             &
                                    ind_wl,WAVE,RREAL,RIMAG, &
                                    use_models,              &
                                    NANG,ANGL,               & ! OUT
                                    GOUT_particles_opt_pixel_wl,   &
                                    GOUT_particles_phmx_pixel_wl,  &
                                    DLSF,KERNELS1,KERNELS2   & 
                                   )
      use mod_par_inv, only : KSHAPE,KIDIM3
      use mod_par_OS,  only : KSD
      use mod_type_DLS
      use mod_alloc_kernels
      !use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------
! IN:
      real,dimension(KIDIM3,KSD), intent(in)  ::  RADIUS, SD
      real,dimension(KSHAPE,KSD), intent(in)  ::  SHD  
      integer,dimension(KSD),     intent(in)  ::  NBIN, NSHAPE1
      real,dimension(KSD),        intent(in)  ::  RREAL, RIMAG !,SSA,EXT
      real,                       intent(in)  ::  WAVE
      integer,                    intent(in)  ::  NSD, ind_wl
      type(iP_flags_for_DLS),     intent(in)  ::  DLSF
      logical,                    intent(in)  ::  use_models
! -----------------------------------------------------------------------
! OUT:
      integer,                        intent(out)    ::  NANG
      real,dimension(KMpar),          intent(out)    ::  ANGL
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_particles_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_particles_phmx_pixel_wl
!	----------------------------------------------------------------------
! INOUT:
      type(kernels_triangle_bin), intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),intent(inout)  ::  KERNELS2
! -----------------------------------------------------------------------
      integer ::  ISD
      real ::  temp
      logical, save :: lfirst = .true.
!	------------------------------------------------------------------------------------------------------
	       
      GOUT_particles_phmx_pixel_wl%ph11(:,:) = 0.
      GOUT_particles_phmx_pixel_wl%ph12(:,:) = 0.
      GOUT_particles_phmx_pixel_wl%ph22(:,:) = 0.
      GOUT_particles_phmx_pixel_wl%ph33(:,:) = 0.
      GOUT_particles_phmx_pixel_wl%ph34(:,:) = 0.
      GOUT_particles_phmx_pixel_wl%ph44(:,:) = 0.
      GOUT_particles_opt_pixel_wl%ext(:) = 0.
      GOUT_particles_opt_pixel_wl%ssa(:) = 0.

#ifdef GRASP_MODELS
      if(use_models) then
        do ISD=1,NSD
            call phmx_models_wl( DLSF, ind_wl, WAVE,                     &
                              NBIN(ISD), SD(1:NBIN(ISD),ISD),            &
                              NANG, ANGL,                                &
                              GOUT_particles_opt_pixel_wl%ssa(ISD),      &
                              GOUT_particles_opt_pixel_wl%ext(ISD),      &
                              GOUT_particles_phmx_pixel_wl%ph11(:,ISD),  &
                              GOUT_particles_phmx_pixel_wl%ph12(:,ISD),  &
                              GOUT_particles_phmx_pixel_wl%ph22(:,ISD),  &
                              GOUT_particles_phmx_pixel_wl%ph33(:,ISD),  &
                              GOUT_particles_phmx_pixel_wl%ph34(:,ISD),  &
                              GOUT_particles_phmx_pixel_wl%ph44(:,ISD),  &
                              KERNELS2)
            if(stop_report%status .eqv. .true.) return
        enddo ! ISD
      endif
#endif
!#else
!#error Code can not be compiled because models directory is not present!
!#endif

      if(.not. use_models) then
        do ISD=1,NSD
!         WRITE(*,*) ' before PHASE_KERNEL keyEL=',DLSF%keyEL
!         WRITE(*,*) RREAL(ISD), RIMAG(ISD),ISD,WAVE,ind_wl,  &
!                  ' RREAL  RIMAG  ISD  WAVE  ind_wl'

! Calculate Optical Characteristics      
          call PHASE_KERNEL_W (                                     & ! IN
                                NBIN(ISD),RADIUS(:,ISD),SD(:,ISD),  &
                                SHD(:,ISD),NSHAPE1(ISD),            &
                                ind_wl,WAVE,RREAL(ISD),RIMAG(ISD),  &
                                NANG,ANGL,DLSF,                     & ! OUT
                                GOUT_particles_opt_pixel_wl%ssa(ISD), &
                                GOUT_particles_opt_pixel_wl%ext(ISD), &
                                !GOUT_particles_opt_pixel_wl%aext(ISD), &
                                GOUT_particles_phmx_pixel_wl%ph11(:,ISD),  &
                                GOUT_particles_phmx_pixel_wl%ph12(:,ISD),  &
                                GOUT_particles_phmx_pixel_wl%ph22(:,ISD),  &
                                GOUT_particles_phmx_pixel_wl%ph33(:,ISD),  &
                                GOUT_particles_phmx_pixel_wl%ph34(:,ISD),  &
                                GOUT_particles_phmx_pixel_wl%ph44(:,ISD),  &
                                KERNELS1,KERNELS2                          &
                              )
          if(stop_report%status .eqv. .true.) return
!      IF(DLSF%keyEL .NE. 4) WRITE(*,*) "NEL.NE.4 !!!"

!      temp = GOUT_main_pixel_wl%ssa(ISD)*GOUT_main_pixel_wl%ext(ISD)
!      if(DLSF%keyEL .gt. 0) GOUT_phmx_pixel_wl%ph11(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph11(1:NANG,ISD)/temp
!      if(DLSF%keyEL .gt. 1) GOUT_phmx_pixel_wl%ph12(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph12(1:NANG,ISD)/temp
!      if(DLSF%keyEL .gt. 2) GOUT_phmx_pixel_wl%ph22(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph22(1:NANG,ISD)/temp
!      if(DLSF%keyEL .gt. 3) GOUT_phmx_pixel_wl%ph33(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph33(1:NANG,ISD)/temp
!      if(DLSF%keyEL .gt. 4) GOUT_phmx_pixel_wl%ph34(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph34(1:NANG,ISD)/temp
!      if(DLSF%keyEL .gt. 5) GOUT_phmx_pixel_wl%ph44(1:NANG,ISD) = GOUT_phmx_pixel_wl%ph44(1:NANG,ISD)/temp       
		   	 
!write(*,*) ISD,GOUT_main_pixel_wl%ssa(ISD),GOUT_main_pixel_wl%ext(ISD),GOUT_phmx_pixel_wl%ph11(NANG,ISD), &
!           ' - ISD,ssa,ext,ph11 - after PHASE_KERNEL_W'
        enddo ! DO ISD=1,NSD
      endif

      return
      end subroutine spheroid_package

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      
      SUBROUTINE PHASE_KERNEL_W (                          & ! IN
                                 NBIN,RADIUS,SD1,          &
                                 RD1,NSHAPE,               &
                                 IW,RWL,RNW,RKW,           &
                                 IANGL,ANG,DLSF,           & ! OUT
                                 SSA,EXT1,                 &
                                 PTP11,PTP12,PTP22,        &
                                 PTP33,PTP34,PTP44,        &
                                 KERNELS1,KERNELS2         &
                                )
	 
!C*****************************************************
!C***This subroutine calculates the 
!C    EXTINCTION, SSA, PHASE FUNCTION using kernel matrices                     
!C
!C   INPUT:
!C*   ISD    I  - number of size distribution component
!C*   IEL    I  - number of needed elements of phase matrix
!C*            =1 - P11 only
!C*            =6 - P11, P12, P22, P33, P34, P44 
!C*   NBIN   I  - the number of the bins in the
!C*                       size distributions
!C*   RADIUS R(KBIN) - the radii corresponding to  
!C*                the bins in the size distributions
!C*               size distributions
!C*   SD     D(NBIN) - size distribution
!C*                    Usually dV/dlnR (it can be also for 
!C*                    radii, numbers, areas,
!C*                    but should agree with KERNEL)
!C*   WAV    R       - wavelength
!C*   RREAL  R       - real part of the refractive index
!C*   RIMAG  R       - imaginary part of the refractive index
!C***************************************************** 
!C*  IANGL        - the number of angles used for 
!C*                   phase function
!C*   ANG    R(IANG) - the values of the angles for the phase 
!C*                        function modeling;
!C*   SSA    R       - single scattering albedo
!C*   EXT    R       - extinction
!C*   PTPii     R(IANG)- phase matrix (UNNORMALIZED !!!)
!C*****************************************************
!C
!C
!C***************************************************************** c
!C***   Definition of the parameters for OPTCHAR1               *** c
!C***************************************************************** c
!c ** 11/29/05 Size Distribution option: SD table or
!c ** LogNormal parameters (only KN=-KN (logarithmic intervals) &
!c ** IA=1 - trapezoidal approximation case)
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 10/09/03 don't calculate f34 and f44 if key_f344=1
!c ** 08/27/03 smooth f33 and f44 for 40< angle <60 
!c ** 05/05/03 this version can be used to retrieve an axis ratio or
!c **          axis ratio distribution
!c **************************************************************** c
!c **   02/28/03                                                 ** c
!c **   Program calculates optical characteristics for given     ** c
!c **   size distribution, refractive index, axis ratio          ** c
!c **   distribution and wavelength                              ** c
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c 
!c **   key  = 1 - create fixed kernels (for fixed axis          ** c 
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **          2 - read fixed kernels from 'Rke...fix...'  files ** c
!c **          3 - create fixed kernels but don't save them      ** c
!c **          4 - don't create fixed kernels, calculate         ** c
!c **              opt.characteristics from original kernels     ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   key_f11= 0 - calculate scattering matrix                 ** c
!c **            1 - calculate only phase function               ** c
!c **   key_f344=0 - calculate scattering matrix                 ** c
!c **            1 - don't calculate f34 and f44                 ** c
!c **   key_org=0 - read original kernels simultaneously make    ** c 
!c **               angle interpolation and calculate opt.char.  ** c
!c **           1 -  -"-, save new kernels in                    ** c
!c **                /distname_N/ directory, STOP                ** c
!c **   key_fx works when key=1                                  ** c
!c **        =0 -  create fixed kernels (for fixed axis          ** c 
!c **              ratio distr.) and save them                   ** c
!c **              into 'Rke...fix...' files and calculate       ** c
!c **              opt.characteristics                           ** c
!c **         1 - save fixed kernels with original kernel format ** c
!c **             in order to be used as input kernels;          ** c
!c **             'Rke...fix' kernels have to be renamed and moved* c
!c **             into directory 'dir_name'(see 'matrix_fixget.f')* c
!c **             The files can not be used if key=2.            ** c
!c **                                                            ** c
!c **   WL   - wavelength                                        ** c
!c **   RN   - real part of the refractive index                 ** c
!c **   RK   - imaginary part of the refractive index            ** c
!c **   rgmin,rgmax &                                            ** c 
!c **   wlmin,wlmax - min,max radii and wlmin,wlmax wavelengths  ** c
!c **                 that are used to recalculate grid radii for** c
!c **                 fixed kernels. New input file              ** c
!c **                'grid1.dat.new' will be created if key=1    ** c
!c **                 or key_org=1. Use key_grid1 to choose      ** c
!c **                 'grid1.dat' or 'grid1.dat.new' will be read** c
!c **                 for further calculations                   ** c  
!c **   key_SD=0 - read Size Distribution table dV/dlnR          ** c
!c **         =1 - calculate Size Distribution for grid radii    ** c
!c **              using Log Normal function                     ** c
!c **   ID    - dimension of d(...)/dlnR or d(...)/dR            ** c
!c **       = 0 - number                                         ** c
!c **       = 1 - radius                                         ** c
!c **       = 2 - area                                           ** c
!c **       = 3 - volume                                         ** c
!c **   NMD   - number of modes (up to 2)                        ** c
!c **   KN   - number of grid radii                              ** c
!c **   grid(KN) - grid radii                                    ** c
!c **   SD(KN)   - size distribution for grid radii              ** c
!c **   (CM(i),SM(i),RMM(i),i=1,NMD) - size distribution         ** c
!c **   function (LogNormal) parameters:                         ** c
!c **                         CM - concentration                 ** c
!c **                         SM - standard deviation            ** c
!c **                         RMM - median radius                ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **   KR  - number of axis ratios                              ** c
!c **   R(KR)  - grid axis ratios                                ** c
!c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
!c **   KM   - number of scattering angles                       ** c
!c **   ANGLE(KM) - scattering angles                            ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   ext     - extinction                                     ** c
!c **   albedo  - albedo                                         ** c
!c **   f... - scattering matrix elements                        ** c
!c **************************************************************** c

      USE mod_par_DLS,     only : KMpar,KNpar	  
      USE mod_par_DLS_bin, only : KWLpar,KCpar
      USE mod_type_DLS
      USE mod_par_inv,     only : KPARS,KSHAPE,KIDIM3
      USE mod_index_cloud	  
      USE mod_retr_settings_derived_type
      USE mod_alloc_kernels
      USE mod_stop_report

      implicit none
!---------------------------------------
! IN:	  
      INTEGER, INTENT(IN)                     :: NBIN,NSHAPE,IW
      !REAL,    INTENT(IN),  DIMENSION(KPARS)  :: RADIUS,SD1
      REAL,    INTENT(IN),  DIMENSION(KIDIM3) :: RADIUS,SD1

      REAL,    INTENT(IN),  DIMENSION(KSHAPE) :: RD1
      REAL,    INTENT(IN)                     :: RWL,RNW,RKW
      TYPE(iP_flags_for_DLS),INTENT(IN)       :: DLSF

! OUT :	  
      INTEGER, INTENT(OUT)                    :: IANGL
      REAL,    INTENT(OUT), DIMENSION(KMpar)  :: ANG
      REAL,    INTENT(INOUT)                    :: SSA,EXT1 !,ABS
      REAL,    INTENT(INOUT), DIMENSION(KMpar)  :: PTP11,PTP12,PTP22,    &
                                                   PTP33,PTP34,PTP44
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2

!---------------------------------------
! LOCAL:
      integer                     :: I,IW1 
      real,  dimension(1,KNpar)   :: RRR,AR
      real,  dimension(1)         :: AC
      real,  dimension(1,1)       :: CM,SM,RMM
      type(DLS_CONFIG), save      :: PHASE_MATRIX
!---------------------------------------
      integer,                    save :: key,keySUB,keyLS,key_RD1,key_RD, &
                                          KN,KM,KR,KC,                     &
                                          NWL,LB,LE
      real,                       save :: xnmin,xnmax,xkmin,xkmax,  &		                                      
                                          pomin,pomax
      real,   dimension(KWLpar),  save :: WL
      real,   dimension(KMpar),   save :: ANG_save=0.0
      real*4, dimension(KSHAPE),  save :: R	  
      logical,                    save :: lfirst = .true.
      real,   dimension(KCpar)         :: RC2  
      real,   dimension(KCpar,KNpar)   :: SD_bin
      real                             :: tiny = 1e-2
      logical :: status
!	------------------------------------------------------------------------------------------------------
!
! ** READ DLS_bin INPUT
!       

!write(*,*) 'in phase_kernel: distname_O=',distname_O 
      if(lfirst) then
        if(DLSF%IWL .eq. 0) then
          call DLS_read_input_bin(DLSF%internal_file_path,DLSF%distname_O,PHASE_MATRIX)
        elseif(DLSF%IWL .eq. 1) then
          call DLS_read_input_bin(DLSF%internal_file_path,DLSF%distname_N,PHASE_MATRIX)   
        else 
          write(tmp_message,'(a,i0,a)') 'DLSF%IWL = ',DLSF%IWL,' - value is not valid'
          G_ERROR(trim(tmp_message))
        endif !
        if(stop_report%status .eqv. .true.) return
        if(PHASE_MATRIX%DLSIN%key .ne. 2) then
          write(tmp_message,'(a,i0,a)') 'key = ',PHASE_MATRIX%DLSIN%key,' .ne. 2'
          G_ERROR(trim(tmp_message))
        endif ! key.ne.2

      !if(DLSF%keyEL .ne. PHASE_MATRIX%DLSIN%keyEL) then
         !write(*,'(2(a,i2))') 'STOP in PHASE_KERNEL_W keyEL=',PHASE_MATRIX%DLSIN%keyEL,  &
         !           ' .ne. DLSF%keyEL=',DLSF%keyEL
         !stop
      !endif ! keyEL

      keySUB    = PHASE_MATRIX%DLSIN%keySUB
      keyLS     = PHASE_MATRIX%DLSIN%keyLS
      key_RD1   = PHASE_MATRIX%DLSIN%key_RD1
      key_RD    = PHASE_MATRIX%DLSIN%key_RD
      KN        = PHASE_MATRIX%DLSIN%KN
      KM        = PHASE_MATRIX%DLSIN%KM
      KC        = PHASE_MATRIX%DLSIN%KC      	  
      xnmin     = PHASE_MATRIX%DLSIN%xnmin
      xnmax     = PHASE_MATRIX%DLSIN%xnmax
      xkmin     = PHASE_MATRIX%DLSIN%xkmin
      xkmax     = PHASE_MATRIX%DLSIN%xkmax
      pomin     = PHASE_MATRIX%DLSIN%pomin
      pomax     = PHASE_MATRIX%DLSIN%pomax
      ANG_save(1:KM) = PHASE_MATRIX%DLSIN%ANGLE(1:KM)
      R(1:NSHAPE)    = PHASE_MATRIX%DLSIN%R(1:NSHAPE)
		
      if(DLSF%IWL.EQ.0) then  ! bi-commponent aerosol with lognormal SDs
         key = 0  	
         NWL=1
      else
         key = 2
         NWL = PHASE_MATRIX%DLSIN%NWL
         if(NWL .gt. KWLpar) then
           write(tmp_message,'(2(a,i0),a)') &
           'Number of wavelengths for spectral kernels NWL = ',NWL, &
           'can not be larger than KWLpar = ',KWLpar,' in mod_par_DLS_bin'
           G_ERROR(trim(tmp_message))
         endif
         WL(1:NWL) = PHASE_MATRIX%DLSIN%WL_bin(1:NWL)
      endif
		  lfirst=.false.
      endif ! lfirst
	         
! ** OPTICAL CARACTERISTIC Subroutine Input preparation for each OPTCHAR_bin CALL
      !WRITE(*,*) 'SD1: NBIN=',NBIN
      !WRITE(*,'(10e16.5)') SD1(1:NBIN) 

      IF(DLSF%IWL.EQ.0) THEN  
         IW1=1
         WL(IW1)=RWL
         RC2(1:KC)=1.e-30   ! is not used with key=0 in OPTCHAR_bin
         if(DLSF%key.eq.0) then                           
! Triangle bins
            SD_bin(1,1:NBIN) = SD1(1:NBIN) 
         else if(DLSF%key.eq.3) then                
! Bi-modal lognormal distribution
            CM (1,1) = SD1(3)
            SM (1,1) = SD1(2)
            RMM(1,1) = SD1(1)
!             write(*,*) CM(1,1),SM(1,1),RMM(1,1),' CM(1,1),SM(1,1),RMM(1,1) in PHASE_KERNEL_W'
!             write(*,*) KN,NBIN,RADIUS(1),RADIUS(NBIN),' KN,NBIN,RADIUS(1),RADIUS(NBIN) in PHASE_KERNEL_W'
            call SIZEDISDN (                           &
                              0,                       & ! IPRI
                             -KN,1,3,1,1,              & ! KN,ID,NSD,NMD
                              CM (1,1),                &
                              SM (1,1),                &
                              RMM(1,1),                &
                              RADIUS(1),RADIUS(NBIN),  &
                              RRR(1,:),AR(1,:),AC(1)   &
                           )
            SD_bin(1,1:KN) = AR(1,1:KN)
         endif ! DLSF%key		
      ELSE ! DLSF%IWL .NE. 0
! Lognormal bins
         KC = NBIN
         IW1=IW
        !do i = 1,NWL
        !write(*,*) "WL(",i,") = ",WL(i)
        !enddo
        !write(*,*) "RWL = ",RWL
        !write(*,*)"In SUBROUTINE forw_IMAGE_I "
         call R_subset_index (tiny,            & ! IN
                              NWL, WL(1:NWL),  &
                              1, RWL,          &
                              IW1, status      & ! OUT

                             )
          if(.not. status) then
            write(tmp_message,'(a,f7.4,a,a,20f7.4)') &
            'Can not find index for wavelength, wl = ',RWL, &
            NEW_LINE('A'),'WL: ',WL(1:NWL)
            G_ERROR(trim(tmp_message))
          endif

         RC2(1:KC)=SD1(1:KC)
         DO I=1,KN
         SD_bin(1:KC,I)=1.e-30 ! is not used with key=2 in OPTCHAR_bin
         ENDDO ! I
      ENDIF ! DLSF%IWL


      if(DLSF%key.eq.3) then
         PHASE_MATRIX%DLSIN%grid(1:KN) = RRR(1,1:KN)
      else  
         KN = NBIN
         PHASE_MATRIX%DLSIN%grid(1:KN) = RADIUS(1:KN)
      endif ! DLSF%key.eq.3
	  
      KR=NSHAPE
      PHASE_MATRIX%DLSIN%RD(1:KR)= RD1(1:KR)
      PHASE_MATRIX%DLSIN%R(1:KR) = R(1:KR)

      PHASE_MATRIX%DLSIN%RN_bin(IW1)=RNW
      PHASE_MATRIX%DLSIN%RK_bin(IW1)=RKW
	  
      LB=IW1
      LE=IW1

      IANGL=KM
      ANG(1:KM)=ANG_save(1:KM)
goto 111
      WRITE(*,*) 'before OPTCHAR_bin:'
      WRITE(*,*) LB,LE,IW1,' - LB,LE,IW1'
      WRITE(*,*) RC2(1:KC),' - RC2'
      WRITE(*,*) PHASE_MATRIX%DLSIN%RN_bin(IW1),PHASE_MATRIX%DLSIN%RK_bin(IW1),' - RN(IW),RK(IW)'
      WRITE(*,*) KC,NWL,WL(1:NWL),' - KC,NWL,WL(1:NWL)'
      WRITE(*,*) key,DLSF%keyEL,keySUB,keyLS,key_RD1,key_RD,' - key, DLSF%keyEL, keySUB, keyLS, key_RD1, key_RD'
      WRITE(*,*) 'ANGLE: KM=',KM
      WRITE(*,'(10f16.5)') ANG(1:KM)
      WRITE(*,*) 'SD_bin: KN=',KN
      WRITE(*,'(10e16.5)') SD_bin(1,1:KN) 
      WRITE(*,*) 'grid: '
      WRITE(*,'(10e16.5)') PHASE_MATRIX%DLSIN%grid(1:KN)
      WRITE(*,*) 'R: '
      WRITE(*,'(10f16.5)') PHASE_MATRIX%DLSIN%R(1:KR)
      WRITE(*,*) 'RD: KR=',KR
      WRITE(*,'(10e16.5)') PHASE_MATRIX%DLSIN%RD(1:KR)
      WRITE(*,*) pomin,pomax,xnmin,xnmax,xkmin,xkmax,'  - pomin, pomax,xnmin, xnmax, xkmin, xkmax'      
      WRITE(*,*) 'XGRID1: '
      WRITE(*,*) 'KN1,KM1,NRATN: ',PHASE_MATRIX%XGRID1%KN1, &
                                   PHASE_MATRIX%XGRID1%KM1, &
                                   PHASE_MATRIX%XGRID1%NRATN
      WRITE(*,*) 'XWL: ',PHASE_MATRIX%XGRID1%XWL
      WRITE(*,*) 'ANGLE1: '
      WRITE(*,'(10e16.5)') PHASE_MATRIX%XGRID1%ANGLE1(1:PHASE_MATRIX%XGRID1%KM1)
      WRITE(*,*) 'grid1: '                                         
      WRITE(*,'(10e16.5)') PHASE_MATRIX%XGRID1%grid1(1:PHASE_MATRIX%XGRID1%KN1)
      WRITE(*,*) 'RATIO: '
      WRITE(*,'(10e16.5)') PHASE_MATRIX%XGRID1%RATIO
      WRITE(*,'(a,/)') 'before OPTCHAR_bin' 
111 continue
!
! ** COMPUTE OPTICAL CARACTERISTICS
!
!  For Subroutine OPTCHAR_bin
!   IN:                 key,key_RD,keyEL,keySUB                   
!                      ,NWL,WL_bin,RN_bin,RK_bin,KN,grid,SD_bin   
!                      ,KR,R,RD,KM,ANGLE                         
!                      ,pomin,pomax,distname_O,distname_N
!                      ,KC,RC2,LB,LE                   
!                      ,xnmin,xnmax,xkmin,xkmax,key_RD1
!
!
!   OUT :              xext_bin,xabs_bin,xsca_bin,albedo_bin     
!                     ,f11_bin,f12_bin,f22_bin                  
!					            ,f33_bin,f34_bin,f44_bin                   
!					           !

      PHASE_MATRIX%KERNEL%internal_file_path = DLSF%internal_file_path
      PHASE_MATRIX%KERNEL%external_file_path = DLSF%external_file_path

      PHASE_MATRIX%KERNEL%distname_O = DLSF%distname_O
      PHASE_MATRIX%KERNEL%distname_N = DLSF%distname_N
      !WRITE(0,*) 'KERNEL: '
      !WRITE(0,*) 'internal_file_path=',trim(PHASE_MATRIX%KERNEL%internal_file_path)      
      !WRITE(0,*) 'external_file_path=',trim(PHASE_MATRIX%KERNEL%external_file_path)
      !WRITE(0,*) 'distname_O=',trim(PHASE_MATRIX%KERNEL%distname_O)
      !WRITE(0,*) 'distname_N=',trim(PHASE_MATRIX%KERNEL%distname_N)
       
      CALL OPTCHAR_bin (                                                             &
                        key, DLSF%keyEL, keySUB, keyLS, key_RD1, key_RD,             &						 
                        NWL,WL,                                                      &
                        PHASE_MATRIX%DLSIN%RN_bin,   PHASE_MATRIX%DLSIN%RK_bin,      &
						 
                        KN,                          PHASE_MATRIX%DLSIN%grid,        &
                        SD_bin,                                                      &
						 
                        KR,                          PHASE_MATRIX%DLSIN%R,           &
                                                     PHASE_MATRIX%DLSIN%RD,          &
						 
                        KM, ANG,                                                     &

                        PHASE_MATRIX%DLSIN%xext_bin, PHASE_MATRIX%DLSIN%xabs_bin,    &
                        PHASE_MATRIX%DLSIN%xsca_bin, PHASE_MATRIX%DLSIN%albedo_bin,  & 

                        PHASE_MATRIX%DLSIN%f11_bin,  PHASE_MATRIX%DLSIN%f12_bin,     &
                        PHASE_MATRIX%DLSIN%f22_bin,  PHASE_MATRIX%DLSIN%f33_bin,     &
                        PHASE_MATRIX%DLSIN%f34_bin,  PHASE_MATRIX%DLSIN%f44_bin,     &
						 
                        KC, RC2,                                                     &
                        LB, LE,                                                      &
                        pomin, pomax,                                                &
                        xnmin, xnmax, xkmin, xkmax,                                  &
                        PHASE_MATRIX%KERNEL,         PHASE_MATRIX%XGRID1,            &
                        KERNELS1,KERNELS2                                            &
                       )

     IF(DLSF%keyEL.GT.0) PTP11(1:KM)=PHASE_MATRIX%DLSIN%f11_bin(1:KM,IW1)
		 IF(DLSF%keyEL.GT.1) PTP12(1:KM)=PHASE_MATRIX%DLSIN%f12_bin(1:KM,IW1)
		 IF(DLSF%keyEL.GT.2) PTP22(1:KM)=PHASE_MATRIX%DLSIN%f22_bin(1:KM,IW1)
		 IF(DLSF%keyEL.GT.3) PTP33(1:KM)=PHASE_MATRIX%DLSIN%f33_bin(1:KM,IW1)
		 IF(DLSF%keyEL.GT.4) PTP34(1:KM)=PHASE_MATRIX%DLSIN%f34_bin(1:KM,IW1)
		 IF(DLSF%keyEL.GT.5) PTP44(1:KM)=PHASE_MATRIX%DLSIN%f44_bin(1:KM,IW1)

                           EXT1 = PHASE_MATRIX%DLSIN%xext_bin(IW1)
                           !ABS  = PHASE_MATRIX%DLSIN%xabs_bin(IW1)
                           SSA  = PHASE_MATRIX%DLSIN%albedo_bin(IW1)
                           !WRITE(*,*) SSA,EXT1,' - SSA,EXT in sub PHASE_KERNEL_W'

      RETURN 
      END SUBROUTINE PHASE_KERNEL_W

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

