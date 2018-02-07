! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **


!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
#include "../constants_set/mod_globals.inc"
      subroutine validator_constants ( RIN )

      use mod_par_inv
      use mod_par_OS
      use mod_par_DLS_bin, only : KWLpar, KCpar
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) ::  RIN
      integer ::  i, j, IDIM1, IDIM2, par_type
! -----------------------------------------------------------------------------------------
      if(KPARS .gt. KDIF) then
        write(tmp_message,'(2(a,i0),3x,a,a,a)' ) &
        'KPARS = ',KPARS,' .gt. KDIF = ',KDIF,'KDIF is used in smoothness subroutines', &
        NEW_LINE('A'), &
        'Check inversion constants: KDIF = max(KPARS,KITIME,KIX,KIY).'
        G_ERROR(trim(tmp_message))
      endif
      if (KVERTM .gt. KNBVM) then
        write(tmp_message,'(3(a,i0))') &
        'KVERTM = ',KVERTM,' .gt. KNBVM = ',KNBVM, &
        ' not valid combination (KNBVM = max(KVERTM,NBVM)) NBVM = ',NBVM
        G_ERROR(trim(tmp_message))
      endif

      if( RIN%use_models .eqv. .true. ) then
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
            if(par_type .eq. par_type_SD_LB) then
              do IDIM2=1,RIN%NDIM%n2(IDIM1)
                if(RIN%NDIM%n3(IDIM2,IDIM1) .gt. KCpar) then
                  write(tmp_message,'(2(a,i0,3x),a,2(a,i0,3x),a)') &
                  'IDIM1 = ',IDIM1,'IDIM2 = ',IDIM2,NEW_LINE('A'), &
                  'Number of models = ',RIN%NDIM%n3(IDIM2,IDIM1), &
                  'can not be larger than constant KCpar = ',KCpar,' in mod_par_DLS_bin.'
                  G_ERROR(trim(tmp_message))
                endif
              enddo
              if(RIN%NW .gt. KWLpar) then
                write(tmp_message,'(2(a,i0,2x),a)') &
                'Number of wavelengths = ',RIN%NW, &
                'in settings can not be larger than constant KWLpar = ',KWLpar, &
                'in mod_par_DLS_bin.f90'
                G_ERROR(trim(tmp_message))
              endif
            endif
        enddo ! IDIM1
      endif

      return
      end subroutine validator_constants

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine validator_SOS_RT_parameters ( RIN )

      use mod_par_OS
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) ::  RIN
      !type(stop_report_type),     intent(inout) ::  stop_report
      integer ::  i,j
! -----------------------------------------------------------------------------------------
      if(RIN%NLYRS+1 .gt. KNT) then
        write(tmp_message,'(2(a,i0),3a)') &
        'RIN%NLYRS+1 = ',RIN%NLYRS+1,' .gt. KNT = ',KNT,'  in validator_SOS_RT_parameters', &
        NEW_LINE('A'), &
        'Check settings and parameters in both Settings and mod_par_OS'
        G_ERROR(trim(tmp_message))
      endif

      return
      end subroutine validator_SOS_RT_parameters

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

      subroutine validator_inversion_settings ( RIN )

      use mod_par_inv
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) ::  RIN
      !type(stop_report_type),     intent(inout) ::  stop_report
      integer ::  i, j, nwmax1, nwmax2, IDIM1, IDIM2, par_type, nbins
      logical ::  lpresent
! -----------------------------------------------------------------------------------------
! Validate KNSING,KNSINGF - number of parameters driving forwar model and 
!                        number of retrieved parameters 
      if(RIN%KNSING .gt. KPARS) then
        write(tmp_message,'(2(a,i4),2a)') 'KNSING = ',RIN%KNSING,' .GT. KPARS = ',KPARS, &
        NEW_LINE('A'), &
        'KPARS in module mod_par_inv'
        G_ERROR(trim(tmp_message))
      endif ! KNSING.GT.KPARS
      if(RIN%KNSINGF .gt. RIN%KNSING) then
        write(tmp_message,'(2(a,i0))') 'KNSINGF = ',RIN%KNSINGF,' .GT. KNSING = ',RIN%KNSING
        G_ERROR(trim(tmp_message))
      endif ! KNSING .GT. KPARS
! Validate retrieved parametres order 
      do IDIM1=2,RIN%NDIM%n1
        if(RIN%NDIM%par_retr(IDIM1)) then
          if(.not. RIN%NDIM%par_retr(IDIM1-1)) then
            write(tmp_message,'(2(a,i0))') 'Retrieved characteristic # ',IDIM1, &
            ' can not follow not retrived characteristic # ',IDIM1-1
            G_ERROR(trim(tmp_message))
          endif ! .not. RIN%NDIM%par_retr(IDIM1-1)
        endif ! RIN%NDIM%par_retr(IDIM1)
      enddo
! Validate NW - number of wavelengths for refractive index and surface
      nwmax1  = 0
      nwmax2  = 0
      do IDIM1=1,RIN%NDIM%n1 
      if(RIN%NDIM%par_type(IDIM1) .eq. par_type_RERI_spect) then  
        do IDIM2=1,RIN%NDIM%n2(IDIM1)
          nwmax1 = max( nwmax1,RIN%NDIM%n3(IDIM2,IDIM1) )
        enddo
      endif ! RIN%NDIM%par_type(IDIM1) .eq. par_type_RERI_spect
      if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and. & 
      RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF_water_end) then
        do IDIM2=1,RIN%NDIM%n2(IDIM1)
          nwmax2 = max( nwmax2,RIN%NDIM%n3(IDIM2,IDIM1) )
        enddo       
      endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and.
      enddo ! IDIM1
      if(nwmax1 .ne. 0 .and. nwmax2 .ne. 0) then
        if(nwmax1 .ne. nwmax2) then
          write(tmp_message,'(2(a,i0))') 'NW for refr.index = ',nwmax1,' .ne. NW for surface = ',nwmax2
          G_ERROR(trim(tmp_message))
        endif ! nwmax1 .ne. nwmax2
        if(nwmax1 .gt. KW) then
          write(tmp_message,'(2(a,i0))') 'NW = ',nwmax1,'  KW = ',KW,'  KW in module mod_par_inv'
          G_ERROR(trim(tmp_message))
        endif ! nwmax1 .gt. KW
      endif ! nwmax1 .ne. 0 .and.

!! Check keyEL and presence of polarized surface
!      do IDIM1=1,RIN%NDIM%n1
!        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF2_land_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF2_land_end) then
!          if(RIN%DLSF%keyEL .lt. 4) then
!            write(tmp_message,'(a,i4,2a)') 'keyEL=',RIN%DLSF%keyEL, &
!            NEW_LINE('A'), &
!            'Number of elements_of_phase_matrix has to be set to 4 if polarized surface is present.'
!            G_ERROR(trim(tmp_message))
!          endif
!        endif ! RIN%NDIM%par_type(IDIM1) .gt.
!      enddo ! IDIM1

! Validate IO - difference order for single pixel a priori estimate constraints
      if( any(RIN%SPCA%IO(:,:,:) .ne. 0) ) then
        write(tmp_message,'(a)') 'RIN%SPCA%IO .ne. 0', &
        ' IO is difference order for single pixel a priori estimate constraints'
        G_ERROR(trim(tmp_message))
      endif ! RIN%iPOBS .lt. 1 .or.

! Validate iPOBS - fit polarization observation option
      if(RIN%iPOBS .lt. 1 .or. RIN%iPOBS .gt. 5) then
        write(tmp_message,'(a,i0,a)') 'iPOBS = ',RIN%iPOBS,'  - Unknown value ( 1 <= iPOBS <= 5 )'
        G_ERROR(trim(tmp_message))
      endif ! RIN%iPOBS .lt. 1 .or.

! Validate edges parameters
      if( KIEDGE .lt. RIN%edges%nx .or.  &
          KIEDGE .lt. RIN%edges%ny .or.  &
          KIEDGE .lt. RIN%edges%nt )     &
      then
        write(tmp_message,'(6(a,i0,a),3a)') &
        'KIEDGE = ',KIEDGE,'  RIN%edges_max%nx = ',RIN%edges%nx, &
        NEW_LINE('A'), &
        'KIEDGE = ',KIEDGE,'  RIN%edges_max%ny = ',RIN%edges%ny, &
        NEW_LINE('A'), &
        'KIEDGE = ',KIEDGE,'  RIN%edges_max%nt = ',RIN%edges%nt, &
        NEW_LINE('A'), &
        'Constant KIEDGE can not be smaller than', &
        NEW_LINE('A'), &
        'settings values for RIN%edges%nx,RIN%edges%ny,RIN%edges%nt'
        G_ERROR(trim(tmp_message))
      endif

! Validate number of chemical aerosol component
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .eq. par_type_CXRI_nmix) then
          if(RIN%NSD .eq. 1) then
            if(RIN%NDIM%n3(1,IDIM1) .ne. 4) then
              write(tmp_message,'(a,i0,a)') &
              'Number of chem. components nfract = ',RIN%NDIM%n3(1,IDIM1),' .ne. 4 for one component aerosol.'
              G_ERROR(trim(tmp_message))
            endif
          elseif(RIN%NSD .eq. 2) then
            if(RIN%NDIM%n3(1,IDIM1) .ne. 4) then
              write(tmp_message,'(a,i0,a)') &
              'Number of chem. components nfract = ',RIN%NDIM%n3(1,IDIM1),' .ne. 4 for fine mode.'
              G_ERROR(trim(tmp_message))
            endif
            if(RIN%NDIM%n3(2,IDIM1) .ne. 3) then
              write(tmp_message,'(a,i0,a)') &
              'Number of chem. components nfract = ',RIN%NDIM%n3(2,IDIM1),' .ne. 3 for coarse mode.'
              G_ERROR(trim(tmp_message))
            endif
          endif
        exit
        endif ! par_type(IDIM1) .eq. par_type_CXRI_nmix
      enddo ! IDIM1

! Validate number of precomputed lognormal bins
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .eq. par_type_SD_LB) then
          if(RIN%NSD .eq. 2) then
            nbins = sum( RIN%NDIM%n3(1:RIN%NSD,IDIM1) )
            if(nbins .gt. KIDIM3) then
              write(tmp_message, '(2(a,i0,2x),a)') &
              'Number of precomputed lognormal bins nbins = ',nbins, &
              'can not be > than constant KIDIM3 = ',KIDIM3,'in module mod_par_inv'
              G_ERROR(trim(tmp_message))
            endif
          endif
        exit
        endif ! par_type(IDIM1) .eq. par_type_CXRI_nmix
      enddo ! IDIM1

! Validate normal system solver
      if(RIN%IMQ .eq. 2) then
        write(tmp_message,'(a)') 'Normal system solver option singular_value_decomposition', &
        ' is currently not supported.'
        G_ERROR(trim(tmp_message))
      endif

! Validate product flags
!      if(RIN%NSD .eq. 1) then
      if(RIN%products%aerosol%pm .or. RIN%products%aerosol%types) then
      if(.not. RIN%products%aerosol%sd2m_mph .and. &
         .not. RIN%products%aerosol%sd2m_ext) then
        write(tmp_message,'(a,a,2(a,es11.4),a,a,a,2(a,es11.4))') &
        'To provide these products:', &
        NEW_LINE('A'), &
        'products.aerosol.pm =',RIN%products%aerosol%pm, &
        '  products.aerosol.types =',RIN%products%aerosol%types, &
        NEW_LINE('A'), &
        'one of the following product options must be .true.', &
        NEW_LINE('A'), &
        'products.aerosol.sd2m_mph =',RIN%products%aerosol%sd2m_mph, &
        '  products.aerosol.sd2m_ext =',RIN%products%aerosol%sd2m_ext
        G_ERROR(trim(tmp_message))
      endif
      endif
!      endif

      if(RIN%products%forcing%forcing .or. RIN%products%forcing%bbflux) then
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .eq. par_type_SD_LB) then
            write(tmp_message,'(3a)') &
            'Radiative forcing product (products.forcing) can not be provided', &
            NEW_LINE('A'), &
            'for precomputed lognormal bins. Choose another size distribution model.'
            G_ERROR(trim(tmp_message))
          endif ! par_type .gt. par_type_SD_beg .and.
        enddo ! IDIM1
      endif
      if(RIN%products%aerosol%lidar) then
      if(RIN%DLSF%keyEL .lt. 1) then
        write(tmp_message,'(3a)') &
        'Lidar product (lidar ratio) was requested in settings (products.aerosol.lidar = .true.).', &
        NEW_LINE('A'), &
        'Lidar ratio can be provided if', &
        'kernel.phase_matrix_package.elements_of_phase_matrix set => 1 .'
        G_ERROR(trim(tmp_message))
      endif
      endif

! Validate presence of vertical profile parameter (standard deviation)
! if vertical profile is gaussian
      if( any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_AVP_par_height) ) then
      if(RIN%aer_prof_type .eq. 1) then
        lpresent = .false.
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .eq. par_type_AVP_par_std) then
            lpresent = .true.
          exit
          endif ! RIN%NDIM%par_type(IDIM1) .eq.
        enddo ! IDIM1
        if(.not. lpresent) then
          write(tmp_message,'(3a)') 'Aerosol vertical profile standard deviation characteristic must be', &
          NEW_LINE('A'), &
          'provided in settings along with gaussian aerosol vertical profile type.'
          G_ERROR(trim(tmp_message))
        endif
      endif
      endif

! Validate presence of concentration characteristic for Lognormal SD model
      if( any(RIN%NDIM%par_type(1:RIN%NDIM%n1) .eq. par_type_SD_LN) ) then
        lpresent = .false.
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .eq. par_type_Cv) then
            lpresent = .true.
          exit
          endif ! par_type .eq.
        enddo ! IDIM1
        if(.not. lpresent) then
          write(tmp_message,'(3a)') 'For Lognormal size distribution model concentration', &
          NEW_LINE('A'), &
          'characteristic must to be provided in settings.'
          G_ERROR(trim(tmp_message))
        endif
      endif

! Validate model of size disrtribution and use_models flag
      if( RIN%use_models .eqv. .true. ) then
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .gt. par_type_SD_beg  .and. par_type .lt. par_type_SD_end) then
            if(par_type .ne. par_type_SD_LB) then
              write(tmp_message,'(a)') 'Phase matrix models work only with lognormal bin size distribution model.'
              G_ERROR(trim(tmp_message))
            endif
          endif ! par_type .eq.
        enddo ! IDIM1
        !if(RIN%NSD .gt. 1) then
        !  write(tmp_message,'(a))') 'Phase matrix lut is in use. Number of aerosol modes (NSD) can not be > 1.'
        !  G_ERROR(trim(tmp_message))
        !endif
      endif

      return
      end subroutine validator_inversion_settings

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine validator_sdata ( RIN, segment )

      use mod_sdata
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) :: RIN
      type(segment_data),         intent(in) :: segment
      !type(stop_report_type),     intent(inout) ::  stop_report

      integer :: ipix, iw, ip, ip1, npixels, meas_type, meas_type1, nmeas_type
      logical :: lsurf, SvR_meas_present
! -----------------------------------------------------------------------------------------

      npixels = segment%npixels 
! validate surface parameters
      lsurf = SvR_meas_present(segment)

      if(lsurf) then
        if(any(NINT(segment%pixels(1:npixels)%land_percent) .ne. 100) .and. RIN%isurf_water .eq. -999) then
          write(tmp_message,'(3a)') 'Water surface parameters are not provided for pixels with water surface.', &
          NEW_LINE('A'), &
          'Inversion settings and SDATA are inconsistent.'
          G_ERROR(trim(tmp_message))
        endif
        if(any(NINT(segment%pixels(1:npixels)%land_percent) .eq. 100) .and. RIN%isurf_land(1) .eq. -999) then
          write(tmp_message,'(3a)') 'Land surface parameters are not provided for pixels with land surface.', &
          NEW_LINE('A'), &
          'Inversion settings and SDATA are inconsistent.'
          G_ERROR(trim(tmp_message))
        endif
      endif ! lsurf
! validate sets of angles for radiance and polarization
!   loop_pix1: do ipix=1,npixels
!      do iw=1,segment%pixels(ipix)%nwl
!      nmeas_type = segment%pixels(ipix)%meas(iw)%NIP
!      do ip=1,nmeas_type
!        meas_type = segment%pixels(ipix)%meas(iw)%meas_type(ip)
!        if(meas_type .ge. meas_type_SvR_beg .and. meas_type .ge. meas_type_SvR_end) then
!          if(nmeas_type-ip .gt. 0) then
!          meas_type1 = segment%pixels(ipix)%meas(iw)%meas_type(ip+1)
!          if(meas_type1 .ge. meas_type_SvR_beg .and. meas_type1 .ge. meas_type_SvR_end) then

!          endif
!          endif ! nmeas_type-ip .gt. 0
!            exit loop_pix1
!        endif ! meas_type .ge. meas_type_lid_beg .and.
!      enddo ! ip
!      enddo ! iw
!    enddo loop_pix1
! validate sets of angles for p11 and p12

! validate measurement types dublicated in single pixel
      loop_pix: do ipix=1,npixels
      do iw=1,segment%pixels(ipix)%nwl
      nmeas_type = segment%pixels(ipix)%meas(iw)%NIP
      if(nmeas_type .gt. 1) then
        do ip=1,nmeas_type-1
          meas_type = segment%pixels(ipix)%meas(iw)%meas_type(ip)
          do ip1=ip+1,nmeas_type
            if(meas_type .eq. segment%pixels(ipix)%meas(iw)%meas_type(ip1)) then
              write(tmp_message,'(4(a,i0,3x),a,2(a,i0,3x),4a)') 'ipix = ',ipix,'iw = ',iw, &
              'ip = ',ip,'ip1 = ',ip1, &
              NEW_LINE('A'), &
              'meas_type(ip) = ',meas_type, &
              'meas_type(ip1) = ',segment%pixels(ipix)%meas(iw)%meas_type(ip1), &
              NEW_LINE('A'), &
              'Measurements of the same meas type can not be present', &
              NEW_LINE('A'), &
              'for single wavelength in single pixel.'
              G_ERROR(trim(tmp_message))
            endif
          enddo
        enddo ! ip
      endif
      enddo ! iw
      enddo loop_pix

      return
      end subroutine validator_sdata

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine validator_settings_sdata ( RIN, segment )

      use mod_sdata
      use mod_retr_settings_derived_type
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in) :: RIN
      type(segment_data),         intent(in) :: segment
      !type(stop_report_type),     intent(inout) ::  stop_report

      integer :: IDIM1
      integer :: par_type
      logical :: lSvR, SvR_meas_present, lsurf, lvd
! -----------------------------------------------------------------------------------------
      lsurf = .false.
      lvd = .false.

! search if Stokes parameteres are present as measurements
      lSvR = SvR_meas_present(segment)
      if(lSvR) then
! validate if surface characteristics are present in settings
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .gt. par_type_surface_beg .and. par_type .lt. par_type_surface_end) then
            lsurf = .true.
          exit
          endif ! RIN%NDIM%par_type(IDIM1) .eq.
        enddo ! IDIM1
        if(.not. lsurf) then
          write(tmp_message,'(3a)') 'Surface characteristics are not provided for radiative transfer.', &
          NEW_LINE('A'), &
          'Inversion settings and SDATA are inconsistent.'
          G_ERROR(trim(tmp_message))
        endif
! validate if vertical distribution characteristics are present in settings
        do IDIM1=1,RIN%NDIM%n1
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .gt. par_type_AVP_beg .and. par_type .lt. par_type_AVP_end) then
            lvd = .true.
          exit
          endif ! RIN%NDIM%par_type(IDIM1) .eq.
        enddo ! IDIM1
        if(.not. lvd) then
          write(tmp_message,'(3a)') &
          'Vertical distribution characteristics are not provided for radiative transfer.', &
          NEW_LINE('A'), &
          'Inversion settings and SDATA are inconsistent.'
          G_ERROR(trim(tmp_message))
        endif
      endif ! lSvR

      return
      end subroutine validator_settings_sdata

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      logical function SvR_meas_present ( segment )

      use mod_sdata
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------------------------
      type(segment_data), intent(in) :: segment

      integer :: ipix, iw, ip, npixels, meas_type, nmeas_type
! -----------------------------------------------------------------------------------------
      SvR_meas_present = .false.

      npixels = segment%npixels 
! search if Stokes parameteres are present as measurements
      loop_pix: do ipix=1,npixels
      do iw=1,segment%pixels(ipix)%nwl
      nmeas_type = segment%pixels(ipix)%meas(iw)%NIP
      do ip=1,nmeas_type
        meas_type = segment%pixels(ipix)%meas(iw)%meas_type(ip)
        if(meas_type .gt. meas_type_SvR_beg .and. meas_type .lt. meas_type_SvR_end) then
            SvR_meas_present = .true.
            return
        endif ! meas_type .gt. meas_type_lid_beg .and.
      enddo ! ip
      enddo ! iw
      enddo loop_pix

      return
      end function SvR_meas_present

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



