        !> @file vertical_distr.f90
        !> Contains routines related to computation of discret vertical
        !> distributions for given number of altitudes to be used in 
        !> in Successive Order of Scattering (SOS) radiative transfer routine
        !>

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine prepares input and calls routine computing discret
        !> @brief verical distributions for all assumed atmospheric components
        !>
        !> @param[in]  RIN - input settings structure
        !> @param[in]  HGR_km - ground height above sea level
        !> @param[in]  HMAX_atm_km - maximum height of atmosphere
        !> @param[in]  NHVP_retr - number of altitudes for retrieved vertical profile
        !> @param[in]  HVP_retr_km - altitudes for retrieved vertical profile
        !> @param[in]  H0 - parameters of vertical distribution or distribution itself
        !> @param[in]  sigma_aerosol - standard deviation for normal distribution
        !> @param[in]  ifgas - defines presence of gas absorption in sdata structure
        !> @param[in]  gaspar - gas absorption value from sdata structure
        !> @param[out] DISCRVD - discret vertical distribution structure

      subroutine RT_discret_vertical_distr ( RIN, HGR_km, HMAX_atm_km,     &
                                             NHVP_retr, HVP_retr_km, H0,   &
                                             sigma_aerosol, ifgas, gaspar, &
                                             DISCRVD )

      use mod_retr_settings_derived_type
      use mod_vertical_distr_derived_type
      use mod_par_inv, only : KVERTM
      use mod_par_OS, only : KSD

      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings), intent(in) :: RIN
      real, intent(in) :: HGR_km, HMAX_atm_km
      integer, intent(in) :: NHVP_retr
      real,dimension(KVERTM), intent(in) :: HVP_retr_km
      real,dimension(KVERTM,KSD), intent(in) :: H0 ! AL contains parameters of aerosol vertical
                                                   ! distribution or distribution itself
      real,dimension(KSD), intent(in) :: sigma_aerosol
      integer, intent(in) :: ifgas
      real, intent(in) :: gaspar
      type(discret_vertical_distribution), intent(inout) :: DISCRVD
! -----------------------------------------------------------------------------------------
      character(len=20), dimension(NMM+NMG) :: distr_type
      real, dimension(NMM+NMG) :: hm_km, sigma_km
      integer, dimension(NMM+NMG) :: nhinp
      real, dimension(KVERTM,NMM+NMG) :: hinp_km
      real, dimension(KVERTM,NMM+NMG) :: vdinp
      integer :: nhinp_max
! -----------------------------------------------------------------------------------------
      integer :: natm, naer, ncld, nmol, ngas
      integer :: i, ibeg, iend
! -----------------------------------------------------------------------------------------
      nhinp_max = KVERTM

      distr_type(:) = ''
      hm_km(:) = -999.0
      sigma_km(:) = -999.0
      nhinp(:) = -999
      hinp_km(:,:) = -999.0
      vdinp(:,:) = -999.0
      natm = 0
      naer = 0
      ncld = 0
      nmol = 0
      ngas = 0
! Aerosol
      naer = RIN%NSD
      natm = natm + naer
      if(naer .gt. 0) then
        ibeg = 1
        iend = natm
        if(NHVP_retr .eq. 1) then
          select case(RIN%aer_prof_type)
          case(0) ! exponential
          ! exponential
            distr_type(ibeg:iend) = 'exponential'
            hm_km(ibeg:iend) = H0(1,1:naer)*0.001  ! km
          case(1)
          ! gaussian
            distr_type(ibeg:iend) = 'gaussian'
            hm_km(ibeg:iend) = H0(1,1:naer)*0.001  ! km
            sigma_km(ibeg:iend) = sigma_aerosol(1:naer)*0.001 ! km   ! 530.33m -> 0.53033km
          end select
        else
          ! lookup table
          distr_type(ibeg:iend) = 'lut'
          nhinp(ibeg:iend) = NHVP_retr
          vdinp(1:NHVP_retr,ibeg:iend) = H0(1:NHVP_retr,1:naer)
          do i=ibeg,iend
          hinp_km(1:NHVP_retr,i) = HVP_retr_km(1:NHVP_retr)
          enddo
        endif
      endif

! Clouds
!      ncld = RIN%NSD_cloud
!      natm = natm + ncld
!      if(ncld .gt. 0) then
!        ibeg = naer + 1
!        iend = natm
!        if(NHVP_retr .eq. 1) then
!          select case(RIN%cld_prof_type)
!          case(1) ! gaussian
!            distr_type(ibeg:iend) = 'gaussian'
!            hm_km(ibeg:iend) = H0_cloud(1,1:ncld)*0.001  ! km
!            sigma_km(ibeg:iend) = sigma_cloud(1:ncld)*0.001 ! km   ! 530.33m -> 0.53033km
!          case(2) ! exponential
!            distr_type(ibeg:iend) = 'exponential'
!            hm_km(ibeg:iend) = H0_cloud(1,1:ncld)*0.001  ! km
!          end select
!        else
!          distr_type(ibeg:iend) = 'lut'
!          nhinp(ibeg:iend) = NHVP_retr
!          hinp_km(1:NHVP_retr,ibeg:iend) = HVP_retr_km(1:NHVP_retr)
!          vdinp(1:NHVP_retr,ibeg:iend) = H0_cloud(1:NHVP_retr,1:ncld)
!        endif
!      endif

! Molecular scattering
      nmol = 1
      natm = natm + nmol
        ibeg = naer + ncld + 1
        iend = natm
        select case(RIN%mol_prof_type)
        case(0) ! exponential
            distr_type(ibeg:iend) = 'exponential'
            hm_km(ibeg:iend) = 8.0
        case(1) ! standard atmosphere
            distr_type(ibeg:iend) = 'stdatm'
        end select

! Gas absorption ! do not delete, future development
!      ngas = 0
!      gasabs = 0.0
!      if(ifgas .eq. 1) then
!            ngas = 1
!            natm = natm + ngas
!            ibeg = naer + ncld + nmol + 1
!            iend = naer + ncld + nmol + ngas
!            gasabs = gaspar
!            distr_type(ibeg:iend) = 'stdatm'
!      endif

      DISCRVD%natm = natm
      DISCRVD%naer = naer
      DISCRVD%ncld = ncld
      DISCRVD%nmol = nmol
      DISCRVD%ngas = ngas

      !write(*,'(5(a,i0,2x))') 'natm = ',natm,'naer = ',naer, &
      !                     'ncld = ',ncld,'nmol = ',nmol,'ngas = ',ngas

      !if(RIN%IPRI_verbose) then
      !  if(nhinp_max .lt. maxval(nhinp(1:natm))) then
      !    write(*,'(2(a,i0,x))') 'nhinp_max = ',nhinp_max, &
      !    'must be greater or equal to maxval(nhinp(1:natm)) = ',maxval(nhinp(1:natm))
      !    write(*,'(a)') 'Execution has to be terminated.'
      !    stop 'stop in RT_discret_vertical_distr'
      !  endif
      !endif
      ! nhinp_max = KVERTM
!write(*,*) 'RT_discret_vertical_distr begin '
      call discret_vertical_distr ( HGR_km, HMAX_atm_km, distr_type,  &
                                    hm_km, sigma_km,                  &
                                    nhinp_max, nhinp(:),              &
                                    hinp_km(1:nhinp_max,:),           &
                                    vdinp(1:nhinp_max,:),             &
                                    DISCRVD )

!write(*,*) 'RT_discret_vertical_distr end '
      return
      end subroutine RT_discret_vertical_distr

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine computes discret verical distributions for
        !> @brief all assumed atmospheric components
        !>
        !> @param[in]  HGR_km - ground height above sea level in km
        !> @param[in]  HMAX_atm_km - maximum height of atmosphere in km
        !> @param[in]  distr_type - extinction vertical distributions
        !> @param[in]  hm_km - mean/scale height for normal distribution
        !> @param[in]  sigma_km - standard deviation for normal distribution
        !> @param[in]  nhinp - number of altitudes for input LUT profile
        !> @param[in]  hinp_km - altitudes for input LUT profile
        !> @param[in]  vdinp - input LUT profile
        !> @param[out] DISCRVD - discret vertical distribution structure

      subroutine discret_vertical_distr ( HGR_km, HMAX_atm_km, distr_type,  &
                                          hm_km, sigma_km,                  & ! model (input)
                                          nhinp_max, nhinp, hinp_km, vdinp, & ! lut   (input)
                                          DISCRVD )

      use mod_vertical_distr_derived_type
      use mod_stop_report

      implicit none
!	----------------------------------------------------------------------
      real, intent(in) :: HGR_km, HMAX_atm_km
      character(len=20), dimension(NMM+NMG) :: distr_type
      real, intent(in), dimension(NMM+NMG) :: hm_km, sigma_km
      integer, intent(in) :: nhinp_max
      integer, intent(in), dimension(NMM+NMG) :: nhinp
      real, intent(in), dimension(nhinp_max,NMM+NMG) :: hinp_km
      real, intent(in), dimension(nhinp_max,NMM+NMG) :: vdinp
      type(discret_vertical_distribution), intent(inout) :: DISCRVD
!	----------------------------------------------------------------------
      integer :: i, i1, ilut
      real,dimension(KVERT_WD+2) :: h_km_temp
      integer :: natm
      logical :: lut_present
!	----------------------------------------------------------------------
      lut_present = .false.
! distr_type: const / exponential / gaussian / lut / stdatm
      natm = DISCRVD%natm

!XH   if we have LUT of vertical profile, the computation will be based on the grid points in the LUT

!XH      find the index for the component using the LUT
!XH      search from the end to prevent conflicts when absorbing gas profile and lidar profile coexist
      ilut = -999
      do i1=natm,1,-1
        if(trim(distr_type(i1)) .eq. 'lut') then
          lut_present = .true.
          ilut = i1
          exit
        endif
      end do
      if (lut_present) then
!XH     grid points before regriding
        !write(*,'(2(a,f10.4))') 'HGR_km =',HGR_km,'  HMAX_atm_km =',HMAX_atm_km
        !write(*,'(a)') 'in discret_vertical_distr before grid_altitudes_LUT'
        !write(*,'(a,i0,a,10(i0,2x))') 'natm = ',natm,'  nhinp: ',nhinp(1:natm)
        !do i1=1,natm
        !do i=1,nhinp(i1)
          !write(*,'(2(2x,i0),f12.4,e12.4,2x,a)')  i1, i, hinp_km(i,i1), vdinp(i,i1), &
          !                                     '- i1, i, hinp_km(i,i1), vdinp(i,i1)'
        !end do
        !end do

        ! creating altitude grid from lidar measurements adding additional points at ground
        ! level (HGR_km) and at top pf the atmosphere HMAX_atm_km) km
        if(nhinp(ilut) .gt. KVERT_WD) then
          write(*,'(2(a,i0))') 'Number of altitudes for given vertical distribution lut = ',nhinp(ilut), &
          '  ilut = ',ilut
          write(*,'(a,i0,a)') 'can not be larger than parameter KVERT_WD = ',KVERT_WD, &
                              ' for discret vertical distribution arrays.'
          stop 'stop in discret_vertical_distr'
        endif
        call grid_altitudes_LUT ( HGR_km, HMAX_atm_km,                      & ! IN
                                  nhinp(ilut), hinp_km(1:nhinp(ilut),ilut), &
                                  DISCRVD%nh,  h_km_temp(1:nhinp(ilut)+2)   & ! OUT
                                )
        if(nhinp(ilut) .eq. KVERT_WD .and. nhinp(ilut) .lt. DISCRVD%nh) then
          write(*,'(2(a,i0))') 'in discret_vertical_distr  nhinp = ',nhinp(ilut), &
                               ' .le. DISCRVD%nh = ',DISCRVD%nh
          write(*,'(a)') 'DISCRVD%h_km() and h_km_temp() array size conflict.'
          stop 'stop in discret_vertical_distr'
        endif

        DISCRVD%h_km(1:DISCRVD%nh) = h_km_temp(1:DISCRVD%nh)
      else
         ! creating a full-scale altitude grid from ground level (HGR_km) to the top of
         ! the atmosphere (HMAX_atm_km)
         DISCRVD%nh = KVERT_WD
         call grid_altitudes ( HGR_km, HMAX_atm_km,        & ! IN
                               DISCRVD%nh,                 &
                               DISCRVD%h_km(1:DISCRVD%nh)  & ! OUT
                             )
      end if
      !write(*,'(a)') 'in discret_vertical_distr altitudes for discret vertical distributions'
      !do i=1,DISCRVD%nh
      !  write(*,'(2x,i0,f12.4,2x,a)') i, DISCRVD%h_km(i),'- i, DISCRVD%h_km(i)'
      !end do

      do i1=1,natm
         ! calculate vertical distributions for each atmospheric compinent 
         ! depending on distribution type
         call discrvd_single_atm_component (  distr_type(i1), hm_km(i1), sigma_km(i1),       &
                                              nhinp(i1), hinp_km(1:nhinp(i1),i1), vdinp(1:nhinp(i1),i1), &
                                              DISCRVD%nh, DISCRVD%h_km(1:DISCRVD%nh),        &
                                              DISCRVD%norm(i1), DISCRVD%val(1:DISCRVD%nh,i1) & ! OUT
                                           )
         !write(*,*) 'in discret_vertical_distr:'
         !write(*,*) i1,DISCRVD%norm(i1),trim(distr_type(i1)),'  - atm component, vd norm, distr type'
         !write(*,'(10f16.7)') DISCRVD%h_km(1:DISCRVD%nh)
         !write(*,'(10e16.7)') DISCRVD%val(1:DISCRVD%nh,i1)
      end do ! i1

      return
      end subroutine discret_vertical_distr

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine computes verical distributions at discret
        !> @brief altitudes for given single atmospheric component
        !>
        !> @param[in]  distr_type - type of vertical distribution
        !> @param[in]  hm_km - mean/scale height for normal/exponential distribution
        !> @param[in]  sigma_km - standard deviation for normal distribution
        !> @param[in]  nhinp - number of altitudes for input LUT distribution
        !> @param[in]  hinp_km - altitudes for input LUT distribution
        !> @param[in]  vdinp - input LUT distribution
        !> @param[in]  nhout - number of altitudes for discret distribution
        !> @param[in]  hout_km - altitudes for discret distribution
        !> @param[out] xnorm - norm of output distribution
        !> @param[out] vdout - discret vertical distribution

      subroutine discrvd_single_atm_component ( distr_type, hm_km, sigma_km, & ! IN
                                                nhinp, hinp_km, vdinp,    &
                                                nhout, hout_km,           &
                                                xnorm, vdout              & ! OUT
                                              )
      use mod_intrpl_linear
      use mod_molecular_scattering, only: std_atm_density
      use mod_stop_report

      implicit none 
!	------------------------------------------------------------------------------------------------------	  
      real, intent(in) :: hm_km, sigma_km
      character(len=20) :: distr_type
      integer, intent(in) :: nhinp
      real, dimension(nhinp), intent(in) :: hinp_km
      real, dimension(nhinp), intent(in) :: vdinp
      integer, intent(in) :: nhout
      real, dimension(nhout), intent(in) :: hout_km
      real, intent(out) :: xnorm
      real, dimension(nhout), intent(out) :: vdout
!	------------------------------------------------------------------------------------------------------
      double precision, dimension(nhout) :: vdtemp
      real, dimension(nhout) :: y
      real, parameter :: tiny = 0.0 !1e-25
      real :: ea, pi
      integer :: i
      logical :: lstop
!	------------------------------------------------------------------------------------------------------
! Calculate vertical distribution (vdout) at nhout given grid altitudes (hout_km)
! NOTE: altitudes are in descending order
!       input : hinp_km   hinp_km(1) > hinp_km(nhinp)
!       output: hout_km   hout_km(1) > hout_km(nhout)
!	------------------------------------------------------------------------------------------------------
! distr_type: const / exponential / gaussian / lut / stdatm
      lstop = .false.

      pi = acos(-1e0)
      vdtemp(:) = 0.0 ! vertical distribution initialization

      !write(*,'(2(a,2x),a)') 'in discrvd_single_atm_component', &
      !trim(distr_type),'profile before modification:'
      !do i=1,nhinp
      !   write(*,'(2x,i0,f12.4,e12.4)') i, hinp_km(i), vdinp(i)
      !end do

      select case(trim(distr_type))
      case('const')
          vdtemp(1:nhout) = 1.
      case('gaussian')
          ! p(h) = exp( -(hg-hm)**2/(2.*sigma*sigma) + (hg-h )**2/(2.*sigma*sigma) ) =
          !      = exp( (hg-h)(hg+h-2hm)/(2.*sigma*sigma) )
          !
          !  hg - hout_km(nhout)   ground hight over sea level
          !  h  - hout_km          current hight
          !  hm - hm_km            mean hight
          ea = 2. * sigma_km * sigma_km
          y(1:nhout) = ( hout_km(nhout)-hout_km(1:nhout) ) * &
                       ( hout_km(nhout)+hout_km(1:nhout)-2.*hm_km )
          vdtemp(1:nhout) = exp ( y(1:nhout)/ea )
      case('exponential')
            ! exponential profile is divided by hm_km so the profile will be normalized
            ! vdtemp(1:nhout) = dexp( dble(hout_km(nhout)-hout_km(1:nhout))/hm_km )/dble(hm_km)
            ! /dble(hm_km) can be deleted because the profile will be normalized
            vdtemp(1:nhout) = dexp( dble(hout_km(nhout)-hout_km(1:nhout))/hm_km )
      case('lut')
         do i=1,nhout
            if (hout_km(i) .gt. hinp_km(1)) then
!XH            ! points higher than maximum height in LUT are assumed 0
               vdtemp(i) = 0.0
               !write(*,*) i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i),
               !           '  1: i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i)'
            else if(hout_km(i) .le. hinp_km(nhinp)) then
!XH            ! points lower than minimum height in LUT assumed to be the same
               vdtemp(i) = vdinp(nhinp)
               !write(*,*) i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i),
               !           '  2: i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i)'
            else
!XH            ! otherwise interpolate to LUT grid
               vdtemp(i) = LINEAR(hinp_km(1:nhinp),vdinp(1:nhinp),nhinp,hout_km(i))
               !write(*,*) i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i),  &
               !           '  3: i0,hinp_km(i0),vdinp(i0),i,hout_km(i),vdtemp(i)'
            end if ! hout_km(i .gt. hinp_km(1)
         enddo ! i=1,nhout
      case('stdatm')
          ! used only for molecular scattering vertical distr.
          call std_atm_density(nhout,hout_km,vdtemp)
      case default
          write(*,'(3a)') 'distr_type = ',trim(distr_type),'  - unknown value'
          call print_crash_report()
          lstop = .true.
      end select ! distr_type
      !write(*,'(2(a,2x),a)') 'in discrvd_single_atm_component modified ',trim(distr_type), &
      !                       'discret profile before normalization:'
      !do i=1,nhout
      !   write(*,'(2x,i0,f12.4,e12.4)') i, hout_km(i), vdtemp(i)
      !end do

! calculate the norm of the profile
        xnorm = 0.0
        do i=2,nhout
          xnorm = xnorm + 0.5*(vdtemp(i-1)+vdtemp(i))*(hout_km(i-1)-hout_km(i))
          !write(*,*) '*** i=',i,'  hout_km(i) =',hout_km(i),'  xnorm =',xnorm
        enddo ! i
        if(xnorm .le. tiny) then
          write(*,'(2a,e12.4)') 'Execution has to be terminated.', &
          ' In discrvd_single_atm_component xnorm =', xnorm
          do i=1,nhout
            write(*,'(2x,i0,f12.4,e12.4,2x,a)') i, hout_km(i), vdtemp(i), &
            '- i, hout_km, vdtemp'
          enddo
          do i=1,nhinp
            write(*,'(2x,i0,f12.4,2x,a)') i, hinp_km(i), &
            '- i, hout_km'
          enddo
          call print_crash_report()
          lstop = .true.
        endif

      if(lstop) &
      stop 'stop in discrvd_single_atm_component'

      !vdout(:) = vdtemp(:)/xnorm
      vdout(:) = vdtemp(:)

      !write(*,'(2(a,2x),a)') 'in discrvd_single_atm_component normalized discret', &
      !                        trim(distr_type),'profile:'
      !write(*,'(a,f10.6)') 'xnorm =',xnorm
      !do i=1,nhout
      !  write(*,'(2x,i0,f12.4,e12.4)') i, hout_km(i), vdout(i)
      !enddo ! i

      return
      end subroutine discrvd_single_atm_component

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine computes given number of equdistant in natural 
        !> @brief logarithm altitudes from top of the atmosphere altitude
        !> @brief to ground height above sea level
        !>
        !> @param[in]  HGR_km - ground height above sea level in km
        !> @param[in]  HMAX_atm_km - top of the atmosphere altitude in km
        !> @param[in]  nhout - number of altitudes
        !> @param[out] hout_km - output altitudes in descending order

      subroutine grid_altitudes ( HGR_km, HMAX_atm_km, &
                                  nhout, hout_km       &
                                )
      implicit none
!	------------------------------------------------------------------------------------------------------      
      real, intent(in) :: HGR_km, HMAX_atm_km
      integer, intent(in) :: nhout
      real, dimension(nhout), intent(out) :: hout_km  ! hout_km(1) > hout_km(nhout)
!	------------------------------------------------------------------------------------------------------
      real :: shift_km, DH_km
      integer :: i
!	------------------------------------------------------------------------------------------------------
      hout_km(:)  = 0.0
      shift_km = 0.0
      if(HGR_km .le. 0.0) then
        shift_km = 1.0
        if(abs(HGR_km) .ge. shift_km) then
          write(*,'(2(a,f12.4),a)') '!!! WARNING !!! abs(HGR_km)=',abs(HGR_km), &
          ' .ge. shift_km=',shift_km,'  in grid_altitudes'
        endif
      endif
            
      hout_km(1) = HMAX_atm_km
      hout_km(nhout) = HGR_km
      DH_km = (log(HMAX_atm_km+shift_km)-log(HGR_km+shift_km))/(nhout-1)
      do i=2,nhout-1
        hout_km(i) = exp(log(hout_km(i-1)+shift_km)-DH_km)-shift_km
      enddo ! i	  	  	  
!      write(*,'(a)') 'in grid_altitudes hout_km:'
!      write(*,'(10f12.4)') hout_km(1:nhout)

      return
      end subroutine grid_altitudes

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine returns altitudes including top of the atmosphere
        !> @brief altitude and ground height above sea level
        !>
        !> @param[in]  HGR_km      - ground height above sea level in km
        !> @param[in]  HMAX_atm_km - top of the atmosphere altitude in km
        !> @param[in]  nhinp       - number of input altitudes
        !> @param[in]  hinp_km     - input altitudes in descending order
        !> @param[out] nhout       - number of output altitudes
        !> @param[out] hout_km     - output altitudes in descending order

      subroutine grid_altitudes_LUT ( HGR_km, HMAX_atm_km, nhinp, hinp_km, &
                                      nhout, hout_km  &
                                    )
      implicit none
!	------------------------------------------------------------------------------------------------------      
      real, intent(in) :: HGR_km, HMAX_atm_km
      integer, intent(in) :: nhinp
      real, dimension(nhinp), intent(in) :: hinp_km ! hinp_km(1) > hinp_km(nhinp)
      integer, intent(out) :: nhout
      real, dimension(nhinp+2), intent(out) :: hout_km  ! hout_km(1) > hout_km(nhout)
!	------------------------------------------------------------------------------------------------------
      nhout = nhinp
      hout_km(1:nhinp) = hinp_km(1:nhinp)

      if(hinp_km(1) .lt. HMAX_atm_km) then
        nhout = nhout+1
        hout_km(2:nhout) = hinp_km(1:nhinp)
        hout_km(1) = HMAX_atm_km
      endif
      if(hinp_km(nhinp) .gt. HGR_km) then
        nhout = nhout+1
        hout_km(nhout) = HGR_km
      endif
!      write(*,'(a)') 'in grid_altitudes_LUT hinp_km:'
!      write(*,'(10f12.4)') hinp_km(1:nhinp)
!      write(*,'(a)') 'in grid_altitudes_LUT hout_km:'
!      write(*,'(10f12.4)') hout_km(1:nhout)

      return
      end subroutine grid_altitudes_LUT
            
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



