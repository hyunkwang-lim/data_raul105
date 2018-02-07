! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **
#include "../../../constants_set/mod_globals.inc"
      subroutine lidar_garrlic ( HOBS_km, HGR_km, HMAX_atm_km, & ! IN
                                 NHVP_fit, HVP_fit_km,         &
                                 NHVP_retr, HVP_retr_km,       &
                                 AVP_retr, CL,                 &
                                 NSD, EXTA, LRA, EXTM,         &
                                 IFMP, MPROF,                  &
                                 mol_prof_type, aer_prof_type, & ! AL
                                 meas_type, meas               & ! INOUT
                               )

      use mod_par_inv, only : KVERTM, KNBVM ! KVERTM = KVERTM + 2
      use mod_par_OS,  only : KSD, NBVM
      use mod_sdata,    only : meas_type_LS,meas_type_DP,meas_type_RL
      use mod_molecular_scattering, only : std_atm_density
      use mod_os
      use mod_intrpl_linear
      use mod_stop_report

      implicit none
!-----------------------------------------------------------------------------------------
! IN
      real,                       intent(in)    :: HOBS_km, HGR_km, HMAX_atm_km
      integer,                    intent(in)    :: NHVP_fit, NHVP_retr
      real,dimension(KVERTM),     intent(in)    :: HVP_fit_km, MPROF
      real,dimension(KVERTM),     intent(in)    :: HVP_retr_km
      real,dimension(KVERTM,KSD), intent(in)    :: AVP_retr
      real,                       intent(in)    :: CL, EXTM
      real,dimension(KSD),        intent(in)    :: EXTA, LRA
      integer,                    intent(in)    :: NSD, IFMP
      integer,                    intent(in)    :: meas_type
      integer,                    intent(in)    :: mol_prof_type, aer_prof_type
!-----------------------------------------------------------------------------------------
! OUT
      real,dimension(KVERTM),     intent(out)    :: meas
!-----------------------------------------------------------------------------------------
! LOCAL
      integer :: isd, ivert
      integer :: NH
      real,dimension(KVERTM) :: H_km, prof_temp
      real :: xnorm
      integer :: NH_temp
      real,dimension(KVERTM) :: HVP_temp_km
      real,dimension(KVERTM,KSD) :: AVP_temp_norm
      integer :: IFMP1
      real,dimension(KVERTM) :: MVP_norm
      real :: hm_km, sigma
      character(len=20) :: distr_type
!-----------------------------------------------------------------------------------------
! NOTE: Altitude arrays are in descending order

      if(NHVP_retr .lt. NHVP_fit) then
         write(tmp_message,'(2(a,i0),a)') &
         'NHVP_retr = ',NHVP_retr,' .lt. NHVP_fit = ',NHVP_fit,'  not supported'
         G_ERROR(trim(tmp_message))
      endif

      NH = 0
      H_km(:) = 0.0
      prof_temp(:) = 0.0
      xnorm = 0.0
      NH_temp = 0
      HVP_temp_km(:) = 0.0
      AVP_temp_norm(:,:) = 0.0
      hm_km = 0.0
      sigma = 0.0
      distr_type = ''

! Number of altitudes and altitude values for lidar signal computation
      !   include HGR_km altitude in array HVP_temp_km of altitudes for 
      !   normalization of vertical distribution (if it is not present
      !   in HVP_fit_km)
      NH_temp = NHVP_fit
      HVP_temp_km(1:NHVP_fit) = HVP_fit_km(1:NHVP_fit)
      if(HVP_fit_km(NHVP_fit) .ne. HGR_km) then
        NH_temp = NH_temp + 1
        HVP_temp_km(NHVP_fit+1) = HGR_km
      endif

! Aerosol vertical distribution
      !   include HGR_km and HMAX_atm_km altitudes in array H_km of
      !   altitudes for normalization of vertical distribution (if they 
      !   are not present in HVP_fit_km)
      call grid_altitudes_LUT ( HGR_km, HMAX_atm_km,                 & ! IN
                                NHVP_retr, HVP_retr_km(1:NHVP_retr), &
                                NH, H_km(1:NHVP_retr+2)              & ! OUT
                              )
      !   normalization of vertical distribution from ground (HGR_km) level
      !   to top of atmosphere (HMAX_atm_km)
      distr_type = 'lut'
      do isd=1,NSD
        call discrvd_single_atm_component ( distr_type, hm_km, sigma,            & ! IN
                                            NHVP_retr, HVP_retr_km(1:NHVP_retr), &
                                            AVP_retr(1:NHVP_retr,isd),           &
                                            NH, H_km(1:NH),                      &
                                            xnorm, prof_temp(1:NH)               & ! OUT
                                           )
        if ( stop_report%status ) return
        prof_temp(1:NH) = prof_temp(1:NH)/xnorm
        !   normalized aerosol vertical distribution at altidudes of
        !   lidar signal measurements + at HGR_km
        do ivert=1,NH_temp
          AVP_temp_norm(ivert,isd) = &
          LINEAR ( H_km(1:NH), prof_temp(1:NH), NH, HVP_temp_km(ivert) )
        enddo ! ivert
        !write(*,*) isd,xnorm,'  - isd, xnorm'
        !write(*,'(a,i0,a)') 'isd = ',isd,'  HVP_temp_km:'
        !write(*,'(10f14.5)') HVP_temp_km(1:NH_temp)
        !write(*,'(a,i0,a)') 'isd = ',isd,'  AVP_norm:'
        !write(*,'(10e14.5)') AVP_temp_norm(1:NH_temp,isd)
      enddo ! isd

! Molecular vertical distribution/profile
      !IFMP1 = IFMP
      !   TODO for IFMP1 = IFMP :
      !AL add extrapolation at ground level altitude of molecular profile loaded from SData
      !AL change discr_vertical_distribution_single to have a switch between nearest neighbour (hardcoded now)
      !AL and other extrapolation methods, for example linear
      IFMP1 = 0

      select case(IFMP1)
      case(1) ! from SDATA file
      ! NOTE: Development needed
        !   linear extrapolation of molecular backscatter profile provided
        !   with lidar signal measurements
        !MPROF(NHVP_fit+1) = MPROF(NHVP_fit) +   &
        !                    (HVP_retr_km(NHVP_fit+1)-HVP_retr_km(NHVP_fit)) /  &
        !                    (HVP_retr_km(NHVP_fit-1)-HVP_retr_km(NHVP_fit)) *  &
        !                    (MPROF(NHVP_fit-1)-MPROF(NHVP_fit))
      case(0) ! molecular scattering vertical distribution
        !   include HGR_km and HMAX_atm_km altitudes in array H_km of
        !   altitudes for normalization of vertical distribution (if they
        !   are not present in HVP_temp_km)
        call grid_altitudes_LUT ( HGR_km, HMAX_atm_km,             & ! IN
                                  NH_temp, HVP_temp_km(1:NH_temp), &
                                  NH, H_km(1:NH_temp+2)            & ! OUT
                                )
        select case(mol_prof_type)
        case(1)
          distr_type = 'exponential'
          hm_km = 8.0
        case(2)
          distr_type = 'stdatm'
        end select
        !   vertical distribution from ground (HGR_km) level
        !   to top of atmosphere (HMAX_atm_km)
        call discrvd_single_atm_component ( distr_type, hm_km, sigma, & ! IN
                                            NH, H_km(1:NH),           &
                                            prof_temp(1:NH),          &
                                            NH, H_km(1:NH),           &
                                            xnorm, prof_temp(1:NH)    & ! OUT
                                           )
        if ( stop_report%status ) return
        prof_temp(1:NH) = prof_temp(1:NH)/xnorm
        do ivert=1,NH_temp
          MVP_norm(ivert) = &
          LINEAR ( H_km(1:NH), prof_temp(1:NH), NH, HVP_temp_km(ivert) )
        enddo ! ivert
        !write(*,*) xnorm,'  - xnorm'
        !write(*,'(a,i0,a)') 'isd = ',isd,'  HVP_temp_km:'
        !write(*,'(10f14.5)') HVP_temp_km(1:NH_temp)
        !write(*,'(a,i0,a)') 'isd = ',isd,'  MVP_norm:'
        !write(*,'(10e14.5)') MVP_norm(1:NH_temp)
      end select

! Lidar signal computation
      meas(:)=0.0
      if(meas_type .eq. meas_type_LS) then
        call lidar_signal_elastic ( HOBS_km, NH_temp, HVP_temp_km,  &
                                    NSD, AVP_temp_norm, EXTA, LRA,  &
                                    MVP_norm, EXTM,                 &
                                    meas(:) )
        if ( stop_report%status ) return
        meas(1:NHVP_fit) = meas(1:NHVP_fit)*CL*0.001 ! units coefficient 0.001 meters to kilometers
        !do ivert=1,NHVP_fit
        !  write(*,*) ivert,HVP_fit_km(ivert),meas(ivert),'  ivert,HVP_fit_km, meas'
        !enddo
      elseif(meas_type .eq. meas_type_DP) then
        stop 'stop in lidar_garrlic: Depolarization is under development'
      endif

      return
      end subroutine lidar_garrlic

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



