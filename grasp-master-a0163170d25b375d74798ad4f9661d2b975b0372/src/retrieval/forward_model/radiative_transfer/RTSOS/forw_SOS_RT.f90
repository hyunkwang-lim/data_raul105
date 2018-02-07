! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

      module mod_rt

      use mod_stop_report

      contains
!
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine get_radiative_transfer_SOS_flags_surf ( RIN, iBRDF_land, iBPDF_land, iBRM_water )

      use mod_retr_settings_derived_type
      
      implicit none
! -----------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in)   ::  RIN
      integer,                    intent(out)  ::  iBRDF_land, iBPDF_land, iBRM_water
      integer   :: IDIM1,par_type	  
! -----------------------------------------------------------------------------------------
! The parameters for radiative transfer module 
!   iBRDF    
!   iBPDF 
! -----------------------------------------------------------------------------------------
! Surface model options

      iBRDF_land = -1
      iBPDF_land = -1
      iBRM_water = -1
      
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_land_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF1_land_end) then
          par_type = RIN%NDIM%par_type(IDIM1)	    	  
          select case (par_type)
          case(par_type_SURF1_land_Ross_Li_BRDF)
            iBRDF_land = 1
          case(par_type_SURF1_land_RPV_BRDF)
            iBRDF_land = 0
          case(par_type_SURF1_land_Litvinov_fast)
            iBRDF_land = 2
!            iBPDF_land = 9	
          case(par_type_SURF1_land_Litvinov)
            iBRDF_land = 9
            iBPDF_land = 9
          end select
        endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_beg .and.

        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF2_land_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF2_land_end) then
          par_type = RIN%NDIM%par_type(IDIM1)        
          select case(par_type)
          case(par_type_SURF2_land_Maignan_Breon)
            iBPDF_land = 0
          case(par_type_SURF2_land_Litvinov)  
            iBPDF_land = 1
          end select
        endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF2_beg .and.

        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF_water_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SURF_water_end) then
          par_type = RIN%NDIM%par_type(IDIM1) 
 	        select case (par_type)
		      case(par_type_SURF_water_Cox_Munk_iso)
		        iBRM_water = 0
		      case(par_type_SURF_water_Cox_Munk_ani)
            iBRM_water = 1
		      case(par_type_SURF_water_Litvinov)
            iBRM_water = 9		  
		      end select	
		    endif ! RIN%NDIM%par_type(IDIM1) .gt. par_type_SURF1_beg .and.
      enddo ! IDIM1

      !write(*,*) 'iBRDF_land,iBPDF_land,iBRM_water:' 
      !write(*,'(3i8,a)') iBRDF_land,iBPDF_land,iBRM_water,'  in get_radiative_transfer_SOS_flags_surf'
     	   	  
      return
      end subroutine get_radiative_transfer_SOS_flags_surf

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      SUBROUTINE forw_SOS_RT (  iFlux,                                   & ! IN
                                RIN,                                     &
                                IW, WAVE, ind_wl, IP,                    &
                                OSHP,                                    &
                                surf_land_par_num, surf_land_par_vect,   &
                                surf_water_par_num, surf_water_par_vect, &
                                tau_mol,                                 &
                                HOBS_km, HGR_km, HMAX_atm_km,            &
                                DISCRVD,                                 &
                                laerosol, lsurface,                      &
                                NANG, ANGL,                              &
                                abs_data_forw_im,                        &
                                pixel_fit,                               &
                                 
                                GOUT_aerosol_opt_pixel_wl,               & ! INOUT
                                GOUT_aerosol_phmx_pixel_wl,              &
                                GOUT_clouds_opt_pixel_wl,                &
                                GOUT_clouds_phmx_pixel_wl,               &
                                GOUT_surface_pixel_wl,                   &
                                GOUT_bbflux_pixel,                       &

                                NBV_comb,                                & ! OUT
                                SLout_comb, SQout_comb,                  &
                                SUout_comb, SLPout_comb                  &
                              )
      use mod_os
      use mod_par_DLS,   only : KMpar
      use mod_par_inv,   only : KPARS, KSHAPE, KBF, KIP, KNBVM
      use mod_par_OS,    only : NMG, KSD, NBVM, NMM, KNT, NG0T, KVERT_WD
      use mod_vertical_distr_derived_type

!XH   modules related to gas absorption
      use mod_abs_kd

      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_sdata
      use mod_inversion_utils, only : check_nan
      use mod_intrpl_linear
      use mod_alloc_kernels
      use mod_stop_report
      use mod_vertical_distr_derived_type

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN:
!XH   switch for broadband flux calculation
      logical,                    intent(in)  ::  iFlux
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  IW,ind_wl,IP 
      real,                       intent(in)  ::  HOBS_km,HGR_km,HMAX_atm_km                            
      real,                       intent(in)  ::  WAVE
      type(OSH_par),              intent(in)  ::  OSHP
      real,                       intent(in)  ::  tau_mol
      type(discret_vertical_distribution), intent(inout) :: DISCRVD
!XH   data of gas absorption
      type (DATA_ABS),            intent(in)  ::  abs_data_forw_im
      type(pixel),                intent(in)  ::  pixel_fit
      integer,                    intent(in)  ::  NANG
      real,dimension(KMpar),      intent(in)  ::  ANGL
      integer,dimension(2),       intent(in)  ::  surf_land_par_num
      real, dimension(2*KBF),     intent(in)  ::  surf_land_par_vect
      integer,                    intent(in)  ::  surf_water_par_num
      real, dimension(KBF),       intent(in)  ::  surf_water_par_vect
!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_aerosol_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_aerosol_phmx_pixel_wl
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_clouds_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_clouds_phmx_pixel_wl
      type(output_pixel_surface_wl),  intent(inout)  ::  GOUT_surface_pixel_wl
      type(output_pixel_bbflux),      intent(inout)  ::  GOUT_bbflux_pixel
!	------------------------------------------------------------------------------------------------------
! OUT:
      integer,                 intent(out) ::  NBV_comb
      real, dimension(2*NBVM), intent(out) ::  SLout_comb, SQout_comb, SUout_comb, SLPout_comb
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer	                      ::  ISD, IV, i1, i2, i3, ii
      real,dimension(2*NBVM)        ::  ASCAT
!	------------------------------------------------------------------------------------------------------
      logical                       ::	status_funct 
!	------------------------------------------------------------------------------------------------------
      integer                       ::  NBV
      real, dimension(NBVM,KIP)     ::  vis,fiv
      real, dimension(NMM+NMG)      ::  EXT_os, SSA_os, EXT0_os
      integer                       ::  ITRONC, iop, NQDR, NT1
      real                          ::  sza, gteta, UFG, DFG,gas
!	------------------------------------------------------------------------------------------------------   
      integer                         ::  natm, naer, ncld, nmol
      integer, dimension(NMM+NMG)     ::  nhinp
      real, dimension(KVERT_WD,NMM+NMG) ::  hinp_km
      real, dimension(KVERT_WD,NMM+NMG) ::  vdinp
      integer :: nhinp_max
      logical ::  laerosol, lsurface
      integer ::  nmeas_type
      integer ::  IDIM1,par_type
!	------------------------------------------------------------------------------------------------------   
! gas absorption
      integer                       ::  ngas, ngas_kd
      integer                       ::  curr_wl_index
      integer,dimension(NCORPS)     ::  ns
      logical                       ::	abs_run
      real,dimension(KNT)           ::  UFX_tmp, DFX_tmp, UFX0_tmp, DFX0_tmp
      character(len=20), dimension(NMM+NMG) :: distr_type
!	------------------------------------------------------------------------------------------------------
      real,dimension(2*NBVM)        ::  SQout_tmp,SUout_tmp,SLPout_tmp,SLout_tmp
      real                          ::  DFG_tmp,UFG_tmp
      real, dimension(2*NBVM)       ::  vis_comb,fiv_comb  
      logical                       ::  lexist
      character(len=255)            ::  external_file_path
      integer                       ::  iBRDF_land,iBPDF_land,iBRM_water
!	------------------------------------------------------------------------------------------------------
! timer
      real,save                     ::  xtime=0.0
      real                          ::  time_start,time_finish
!	------------------------------------------------------------------------------------------------------
!     IW     - index of wave length in array of wave length for current pixel
!     WAVE   - value of wave length for index IW from array of wave length for current pixel
!     ind_wl - index of wave length in general array of wave length for inversion
!     write(*,*) 'in forw_IMAGE_I_IW: iw =',iw,WAVE,ind_wl
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
      nhinp_max = KVERT_WD

      distr_type(:) = 'lut'

      natm = DISCRVD%natm
      naer = DISCRVD%naer
      ncld = DISCRVD%ncld
      nmol = DISCRVD%nmol

      iop = RIN%ipplane
      EXT_os(:) = 0.0
      SSA_os(:) = 0.0
!	------------------------------------------------------------------------------------------------------
! Extinction and single scattering albedo as an input for SOS_RT
!     aerosol
      do ISD=1,naer
        if(RIN%IPRI_verbose) then
        if(GOUT_aerosol_opt_pixel_wl%EXT(ISD) .ge. 70.0) then
          write(*,'(a)') '!!! AOD is too high before call radiative_transfer_SOS !!!'
          write(*,'(a,i0,a,f11.4)') 'isd = ',isd,'  aod = ',GOUT_aerosol_opt_pixel_wl%EXT(ISD)
          call print_crash_report()
          !stop 'stop in forw_SOS_RT'
        endif
        endif
        EXT_os(ISD) = GOUT_aerosol_opt_pixel_wl%EXT(ISD)
        SSA_os(ISD) = GOUT_aerosol_opt_pixel_wl%SSA(ISD)
      enddo ! ISD
!      write(*,*) 'before OS_H: EXT_os(1:naer)  ',EXT_os(1:naer)
!      write(*,*) 'before OS_H: SSA_os(1:naer)  ',SSA_os(1:naer)
!      print *, wave, sum(ext_os(1:naer))
!	------------------------------------------------------------------------------------------------------
!     molecular
      EXT_os(naer+ncld+nmol) = tau_mol
      SSA_os(naer+ncld+nmol) = 1.0
!	------------------------------------------------------------------------------------------------------
!     gas k distribution computation
      ngas_kd = 0
      abs_run = .false.
      if (iFlux) then
!XH      checking if we need absorption calculation for this wavelength
         curr_wl_index = 0
         do ii = 1, abs_data_forw_im%NWL
            if (abs_data_forw_im%ABS_GS_WL(ii)%WVL .eq. WAVE) then
!XH            considering gas absorption for broadband flux calculation
              ngas_kd = 1
              abs_run = .true.
              if (ngas_kd .gt. NMG) then
                write(*,'(2(a,i0))') 'ngas_kd = ',ngas_kd,'  NMG = ',NMG
                write(*,'(a)') 'Number of gas component > NMG parameter.'
                stop 'stop in forw_SOS_RT'
              end if ! ngas_kd .gt. NMG
              curr_wl_index = ii
            exit
            end if
         end do ! ii
      end if
!	------------------------------------------------------------------------------------------------------
!     gas provided in sdata structure
      if(.not. abs_run) then
        if (pixel_fit%IFGAS .eq. 1) then
          gas = pixel_fit%meas(IW)%gaspar
          SSA_os(1) = EXT_os(1)*SSA_os(1)/(EXT_os(1)+gas)
          EXT_os(1) = EXT_os(1)+gas
        end if
      end if
!	------------------------------------------------------------------------------------------------------

      NT1=RIN%NLYRS+1
  
!     Ground based observations (SunPhotometer) or Plane observations upward
      NQDR = NG0T
      if (pixel_fit%HOBS .gt. pixel_fit%MASL) then
!        Satellite or Plane observations downward
         NQDR = OSHP%NG
      else if(pixel_fit%HOBS .lt. pixel_fit%MASL) then
         write(*,*) 'HOBS_km=',HOBS_km,' .lt.  HGR_km=',HGR_km
         stop 'stop in forw_SOS_RT'
      end if
!      NQDR=51

      ITRONC = 1
      if (.not. RIN%ITRONC .or. OSHP%IMSC .eq. 1) ITRONC = 0

      if (iFlux) then
         sza = abs_data_forw_im%SZA
         NBV_COMB = 1
         vis_comb(1) = 0.0
         fiv_comb(1) = 0.0
      else
!         write(*,*) 'in forw_SOS_RT  ip=',ip
         call get_pixel_geom ( IW,IP, pixel_fit, NBV, sza, vis(:,IP), fiv(:,IP) )
!         write(*,*) '11: vis:'
!         write(*,'(10f16.4)') vis(1:NBV,ip)
!         write(*,*) '11: fiv:'
!         write(*,'(10f16.4)') fiv(1:NBV,ip)

         NBV_comb = NBV
         vis_comb(1:NBV) = vis(1:NBV,IP)
         fiv_comb(1:NBV) = fiv(1:NBV,IP)
         
         nmeas_type = pixel_fit%meas(IW)%NIP
         
         if (nmeas_type-IP .gt. 0) then
            call get_pixel_geom (IW,IP+1,      &
                                 pixel_fit,    &
                                 NBV,          &
                                 sza,          &
                                 vis(:,IP+1),  &
                                 fiv(:,IP+1)   &
                                )
!            write(*,*) '12: vis:'
!            write(*,'(10f16.4)') vis(1:NBV,ip+1)

            if ( (NBV .ne. NBV_comb) .or.   &
                 (NBV .eq. NBV_comb .and.   &
                 (any(vis(1:NBV,IP+1) .ne. vis_comb(1:NBV)) .or.   &
                  any(fiv(1:NBV,IP+1) .ne. fiv_comb(1:NBV))) ) ) then
               vis_comb(NBV_comb+1:NBV_comb+NBV) = vis(1:NBV,IP+1)
               fiv_comb(NBV_comb+1:NBV_comb+NBV) = fiv(1:NBV,IP+1)
               NBV_comb = NBV_comb+NBV
            end if ! NBV .ne. NBV_comb .or.
         end if ! nmeas_type-IP .gt. 0
      end if ! iFlux

!      write(*,*) 'sza=',sza
!      write(*,*) '1: vis:'
!      write(*,'(10f16.4)') vis_comb(1:NBV_comb)
!      write(*,*) '1: fiv:'
!      write(*,'(10f16.4)') fiv_comb(1:NBV_comb)

!      write(*,'(a)') 'Vertical distributions:'
!      do iv=1,DISTRVD%nh
!         write(*,'(a,i0,f12.4,10e12.4)') 'iv, h_km(iv), VD(iv,1:natm)  ', &
!         iv,DISTRVD%h_km(iv),DISTRVD%val(iv,1:natm)
!      end do ! iv
!	------------------------------------------------------------------------------------------------------
      call get_radiative_transfer_SOS_flags_surf ( RIN,iBRDF_land,iBPDF_land,iBRM_water )
!      stop 'stop in forw_SOS_RT: development needed'
      ns=1
      if (abs_run) then
         ns = abs_data_forw_im%ABS_GS_WL(curr_wl_index)%NEXP
         SLout_tmp  = 0.0
         SQout_tmp  = 0.0
         SUout_tmp  = 0.0
         UFG_tmp = 0.0
         DFG_tmp = 0.0
         UFX_tmp = 0.0
         DFX_tmp = 0.0
         UFX0_tmp= 0.0
         DFX0_tmp= 0.0
      end if ! abs_run
      DO i1=1,ns(1)
      DO i2=1,ns(2)
      DO i3=1,ns(3)
         !IF((RIN%IPRI_additional_info) .and. (IW .eq. 1))  THEN ! "write_rad_transf_input" does not print OS_H input for GAS
         IF (RIN%IPRI_additional_info)  THEN ! "write_rad_transf_input" does not print OS_H input for GAS
            CALL write_rad_transf_input ( RIN%DLSF%external_file_path, RIN%IPRI_additional_info,    &
                                          OSHP%IMSC, OSHP%NG, OSHP%NN, OSHP%NF, RIN%eps_err,        &
                                          iBRDF_land, iBPDF_land, iBRM_water, pixel_fit%land_percent, &
                                          NQDR, ITRONC, iop, NT1,                                   &
                                          sza, NBV_comb, vis_comb(:), fiv_comb(:),                  &
                                          IW, WAVE,                                                 &
                                          surf_land_par_num, surf_land_par_vect,                    &
                                          surf_water_par_num, surf_water_par_vect,                  &
                                          EXT_os, SSA_os,                                           &
                                          NANG, ANGL(1:NANG),                                       &
                                          HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,                    &
                                          laerosol, lsurface ,           	                          &
                                          GOUT_aerosol_phmx_pixel_wl%ph11(1:NANG,:),                &
                                          GOUT_aerosol_phmx_pixel_wl%ph12(1:NANG,:),                &
                                          GOUT_aerosol_phmx_pixel_wl%ph22(1:NANG,:),                &
                                          GOUT_aerosol_phmx_pixel_wl%ph33(1:NANG,:)                 &
                                        )
         ENDIF ! RIN%IPRI_additional_info

         if (abs_run) then
            nhinp(1:natm) = DISCRVD%nh
            do isd=1,natm
              hinp_km(:,isd) = DISCRVD%h_km(:)
              vdinp(:,isd) = DISCRVD%val(:,isd)
            enddo
            ! gas segment profile
            call ABS_GAS_PROFILE (  abs_data_forw_im,                  & ! IN
                                    curr_wl_index,                     &
                                    i1, i2, i3,                        &
                                    nhinp_max,                         &
                                    nhinp(natm+ngas_kd),               &  ! OUT
                                    hinp_km(1:nhinp_max,natm+ngas_kd), &
                                    vdinp(1:nhinp_max,natm+ngas_kd),   &
                                    EXT_os(natm+ngas_kd)               &
                                 )
!XH         do not consider atmospheric profile for zero optical thickness
            if (EXT_os(natm+ngas_kd) .eq. 0.0) then
              ngas_kd = 0
            endif
            DISCRVD%natm = natm + ngas_kd
            DISCRVD%ngas = ngas_kd

            ! Recompute all discret profiles of atmospheric components at altitudes for fluxes
            ! above nhinp_max = KVERT_WD
            call discret_vertical_distr ( HGR_km, HMAX_atm_km, distr_type, &
                                          0.0, 0.0, nhinp_max, nhinp(:),   &
                                          hinp_km(1:nhinp_max,:),          &
                                          vdinp(1:nhinp_max,:),            &
                                          DISCRVD )
         end if ! abs_run

         call CPU_time(time_start)

!        Radiative Transfer
!        NBV           - number of observation angles
!        vis(1:NBV,IP) - zenith  observation angles
!        fiv(1:NBV,IP) - azimuth observation angles
!        sza           - Solar Zenith Angle tetas
!         write(*,*) 'before OS_HERMAN: HOBS,MASL  ',pixel_fit%HOBS,pixel_fit%MASL

!write(*,*) 'before radiative_transfer_SOS keyEl=',RIN%DLSF%keyEl
         if (abs_run) then
!XH         calculate broadband flux without aerosol
            EXT0_os = EXT_os
            EXT0_os(1:naer) = 0.0
            CALL radiative_transfer_SOS (  RIN,                                                          & ! IN
                                           IW,OSHP%IMSC,OSHP%NG,OSHP%NN,OSHP%NF,                         &
                                           iBRDF_land,iBPDF_land,iBRM_water,pixel_fit%land_percent,      &
                                           NQDR,ITRONC,iop, NT1,                                         &
                                           sza, NBV_comb, vis_comb(:), fiv_comb(:),                      &
                                           WAVE,                                                         &
                                           surf_land_par_num(1),surf_land_par_num(2),surf_land_par_vect, &
                                           surf_water_par_num,surf_water_par_vect,                       &
                                           !NM, ngas,                                                     &
                                           EXT0_os, SSA_os,                                              &
                                           KMpar, NANG, ANGL,                                            &
                                           HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,                        &
                                           !KVERTM,NH0,H01,PROF0,                                         &
                                           laerosol,lsurface,                                            &
                                           !RIN%mol_prof_type, RIN%aer_prof_type,                         & ! AL
                                           GOUT_aerosol_phmx_pixel_wl%ph11(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph12(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph22(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph33(:,:),                         &
                                           ASCAT,                                                        &
                                           SLout_comb,SQout_comb,                                        & ! OUT
                                           SUout_comb,SLPout_comb,                                       &
                                           UFG,DFG,                                                      &
                                           GOUT_bbflux_pixel%NHLV,                                       &
                                           GOUT_bbflux_pixel%HLV,                                        &
                                           GOUT_bbflux_pixel%BBUFX0,                                     &
                                           GOUT_bbflux_pixel%BBDFX0                                      &
                                        )
            if ( stop_report%status ) return
            if (RIN%DLSF%keyEl .eq. 1) then
               SQout_comb=0.0
               SUout_comb=0.0
               SLPout_comb=0.0
            end if
!XH         calculate broadband flux with aerosol
!XH         interface to call scalar or vector RTF subroutine
            CALL radiative_transfer_SOS (  RIN,                                                          & ! IN
                                           IW,OSHP%IMSC,OSHP%NG,OSHP%NN,OSHP%NF,                         &
                                           iBRDF_land,iBPDF_land,iBRM_water,pixel_fit%land_percent,      &
                                           NQDR, ITRONC, iop, NT1,                                       &
                                           sza, NBV_comb, vis_comb(:), fiv_comb(:),                      &
                                           WAVE,                                                         &
                                           surf_land_par_num(1),surf_land_par_num(2),surf_land_par_vect, &
                                           surf_water_par_num,surf_water_par_vect,                       &
                                           !NM, ngas,                                                     &
                                           EXT_os, SSA_os,                                               &
                                           KMpar, NANG, ANGL,                                            &
                                           HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,                        &
                                           !KVERTM,NH0,H01,PROF0,                                         &
                                           laerosol,lsurface,                                            &
                                           !RIN%mol_prof_type, RIN%aer_prof_type,                         & ! AL
                                           GOUT_aerosol_phmx_pixel_wl%ph11(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph12(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph22(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph33(:,:),                         &
                                           ASCAT,                                                        &
                                           SLout_comb,SQout_comb,                                        & ! OUT
                                           SUout_comb,SLPout_comb,                                       &
                                           UFG,DFG,                                                      &
                                           GOUT_bbflux_pixel%NHLV,                                       &
                                           GOUT_bbflux_pixel%HLV,                                        &
                                           GOUT_bbflux_pixel%BBUFXA,                                     &
                                           GOUT_bbflux_pixel%BBDFXA                                      &
                                        )
            if ( stop_report%status ) return
            if (RIN%DLSF%keyEl .eq. 1) then
               SQout_comb=0.0
               SUout_comb=0.0
               SLPout_comb=0.0
            end if
         else
!XH         without flux calculation
!XH         interface to call scalar or vector RTF subroutine
            CALL radiative_transfer_SOS (  RIN,                                                          & ! IN
                                           IW,OSHP%IMSC,OSHP%NG,OSHP%NN,OSHP%NF,                         &
                                           iBRDF_land,iBPDF_land,iBRM_water,pixel_fit%land_percent,      &
                                           NQDR, ITRONC, iop, NT1,                                       &
                                           sza, NBV_comb, vis_comb(:), fiv_comb(:),                      &
                                           WAVE,                                                         &
                                           surf_land_par_num(1),surf_land_par_num(2),surf_land_par_vect, &
                                           surf_water_par_num,surf_water_par_vect,                       &
                                           !NM, ngas,                                                     &
                                           EXT_os, SSA_os,                                               &
                                           KMpar, NANG, ANGL,                                            &
                                           HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,                        &
                                           !KVERTM,NH0,H01,PROF0,                                         &
                                           laerosol,lsurface,                                            &
                                           !RIN%mol_prof_type, RIN%aer_prof_type,                         & ! AL
                                           GOUT_aerosol_phmx_pixel_wl%ph11(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph12(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph22(:,:),                         &
                                           GOUT_aerosol_phmx_pixel_wl%ph33(:,:),                         &
                                           ASCAT,                                                        &
                                           SLout_comb,SQout_comb,                                        & ! OUT
                                           SUout_comb,SLPout_comb,                                       &
                                           UFG,DFG                                                       &
                                        )
            if ( stop_report%status ) return
            if (RIN%DLSF%keyEl .eq. 1) then
               SQout_comb=0.0
               SUout_comb=0.0
               SLPout_comb=0.0
            end if
         end if !abs_run

!XH      integration over the absorbing band to get the radiance with K-distribution for vector case
         if (abs_run) then
!XH         set number of gas components back
            ngas_kd = 1
            call ABS_GAS_INT(                                           &
                               curr_wl_index,i1,i2,i3,                  & ! IN
                               abs_data_forw_im,                        &
                               NBV=NBV_comb,                            &
                               SL=SLout_comb,                           &
                               SQ=SQout_comb,                           &
                               SU=SUout_comb,                           &
                               SL_TMP=SLout_tmp,                        &
                               SQ_TMP=SQout_tmp,                        &
                               SU_TMP=SUout_tmp,                        &
                               UFG=UFG,                                 &
                               DFG=DFG,                                 &
                               UFG_TMP=UFG_tmp,                         &
                               DFG_TMP=DFG_tmp,                         &
                               NLV=GOUT_bbflux_pixel%NHLV,              &
                               UFX0=GOUT_bbflux_pixel%BBUFX0,           &
                               DFX0=GOUT_bbflux_pixel%BBDFX0,           &
                               UFX =GOUT_bbflux_pixel%BBUFXA,           &
                               DFX =GOUT_bbflux_pixel%BBDFXA,           &

                               UFX0_TMP=UFX0_tmp,                       &
                               DFX0_TMP=DFX0_tmp,                       &
                               UFX_TMP =UFX_tmp,                        &
                               DFX_TMP =DFX_tmp                         & ! OUT
                            )
         end if ! abs_run

      END DO ! i3
      END DO ! i2
      END DO ! i1

      if (abs_run) then
         SLout_comb(1:NBV_comb)  = SLout_tmp(1:NBV_comb)
         SQout_comb(1:NBV_comb)  = SQout_tmp(1:NBV_comb)
         SUout_comb(1:NBV_comb)  = SUout_tmp(1:NBV_comb)
         SLPout_comb(1:NBV_comb) = SQRT( SQout_comb(1:NBV_comb)**2      &
                                        +SUout_comb(1:NBV_comb)**2)
         UFG = UFG_tmp
         DFG = DFG_tmp
         GOUT_bbflux_pixel%BBUFX0=UFX0_tmp
         GOUT_bbflux_pixel%BBDFX0=DFX0_tmp
         GOUT_bbflux_pixel%BBUFXA=UFX_tmp
         GOUT_bbflux_pixel%BBDFXA=DFX_tmp
      end if ! abs_run

      call CPU_time(time_finish)
      xtime=xtime+time_finish-time_start

      GOUT_surface_pixel_wl%salbedo = 0.0
      IF (OSHP%IMSC .NE. 1) THEN
         GOUT_surface_pixel_wl%salbedo = UFG/DFG
!         WRITE(*,*) WAVE,UFG/DFG,'   SURF_ALBEDO'
      END IF ! OSHP%IMSC .NE. 1

      lexist=.false.
      if (RIN%IPRI_verbose) then
         DO IV=1,NBV_comb
            status_funct = check_nan(1,SLout_comb(IV:IV))
            if (.not. status_funct) then
               write(*,*) 'in forw_SOS_RT :'
               write(*,'(i5,e13.5,a)') IV,SLout_comb(IV),'  - IV,SLout_comb(IV)'
            end if
            IF (SLout_comb(IV) .LE. 0. .OR. (.not. status_funct)) THEN
               !SLout_comb(IV) = 1e+1
               lexist = .true.
            END IF
            status_funct = check_nan(1,SQout_comb(IV:IV))
            if (.not. status_funct) then
               write(*,*) 'in forw_SOS_RT :'
               write(*,'(i5,e13.5,a)') IV,SQout_comb(IV),'  - IV,SQout_comb(IV)'
            end if
            IF (.not. status_funct) THEN
               !SQout_comb(IV) = 1e-5
               lexist = .true.
            END IF
            status_funct = check_nan(1,SUout_comb(IV:IV))
            if (.not. status_funct) then
               write(*,*) 'in forw_SOS_RT :'
               write(*,'(i5,e13.5,a)') IV,SUout_comb(IV),'  - IV,SUout_comb(IV)'
            end if
            IF (.not. status_funct) THEN
               !SUout_comb(IV) = 1e-5
               lexist = .true.
            END IF
         END DO ! IV
      end if

      IF(RIN%IPRI_additional_info .and. lexist) THEN
        write(*,'(3a8,4(a11),a11)')'thd','vis','phi','plh','pqh','puh','pph',' wl'
        DO IV=1,NBV_comb
          write (*,'(4f9.3,4f11.6,f12.6)')       &
                ASCAT(iv),vis_comb(iv),fiv_comb(iv),sza,      &
                SLout_comb(IV),SQout_comb(IV),SUout_comb(IV),SLPout_comb(IV)/SLout_comb(IV),WAVE
        END DO ! IV
        !call print_crash_report()
        !stop 'stop: NaN in forw_SOS_RT'
      END IF ! .not. RIN%IPRI_additional_info .and
!
      RETURN
      END SUBROUTINE forw_SOS_RT
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      end module mod_rt
