! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains :
   
! subroutine forward_model_pixel_wl
! subroutine forward_model_pixel
! subroutine forw_phase_matrix
! subroutine forw_lidar_signal
! subroutine forw_radiative_transfer
! subroutine forward_model_pixel_PHMX

      module mod_forward_model

      use mod_stop_report

      contains
!
      subroutine forward_model_pixel_wl (  RIN,ipix,                &  ! IN
                                           IW,WAVE,ind_wl,lresult,  &
                                           NBIN,RADIUS,SD,          &
                                           SHD,NSHAPE1,             &
                                           RREAL,RIMAG,             &
                                           OSHP,                    &
                                           BRF_land,BRP_land,BRM_water,tau_mol,  &
                                           HOBS_km,HGR_km,HMAX_atm_km, &
                                           NHVP_fit,HVP_fit_km,     &
                                           NHVP_retr,HVP_retr_km,   &
                                           H0,sigma_aerosol,CL,     &
                                           abs_data_forw_im,        &
                                           pixel_fit,               & ! INOUT
                                           GOUT_aerosol,            &
                                           GOUT_clouds,             &
                                           GOUT_surface,            &
                                           GOUT_bbflux_pixel,       &
                                           NANG,ANGL,               &  ! OUT
                                           KERNELS1,KERNELS2        &
                                        )

      use mod_par_DLS,   only : KMpar
      use mod_par_inv,   only : KPARS,KSHAPE,KBF,KIDIM3,KVERTM,KIP,KNBVM
      use mod_par_OS,    only : NMG,KSD,NBVM,NMM,KNT,NG0T
!XH   module related to gas absorption
      use mod_abs_kd

      use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_alloc_kernels
	  	  
      IMPLICIT NONE
!	------------------------------------------------------------------------------------------------------
! IN:
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  ipix,IW,ind_wl 
      real,                       intent(in)  ::  HOBS_km,HGR_km,HMAX_atm_km                            
      integer,                    intent(in)  ::  NHVP_fit,NHVP_retr
      real,dimension(KVERTM),     intent(in)  ::  HVP_fit_km,HVP_retr_km 
      real,                       intent(in)  ::  WAVE
      integer,dimension(KSD),     intent(in)  ::  NSHAPE1  
      real,dimension(KBF),        intent(in)  ::  BRF_land,BRP_land,BRM_water
      type(OSH_par),              intent(in)  ::  OSHP		
      real,                       intent(in)  ::  tau_mol
      real,                       intent(in)  ::  CL
      real,dimension(KVERTM,KSD), intent(inout) ::  H0 ! AL contains parameters of aerosol vertical 
                                                       ! distribution or distribution itself
      real,dimension(KSD),        intent(in)  ::  RREAL, RIMAG, sigma_aerosol
!XH   data of gas absorbing band
      type (DATA_ABS),            intent(in)  ::  abs_data_forw_im

      real,dimension(KIDIM3,KSD), intent(in)  ::  RADIUS,SD	  
      real,dimension(KSHAPE,KSD), intent(in)  ::  SHD  
      integer,dimension(KSD),     intent(in)  ::  NBIN
      logical,                    intent(in)  ::  lresult
!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(output_segment_particles),  intent(inout)  ::  GOUT_aerosol
      type(output_segment_particles),  intent(inout)  ::  GOUT_clouds
      type(output_segment_surface),    intent(inout)  ::  GOUT_surface
      type(output_pixel_bbflux),       intent(inout), optional  ::  GOUT_bbflux_pixel
      type(pixel),                 intent(inout)  ::  pixel_fit
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2
!	------------------------------------------------------------------------------------------------------
! OUT:
      integer,                        intent(out)           ::  NANG
      real,dimension(KMpar),          intent(out)           ::  ANGL
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer                       ::  ISD,IP,IAN
      real, dimension(KNBVM)        ::  meas
      logical                       ::  laerosol, lsurface
      integer                       ::  nmeas_type, meas_type
      real                          ::  sca
      integer                       ::  NBV_comb
      real, dimension(2*NBVM)       ::  SQout_comb,SUout_comb,SLPout_comb,SLout_comb
      logical                       ::  SvR_first
!XH   switch for broadband flux calculation
      logical                       ::  iFlux
!	------------------------------------------------------------------------------------------------------
! SAVE:
      real,                       save  ::  WAVE_save = 0.0
      real,dimension(KIDIM3,KSD), save  ::  SD_save  = 0.0
      real,dimension(KSHAPE,KSD), save  ::  SHD_save = 0.0
      real,dimension(KVERTM,KSD), save  ::  H0_save = 0.0
      real,dimension(KSD),        save  ::  RREAL_save = 0.0 	
      real,dimension(KSD),        save  ::  RIMAG_save = 0.0
      real,dimension(KBF),        save  ::  BRF_land_save  = 0.0
      real,dimension(KBF),        save  ::  BRP_land_save  = 0.0
      real,dimension(KBF),        save  ::  BRM_water_save = 0.0
      real,dimension(KSD),        save  ::  sigma_aerosol_save = 0.0
      integer,                    save  ::  IMSC_save = -1
!	------------------------------------------------------------------------------------------------------
!     IW     - index of wave length in array of wave length for current pixel
!     WAVE   - value of wave length for index IW from array of wave length for current pixel
!     ind_wl - index of wave length in general array of wave length for inversion
!     write(*,*) 'in forward_model_pixel_wl: iw =',iw,WAVE,ind_wl
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
! lresult=.true. only single scattering properties are calculated
      if(.not. lresult) then
        laerosol = .false.
        lsurface = .false.
        if (WAVE .ne. WAVE_save) then
          laerosol = .true.
          lsurface = .true.
        else
          if (OSHP%IMSC .ne. IMSC_save) then
            laerosol = .true.
            lsurface = .true.
          end if
! aerosol
          if (.not. laerosol) then
            if(ANY(SD(:,:)  .ne. SD_save(:,:)))   &
            laerosol = .true.
          endif
          if (.not. laerosol) then
            if(ANY(SHD(:,:) .ne. SHD_save(:,:)))  &
            laerosol = .true.
          endif
          if (.not. laerosol) then
            if(ANY(H0(:,:)  .ne. H0_save(:,:)))   &
            laerosol = .true.
          endif
          if (.not. laerosol) then
            if(ANY(RREAL(:) .ne. RREAL_save(:)))  &
            laerosol = .true.
          endif
          if (.not. laerosol) then
            if(ANY(RIMAG(:) .ne. RIMAG_save(:)))  &
            laerosol = .true.
          endif
          if (.not. laerosol) then
            if(ANY(sigma_aerosol(:) .ne. sigma_aerosol_save(:)))  &
            laerosol = .true.
          endif
! surface
          if (.not. lsurface) then
            if(ANY(BRF_land(:)  .ne. BRF_land_save(:)))  &
            lsurface = .true.
          endif
          if (.not. lsurface) then
            if(ANY(BRP_land(:)  .ne. BRP_land_save(:)))  &
            lsurface = .true.
          endif
          if (.not. lsurface) then
            if(ANY(BRM_water(:) .ne. BRM_water_save(:))) &
            lsurface = .true.
          endif
        endif ! WAVE .eq. WAVE_save
        if (laerosol) then
          SD_save(:,:)  = SD(:,:)
          SHD_save(:,:) = SHD(:,:)
          H0_save(:,:)  = H0(:,:)
          RREAL_save(:) = RREAL(:)
          RIMAG_save(:) = RIMAG(:)
          sigma_aerosol_save(:) = sigma_aerosol(:)
        endif
        if (lsurface) then
          BRF_land_save(:)  = BRF_land(:)
          BRP_land_save(:)  = BRP_land(:)
          BRM_water_save(:) = BRM_water(:)
        endif
        !laerosol=.true.
        !lsurface=.true.
!       write(*,*) 'laerosol=',laerosol,'  lsurface=',lsurface
        WAVE_save = WAVE
        IMSC_save = OSHP%IMSC
      else
        !only single scattering properties are calculated
        laerosol = .true.
      endif ! lresult
!	------------------------------------------------------------------------------------------------------
!  Aerosol phase matrix, extinction and single scattering albedo
      if(RIN%NSD .gt. 0 .and. laerosol) then
        call forw_single_scattering_particle_properties ( RIN, RIN%NSD, NBIN, RADIUS, SD, & ! IN
                                                          SHD, NSHAPE1,               &
                                                          ind_wl, WAVE, RREAL, RIMAG, &
                                                          NANG, ANGL, ipix,           & ! OUT
                                                          GOUT_aerosol,               &
                                                          KERNELS1, KERNELS2          &
                                                        )
        if ( stop_report%status ) return
        if(RIN%IPRI_additional_info .and. OSHP%IMSC .eq. -2) then
          do ISD=1,RIN%NSD
          write(*,*) RREAL(ISD),RIMAG(ISD),ISD,IW,WAVE,ind_wl,'  RREAL RIMAG ISD IW WAVE ind_wl in forw_IMAGE_I_IW'		 
          write(*,*) ISD,GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl)%ssa(ISD),  & 
                         GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl)%ext(ISD),'  ISD SSA EXT - aerosol'
          !do IAN=1,NANG
          !  write(*,'(i4,5e14.5)') ISD,ANGL(IAN),  &
          !  GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl)%ph11(IAN,ISD),GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl)%ph12(IAN,ISD),  &
          !  GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl)%ph22(IAN,ISD),GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl)%ph33(IAN,ISD)
          !enddo ! IAN
          enddo ! ISD
        endif ! IPRI
      endif ! RIN%NSD .gt. 0 .and. laerosol
!	------------------------------------------------------------------------------------------------------
      if (present(GOUT_bbflux_pixel)) then
!XH      for broadband flux calculation
         iFlux = .true.
         call forw_radiative_transfer (   iFlux,                                    & ! IN
                                          RIN,                                      &
                                          IW,WAVE,ind_wl,IP,                        &
                                          OSHP,                                     &
                                          BRF_land,BRP_land,BRM_water,tau_mol,      &
                                          HOBS_km,HGR_km,HMAX_atm_km,               &
                                          NHVP_retr,HVP_retr_km,H0,sigma_aerosol,   &
                                          laerosol,lsurface,                        &
                                          NANG,ANGL,                                &
                                          abs_data_forw_im,                         &
                                          pixel_fit,                                &

                                          GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl),  & ! INOUT
                                          GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl), &
                                          GOUT_clouds%opt%pixel(ipix)%wl(ind_wl),   &
                                          GOUT_clouds%phmx%pixel(ipix)%wl(ind_wl),  &
                                          GOUT_surface%pixel(ipix)%wl(ind_wl),      &
                                          GOUT_bbflux_pixel,                        &

                                          NBV_comb,                                 & ! OUT
                                          SLout_comb,SQout_comb,                    &
                                          SUout_comb,SLPout_comb                    &
                                      )
         return
      end if
!	------------------------------------------------------------------------------------------------------
      if (.not. lresult) then
         nmeas_type = pixel_fit%meas(IW)%NIP
         !write(*,*) 'nmeas_type =',nmeas_type
         SvR_first = .true.
         LOOP_meas_type: do IP=1,nmeas_type
            meas_type = pixel_fit%meas(IW)%meas_type(IP)
            !write(*,*) 'ip =',ip,'  meas_type =',meas_type
!	------------------------------------------------------------------------------------------------------
!           Single scattering optical properties
            if (meas_type .gt. meas_type_tau_beg .and. meas_type .lt. meas_type_phm_end)  then
               call set_pixel_phase_matr_fit (                                            &
                                                RIN%NSD,tau_mol,                          &
                                                NANG,ANGL,                                &
                                                GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl),  &
                                                GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl), &
                                                meas_type,IW,IP,RIN%iPOBS,                &
                                                pixel_fit                                 &
                                             )
              if ( stop_report%status ) return
!	------------------------------------------------------------------------------------------------------
!           Lidar signal and if developped depolarization
!           This part has to be developed for depolarization and Roman lidar
            else if (meas_type .gt. meas_type_lid_beg .and. meas_type .lt. meas_type_lid_end) then
               call forw_lidar_signal (                                                   &
                                         HOBS_km,HGR_km,HMAX_atm_km,                      & ! IN
                                         NHVP_fit,HVP_fit_km,                             &
                                         NHVP_retr,HVP_retr_km,                           &
                                         H0,CL,                                           &
                                         RIN%NSD,                                         &
                                         GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl)%ext(:),  &
                                         GOUT_aerosol%lidar%pixel(ipix)%wl(ind_wl)%lr(:), &
                                         tau_mol,                                         &
                                         pixel_fit%meas(IW)%IFMP(IP),                     &
                                         pixel_fit%meas(IW)%MPROF(:,IP),                  &
                                         RIN%mol_prof_type,                               & ! AL
                                         RIN%aer_prof_type,                               & ! AL
                                         meas_type,meas(:)                                & ! INOUT
                                      )
               if ( stop_report%status ) return
               call set_pixel_lidar_signal_fit ( NHVP_fit,meas(1:NHVP_fit), &
                                                 meas_type,IW,IP,  &
                                                 pixel_fit         &
                                               )
!	------------------------------------------------------------------------------------------------------
!           Radiative transfer accounting for multiple scattering
            else if(meas_type .gt. meas_type_SvR_beg .and. meas_type .lt. meas_type_SvR_end) then
               if (SvR_first) then
                  ! Radiative transfer routine is called only once for first Stokes vector measurement
                  ! type because its output is a Stokes vector (I,Q,U) in one routine call.
                  iFlux = .false.
                  call forw_radiative_transfer (   iFlux,                                    & ! IN
                                                   RIN,                                      &
                                                   IW,WAVE,ind_wl,IP,                        &
                                                   OSHP,                                     &
                                                   BRF_land,BRP_land,BRM_water,tau_mol,      &
                                                   HOBS_km,HGR_km,HMAX_atm_km,               &
                                                   NHVP_retr,HVP_retr_km,H0,sigma_aerosol,   &
                                                   laerosol,lsurface,                        &
                                                   NANG,ANGL,                                &
                                                   abs_data_forw_im,                         &
                                                   pixel_fit,                                &

                                                   GOUT_aerosol%opt%pixel(ipix)%wl(ind_wl),  & ! INOUT
                                                   GOUT_aerosol%phmx%pixel(ipix)%wl(ind_wl), &
                                                   GOUT_clouds%opt%pixel(ipix)%wl(ind_wl),   &
                                                   GOUT_clouds%phmx%pixel(ipix)%wl(ind_wl),  &
                                                   GOUT_surface%pixel(ipix)%wl(ind_wl),      &
                                                   GOUT_bbflux_pixel,                        &

                                                   NBV_comb,                                 & ! OUT
                                                   SLout_comb,SQout_comb,                    &
                                                   SUout_comb,SLPout_comb                    &
                                               )
                  if ( stop_report%status ) return
                  call set_pixel_Stokes_vec_fit (  IW,IP,nmeas_type,RIN%iPOBS, &
                                                   NBV_comb,                   &
                                                   SLout_comb,SQout_comb,      &
                                                   SUout_comb,SLPout_comb,     &
                                                   pixel_fit                   &
                                                )
                  SvR_first = .false.
               end if
            end if ! meas_type .gt. meas_type_tau_beg .and.
!	------------------------------------------------------------------------------------------------------
         end do LOOP_meas_type
      endif ! lresult
!
      RETURN
      END SUBROUTINE forward_model_pixel_wl

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! forw_phase_matrix
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine forw_phase_matrix (                         &
                                    NSD,NBIN,RADIUS,SD,      & ! IN 
                                    SHD,NSHAPE1,             &
                                    ind_wl,WAVE,RREAL,RIMAG, &
                                    use_models,              &
                                    NANG,ANGL,ipix,          & ! OUT
                                    GOUT_particles_opt_pixel_wl,  &
                                    GOUT_particles_phmx_pixel_wl, &
                                    DLSF,KERNELS1,KERNELS2   & 
                                   )
                                   
      use mod_par_inv, only : KSHAPE,KIDIM3
      use mod_par_OS,  only : KSD
      !use mod_type_DLS
      use mod_alloc_kernels
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      !use mod_index_cloud
            
      implicit none
! -----------------------------------------------------------------------
! IN:
      real,dimension(KIDIM3,KSD), intent(in)  ::  RADIUS,SD	  
      real,dimension(KSHAPE,KSD), intent(in)  ::  SHD  
      integer,dimension(KSD),     intent(in)  ::  NBIN,NSHAPE1
      real,dimension(KSD),        intent(in)  ::  RREAL,RIMAG
      real,                       intent(in)  ::  WAVE
      integer,                    intent(in)  ::  NSD,ind_wl,ipix
      type(iP_flags_for_DLS),     intent(in)  ::  DLSF
      logical,                    intent(in)  ::  use_models
! -----------------------------------------------------------------------
! OUT:
      integer,                        intent(out) ::  NANG
      real,dimension(KMpar),          intent(out) ::  ANGL
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_particles_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_particles_phmx_pixel_wl

      !type(output_pixel_main_wl),     intent(inout)  ::  GOUT_main_pixel_wl
      !type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_phmx_pixel_wl
!	----------------------------------------------------------------------
! INOUT:
      type(kernels_triangle_bin), intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),intent(inout)  ::  KERNELS2
! -----------------------------------------------------------------------

#if defined(SPHEROID)
!      if(.not. present(KERNELS1) .or. .not. present(KERNELS2)) then
!         write(*,*) 'KERNELS1,KERNELS2 are not present in forw_phase_matrix'
!         stop 'stop in forw_phase_matrix'
!      endif
      call spheroid_package (                          & ! IN
                              NSD,NBIN,RADIUS,SD,      & 
                              SHD,NSHAPE1,             &
                              ind_wl,WAVE,RREAL,RIMAG, &
                              use_models,              &
                              NANG,ANGL,               & ! OUT
                              GOUT_particles_opt_pixel_wl,  &
                              GOUT_particles_phmx_pixel_wl, &
                              DLSF,KERNELS1,KERNELS2   &
                            )

!#elif defined(ANOTHER)
!      call another_model()
#else
#error No PHASE_MATRIX configured!
#endif

      return
      end subroutine forw_phase_matrix
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! forw_lidar_signal
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine forw_lidar_signal (                                &
                                    HOBS_km,HGR_km,HMAX_atm_km,     & ! IN
                                    NHVP_fit,HVP_fit_km,            &
                                    NHVP_retr,HVP_retr_km,          &
                                    H0,CL,                          &
                                    NSD,EXTA,LRA,tau_mol,           & 																																			
                                    IFMP,MPROF,                     &
                                    mol_prof_type, aer_prof_type,   & !AL
                                    meas_type,meas                  & ! INOUT									
                                   )

      use mod_par_inv,   only : KVERTM
      use mod_par_OS,    only : KSD  
	  	  
      implicit none
!	------------------------------------------------------------------------------------------------------
! IN:
      integer,                    intent(in)  ::  NSD,IFMP
      real,                       intent(in)  ::  HOBS_km,HGR_km,HMAX_atm_km                            
      integer,                    intent(in)  ::  NHVP_fit,NHVP_retr
      real,dimension(KVERTM),     intent(in)  ::  HVP_fit_km,HVP_retr_km 
      real,                       intent(in)  ::  CL,tau_mol
      real,dimension(KSD),        intent(in)  ::  EXTA,LRA
      real,dimension(KVERTM,KSD), intent(in)  ::  H0 ! contains parameters of aerosol vertical 
                                                     ! distribution or distribution itself
      real,dimension(KVERTM),     intent(in)  ::  MPROF
      integer,                    intent(in)  ::  meas_type
      integer,                    intent(in)  ::  mol_prof_type, aer_prof_type
!	------------------------------------------------------------------------------------------------------
! INOUT:
      real,dimension(KVERTM),     intent(out) ::  meas
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------

#if defined(GARRLIC)
      call lidar_garrlic (                                  &
                            HOBS_km,HGR_km,HMAX_atm_km,     & ! IN
                            NHVP_fit,HVP_fit_km,            &
                            NHVP_retr,HVP_retr_km,          &
                            H0,CL,                          &
                            NSD,EXTA,LRA,tau_mol,           & 																																			
                            IFMP,MPROF,                     &
                            mol_prof_type, aer_prof_type,   &
                            meas_type,meas                  & ! INOUT
                         )

!#elif defined(ANOTHER)
!      call another_model()
#else
!#error No LIDAR configured!
#warning No LIDAR configured!
#endif

      return
      end subroutine forw_lidar_signal
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! forw_radiative_transfer
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine forw_radiative_transfer (                                          &
                                             iFlux,                                 &
                                             RIN,                                   & ! IN
                                             IW,WAVE,ind_wl,IP,                     &
                                             OSHP,                                  &
                                             BRF_land,BRP_land,BRM_water,tau_mol,   &
                                             HOBS_km,HGR_km,HMAX_atm_km,            &
                                             NHVP_retr,HVP_retr_km,H0,sigma_aerosol,&
                                             laerosol,lsurface,                     &
                                             NANG,ANGL,                             &
                                             abs_data_forw_im,                      &
                                             pixel_fit,                             &
                                 
                                             GOUT_aerosol_opt_pixel_wl,             & ! INOUT
                                             GOUT_aerosol_phmx_pixel_wl,            &
                                             GOUT_clouds_opt_pixel_wl,              &
                                             GOUT_clouds_phmx_pixel_wl,             &
                                             GOUT_surface_pixel_wl,                 &
                                             GOUT_bbflux_pixel,                     &

                                             NBV_comb,                              & ! OUT
                                             SLout_comb,SQout_comb,                 &
                                             SUout_comb,SLPout_comb                 &
                                         )
 
      use mod_os
      use mod_rt
      use mod_par_DLS,   only : KMpar
      use mod_par_inv,   only : KPARS,KSHAPE,KBF,KIDIM3,KVERTM,KIP,KNBVM
      use mod_par_OS,    only : NMG,KSD,NBVM,NMM,KNT,NG0T
!XH   module related to gas absorption
      use mod_abs_kd

!      use mod_index_cloud
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_sdata_meas_type
      use mod_sdata_derived_type
!      use mod_atm_lid
      use mod_alloc_kernels
      use mod_vertical_distr_derived_type

      IMPLICIT NONE
!	------------------------------------------------------------------------------------------------------
! IN:
!XH   switch for broadband flux calculation
      logical,                    intent(in)  ::  iFlux
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  IW,ind_wl,IP 
      real,                       intent(in)  ::  HOBS_km,HGR_km,HMAX_atm_km                            
      integer,                    intent(in)  ::  NHVP_retr
      real,dimension(KVERTM),     intent(in)  ::  HVP_retr_km 
      real,                       intent(in)  ::  WAVE
      real,dimension(KBF),        intent(in)  ::  BRF_land,BRP_land,BRM_water
      type(OSH_par),              intent(in)  ::  OSHP		
      real,                       intent(in)  ::  tau_mol
      real,dimension(KVERTM,KSD), intent(inout)  ::  H0 ! AL contains parameters of aerosol vertical
                                                        ! distribution or distribution itself
      real,dimension(KSD),        intent(in)     ::  sigma_aerosol
      logical,                    intent(in) ::  laerosol, lsurface
!XH   data of o2 absorption, do not delete
      type (DATA_ABS),            intent(in)  ::  abs_data_forw_im

      integer,                    intent(in)    ::  NANG
      real,dimension(KMpar),      intent(in)    ::  ANGL
      integer,                    intent(inout) ::  NBV_comb
      real, dimension(2*NBVM),    intent(inout) ::  SQout_comb,SUout_comb, &
                                                    SLPout_comb,SLout_comb
!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_aerosol_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_aerosol_phmx_pixel_wl
      type(output_pixel_opt_wl),      intent(inout)  ::  GOUT_clouds_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(inout)  ::  GOUT_clouds_phmx_pixel_wl
      type(output_pixel_surface_wl),  intent(inout)  ::  GOUT_surface_pixel_wl
      type(pixel),                    intent(inout)  ::  pixel_fit
      type(output_pixel_bbflux),      intent(inout)  ::  GOUT_bbflux_pixel
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      integer,dimension(2)          ::  surf_land_par_num
      real, dimension(2*KBF)        ::  surf_land_par_vect
      integer                       ::  surf_water_par_num
      real, dimension(KBF)          ::  surf_water_par_vect
      integer                       ::  IDIM1, par_type, ns3
      type(discret_vertical_distribution) :: DISCRVD
!	------------------------------------------------------------------------------------------------------

#if defined(OSH)

      ! Vector of surface parametrs
      ! surf_par_num(1) - total number of surface parameters (BRDF+BPDF)
      ! surf_par_num(2) - number of SURF1 (BRDF) surface parameters
! ask PL
!tl Temporary solution if water parameters are not provided.
!tl if surf_water_par_num = 0 or surf_land_par_num = 0, there is a memory problem 
!tl in call developpe_ocean_land: x_vect_water(1:n_par_water)
!tl                               x_vect_land(1:n_par_land)
!tl      surf_land_par_num(:)   = 0
      surf_land_par_num(:)   = 1
      surf_land_par_vect(:)  = 0.0 
!tl      surf_water_par_num     = 0
      surf_water_par_num     = 1

      surf_water_par_vect(:) = 0.0 

      do IDIM1=1,RIN%NDIM%n1
         par_type = RIN%NDIM%par_type(IDIM1)
         if (par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF1_land_end) then
            surf_land_par_num(1:2) = RIN%NDIM%n2(IDIM1)
            surf_land_par_vect(1:RIN%NDIM%n2(IDIM1)) = BRF_land(1:RIN%NDIM%n2(IDIM1))
          if(par_type .eq. par_type_SURF1_land_RPV_BRDF) then
            ns3 = 3
            surf_land_par_vect(ns3) = surf_land_par_vect(ns3) - 1.0
          endif
         else if (par_type .gt. par_type_SURF2_land_beg .and. par_type .lt. par_type_SURF2_land_end) then
            surf_land_par_vect(surf_land_par_num(1)+1:surf_land_par_num(1)+RIN%NDIM%n2(IDIM1)) = BRP_land(1:RIN%NDIM%n2(IDIM1))
            surf_land_par_num(1) = surf_land_par_num(1) + RIN%NDIM%n2(IDIM1)
         else if (par_type .gt. par_type_SURF_water_beg .and. par_type .lt. par_type_SURF_water_end) then
            surf_water_par_num = RIN%NDIM%n2(IDIM1)
            surf_water_par_vect(1:RIN%NDIM%n2(IDIM1)) = BRM_water(1:RIN%NDIM%n2(IDIM1))
         end if ! par_type .gt.
      end do ! IDIM1
      !write(*,*) 'surf_model:    ',surf_model(1:2)
      !write(*,*) 'surf_par_vect: ',surf_par_vect(1:sum(surf_par_num(1:2)))
      call RT_discret_vertical_distr( RIN, HGR_km, HMAX_atm_km, &
                                      NHVP_retr, HVP_retr_km, H0, sigma_aerosol,  &
                                      pixel_fit%ifgas, pixel_fit%meas(IW)%gaspar, &
                                      DISCRVD )
      if ( stop_report%status ) return
      call forw_SOS_RT (iFlux,                                   & ! IN
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

!#elif defined(ANOTHER)
!      call another_model()
#else
!#error No RAD_TRANSF configured!
#warning No RAD_TRANSF configured!
#endif
!
      return
      end subroutine forw_radiative_transfer

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE forward_model_pixel (                      &
                                        RIN,OSHP,ipix,      & ! IN
                                        IWb,IWe,IWW,lresult,&
                                        tau_mol,            &
                                        NHVP_meas,HVP_meas, &
                                        nwl,wl,ind_wl,AP,   & ! INOUT
                                        pixel_fit,          &
                                        pixel_vec_fit,      &
                                        GOUT_aerosol,       &
                                        GOUT_clouds,        &
                                        GOUT_surface,       &

                                        NANG,ANGL,          & ! OUT
                                        KERNELS1,KERNELS2   &
                                     )
      use mod_os
      use mod_par_DLS, only : KMpar
      use mod_par_inv, only : KW,KMESS,KPARS,KSHAPE,KBF,  &
                              KIDIM2,KIDIM3,KVERTM
      use mod_par_OS, only : NMG,KSD,NBVM,NMM,KNT,NG0,HMAX_atm
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_sdata, only : set_pixel_meas_vector_fs
      use mod_alloc_kernels
!XH   modules related to gas absorption
      use mod_abs_kd

      use mod_forward_model_characteristics
!      use mod_index_clouds
!      use mod_ip_utility, only : check_nan

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      type(retr_input_settings),  intent(in)  ::  RIN
      type(OSH_par),              intent(in)  ::  OSHP
      integer,                    intent(in)  ::  IWW,IWb,IWe,ipix		                                            
      real,dimension(KW),         intent(in)  ::  tau_mol
      integer,                    intent(in)  ::  NHVP_meas ! number of heights for lidar measurements
      real,dimension(KVERTM),     intent(in)  ::  HVP_meas  ! heights for lidar measurements
      integer,                    intent(in)  ::  nwl
      !real,dimension(KW),         intent(in)  ::  WAVE
      !real,dimension(KWM),        intent(in)  ::  wl
      !integer,dimension(KWM),     intent(in)  ::  ind_wl
      real,dimension(KW),         intent(in)  ::  wl
      integer,dimension(KW),      intent(in)  ::  ind_wl
      real,dimension(KPARS),      intent(in)  ::  AP
!     lresult=.true. only single scattering properties are calculated for parameters retrieved at all wavelengths
      logical,                    intent(in)  ::  lresult
!	------------------------------------------------------------------------------------------------------
! OUT :
      type(pixel_vector),         intent(inout)::  pixel_vec_fit      
      integer,                    intent(out)  ::  NANG
      real,dimension(KMpar),      intent(out)  ::  ANGL
!	------------------------------------------------------------------------------------------------------
! INOUT :
      type(output_segment_particles),intent(inout)  :: GOUT_aerosol
      type(output_segment_particles),intent(inout)  :: GOUT_clouds
      type(output_segment_surface),  intent(inout)  :: GOUT_surface
!      type(output_segment_forcing),  intent(inout)  :: GOUT_forcing

      type(pixel),                intent(inout)  ::  pixel_fit
      type(kernels_triangle_bin), intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),intent(inout)  ::  KERNELS2
!	------------------------------------------------------------------------------------------------------
! LOCAL :
      real                       ::  HOBS_km      ! height of observations
      real                       ::  HGR_km       ! height above sea level
      real                       ::  HMAX_atm_km  ! top of atmosphere
      real,dimension(KVERTM)     ::  HVP_meas_km      
      real,dimension(KPARS)      ::	 APSING
!	------------------------------------------------------------------------------------------------------
      !integer,dimension(KSD)               ::  NBIN,NBIN_cloud
      !real,dimension(KIDIM3,KIDIM2)        ::  RIMAGL,RREALL,RIMAGL_cloud,RREALL_cloud  
      !real,dimension(KIDIM3,KIDIM2)        ::  BRF1_land,BRP1_land,BRM1_water
      !real,dimension(KVERTM,KSD)           ::  H0,H0_cloud
      !real,dimension(KSD)                  ::  C0,C0_cloud
      !real,dimension(KSD)                  ::  sigma_aerosol,sigma_cloud
      !real,dimension(KW)                   ::  CL 
      type(forward_model_characteristics_particles)  :: forw_aerosol
      type(forward_model_characteristics_particles)  :: forw_clouds
      type(forward_model_characteristics_surface)    :: forw_surface
!	------------------------------------------------------------------------------------------------------
      integer	::	II,ISD,IW !, IW1, IP
      integer	::	IDIM1 !, IDIM2, IDIM3
!      integer,dimension(KSD)	::	 NSHAPE
!      integer                 ::  NHVP_retr
!      real,dimension(KVERTM)  ::  HVP_retr_km
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
      real,dimension(KBF)	                ::	BRF_land,BRP_land,BRM_water
      !!real,dimension(KPARS,KSD)	        ::	RADIUS,SD
      !real,dimension(KIDIM3,KSD)	        ::	SD,RADIUS
      !real,dimension(KIDIM3,KSD)	        ::	SD_cloud,RADIUS_cloud
      !real,dimension(KSHAPE,KSD)	        ::	RATIO,SHD,RATIO_cloud,SHD_cloud
      !real,dimension(KSD)	                ::	fcloud
!tl      real,dimension(3,KW),save	        ::	REALM,RIMAGM
!	------------------------------------------------------------------------------------------------------
!      real (selected_real_kind(p=15)),dimension(RIN%NW) :: RREALMM,RIMAGMM
!      real (selected_real_kind(p=15)),dimension(RIN%NW) :: WAVELM
!      real (selected_real_kind(p=15))	::	rh_mxtr       ! Valid range: 0--1
!      real (selected_real_kind(p=15))	::	fract_inslbl  ! Volume fraction of insoluble inclusions
!      real (selected_real_kind(p=15))	::	fract_inslbl1 ! Volume fraction of insoluble inclusions 
!      real (selected_real_kind(p=15))	::	fract_soot    ! Volume fraction of mixture that is composed of soot (aka BC).
!      real (selected_real_kind(p=15))	::	fract_iron    ! Volume fraction of mixture that is composed of iron (aka Hematite).
!      real (selected_real_kind(p=15))	::	fract_slbl    ! Volume fraction of soluble inclusions 
!      real (selected_real_kind(p=15))	::	fract_wtr     ! Volume fraction of mixture that is composed of water.
!      character(12)                    :: instrument
!	------------------------------------------------------------------------------------------------------
      integer	::	IB
!	------------------------------------------------------------------------------------------------------
      real,dimension(KSD)	::	RREAL, RIMAG 
      real,dimension(KSD)	::	RREAL_clouds, RIMAG_clouds 

!      real	   ::	AA, AAT, BB
      !integer	::	par_type,par_type_SD 

      !integer,dimension(KSD)	::	NSHAPE1,NSHAPE1_clouds
!      character*15,save :: ref_ind_flag
      !integer           :: par_type_SD,par_type_SHD,par_type_RI,par_type_SURF
!	------------------------------------------------------------------------------------------------------
      logical		::	status_funct 
      integer, parameter  ::  KANG=2*NG0+1
!	------------------------------------------------------------------------------------------------------
!     variables related to gas absorption
      type(DATA_ABS)  ::  abs_data_forw_im
!	------------------------------------------------------------------------------------------------------
      !RADIUS(:,:) = 0.0
      !SD(:,:)     = 0.0
      !RATIO(:,:)  = 0.0
      !SHD(:,:)    = 0.0
      ANGL(:)     = 0.0
      !H0(:,:)     = 0.0
      !RREAL(:)    = 0.0
      !RIMAG(:)    = 0.0
      
      !BRF(:)    = 0.0
      !BRP(:)    = 0.0
      !BRF1_land(:,:) = 0.0
      !BRP1_land(:,:) = 0.0
!	------------------------------------------------------------------------------------------------------
      HOBS_km     = pixel_fit%HOBS * 0.001
      HGR_km      = pixel_fit%MASL * 0.001
      HMAX_atm_km = HMAX_atm       * 0.001

      if (NHVP_meas .gt. 1) then
         HVP_meas_km(1:NHVP_meas) = HVP_meas(1:NHVP_meas) * 0.001
      end if
!	------------------------------------------------------------------------------------------------------
      select case(RIN%KL) 
      case(1)
        APSING(1:RIN%KNSING) = EXP(AP(1:RIN%KNSING))
      case(0) 
        APSING(1:RIN%KNSING) =     AP(1:RIN%KNSING)
      end select
!tl      DO II=1,RIN%KNSING
!tl      IF(APSING(II) .LT. RIN%APSMIN(II))  APSING(II) = RIN%APSMIN(II)
!tl      IF(APSING(II) .GT. RIN%APSMAX(II))  APSING(II) = RIN%APSMAX(II)
!tl      ENDDO ! II
	   			   
!           write(*,*) ' 3: forward_model_pixel APSMIN: '
!           write(*,'(10e13.4)') RIN%APSMIN(1:RIN%KNSING)
!           write(*,*) ' 3: forward_model_pixel APSMAX: '
!           write(*,'(10e13.4)') RIN%APSMAX(1:RIN%KNSING)

!al      write(*,*) 'in forward_model_pixel AP,APSING: '
!al      do ii=1,RIN%KNSING
!al      write(*,*) ii,exp(AP(ii)),APSING(ii)
!al      enddo 

! Unpack parameter vector AP (driving forward model)
      call unpack_parameter_vector_ap(  RIN, APSING, forw_aerosol, forw_clouds, forw_surface, &
                                        pixel_fit=pixel_fit, &      
                                        HGR_km=HGR_km, NHVP_meas=NHVP_meas, HVP_meas_km=HVP_meas_km )
goto 111 ! test

Write(*,*) 'Values after unpack_parameter_vector_ap '
! Aerosol
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'RADIUS: isd=',isd,'  C0=',forw_aerosol%C0(isd)
  write(*,'(10e14.6)') forw_aerosol%RADIUS(1:forw_aerosol%NBIN(isd),isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'SD:     isd=',isd
  write(*,'(10e14.6)') forw_aerosol%SD(1:forw_aerosol%NBIN(isd),isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'NSHAPE1: isd=',isd
  write(*,'(2i5)') forw_aerosol%NSHAPE(isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'RATIO:     isd=',isd
  write(*,'(10e14.6)') forw_aerosol%RATIO(1:forw_aerosol%NSHAPE(isd),isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'SHD:     isd=',isd
  write(*,'(10e14.6)') forw_aerosol%SHD(1:forw_aerosol%NSHAPE(isd),isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'H0: isd=',isd,'  NHVP_retr=',forw_aerosol%NHVP_retr
  write(*,'(10e14.6)') forw_aerosol%H0(1:forw_aerosol%NHVP_retr,isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'prifile std: isd=',isd
  write(*,'(10e14.6)') forw_aerosol%sigma(isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'RREALL: isd=',isd
  write(*,'(10e14.6)') forw_aerosol%RREAL(1:RIN%NW,isd)
enddo
write(*,*) '******** in forward_model_pixel'
do isd=1,RIN%NSD
  write(*,*) 'RIMAGL: isd=',isd
  write(*,'(10e14.6)') forw_aerosol%RIMAG(1:RIN%NW,isd)
enddo

111 continue

goto 222

! Surface
write(*,*) '******** in forward_model_pixel'
do ii=1,4
  write(*,*) 'BRF1_land:'
  write(*,'(10e14.6)') forw_surface%BRF_land(1:RIN%NW,ii)
enddo
write(*,*) '******** in forward_model_pixel'
do ii=1,4
  write(*,*) 'BRP1_land:'
  write(*,'(10e14.6)') forw_surface%BRP_land(1:RIN%NW,ii)
enddo
do ii=1,4
  write(*,*) 'BRM1_water:'
  write(*,'(10e14.6)') forw_surface%BRM_water(1:RIN%NW,ii)
enddo

222 continue

goto 333 ! test

!if(RIN%NSD_clouds .ne. 0) then
!! Clouds
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'RADIUS_clouds: isd=',isd,'  C0_clouds=',forw_clouds%C0(isd)
!  write(*,'(10e14.6)') forw_clouds%RADIUS(1:forw_clouds%NBIN(isd),isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'SD_clouds:     isd=',isd
!  write(*,'(10e14.6)') forw_clouds%SD(1:forw_clouds%NBIN(isd),isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'NSHAPE1_clouds: isd=',isd
!  write(*,'(2i5)') forw_clouds%NSHAPE(isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'RATIO_clouds:     isd=',isd
!  write(*,'(10e14.6)') forw_clouds%RATIO(1:forw_clouds%NSHAPE(isd),isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'SHD_clouds:     isd=',isd
!  write(*,'(10e14.6)') forw_clouds%SHD(1:forw_clouds%NSHAPE(isd),isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'H0_clouds:     isd=',isd
!  write(*,'(10e14.6)') forw_clouds%H0(1:forw_clouds%NHVP_retr,isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'RREALL_clouds: isd=',isd
!  write(*,'(10e14.6)') forw_clouds%RREAL(1:RIN%NW,isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'RIMAGL_clouds: isd=',isd
!  write(*,'(10e14.6)') forw_clouds%RIMAG(1:RIN%NW,isd)
!enddo
!write(*,*) '******** in forward_model_pixel'
!do isd=1,RIN%NSD_clouds
!  write(*,*) 'fract_clouds:     isd=',isd
!  write(*,'(10e14.6)') forw_clouds%fclouds(isd)
!enddo
!endif ! RIN%NSD_cloud .ne. 0

333 continue


!     Calculate Optical Characteristics
      LOOP_WL:  DO IW=IWb,IWe
!tl      IF(IWW .NE. 0 .AND. IWW .NE. ind_wl(IW)) CYCLE LOOP_WL 
! Surface parameters at IW-th wave length 
      call get_SURF_wl( RIN,IW,ind_wl(IW),  & 
                        forw_surface%BRF_land,forw_surface%BRP_land,forw_surface%BRM_water,  & 
                        BRF_land,BRP_land,BRM_water )
      
! Complex refractive index at IW-th wave length
      if(RIN%NSD .ne. 0) &
      call get_REFI_wl( RIN,GOUT_aerosol%chem%pixel(ipix),IW,ind_wl(IW),  &
                        forw_aerosol%RREAL,forw_aerosol%RIMAG,RREAL,RIMAG )
!      if(RIN%NSD_clouds .ne. 0) &
!      call get_REFI_wl( RIN,GOUT_clouds%chem%pixel(ipix),IW,ind_wl(IW),ind_wl(1),  &
!                        forw_clouds%RREAL,forw_clouds%RIMAG,RREAL_clouds,RIMAG_clouds)
!write(*,*) 'in forward_model_pixel: iw     =',iw,'  WAVE(IW)    =',RIN%WAVE(IW)
!write(*,*) 'in forward_model_pixel: iw     =',iw,'  wl(IW)      =',wl(IW)									
!write(*,*) 'in forward_model_pixel: RREAL:',RREAL(1:RIN%NSD)
!write(*,*) 'in forward_model_pixel: RIMAG:',RIMAG(1:RIN%NSD)

!      call forward_model_pixel_wl ( RIN,ipix,                         & ! IN
!                                    IW,wl(IW),ind_wl(IW),lresult,     &
!                                    forw_aerosol%NBIN,forw_aerosol%RADIUS,forw_aerosol%SD, &
!                                    forw_aerosol%SHD,forw_aerosol%NSHAPE,            &
!                                    RREAL,RIMAG,                                     &
!                                    forw_clouds%NBIN,forw_clouds%RADIUS,forw_clouds%SD, &
!                                    forw_clouds%SHD,forw_clouds%NSHAPE, &
!                                    RREAL_clouds,RIMAG_clouds,          &
!                                    OSHP,                             &
!                                    BRF_land,BRP_land,BRM_water,tau_mol(ind_wl(IW)), &
!                                    HOBS_km,HGR_km,HMAX_atm_km,       &
!                                    NHVP_meas,HVP_meas_km,            &
!                                    forw_aerosol%NHVP_retr,forw_aerosol%HVP_retr_km,         &
!                                    forw_aerosol%H0,forw_aerosol%sigma,forw_aerosol%CL(IW),  &
!                                    forw_clouds%H0,forw_clouds%sigma,forw_clouds%fclouds,  &
!                                    abs_data_forw_im,               &
!                                    pixel_fit,                      & ! INOUT
!                                    GOUT_aerosol,                   &
!                                    GOUT_clouds,                    &
!                                    GOUT_surface,                   &
!
!                                    NANG=NANG,ANGL=ANGL,            & ! OUT
!                                    KERNELS1=KERNELS1,              &
!                                    KERNELS2=KERNELS2               &
!                                  )

         CALL forward_model_pixel_wl (  RIN,ipix,                                        & ! IN
                                        IW,wl(IW),ind_wl(IW),lresult,                    &
                                        forw_aerosol%NBIN,forw_aerosol%RADIUS,forw_aerosol%SD, &
                                        forw_aerosol%SHD,forw_aerosol%NSHAPE,            &
                                        RREAL,RIMAG,                                     &
                                        OSHP,                                            &
                                        BRF_land,BRP_land,BRM_water,tau_mol(ind_wl(IW)), &
                                        HOBS_km,HGR_km,HMAX_atm_km,                      &
                                        NHVP_meas,HVP_meas_km,                           &
                                        forw_aerosol%NHVP_retr,forw_aerosol%HVP_retr_km,         &
                                        forw_aerosol%H0,forw_aerosol%sigma,forw_aerosol%CL(IW),  &
                                        abs_data_forw_im,                                &
                                        pixel_fit,                                       & ! INOUT
                                        GOUT_aerosol,                                    &
                                        GOUT_clouds,                                     &
                                        GOUT_surface,                                    &

                                        NANG=NANG,ANGL=ANGL,                             & ! OUT
                                        KERNELS1=KERNELS1,                               &
                                        KERNELS2=KERNELS2                                &
                                     )
        if ( stop_report%status ) return
!write(*,*) 'after  forward_model_pixel: iw =',iw,'  iwb =',iwb,'  iwe =',iwe
      ENDDO LOOP_WL  ! IW
!	-------------------------------------------------------------------------------------
      if (.not. lresult)  &
         CALL set_pixel_meas_vector_FS (  RIN,IWb,IWe,ipix, & ! IN
                                          pixel_fit,        & ! INOUT
                                          pixel_vec_fit     &
                                       )
		 
!      IB=SUM(pixel_vec_fit%nFS(1:IWb))-pixel_vec_fit%nFS(IWb) ! number of elements in FPS before IWb-th wavelength 
!      II=SUM(pixel_vec_fit%nFS(IWb:IWe))     ! number of elements in FPS for wavelengths from IWb to IWe
!      do iw=1,pixel_vec_fit%KMIMAGE
!         write(*,*) iw,FPS(iw),'  - iw,pixel_vec_fit%FS(iw)'
!      enddo ! iw
!      write(*,*) IB,II,pixel_vec_fit%KMIMAGE,'  IB,II,KMIMAGE'

!      status_funct = check_nan(II,pixel_vec_fit%FS(IB+1:IB+II))
!      if(.not. status_funct) then
!		   write(*,*) 'NaN: after set_pixel_meas_vector_FS in forward_model_pixel'
!		   stop 'stop in forward_model_pixel'
!      endif
!      stop 'stop test after check_nan in forward_model_pixel_wl'

!	-------------------------------------------------------------------------------------

!tl      IF(RIN%KL .EQ. 1) AP(1:RIN%KNSING) = LOG(APSING(1:RIN%KNSING))
!tl      IF(RIN%KL .EQ. 0) AP(1:RIN%KNSING) =     APSING(1:RIN%KNSING)
!
      RETURN
      END SUBROUTINE forward_model_pixel

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine forw_single_scattering_particle_properties ( RIN,NSD,NBIN,RADIUS,SD,      & ! IN
                                                              SHD,NSHAPE1,                 &
                                                              IW,WAVE,RREAL,RIMAG,         &
                                                              NANG,ANGL,ipix,              & ! OUT
                                                              GOUT_particles,              &
                                                              KERNELS1,KERNELS2            & 
                                                            )
      use mod_par_DLS,   only : KMpar
      use mod_par_inv,   only : KSHAPE,KIDIM3 
      use mod_par_OS,    only : KSD
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
      use mod_alloc_arrays
      	  	  
      implicit none
!	------------------------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(in)  ::  RIN
      integer,                    intent(in)  ::  ipix,IW,NSD
      real,                       intent(in)  ::  WAVE
      integer,dimension(KSD),     intent(in)  ::  NSHAPE1  
      real,dimension(KSD),        intent(in)  ::  RREAL,RIMAG
      real,dimension(KIDIM3,KSD), intent(in)  ::  RADIUS,SD	  
      real,dimension(KSHAPE,KSD), intent(in)  ::  SHD  
      integer,dimension(KSD),     intent(in)  ::  NBIN
!	------------------------------------------------------------------------------------------------------
      type(output_segment_particles),intent(inout)  ::  GOUT_particles
      type(kernels_triangle_bin),    intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),   intent(inout)  ::  KERNELS2
!	------------------------------------------------------------------------------------------------------
      integer,                    intent(out)  ::  NANG
      real,dimension(KMpar),      intent(out)  ::  ANGL      
!	------------------------------------------------------------------------------------------------------
      integer	                      ::	ISD
      real                          ::  sca
      logical                       ::  use_models
!	------------------------------------------------------------------------------------------------------  
!	------------------------------------------------------------------------------------------------------
      use_models = RIN%use_models

      call forw_phase_matrix (                               &
                                NSD,NBIN,RADIUS,SD,          & ! IN
                                SHD,NSHAPE1,                 &
                                IW,WAVE,RREAL,RIMAG,         &
                                use_models,                  &
                                NANG,ANGL,ipix,              & ! OUT
                                GOUT_particles%opt%pixel(ipix)%wl(IW),  &
                                GOUT_particles%phmx%pixel(ipix)%wl(IW), &
                                RIN%DLSF,KERNELS1,KERNELS2   & 
                              )
      if ( stop_report%status ) return
      ! Set optical properties for output
      do ISD=1,NSD
        GOUT_particles%rind%pixel(ipix)%wl(IW)%mreal(ISD) = RREAL(ISD)
        GOUT_particles%rind%pixel(ipix)%wl(IW)%mimag(ISD) = RIMAG(ISD)

        sca = GOUT_particles%opt%pixel(ipix)%wl(IW)%ssa(ISD)*GOUT_particles%opt%pixel(ipix)%wl(IW)%ext(ISD)
        GOUT_particles%opt%pixel(ipix)%wl(IW)%aext(ISD) = &
                    GOUT_particles%opt%pixel(ipix)%wl(IW)%ext(ISD) * &
                    (1.-GOUT_particles%opt%pixel(ipix)%wl(IW)%ssa(ISD))
        if(RIN%DLSF%keyEL .gt. 0) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph11(1:NANG,ISD) =  &
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph11(1:NANG,ISD)/sca
        if(RIN%DLSF%keyEL .gt. 1) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph12(1:NANG,ISD) =  & 
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph12(1:NANG,ISD)/sca
        if(RIN%DLSF%keyEL .gt. 2) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph22(1:NANG,ISD) =  & 
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph22(1:NANG,ISD)/sca
        if(RIN%DLSF%keyEL .gt. 3) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph33(1:NANG,ISD) =  & 
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph33(1:NANG,ISD)/sca
        if(RIN%DLSF%keyEL .gt. 4) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph34(1:NANG,ISD) =  & 
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph34(1:NANG,ISD)/sca
        if(RIN%DLSF%keyEL .gt. 5) GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph44(1:NANG,ISD) =  & 
                                  GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph44(1:NANG,ISD)/sca       
        GOUT_particles%rind%pixel(ipix)%wl(IW)%mreal(ISD) = RREAL(ISD)
        GOUT_particles%rind%pixel(ipix)%wl(IW)%mimag(ISD) = RIMAG(ISD)
      enddo ! ISD
      if(RIN%DLSF%keyEL .gt. 0) &
        GOUT_particles%lidar%pixel(ipix)%wl(IW)%lr(1:NSD)= 4.*3.14159265/  &
        ( GOUT_particles%opt%pixel(ipix)%wl(IW)%ssa(1:NSD)*  &
          GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph11(NANG,1:NSD) )
      if(RIN%DLSF%keyEL .gt. 2) &
        GOUT_particles%lidar%pixel(ipix)%wl(IW)%ldpr(1:NSD) =          &
        (1.-GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph22(NANG,1:NSD))/  &
        (1.+GOUT_particles%phmx%pixel(ipix)%wl(IW)%ph22(NANG,1:NSD))*100.

      return
      end subroutine forw_single_scattering_particle_properties

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      end module mod_forward_model
