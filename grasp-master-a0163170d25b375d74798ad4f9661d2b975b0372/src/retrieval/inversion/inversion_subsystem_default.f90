module inversion_subsystem
    use mod_forward_model
    use mod_sdata_derived_type
    use mod_retr_settings_derived_type
    use mod_retr_general_output_derived_type
    use mod_alloc_kernels
    use mod_par_DLS,                                   only : KMPAR
    use mod_par_inv,                                   only : KW, KPARS, KIMAGE
    use mod_sdata,                                     only : get_pixel_wl
    use mod_stop_report

    implicit none

#ifdef MUMPS
include 'mpif.h'
#endif

  contains


    subroutine inversion_init(                         &
            input_settings                             &
        )

        type(retr_input_settings),intent(in) :: input_settings
        integer                              :: IERR

#ifdef MUMPS
        call MPI_INIT(IERR)
#endif
    end subroutine inversion_init


    subroutine inversion_finalize()
        integer :: IERR

#ifdef MUMPS
        call MPI_FINALIZE(IERR)
#endif
    end subroutine inversion_finalize


    subroutine inversion_forward_model(                &
            ipixstart,                                 &
            ipixstop,                                  &
            RIN,                                       &
            OSHP,                                      &
            IWW,                                       &
            lresult,                                   &
            tau_mol,                                   &
            NHVP_meas,                                 &
            HVP_meas,                                  &
            AP,                                        &
            pixel_fit,                                 &
            pixel_vec_fit,                             &
            GOUT_aerosol,                              &
            GOUT_clouds,                               &
            GOUT_surface,                              &
            NANG,                                      &
            ANGL,                                      &
            KERNELS1,                                  &
            KERNELS2                                   &
        )

        integer,                              intent(in)       :: ipixstart
        integer,                              intent(in)       :: ipixstop

        type(retr_input_settings),            intent(in)       :: RIN
        type(OSH_par),                        intent(in)       :: OSHP
        integer,                              intent(in)       :: IWW
        logical,                              intent(in)       :: lresult
        real,dimension(KW,KIMAGE),            intent(in)       :: tau_mol
        integer,                              intent(in)       :: NHVP_meas
        real,dimension(KVERTM),               intent(in)       :: HVP_meas
        real,dimension(KPARS,KIMAGE),         intent(in)       :: AP
        type(pixel),dimension(KIMAGE),        intent(inout)    :: pixel_fit
        type(pixel_vector),dimension(KIMAGE), intent(inout)    :: pixel_vec_fit
        type(output_segment_particles),       intent(inout)    :: GOUT_aerosol
        type(output_segment_particles),       intent(inout)    :: GOUT_clouds
        type(output_segment_surface),         intent(inout)    :: GOUT_surface
        integer,                              intent(out)      :: NANG
        real,dimension(KMPAR),                intent(out)      :: ANGL
        type(kernels_triangle_bin),           intent(inout)    :: KERNELS1
        type(kernels_lognormal_bin),          intent(inout)    :: KERNELS2

        integer                                                :: IWb
        integer                                                :: IWe
        integer,dimension(KIMAGE)                              :: nwl
        real,dimension(KW,KIMAGE)                              :: wl
        integer,dimension(KW,KIMAGE)                           :: ind_wl
        integer                                                :: ipix
        integer                                                :: i


        if(.not. lresult) then
          do ipix=ipixstart,ipixstop
            call get_pixel_wl(                    &
                                pixel_fit(ipix),  &
                                nwl(ipix),        &
                                wl(:,ipix),       &
                                ind_wl(:,ipix)    &
                             )
          enddo
        else
          do ipix=ipixstart,ipixstop
            nwl(ipix) = RIN%NW
            wl(1:RIN%NW,ipix) = RIN%WAVE(1:RIN%NW)
            do i=1,RIN%NW
              ind_wl(i,ipix) = i
            enddo
          enddo
        endif

        IWb = 1
        do ipix=ipixstart,ipixstop
            IWe = nwl(ipix)
            call forward_model_pixel(                               &
                    RIN,                                            &
                    OSHP,                                           &
                    ipix,                                           &
                    IWb,                                            &
                    IWe,                                            &
                    0,                                              &
                    lresult,                                        &
                    tau_mol(:,ipix),                                &
                    NHVP_meas,                                      &
                    HVP_meas,                                       &
                    nwl(ipix),                                      &
                    wl(:,ipix),                                     &
                    ind_wl(:,ipix),                                 &
                    AP(:,ipix),                                     &
                    pixel_fit(ipix),                                &
                    pixel_vec_fit(ipix),                            &
                    GOUT_aerosol,                                   &
                    GOUT_clouds,                                    &
                    GOUT_surface,                                   &
                    NANG,                                           &
                    ANGL,                                           &
                    KERNELS1,                                       &
                    KERNELS2                                        &
                )
        if ( stop_report%status ) return
        enddo
    end subroutine inversion_forward_model


    subroutine inversion_jacobian_matrix(                   &
            iu_main_output,                                 &
            lp,                                             &
            ipixstart,                                      &
            ipixstop,                                       &
            RIN,                                            &
            lresult,                                        &
            tau_mol,                                        &
            NHVP_meas,                                      &
            HVP_meas,                                       &
            AP,                                             &
            ALS,                                            &
            pixel_fit,                                      &
            segment_vec_fit,                                &
            GOUT_aerosol,                                   &
            GOUT_clouds,                                    &
            GOUT_surface,                                   &
            NANG,                                           &
            ANGL,                                           &
            KERNELS1,                                       &
            KERNELS2,                                       &
            APMIN,                                          &
            APMAX,                                          &
            dermask,                                        &
            US                                              &
        )

        use mod_covariance_matrix,                         only : delta_derivatives

        integer,                              intent(in)       :: iu_main_output
        integer,                              intent(in)       :: lp
        integer,                              intent(in)       :: ipixstart
        integer,                              intent(in)       :: ipixstop

        type(retr_input_settings),            intent(in)       :: RIN
        logical,                              intent(in)       :: lresult
        real,dimension(KW,KIMAGE),            intent(in)       :: tau_mol
        integer,                              intent(in)       :: NHVP_meas
        real,dimension(KVERTM),               intent(in)       :: HVP_meas
        real,dimension(KPARS,KIMAGE),         intent(in)       :: AP
        real,dimension(KIMAGE),               intent(in)       :: ALS
        type(pixel),dimension(KIMAGE),        intent(in)       :: pixel_fit
        type(pixel_vector),dimension(KIMAGE), intent(in)       :: segment_vec_fit
        type(output_segment_particles),       intent(inout)    :: GOUT_aerosol
        type(output_segment_particles),       intent(inout)    :: GOUT_clouds
        type(output_segment_surface),         intent(inout)    :: GOUT_surface
        integer,                              intent(out)      :: NANG
        real,dimension(KMPAR),                intent(out)      :: ANGL
        type(kernels_triangle_bin),           intent(inout)    :: KERNELS1
        type(kernels_lognormal_bin),          intent(inout)    :: KERNELS2
        real,dimension(KPARS),                intent(in)       :: APMIN
        real,dimension(KPARS),                intent(in)       :: APMAX
        logical,dimension(KPARS,KIMAGE),      intent(in)       :: dermask
        real,dimension(KMESS,KPARS,KIMAGE),   intent(inout)    :: US

        integer                                                :: ipix
        integer                                                :: I
        integer                                                :: J
        integer                                                :: IW
        integer                                                :: IWb
        integer                                                :: IWe
        integer                                                :: J1
        real                                                   :: DL1
        real,dimension(KPARS)                                  :: AP_temp
        type(pixel)                                            :: pixel_fit_temp
        type(pixel_vector)                                     :: pixel_vec_fit_temp
        type(pixel_vector),dimension(KIMAGE)                   :: pixel_vec_fit_deriv
        integer                                                :: nwl
        real,dimension(KW)                                     :: wl
        integer,dimension(KW)                                  :: ind_wl

        pixel_vec_fit_deriv(ipixstart:ipixstop) = segment_vec_fit(ipixstart:ipixstop)
        do ipix=ipixstart,ipixstop
            ! Jacobian matrix for each pixel (Kp);  US(:,:,ipix)
            AP_temp(1:RIN%KNSING) = AP(1:RIN%KNSING,ipix)
            pixel_fit_temp = pixel_fit(ipix)
            pixel_vec_fit_temp = pixel_vec_fit_deriv(ipix)
            call get_pixel_wl(                    &
                                pixel_fit(ipix),  &
                                nwl,              &
                                wl,               &
                                ind_wl            &
                             )
            LOOP_WL: do IW=1,nwl
                J1=SUM(segment_vec_fit(ipix)%nFS(1:IW))-segment_vec_fit(ipix)%nFS(IW) ! number of elements in FPS before IW-th wavelength
                IWb=IW
                IWe=IW
                LOOP_parameters: do I=1,RIN%KNSINGF
                    if(RIN%IWW_SINGL(I) .ne. 0 .and. RIN%IWW_SINGL(I) .ne. ind_wl(IW)) then
                        cycle LOOP_parameters
                    endif
                    if(.not. dermask(I,ipix)) then
                        cycle LOOP_parameters
                    endif
                    DL1 = delta_derivatives(ALS(ipix),APMIN(I),APMAX(I),AP(I,ipix),RIN%DL)
                    AP_temp(I) = AP(I,ipix)+DL1
                    if(I .eq. -24) then
                        write(*,'(3(3x,i0),3e16.7,3x,a)') lp, ipix, i, exp(AP_temp(i)), exp(AP(i,ipix)), DL1, &
                                                         'lp, ipix, i, AP_temp(i), AP(i), DL1  - in inversion_jacobian_matrix'
                    endif

                    call forward_model_pixel(          &
                            RIN,                       &
                            RIN%OSHD,                  &
                            ipix,                      &
                            IWb,                       &
                            IWe,                       &
                            RIN%IWW_SINGL(I),          &
                            lresult,                   &
                            tau_mol(:,ipix),           &                                                          
                            NHVP_meas,                 &
                            HVP_meas,                  &
                            nwl,                       &
                            wl,                        &
                            ind_wl,                    &
                            AP_temp,                   &
                            pixel_fit_temp,            &
                            pixel_vec_fit_temp,        & 
                            GOUT_aerosol,              &
                            GOUT_clouds,               &         
                            GOUT_surface,              &
                            NANG,                      &
                            ANGL,                      &
                            KERNELS1,                  &
                            KERNELS2                   &
                        )
                    if ( stop_report%status ) return
                    do J=J1+1,J1+segment_vec_fit(ipix)%nFS(IW)
                        US(J,I,ipix)=(pixel_vec_fit_temp%FS(J)-pixel_vec_fit_deriv(ipix)%FS(J))/DL1
                    enddo ! J
                    if(ipix .eq. -7) then  !(I .eq. -24)
                        do J=J1+1,J1+segment_vec_fit(ipix)%nFS(IW)
                            write(*,'(2i5,7e18.8,2a)') I,J,US(J,I,ipix),DL1,RIN%DL, &
                                pixel_vec_fit_temp%FS(J),pixel_vec_fit_deriv(ipix)%FS(J), exp(pixel_vec_fit_temp%FS(J)), &
                                exp(pixel_vec_fit_deriv(ipix)%FS(J)),'  I,J,US(J,I),DL1,DL,FPS1(J),', &
                                'FPS(J),exp(FPS1(J)),exp(FPS(J))'
                        enddo
                     endif ! I .eq. -1
                    AP_temp(I) = AP(I,ipix)
                enddo LOOP_parameters 
            enddo LOOP_WL
            if(RIN%IPRI_verbose) then
              if(RIN%IPRI_additional_info .or. RIN%IPFP%INVSING .gt. 0) then
                write(iu_main_output,'(10x,a,i0,3x,a,i0)') 'pixel # ',ipix, &
                'AFTER Jacobian Matrix     for iteration # ',lp+1
              endif
            endif
        enddo
    end subroutine inversion_jacobian_matrix

end module inversion_subsystem
