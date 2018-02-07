! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../constants_set/mod_globals.inc"
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine OPT_ERR_estimates (iu_main_output,         & ! IN
                                    sparse_solver_name,     &
                                    UF, UFS, nnz_err, UFNZ,       &
                                    QS, QIM, AGENP, ALS, KM1_pix, &
                                    APMIN, APMAX,                 &
                                    RIN, INVSING, segment_meas,   &
                                    errest_pix, AP,               &
                                    tau_mol, NHVP_meas, HVP_meas, &
                                    GOUT,                   &
                                    solver_timer,           &
                                    KERNELS1,KERNELS2       &
                                   )

      use mod_par_inv, only : KPARS,KIMAGE,KPAR,KW
      use mod_par_OS,  only : KSD
      use mod_par_DLS, only : KMpar
      use mod_fisher_matrix_ccs    
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type 
      use mod_alloc_kernels
      use mod_sdata_derived_type
      use mod_forward_model
      use mod_covariance_matrix, only : delta_derivatives
      use mod_stop_report

      !use iso_c_binding ! c_intptr_t
        
      implicit none
!	----------------------------------------------------------------------
! IN :
      type(retr_input_settings),         intent(in)  ::  RIN
      type(segment_data),                intent(in)  ::  segment_meas
      integer,                           intent(in)  ::  iu_main_output
      character(*),                      intent(in)  ::  sparse_solver_name
      integer,                           intent(in)  ::  INVSING
      integer,                           intent(in)  ::  nnz_err
      type(nonzero),                     intent(inout)  ::  UFNZ    
      real,                              intent(in)  ::  AGENP 
      real,dimension(KIMAGE),            intent(in)  ::  ALS
      integer,dimension(KIMAGE),         intent(in)  ::  KM1_pix
      
      real,dimension(KPARS,KPARS,KIMAGE),intent(inout)  ::  UFS     
      real,dimension(KPARS,KIMAGE),      intent(in)  ::  QS
      real,dimension(KPAR),              intent(in)  ::  QIM    
      real,dimension(KPAR,KPAR),         intent(in)  ::  UF
      logical, dimension(KIMAGE),        intent(in)  ::  errest_pix
      real,dimension(KPARS),             intent(in)  ::  APMIN,APMAX
      real,dimension(KW,KIMAGE),         intent(in)  ::  tau_mol
      integer,                           intent(in)  ::  NHVP_meas
      real,dimension(KVERTM),            intent(in)  ::  HVP_meas  ! heights for lidar measurements
!	----------------------------------------------------------------------
! INOUT :
      real,                              intent(inout)  ::  solver_timer
      real,dimension(KPARS,KIMAGE),      intent(inout)  ::  AP
      type(output_segment_general),      intent(inout)  ::  GOUT
      type(kernels_triangle_bin),        intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin),       intent(inout)  ::  KERNELS2
!	-------------------------------------------------------------------------------------
! LOCAL :
      type(output_segment_particles)  ::  GOUT_aerosol_temp
      type(output_segment_particles)  ::  GOUT_clouds_temp		
      type(output_segment_surface)    ::  GOUT_surface_temp
      type(output_segment_products_particles)  ::  GOUT_products_aerosol_temp
      real,dimension(KPARS)		        ::  AP_temp
      real,dimension(:,:),allocatable ::  UW0
      real                            ::  DL1
      integer                ::  IWW
      !real*4                 ::  time_start,time_finish
      real, dimension(KPARS) ::  FFS, QS1
      real*8                               ::  TEMP
      real*8, dimension(:,:), allocatable  ::  ATEMP
      real*8, dimension(KPAR) ::  b
      real,   dimension(KPAR) ::  b1
      logical ::  ltest, lresult
      integer ::  alloc_stat
      integer ::  i, i1, j, j1, ii, ii1, IOF, ISD, IW, ipix
      integer ::  IWb, IWe
      real    ::  tiny, ERR_temp, BIAS_temp
      integer            ::  npixels, NW, IMQ, KNSING, KNSINGF
      real,dimension(KW) ::  WAVE
      integer  ::  nnz_pix
      integer  ::  KN, KNF, NOF, NSD1
! NOF  - number of optical functions for error estimates
! NSD1 - number of SD modes + 1(total)
      integer :: NANG
      real,dimension(KMPAR) :: ANGL
      type(pixel) :: pixel_fit_temp
      type(pixel_vector) :: pixel_vec_fit_temp
      integer,dimension(KW) :: ind_wl
!	-------------------------------------------------------------------------------------
      lresult=.true.

      npixels = segment_meas%npixels
      NW     = RIN%nw
      WAVE(:)= RIN%wave(:)
      KNSING = RIN%KNSING
      KNSINGF= RIN%KNSINGF
      IMQ    = RIN%IMQ
      do iw=1,NW
        ind_wl(iw) = iw
      enddo
!write(*,*) 'QS:'
!write(*,'(10e14.5)') QS(1:KNSING,1)      
!write(*,*) 'UFS(i,i):'
!write(*,'(10e14.5)') (UFS(ii,ii,1),ii=1,KNSING)
!write(*,*) 'ALS:'
!write(*,'(10e14.5)') ALS(1:npixels)

      if( RIN%products%errest%aerosol%lidar .and. (RIN%DLSF%keyEL .gt. 0)) then
         NOF=3 ! ext, ssa, lr
      else
         NOF=2 ! ext and ssa
      endif !

      select case(RIN%NSD)
      case(1)
         NSD1=RIN%NSD
      case(2)
         NSD1=RIN%NSD+1  ! fine, coarse, total    
      case default
         write(tmp_message,'(a,i0,a)') &
         'in sub set_UW0 :  NSD = ',RIN%NSD,' - is not supported'
         G_ERROR(trim(tmp_message))
      end select

      allocate(UW0(NW*(RIN%NSD+1)*(RIN%NSD+1),RIN%KNSINGF),stat=alloc_stat)
      if (alloc_stat /= 0) then
        write(tmp_message,'(a)') 'error while trying to allocate UW0'
        G_ERROR(trim(tmp_message))
      endif

      !call CPU_TIME(time_start)
      KNF = KNSINGF*npixels      
      ltest = .false.      
      if(ltest) then
        allocate(ATEMP(KNF,KNF),stat=alloc_stat)
        if (alloc_stat /= 0) then
          write(tmp_message,'(a)') 'error while trying to allocate ATEMP'
          G_ERROR(trim(tmp_message))
        endif
          ATEMP(:,:) = 0.0
          do i=1,nnz_err
            ATEMP(UFNZ%row(i),UFNZ%col(i)) = UFNZ%val(i)
          enddo ! i
      endif ! ltest

      IWb = 1
      IWe = NW    

      do ipix=1,npixels
!C*** FORWARD CALCULATIONS for EACH pixel:  *** 
!         write(*,*) 'before forward_model_pixel_PHMX in opt_err_estimates :'
!         write(*,*) IWb,IWe,ipix,'  - IWb,IWe,ipix'

      if(INVSING .lt. 2) then
        if(.not. errest_pix(ipix)) cycle
        if(RIN%IPRI_verbose .or. RIN%IPRI_additional_info) then
          do ii=1,RIN%KNSINGF
            if(UFS(ii,ii,ipix) .eq. 0.0) then
              !UFS(ii,ii,ipix) = 1e-30
              !if(RIN%IPRI_verbose)  &
              write(tmp_message,'(a,2i4,a)') &
              'Warning OPT_ERR_estimates: UFS diagonal element is ZERO for', &
              ipix,ii,'  ipix,ii'
              G_ERROR(trim(tmp_message))
            endif
          enddo
        endif ! IPRI_verbose .or.
        call UF_nonzero_one_pixel ( KNSINGF,UFS(:,:,ipix), & ! IN
                                    UFNZ,nnz_pix           & ! OUT
                                  )
      if ( stop_report%status ) return
      endif ! INVSING .lt. 2

      pixel_fit_temp = GOUT%retrieval%fit%segment_fit%pixels(ipix)

      UW0(:,:) = 0.0
      AP_temp(1:RIN%KNSING) = AP(1:RIN%KNSING,ipix)

LOOP_WL: do IW=1,NW
      IWb = IW
      IWe = IW	  	        	  
LOOP_parameters: do I=1,RIN%KNSINGF                  
      IWW = RIN%IWW_SINGL(I)
      IF(IWW .NE. 0 .AND. IWW .NE. IW) CYCLE LOOP_parameters
      DL1 = delta_derivatives(ALS(ipix),APMIN(I),APMAX(I),AP(I,ipix),RIN%DL)
!      DL1 = delta_derivatives(ALS(ipix)*KM1_pix(ipix),APMIN(I),APMAX(I),AP(I,ipix),RIN%DL)
      AP_temp(I) = AP(I,ipix)+DL1                                
!!      call MIN_AP_MAX(I,I,APMIN,APMAX,AP_temp)
                    call forward_model_pixel(          &
                            RIN,                       &
                            RIN%OSHD,                  &
                            ipix,                      &
                            IWb,                       &
                            IWe,                       &
                            IWW,                       &
                            lresult,                   &
                            tau_mol(:,ipix),           &                                                          
                            NHVP_meas,                 &
                            HVP_meas,                  &
                            NW,                        &
                            WAVE,                      &
                            ind_wl,                    &
                            AP_temp,                   &
                            pixel_fit_temp,            &
                            pixel_vec_fit_temp,        & 
                            GOUT_aerosol_temp,         &
                            GOUT_clouds_temp,          &
                            GOUT_surface_temp,         &
                            NANG,                      &
                            ANGL,                      &
                            KERNELS1,                  &
                            KERNELS2                   &
                        )

                    call total_single_scattering_properties (        &
                                    iu_main_output, ipix, NANG,      &
                                    RIN%DLSF%keyEL, RIN%NSD, RIN%products%aerosol, &
                                    IW, RIN%WAVE(IW),                &
                                    GOUT_aerosol_temp,               &
                                    GOUT_products_aerosol_temp       &
                                  )

      AP_temp(I) = AP(I,ipix)

      call set_UW0 ( IW,RIN%KL,                    & ! IN
                     NW,RIN%NSD,NOF,NSD1,          & 
                     GOUT%aerosol%opt%pixel(ipix)%wl(IW),        &
                     GOUT%aerosol%lidar%pixel(ipix)%wl(IW),      &
                     GOUT_aerosol_temp%opt%pixel(ipix)%wl(IW),   &
                     GOUT_aerosol_temp%lidar%pixel(ipix)%wl(IW), &
                     DL1,                          &
                     UW0(:,I)                      & ! INOUT
                   )						
enddo LOOP_parameters ! I
enddo LOOP_WL ! IW
	  
!      write(*,*) 'ipix=',ipix,'  AP:'
!      write(*,'(10e14.5)') AP(1:RIN%KNSING,ipix)
!      write(*,*) '1: UW0'
!      write(*,'(10e14.5)') UW0(1:NW*(RIN%NSD+1)*(RIN%NSD+1),1:RIN%KNSINGF)

! opt_funct - ext, ssa, lr         
LOOP_opt_funct: do IOF=1,NOF            

        do ISD=1,NSD1
          do IW=1,NW 
            j1=(IOF-1)*NSD1*NW+(ISD-1)*NW+IW

            select case(INVSING)
            case(2) 
! Multi pixel inversion               
              b(:) = 0.0d0
              do i=1,KNSINGF
                j = KNSINGF*(ipix-1)+i
                b(j) = UW0(j1,i)
                !write(*,*) ipix,i,j,j1,b(j),UW0(j1,i),'  ipix,i,j,j1,b(j),UW0(j1,i)'
              enddo ! i	  
              if(IMQ .eq. 2) then
                call ITERQ (  IMQ,KNF,               & ! IN
                              real(UF(1:KNF,1:KNF)), &
                              real(b(1:KNF)),        &
                              b1(1:KNF)              & ! OUT
                           )
              else if(IMQ .eq. 3 .or. IMQ .eq. 1) then
                ! write(iu_main_output,*) 'in OPT_ERR_estimates, before sparse_solver','  nnz=',nnz_err  
                call sparse_solver (  iu_main_output,     & ! IN
                                      KNF, nnz_err, UFNZ, &
                                      b,                  & ! INOUT
                                      solver_timer        &
                                   )
                ! write(*,*) 'after  solver_SuperLU: j=',j
                b1(1:KNF) = b(1:KNF)
              endif ! IMQ .eq. 2
              ERR_temp  = 0.0
              BIAS_temp = 0.0
              do i=1,KNSINGF
                j = KNSINGF*(ipix-1)+i
                ERR_temp  = ERR_temp +UW0(j1,i)*b1(j)
                BIAS_temp = BIAS_temp+UW0(j1,i)*QIM(j)
                !write(*,*) ipix,i,j,j1,b1(j),UW0(j1,i,ipix),'  ipix,i,j,j1,b1(j),UW0(j1,i,ipix)'
              enddo ! i
              ERR_temp  = sqrt(ERR_temp*AGENP)

            case(0,1)
! Single pixel or Multi pixel inversion + Single pixel

              FFS(:) = 0.0
              do i=1,KNSINGF
                FFS(i) = UW0(j1,i)         
              enddo ! i
              if(IMQ .eq. 2) then
!write(*,*) 'IOF=',IOF,'  ISD=',ISD,'  IW=',IW,'  before ITERQ, NOF=',NOF,'  NSD1=',NSD1
                call ITERQ (  IMQ,KNSINGF,                   & ! IN
                              UFS(1:KNSINGF,1:KNSINGF,ipix), &
                              FFS(1:KNSINGF),                &
                              QS1(1:KNSINGF)                 & ! OUT
                            )
              else if(IMQ .eq. 3) then
                b(1:KNSINGF)=FFS(1:KNSINGF)
                !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
                !'UF inversion, before sparse_solver ',trim(sparse_solver_name),'  nnz_pix=',nnz_pix
!write(*,*) 'IOF=',IOF,'  ISD=',ISD,'  IW=',IW,'  before sparse_solver, NOF=',NOF,'  NSD1=',NSD1
                call sparse_solver (  iu_main_output,       & ! IN
                                      KNSINGF,nnz_pix,UFNZ, &
                                      b,                    & ! INOUT
                                      solver_timer          &
                                    )
                !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
                !'UF inversion, after  sparse_solver ',trim(sparse_solver_name)
                QS1(1:KNSINGF) = b(1:KNSINGF)
!write(*,*) 'after  sparse_solver'
              endif ! IMQ .eq. 2
              ERR_temp  = 0.0
              BIAS_temp = 0.0
              do i=1,KNSINGF
                ERR_temp  = ERR_temp +UW0(j1,i)*QS1(i)
                BIAS_temp = BIAS_temp+UW0(j1,i)*QS(i,ipix)
              enddo ! i
              ERR_temp  = sqrt(ERR_temp*ALS(ipix))
!write(*,*) 'iw=',iw,'  IOF=',IOF,'  ERR_temp=',ERR_temp,'  BIAS_temp=',BIAS_temp
            end select ! INVSING

            call set_opt_err (  IOF,IW,ISD,ipix,    & ! IN
                                ERR_temp,BIAS_temp, &
                                GOUT%errest%aerosol & ! INOUT
                             )      
            if ( stop_report%status ) return
          enddo ! IW
        enddo ! ISD

enddo LOOP_opt_funct

      enddo ! ipix  
      
      if(ltest) then
         deallocate(ATEMP,stat=alloc_stat)
         if (alloc_stat /= 0) then
           write(tmp_message,'(a)') 'error while trying to deallocate ATEMP'
           G_ERROR(trim(tmp_message))
         endif
      endif ! ltest

      !call CPU_TIME(time_finish)

      deallocate(UW0,stat=alloc_stat)
      if (alloc_stat /= 0) then
        write(tmp_message,'(a)') 'error while trying to deallocate UW0'
        G_ERROR(trim(tmp_message))
      endif

      !write(iu_main_output,*) 'OPT_ERR_estimates: CPU_TIME(sec)=',time_finish-time_start
      
      return
      end subroutine OPT_ERR_estimates

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_UW0 ( IW,KL,                & ! IN
                           NW,NSD,NOF,NSD1,      & 
                           GOUT_aerosol_opt_pixel_wl,         &
                           GOUT_aerosol_lidar_pixel_wl,       &
                           GOUT_aerosol_opt_pixel_wl_temp,    &
                           GOUT_aerosol_lidar_pixel_wl_temp,  &
                           DL1,                  &
                           UW0                   & ! INOUT
                         )

      use mod_par_inv, only : KW,KPARS
      use mod_retr_general_output_derived_type

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      integer,                          intent(in)  ::  IW,KL,NW,NSD,NOF,NSD1
      real,                             intent(in)  ::  DL1
      type(output_pixel_opt_wl),        intent(in)  ::  GOUT_aerosol_opt_pixel_wl
      type(output_pixel_opt_wl),        intent(in)  ::  GOUT_aerosol_opt_pixel_wl_temp
      type(output_pixel_lidar_ratio_wl),intent(in)  ::  GOUT_aerosol_lidar_pixel_wl
      type(output_pixel_lidar_ratio_wl),intent(in)  ::  GOUT_aerosol_lidar_pixel_wl_temp     
!	------------------------------------------------------------------------------------------------------
! INOUT :
      real,dimension(NW*(NSD+1)*(NSD+1)), intent(inout)  ::  UW0
!	------------------------------------------------------------------------------------------------------
! LOCAL : 
      integer  ::  J,ISD,IOF

! Error estimates: indexes for UW0 matrix 
!      integer,parameter	::	index_ext = 1
!      integer,parameter	::	index_ssa = 2
!      integer,parameter	::	index_lr  = 3

! NOF  - number of optical functions for error estimates
! NSD1 - number of SD mode + 1(total) 
!	------------------------------------------------------------------------------------------------------
      
      do IOF=1,NOF            

         do ISD=1,NSD
         J=(IOF-1)*NSD1*NW+(ISD-1)*NW+IW
         select case(IOF)
         case(index_ext) 
            if(KL .eq. 1) then
               UW0(J) = (log(GOUT_aerosol_opt_pixel_wl%ext(ISD))-log(GOUT_aerosol_opt_pixel_wl_temp%ext(ISD)))/DL1
            else
               UW0(J) = (GOUT_aerosol_opt_pixel_wl%ext(ISD)-GOUT_aerosol_opt_pixel_wl_temp%ext(ISD))/DL1            
            endif
         case(index_ssa) 
            if(KL .eq. 1) then
               UW0(J) = (log(GOUT_aerosol_opt_pixel_wl%ssa(ISD))-log(GOUT_aerosol_opt_pixel_wl_temp%ssa(ISD)))/DL1
            else
               UW0(J) = (GOUT_aerosol_opt_pixel_wl%ssa(ISD)-GOUT_aerosol_opt_pixel_wl_temp%ssa(ISD))/DL1            
            endif
         case(index_lr) 
            if(KL .eq. 1) then
               UW0(J) = (log(GOUT_aerosol_lidar_pixel_wl%lr(ISD))-log(GOUT_aerosol_lidar_pixel_wl_temp%lr(ISD)))/DL1
            else
               UW0(J) = (GOUT_aerosol_lidar_pixel_wl%lr(ISD)-GOUT_aerosol_lidar_pixel_wl_temp%lr(ISD))/DL1            
            endif
         end select
         enddo ! ISD
                  
         if(NSD1 .eq. 3) then
         J=(IOF-1)*NSD1*NW+(NSD1-1)*NW+IW

         select case(IOF)
         case(index_ext) 
            if(KL .eq. 1) then
               UW0(J) = (log(GOUT_aerosol_opt_pixel_wl%extt)-log(GOUT_aerosol_opt_pixel_wl_temp%extt))/DL1
            else
               UW0(J) = (GOUT_aerosol_opt_pixel_wl%extt-GOUT_aerosol_opt_pixel_wl_temp%extt)/DL1
            endif
         case(index_ssa)
            if(KL .eq. 1) then         
               UW0(J) = (log(GOUT_aerosol_opt_pixel_wl%ssat)-log(GOUT_aerosol_opt_pixel_wl_temp%ssat))/DL1
            else
               UW0(J) = (GOUT_aerosol_opt_pixel_wl%ssat-GOUT_aerosol_opt_pixel_wl_temp%ssat)/DL1            
            endif
         case(index_lr)
            if(KL .eq. 1) then
               UW0(J) = (log(GOUT_aerosol_lidar_pixel_wl%lrt)-log(GOUT_aerosol_lidar_pixel_wl_temp%lrt))/DL1 
            else
               UW0(J) = (GOUT_aerosol_lidar_pixel_wl%lrt-GOUT_aerosol_lidar_pixel_wl_temp%lrt)/DL1
            endif        
         end select
         endif !  NSD1 .eq. 3  

      enddo ! IOF

      return
      end subroutine set_UW0
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_opt_err (IOF,IW,ISD,ipix,          & ! IN
                              ERR_temp,BIAS_temp,       &
                              GOUT_errest_aerosol       & ! INOUT
                             )

      use mod_retr_general_output_derived_type
      use mod_stop_report

      implicit none
!	------------------------------------------------------------------------------------------------------
! IN :
      integer,                   intent(in)  ::  IOF,IW,ISD,ipix
      real,                      intent(in)  ::  ERR_temp,BIAS_temp     
!	------------------------------------------------------------------------------------------------------
! INOUT :
      type(output_segment_err_estim_particles), intent(inout)  ::  GOUT_errest_aerosol

!	------------------------------------------------------------------------------------------------------
! Error estimates: indexes for UW0 matrix 
!      integer,parameter	::	index_ext = 1
!      integer,parameter	::	index_ssa = 2
!      integer,parameter	::	index_lr  = 3

! NOF  - number of optical functions for error estimates
!	------------------------------------------------------------------------------------------------------      

      if(ISD .eq. 1 .or. ISD .eq. 2) then

         select case(IOF)
         case(index_ext)
            GOUT_errest_aerosol%opt%pixel(ipix)%ERR_ext(IW,ISD)  = ERR_temp
            GOUT_errest_aerosol%opt%pixel(ipix)%BIAS_ext(IW,ISD) = BIAS_temp
         case(index_ssa) 
            GOUT_errest_aerosol%opt%pixel(ipix)%ERR_ssa(IW,ISD)  = ERR_temp
            GOUT_errest_aerosol%opt%pixel(ipix)%BIAS_ssa(IW,ISD) = BIAS_temp
         case(index_lr)  
            GOUT_errest_aerosol%lidar%pixel(ipix)%ERR_lr(IW,ISD)  = ERR_temp
            GOUT_errest_aerosol%lidar%pixel(ipix)%BIAS_lr(IW,ISD) = BIAS_temp
         case default
            write(tmp_message,'(a,i0,a)') &
            'Index of Optical Function (IOF) = ',IOF,' value is not valid'
            G_ERROR(trim(tmp_message))
         end select

      elseif(ISD .eq. 3) then
 
         select case(IOF)
         case(index_ext)
            GOUT_errest_aerosol%opt%pixel(ipix)%ERR_extt(IW)  = ERR_temp
            GOUT_errest_aerosol%opt%pixel(ipix)%BIAS_extt(IW) = BIAS_temp
         case(index_ssa) 
            GOUT_errest_aerosol%opt%pixel(ipix)%ERR_ssat(IW)  = ERR_temp
            GOUT_errest_aerosol%opt%pixel(ipix)%BIAS_ssat(IW) = BIAS_temp
         case(index_lr)  
            GOUT_errest_aerosol%lidar%pixel(ipix)%ERR_lrt(IW)  = ERR_temp
            GOUT_errest_aerosol%lidar%pixel(ipix)%BIAS_lrt(IW) = BIAS_temp
         case default
            write(tmp_message,'(a,i0,a)') &
            'Index of Optical Function (IOF) = ',IOF,' value is not valid'
            G_ERROR(trim(tmp_message))
         end select
      
      else

         write(tmp_message,'(a,i0,a)') 'ISD = ',ISD,' value is not valid'
         G_ERROR(trim(tmp_message))

      endif ! ISD .lt. 3
                        
      return
      end subroutine set_opt_err

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! subroutine forward_model_pixel_PHMX moved to forw_model.f90 file
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine PAR_ERR_estimates (iu_main_output,        & ! IN
                                    sparse_solver_name,    &     
                                    IP2,                   &
                                    UF,UFS,nnz_err,        & 
                                    QS,QIM,AGENP,ALS,      &
                                    RIN,INVSING,segment_meas, &
                                    errest_pix,            &                                 
                                    UFNZ,GOUT_errest_par,  & ! INOUT 
                                    solver_timer           &                                                                                               
                                   )

      use mod_par_inv, only : KPARS,KIMAGE,KPAR
      use mod_fisher_matrix_ccs    
      use mod_retr_general_output_derived_type 
      !use iso_c_binding ! c_intptr_t
      use mod_time_utils
      use mod_sdata_derived_type
      use mod_stop_report

      implicit none
!	-------------------------------------------------------------------------------------
! IN :
      type(retr_input_settings), intent(in)  ::  RIN
      type(segment_data),        intent(in)  ::  segment_meas
      integer,                   intent(in)  ::  iu_main_output
      character(*),              intent(in)  ::  sparse_solver_name
      integer,                   intent(in)  ::  IP2,INVSING
      integer,                   intent(in)  ::  nnz_err
      real,dimension(KPAR),      intent(in)  ::  QIM    
      real,                      intent(in)  ::  AGENP 
      real,dimension(KIMAGE),            intent(in)  ::  ALS
      real,dimension(KPARS,KIMAGE),      intent(in)  ::  QS
      real,dimension(KPARS,KPARS,KIMAGE),intent(in)  ::  UFS     
      real,dimension(KPAR,KPAR),         intent(in)  ::  UF
      logical, dimension(KIMAGE),        intent(in)  ::  errest_pix
!	-------------------------------------------------------------------------------------
! INOUT :
      real,                              intent(inout)  ::  solver_timer
      type(nonzero),                     intent(inout)  ::  UFNZ    
      type(output_segment_err_estim_par),intent(inout)  ::  GOUT_errest_par		
!	-------------------------------------------------------------------------------------
! LOCAL :
      logical                              ::  IPRI_additional_info
      integer                              ::  npixels, KNSING, KNSINGF
      integer                              ::  IMQ, npix
      !real*4                               ::  time_start, time_finish
      real, dimension(KPARS)               ::  FFS,QS1
      integer                              ::  nnz_pix
      real*8                               ::  TEMP
      real*8, dimension(:,:), allocatable  ::  ATEMP 

      real*8, dimension(KPAR)              ::  b	
      real,   dimension(KPAR)              ::  b1	

      logical                              ::  ltest
      integer                              ::  alloc_stat
      integer                              ::  i, i1, j, ii, ii1, KN, KNF, ipix
      real                                 ::  tiny
      character(len=12),dimension(KIMAGE)  ::  pdate, ptime
!	-------------------------------------------------------------------------------------
      !call CPU_TIME(time_start)
      npixels = segment_meas%npixels
      KNSING  = RIN%KNSING
      KNSINGF = RIN%KNSINGF
      KN      = KNSING*npixels
      KNF     = KNSINGF*npixels

      IPRI_additional_info   = RIN%IPRI_additional_info
      IMQ     = RIN%IMQ

      pdate(:) = ' '
      ptime(:) = ' '
      do j=1,npixels
         call convert_time_to_string(segment_meas%pixels(j)%t, "%F", pdate(j))
         call convert_time_to_string(segment_meas%pixels(j)%t, "%X", ptime(j))
      enddo 

      select case(INVSING)
      case(2)
! Multi pixel inversion

        ltest = .false.
        if(ltest) then
          allocate(ATEMP(KNF,KNF),stat=alloc_stat)
          if (alloc_stat /= 0) then
            write(tmp_message,'(a)') &
            'error while trying to allocate ATEMP'
            G_ERROR(trim(tmp_message))
          endif
          ATEMP(:,:) = 0.0
          do i=1,nnz_err
            ATEMP(UFNZ%row(i),UFNZ%col(i)) = UFNZ%val(i)
          enddo ! i
        endif ! ltest
!write(*,*) 'IP2=',IP2,'  KNF=',KNF,'  KNSINGF=',KNSINGF,'  npixels=',npixels

        if(IMQ .eq. 2) then ! SVD method
          do ipix=1,npixels
          do i=1,KNSINGF
            j = KNSINGF*(ipix-1)+i
            b(:) = 0.0d0
            b(j) = 1.0d0
            call ITERQ (  IMQ, KNF,             & ! IN
                          real(UF(1:KNF,1:KNF)),&
                          real(b(1:KNF)),       &
                          b1(1:KNF)             & ! OUT
                       )
            !call ITERQ ( !IMQ,KNF,                 & ! IN
                          !real(ATEMP(1:KNF,1:KNF)),&
                          !real(b(1:KNF)),          &
                          !b1(1:KNF)                & ! OUT
                        !)
            GOUT_errest_par%pixel(ipix)%ERRP(i)  = sqrt( b1(j)*AGENP )
            GOUT_errest_par%pixel(ipix)%BIASP(i) = QIM(j)

            if(ltest) then
              b(1:KNF) = b1(1:KNF)
              !write(*,*) 'col=',j
              !write(*,'(10e14.5)') (b(i1),i1=1,KNF)
              tiny = 1e-4
              do ii=1,KNF
                TEMP = 0.0d0
                do ii1=1,KNF
                  TEMP=TEMP+ATEMP(ii,ii1)*b(ii1)
                enddo ! ii1
                if(TEMP .gt. tiny) then
                  write(*,'(a,i5,a,e14.5)') 'i=',ii,'  TEMP=',TEMP
                endif ! TEMP .gt. tiny
              enddo ! ii
            endif ! ltest         
          enddo ! i
          enddo ! ipix	               
          if(ltest) then
              deallocate(ATEMP,stat=alloc_stat)
              if (alloc_stat /= 0) then
                write(tmp_message,'(a)') &
                'error while trying to deallocate ATEMP'
                G_ERROR(trim(tmp_message))
              endif
          endif ! ltest

        else if(IMQ .eq. 3 .or. IMQ .eq. 1) then

          if(IP2 .gt. 0) &
          call compress_zeros ( KNSINGF,npixels, & ! IN
                                nnz_err,         &
                                UFNZ             & ! INOUT
                              )
          !write(iu_main_output,*) 'in PAR_ERR_estimates, before sparse_solver', &
          !                        '  nnz=',nnz_err'  
          npix = npixels
          call sparse_solver (  iu_main_output,     & ! IN
                                KNF, nnz_err, UFNZ, &
                                b,                  & ! INOUT
                                solver_timer,       & ! OUT
                                npix, KNSINGF       & ! OPTIONAL
                             )           
          if ( stop_report%status ) return
          b1(1:KNF) = b(1:KNF)
          do ipix=1,npixels
          do i=1,KNSINGF
            j = KNSINGF*(ipix-1)+i
            ! write(*,*) 'after  solver_SuperLU: j=',j
            GOUT_errest_par%pixel(ipix)%ERRP(i)  = sqrt( b1(j)*AGENP )
            GOUT_errest_par%pixel(ipix)%BIASP(i) = QIM(j)
            !write(*,*) 'col=',j
            !write(*,'(10e14.5)') (b1(i1),i1=1,KNF)
          enddo ! i
          enddo ! ipix	               
        endif ! IMQ .eq. 2

      case(0,1)
! Pixel by pixel inversion

        do ipix=1,npixels
        if(.not. errest_pix(ipix)) cycle
          if(IMQ .eq. 2) then ! SVD method
            do i=1,KNSINGF
              FFS(:) = 0.0
              FFS(i) = 1.0
              call ITERQ (IMQ,KNSINGF,               & ! IN
                          UFS(1:KNSINGF,1:KNSINGF,ipix), &
                          FFS(1:KNSINGF),            &
                          QS1(1:KNSINGF)             & ! OUT
                        )
              GOUT_errest_par%pixel(ipix)%ERRP(i)  = sqrt( QS1(i)*ALS(ipix) )
              GOUT_errest_par%pixel(ipix)%BIASP(i) = QS(i,ipix)
! OD 2015-10-09          GOUT_errest_par%pixel(ipix)%BIASP(i) = sqrt(QS(i,ipix)*QS(i,ipix)*ALS(ipix))
            enddo ! i
          else if(IMQ .eq. 3) then ! sparse matrix solver
            call UF_nonzero_one_pixel (KNSINGF, UFS(:,:,ipix), UFNZ,nnz_pix)
            if ( stop_report%status ) return
            !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
            !'UF inversion, before sparse_solver ',trim(sparse_solver_name),'  nnz_pix=',nnz_pix
            npix = 1
            call sparse_solver (  iu_main_output,         & ! IN
                                  KNSINGF, nnz_pix, UFNZ, &
                                  b,                      & ! INOUT
                                  solver_timer,           &
                                  npix, KNSINGF           & ! OPTIONAL
                                 )
            if ( stop_report%status ) return
              !if(RIN%IPRI_verbose)  write(iu_main_output,*)  &
              !'UF inversion, after  sparse_solver ',trim(sparse_solver_name)         
            QS1(1:KNSINGF) = b(1:KNSINGF)
            do i=1,KNSINGF
              GOUT_errest_par%pixel(ipix)%ERRP(i)  = sqrt( QS1(i)*ALS(ipix) )
              GOUT_errest_par%pixel(ipix)%BIASP(i) = QS(i,ipix)
            enddo ! i
          endif ! IMQ .eq. 2
        enddo ! ipix

      end select ! INVSING

      !call CPU_TIME(time_finish)

      !write(iu_main_output,*) 'PAR_ERR_estimates : CPU_TIME(sec)=',time_finish-time_start
      
      return
      end subroutine PAR_ERR_estimates
            
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss



