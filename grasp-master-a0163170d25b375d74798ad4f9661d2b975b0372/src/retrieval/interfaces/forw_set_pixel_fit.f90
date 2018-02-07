! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../constants_set/mod_globals.inc"
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine set_pixel_phase_matr_fit (                       &
                                             NSD,tau_mol,         & 
                                             NANG,ANGL,           &
                                             GOUT_aerosol_opt_pixel_wl,  &
                                             GOUT_aerosol_phmx_pixel_wl, & ! INOUT									
                                             meas_type,IW,IP,iPOBS,      &
                                             pixel_fit            &
                                           )
      
      use mod_par_OS,   only : NBVM
      use mod_inversion_utils
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_sdata, only : get_pixel_geom, set_pixel_meas
      use mod_intrpl_linear
      use mod_retr_general_output_derived_type
      use mod_stop_report

      implicit none
!	------------------------------------------------------------------------------------------------------  
! IN:
      real,                           intent(in)  ::  tau_mol                            
      integer,                        intent(in)  ::  NSD,NANG,meas_type,IW,IP,iPOBS
      real,dimension(NANG),           intent(in)  ::  ANGL
      type(output_pixel_opt_wl),      intent(in)  ::  GOUT_aerosol_opt_pixel_wl
      type(output_pixel_ph_matrix_wl),intent(in)  ::  GOUT_aerosol_phmx_pixel_wl
!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(pixel),                    intent(inout)  ::  pixel_fit
!	------------------------------------------------------------------------------------------------------
      integer                ::  nsca_ang, ISD, iv
      real                   ::  sza, p11
      real, dimension(NBVM)  ::  vis, fiv, scat_angl, temp
      real, dimension(NBVM)  ::  meas
      real, dimension(NSD)   ::  sca
      logical                ::  status_funct
!	------------------------------------------------------------------------------------------------------
         vis(:) = 0.0
         call get_pixel_geom (IW,IP,		    &
                              pixel_fit,    &
                              nsca_ang,sza,	&
                              vis(:),fiv(:) &
                             )

         if(nsca_ang .gt. NBVM) then
            write(tmp_message,'(2(a,i0))') &
            'nsca_ang = ',nsca_ang,' .gt. NBVM = ',NBVM
            G_ERROR(trim(tmp_message))
         endif

         do iv=1,nsca_ang
            status_funct = geom2scat_angl(sza,vis(iv),fiv(iv),scat_angl(iv))
            !write(*,*) iv,NBVM,scat_angl(iv),'  iv,NBVM,scat_angl(iv)'
         enddo ! iv

         temp = 0.0
         meas(1:nsca_ang) = 0.0   
         select case(meas_type)
         case(meas_type_aod)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + GOUT_aerosol_opt_pixel_wl%EXT(ISD)
            enddo  ! ISD
            enddo ! iv
         case(meas_type_tod)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + GOUT_aerosol_opt_pixel_wl%EXT(ISD)
            enddo  ! ISD
            meas(iv) = meas(iv) + tau_mol
            if(pixel_fit%IFGAS .eq. 1) then
              meas(iv) = meas(iv) + pixel_fit%meas(IW)%gaspar
            endif ! IFGAS .eq. 1
            enddo ! iv
         case(meas_type_p11)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + sca(ISD)*LINEAR_LN(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph11(1:NANG,ISD),NANG,scat_angl(iv))
            enddo ! ISD
            enddo ! iv
         case(meas_type_p12)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            select case( iPOBS )
            case(1,2)
              do iv=1,nsca_ang
              do ISD=1,NSD
              meas(iv) = meas(iv) + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph12(1:NANG,ISD),NANG,scat_angl(iv))
              enddo ! ISD
              enddo ! iv
            case(5)
              do iv=1,nsca_ang
              p11 = 0.0
              do ISD=1,NSD
              meas(iv) = meas(iv) + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph12(1:NANG,ISD),NANG,scat_angl(iv))
              p11 = p11 + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph11(1:NANG,ISD),NANG,scat_angl(iv))
              enddo ! ISD
              meas(iv) = -meas(iv)/p11
              enddo ! iv
            end select
         case(meas_type_p22)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + sca(ISD)*LINEAR_LN(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph22(1:NANG,ISD),NANG,scat_angl(iv))
            enddo ! ISD
            enddo ! iv
         case(meas_type_p33)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph33(1:NANG,ISD),NANG,scat_angl(iv))
            enddo ! ISD
            enddo ! iv
         case(meas_type_p34)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph34(1:NANG,ISD),NANG,scat_angl(iv))
            enddo ! ISD
            enddo ! iv
         case(meas_type_p44)
            sca(1:NSD) = GOUT_aerosol_opt_pixel_wl%ssa(1:NSD)*GOUT_aerosol_opt_pixel_wl%ext(1:NSD)
            do iv=1,nsca_ang
            do ISD=1,NSD
            meas(iv) = meas(iv) + sca(ISD)*LINEAR(ANGL(1:NANG),GOUT_aerosol_phmx_pixel_wl%ph44(1:NANG,ISD),NANG,scat_angl(iv))
            enddo ! ISD
            enddo ! iv
         case default
            write(tmp_message,'(a,i0,a)') &
            'meas_type = ',meas_type,' - unknown value'
            G_ERROR(trim(tmp_message))
         end select

         call set_pixel_meas (                  &
                              IW,               & ! IN
                              nsca_ang,         &
                              meas_type,        &
                              1,                & ! 1 - set pixel with meas 
                              meas(1:nsca_ang),	& ! INOUT
                              pixel_fit         &
                             )

      return
      end subroutine set_pixel_phase_matr_fit
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_pixel_lidar_signal_fit (                  &
                                              nhvp,meas,       &
                                              meas_type,IW,IP, &
                                              pixel_fit        &
                                            )

      use mod_par_inv,   only :  KNBVM !,KVERTM
      use mod_sdata_derived_type
      use mod_sdata, only : set_pixel_meas
                  
      implicit none
!	------------------------------------------------------------------------------------------------------  
! IN:
      integer,               intent(in)     ::  meas_type,IW,IP,nhvp
      real,dimension(nhvp), intent(inout)  ::  meas
      !real,dimension(KNBVM), intent(inout)  ::  meas
      !real,dimension(KVERTM),intent(inout)  ::  meas

!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(pixel),           intent(inout)  ::  pixel_fit
!	------------------------------------------------------------------------------------------------------

        call set_pixel_meas (                &
                              IW,            & ! IN
                              nhvp,          &
                              meas_type,     &
                              1,             & ! 1 - set pixel with meas
                              meas(1:nhvp),	 & ! INOUT
                              pixel_fit      &
                             )

      return
      end subroutine set_pixel_lidar_signal_fit


! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine set_pixel_Stokes_vec_fit (  IW,IP,nmeas_type,iPOBS,  &
                                             NBV_comb,                &
                                             SLout_comb,SQout_comb,   &
                                             SUout_comb,SLPout_comb,  &
                                             pixel_fit                &
                                          )
                                             

      use mod_par_OS,  only : NBVM
      use mod_sdata_meas_type
      use mod_sdata_derived_type
      use mod_sdata, only : set_pixel_meas
            
      implicit none
!	------------------------------------------------------------------------------------------------------  
! IN:
      integer,                intent(in)  ::  IW,IP,NBV_comb,nmeas_type,iPOBS
      real,dimension(2*NBVM), intent(in)  ::  SLout_comb,SQout_comb,  &
                                              SUout_comb,SLPout_comb   

!	------------------------------------------------------------------------------------------------------
! INOUT:
      type(pixel),            intent(inout)  ::  pixel_fit
!	------------------------------------------------------------------------------------------------------
! LOCAL:
      real,dimension(NBVM)   ::  SQout,SUout,SLPout,SLout
      integer                ::  ind,NBV
!	------------------------------------------------------------------------------------------------------
         ind = 1 ! set pixel with meas 
         if(pixel_fit%meas(IW)%meas_type(IP) .eq. meas_type_I) &
         
            NBV = pixel_fit%meas(IW)%NBVM(IP)
            SLout(1:NBV) =  SLout_comb(1:NBV)       
            call set_pixel_meas (               & ! IN
                                 iw,            &
                                 NBV,   &
                                 meas_type_I,   &
                                 ind,           &
                                 SLout(1:NBV),	& ! INOUT
                                 pixel_fit      &
                                )
            
         if(nmeas_type-IP .gt. 0) then
            if(pixel_fit%meas(IW)%meas_type(IP+1) .eq. meas_type_Q) then
            NBV = pixel_fit%meas(IW)%NBVM(IP+1)
            SQout(1:NBV) =  SQout_comb(NBV_comb-NBV+1:NBV_comb)       
            call set_pixel_meas (               & ! IN
                                 iw,            &
                                 NBV,           &
                                 meas_type_Q,   &
                                 ind,           &
                                 SQout(1:NBV),	& ! INOUT
                                 pixel_fit      &
                                )
            SUout(1:NBV) =  SUout_comb(NBV_comb-NBV+1:NBV_comb)       
            call set_pixel_meas (               & ! IN
                                 iw,            &
                                 NBV,           &
                                 meas_type_U,   &
                                 ind,           &
                                 SUout(1:NBV),  & ! INOUT
                                 pixel_fit      &
                                )

            else if(pixel_fit%meas(IW)%meas_type(IP+1) .eq. meas_type_P) then
            NBV = pixel_fit%meas(IW)%NBVM(IP+1)
            select case(iPOBS)
            case(3,4)
              SLPout(1:NBV) =  abs(SLPout_comb(NBV_comb-NBV+1:NBV_comb))
            case(5)
              SLPout(1:NBV) =  abs(SLPout_comb(NBV_comb-NBV+1:NBV_comb))/  & 
                                    SLout_comb(NBV_comb-NBV+1:NBV_comb)
            case default
              write(*,*) 'RIN%iPOBS=',iPOBS,  & 
                ' value is not valid to fit [degree of] linear polarization'
              write(*,*) 'stop in set_pixel_Stokes_vec_fit'
              stop
            end select
            call set_pixel_meas (               & ! IN
                                 iw,            &
                                 NBV,           &
                                 meas_type_P,   &
                                 ind,           &
                                 SLPout(1:NBV),	& ! INOUT
                                 pixel_fit      &
                                )         
            endif ! pixel_fit%meas(IW)%meas_type(IP+1) .eq. meas_type_Q   
         endif ! nmeas_type-IP .gt. 0 

      return 
      end subroutine set_pixel_Stokes_vec_fit
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      
