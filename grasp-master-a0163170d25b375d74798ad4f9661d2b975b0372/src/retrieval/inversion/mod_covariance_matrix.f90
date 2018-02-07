! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! subroutine cmatrix
! subroutine covariance_matrix_segment
! subroutine covariance_matrix_pixel
module mod_covariance_matrix

contains
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine cmatrix (  KL,INN,DNN,     & ! IN 
                            IFCOV,CMTRX,FS, & 
                            CSN             & ! OUT	
                         )
      implicit none
! -----------------------------------------------------------------
! IN :
      integer, intent(in) :: KL,INN,IFCOV	
      real,    intent(in) :: DNN,CMTRX,FS
! -----------------------------------------------------------------
! OUT :
      real,    intent(out) :: CSN	  	  
! -----------------------------------------------------------------
! -----------------------------------------------------------------
! Accounting for different accuracy levels in the data     
! INN - EQ.1.THEN error is assumed absolute   
!     - EQ.0 THEN error is assumed relative    
! DNN - variation of the noise of the I-th source  
! -----------------------------------------------------------------
            
! Covariance matrix element

      CSN = DNN*DNN

      SELECT CASE(KL)
      CASE(1)
         IF(INN .eq. 1) CSN = CSN/(FS*FS)        
         !write(*,*) DNN,INN,FS,'  - DNN,INN,FS'
      CASE(0)
         IF(INN .eq. 0) CSN = CSN*FS*FS
      END SELECT ! KL

      if(IFCOV .eq. 1) then
         CSN = CSN*CMTRX
         !write(*,*) 'in cmatrix: CSN=',CSN,'  CMTRX=',CMTRX
         !stop
      endif ! IFCOV .eq. 1
         
   return
   end subroutine cmatrix

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss 

! Accounting for different accuracy levels in the data     
! INOISE  - the number of different noise sources              
! SGMS(I) - std of noise in i -th source                      
! INN(I)  - EQ.1.THEN error is absolute with                  
!         - EQ.0 THEN error assumed relative
! DNN(I)  - variation of the noise of the I-th source
! IKI(I)   - total number of measurements of i-th source       
! IKI_shift(I) - presence of shift for measurements of i-th source       
! KNOISE(1,K) - specific numbers affected by i-th source        
! +
! Covariance matrix 

      subroutine covariance_matrix_segment (RIN,INVSING,           & ! IN
                                            segment_meas,          & 
                                            segment_vec_meas,      & ! INOUT 
                                            MNOISEI,               & ! IN
                                            IKI,IKI_shift,KNOISEI, & ! OUT
                                            CS,ARES2,              &
                                            GOUT_segment_residual  &       
                                           ) 	  
      use mod_sdata_derived_type
      use mod_par_inv, only      : KMESS,KKNOISE,KIP,KWM,KIMAGE
      use mod_retr_settings_derived_type
      use mod_retr_general_output_derived_type
            
      implicit none
!	----------------------------------------------------------------------------
! IN :	
      type(retr_input_settings),             intent(in) ::  RIN
      integer,                               intent(in) ::  INVSING
      type(segment_data),                    intent(in) ::  segment_meas
      integer,dimension(KIP,KWM,KIMAGE),     intent(in) ::  MNOISEI   
      type(pixel_vector),dimension(KIMAGE),  intent(in) ::  segment_vec_meas
      type(output_segment_residual),optional,intent(in) ::  GOUT_segment_residual   
!	----------------------------------------------------------------------------
! OUT :
      integer,dimension(KKNOISE,KIMAGE),      intent(out) ::  IKI       ! number of meas for i-th noise 
      logical,dimension(KKNOISE,KIMAGE),      intent(out) ::  IKI_shift ! presence of SHIFT for i-th noise 
      integer,dimension(KMESS,KKNOISE,KIMAGE),intent(out) ::  KNOISEI
      real,dimension(KMESS,KIMAGE),           intent(out) ::  CS
      real,dimension(KIMAGE),                 intent(out) ::  ARES2
!	-----------------------------------------------------------------------------
! LOCAL :
      real                                 ::  CS1
      integer                              ::  ipix,npixels,KMIMAGE,JJS
      type(pixel_vector)                   ::  segment_vec_meas_temp
!	------------------------------------------------------------------------------------------
      npixels = segment_meas%npixels

      ARES2(:) = 0.0

      do ipix=1,npixels
! Calculate covariance matrix
        JJS = SUM(segment_vec_meas(ipix)%nFS(1:segment_meas%pixels(ipix)%nwl))      
        segment_vec_meas_temp = segment_vec_meas(ipix)
        if(RIN%KL .eq. 1) &
        segment_vec_meas_temp%FS(1:JJS) = EXP(segment_vec_meas(ipix)%FS(1:JJS))
        if(present(GOUT_segment_residual)) then
        call covariance_matrix_pixel (RIN,ipix,INVSING,          & ! IN
                                      segment_meas%pixels(ipix), & 
                                      segment_vec_meas_temp,     &
                                      MNOISEI(:,:,ipix),         &
                                      IKI(:,ipix),               & ! OUT 
                                      IKI_shift(:,ipix),         &
                                      KNOISEI(:,:,ipix),         &
                                      CS(:,ipix),CS1,            &
                                      ARES2(ipix),               &
                                      GOUT_segment_residual%pixel(ipix)%resa(:),  &
                                      GOUT_segment_residual%pixel(ipix)%resr(:)   &
                                     )
        else
          call covariance_matrix_pixel (RIN,ipix,INVSING,          & ! IN
                                        segment_meas%pixels(ipix), & 
                                        segment_vec_meas_temp,     &
                                        MNOISEI(:,:,ipix),         &
                                        IKI(:,ipix),               & ! OUT 
                                        IKI_shift(:,ipix),         &
                                        KNOISEI(:,:,ipix),         &
                                        CS(:,ipix),CS1,            &
                                        ARES2(ipix)                &
                                       )
        endif
      enddo ! ipix

      if(INVSING .ne. 0) then ! multi pixel inversion
        do ipix=1,npixels
          KMIMAGE = segment_vec_meas(ipix)%KMIMAGE      
          CS(1:KMIMAGE,ipix) = CS1/CS(1:KMIMAGE,ipix)
        enddo ! ipix
! ARES - parameters adjusting the weight of a priori constraints during iterations
! if RESIDIAL > ARES, the weight of  a priori constraints is increased, once  
! RESIDIAL < ARES the weight of  a priori constraints is fixed
        ARES2(1:npixels) = ARES2(npixels)
      endif ! INVSING .ne. 0
      
      !if(.not. RIN%IPRI_silent_mode .and. RIN%IPRI2) then
      !   do ipix=1,npixels
      !      KMIMAGE=segment_vec_meas(ipix)%KMIMAGE            
      !      write(*,*)
      !      write(*,*) 'CS:  ipix=',ipix,'  CS1=',CS1,'  ARES2=',ARES2(ipix) 
      !      write(*,'(10e16.5)') CS(1:KMIMAGE,ipix)
      !   enddo ! ipix
      !   !write(0,*) 'in inversion: ipix=',ipix,'  IKI: ',IKI(1:RIN%NOISE%INOISE,ipix) 
      !endif ! .not. RIN%IPRI_silent_mode
  
      return
      end subroutine covariance_matrix_segment

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine covariance_matrix_pixel (RIN,ipix,INVSING,       & ! IN
                                          pixel_cont,             & 
                                          pixel_vec,              &
                                          MNOISEI,                &
                                          IKI,IKI_shift,KNOISEI,  & ! OUT 
                                          CS,CS1,                 &
                                          ARES2,                  &
                                          resa,resr               &
                                         )
      use mod_sdata
      use mod_par_inv, only      : KMESS,KKNOISE,KIP,KWM   									  
      use mod_retr_settings_derived_type
      
      implicit none
!	----------------------------------------------------------------------------
! IN :	
      type(retr_input_settings), intent(in) :: RIN
      integer,                   intent(in) :: ipix,INVSING
      type(pixel),               intent(in) :: pixel_cont
      type(pixel_vector),        intent(in) :: pixel_vec
      integer,dimension(KIP,KWM),intent(in) :: MNOISEI
      real,dimension(KKNOISE),optional,intent(in) :: resa,resr 
!	----------------------------------------------------------------------------
! OUT :
      integer,dimension(KKNOISE),      intent(out)   :: IKI       ! number of meas for i-th noise 
      logical,dimension(KKNOISE),      intent(out)   :: IKI_shift ! presence of SHIFT for i-th noise 
      integer,dimension(KMESS,KKNOISE),intent(out)   :: KNOISEI
      real,dimension(KMESS),           intent(out)   :: CS
      real,                            intent(inout) :: CS1
      real,                            intent(inout) :: ARES2		
!	-----------------------------------------------------------------------------
! LOCAL :
      integer                  ::  NWL, NIP, meas_type, iMN
      integer                  ::  JJS, IW, IV, IP, JJS1
      real                     ::  meas, I, Q2, U2, XIKN 
      real,save                ::  XIK1 
      integer                  ::  IFCOV,iw_1,IP_1,JJS_1
      real                     ::  CMTRX
      integer                  ::  KMIMAGE,iMN_min
!      real,dimension(KKNOISE)  ::  DNN
      real                     ::  DNN
!	-----------------------------------------------------------------------------
! iPOBS = 
!           1    I, Q, U
!           2    I, Q/I, U/I                  
!           3    I, P    or sqrt(Q*Q+U*U)
!           4    I, P/I  or sqrt(Q*Q+U*U)/I
!           5    I, P/I meas
! Types of measurements :
! 1 - tod(wl)        - total optical depth
! 2 - aod(wl)        - aerosol optical depth
! 3 - p11(angle,wl)  - phase matrix element p11
! 4 - p12(angle,wl)  - phase matrix element p12
! 5 - p22(angle,wl)  - phase matrix element p22
! 6 - p33(angle,wl)  - phase matrix element p33
! 7 - p34(angle,wl)  - phase matrix element p34
! 8 - p44(angle,wl)  - phase matrix element p44

! 9  - LS(height,wl)  - lidar signal
! 10  - DP(height,wl)  - depolarization ratio
! 11 - RL(height,wl)  - Raman lidar signal

! 12 - I(angle,wl)    - Stokes parameter I
! 13 - Q(angle,wl)    - Stokes parameter Q
! 14 - U(angle,wl)    - Stokes parameter U
! 15 - P(angle,wl)    - polarization sqrt(Q*Q+U*U)

! NOISE
! MNOISEI(1:NIP,iw)
!	-----------------------------------------------------------------------------      
        IKI(:)       = 0
        IKI_shift(:) = .false.
        CS(:)        = 0.0
        KNOISEI(:,:) = 0
        KMIMAGE      = 0

        if(INVSING .eq. 0 .or. (INVSING .gt. 0 .and. ipix .eq. 1)) then
! Search for index (iMN) of first measerement with first noise (MNOISEI(ip,iw)=1) 
          JJS_1 = 0
          iw_1  = 0
          ip_1  = 0
          if(INVSING .eq. 0) then ! single pixel inversion
            iMN=9999
            do iw=1,pixel_cont%NWL
            do ip=1,pixel_cont%meas(iw)%NIP
               if(MNOISEI(ip,iw) .lt. iMN) &
               iMN=MNOISEI(ip,iw)
            enddo ! ip
            enddo ! iw
          else                    ! multi  pixel inversion
            iMN=1
          endif ! INVSING .eq. 0

          loop_iw: do iw=1,pixel_cont%NWL
          do ip=1,pixel_cont%meas(iw)%NIP
            if(MNOISEI(IP,iw) .gt. iMN) cycle
            iw_1 = iw
            ip_1 = ip
!write(*,*) 'iw=',iw,'  IP=',IP,'  iw_1=',iw_1,'  IP_1=',IP_1,'  iMN=',iMN,'  in covariance_matrix_pixel'
            exit loop_iw
          enddo ! ip
          enddo loop_iw
!write(*,*) "Search for index (iMN): ipix=",ipix
        endif ! INVSING .eq. 0 .or.
! Define IKI and KNOISE 
      JJS = 0

      NWL = pixel_cont%NWL
LOOP_WL : do iw=1,NWL
      NIP  = pixel_cont%meas(iw)%NIP
LOOP_meas_type : do IP=1,NIP
         IFCOV     = pixel_cont%meas(iw)%IFCOV(ip)
         meas_type = pixel_cont%meas(iw)%meas_type(IP)
         if(RIN%iPOBS .ge. 3 .and. meas_type .eq. meas_type_U) &
         exit LOOP_meas_type
         iMN = MNOISEI(IP,iw)  
         if(iMN .lt. 1 .or. iMN .gt. RIN%NOISE%INOISE) then
            write(*,*) 'iMN=',iMN,'  INOISE=',RIN%NOISE%INOISE
            write(*,*) 'Problem with assigned noise indices ( 1<iMN<INOISE )'
            write(*,*) 'Please check noise definition in settings file. This definition has to match with data to be processed'
            stop 'stop in covariance_matrix_pixel'         
         endif !
         if(                                   &
            (meas_type .ge. meas_type_Q   .and. meas_type .le. meas_type_U .and. RIN%iPOBS .le. 2) .or.  &
            (meas_type .eq. meas_type_p12) .or.                                   &
            (meas_type .ge. meas_type_p33 .and. meas_type .le. meas_type_p44)     &
           ) IKI_shift(iMN) = .true.

         do IV=1,pixel_cont%meas(iw)%NBVM(IP)
         JJS=JJS+1              
         !iMN = MNOISEI(IP,iw)
         !if(iMN .lt. 1 .or. iMN .gt. RIN%NOISE%INOISE) then
            !write(*,*) 'iMN=',iMN,'  INOISE=',RIN%NOISE%INOISE
            !write(*,*) 'Problem with assigned noise indices ( 1<iMN<INOISE )'
            !write(*,*) 'Please check noise definition in settings file. This definition has to match with data to be processed'
            !stop 'stop in covariance_matrix_pixel'
         !endif !
         IKI(iMN) = IKI(iMN)+1
         KNOISEI(IKI(iMN),iMN) = JJS
         !if(                                   &
         !   (meas_type .ge. meas_type_Q   .and. meas_type .le. meas_type_U)   .or.  &
         !   (meas_type .eq. meas_type_p12) .or.                                     &
         !   (meas_type .ge. meas_type_p33 .and. meas_type .le. meas_type_p44)       &
         !  ) IKI_shift(iMN) = .true.

         if(meas_type .eq. meas_type_Q .and. RIN%iPOBS .gt. 2) then
         CMTRX = pixel_cont%meas(iw)%CMTRX(iv,ip  ) +  &
                 pixel_cont%meas(iw)%CMTRX(iv,ip+1)      
         else
         CMTRX = pixel_cont%meas(iw)%CMTRX(iv,ip)
         endif ! meas_type .eq. meas_type_Q .and.

         if( present(resa) .and. present(resr) ) then
          if(RIN%NOISE%INN(iMN) .eq. 1) then
            DNN = resa(iMN)
!            write(*,*) 'in cmatrix: ipix=',ipix,'  DNN=',DNN,'  INN(iMN)=',RIN%NOISE%INN(iMN)
          elseif(RIN%NOISE%INN(iMN) .eq. 0) then
            DNN = resr(iMN)
!            write(*,*) 'in cmatrix: ipix=',ipix,'  DNN=',DNN,'  INN(iMN)=',RIN%NOISE%INN(iMN)
          endif ! RIN%NOISE%INN(iMN) .eq. 1
         else
            DNN = RIN%NOISE%DNN(iMN)
         endif ! present
         call cmatrix ( RIN%KL,             & ! IN
                        RIN%NOISE%INN(iMN), & 
                        DNN,                &  
                        IFCOV,CMTRX,        & 
                        pixel_vec%FS(JJS),  & 
                        CS(JJS)             & ! OUT	
                      )
         if (CS(JJS) .EQ. 0) then
            write(*,*) 'zero in diagonal element of covariance matrix #', JJS
         endif

         if(INVSING .eq. 0 .or. (INVSING .gt. 0 .and. ipix .eq. 1)) then
            if(iw_1 .eq. iw .and. IP_1 .eq. IP .and. IV .eq. 1) JJS_1 = JJS
         endif ! INVSING .eq. 0 .or.
         enddo ! IV 
enddo LOOP_meas_type
enddo LOOP_WL

      IF(INVSING .eq. 0 .or. (INVSING .gt. 0 .and. ipix .eq. 1)) THEN
        iMN = MNOISEI(IP_1,iw_1)  
        XIK1 = IKI(iMN)
        if( present(resa) .and. present(resr) ) then
          if(RIN%NOISE%INN(iMN) .eq. 1) then
            DNN = resa(iMN)
          elseif(RIN%NOISE%INN(iMN) .eq. 0) then
            DNN = resr(iMN)
          endif ! RIN%NOISE%INN(iMN) .eq. 1
        else
          DNN = RIN%NOISE%DNN(iMN)
        endif ! present

        CS1 = DNN*DNN*XIK1
        IF(pixel_cont%meas(iw_1)%IFCOV(IP_1) .eq. 1) CS1 = CS1*pixel_cont%meas(iw_1)%CMTRX(1,IP_1) 
        SELECT CASE(RIN%KL)
        CASE(1)
          IF(RIN%NOISE%INN(iMN) .eq. 1) CS1 = CS1/(pixel_vec%FS(JJS_1)*pixel_vec%FS(JJS_1))
        CASE(0)
          IF(RIN%NOISE%INN(iMN) .eq. 0) CS1 = CS1*pixel_vec%FS(JJS_1)*pixel_vec%FS(JJS_1)
        END SELECT ! KL
! write(*,*) '***  ',iw_1,IP_1,JJS_1,iMN,XIK1,CS1,pixel_vec%FS(JJS_1), &
! '  -  iw_1,IP_1,JJS_1,iMN,XIK1,CS1,pixel_vec%FS(JJS_1)'
      ENDIF ! INVSING .eq. 0 .or.
      ARES2 = CS1/XIK1 
                 
      do iMN=1,RIN%NOISE%INOISE
        XIKN = IKI(iMN)       
        do iv=1,IKI(iMN)
          CS(KNOISEI(iv,iMN)) = CS(KNOISEI(iv,iMN))*XIKN
          if(INVSING .eq. 0) CS(KNOISEI(iv,iMN)) = CS1/CS(KNOISEI(iv,iMN))         
! write(*,*) iMN,XIK1,XIKN,CS(KNOISEI(iv,iMN)),'  - iMN,XIK1,XIK,CS'

        enddo ! iv 
      enddo ! iMN
                          
      return
      end subroutine covariance_matrix_pixel

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      real function delta_derivatives ( ALS,APMIN,APMAX,AP,DL )
                 
      implicit none
! ----------------------------------------------------------------	  
      real,intent(in)    ::  ALS,APMAX,APMIN,AP,DL
! ----------------------------------------------------------------	  
      real  ::  AFOR,ALS1,DELTA,DL1
! ----------------------------------------------------------------	  
      DELTA = APMAX-APMIN
      ALS1 = SQRT(ALS)
      if(ALS1 .gt. 100.) ALS1 = 100.       
      if( DELTA .gt. 0.5 ) then
        AFOR = 10.0*DELTA
      else
        AFOR = 1.0
      endif      
      DL1 = DL*DELTA*AFOR*ALS1*0.5
      if( DL1 .ge. DELTA ) DL1 = DELTA/3.0
      if(DL1 .lt. DL) DL1 = DL

      if( AP+DL1 .gt. APMIN .and. AP+DL1 .le. APMAX ) then
        delta_derivatives = DL1
        return
      elseif( AP+DL1 .gt. APMAX ) then
        if( AP-DL1 .ge. APMIN ) then
          DL1 = -DL1
        else 
          if( APMAX-AP .ge. AP-APMIN) then
            DL1 = APMAX-AP
          else
            DL1 = APMIN-AP
          endif
        endif
      elseif( AP+DL1 .le. APMIN ) then
        write(*,*) 'AP+DL1=',AP+DL1,' MIN=',APMIN
        stop 'stop in delta_derivatives: AP+DL1 .le. APMIN'
      endif
      delta_derivatives = DL1

      end function delta_derivatives 

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine set_surf_derivative_mask ( RIN, segment, dermask )

      use mod_par_inv, only : KPARS, KIMAGE
      use mod_retr_settings_derived_type
      use mod_sdata_derived_type

      implicit none
! -----------------------------------------------------------------
      type(retr_input_settings), intent(in) :: RIN
      type(segment_data), intent(in) :: segment
      logical, dimension(KPARS,KIMAGE), intent(out) :: dermask
! -----------------------------------------------------------------
      integer :: meas_type, par_type, npixels
      integer :: iw, ip, ipix, ibeg, iend
! lland  - presence of retrieved land surface parameters
! lwater - presence of retrieved water surface parameters
      logical :: lland = .false., lwater = .false.
! IDIM1 - chatacteristic number
      integer :: IDIM1, IDIM1_land1, IDIM1_land2, IDIM1_water
! npar - number of surface parameters
      integer :: npar_land1, npar_land2, npar_water
! -----------------------------------------------------------------
      npixels = segment%npixels

      do ipix=1,npixels
        dermask(1:RIN%KNSING,ipix) = .true.
      enddo
      if(RIN%KNSING .gt. RIN%KNSINGF) then
        do ipix=1,npixels
          dermask(RIN%KNSINGF+1:RIN%KNSING,ipix) = .false.
        enddo
      endif

      IDIM1_land1 = 0
      IDIM1_land2 = 0
      IDIM1_water = 0

      do IDIM1=1,RIN%NDIM%n1
        par_type = RIN%NDIM%par_type(IDIM1)
        if(par_type .gt. par_type_SURF1_land_beg .and. par_type .lt. par_type_SURF1_land_end) then
          if(RIN%NDIM%par_retr(IDIM1)) then
            npar_land1 = sum(RIN%NDIM%n3(1:RIN%NDIM%n2(IDIM1),IDIM1))
            IDIM1_land1 = IDIM1
            lland = .true.
          endif
        elseif(par_type .gt. par_type_SURF2_land_beg .and. par_type .lt. par_type_SURF2_land_end) then
          if(RIN%NDIM%par_retr(IDIM1)) then
            npar_land2 = sum(RIN%NDIM%n3(1:RIN%NDIM%n2(IDIM1),IDIM1))
            IDIM1_land2 = IDIM1
            lland = .true.
          endif
        elseif(par_type .gt. par_type_SURF_water_beg .and. par_type .lt. par_type_SURF_water_end) then
          if(RIN%NDIM%par_retr(IDIM1)) then
            npar_water = sum(RIN%NDIM%n3(1:RIN%NDIM%n2(IDIM1),IDIM1))
            IDIM1_water = IDIM1
            lwater = .true.
          endif
        endif
      enddo ! IDIM1

      if(lland .eqv. .true.) then
        do ipix=1,npixels
          if(NINT(segment%pixels(ipix)%land_percent) .eq. 0) then
            if(IDIM1_land1 .gt. 0) then
              ibeg = RIN%NDIM%ISTARSING(1,IDIM1_land1)
              iend = ibeg + npar_land1 - 1
              dermask(ibeg:iend,ipix) = .false.
            endif
            if(IDIM1_land2 .gt. 0) then
              ibeg = RIN%NDIM%ISTARSING(1,IDIM1_land2)
              iend = ibeg + npar_land2 - 1
              dermask(ibeg:iend,ipix) = .false.
            endif
          endif
        enddo ! ipix
      endif ! lland .eqv. .true.

      if(lwater .eqv. .true.) then
        do ipix=1,npixels
          if(NINT(segment%pixels(ipix)%land_percent) .eq. 100) then
            ibeg = RIN%NDIM%ISTARSING(1,IDIM1_water)
            iend = ibeg + npar_water - 1
            dermask(ibeg:iend,ipix) = .false.
          endif
        enddo ! ipix
      endif ! lwater .eqv. .true.

      return
      end subroutine set_surf_derivative_mask

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
end module mod_covariance_matrix
