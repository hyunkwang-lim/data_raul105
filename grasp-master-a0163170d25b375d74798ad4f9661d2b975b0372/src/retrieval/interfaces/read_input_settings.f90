! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **


! file contains :
! subroutine read_input_settings

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! with description by OD

! read Retrieval setting INput file (RIN structure)

      subroutine read_input_settings (inversion_input_settings_file, & ! IN
                                      main_output_file,            & 
                                      plotting_output_file,        &
                                      sdata_sim_file,              &
                                      RIN,iu_main_output           & ! OUT
                                     )

      use mod_par_inv
      use mod_retr_settings_derived_type
     
      implicit none

      character(*), intent(in)       ::  inversion_input_settings_file
      character(*), intent(in)       ::  main_output_file
      character(*), intent(in)       ::  plotting_output_file
      character(*), intent(in)       ::  sdata_sim_file
! --------------------------------------------------------------------------------------
      type(retr_input_settings),  intent(out) ::  RIN
      integer,                    intent(out) ::  iu_main_output
! --------------------------------------------------------------------------------------
! LOCAL :
      integer                         ::  I, II, IW, KSTAR, keyBIN,   &
                                          IDIM1, IDIM2, IDIM3,        &
                                          iu_iP_input
      character(255)                  ::  msg 
      character(255)                  ::  distname_O,distname_N
      real,dimension(KIDIM2)          ::  RMIN,RMAX
      integer                         ::  par_type
! ---------------------------------------------------------------------------------------

      if(KIDIM3 .gt. 22 .and. KVERTM .lt. KIDIM3) then
         write(*,*) 'KVERTM=',KVERTM,' <  KIDIM3=',KIDIM3
         stop 'stop in read_input_settings: Check parameters in mod_par_inv'
      endif ! KIDIM3 .gt. 22 .and.
	
      RIN%main_output_file     = main_output_file
      RIN%plotting_output_file = plotting_output_file
      RIN%sdata_sim_file       = sdata_sim_file
      
!write(*,*) 'inversion_input_settings_file=',inversion_input_settings_file
! iu_iP_input=5 read Retrieval INput settings from STDIN	  	  
      if (inversion_input_settings_file .ne. '-' ) then	  
         open (newunit=iu_iP_input, FILE=inversion_input_settings_file, STATUS='old')
      else
         iu_iP_input = 5      
      end if ! inversion_input_settings_file .ne. '-'
	  
!****************************************************************************************
!*** The parameters for DLS_bin code
!****************************************************************************************
! for DLS_bin code
! DLSF (DLS Flags)   
!   keyEL        - number of scattering matrix elements in inversion
!**************************************************************************************** 

! Check of parameters for arrays

      call retr_input_initialization ( RIN )
      
      read (iu_iP_input,*) RIN%DLSF%keyEL 

      if (RIN%main_output_file .eq. '-') then
         iu_main_output = 6
      else
! open main output file
         open (newunit=iu_main_output, FILE=trim(RIN%main_output_file), STATUS='replace',recl=20000)
      endif ! RIN%main_output_file .eq. '-'

!	----------------------------------------------------------

!****************************************************************************************
!***  The parameters that define GENERAL RETRIEVAL !!!
!****************************************************************************************
!*** KNSING  - total number of parameters driving forward model
!*** KNSINGF - number of retrieved pameters (i.e. KNSING-KNSINGF parameters are fixed 
!***           to the values given by initial guess
!*** KL    = 0 - minimization of absolute errors
!***       = 1 - minimization of log
!*** ISTOP = 0 - the full retrieval is implemented for the satellite data read from 
!***             "SDATA_NEW.dat".
!***       = 1 - only single forward run is implemented (if INPUT=1, the program simulates 
!***             the satellite data  using aerosol and surface propertes given in 
!***             "input_iguess.dat" file.
!***             These simulated satelite data are writted in "SDATA_NEW_S.dat" file.
!***       = 2 - the full retrieval is implemented for the simulated satellite data read 
!***             from "SDATA_NEW.dat". 
!*** IPRI_additional_info = 0 - no additional print
!***                      = 1 - additional print for debagging
!*** IPRI_output_retr = 0 print retr. results after last iteration 
!***                  = 1 print retr. results after each iteration
!****************************************************************************************

      read (iu_iP_input,*) RIN%KL,RIN%ISTOP,RIN%IPRI_additional_info,RIN%IPRI_iter_result

      !if(RIN%IPRI_additional_info .and. .not. RIN%IPRI_verbose) then 
      !  write(*,*) 'IPRI_additional_info can not be ',RIN%IPRI_additional_info,' for IPRI_silent_mode=',RIN%IPRI_verbose
      !  stop 'stop in read_input_settings'
      !endif

      !if(.not. RIN%IPRI_verbose) RIN%IPRI_additional_info = .false.
      
      if(RIN%IPRI_additional_info) then
         write(iu_main_output,*) 'KL, ISTOP, IPRI_additional_info :'
         write(iu_main_output,'(i5,2l9)') RIN%KL,RIN%ISTOP,RIN%IPRI_additional_info
      endif

!****************************************************************************************
!***   NSD   - number of aerosol component
!***   
!***   NLYRS - number of layers in radiance transfer calculations
!***   
!***   iPOBS = 0 - to fit Radiance I and Degree of polarization sqrt((Q*Q+U*U)/I)
!***         = 1 - to fit I, Q/I, U/I 
!***         = 2 - to fit I, Q, U
!***         = 3 - to fit I,sqrt(Q*Q+U*U)
!***   ipplane = 0 - calculated polarization is projected to plane of scattering
!***           = 1 - calculated polarization is projected to meridian plane 
!***   
!C**   SHIFT   - is value (usually 1) that is added to negative
!C**             parameters in order to be able use log transformations
!*** 
!*** parameters for calculating sun/sky radiances  (OS_HERMAN)
!***   iBRDF =  - to be described by PL
!***          
!***   iBPDF =  - to be described by PL
!***   
!***   isurf = 0 - OS_HERMAN w/o changes by PL
!***         = 1 - land  surface
!***         = 2 - ocean surface 
!***   ITRONC = .true.  - phase function truncation in radiative transfer
!***          = .false. - no phase function truncation in radiative transfer
!***   eps_err - threshold for convergence of radiance calculation in radiative transfer
!****************************************************************************************
      read (iu_iP_input,*) RIN%NSD,RIN%NLYRS,RIN%ipplane,RIN%SHIFT,  &
                           RIN%iPOBS,RIN%iIOBS,                      &
                           RIN%ITRONC,RIN%eps_err,RIN%mol_prof_type, RIN%aer_prof_type

      if(RIN%IPRI_additional_info) then
        write(iu_main_output,'(2a)')  &
        'NSD,NLYRS,ipplane,SHIFT,iPOBS,iIOBS,ITRONC,eps_er,mol_prof_type, aer_prof_type:'
        write(iu_main_output,'(3i5,e12.4,2i5,l5,e12.4,i5)')    &
                            RIN%NSD,RIN%NLYRS,RIN%ipplane,RIN%SHIFT,  &
                            RIN%iPOBS,RIN%iIOBS,                      &
                            RIN%ITRONC,RIN%eps_err, RIN%mol_prof_type, RIN%aer_prof_type
      endif
!C****************************************************************************************
!C***  The parameters define LIGHT SCATTERING REGIMEs 
!C***   of BOTH modeling and retrieving      !!!
!C****************************************************************************************
!  OSHF (OS_HERMAN Forward)     - for signal
!  OSHD (OS_HERMAN Derivatives) - for sim. matrix
!C***  
!C***  IMSC = 0 - multiple scattering regime for signal
!C***       = 1 - single scattering regime for signal
!C***  NG1 - number of terms in expansion of scattering matrix
!C***  
!C***  NG2 - number of terms in angle integration         
!C***  
!C***  NF  - to be defined by PL
!*****************************************************************************************

      read (iu_iP_input,*) RIN%OSHF%IMSC,RIN%OSHF%NG,RIN%OSHF%NN,RIN%OSHF%NF,   &
                           RIN%OSHD%IMSC,RIN%OSHD%NG,RIN%OSHD%NN,RIN%OSHD%NF 

      if(RIN%IPRI_additional_info) then
      write(iu_main_output,*)                       &
						  'OSHF%IMSC,OSHF%NG1,OSHF%NG2,OSHF%NF',   &
						  'OSHD%IMSC,OSHD%NG1,OSHD%NG2,OSHD%NF'

      write(iu_main_output,'(10i5)')                                        &
                           RIN%OSHF%IMSC,RIN%OSHF%NG,RIN%OSHF%NN,RIN%OSHF%NF,  &
                           RIN%OSHD%IMSC,RIN%OSHD%NG,RIN%OSHD%NN,RIN%OSHD%NF
      endif 

!C****************************************************************************************
!C***  The the atmospheric characteristics 
!C***  that are retrieved for EACH PIXES    
!C****************************************************************************************
!   NDIM1,( NDIM2(IDIM1),(NDIM3(IDIM2,IDIM1),IDIM2=1,NDIM2(IDIM1)),IDIM1=1,NDIM1 )
!
!  SD(nbin) normalized
!  RE part of Refractive index n(wl) / parameters driving Chemistry model
!  IM part of Refractive index k(wl)
!  Sphericity / Shape distribution
!  Aerosol Profile
!  Aerosol Concentrations 
!  BRDF  
!  BPDF
!  Lidar calibtation coefficient
!**** 
!  NDIM
!C*** NDIM1 - number of different retrieved characteristics
!C*** NDIM(KCHAR) - "dimension of each characteristic
!C*** NCHPAR(KCHAR,KDIM) - number of parameter describing 
!C***                      each characteristic 
!C*** IDIM1       - characteristic
!C**** NDIM2(IDIM1)       -  number of aerosol components 
!C***** NDIM3(IDIM2,IDIM1)  - number of parameters used for characteristic 
!C***      size distribution  (e.g. number of size bins) 
!C***
!C******   WL - number of parameters used for 
!C***     describing real part of ref. index
!C***     (e.g. usually it is number of wavelengths) or 4 chemistry parameters
!C***
!C******  WL - number of parameters used for 
!C***     describing imaginary part of ref. index
!C***     (e.g. usually it is number of wavelengths)
!C***
!C******  number of parameters used for 
!C***      shape distribution  (e.g. number of axis ratio bins)  
!C***
!C******  number of parameters used for height 
!C***
!C******  number of parameters used for aerosol concentration 
!C***
!C******  number of parameters used for 
!C***        describing each component of BRDF of the surface
!C***     (e.g. it can be number of wavelengths)
!C***     
!C******  number of parameters used for 
!C***        describing each component of BPDF of the surface
!C***     (e.g. it can be number of wavelengths)
!C****************************************************************************************

      read (iu_iP_input,*)   RIN%NDIM%n1,(RIN%NDIM%par_type(II),RIN%NDIM%par_retr(II),II=1,RIN%NDIM%n1)
      if(RIN%IPRI_additional_info) WRITE(iu_main_output,*)  &
                    RIN%NDIM%n1,(II,RIN%NDIM%par_type(II),RIN%NDIM%par_retr(II),II=1,RIN%NDIM%n1),  &
                    '   - NDIM1,II,par_type(II),par_retr(II)' 

      if(RIN%NDIM%n1 .gt. KIDIM1) then
        write(*,'(2(a,i4))') 'NDIM1=',RIN%NDIM%n1,'  KIDIM1=',KIDIM1
        write(*,*) 'KW in module mod_par_inv'
        stop 'stop in read_input_settings'
      endif

      do IDIM1=1,RIN%NDIM%n1        
        read (iu_iP_input,*) RIN%NDIM%n2(IDIM1)
        if(RIN%IPRI_additional_info) THEN
          write(iu_main_output,*) 'IDIM1,NDIM2(IDIM1):'
          write(iu_main_output,'(2i5)') IDIM1,RIN%NDIM%n2(IDIM1)
        endif
        read (iu_iP_input,*) (RIN%NDIM%n3(IDIM2,IDIM1),IDIM2=1,RIN%NDIM%n2(IDIM1))
        do IDIM2=1,RIN%NDIM%n2(IDIM1)
          if(RIN%NDIM%n3(IDIM2,IDIM1) .GT. KIDIM3) THEN
            write(*,'(4(a,i4))') 'IDIM1=',IDIM1,'  IDIM2=',IDIM2,'  KIDIM3=',KIDIM3,'  NDIM3=',RIN%NDIM%n3(IDIM2,IDIM1)
            write(*,*) 'KIDIM3 in module mod_par_inv'
            stop 'STOP in read_input_settings'
          endif
          if(RIN%IPRI_additional_info) then
            write(iu_main_output,*) 'IDIM1,IDIM2,NDIM3(IDIM2,IDIM1): ' 
            write(iu_main_output,'(3i5)') IDIM1,IDIM2,RIN%NDIM%n3(IDIM2,IDIM1) 
          endif
        enddo ! IDIM2  
      enddo ! IDIM1

! Set Size Distribution flags for spheroid package (DLS)
      call set_RIN_nsph_package_flags_SD ( RIN,iu_main_output )
      
!C****************************************************************************************
!C***  ISTARSING(IDIM1,IDIM2,IDIM3)
!C***  Definition of the places of the retrieved parameters
!C***  in the vector of unknown AP for single pixel
!C****************************************************************************************
      KSTAR=1
      do IDIM1=1,RIN%NDIM%n1
        do IDIM2=1,RIN%NDIM%n2(IDIM1) 
          RIN%NDIM%ISTARSING(IDIM2,IDIM1) = KSTAR
          KSTAR = KSTAR+RIN%NDIM%n3(IDIM2,IDIM1)
        enddo   
      enddo
!C****************************************************************************************
!C*** INOISE  - the number of different noise sources             
!C*** SGMS(I) - std of ADDED synthetic rundom noise to the data where i -th source  of 
!C***           noise is expected, if SGMS(I)=0 - no synthetic random noise is added                    
!C*** INN(I)  - EQ.1.THEN error is absolute with  (for both definiton of cov matrix used 
!C***           in fitting and in synthetic noise                
!C***         - EQ.0 THEN error assumed relative
!C*** DNN(I)  - variation of the noise of the I-th source (this defines the matrix used 
!C***           in fitting) 
!C****************************************************************************************
      read (iu_iP_input,*)   RIN%NOISE%INOISE		
      do I=1,RIN%NOISE%INOISE
        read (iu_iP_input,*)  RIN%NOISE%SGMS(I),RIN%NOISE%INN(I),RIN%NOISE%DNN(I) ,          &
                              RIN%NOISE%NMT(I),(RIN%NOISE%MT(II,I),RIN%NOISE%NWLP(II,I),     & 
                              RIN%NOISE%IWLP(1:RIN%NOISE%NWLP(II,I),II,I),II=1,RIN%NOISE%NMT(I))       
        if(RIN%IPRI_additional_info) then
          if(I .EQ. 1) write (iu_main_output,*) 'NOISE PARAMETERS:' 
          write (iu_main_output,'(2(i5,e15.5))')     &
                           I,RIN%NOISE%SGMS(I),RIN%NOISE%INN(I),RIN%NOISE%DNN(I) 
        endif ! RIN%IPRI_additional_info
      enddo ! I
!C****************************************************************************************
!C***  IBIN - index determining the binning of SD:
!C***         = -1 equal in logarithms
!C***         =  1 equal in absolute scale
!C****************************************************************************************
      read (iu_iP_input,*)   RIN%IBIN
      if(RIN%IPRI_additional_info) write(iu_main_output,'(a11,2i5)') 'IBIN: ',RIN%IBIN
!C****************************************************************************************
!C assigning specifications for retrieved characteristics
!C e.g. sizes, aspect ratios, etc.
!C**************************************************************************************** 
      do IDIM1=1,RIN%NDIM%n1 
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SD_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SD_end) then
          par_type = RIN%NDIM%par_type(IDIM1)
          if(par_type .eq. par_type_SD_TB) then
            do IDIM2=1,RIN%NDIM%n2(IDIM1)
              RIN%NBIN(IDIM2)=RIN%NDIM%n3(IDIM2,IDIM1)
              read(iu_iP_input,*) RIN%RMIN(IDIM2),RIN%RMAX(IDIM2)
              if(RIN%IPRI_additional_info)  write(iu_main_output,*)    &
              RMIN(IDIM2),RMAX(IDIM2),IDIM2,RIN%NDIM%n3(IDIM2,IDIM1),' RMIN(IDIM2),RMAX(IDIM2),IDIM2,NDIM3' 
            enddo ! IDIM2
          else 
            do IDIM2=1,RIN%NDIM%n2(IDIM1)
              RIN%NBIN(IDIM2)=RIN%NDIM%n3(IDIM2,IDIM1)
              read(iu_iP_input,*) RIN%RMIN(IDIM2),RIN%RMAX(IDIM2),  &
                                 (RIN%RADIUS1(IDIM3,IDIM2),IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1))  
              if(RIN%IPRI_additional_info) then
                do IDIM3=1,RIN%NDIM%n3(IDIM2,1)
                    write(iu_main_output,'(e12.4,i5,a8)') RIN%RADIUS1(IDIM3,IDIM2),IDIM3,'  RADIUS'
                enddo ! IDIM3
              endif  ! RIN%IPRI_additional_info
            enddo ! IDIM2
          endif ! par_type_SD .eq. par_type_SD_TB
          exit ! IDIM1 loop 
        endif
      enddo ! IDIM1
      do IDIM1=1,RIN%NDIM%n1
        if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SHD_beg .and. RIN%NDIM%par_type(IDIM1) .lt. par_type_SHD_end) then
          do IDIM2=1,RIN%NDIM%n2(IDIM1) 
            if(RIN%NDIM%par_type(IDIM1) .gt. par_type_SHD_distr) then
              read(*,*) (RIN%RATIO1(IDIM3,IDIM2),IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1))
            else
              RIN%RATIO1(1:RIN%NDIM%n3(IDIM2,IDIM1),IDIM2)=1.0
            endif
!           write(*,*) 'IDIM2 IDIM3  RATIO1'
!           write(*,*) IDIM2, IDIM3, RIN%RATIO1(1:RIN%NDIM%n3(IDIM2,IDIM1),IDIM2)
          enddo  ! IDIM2   
        exit ! IDIM1 loop 
        endif
      enddo  ! IDIM1
!C*****************************************************************************************
!C***  Parameters for ITERQ.f (details in "iterqP" subr.)
!C***  IMQ - parameter defining the inversion in Q-image iterations
!C***  IMQ=1 - linear iterations
!C***  IMQ=2 - Q-iterations based on SVD inversion of for each IMAGE (not supported)
!C***  IMQ=3 - sparce matrix solver
!C***  INVSING - determines if inversion is done for each image, or for all images at once
!C***  INVSING = 0  only PIXEL by PIXEL inversion (inversion is done for each pixel) 
!C***  INVSING = 1  first, the MULTI-PIXEL inversion (inversion is done for whole group of 
!C***               pixels at once) done, 
!C***          = 2  THEN once residual stops in MULTI-PIXEL fitting,  PIXEL by PIXEL 
!C***               inversion starts with initial guess equal the solution of 
!C***               MULTI-PIXEL fitting
!C*** 
!C***  TTEST - the parameter defining if retrieved parameter is variable in time or not:
!C***          IF(DT.EQ.0.OR.TTEST.LE.GT(i) /DT) then the parameter which has GT(i) is 
!C***    	     constant with time, i.e it is constant in different pixels with the same X and Y
!C***          coordinates but different T (GT(i) - Lagrange parameter for a priori restrictions 
!C***          of the smoothness with time)
!C***  XTEST - the parameter defining if retrieved parameter is variable along X 
!C***  YTEST - the parameter defining if retrieved parameter is variable along Y
!C***
!C*****************************************************************************************
      read (iu_iP_input,*)
      read (iu_iP_input,*)  RIN%IMQ,   &
         RIN%IPFP%INVSING,RIN%IPFP%TXY_group,RIN%IPFP%TTEST,RIN%IPFP%XTEST,RIN%IPFP%YTEST   
      if(RIN%IPRI_additional_info)  then
         write(iu_main_output,*) 'IMQ,INVSING,TXY_group,TTEST,XTEST,YTEST:'
         write(iu_main_output,'(2i5,l10,3e15.5,l10)') RIN%IMQ,   & 
         RIN%IPFP%INVSING,RIN%IPFP%TXY_group,RIN%IPFP%TTEST,RIN%IPFP%XTEST,RIN%IPFP%YTEST
      endif ! RIN%IPRI_additional_info
!C*****************************************************************************************
!C*** PARAMETERS FOR Levenb.-Marquardt correction:
!C i.e. the Lagrange parameters are adjusted with
!C  each p-th iterations:
!C*********
!C***  IPSTOP - the maximum number of IP iterations
!C***           where  Levenb.-Marquardt cor. applied
!C*** ...MIN - minimum 
!C*** ...MAX - maximum
!C***          - the variability limits of retrieved 
!C***    parameters they will be used for calculating
!C***    the strength of Lev.-Marq. correction  
!C***
!C*************************

      read(iu_iP_input,*)
      read(iu_iP_input,*)  RIN%IPSTOP

!C*****************************************************************************************
!C***       Set a priori constraints
!C*****************************************************************************************
!C***         A PRIORI Estimates (single pixel)
!C*****************************************************************************************
!C***  The corresponding Lagrange parameter is used as
!C***  the one assumed for a priori SD
!C***  A priori values for SD ends are taken from initial guess
!C*****************************************************************************************
!C  SPCA (Single Pixel Constraints, A priori) :
!C***  IO(NPAR) - order of correlations for a priori estimates
!C***  GSM(NPAR)- Lagrange multiplayer for a priori estimates
!C*****************************************************************************************
      if(RIN%IPRI_additional_info)  WRITE(iu_main_output,*) 'read single pixel a priori constraints'
      read(iu_iP_input,*)
      do IDIM1=1,RIN%NDIM%n1
        do IDIM2=1,RIN%NDIM%n2(IDIM1) 
          do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            read(iu_iP_input,*) RIN%SPCA%IO(IDIM3,IDIM2,IDIM1),RIN%SPCA%GSM(IDIM3,IDIM2,IDIM1)
            if(RIN%IPRI_additional_info) then
              write(iu_main_output,'(4i5,e15.5,2a)')  IDIM3,IDIM2,IDIM1, &
              RIN%SPCA%IO(IDIM3,IDIM2,IDIM1),RIN%SPCA%GSM(IDIM3,IDIM2,IDIM1), &
              '  IDIM3,IDIM2,IDIM1,',' IO0(IDIM3,IDIM2,IDIM1),GSM0(IDIM3,IDIM2,IDIM1)'
            endif
          enddo ! IDIM3
        enddo ! IDIM2
     enddo ! IDIM1

!C*****************************************************************************************
!C***         Smoothness constraints (single pixel)                 
!C*****************************************************************************************
!C SPCS - Single Pixel Constraints (Smoothness) :
!C***  IO(NPAR) - order of smoothness for I-th parameter depend.
!C***  GSM(NPAR)- Lagrange multiplayer for  smoothness 
!C***             for I-th parameter dependence
!C*****************************************************************************************
      if(RIN%IPRI_additional_info)  WRITE(iu_main_output,*) 'read single pixel smoothness constraints'
      read(iu_iP_input,*)
      do IDIM1=1,RIN%NDIM%n1
        do IDIM2=1,RIN%NDIM%n2(IDIM1)
            read(iu_iP_input,*) RIN%SPCS%IO(IDIM2,IDIM1),RIN%SPCS%GSM(IDIM2,IDIM1)
            if(RIN%IPRI_additional_info) WRITE(iu_main_output,'(3i5,e15.5,2a)') IDIM2,IDIM1, &
                          RIN%SPCS%IO(IDIM2,IDIM1),RIN%SPCS%GSM(IDIM2,IDIM1), &
                          '  IDIM2,IDIM1,',' IO(IDIM2,IDIM1),GSM(IDIM2,IDIM1)'
        enddo ! IDIM2
      enddo ! IDIM1

!C*****************************************************************************************
!C***         Smoothness constraints (multi pixel)                
!C*****************************************************************************************
!C MPCS - Multi Pixel Constraints (Smoothness) 
!C***  IOT(...) - order of smoothness of time dependence
!C***  GSMT(NPAR)- Lagrange multiplayer for  smoothness 
!C***             of time dependence
!C***  IOX(...) - order of smoothness of X dependence
!C***  GSMX(NPAR)- Lagrange multiplayer for  smoothness 
!C***             of X dependence
!C***  IOY(...) - order of smoothness of Y dependence
!C***  GSMY(NPAR)- Lagrange multiplayer for  smoothness 
!C***             of Y dependence
!C*****************************************************************************************
      if(RIN%IPRI_additional_info)  WRITE(iu_main_output,*) 'read multi pixel smoothness constraints'
      read(iu_iP_input,*)
      do IDIM1=1,RIN%NDIM%n1
        do IDIM2=1,RIN%NDIM%n2(IDIM1) 
          read(iu_iP_input,*) RIN%MPCS%IOT(IDIM2,IDIM1),RIN%MPCS%GSMT(IDIM2,IDIM1), &
                              RIN%MPCS%IOX(IDIM2,IDIM1),RIN%MPCS%GSMX(IDIM2,IDIM1), &
                              RIN%MPCS%IOY(IDIM2,IDIM1),RIN%MPCS%GSMY(IDIM2,IDIM1) 
          if(RIN%IPRI_additional_info) write(iu_main_output,'(2i5,3(i5,e15.5),a46)') IDIM1,IDIM2, & 
                        RIN%MPCS%IOT(IDIM2,IDIM1),RIN%MPCS%GSMT(IDIM2,IDIM1), &
                        RIN%MPCS%IOX(IDIM2,IDIM1),RIN%MPCS%GSMX(IDIM2,IDIM1), &
                        RIN%MPCS%IOY(IDIM2,IDIM1),RIN%MPCS%GSMY(IDIM2,IDIM1), &
                        '  IDIM2,IDIM1,IOT,GSMT,IOX,GSMX,IOY,GSMY' 
       enddo ! IDIM2  
      enddo ! IDIM1

!C*****************************************************************************************
!C***  EPSP - for stopping p - iterations
!C***  EPSQ - for stopping q - iterations in ITERQ
!C***  DL   - for calc. derivatives, (see FMATRIX)
!C*****************************************************************************************

      read (iu_iP_input,*) RIN%MAXP,RIN%EPSP,RIN%EPSQ,RIN%DL
      if(RIN%IPRI_additional_info) then
         write(iu_main_output,*) 'MAXP,EPSP,EPSQ,DL:'
         write(iu_main_output,'(i5,3e15.5)') RIN%MAXP,RIN%EPSP,RIN%EPSQ,RIN%DL
      endif

!C*****************************************************************************************
!C***  Read initial guess (it is the same for each SINGLE PIXEL)
!C***  for the solution (absolute values)
!C*****************************************************************************************

      if(RIN%IPRI_additional_info) WRITE(iu_main_output,*) 'read AP initial guess'
      read (iu_iP_input,*)
	  
      do IDIM1=1,RIN%NDIM%n1
        do IDIM2=1,RIN%NDIM%n2(IDIM1) 
          do IDIM3=1,RIN%NDIM%n3(IDIM2,IDIM1)
            II= RIN%NDIM%ISTARSING(IDIM2,IDIM1)+IDIM3-1
            read(iu_iP_input,*)                            &
            RIN%APSING(II),RIN%APSMIN(II),RIN%APSMAX(II),RIN%IWW_SINGL(II)
            if(RIN%IPRI_additional_info) then 
              write(iu_main_output,'(4i5,3e15.5,a28)') II,IDIM1,IDIM2,IDIM3, &
              RIN%APSING(II),RIN%APSMIN(II),RIN%APSMAX(II),'  II,IDIM1,IDIM2,IDIM3:'
            endif
! ask Pavel, David
            !if(RIN%par_type(IDIM1) .eq. par_type_SURF1_RPV_BRDF) then 
              !if(IDIM2 .eq. 3 .and. RIN%APSING(II) .lt. 0.0) then
                !RIN%APSING(II) = RIN%APSING(II) + 1.0
                !RIN%APSMIN(II) = RIN%APSMIN(II) + 1.0
                !RIN%APSMAX(II) = RIN%APSMAX(II) + 1.0
              !endif
            !endif ! RIN%par_type(IDIM1) .eq. par_type_SURF1_RPV_BRDF
          enddo ! IDIM3
        enddo ! IDIM2  
      enddo ! IDIM1

!C*****************************************************************************************
!C    nx  -   number of pixels for X coordinate can be with data for X edge segments   
!C    ny  -   number of pixels for Y coordinate can be with data for Y edge segments  
!C    nt  -   number of pixels for T coordinate can be with data for T edge segments 
!C*****************************************************************************************
! Sizes of segment edges
      read(iu_iP_input,*)
      read(iu_iP_input,*) RIN%edges%nx,RIN%edges%ny,RIN%edges%nt
!C*****************************************************************************************
!C    INPUT  -   determines the initial guess:
!C            = F  INITIAL GUESS   for each pixel is equal the the values given above 
!C                 in INPUT... file 
!C            = T  INITIAL GUESS is read from file input_iguess.dat,
!C                 then INITIAL GUESS   for each pixel can be different
!C
!C*****************************************************************************************
! General output product
      read(iu_iP_input,*)
      ! Retrieval
      read(iu_iP_input,*) RIN%products%retrieval%res,RIN%products%retrieval%par,RIN%products%retrieval%fit  
      ! Aerosol
      read(iu_iP_input,*) RIN%products%aerosol%opt,RIN%products%aerosol%rind,RIN%products%aerosol%chem,  &
                          RIN%products%aerosol%phmx,RIN%products%aerosol%lidar,  &
                          RIN%products%aerosol%sd2m_mph,RIN%products%aerosol%sd2m_ext  
      ! Clouds
      read(iu_iP_input,*) RIN%products%clouds%opt,RIN%products%clouds%rind,RIN%products%clouds%chem,  &
                          RIN%products%clouds%phmx,RIN%products%clouds%lidar,  &
                          RIN%products%clouds%sd2m_mph,RIN%products%clouds%sd2m_ext  
      ! Surface 
      read(iu_iP_input,*) RIN%products%surface%surf  
      ! Radiative Forcing
      read(iu_iP_input,*) RIN%products%forcing%bbflux,RIN%products%forcing%forcing  
      ! Error estimates
      read(iu_iP_input,*) RIN%products%errest%par,  &
                          RIN%products%errest%aerosol%opt,RIN%products%errest%aerosol%lidar,  &
                          RIN%products%errest%clouds%opt, RIN%products%errest%clouds%lidar                            

      read(iu_iP_input,*) RIN%INPUT,RIN%Aexp_iwl(1:2),RIN%ndvi_iwl(1:2)

!C*****************************************************************************************
!C  for DLS package 
!C  DLSF :
!C    distname_O  - KERNEL directory name for 22/16 bins and for lognormal SD
!C    distname_N  - KERNEL directory name for precalculated lognormal bins
!C*****************************************************************************************

      read(iu_iP_input,'(a)') RIN%DLSF%internal_file_path
      read(iu_iP_input,'(a)') RIN%DLSF%external_file_path
!      write(*,*) 'in read_input_settings.f90 internal_file_path=',TRIM(RIN%DLSF%internal_file_path)
!      write(*,*) 'in read_input_settings.f90 external_file_path=',TRIM(RIN%DLSF%external_file_path)
      read(iu_iP_input,'(a)') RIN%DLSF%distname_O
      read(iu_iP_input,'(a)') RIN%DLSF%distname_N
!      write(*,*) 'in read_input_settings.f90 distname_O=',TRIM(RIN%DLSF%distname_O)
!      write(*,*) 'in read_input_settings.f90 distname_N=',TRIM(RIN%DLSF%distname_N)

!tl      if(iu_iP_input .ne. 5) close(iu_iP_input) 

      if(((RIN%DLSF%key .eq. 0.or.RIN%DLSF%key .eq. 3) .and. (RIN%DLSF%IWL .ne. 0)) .or.  &
	      ((RIN%DLSF%key .eq. 2).and.(RIN%DLSF%IWL .ne. 1))) then
        write(*,*) 'key=',RIN%DLSF%key,' IWL=',RIN%DLSF%IWL
        write(*,*) '((key.eq.0/3.and.IWL.ne.0).or.(key.eq.2.and.IWL.ne.1))'
        stop 'stop in read_input_settings'
      endif ! key and RIN%DLSF%IWL
!stop ''
      return
      end subroutine read_input_settings

!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
