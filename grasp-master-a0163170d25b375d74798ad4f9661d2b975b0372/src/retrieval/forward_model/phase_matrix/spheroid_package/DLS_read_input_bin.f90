! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../../../constants_set/mod_globals.inc"
      SUBROUTINE DLS_read_input_bin(file_path_git_track,distname_O,PHASE_MATRIX) 

      USE mod_par_DLS,     only : KMD,KNpar,KMpar,KRpar,rootdir						  
      USE mod_par_DLS_bin, only : KCpar,KWLpar           	 
      use mod_type_DLS
      use mod_stop_report

      implicit none
! --------------------------------------------------------------------------------
      character(len=255), intent(in)   ::  file_path_git_track
      character(len=255), intent(in)   ::  distname_O
      type(DLS_CONFIG),   intent(out)  ::  PHASE_MATRIX

      integer,parameter :: KSD=1 ! number of aerosol component, needs to be developed for KSD>1
      real,   dimension(KMD,KCpar)    :: RC1,XRM,XSM
      real,   dimension(KMD,KSD)      :: RCR
   
      real,  dimension(KMD,KSD)     :: CM,SM,RMM
      real,  dimension(KSD,KNpar)   :: RRR,AR,xgrid
      real,  dimension(KSD)         :: AC
      integer   :: ID,LB1,LE1,key_SD,NMD
      real      :: rgmin,rgmax,wlmin,wlmax
      integer   :: i,i1,i2,j,NSD
      real      :: RCsum,pi2

      CHARACTER(len=255) :: full_name
! --------------------------------------------------------------------------------
!        integer                   :: key,keyEL,keySUB,keyLS,key_RD1,key_RD, &    
!				                     KN,KM,KR          
!        real*4, dimension (KRpar) :: R
!		real,   dimension (KNpar) :: grid
!		real,   dimension (KRpar) :: RD
!		real,   dimension (KMpar) :: ANGLE 
!		real 	                  :: pomin,pomax
! --------------------------------------------------------------------------------
!   integer                         :: KC,NWL,LB,LE                        
																		
!   real,   dimension(KCpar)        :: RC2  
!   real,   dimension(KMpar,KWLpar) :: f11_bin,f12_bin,f22_bin    &
!									 ,f33_bin,f34_bin,f44_bin
!   real,   dimension(KWLpar)       :: WL_bin,RN_bin,RK_bin       & 
!									 ,xext_bin,xabs_bin,xsca_bin,albedo_bin
!   real,   dimension(KCpar,KNpar)  :: SD_bin
!   real*4                          :: xnmin,xnmax,xkmin,xkmax     
! --------------------------------------------------------------------------------
      !write(*,*) 'in read dls input: distname_O=',TRIM(distname_O)

      NSD=KSD
!
! ** READ INPUT
!       

!tl      write(*,*) 'in DLS_read_input: distname_O=',TRIM(distname_O)

      full_name=trim(file_path_git_track)//trim(distname_O)//'input_bin.dat'
!tl      write(*,*) 'file_name=',TRIM(full_name)

      open(10,file=trim(full_name),status='old')
      read(10,*) PHASE_MATRIX%DLSIN%key,    PHASE_MATRIX%DLSIN%keyEL, &
                 PHASE_MATRIX%DLSIN%keySUB, PHASE_MATRIX%DLSIN%keyLS, &
                 PHASE_MATRIX%DLSIN%key_RD1

! **   key_RD =1 - volume mixture of spheroids                  ** 
! **           2 - surface area  mixture of spheroids           ** 
      PHASE_MATRIX%DLSIN%key_RD=1
      if(PHASE_MATRIX%DLSIN%keySUB.eq.0) then
       write(*,*) 'key,keyEL,keySUB,keyLS,key_RD1,key_RD:'
	     write(*,*) PHASE_MATRIX%DLSIN%key,   PHASE_MATRIX%DLSIN%keyEL,   PHASE_MATRIX%DLSIN%keySUB,  &
	                PHASE_MATRIX%DLSIN%keyLS, PHASE_MATRIX%DLSIN%key_RD1, PHASE_MATRIX%DLSIN%key_RD
      endif ! keySUB
      read(10,*) rgmin,rgmax,wlmin,wlmax
      read(10,*) PHASE_MATRIX%DLSIN%xnmin, PHASE_MATRIX%DLSIN%xnmax,   &
	               PHASE_MATRIX%DLSIN%xkmin, PHASE_MATRIX%DLSIN%xkmax
      read(10,*) PHASE_MATRIX%DLSIN%NWL, LB1, LE1
	if(PHASE_MATRIX%DLSIN%key.eq.0.and.PHASE_MATRIX%DLSIN%NWL.ne.1) then
      write(tmp_message,'(2(a,i0))') '(key.eq.0.and.NWL.ne.1) !!!'
      G_ERROR(trim(tmp_message))
  endif
      if(PHASE_MATRIX%DLSIN%key.gt.0) then
	if(LB1.gt.LE1.or.LE1.gt.PHASE_MATRIX%DLSIN%NWL.or.LB1.gt.PHASE_MATRIX%DLSIN%NWL) then
      write(tmp_message,'(a)') 'check LB1, LE1 and NWL in input_bin.dat'
      G_ERROR(trim(tmp_message))
	endif
      endif
      read(10,*) (PHASE_MATRIX%DLSIN%WL_bin(i),i=1,PHASE_MATRIX%DLSIN%NWL)
      read(10,*) (PHASE_MATRIX%DLSIN%RN_bin(i),i=1,PHASE_MATRIX%DLSIN%NWL)
      read(10,*) (PHASE_MATRIX%DLSIN%RK_bin(i),i=1,PHASE_MATRIX%DLSIN%NWL)

    if(PHASE_MATRIX%DLSIN%keySUB.eq.0) then
      write(*,*) 'NWL=',PHASE_MATRIX%DLSIN%NWL
      write(*,*) 'WL=', PHASE_MATRIX%DLSIN%WL_bin(1:PHASE_MATRIX%DLSIN%NWL)
      write(*,*) 'RN=', PHASE_MATRIX%DLSIN%RN_bin(1:PHASE_MATRIX%DLSIN%NWL)
      write(*,*) 'RK=', PHASE_MATRIX%DLSIN%RK_bin(1:PHASE_MATRIX%DLSIN%NWL)
    endif ! keySUB
	
	read(10,*) key_SD,ID,NMD
	if(NMD.gt.1) then
    write(tmp_message,'(a)') 'NMD.gt.1 - Code does not support!!!'
    G_ERROR(trim(tmp_message))
	endif ! NMD
	if(key_SD.eq.0) then
    write(tmp_message,'(a)') 'key_SD.eq.0 - DOES NOT WORK for lognormal bins!!!'
    G_ERROR(trim(tmp_message))
      read(10,*)
      read(10,*) PHASE_MATRIX%DLSIN%KN
      if(PHASE_MATRIX%DLSIN%keySUB.eq.0) write(*,*) 'KN=',PHASE_MATRIX%DLSIN%KN
      do i=1,PHASE_MATRIX%DLSIN%KN
      read(10,*) PHASE_MATRIX%DLSIN%grid(i), PHASE_MATRIX%DLSIN%SD_bin(1,i)
      enddo ! i
	else ! key_SD=1
	read(10,*) (CM(i,1),SM(i,1),RMM(i,1),i=1,NMD)

      read(10,*) PHASE_MATRIX%DLSIN%KC,PHASE_MATRIX%DLSIN%RC2(1:PHASE_MATRIX%DLSIN%KC)
	 if(PHASE_MATRIX%DLSIN%key.ne.0) then
	 if(PHASE_MATRIX%DLSIN%key.eq.1) PHASE_MATRIX%DLSIN%RC2(1:PHASE_MATRIX%DLSIN%KC)=1.
       if(PHASE_MATRIX%DLSIN%keySUB.eq.0) then
	     write(*,*) 'RC2:  '
       write(*,'(5e12.4)') PHASE_MATRIX%DLSIN%RC2(1:PHASE_MATRIX%DLSIN%KC)
       endif ! keySUB
	 endif ! key.ne.0
	do i1=1,PHASE_MATRIX%DLSIN%KC
      read(10,*) (RC1(i,i1),XRM(i,i1),XSM(i,i1), i=1,NMD)
	enddo ! i1
	if(PHASE_MATRIX%DLSIN%key.eq.1) then
	  write(*,*) 'i1, RC1, XRM, XSM:'
	  do i1=1,PHASE_MATRIX%DLSIN%KC
	  write(*,*) i1,(RC1(i2,i1),XRM(i2,i1),XSM(i2,i1),i2=1,NMD)
	  enddo ! i1
	else if(PHASE_MATRIX%DLSIN%key.eq.2) then
	  if(PHASE_MATRIX%DLSIN%keySUB.eq.0) then
	  write(*,*) 'i1, RC2, XRM, XSM:'
	  do i1=1,PHASE_MATRIX%DLSIN%KC
	  write(*,*) i1,PHASE_MATRIX%DLSIN%RC2(i1),(XRM(i2,i1),XSM(i2,i1),i2=1,NMD)
	  enddo ! i1
	  endif ! keySUB	
	endif

	read(10,*) PHASE_MATRIX%DLSIN%KN
      if(PHASE_MATRIX%DLSIN%keySUB.eq.0) write(*,*) 'KN=',PHASE_MATRIX%DLSIN%KN
      do i=1,PHASE_MATRIX%DLSIN%KN
      read(10,*) PHASE_MATRIX%DLSIN%grid(i)
      enddo ! i
	endif ! key_SD
    xgrid(1,1:PHASE_MATRIX%DLSIN%KN)=PHASE_MATRIX%DLSIN%grid(1:PHASE_MATRIX%DLSIN%KN)

	if(PHASE_MATRIX%DLSIN%KN.gt.KNpar) then
    write(tmp_message,'(2(a,i0),3a)') 'in input.dat KN = ',PHASE_MATRIX%DLSIN%KN, &
    ' .ne. KNpar = ',KNpar,' in mod_par_DLS', &
    NEW_LINE('A'),'KN should be < or = KNpar'
    G_ERROR(trim(tmp_message))
	endif ! KM&KMpar

!      write(*,*) 'after grid'
      read(10,'(a)') PHASE_MATRIX%KERNEL%distname_O
	    read(10,'(a)') PHASE_MATRIX%KERNEL%distname_N
      if(PHASE_MATRIX%DLSIN%keySUB.eq.0) then 
	    write(*,*) 'distname_O = ',TRIM(PHASE_MATRIX%KERNEL%distname_O)
	    write(*,*) 'distname_N = ',TRIM(PHASE_MATRIX%KERNEL%distname_N)
      endif ! keySUB
	  
      read(10,*) PHASE_MATRIX%DLSIN%KR
      do i=1,PHASE_MATRIX%DLSIN%KR
      read(10,*) PHASE_MATRIX%DLSIN%R(i), PHASE_MATRIX%DLSIN%RD(i)
      enddo ! i
      read(10,*) PHASE_MATRIX%DLSIN%KM
	if(PHASE_MATRIX%DLSIN%KM.gt.KMpar) then
    write(tmp_message,'(2(a,i0),3a)') 'in input.dat KM = ', &
    PHASE_MATRIX%DLSIN%KM,' .gt. Kmpar = ',KMpar,' in mod_par_DLS', &
    NEW_LINE('A'),'KM should be < or = KMpar'
    G_ERROR(trim(tmp_message))
	endif ! KM&KMpar
      do j=1,PHASE_MATRIX%DLSIN%KM
      read(10,*) PHASE_MATRIX%DLSIN%ANGLE(j)
      enddo ! j

! READ common names for Kernel files
      READ(10,*) PHASE_MATRIX%XGRID1%NRATN
      IF(PHASE_MATRIX%XGRID1%NRATN.gt.KR1par) THEN
        write(tmp_message,'(a)') 'NRATN .gt. KR1par !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      DO i=1,PHASE_MATRIX%XGRID1%NRATN
      READ(10,*) PHASE_MATRIX%KERNEL%comm_name(i)
      ENDDO ! i

      close(10)

! READ grid file

      !full_name=TRIM(rootdir)//TRIM(PHASE_MATRIX%KERNEL%distname_O)//'/grid1.dat'
      full_name=trim(file_path_git_track)//trim(distname_O)//'grid1.dat'
!tl      write(*,*) 'file=',TRIM(full_name)
      OPEN (10,file=full_name)
      READ (10,*) PHASE_MATRIX%XGRID1%KN1, PHASE_MATRIX%XGRID1%XWL

      if(PHASE_MATRIX%XGRID1%KN1.gt.KN1par) then
        write(tmp_message,'(2(a,i0),2a)') 'in GET_MATRIX: KN1 = ', &
        PHASE_MATRIX%XGRID1%KN1,'  KN1par = ',KN1par, &
        NEW_LINE('A'),' !!! KN1.ne.KN1par'
        G_ERROR(trim(tmp_message))
      endif
      DO I=1,PHASE_MATRIX%XGRID1%KN1
      READ (10,*) PHASE_MATRIX%XGRID1%grid1(I)
      ENDDO ! I 
      READ (10,*) PHASE_MATRIX%XGRID1%KM1
      if(PHASE_MATRIX%XGRID1%KM1.gt.KM1par) then
        write(tmp_message,'(2(a,i0),2a)') 'KN1 = ', &
        PHASE_MATRIX%XGRID1%KN1,'  KN1par = ',KN1par, &
        NEW_LINE('A'),' !!! KM1.ne.KM1par'
        G_ERROR(trim(tmp_message))
      endif
      DO J=1,PHASE_MATRIX%XGRID1%KM1
      READ (10,*) PHASE_MATRIX%XGRID1%ANGLE1(J)
      ENDDO ! J
      CLOSE(10)
	IF((PHASE_MATRIX%DLSIN%key.eq.0).and.(PHASE_MATRIX%XGRID1%KM1.ne.PHASE_MATRIX%DLSIN%KM)) THEN
      write(tmp_message,'(a)') '(key .eq. 0) .and. (KM1 .ne. KM)'
      G_ERROR(trim(tmp_message))
	ENDIF


! **   key_grid1 read grid radii and scat.angles which were used** c
! **             for kernel look up tables or fixed kernels     ** c
! **          =0 - 'grid1.dat'                                  ** c
! **           1 - 'grid1.dat.fix'                              ** c

!      if(key.eq.2) then
!	 key_grid1=1
!	else 
!	 key_grid1=0
!	endif

      IF(PHASE_MATRIX%DLSIN%key.EQ.0) PHASE_MATRIX%DLSIN%KC=1
!
! ** CALCULATE SD if key_SD=1
!      
	  IF(PHASE_MATRIX%DLSIN%keySUB.EQ.0) THEN
      if(key_SD.eq.1) then
      OPEN (7,FILE='LNPAR.dat',status='unknown')
      OPEN (9,FILE='SizeDis.dat',status='unknown')
      DO i1=1,PHASE_MATRIX%DLSIN%KC
      if(PHASE_MATRIX%DLSIN%key.eq.0) then
       call SIZEDISDN ( 1,                                   & ! IPRI
                        -PHASE_MATRIX%DLSIN%KN,1,ID,NSD,NMD, &
                        CM (1:NMD,1:NSD),                    &
                        SM (1:NMD,1:NSD),                    &
                        RMM(1:NMD,1:NSD),                    &
                        xgrid(1:NSD,1),xgrid(1:NSD,PHASE_MATRIX%DLSIN%KN), &
                        RRR(1:NSD,:),AR(1:NSD,:),AC(1:NSD)   &
                      )
!       write(*,*) 'CM=',CM(1:NMD,1),' AC=',AC(1)
      else ! key=1 or 2
	   RCsum=sum(RC1(1:NMD,i1))
	   RCR(1:NMD,1)=RC1(1:NMD,i1)/RCsum

	   call SIZEDISDN ( 1,                          & ! IPRI
					   -PHASE_MATRIX%DLSIN%KN,1,ID,NSD,NMD, &
						RCR(1:NMD,1),                         &
						XSM(1:NMD,i1),                        &
						XRM(1:NMD,i1),                        &
						xgrid(1:NSD,1),xgrid(1:NSD,PHASE_MATRIX%DLSIN%KN), &
						RRR(1:NSD,:),AR(1:NSD,:),AC(1:NSD)    &
                    )
	 
!         write(*,*) i1,' NMD=',NMD,' CM=',RCR(1:NMD,1),' AC=',AC(1)
	  endif ! key
	  do i=1,PHASE_MATRIX%DLSIN%KN
!	    write(*,13) i,grid(i),RRR(1,i),SD_bin(i1,i),AR(1,i) 
	  PHASE_MATRIX%DLSIN%SD_bin(i1,i)=AR(1,i)
!	    write(*,13) i,RRR(1,i),SD_bin(i1,i)
	  enddo ! i   
      ENDDO ! i1  

      CLOSE(9)
	  CLOSE(7)	  
	  endif ! key_SD
!
! ** CHECK INPUT
!
      write(*,*) 'CHECK INPUT in OPTCHAR1:'
      write(*,*) 'key, keyEL, keySUB, keyLS, key_RD1, key_RD'
      write(*,'(i3,4i6)') PHASE_MATRIX%DLSIN%key,   PHASE_MATRIX%DLSIN%keyEL,   PHASE_MATRIX%DLSIN%keySUB,  &
	                      PHASE_MATRIX%DLSIN%keyLS, PHASE_MATRIX%DLSIN%key_RD1, PHASE_MATRIX%DLSIN%key_RD
      write(*,*) 'rgmin,rgmax,wlmin,wlmax:'
      write(*,'(2e12.5,2f7.4)')  rgmin,rgmax,wlmin,wlmax
      write(*,*) 'WL:'
      write(*,'(10f7.4)')  PHASE_MATRIX%DLSIN%WL_bin(1:PHASE_MATRIX%DLSIN%NWL)
	    write(*,*) 'RN:'
      write(*,'(10e12.4)') PHASE_MATRIX%DLSIN%RN_bin(1:PHASE_MATRIX%DLSIN%NWL)
	    write(*,*) 'RK:'
      write(*,'(10e12.4)') PHASE_MATRIX%DLSIN%RK_bin(1:PHASE_MATRIX%DLSIN%NWL)
      write(*,*) 'size distribution:'
	  if(PHASE_MATRIX%DLSIN%key .ne. 0) write(*,'(a,10i12)')  &
              'lognormal bin N:',(i1,i1=1,PHASE_MATRIX%DLSIN%KC)
	  do i=1,PHASE_MATRIX%DLSIN%KN
	  write(*,13) i,PHASE_MATRIX%DLSIN%grid(i),  &
	  (PHASE_MATRIX%DLSIN%SD_bin(i1,i),i1=1,PHASE_MATRIX%DLSIN%KC)
	  enddo ! i
	
      write(*,*) 'axis ratio distribution:'
      do i=1,PHASE_MATRIX%DLSIN%KR
      write(*,13) i,PHASE_MATRIX%DLSIN%R(i), PHASE_MATRIX%DLSIN%RD(i)
      enddo ! i
      write(*,*) 'SCATTERING ANGLES: KM=',PHASE_MATRIX%DLSIN%KM
      write(*,14) PHASE_MATRIX%DLSIN%ANGLE(1:PHASE_MATRIX%DLSIN%KM)
      if(PHASE_MATRIX%DLSIN%keyEL.eq.1)   &
            write(*,*) 'F11 will be calculated'
      if(PHASE_MATRIX%DLSIN%keyEL.eq.2)   &
            write(*,*) 'F11 F12 will be calculated'
      if(PHASE_MATRIX%DLSIN%keyEL.eq.3)   &
			write(*,*) 'F11 F12 F22 will be calculated'
      if(PHASE_MATRIX%DLSIN%keyEL.eq.4)   &
            write(*,*) 'F11 F12 F22 F33 will be calculated'
      if(PHASE_MATRIX%DLSIN%keyEL.eq.5)   &
            write(*,*) 'F11 F12 F22 F33 F34 will be calculated'
      if(PHASE_MATRIX%DLSIN%keyEL.eq.6)   &
            write(*,*) 'F11 F12 F22 F33 F34 F44 will be calculated'

	  ENDIF ! keySUB=0
!
! ** rgrid min & max calculation for fixed or NEWORG kernels
!
	pi2=2.*ACOS(-1.)
	PHASE_MATRIX%DLSIN%pomin=pi2*rgmin/wlmax
	PHASE_MATRIX%DLSIN%pomax=pi2*rgmax/wlmin  

	PHASE_MATRIX%DLSIN%LB=1
	PHASE_MATRIX%DLSIN%LE=PHASE_MATRIX%DLSIN%NWL
	if(PHASE_MATRIX%DLSIN%key.eq.2) then
	 PHASE_MATRIX%DLSIN%LB=LB1
	 PHASE_MATRIX%DLSIN%LE=LE1
	endif
	
13    FORMAT(i4,22E12.4)
14    FORMAT(7F8.2)

      RETURN
      END SUBROUTINE DLS_read_input_bin


