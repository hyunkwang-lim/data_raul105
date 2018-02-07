! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../../../constants_set/mod_globals.inc"

      SUBROUTINE MATRIX_FIX_bin(key,key_RD,keyEL,keySUB             &
                               ,KR                                  &
                               ,KN,grid,KM,ANGLE                    &
                               ,WAVEL,KRE,KIM,ARE,AIM,pomin,pomax   &
                               ,KERNEL,XGRID1                       &
                               ,NWL,WL,KC,SD                        &
                               ,xnmin,xnmax,xkmin,xkmax             & 
                               ,KERNELS1,KERNELS2)
!c ** ANGLE SPLINE interpolation version
!c ** version for Oleg
!c **
!c ** rootdir is defined in "optchar.par"
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 05/05/03 this version can be used to retrieve an aspect ratio
!c ** 13/10/03 IF(ANGLE<=40.)ln(f33)-interpol.
!c **          IF(ANGLE<=50.)ln(f44)-interpol.
!c **************************************************************** c
!c **   Subroutine gets original and calculates fixed kernel     ** c
!c **   matrices for given                                       ** c 
!c **   aspect ratio distribution and scattering angles          ** c
!c **                                                            ** c
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c
!c **   key  = 0 - calculate opt.characteristics from original   ** c
!c **              kernels 'Rke...'                              ** c
!c **          1 - read original kernels, create WL kernels      ** c
!c **              and save them in dirname_N  and calculate     ** c
!c **              opt.characteristics                           ** c
!c **          2 - read WL kernels and calculate                 ** c
!c **              opt.characteristics                           ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   keyEL=  1 - calculate F11                                ** c
!c **           2 - F11,F12                                      ** c
!c **           3 - F11,F12,F22                                  ** c
!c **           4 - F11,F12,F22,F33                              ** c
!c **           5 - F11,F12,F22,F33,F34                          ** c
!c **           6 - F11,F12,F22,F33,F34,F44                      ** c
!c **   KR  - number of aspect ratios                            ** c
!c **   R(KR)  - grid ratios                                     ** c
!c **   RD(KR) - aspect ratio distribution for grid aspect ratios** c
!c **   KM   - number of scattering angles for fixed kernels     ** c
!c **   ANGLE(KM)- scatering angles for fixed kernels            ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   UF...(KMpar,KN1par,KIMpar,KREpar) - kernels for          ** c
!c **           given aspect ratio distribution                  ** c
!c **                                                            ** c
!c **   UFEA(2,KN1par,KIMpar,KREpar) - extinction and absorption ** c
!c **                                                    kernels ** c
!c **              1 - extinction                                ** c
!c **              2 - absorption                                ** c
!c **   KN1 - number of grid radii from original or fixed kernels** c
!c **   grid1(KN1) - grid radii                                  ** c
!c **   WAVEL - wavelength from original or fixed kernels        ** c
!c **   KRE   - number of real parts of refr.ind.                ** c
!c **   KIM   - number of imaginary parts of refr.ind.           ** c
!c **   ARE(KRE) - real parts of refr.ind.                       ** c
!c **   AIM(KIM) - imaginary parts of refr.ind.                  ** c
!c **************************************************************** c
!c **************************************************************** c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.dat'                                  ** c
!c **           1 - 'grid1.dat.fix'                              ** c

      !use alloc1
      !use alloc2
      use mod_alloc_kernels 
      USE mod_par_DLS
      USE mod_par_DLS_bin
      use mod_intrpl_spline
      use mod_type_DLS
      use mod_stop_report
!c -------------------------------------------------------------------
      implicit none
!c -------------------------------------------------------------------	  
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2

!      real      :: RD(KRpar)
      integer   :: NEL(0:6)
      character(len=2) NELC(0:6)
!      real*4    :: R(KRpar)
      real      :: grid(KNpar)     &
                  ,ANGLE(KMpar)    &
                  ,ANGLE2(KM1par)  &
                  ,RAB1(KM1par)    &
                  ,RAB2(KM1par,KN1par,KIMpar,KREpar)
      real      :: ARE(KREpar), AIM(KIMpar)  &
                  ,SD(KCpar,KNpar)
      real      :: WAVEL,WL(KWLpar), dlnr
      integer   :: key,key_RD,keyEL,keySUB
      integer   :: KM1,KN1,KN,KKEL,KNWL
      CHARACTER(255)  name,dir_name_O,dir_name_N,  &
                     distname_O,distname_N
      CHARACTER(255) full_name
      type(DLS_kernel_location) :: KERNEL
      type(DLS_kernel_grid1)    :: XGRID1
      CHARACTER(60)  comm_name1	  
      CHARACTER(15)  NBINC,NWLC
	  
      integer,save :: NDP=0

      DOUBLE PRECISION, DIMENSION(KM1par) :: XXS2, YYS2, b, c, d
      DOUBLE PRECISION :: XARG,YFIT
      real(4) :: xnmin,xnmax,xkmin,xkmax

!c ----- Timer ------
      real(4) tarray(2),T_RFM,T_CFM
      real(4) dtime
!c -------------------------------------------------------------------
      integer :: KR,KM,NWL,KC,nn1,nn2,nn3,nn,NDPP,ind
      integer :: KRE,KIM,I,I1,L,J,KNN,KMM,  &
	              IRE,IIM,IRATN,KEL,NRATN,NNEL,ierr
      real    :: pomin,pomax,coeff,xx,WAVEL1,RREMIN,RREMAX,  &
	              rmin,rmax,RIMMIN,RIMMAX,rgmin,rgmax,PI2
!c -------------------------------------------------------------------
    
      NRATN      = XGRID1%NRATN
      
      KM1 = XGRID1%KM1
      KN1 = XGRID1%KN1
	  
      !write(*,*)'distname_O=',trim(KERNEL%internal_file_path)
      !write(*,*)'distname_N=',trim(KERNEL%external_file_path)
      !write(*,*)'distname_O=',trim(KERNEL%distname_O)
      !write(*,*)'distname_N=',trim(KERNEL%distname_N)

      write(NBINC,*) KC
      write(NWLC,*)  NWL
!       write(*,*) 'NBINC : ',trim(adjustl(NBINC))
!       write(*,*) 'NWLC  : ',trim(adjustl(NWLC))	         		

      dir_name_O=trim(KERNEL%internal_file_path)//trim(KERNEL%distname_O) !//'/'
      dir_name_N=trim(KERNEL%internal_file_path)//trim(KERNEL%distname_N) !//'/'

!tl      write(*,*)'internal_file_path=',KERNEL%internal_file_path
!tl      write(*,*)'dir_name_O=',TRIM(dir_name_O)
!tl      write(*,*)'dir_name_N=',TRIM(dir_name_N)

      PI2=2.*ACOS(-1.)
      NEL  = (/0, 11, 12, 22, 33, 34, 44/)
      NELC = (/'00','11','12','22','33','34','44'/)
      T_RFM=0.
      T_CFM=0.
!c
!c ** Redefine grids for NEW WL kernels
!c
      NDPP=0
      if(key .eq. 0) then
      rgmin=pomin*XGRID1%XWL/pi2
      rgmax=pomax*XGRID1%XWL/pi2 
!c
      if(rgmin .lt. XGRID1%grid1(1)) then
       write(tmp_message,'(a,es11.4,a,a,2(a,es11.4))') &
       'XWL =',XGRID1%XWL, &
      ' IS wavelength in kernels equal to XWL?', &
       NEW_LINE('A'), &
       'check input.dat: rgmin =',rgmin,' < grid1(1) =',XGRID1%grid1(1)
       G_ERROR(trim(tmp_message))
      endif
      if(rgmax .gt. XGRID1%grid1(KN1)) then
       write(tmp_message,'(a,es11.4,a,a,2(a,es11.4))') &
       'XWL =',XGRID1%XWL,' IS wavelength in kernels equal to XWL?', &
        NEW_LINE('A'), &
       'check input.dat: rgmax=',rgmax,' > grid1(KN1) =',XGRID1%grid1(KN1)
       G_ERROR(trim(tmp_message))
      endif

      GOTO 700
      ind=0
      do i=1,KN1-1
!c	write(*,*) 'i=',i,' rgmin=',rgmin,' grid1(i+1)=',XGRID1%grid1(i+1)
      if(rgmin .le. XGRID1%grid1(i+1)) then
!c	xr1=XGRID1%grid1(i)
      nn1=i
      ind=ind+1
      endif
      if(ind .eq. 1) EXIT
      enddo ! i
	 
      ind=0
      do i=KN1,2,-1
!c	write(*,*) 'i=',i,' rgmax=',rgmax,' grid1(i-1)=',XGRID1%grid1(i-1)
      if(rgmax.ge.XGRID1%grid1(i-1)) then
!c	xr2=XGRID1%grid1(i)
      nn2=i
      ind=ind+1
      endif
      if(ind.eq.1) EXIT
      enddo ! i

      if(rgmin.eq.XGRID1%grid1(1)) then
!c	xr1=XGRID1%grid1(1)
      nn1=1
      endif ! rgmin 
      if(rgmax.eq.XGRID1%grid1(KN1)) then
!c	xr2=XGRID1%grid1(KN1)
      nn2=KN1
      endif ! rgmax
      nn3=nn2-nn1+1
!c	write(*,*) 'nn1=',nn1,' nn2=',nn2,' nn3=',nn3
!c	write(*,*) 'xr1=',xr1,' xr2=',xr2,' KN1=',KN1

!cl	STOP 'TEST STOP'
	 
700   CONTINUE
      endif ! key=0

      NNEL=keyEL+1
      KKEL=6-keyEL
!c
!c ** Read KC WL original kernels UO...
!c
      IF (key .eq. 2) THEN
	  	  
      DO IRATN=1,NRATN
!c **
!c ** READ UO11, UO12, UO22, UO33, UO34, UO44 matrices
!c **
      comm_name1 = trim(KERNEL%comm_name(IRATN))//'_bin'//  &
                   trim(adjustl(NBINC))//'_wl'//trim(adjustl(NWLC))//'_'              

      DO KEL=1,NNEL-1		
      if(KEL .gt. keyEL) CYCLE

      name=trim(comm_name1)//NELC(KEL)//'.txt'
      full_name=TRIM(dir_name_N)//TRIM(name)
!tl      write(*,*) 'name=',name
!tl      write(*,*) 'dir_name_N=',dir_name_N
!tl      write(*,*) 'full_name=',full_name

      OPEN(11,FILE=full_name,status='old')
      READ(11,*) XGRID1%RATIO(iratn)
      READ(11,*) KNN,KNWL
      READ(11,*) KMM
      IF(KNN .ne. KC)   THEN
        write(tmp_message,'(a)') 'KNN .ne. KC !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(KNWL .ne. NWL) THEN
        write(tmp_message,'(a)') 'KNWL .ne. NWL !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(KNWL .gt. KWLpar) THEN
        write(tmp_message,'(a)') 'KNWL .gt. KWLpar !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(KMM .ne. KM1)  THEN
        write(tmp_message,'(a)') 'KMM .ne. KM1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF

      READ(11,*) ANGLE2(1:KMM)
!cl        WRITE(*,*)  'ANGLE2 in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) ANGLE2(:KM)
!cl        WRITE(*,*)  'ANGLE1  in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) XGRID1%ANGLE1(:KM)

!tl      WRITE(*,*) 'READ matrix UO',NEL(KEL)
      READ(11,*) RREMIN, RREMAX
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
      READ(11,*) RIMMIN, RIMMAX
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
      READ(11,*) KRE, KIM
!cl      WRITE(*,*) KRE, KIM,' KRE, KIM'

      IF(KRE .LT. 0) KRE=-KRE
      IF(KIM .LT. 0) KIM=-KIM 
      DO  L=1,NWL
      DO  IRE=1,KRE
      DO  IIM=1,KIM
         READ(11,*) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
         READ(11,*) WL(L),ARE(IRE),AIM(IIM) 
         AIM(IIM)=-AIM(IIM)
         IF(KEL .EQ. 1) then           
            DO  I=1,KC
            READ(11,*) (RAB1(J),J=1,KM)
            DO J=1,KM
               KERNELS2%UO11(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ELSE IF(KEL .EQ. 2) THEN
            DO  I=1,KC
            READ(11,*) RAB1(1:KM)
            if(XGRID1%ANGLE1(1).eq.0.)     RAB1(1)  = 0.
            if(XGRID1%ANGLE1(KM1).eq.180.) RAB1(KM) = 0.
            DO J=1,KM 
               KERNELS2%UO12(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ELSE IF(KEL .EQ. 3) THEN
            DO  I=1,KC
            READ(11,*) (RAB1(J),J=1,KM)
            DO J=1,KM
               KERNELS2%UO22(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ELSE IF(KEL .EQ. 4) THEN
            DO  I=1,KC
            READ(11,*) (RAB1(J),J=1,KM)
            DO J=1,KM
               KERNELS2%UO33(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ELSE IF(KEL .EQ. 5) THEN
            DO  I=1,KC
            READ(11,*) (RAB1(J),J=1,KM)
            DO J=1,KM
               KERNELS2%UO34(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ELSE IF(KEL .EQ. 6) THEN
            DO  I=1,KC
            READ(11,*) (RAB1(J),J=1,KM)
            DO J=1,KM
               KERNELS2%UO44(IRATN,J,I,IIM,IRE,L)=RAB1(J)
            ENDDO ! J
            ENDDO ! I
         ENDIF ! KEL
      ENDDO ! IIM
      ENDDO ! IRE
      ENDDO ! L NWL
      CLOSE (11)
      ENDDO ! KEL	  	
!cl         write(*,*) 'after read IRATN=',IRATN,
!cl     &   ' UO11=',KERNELS2%UO11(IRATN,1,1,1,1,1),
!cl     &   ' UO11=',KERNELS2%UO11(IRATN,KM1,KN1,KIM,KRE,NWL)
!c **
!c ** READ UOEA MATRIX
!c **              
      name=trim(comm_name1)//NELC(0)//'.txt'

      full_name=TRIM(dir_name_N)//TRIM(name)
!tl      write(*,*) 'name=',name,'full_name=',full_name
      OPEN(11,FILE=full_name,status='old')

      READ(11,*) XGRID1%RATIO(iratn)
!cl      WRITE(*,*) RATIO(iratn),' ratio'
      READ(11,*) KNN
      IF(KNN .ne. KC) THEN
        write(tmp_message,'(2i5,3a)') KNN,KC,' KNN,KC', &
        NEW_LINE('A'),'2: KNN.ne.KC !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
!tl      WRITE(*,*) 'READ matrix UO',NEL(0)
      READ(11,*) RREMIN, RREMAX
!cl      WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'  
      READ(11,*) RIMMIN, RIMMAX
!cl      WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
      READ(11,*) KRE, KIM
      IF(KRE .LT. 0) KRE=-KRE
      IF(KIM .LT. 0) KIM=-KIM 
      DO L=1,NWL
      DO IRE=1,KRE
      DO IIM=1,KIM
      READ(11,*) nn,xx
!cl      WRITE(*,*) nn,xx,' KEL,RATIO'
      READ(11,*) WAVEL1,ARE(IRE),AIM(IIM)
      IF(WL(L) .ne. WAVEL1) THEN
        write(tmp_message,'(i5,2es11.4,3a)') &
        L, WL(L), WAVEL1,' L,WL(L),WAVEL1', &
        NEW_LINE('A'),'WL(L).ne.WAVEL1'
        G_ERROR(trim(tmp_message))
      ENDIF
      AIM(IIM)=-AIM(IIM)             

      READ(11,*)           
      READ(11,*) (KERNELS2%UOEA(IRATN,1,I,IIM,IRE,L),I=1,KC)
      READ(11,*)           
      READ(11,*) (KERNELS2%UOEA(IRATN,2,I,IIM,IRE,L),I=1,KC)

      ENDDO ! IIM
      ENDDO ! IRE
      ENDDO ! L NWL

      CLOSE (11)
      ENDDO ! IRATN
      CLOSE(10)
!c END for Read WL original kernels
	
!c **************************************************


!13    FORMAT(3E12.4,I4)                    
!14    FORMAT('II=',i2,' key_RD=',i2,' key_RD1=',i2)
!16    FORMAT(12x,'wl=',f5.2,2x,'n=',f8.5)
!cl      WRITE(*,*) 'Fixed kernel matrices have been read'
!tl      IF(key_RD.EQ.1) WRITE(*,*) 'Volume mixture of spheroids'
!tl      IF(key_RD.EQ.2) WRITE(*,*) 'Surface area mixture of spheroids'
										if(keySUB.eq.0) then
                                            T_RFM=dtime(tarray) !+++
                                            T_RFM=tarray(1)+tarray(2)
      write(6,*) 
      write(6,*)  &
              '------------------ T I M I N G ------------------' 
      WRITE(*,61) T_RFM/60.     
										endif
   61 format('  Read WL original kernels. ........ ',f8.3,' min.')

      RETURN
      ENDIF ! key=2

      IF(key .NE. 2) THEN     
      IF(key .EQ. 1) THEN
!c
!c *** ALLOCATE and INITIALIZE ARRAYS
!c
        ALLOCATE(KERNELS1%UEA(KR1par,2,KN1par,KIMpar,KREpar),    stat=ierr)
        if(ierr/=0) stop 'Can not allocate UEA array' 
        KERNELS1%UEA=0.
        if(keyEL .gt. 0) then
        ALLOCATE(KERNELS1%U11(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
        if(ierr/=0) stop 'Can not allocate U11 array'
        KERNELS1%U11=0.
        endif
        
        if(keyEL .gt. 1) then
          ALLOCATE(KERNELS1%U12(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U12 array'
          KERNELS1%U12=0.
        endif
        if(keyEL .gt. 2) then
          ALLOCATE(KERNELS1%U22(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U22 array'
          KERNELS1%U22=0.
        endif
        if(keyEL .gt. 3) then
          ALLOCATE(KERNELS1%U33(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U33 array'
          KERNELS1%U33=0.
		endif
        if(keyEL .gt. 4) then
          ALLOCATE(KERNELS1%U34(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U34 array'
          KERNELS1%U34=0.
        endif
        if(keyEL .gt. 5) then
          ALLOCATE(KERNELS1%U44(KR1par,KM1par,KN1par,KIMpar,KREpar),stat=ierr)
          if(ierr/=0) stop 'Can not allocate U44 array'
          KERNELS1%U44=0.
        endif 

      ENDIF ! key=1

      IF(NDP .EQ. 0) THEN
      T_CFM=0.
!c
!c ** READ ORIGINAL kernels 
!c              
      DO IRATN=1,NRATN   ! Aspect ratio loop

      comm_name1 = trim(KERNEL%comm_name(IRATN))//'_'

      IF(key .EQ. 0) THEN
!c **
!c ** READ U11, U12, U22, U33, U34, U44 matrices
!c **
      DO KEL=1,NNEL-1
      if(KEL .gt. keyEL) CYCLE
      name=trim(comm_name1)//NELC(KEL)//'.txt'
      full_name=TRIM(dir_name_O)//TRIM(name)
      OPEN(11,FILE=full_name,status='old')
!tl      write(*,*) 'full_name=',full_name
      READ(11,*) rmin,rmax,XGRID1%RATIO(iratn)
      READ(11,*) KNN
      IF(KNN .LT. 0) KNN=-KNN
      READ(11,*) KMM
      IF(KNN .ne. KN1) THEN
        write(tmp_message,'(a)') 'KNN.ne.KN1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(KMM .ne. KM1) THEN
        write(tmp_message,'(a)') '1: KMM.ne.KM1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      READ(11,*) ANGLE2(1:KMM)
!cl        WRITE(*,*)  'ANGLE2 in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) ANGLE2(:KM1)
!cl        WRITE(*,*)  'ANGLE1  in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) XGRID1%ANGLE1(:KM1)

!tl      WRITE(*,*) 'READ matrix U',NEL(KEL)
      READ(11,*) RREMIN, RREMAX
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
      READ(11,*) RIMMIN, RIMMAX
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
      READ(11,*) KRE, KIM
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'

      IF(KRE .LT. 0) KRE=-KRE
      IF(KIM .LT. 0) KIM=-KIM 
      DO  IRE=1,KRE
      DO  IIM=1,KIM
      READ(11,*) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
      READ(11,*) WAVEL,ARE(IRE),AIM(IIM) 
      AIM(IIM)=-AIM(IIM)
     
      IF(KEL .EQ. 1) then           
         DO  I=1,KN1
            READ(11,*) (KERNELS1%U11(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
      ELSE IF(KEL .EQ. 2) THEN
         DO  I=1,KN1
            READ(11,*) (KERNELS1%U12(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
         if(XGRID1%ANGLE1(1)  .eq.  0.) then 
         do I=1,KN1
         KERNELS1%U12(IRATN,1,I,IIM,IRE) = 0.
			enddo ! j
         endif
         if(XGRID1%ANGLE1(KM1).eq.180.) then
         do I=1,KN1
         KERNELS1%U12(IRATN,KM1,I,IIM,IRE) = 0.
			enddo ! j
         endif 
      ELSE IF(KEL .EQ. 3) THEN
         DO  I=1,KN1
         READ(11,*) (KERNELS1%U22(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
      ELSE IF(KEL .EQ. 4) THEN
         DO  I=1,KN1
         READ(11,*) (KERNELS1%U33(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
      ELSE IF(KEL .EQ. 5) THEN
         DO  I=1,KN1
         READ(11,*) (KERNELS1%U34(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
         if(XGRID1%ANGLE1(1)  .eq.  0.) then 
         do I=1,KN1
         KERNELS1%U34(IRATN,1,I,IIM,IRE) = 0.
			enddo ! j
         endif
         if(XGRID1%ANGLE1(KM1).eq.180.) then
         do I=1,KN1
         KERNELS1%U34(IRATN,KM1,I,IIM,IRE) = 0.
			enddo ! j
         endif 
      ELSE IF(KEL .EQ. 6) THEN
         DO  I=1,KN1
         READ(11,*) (KERNELS1%U44(IRATN,J,I,IIM,IRE),J=1,KM1)
         ENDDO ! I
         ENDIF ! KEL
         ENDDO ! IIM
         ENDDO ! IRE
         CLOSE (11)
      ENDDO ! KEL	
!cl         write(*,*) 'after read IRATN=',IRATN,
!cl     &   ' U11=',U11(IRATN,1,1,1,1),
!cl     &   ' U11=',U11(IRATN,KM1,KN1,KIM,KRE)
      ENDIF ! key=0

      IF(key .EQ. 1) THEN
!c **
!c ** READ U11, U12, U22, U33, U34, U44 matrices
!c **
      DO KEL=1,NNEL-1
      if(KEL .gt. keyEL) CYCLE
      name=trim(comm_name1)//NELC(KEL)//'.txt'
      full_name=TRIM(dir_name_O)//TRIM(name)
      OPEN(11,FILE=full_name,status='old')
!tl      write(*,*) 'full_name',full_name
      READ(11,*) rmin,rmax,XGRID1%RATIO(iratn)
      READ(11,*) KNN
      IF(KNN .LT. 0) KNN=-KNN
      READ(11,*) KMM
      IF(KNN .ne. KN1) THEN
      STOP ' in MATRIX_FIX_bin 1: KNN.ne.KN1 !!!'
        write(tmp_message,'(a)') '1: KNN.ne.KN1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(KMM .ne. KM1) THEN
        write(tmp_message,'(a)') '1: KMM.ne.KM1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
      READ(11,*) ANGLE2(1:KMM)
!cl        WRITE(*,*)  'ANGLE2 in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) ANGLE2(:KM1)
!cl        WRITE(*,*)  'ANGLE1  in MATRIX_FIX_bin NEL=',NEL(KEL)
!cl        WRITE(*,11) XGRID1%ANGLE1(:KM1)

!tl      WRITE(*,*) 'READ matrix U',NEL(KEL)
      READ(11,*) RREMIN, RREMAX
!cl        WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'   
      READ(11,*) RIMMIN, RIMMAX
!cl        WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'        
      READ(11,*) KRE, KIM
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'

      IF(KRE .LT. 0) KRE=-KRE
      IF(KIM .LT. 0) KIM=-KIM 
      DO  IRE=1,KRE
      DO  IIM=1,KIM
      READ(11,*) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
      READ(11,*) WAVEL,ARE(IRE),AIM(IIM) 
      AIM(IIM)=-AIM(IIM)
     
      IF(KEL .EQ. 1) then           
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         RAB2(1:KM1,I,IIM,IRE)=RAB1(1:KM1)
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=LOG(RAB1(1:KM1))
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         KERNELS1%U11(IRATN,J,I,IIM,IRE)=EXP(YFIT)
         ENDDO ! J
         ENDDO ! I
      ELSE IF(KEL .EQ. 2) THEN
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         if(XGRID1%ANGLE1(1) .eq. 0.)     RAB1(1)  = 0.
         if(XGRID1%ANGLE1(KM1) .eq. 180.) RAB1(KM1)= 0.
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         if(J .eq. 1 .and. XGRID1%ANGLE1(J) .eq. 0.)       YFIT = 0.
         if(J .eq. KM1 .and. XGRID1%ANGLE1(KM1) .eq. 180.) YFIT = 0.
         KERNELS1%U12(IRATN,J,I,IIM,IRE)=YFIT*KERNELS1%U11(IRATN,J,I,IIM,IRE)
         ENDDO ! J
         ENDDO ! I
      ELSE IF(KEL .EQ. 3) THEN
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=LOG(RAB1(1:KM1))
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         KERNELS1%U22(IRATN,J,I,IIM,IRE)=EXP(YFIT)
         ENDDO ! J
         ENDDO ! I
      ELSE IF(KEL .EQ. 4) THEN
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         KERNELS1%U33(IRATN,J,I,IIM,IRE)=YFIT*KERNELS1%U11(IRATN,J,I,IIM,IRE)
         ENDDO ! J
         ENDDO ! I
      ELSE IF(KEL .EQ. 5) THEN
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         if(XGRID1%ANGLE1(1) .eq. 0.)     RAB1(1)  = 0.
         if(XGRID1%ANGLE1(KM1) .eq. 180.) RAB1(KM1)= 0.
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         if(J .eq. 1 .and. XGRID1%ANGLE1(J) .eq. 0.)     YFIT = 0.
         if(J .eq. KM1 .and. XGRID1%ANGLE1(J) .eq. 180.) YFIT = 0.
         KERNELS1%U34(IRATN,J,I,IIM,IRE)=YFIT*KERNELS1%U11(IRATN,J,I,IIM,IRE)
         ENDDO ! J
         ENDDO ! I
      ELSE IF(KEL .EQ. 6) THEN
         DO  I=1,KN1
         READ(11,*) RAB1(1:KM1)
         XXS2(1:KM1)=XGRID1%ANGLE1(1:KM1)
         YYS2(1:KM1)=RAB1(1:KM1)/RAB2(1:KM1,I,IIM,IRE)
         DO J=1,KM
         XARG=ANGLE(J)
         CALL intrpl_spline(KM1, XXS2(1:KM1), YYS2(1:KM1) ,XARG, YFIT, 1, J, b, c, d)
         KERNELS1%U44(IRATN,J,I,IIM,IRE)=YFIT*KERNELS1%U11(IRATN,J,I,IIM,IRE)
         ENDDO ! J
         ENDDO ! I
      ENDIF ! KEL
      ENDDO ! IIM
      ENDDO ! IRE
      CLOSE (11)
      ENDDO ! KEL	
      
      ENDIF ! key=1
!cl         write(*,*) 'after read IRATN=',IRATN,
!cl     &   ' U11=',U11(IRATN,1,1,1,1),
!cl     &   ' U11=',U11(IRATN,KM1,KN1,KIM,KRE)
!c **
!c ** READ UEA MATRIX
!c **
	  name=trim(comm_name1)//NELC(0)//'.txt'
      full_name=TRIM(dir_name_O)//TRIM(name)
!tl      write(*,*) 'full_name=',full_name
      OPEN(11,FILE=full_name,status='old')

      READ(11,*) rmin,rmax,XGRID1%RATIO(iratn)
!cl      WRITE(*,*) rmin,rmax,RATIO(iratn),' rmin,rmax,ratio'
      READ(11,*) KNN
      IF(KNN .LT. 0) KNN=-KNN
      IF(KNN .ne. KN1) THEN
        write(tmp_message,'(2i5,3a)') KNN,KN1,' KNN,KN1', &
        NEW_LINE('A'),'  2: KNN.ne.KN1 !!!'
        G_ERROR(trim(tmp_message))
      ENDIF
!tl      WRITE(*,*) 'READ matrix U',NEL(0)
      READ(11,*) RREMIN, RREMAX
!cl      WRITE(*,*) RREMIN, RREMAX,' RREMIN,RREMAX'  
      READ(11,*) RIMMIN, RIMMAX
!cl      WRITE(*,*) RIMMIN, RIMMAX,' RIMMIN,RIMMAX'
      READ(11,*) KRE, KIM

      IF(KRE .LT. 0) KRE=-KRE
      IF(KIM .LT. 0) KIM=-KIM 
      DO IRE=1,KRE
      DO IIM=1,KIM
      READ(11,*) nn,xx
!cl      WRITE(*,*) nn,xx,' KEL,RATIO'
      READ(11,*) WAVEL1,ARE(IRE),AIM(IIM)
      IF(keyEL .ge. 1 .and. WAVEL .ne. WAVEL1) THEN
        write(tmp_message,'(2es11.4,3a)') WAVEL,WAVEL1,' WAVEL, WAVEL1', &
        NEW_LINE('A'),'  WAVEL.ne.WAVEL1'
        G_ERROR(trim(tmp_message))
      ELSE
      WAVEL=WAVEL1
      ENDIF
       
      AIM(IIM)=-AIM(IIM)             

      READ(11,*)           
      READ(11,*) (KERNELS1%UEA(IRATN,1,I,IIM,IRE),I=1,KN1)
      READ(11,*)           
      READ(11,*) (KERNELS1%UEA(IRATN,2,I,IIM,IRE),I=1,KN1)
      ENDDO ! IIM
      ENDDO ! IRE
      CLOSE (11)	

      ENDDO ! IRATN   End of Aspect ratio loop
      CLOSE(10)

      IF(keyEL .GT. 0) THEN
!c
!c ** LOG(U...)
!c
      DO   IRE=1,KRE
      DO   IIM=1,KIM
      DO     I=1,KN1
      DO     J=1,KM
      DO IRATN=1,KR

      KERNELS1%U11(IRATN,J,I,IIM,IRE)=               &
      LOG( KERNELS1%U11(IRATN,J,I,IIM,IRE)*WAVEL )

      if(keyEL .gt. 1)                      &
         KERNELS1%U12(IRATN,J,I,IIM,IRE)=      &
         KERNELS1%U12(IRATN,J,I,IIM,IRE)*WAVEL 
      if(keyEL .gt. 2)                      &
         KERNELS1%U22(IRATN,J,I,IIM,IRE)=      &
         LOG( KERNELS1%U22(IRATN,J,I,IIM,IRE)*WAVEL )
      if(keyEL .gt. 3) then
         KERNELS1%U33(IRATN,J,I,IIM,IRE)=  &
         KERNELS1%U33(IRATN,J,I,IIM,IRE)*WAVEL
         IF(ANGLE(J) .LE. 40.)            &
         KERNELS1%U33(IRATN,J,I,IIM,IRE)=           &
         LOG( KERNELS1%U33(IRATN,J,I,IIM,IRE) )
      endif
      if(keyEL .gt. 4)                      &
         KERNELS1%U34(IRATN,J,I,IIM,IRE)=      &
         KERNELS1%U34(IRATN,J,I,IIM,IRE)*WAVEL
      if(keyEL .gt. 5) then
         KERNELS1%U44(IRATN,J,I,IIM,IRE)=  &
         KERNELS1%U44(IRATN,J,I,IIM,IRE)*WAVEL
      IF(ANGLE(J) .LE. 50.)            &
         KERNELS1%U44(IRATN,J,I,IIM,IRE)=           &
         LOG( KERNELS1%U44(IRATN,J,I,IIM,IRE) )
      endif
	  
      ENDDO ! IRATN
      ENDDO ! J
      ENDDO ! I
      ENDDO ! IIM
      ENDDO ! IRE
      ENDIF ! keyEL .GT. 0
      
      IF(key .eq. 1) THEN
!c      coeff=1.
!tl      coeff=1e+3                  ! for AERONET dV/dlnR
      dlnr=((LOG(grid(KN))-LOG(grid(1)))/(KN-1))  &
             /((LOG(XGRID1%grid1(KN1))-LOG(XGRID1%grid1(1)))/(KN1-1))
      do i1=1,KC
      SD(i1,1:KN)=SD(i1,1:KN)*dlnr  !*coeff !tl
      enddo ! i1
      CALL USUO_WL(keyEL                                           &
                  ,KN1,XGRID1%grid1,KM,ANGLE,WAVEL,KRE,KIM,ARE,AIM &
                  ,KN,grid,NWL,WL,KR,KC,SD                         &
                  ,xnmin,xnmax,xkmin,xkmax                         &
                  ,KERNELS1,KERNELS2)   
!c      write(*,*) 'After   USUO_WL'
!c *****************************************
      do i1=1,KC
      SD(i1,1:KN)=SD(i1,1:KN)/dlnr  !/coeff !tl
      enddo ! i1


!c START Save WL original kernels (UO... arrays)

      DO IRATN=1,NRATN
!c **
!c ** Write UO11, UO12, UO22, UO33, UO34, UO44 matrices
!c **

      comm_name1 = trim(KERNEL%comm_name(IRATN))//'_bin'//  &
           trim(adjustl(NBINC))//'_wl'//trim(adjustl(NWLC))//'_'             
!write(*,*) 'comm_name1 = ',comm_name1

      DO KEL=1,NNEL-1
      if(KEL .gt. keyEL) CYCLE
      name=trim(comm_name1)//NELC(KEL)//'.txt'
      full_name=TRIM(dir_name_N)//TRIM(name)
!cl        full_name=TRIM(name)
      OPEN(21,FILE=full_name,status='unknown')
      write(*,*) full_name
      WRITE(21,*) XGRID1%RATIO(iratn),'  RATIO'
      WRITE(21,*) KC,NWL,' number of lognorm.bins and wave lengths'
      WRITE(21,*) KM,'  number of angles'
      WRITE(21,27) ANGLE(1:KM)

      WRITE(*,*) 'WRITE matrix UO',NEL(KEL)
      WRITE(21,*) ARE(1),ARE(KRE),' real refr. indices'
!cl        WRITE(*,*) ARE(1),ARE(KRE),' RREMIN,RREMAX'   
      WRITE(21,*) -AIM(1),-AIM(KIM),' imag refr. indices'
!cl        WRITE(*,*) -AIM(1),-AIM(KIM),' RIMMIN,RIMMAX'        

      WRITE(21,*) KRE, -KIM,' number of intervals for opt. const'
!cl        WRITE(*,*) KRE, -KIM,' KRE, -KIM'

      DO  L=1,NWL
      DO  IRE=1,KRE
      DO  IIM=1,KIM
      if(NDPP .eq. 0) then
      if(XGRID1%XWL.ne.WAVEL) then
        write(tmp_message,'(2(a,es11.4))') 'in grid1.dat: XWL =',XGRID1%XWL, &
        '.NE. in Rke...: WAVEL =',WAVEL
        G_ERROR(trim(tmp_message))
      endif ! XGRID1%XWL
      NDPP=1
      endif ! NDPP
      
      WRITE(21,21) NEL(KEL),XGRID1%RATIO(iratn)
!cl        WRITE(*,*) NEL(KEL),XGRID1%RATIO(iratn),' KEL,RATIO'
      WRITE(21,20) WL(L),ARE(IRE),-AIM(IIM)

      IF(KEL .EQ. 1) then              
         DO  I=1,KC
         WRITE(21,27) KERNELS2%UO11(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ELSE IF(KEL .EQ. 2) THEN
         DO  I=1,KC
         WRITE(21,27) KERNELS2%UO12(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ELSE IF(KEL .EQ. 3) THEN
         DO  I=1,KC       
         WRITE(21,27) KERNELS2%UO22(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ELSE IF(KEL .EQ. 4) THEN
         DO  I=1,KC
         WRITE(21,27) KERNELS2%UO33(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ELSE IF(KEL .EQ. 5) THEN
         DO  I=1,KC
         WRITE(21,27) KERNELS2%UO34(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ELSE IF(KEL .EQ. 6) THEN
         DO  I=1,KC
         WRITE(21,27) KERNELS2%UO44(IRATN,1:KM,I,IIM,IRE,L)
         ENDDO ! I
      ENDIF ! KEL
      ENDDO ! IIM
      ENDDO ! IRE
      ENDDO ! L NWL
      write(*,*) 'KEL=',KEL
      CLOSE (21)
      ENDDO ! KEL	  	
      write(*,*) 'after read IRATN=',IRATN,  &
               ' RATIO=',XGRID1%RATIO(iratn)

!cl     &   ' U11=',U11(IRATN,1,1,1,1),
!cl     &   ' U11=',U11(IRATN,KM1,KN1,KIM,KRE)
!c **
!c ** WRITE UOEA MATRIX
!c **            
      name=trim(comm_name1)//NELC(0)//'.txt'
      full_name=TRIM(dir_name_N)//TRIM(name)
      OPEN(21,FILE=full_name,status='unknown')
      write(*,*) full_name

      WRITE(21,*) XGRID1%RATIO(iratn),'  RATIO'
      WRITE(21,*) KC,NWL,' number of lognorm.bins and wave lengths'
      WRITE(*,*) 'WRITE matrix UO',NEL(0)
      WRITE(*,*) 
      WRITE(21,*) ARE(1),ARE(KRE),' real refr. indices'
!cl        WRITE(*,*) ARE(1),ARE(KRE),' RREMIN,RREMAX'   
      WRITE(21,*) -AIM(1),-AIM(KIM),' imag refr. indices'
!cl        WRITE(*,*) -AIM(1),-AIM(KIM),' RIMMIN,RIMMAX'        
      WRITE(21,*) KRE, -KIM,' number of intervals for opt. const'
!cl        WRITE(*,*) KRE, KIM,' KRE, KIM'

      DO L=1,NWL
      DO IRE=1,KRE
      DO IIM=1,KIM      
         WRITE(21,21) nn,xx
!cl        WRITE(*,*) nn,xx,' KEL,RATIO'
         WRITE(21,20) WL(L),ARE(IRE),-AIM(IIM)    
         WRITE(21,*) ' EXTINCTION (1/km, for d()/dlnr m3/m3*km):'
         WRITE(21,27) (KERNELS2%UOEA(IRATN,1,I,IIM,IRE,L),I=1,KC)    
         WRITE(21,*) ' ABSORPTION (1/km, for dv/dlnr m3/m3*km):'
         WRITE(21,27) (KERNELS2%UOEA(IRATN,2,I,IIM,IRE,L),I=1,KC)
      ENDDO ! IIM
      ENDDO ! IRE
      ENDDO  ! L NWL

      CLOSE (11)
      ENDDO ! IRATN
      CLOSE(10)

   27 FORMAT (7E16.7)
   20 FORMAT (F14.5,2E16.7,'  wavel,rreal,rimag')
   21 FORMAT (I3,F14.5,'  element,ratn')


      NDP=1

      ENDIF ! NDP

      ENDIF ! key .eq. 1

      IF(key .eq. 1) THEN
!c
!c *** DEALLOCATE ARRAYS
!c
      DEALLOCATE(KERNELS1%UEA,stat=ierr)
      if(ierr/=0) stop 'Can not deallocate UEA array'
      if(keyEL .gt. 0) then
         DEALLOCATE(KERNELS1%U11,stat=ierr)
         if(ierr/=0) stop 'Can not deallocate U11 array'
      endif
      if(keyEL .gt. 1) then
        DEALLOCATE(KERNELS1%U12,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U12 array'
      endif
      if(keyEL .gt. 2) then
        DEALLOCATE(KERNELS1%U22,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U22 array'
      endif
      if(keyEL .gt. 3) then
        DEALLOCATE(KERNELS1%U33,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U33 array'
      endif
      if(keyEL .gt. 4) then
        DEALLOCATE(KERNELS1%U34,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U34 array'
      endif
      if(keyEL .gt. 5) then
        DEALLOCATE(KERNELS1%U44,stat=ierr)
        if(ierr/=0) stop 'Can not deallocate U44 array'
      endif 

      ENDIF ! key .eq. 1
                                       
      if(keySUB .eq. 0) then                                       
      write(6,*) 
      write(6,*)  & 
              '------------------ T I M I N G ------------------' 
      WRITE(*,62) T_CFM/60.     
      endif
62    format('  Calcul. fixed kernels ... ',f8.3,' min.')

      ENDIF ! key .NE. 2
    
!10    FORMAT(3E15.7,' rmin, rmax, RATIO')                    
!11    FORMAT(7E15.7)                    
!12    FORMAT(3E15.7,I4,'    wavelength, n, k, NEL')                    
!15    FORMAT(7F12.2)                    
       
      RETURN
      END SUBROUTINE MATRIX_FIX_bin

!c **************************************************************** 

