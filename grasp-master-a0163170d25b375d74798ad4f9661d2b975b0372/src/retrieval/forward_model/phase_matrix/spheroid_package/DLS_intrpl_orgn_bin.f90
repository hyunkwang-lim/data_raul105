! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../../../constants_set/mod_globals.inc"
      SUBROUTINE USMATRIX_bin(key,keyEL,keySUB,keyLS,key_RD1  &
                             ,key_RD                          &
                             ,NWL,WL,RN,RK,KN,grid            &
                             ,KR,R,RD,KM,ANGLE,dlnr,dlnr1     &
                             ,pomin,pomax                     &
                             ,KERNEL,XGRID1                   &
                             ,KC,SD,LB,LE                     &
                             ,xnmin,xnmax,xkmin,xkmax         &
                             ,KERNELS1,KERNELS2  )
                             
! Linear or Spline interpolation for key=0 02/09/2011
! 
! spline or linear : ext,      po, n, log(k)
! spline or linear : log(F11)...F44, po; F11…F44 n, log(k)

!c ** ANGLE & PO SPLINE interpolation version
!c ** version for Oleg
!c **
!c ** 12/04/03 f22 interpolation is logarithmic
!c ** 05/05/03 this version can be used to retrieve an aspect ratio
!c ** 13/10/03 IF(ANGLE<=40.)ln(f33)-interpol.
!c **          IF(ANGLE<=50.)ln(f44)-interpol.
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
!c **   keyEL=  0 - calculate EXT                                ** c
!c **           1 - EXT,F11                                      ** c
!c **           2 - EXT,F11,F12                                  ** c
!c **           3 - EXT,F11,F12,F22                              ** c
!c **           4 - EXT,F11,F12,F22,F33                          ** c
!c **           5 - EXT,F11,F12,F22,F33,F34                      ** c
!c **           6 - EXT,F11,F12,F22,F33,F34,F44                  ** c        
!c **   WL   - wavelength                                        ** c
!c **   RN   - real part of the refractive index                 ** c
!c **   RK   - imaginary part of the refractive index            ** c
!c **   rgmin,rgmax &                                            ** c 
!c **   wlmin,wlmax - min,max radii and wlmin,wlmax wavelengths  ** c
!c **                 that are used to recalculate grid radii for** c
!c **                 fixed kernels. New input file              ** c
!c **                'grid1.dat.new' will be created if key=1    ** c
!c **                 or key_org=1. Use key_grid1 to choose      ** c
!c **                 'grid1.dat' or 'grid1.dat.new' will be read** c
!c **                 for further calculations                   ** c  
!c **   xnmin,xnmax &                                            ** c
!c **   xkmin,xkmax - min,max parts of refractive index to be    ** c
!c **                 used to calculate WL original kernels      ** c
!c **                 (key=1)                                    ** c
!c **   key_SD=0 - read Size Distribution table dV/dlnR          ** c
!c **         =1 - calculate Size Distribution for grid radii    ** c
!c **              using Log Normal function                     ** c
!c **   ID    - dimension of d(...)/dlnR or d(...)/dR            ** c
!c **       = 0 - number                                         ** c
!c **       = 1 - radius                                         ** c
!c **       = 2 - area                                           ** c
!c **       = 3 - volume                                         ** c
!c **   NMD   - number of modes (up to 2)                        ** c
!c **   KN   - number of grid radii                              ** c
!c **   grid(KN) - grid radii                                    ** c
!c **   SD(KN)   - size distribution for grid radii              ** c
!c **   (CM(i),SM(i),RMM(i),i=1,NMD) - size distribution         ** c
!c **   function (LogNormal) parameters:                         ** c
!c **                         CM - concentration                 ** c
!c **                         SM - standard deviation            ** c
!c **                         RMM - median radius                ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_N - new original kernel directory               ** c
!c **                                      name (key_org=1)      ** c
!c **   KR  - number of axis ratios                              ** c
!c **   R(KR)  - grid axis ratios                                ** c
!c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
!c **   KM   - number of scattering angles                       ** c
!c **   ANGLE(KM) - scattering angles                            ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   ext     - extinction                                     ** c
!c **   albedo  - albedo                                         ** c
!c **   f... - scattering matrix elements                        ** c
!c **************************************************************** c
!c      use alloc1

      USE mod_par_DLS
      USE mod_par_DLS_bin
      use mod_intrpl_linear
      use mod_intrpl_spline
      use mod_type_DLS
      use mod_alloc_kernels
      use mod_stop_report

      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2

      integer KN, KR, KM,  &
              key, keyEL, keySUB, keyLS, key_RD1, key_RD
      real WAVEL
      real*4 dlnr, dlnr1, dlnr_grid, dlnr_grid1
      
      real WL(KWLpar),RN(KWLpar), RK(KWLpar)
      real*4 R(KRpar),RATIO(KR1par)

      real grid(KNpar),RD(KRpar)
      real grid1(KN1par), ANGLE(KMpar)
      real US11,US12,US22,US33,US34,US44,USEA
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
      real ARE(KREpar), AIM(KIMpar)
      type(DLS_kernel_location) :: KERNEL
      type(DLS_kernel_grid1)    :: XGRID1
      integer NWL,KC,LB,LE
      real SD(KCpar,KNpar)
!c      real LINEAR
!cl     &, LINEAR_LN
!c ----- Timer ------
      real*4 tarray(2),T_INT,T_INT0
      real*4 dtime
      real time_begin, time_end, T_CPU
      real time_begin1, time_end1, T_CPU1
      real*4 T_INT1,T_INT01

      save WAVEL,KN1,grid1,KRE,KIM,ARE,AIM,NRATN,RATIO &
                                 ,dlnr_grid,dlnr_grid1
      integer,save :: NDP=0
      real*4 xnmin,xnmax,xkmin,xkmax
!  ------------------------------------------------------------------------	
	  
      PI=ACOS(-1.)
      PI2=2.*PI

      IF(NDP.EQ.0) THEN 
      T_INT=0.
      T_CPU=0.
      T_INT1=0.
      T_CPU1=0.		
                           if(keySUB.eq.0) then
                              T_INT01=dtime(tarray) !+++
                              CALL CPU_TIME (time_begin1)
									         endif

      CALL MATRIX_FIX_bin(key,key_RD,keyEL,keySUB,    &
                  KR,                                 &
                  KN,grid,KM,ANGLE,                   &
                  WAVEL,KRE,KIM,ARE,AIM,pomin,pomax,  &
                  KERNEL,XGRID1,                      &
                  NWL,WL,KC,SD,                       &
                  xnmin,xnmax,xkmin,xkmax,            &
                  KERNELS1,KERNELS2)
      if(stop_report%status .eqv. .true.) return

                           if(keySUB.eq.0) then				                       
                              T_INT01=dtime(tarray)
                              T_INT1=T_INT1+T_INT01   !+++ 
                              CALL CPU_TIME (time_end1)
                              T_CPU1=T_CPU1+(time_end1-time_begin1)   !+++ 	  
!                                  WRITE(*,51) T_INT1/60.     
                              WRITE(*,51) T_INT1
                              WRITE(*,52) T_CPU1
                           endif
                           
   51 format('  Read kernels .................... ',f8.3,' sec')
   52 format('  Read kernels CPU_time ........... ',f8.3,' sec')

!c ** in order to adjust Size Distribution (SD) in OPTCHAR subroutine

!	          dlnr=((LOG(grid(KN))-LOG(grid(1)))/(KN-1))
!     &       /((LOG(grid1(KN1))-LOG(grid1(1)))/(KN1-1))
      KN1 = XGRID1%KN1
      grid1(1:KN1) = XGRID1%grid1(1:KN1)
		  dlnr_grid =(LOG(grid(KN))-LOG(grid(1)))/(KN-1)
		  dlnr_grid1=(LOG(grid1(KN1))-LOG(grid1(1)))/(KN1-1)
!      write(*,'(3(a,e14.5))') 'dlnr_grid=',dlnr_grid
!     &                      ,' dlnr_grid1=',dlnr_grid1
!     &                      ,' dlnr_grid/dlnr_grid1='
!     &                      ,dlnr_grid/dlnr_grid1

      NRATN = XGRID1%NRATN
      RATIO(1:NRATN) = XGRID1%RATIO(1:NRATN)
	  
      NDP=1
      ENDIF ! NDP 
	  
	   
		 
      dlnr = dlnr_grid
      dlnr1= dlnr_grid1

!cl	   write(*,*) 'after MATRIX_FIX_bin'
	   if(keySUB .eq. 0) then
	   write(*,*) 'KRE=',KRE,' KIM=',KIM
	   write(*,*) ARE(1),' <ARE< ',ARE(KRE)
	   write(*,*) AIM(1),' <AIM< ',AIM(KIM)
	   endif

!cl	  write(*,*) '3: U11(1,1,1,1,1): ',U11(1,1,1,1,1)
!c
!c ** CHECK OF SIZE PARAMETER (PO) RANGE
!c 
!cl       write(*,*) '********** in usmatrix.f: WAVEL=',WAVEL
									if(keySUB.eq.0) then
									   T_INT0=dtime(tarray) !+++
							           CALL CPU_TIME (time_begin)
									endif

      if(key.lt.2) then
      DO L=LB,LE
      POS =PI2*grid1(1)/WAVEL
      POF =PI2*grid1(KN1)/WAVEL
      PO1S=PI2*grid(1)/WL(L)
      PO1F=PI2*grid(KN)/WL(L)
      IF(PO1S.LT.POS.OR.PO1S.GT.POF) THEN
        WRITE(tmp_message,'(a,a,a,i0,3(3x,a,es11.4))') &
        'PO is out of look-up table in USMATRIX 1', &
        NEW_LINE('A'), &
        'L = ',L,'POS =',POS,'PO1S =',PO1S,'POF =',POF
        G_ERROR(trim(tmp_message))
      ENDIF
      IF(PO1F.LT.POS.OR.PO1F.GT.POF) THEN
        WRITE(tmp_message,'(a,a,a,i0,3(3x,a,es11.4))') &
        'PO is out of look-up table in USMATRIX 2', &
        NEW_LINE('A'), &
        'L = ',L,'POS =',POS,'PO1F=',PO1F,'POF=',POF
        G_ERROR(trim(tmp_message))
      ENDIF
      ENDDO ! NWL
      endif ! key<2

      if(key.eq.0) then
         if(key_RD1.eq.1) then
            call USU_LS_bin(key_RD,keyEL,keySUB,keyLS,KR,R,RD,RATIO,NRATN  &
                  ,KN1,grid1,KM,ANGLE,WAVEL,KRE,KIM,ARE,AIM          &
                  ,KN,grid,NWL,WL,RN,RK                              &
                  ,KERNELS1) 
         else

            call USU_LS_RD_bin(key_RD,keyEL,keySUB,keyLS      &
                  ,KR,R,RD,RATIO                        &
                  ,KN1,grid1,KM,ANGLE,WAVEL,KRE,KIM,ARE,AIM   &
                  ,KN,grid,NWL,WL,RN,RK,KERNELS1) 
         endif ! key_RD1
      else ! key > 0
         if(key_RD1.eq.1) then
            call USU_WL(key_RD,keyEL,keySUB,KR,R,RD,RATIO,NRATN  &
                   ,KM,KRE,KIM,ARE,AIM     &
                   ,RN,RK,KC,LB,LE,KERNELS2)
         else		
            call USU_WL_RD(key_RD,keyEL,keySUB,KR,R,RD,RATIO  &
                   ,KM,KRE,KIM,ARE,AIM        &
                   ,RN,RK,KC,LB,LE,KERNELS2)
            if(keySUB.eq.0) write(*,*) 'key_RD1=',key_RD1
         endif       
      endif ! key
                        if(keySUB.eq.0) then
                           T_INT0=dtime(tarray)
                           T_INT=T_INT+T_INT0   !+++ 
                           CALL CPU_TIME (time_end)
                           T_CPU=T_CPU+(time_end-time_begin)   !+++ 	  
!                              WRITE(*,61) T_INT/60.     
                           WRITE(*,61) T_INT
                           WRITE(*,62) T_CPU
                        endif

   61 format('  Interpolation .................... ',f8.3,' sec')
   62 format('  Interpolation CPU_time ........... ',f8.3,' sec')
!      write(*,*) 'test DLS_intrpl_orgn_bin.f US22(1,1:KN,1) :'
!      write(*,'(5e14.4)') US22(1,1:KN,1)

      RETURN
      END SUBROUTINE USMATRIX_bin


!c **************************************************************** 
      SUBROUTINE USU_LS_bin(key_RD,keyEL,keySUB,keyLS,KR,R,RD,RATIO  &
                        ,NRATN,KN1,grid1,KM,ANGLE,WAVEL              &
                        ,KRE,KIM,ARE,AIM,KN,grid,NWL,WL,RN,RK        &
                        ,KERNELS1)       
      !use alloc1
      USE mod_alloc_kernels
      USE mod_par_DLS
      USE mod_par_DLS_bin
      use mod_intrpl_linear
      use mod_intrpl_spline
	  
      type(kernels_triangle_bin),  intent(in)  ::  KERNELS1
      
      real*4 R(KRpar), RATIO(KR1par), RRATN
      dimension grid(KNpar),RD(KRpar)
      dimension grid1(KN1par), ANGLE(KMpar)
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
      dimension RDc(KR),AA(6),BB(6),AB(6)  &
                       ,CC(6),DD(6),CD(6)  &
                       ,sumUS(6)
      dimension XPO(2), YPO(2), X(2), Y(2) &
               ,ARE(KREpar), AIM(KIMpar)
      integer NWL, keyEL, keySUB, keyLS
      real WL(KWLpar), RN(KWLpar), RK(KWLpar)
!c **  for SPLINE subroutine
      DOUBLE PRECISION :: XARG, YFIT
      DOUBLE PRECISION, DIMENSION(KN1par) :: XXS1, YYS1, b, c, d

      if(NWL.gt.1) then
         write(*,*) 'STOP NWL.gt.1 in USU subroutine'
         STOP
      endif

      write(*,*) 'DLS_intrpl_orgn_bin.f !!! keyLS=',keyLS
	  
      PI=ACOS(-1.)
      PI2=2.*PI
!cl      cinv=WAVEL/WL(1)  ! invariant for dV(..)/dlnr  
       IF(key_RD.eq.2) then
!c 
!c*** RECALCULATE ASPECT RATIO DISTRIBUTION (RDc()=SAREA/VOLUME)
!c*** RDc()=RD()/RDc(); sumRD=sum(RDc())
!c
!c ** OBLATE
        do IR=1,KR
        if(R(IR).lt.1.) then              
        E=SQRT( 1.-R(IR)*R(IR) )
        xa1=LOG((1.+E)/(1.-E))/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+  &
                  0.5*xa1*R(IR)**(4./3.))
!c ** PROLATE            
        elseif(R(IR).gt.1.) then          
        E=SQRT( 1.-1./R(IR)/R(IR) )
        xa2=ASIN(E)/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+  &
                     xa2*R(IR)**(1./3.))
!c ** SPHERE
        elseif(R(IR).eq.1.) then             
        RDc(IR)=3.
        endif ! R()
!c ** WRITE ASPECT RATIO DISTRIBUTION 
!c          write(*,*) 'R=',R(IR),' B=',RDc(IR),
!c     &  ' 1/B=',1./RDc(IR),' RD=',RD(IR)

        enddo ! IR
        RDc(:KR)=RD(:KR)/RDc(:KR)
      ENDIF ! key_RD
      IF(key_RD.eq.1) RDc(:KR)=RD(:KR)

      sumRD=sum(RDc(:KR))

      if(keySUB.eq.0) then
      do IR=1,KR 
      write(*,*) 'R=',R(IR),' RDc=',RDc(IR)
      enddo ! IR
      write(*,*) 'sumRD=',sumRD
      endif

      RL=RN(1)
      RI=RK(1)
      I0=0
      I1=0
      K0=0
      K1=0
      DO I=1,KRE-1
      IF(RL.GE.ARE(I).AND.RL.LE.ARE(I+1)) THEN
       I0=I
       I1=I+1
      ENDIF
      ENDDO
      DO I=1,KIM-1
      IF(RI.GE.AIM(I).AND.RI.LE.AIM(I+1)) THEN
        K0=I
        K1=I+1
      ENDIF
      ENDDO
      IF(RL.LE.ARE(1)) THEN
        if(keySUB.eq.0) then
        IF(RL.LT.ARE(1)) THEN
        WRITE(*,*) 'n=',RN,' is out of the range:',ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(1)
        ENDIF
        endif
        I0=1
        I1=2
        RN(1)=ARE(1)
        RL=ARE(1)
      ENDIF
      IF(RL.GE.ARE(KRE)) THEN
        if(keySUB.eq.0) then
        IF(RL.GT.ARE(KRE)) THEN
        WRITE(*,*) 'n=',RN,' is out of the range:',ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(KRE)
        ENDIF
        endif
        I0=KRE-1
        I1=KRE
        RN(1)=ARE(KRE)
        RL=ARE(KRE)
      ENDIF
      IF(RI.LE.AIM(1)) THEN
        if(keySUB.eq.0) then
        IF(RI.LT.AIM(1)) THEN
        WRITE(*,*) 'k=',RK,' is out of the range:',AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(1)
        ENDIF
        endif
        K0=1
        K1=2
        RK(1)=AIM(1)
        RI=AIM(1)
      ENDIF
      IF(RI.GE.AIM(KIM)) THEN
        if(keySUB.eq.0) then
        IF(RI.GT.AIM(KIM)) THEN
        WRITE(*,*) 'k=',RK,' is out of the range:',AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(KIM)
        ENDIF
		endif
		
		K0=KIM-1
        K1=KIM
        RK(1)=AIM(KIM)
        RI=AIM(KIM)
      ENDIF
!C      write(*,*) I0,I1,K0,K1,' I0,I1,K0,K1'
!C      WRITE(*,*) ARE(I0),ARE(I1),' ARE'
!C      WRITE(*,*) AIM(K0),AIM(K1),' AIM'
      DO I=1,KN
      IP0=0
      IP1=0
      PO=grid(I)*PI2/WL(1)            

      DO IP=1,KN1-1
      PO1=grid1(IP)  *PI2/WAVEL   
      PO2=grid1(IP+1)*PI2/WAVEL 
      IF(PO.GE.PO1.AND.PO.LE.PO2) THEN
      IP0=IP
      IP1=IP+1
      ENDIF
      IF(PO.LE.(grid1(1)*PI2/WAVEL)) THEN
      IP0=1
      IP1=2
      ENDIF
      IF(PO.GE.(grid1(KN1)*PI2/WAVEL)) THEN
      IP0=KN1-1
      IP1=KN1
      ENDIF      
      ENDDO ! IP

      IF(keyLS.EQ.1) THEN
        XPO(1)=grid1(IP0)*PI2/WAVEL
        XPO(2)=grid1(IP1)*PI2/WAVEL	
      ELSE
	    XXS1(1:KN1)=grid1(1:KN1)*PI2/WAVEL
	    XARG=PO
	  ENDIF

      IF(keyEL .EQ. 0) GOTO 100
!c
!c **  SCATTERING MATRIX ELEMENTS
!c
      DO J=1,KM
      sumUS(:6)=0.
      DO IR=1,KR
      IF(NRATN.EQ.1) then
        RRATN=RATIO(1)
        IF(RRATN.NE.R(IR)) THEN
        WRITE(*,*) 'R=',R(IR),' .NE. RATIO=',RATIO(1)
        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(1)
        R(IR)=RRATN
        ENDIF
      ELSE
        RRATN=R(IR)
      ENDIF
!c
!c ** AXIAL RATIO LOOP
!c
      IF(NRATN.NE.1) THEN
      DO IRATN=1,NRATN-1
      IF(RRATN.GE.RATIO(IRATN).AND.RRATN.LE.RATIO(IRATN+1)) THEN
        L0=IRATN
        L1=IRATN+1
      ENDIF
      ENDDO
      IF(RRATN.LE.RATIO(1)) THEN
        IF(RRATN.LT.RATIO(1)) THEN
        WRITE(*,*) 'R=',R(IR),' is out of the range:',RATIO(1),'< R <',RATIO(NRATN)
        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(1)
        ENDIF
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      IF(RRATN.GE.RATIO(NRATN)) THEN
        IF(RRATN.GT.RATIO(NRATN)) THEN
        WRITE(*,*) 'R=',R(IR),' is out of the range:',RATIO(1),'< R <',RATIO(NRATN)
        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(NRATN)
        ENDIF
        L0=NRATN-1
        L1=NRATN
        R(IR)=RATIO(NRATN)
        RRATN=RATIO(NRATN)
      ENDIF
      ELSE
        L0=1
        L1=1
      ENDIF

      X(1)=ARE(I0)
      X(2)=ARE(I1)
!c
!c ** U11     AA(1)
!c
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U11(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U11(L0,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U11(L0,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	    ELSE
	      YYS1(1:KN1)=KERNELS1%U11(L0,J,1:KN1,K0,I0)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	      YYS1(1:KN1)=KERNELS1%U11(L0,J,1:KN1,K0,I1)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	    ENDIF
      AA(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     AA(2)
!c
      if(keyEL .gt. 1) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U12(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U12(L0,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U12(L0,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	    ELSE
	      YYS1(1:KN1)=KERNELS1%U12(L0,J,1:KN1,K0,I0)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	      YYS1(1:KN1)=KERNELS1%U12(L0,J,1:KN1,K0,I1)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	    ENDIF
      AA(2)=LINEAR( X, Y, 2, RL )
	  endif
!c
!c ** U22     AA(3)
!c
      if(keyEL.gt.2) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U22(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U22(L0,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U22(L0,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	    ELSE
	      YYS1(1:KN1)=KERNELS1%U22(L0,J,1:KN1,K0,I0)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
        YYS1(1:KN1)=KERNELS1%U22(L0,J,1:KN1,K0,I1)
        CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
      ENDIF
      AA(3)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 1
!c
!c ** U33     AA(4)
!c
      if(keyEL .gt. 3) then
	  IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U33(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U33(L0,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U33(L0,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(L0,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(L0,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        AA(4)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 3	
!c
!c ** U34     AA(5)
!c
      if(keyEL .gt. 4) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U34(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U34(L0,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U34(L0,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U34(L0,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(L0,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        AA(5)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 4
!c
!c ** U44     AA(6)
!c
      if(keyEL .gt. 5) then
      IF(keyLS .EQ. 1) THEN
	      YPO(1)=KERNELS1%U44(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U44(L0,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U44(L0,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
 	    YYS1(1:KN1)=KERNELS1%U44(L0,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(L0,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))   
        ENDIF ! ANGLE
        AA(6)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 5
!c
!c ** U11     BB(1)
!c
      IF(keyLS .eq. 1) THEN
        YPO(1)=KERNELS1%U11(L0,J,IP0,K1,I0)
		YPO(2)=KERNELS1%U11(L0,J,IP1,K1,I0) 
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U11(L0,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U11(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U11(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        BB(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     BB(2)
!c
      if(keyEL.gt.1) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U12(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U12(L0,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U12(L0,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U12(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U12(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
	    Y(2)=YFIT
	  ENDIF
        BB(2)=LINEAR( X, Y, 2, RL )
	  endif
!c
!c ** U22     BB(3)
!c
      if(keyEL.gt.2) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U22(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U22(L0,J,IP1,K1,I0) 
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U22(L0,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U22(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U22(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        BB(3)=LINEAR( X, Y, 2, RL )
	  endif
!c
!c ** U33     BB(4)
!c
      if(keyEL.gt.3) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U33(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U33(L0,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U33(L0,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        BB(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     BB(5)
!c
      if(keyEL.gt.4) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U34(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U34(L0,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U34(L0,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
      ELSE
	    YYS1(1:KN1)=KERNELS1%U34(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        BB(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     BB(6)
!c
      if(keyEL.gt.5) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U44(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U44(L0,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U44(L0,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U44(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        BB(6)=LINEAR( X, Y, 2, RL )
      endif

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))

      DO II=1,keyEL
      Y(1)=AA(II)
      Y(2)=BB(II)
      AB(II)=LINEAR( X, Y, 2, log(RI) )
      ENDDO ! II
	  
      IF(NRATN.NE.1) THEN
      X(1)=ARE(I0)
      X(2)=ARE(I1)
!c
!c ** U11     CC(1)
!c
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U11(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U11(L1,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U11(L1,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
      ELSE
	    YYS1(1:KN1)=KERNELS1%U11(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U11(L1,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        CC(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     CC(2)
!c
      if(keyEL.gt.1) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U12(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U12(L1,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U12(L1,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U12(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U12(L1,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        CC(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     CC(3)
!c
      if(keyEL.gt.2) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U22(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U22(L1,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U22(L1,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U22(L1,J,1:KN1,K0,I0)
	    key_spln=0
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U22(L1,J,1:KN1,K0,I1)
	    key_spln=0
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        CC(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     CC(4)
!c
      if(keyEL.gt.3) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U33(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U33(L1,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U33(L1,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(L1,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        CC(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     CC(5)
!c
      if(keyEL.gt.4) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U34(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U34(L1,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U34(L1,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U34(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(L1,J,1:KN1,K0,I1)
	    key_spln=0
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        CC(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     CC(6)
!c
      if(keyEL.gt.5) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U44(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U44(L1,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U44(L1,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U44(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(L1,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        CC(6)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U11     DD(1)
!c
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U11(L1,J,IP0,K1,I0)
		YPO(2)=KERNELS1%U11(L1,J,IP1,K1,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U11(L1,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U11(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U11(L1,J,1:KN1,K1,I1)
	    key_spln=0
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        DD(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     DD(2)
!c
      if(keyEL.gt.1) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U12(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U12(L1,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U12(L1,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U12(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U12(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        DD(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     DD(3)
!c
      if(keyEL.gt.2) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U22(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U22(L1,J,IP1,K1,I0) 
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U22(L1,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U22(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U22(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        DD(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     DD(4)
!c
      if(keyEL.gt.3) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U33(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U33(L1,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U33(L1,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
      ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        DD(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     DD(5)
!c
      if(keyEL.gt.4) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U34(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U34(L1,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U34(L1,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
        YYS1(1:KN1)=KERNELS1%U34(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
      ENDIF
        DD(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     DD(6)
!c
      if(keyEL.gt.5) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U44(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U44(L1,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U44(L1,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U44(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))    
        ENDIF ! ANGLE
        DD(6)=LINEAR( X, Y, 2, RL )
      endif

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      DO II=1,keyEL
      Y(1)=CC(II)
      Y(2)=DD(II)
      CD(II)=LINEAR( X, Y, 2, log(RI) )
      ENDDO ! II
	  
      RRATN1=RRATN

      X(1)=RATIO(L0)
      X(2)=RATIO(L1)

!c ** US11
      Y(1)=AB(1)
      Y(2)=CD(1)
      sumUS(1)=sumUS(1)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
!c ** US12
      if(keyEL.gt.1) then
      Y(1)=AB(2)
      Y(2)=CD(2)
      sumUS(2)=sumUS(2)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      endif
!c ** US22
      if(keyEL.gt.2) then
      Y(1)=AB(3)
      Y(2)=CD(3)
      sumUS(3)=sumUS(3)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      endif
!c ** US33
      if(keyEL.gt.3) then    
      Y(1)=AB(4)
      Y(2)=CD(4)
      sumUS(4)=sumUS(4)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      endif
!c ** US34
      if(keyEL.gt.4) then
      Y(1)=AB(5)
      Y(2)=CD(5)
      sumUS(5)=sumUS(5)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      endif
!c ** US44
      if(keyEL.gt.5) then
      Y(1)=AB(6)
      Y(2)=CD(6)
      sumUS(6)=sumUS(6)+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      endif 
	  
      ELSEIF(NRATN.EQ.1) THEN

      sumUS(1)=AB(1)*RDc(IR)
	  if(keyEL.gt.1) sumUS(2)=AB(2)*RDc(IR)        		
	  if(keyEL.gt.2) sumUS(3)=AB(3)*RDc(IR)
	  if(keyEL.gt.3) sumUS(4)=AB(4)*RDc(IR)
	  if(keyEL.gt.4) sumUS(5)=AB(5)*RDc(IR)
	  if(keyEL.gt.5) sumUS(6)=AB(6)*RDc(IR)
      ENDIF ! NRATN

      ENDDO ! IR RD

      US11(J,I,1)=sumUS(1)/WL(1)/sumRD
      if(keyEL.gt.1) US12(J,I,1)=sumUS(2)/WL(1)/sumRD
      if(keyEL.gt.2) US22(J,I,1)=sumUS(3)/WL(1)/sumRD
      if(keyEL.gt.3) US33(J,I,1)=sumUS(4)/WL(1)/sumRD
      if(keyEL.gt.4) US34(J,I,1)=sumUS(5)/WL(1)/sumRD
      if(keyEL.gt.5) US44(J,I,1)=sumUS(6)/WL(1)/sumRD

      ENDDO   ! J KM

100   CONTINUE

!c
!c ** EXTINCTION & ABSORPTION
!c
      DO J=1,2
      sumUSEA=0.
      DO IR=1,KR
      IF(NRATN.EQ.1) then
      RRATN=RATIO(NRATN)
      ELSE
      RRATN=R(IR)
      ENDIF
      IF(NRATN .NE. 1) THEN
      DO IRATN=1,NRATN-1
      IF(RRATN .GE. RATIO(IRATN) .AND. RRATN .LE. RATIO(IRATN+1)) THEN
          L0=IRATN
          L1=IRATN+1
!cl      write(*,*) IRATN,RATIO(IRATN),L0,L1,R(IR),RRATN, 
!cl     & ' IRATN,RATIO(IRATN),L0,L1 R(IR),RRATN',
!cl     & ' for ext & abs'
      ENDIF
      ENDDO

      IF(RRATN .GE. RATIO(NRATN)) THEN
        L0=NRATN-1
        L1=NRATN
        R(IR)=RATIO(NRATN)
        RRATN=RATIO(NRATN)
      ENDIF
      IF(RRATN .LE. RATIO(1)) THEN
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      ELSE
        L0=1
        L1=1
      ENDIF

!cl	XXS1(1:KN1)=grid1(1:KN1)*PI2/WAVEL
!cl	XARG=PO

      RRATN1=RRATN

      X(1)=ARE(I0)
      X(2)=ARE(I1)

      IF(keyLS .EQ. 1) THEN      
        YPO(1)=KERNELS1%UEA(L0,J,IP0,K0,I0)
        YPO(2)=KERNELS1%UEA(L0,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(L0,J,IP0,K0,I1)
        YPO(2)=KERNELS1%UEA(L0,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%UEA(L0,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(L0,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        AA(1)=LINEAR( X, Y, 2, RL )

      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%UEA(L0,J,IP0,K1,I0)
        YPO(2)=KERNELS1%UEA(L0,J,IP1,K1,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(L0,J,IP0,K1,I1)
        YPO(2)=KERNELS1%UEA(L0,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%UEA(L0,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(L0,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        BB(1)=LINEAR( X, Y, 2, RL )      

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      Y(1)=AA(1)
      Y(2)=BB(1)
      AB(1)=LINEAR( X, Y, 2, log(RI) )

      IF(NRATN .NE. 1) THEN    
      X(1)=ARE(I0)
      X(2)=ARE(I1)

      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%UEA(L1,J,IP0,K0,I0)
        YPO(2)=KERNELS1%UEA(L1,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(L1,J,IP0,K0,I1)
        YPO(2)=KERNELS1%UEA(L1,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
        YYS1(1:KN1)=KERNELS1%UEA(L1,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(L1,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        CC(1)=LINEAR( X, Y, 2, RL )

      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%UEA(L1,J,IP0,K1,I0)
        YPO(2)=KERNELS1%UEA(L1,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(L1,J,IP0,K1,I1)
        YPO(2)=KERNELS1%UEA(L1,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%UEA(L1,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(L1,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        DD(1)=LINEAR( X, Y, 2, RL )     

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      Y(1)=CC(1)
      Y(2)=DD(1)
      CD(1)=LINEAR( X, Y, 2, log(RI) )

      X(1)=RATIO(L0)
      X(2)=RATIO(L1)
      Y(1)=AB(1)
      Y(2)=CD(1)
      sumUSEA=sumUSEA+LINEAR( X, Y, 2, RRATN1 )*RDc(IR)
      ELSEIF(NRATN .EQ. 1) THEN
      sumUSEA=AB(1)*RDc(IR)      
      ENDIF
      ENDDO ! IR KR

      USEA(J,I,1)=sumUSEA*WAVEL/WL(1)/sumRD     
      
      ENDDO ! J 2

      ENDDO ! I KN           
      
	RETURN 
	END SUBROUTINE USU_LS_bin

!***************************
!***************************

      SUBROUTINE USU_LS_RD_bin(key_RD,keyEL,keySUB,keyLS       &
                        ,KR,R,RD,RATIO                         &
                        ,KN1,grid1,KM,ANGLE,WAVEL              &
                        ,KRE,KIM,ARE,AIM,KN,grid,NWL,WL,RN,RK  &
                        ,KERNELS1) 

      !use alloc1
      USE mod_par_DLS
      USE mod_par_DLS_bin
      use mod_intrpl_linear
      use mod_intrpl_spline
      use mod_alloc_kernels
	  
      type(kernels_triangle_bin),  intent(in)  ::  KERNELS1

      real*4 R(KRpar), RATIO(KR1par)
      real grid(KNpar),RD(KRpar)
      real grid1(KN1par), ANGLE(KMpar)
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
      real RDc(KR),AA(6),BB(6),AB(6),sumUS(6)
      real XPO(2), YPO(2), X(2), Y(2) &
               ,ARE(KREpar), AIM(KIMpar)
      integer NWL, keyEL, keySUB, keyLS
      real WL(KWLpar), RN(KWLpar), RK(KWLpar)
!c **  for SPLINE subroutine

      DOUBLE PRECISION :: XARG, YFIT
      DOUBLE PRECISION, DIMENSION(KN1par) :: XXS1, YYS1, b, c, d

      if(NWL.gt.1) then
         write(*,*) 'STOP NWL.gt.1 in USU subroutine'
         write(*,*) 'in DLS_intrpl_orgn_bin.f'
         STOP
      endif
	  
      PI=ACOS(-1.)
      PI2=2.*PI
!cl      cinv=WAVEL/WL(1)  ! invariant for dV(..)/dlnr  
      IF(key_RD.eq.2) then
!c 
!c*** RECALCULATE ASPECT RATIO DISTRIBUTION (RDc()=SAREA/VOLUME)
!c*** RDc()=RD()/RDc(); sumRD=sum(RDc())
!c
!c ** OBLATE
        do IR=1,KR
        if(R(IR).lt.1.) then              
        E=SQRT( 1.-R(IR)*R(IR) )
        xa1=LOG((1.+E)/(1.-E))/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+  &
                  0.5*xa1*R(IR)**(4./3.))
!c ** PROLATE            
        elseif(R(IR).gt.1.) then          
        E=SQRT( 1.-1./R(IR)/R(IR) )
        xa2=ASIN(E)/E
        RDc(IR)=1.5*(R(IR)**(-2./3.)+  &
                     xa2*R(IR)**(1./3.))
!c ** SPHERE
        elseif(R(IR).eq.1.) then             
        RDc(IR)=3.
        endif ! R()
!c ** WRITE ASPECT RATIO DISTRIBUTION 
!c          write(*,*) 'R=',R(IR),' B=',RDc(IR),
!c     &  ' 1/B=',1./RDc(IR),' RD=',RD(IR)

        enddo ! IR
        RDc(:KR)=RD(:KR)/RDc(:KR)
      ENDIF ! key_RD
      IF(key_RD.eq.1) RDc(:KR)=RD(:KR)

      sumRD=sum(RDc(:KR))

      if(keySUB.eq.0) then
      do IR=1,KR 
      write(*,*) 'R=',R(IR),' RDc=',RDc(IR)
      enddo ! IR
      write(*,*) 'sumRD=',sumRD
      endif

      RL=RN(1)
      RI=RK(1)
      I0=0
      I1=0
      K0=0
      K1=0
      DO I=1,KRE-1
      IF(RL.GE.ARE(I).AND.RL.LE.ARE(I+1)) THEN
       I0=I
       I1=I+1
      ENDIF
      ENDDO
      DO I=1,KIM-1
      IF(RI.GE.AIM(I).AND.RI.LE.AIM(I+1)) THEN
        K0=I
        K1=I+1
      ENDIF
      ENDDO
      IF(RL.LE.ARE(1)) THEN
        if(keySUB.eq.0) then
        IF(RL.LT.ARE(1)) THEN
        WRITE(*,*) 'n=',RN,' is out of the range:',ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(1)
        ENDIF
        endif
        I0=1
        I1=2
        RN(1)=ARE(1)
        RL=ARE(1)
      ENDIF
      IF(RL.GE.ARE(KRE)) THEN
        if(keySUB.eq.0) then
        IF(RL.GT.ARE(KRE)) THEN
        WRITE(*,*) 'n=',RN,' is out of the range:',ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(KRE)
        ENDIF
        endif
        I0=KRE-1
        I1=KRE
        RN(1)=ARE(KRE)
        RL=ARE(KRE)
      ENDIF
      IF(RI.LE.AIM(1)) THEN
        if(keySUB.eq.0) then
        IF(RI.LT.AIM(1)) THEN
        WRITE(*,*) 'k=',RK,' is out of the range:',AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(1)
        ENDIF
        endif
        K0=1
        K1=2
        RK(1)=AIM(1)
        RI=AIM(1)
      ENDIF
      IF(RI.GE.AIM(KIM)) THEN
        if(keySUB.eq.0) then
        IF(RI.GT.AIM(KIM)) THEN
        WRITE(*,*) 'k=',RK,' is out of the range:',AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(KIM)
        ENDIF
		endif
		
		K0=KIM-1
        K1=KIM
        RK(1)=AIM(KIM)
        RI=AIM(KIM)
      ENDIF
!C      write(*,*) I0,I1,K0,K1,' I0,I1,K0,K1'
!C      WRITE(*,*) ARE(I0),ARE(I1),' ARE'
!C      WRITE(*,*) AIM(K0),AIM(K1),' AIM'
      DO I=1,KN
      IP0=0
      IP1=0
      PO=grid(I)*PI2/WL(1)            

      DO IP=1,KN1-1
      PO1=grid1(IP)  *PI2/WAVEL   
      PO2=grid1(IP+1)*PI2/WAVEL 
      IF(PO.GE.PO1.AND.PO.LE.PO2) THEN
      IP0=IP
      IP1=IP+1
      ENDIF
      IF(PO.LE.(grid1(1)*PI2/WAVEL)) THEN
      IP0=1
      IP1=2
      ENDIF
      IF(PO.GE.(grid1(KN1)*PI2/WAVEL)) THEN
      IP0=KN1-1
      IP1=KN1
      ENDIF      
      ENDDO ! IP

      IF(keyLS.EQ.1) THEN
        XPO(1)=grid1(IP0)*PI2/WAVEL
        XPO(2)=grid1(IP1)*PI2/WAVEL	
      ELSE
	    XXS1(1:KN1)=grid1(1:KN1)*PI2/WAVEL
	    XARG=PO
	  ENDIF

      IF(keyEL .EQ. 0) GOTO 100
!c
!c **  SCATTERING MATRIX ELEMENTS
!c
      DO J=1,KM        ! angle loop
      sumUS(:6)=0.
      DO IR=1,KR       ! axis ratio loop
        IF(RATIO(IR).NE.R(IR)) THEN
        WRITE(*,*) 'key_RD1=2 and R=',R(IR),  &
      ' .NE. RATIO=',RATIO(IR)
        WRITE(*,*) 'STOP in matrix_intrpl_...f'
        STOP
        ENDIF
		
if(I0.eq.0.or.I1.eq.0) write(*,*) I0,I1,'  -I0,I1'
      X(1)=ARE(I0)
      X(2)=ARE(I1)
!c
!c ** U11     AA(1)
!c
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U11(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U11(IR,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U11(IR,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U11(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U11(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        AA(1)=LINEAR( X, Y, 2, RL )

!c
!c ** U12     AA(2)
!c
      if(keyEL .gt. 1) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U12(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U12(IR,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U12(IR,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U12(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U12(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        AA(2)=LINEAR( X, Y, 2, RL )
	  endif ! keyEL .gt. 1
!c
!c ** U22     AA(3)
!c
      if(keyEL .gt. 2) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U22(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U22(IR,J,IP1,K0,I0)
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U22(IR,J,IP1,K0,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U22(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U22(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
       AA(3)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 2
!c
!c ** U33     AA(4)
!c
      if(keyEL .gt. 3) then
	  IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U33(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U33(IR,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U33(IR,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        AA(4)=LINEAR( X, Y, 2, RL )
      endif	! keyEL .gt. 3
!c
!c ** U34     AA(5)
!c
      if(keyEL .gt. 4) then
      IF(keyLS.EQ.1) THEN
        YPO(1)=KERNELS1%U34(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U34(IR,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U34(IR,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U34(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
       Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        AA(5)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 4
!c
!c ** U44     AA(6)
!c
      if(keyEL .gt. 5) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U44(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%U44(IR,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%U44(IR,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
 	    YYS1(1:KN1)=KERNELS1%U44(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
       Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))   
        ENDIF ! ANGLE
        AA(6)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 5
!c
!c ** U11     BB(1)
!c
      IF(keyLS .eq. 1) THEN
        YPO(1)=KERNELS1%U11(IR,J,IP0,K1,I0)
		  YPO(2)=KERNELS1%U11(IR,J,IP1,K1,I0) 
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U11(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U11(IR,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U11(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
       Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U11(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
       Y(2)=EXP(YFIT)
	  ENDIF
        BB(1)=LINEAR( X, Y, 2, RL )

!c
!c ** U12     BB(2)
!c
      if(keyEL .gt. 1) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U12(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U12(IR,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U12(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U12(IR,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U12(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
       Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U12(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
	    Y(2)=YFIT
	  ENDIF
       BB(2)=LINEAR( X, Y, 2, RL )
	  endif ! keyEL .gt. 1
!c
!c ** U22     BB(3)
!c
      if(keyEL .gt. 2) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U22(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U22(IR,J,IP1,K1,I0) 
        Y(1)=EXP(LINEAR( XPO, YPO, 2, PO ))
        YPO(1)=KERNELS1%U22(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U22(IR,J,IP1,K1,I1)
        Y(2)=EXP(LINEAR( XPO, YPO, 2, PO ))
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U22(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=EXP(YFIT)
	    YYS1(1:KN1)=KERNELS1%U22(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=EXP(YFIT)
	  ENDIF
        BB(3)=LINEAR( X, Y, 2, RL )
	  endif ! keyEL .gt. 2
!c
!c ** U33     BB(4)
!c
      if(keyEL .gt. 3) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U33(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U33(IR,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U33(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U33(IR,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U33(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U33(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.40.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        BB(4)=LINEAR( X, Y, 2, RL )
      endif ! (keyEL .gt. 3
!c
!c ** U34     BB(5)
!c
      if(keyEL .gt. 4) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U34(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U34(IR,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U34(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U34(IR,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
      ELSE
	    YYS1(1:KN1)=KERNELS1%U34(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U34(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        BB(5)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 4
!c
!c ** U44     BB(6)
!c
      if(keyEL .gt. 5) then
      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%U44(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%U44(IR,J,IP1,K1,I0) 
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%U44(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%U44(IR,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%U44(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%U44(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        IF(ANGLE(J).LE.50.) THEN
        Y(1)=EXP(Y(1))  
        Y(2)=EXP(Y(2))  
        ENDIF ! ANGLE
        BB(6)=LINEAR( X, Y, 2, RL )
      endif ! keyEL .gt. 5

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))

      DO II=1,keyEL
      Y(1)=AA(II)
      Y(2)=BB(II)
      AB(II)=LINEAR( X, Y, 2, log(RI) )
      sumUS(II)=sumUS(II)+AB(II)*RDc(IR)
      ENDDO ! II

      ENDDO ! IR RD	  

                       US11(J,I,1)=sumUS(1)/WL(1)/sumRD
      if(keyEL .gt. 1) US12(J,I,1)=sumUS(2)/WL(1)/sumRD
      if(keyEL .gt. 2) US22(J,I,1)=sumUS(3)/WL(1)/sumRD
      if(keyEL .gt. 3) US33(J,I,1)=sumUS(4)/WL(1)/sumRD
      if(keyEL .gt. 4) US34(J,I,1)=sumUS(5)/WL(1)/sumRD
      if(keyEL .gt. 5) US44(J,I,1)=sumUS(6)/WL(1)/sumRD

      ENDDO   ! J KM

100   CONTINUE

!c
!c ** EXTINCTION & ABSORPTION
!c
      DO J=1,2
      sumUSEA=0.
      DO IR=1,KR

      X(1)=ARE(I0)
      X(2)=ARE(I1)

      IF(keyLS .EQ. 1) THEN      
        YPO(1)=KERNELS1%UEA(IR,J,IP0,K0,I0)
        YPO(2)=KERNELS1%UEA(IR,J,IP1,K0,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(IR,J,IP0,K0,I1)
        YPO(2)=KERNELS1%UEA(IR,J,IP1,K0,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%UEA(IR,J,1:KN1,K0,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(IR,J,1:KN1,K0,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        AA(1)=LINEAR( X, Y, 2, RL )

      IF(keyLS .EQ. 1) THEN
        YPO(1)=KERNELS1%UEA(IR,J,IP0,K1,I0)
        YPO(2)=KERNELS1%UEA(IR,J,IP1,K1,I0)
        Y(1)=LINEAR( XPO, YPO, 2, PO )
        YPO(1)=KERNELS1%UEA(IR,J,IP0,K1,I1)
        YPO(2)=KERNELS1%UEA(IR,J,IP1,K1,I1)
        Y(2)=LINEAR( XPO, YPO, 2, PO )
	  ELSE
	    YYS1(1:KN1)=KERNELS1%UEA(IR,J,1:KN1,K1,I0)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(1)=YFIT
	    YYS1(1:KN1)=KERNELS1%UEA(IR,J,1:KN1,K1,I1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
        Y(2)=YFIT
	  ENDIF
        BB(1)=LINEAR( X, Y, 2, RL )      

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      Y(1)=AA(1)
      Y(2)=BB(1)
      AB(1)=LINEAR( X, Y, 2, log(RI) )
      sumUSEA=sumUSEA+AB(1)*RDc(IR)
      ENDDO ! IR KR

      USEA(J,I,1)=sumUSEA*WAVEL/WL(1)/sumRD     
      
      ENDDO ! J 2

      ENDDO ! I KN           
      
	RETURN 
	END SUBROUTINE USU_LS_RD_bin

!***************************



