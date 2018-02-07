! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!c **
!c calculate US for opt.char. from original WL kernels  
!c **************************************************************** 
      SUBROUTINE USU_WL(key_RD,keyEL,keySUB,KR,R,RD,  &
                        RATIO,NRATN,                  &
                        KM,KRE,KIM,ARE,AIM,           &
                        RN,RK,KC,LB,LE,KERNELS2) 
      !use alloc2
      USE mod_alloc_kernels
      use mod_intrpl_linear
      USE mod_par_DLS
      USE mod_par_DLS_bin
	  
      type(kernels_lognormal_bin), intent(in)  ::  KERNELS2

      real*4 R(KRpar), RATIO(KR1par), RRATN
      real RD(KRpar)
	   real US11,US12,US22,US33,US34,US44,USEA
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
      real RDc(KR),AA(6),BB(6),AB(6)   &
                       ,CC(6),DD(6),CD(6)   &
                       ,sumUS(6)
      real X(2), Y(2), ARE(KREpar), AIM(KIMpar)
!      real LINEAR
!      real,external	:: LINEAR
      integer KC,LB,LE,keyEL,keySUB
      real RN(KWLpar), RK(KWLpar)

      PI=ACOS(-1.)
      PI2=2.*PI

       IF(key_RD.eq.2) then
 
!c*** RECALCULATE ASPECT RATIO DISTRIBUTION (RDc()=SAREA/VOLUME)
!c*** RDc()=RD()/RDc(); sumRD=sum(RDc())

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
!c          write(*,*) 'R=',R(IR),' B=',RDc(IR),  &
!c        ' 1/B=',1./RDc(IR),' RD=',RD(IR)

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
	  
      DO L=LB,LE ! WL loop start

      RL=RN(L)
      RI=RK(L)
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
        IF(RL.LT.ARE(1).and.keySUB.eq.0) THEN
        WRITE(*,*) 'n=',RN(L),' is out of the range:',  &
                        ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(1)
        ENDIF
        I0=1
        I1=2
        RN(L)=ARE(1)
        RL=ARE(1)
      ENDIF
      IF(RL.GE.ARE(KRE)) THEN
        IF(RL.GT.ARE(KRE).and.keySUB.eq.0) THEN
        WRITE(*,*) 'n=',RN(L),' is out of the range:',  &
                        ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(KRE)
        ENDIF
        I0=KRE-1
        I1=KRE
        RN(L)=ARE(KRE)
        RL=ARE(KRE)
      ENDIF
      IF(RI.LE.AIM(1)) THEN
        IF(RI.LT.AIM(1).and.keySUB.eq.0) THEN
        WRITE(*,*) 'k=',RK(L),' is out of the range:',  &
                        AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(1)
        ENDIF
        K0=1
        K1=2
        RK(L)=AIM(1)
        RI=AIM(1)
      ENDIF
      IF(RI.GE.AIM(KIM)) THEN
        IF(RI.GT.AIM(KIM).and.keySUB.eq.0) THEN
        WRITE(*,*) 'k=',RK(L),' is out of the range:',  &
                        AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(KIM)
        ENDIF
        K0=KIM-1
        K1=KIM
        RK(L)=AIM(KIM)
        RI=AIM(KIM)
      ENDIF
!C      write(*,*) I0,I1,K0,K1,' I0,I1,K0,K1'
!C      WRITE(*,*) ARE(I0),ARE(I1),' ARE'
!C      WRITE(*,*) AIM(K0),AIM(K1),' AIM'
      DO I=1,KC

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
        WRITE(*,*) 'R has been changed',  &
                        R(IR),' => ',RATIO(1)
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
        WRITE(*,*) 'R=',R(IR),' is out of the range:',  &
                        RATIO(1),'< R <',RATIO(NRATN)
        WRITE(*,*) 'R has been changed',R(IR),' => ',RATIO(1)
        ENDIF
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      IF(RRATN.GE.RATIO(NRATN)) THEN
        IF(RRATN.GT.RATIO(NRATN)) THEN
        WRITE(*,*) 'R=',R(IR),' is out of the range:',  &
                        RATIO(1),'< R <',RATIO(NRATN)
        WRITE(*,*) 'R has been changed',R(IR),          &
                                  ' => ',RATIO(NRATN)
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
!cl      write(*,*) L0,L1,R(IR),RRATN,
!cl     &  ' L0,L1 R(IR),RRATN for phase function'
!cl      write(*,*) I0,I1,K0,K1,IP0,IP1,
!cl     &  ' I0,I1,K0,K1,IP0,IP1,'

      X(1)=ARE(I0)
      X(2)=ARE(I1)
!c
!c ** U11     AA(1)
!c
      Y(1)=KERNELS2%UO11(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO11(L0,J,I,K0,I1,L)
      AA(1)=LINEAR( X, Y, 2, RL )
!cl      if(J.eq.1.and.I.eq.1) then
!cl      write(*,*) 'AA(1)=',AA(1)
!cl	endif
!c
!c ** U12     AA(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO12(L0,J,I,K0,I1,L)
      AA(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     AA(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO22(L0,J,I,K0,I1,L)
      AA(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     AA(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO33(L0,J,I,K0,I1,L)
      AA(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     AA(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO34(L0,J,I,K0,I1,L)
      AA(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     AA(6)
!c
      if(keyEL.gt.5) then
	   Y(1)=KERNELS2%UO44(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO44(L0,J,I,K0,I1,L)
      AA(6)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U11     BB(1)
!c
      Y(1)=KERNELS2%UO11(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO11(L0,J,I,K1,I1,L)
      BB(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     BB(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO12(L0,J,I,K1,I1,L)
      BB(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     BB(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO22(L0,J,I,K1,I1,L)
      BB(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     BB(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO33(L0,J,I,K1,I1,L)
      BB(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     BB(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO34(L0,J,I,K1,I1,L)
      BB(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     BB(6)
!c
      if(keyEL.gt.5) then
      Y(1)=KERNELS2%UO44(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO44(L0,J,I,K1,I1,L)
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
      Y(1)=KERNELS2%UO11(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO11(L1,J,I,K0,I1,L)
      CC(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     CC(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO12(L1,J,I,K0,I1,L)
      CC(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     CC(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO22(L1,J,I,K0,I1,L)
      CC(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     CC(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO33(L1,J,I,K0,I1,L)
      CC(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     CC(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO34(L1,J,I,K0,I1,L)
      CC(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     CC(6)
!c
      if(keyEL.gt.5) then
      Y(1)=KERNELS2%UO44(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO44(L1,J,I,K0,I1,L)
      CC(6)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U11     DD(1)
!c
      Y(1)=KERNELS2%UO11(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO11(L1,J,I,K1,I1,L)
      DD(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     DD(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO12(L1,J,I,K1,I1,L)
      DD(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     DD(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO22(L1,J,I,K1,I1,L)
      DD(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     DD(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO33(L1,J,I,K1,I1,L)
      DD(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     DD(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO34(L1,J,I,K1,I1,L)
      DD(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     DD(6)
!c
      if(keyEL.gt.5) then
      Y(1)=KERNELS2%UO44(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO44(L1,J,I,K1,I1,L)
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

                     US11(J,I,L)=sumUS(1)/sumRD
      if(keyEL.gt.1) US12(J,I,L)=sumUS(2)/sumRD
      if(keyEL.gt.2) US22(J,I,L)=sumUS(3)/sumRD
      if(keyEL.gt.3) US33(J,I,L)=sumUS(4)/sumRD
      if(keyEL.gt.4) US34(J,I,L)=sumUS(5)/sumRD
      if(keyEL.gt.5) US44(J,I,L)=sumUS(6)/sumRD

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
      IF(NRATN.NE.1) THEN
      DO IRATN=1,NRATN-1
      IF(RRATN.GE.RATIO(IRATN).AND.RRATN.LE.RATIO(IRATN+1)) THEN
          L0=IRATN
          L1=IRATN+1
      ENDIF
      ENDDO

      IF(RRATN.GE.RATIO(NRATN)) THEN
        L0=NRATN-1
        L1=NRATN
        R(IR)=RATIO(NRATN)
        RRATN=RATIO(NRATN)
      ENDIF
      IF(RRATN.LE.RATIO(1)) THEN
        L0=1
        L1=2
        R(IR)=RATIO(1)
        RRATN=RATIO(1)
      ENDIF
      ELSE
        L0=1
        L1=1
      ENDIF

      RRATN1=RRATN

      X(1)=ARE(I0)
      X(2)=ARE(I1)
      
      Y(1)=KERNELS2%UOEA(L0,J,I,K0,I0,L)
      Y(2)=KERNELS2%UOEA(L0,J,I,K0,I1,L)
      AA(1)=LINEAR( X, Y, 2, RL )

      Y(1)=KERNELS2%UOEA(L0,J,I,K1,I0,L)
      Y(2)=KERNELS2%UOEA(L0,J,I,K1,I1,L)
      BB(1)=LINEAR( X, Y, 2, RL )

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      Y(1)=AA(1)
      Y(2)=BB(1)
      AB(1)=LINEAR( X, Y, 2, log(RI) )

      IF(NRATN.NE.1) THEN    
      X(1)=ARE(I0)
      X(2)=ARE(I1)

      Y(1)=KERNELS2%UOEA(L1,J,I,K0,I0,L)
      Y(2)=KERNELS2%UOEA(L1,J,I,K0,I1,L)
      CC(1)=LINEAR( X, Y, 2, RL )

      Y(1)=KERNELS2%UOEA(L1,J,I,K1,I0,L)
      Y(2)=KERNELS2%UOEA(L1,J,I,K1,I1,L)
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
      ELSEIF(NRATN.EQ.1) THEN
      sumUSEA=AB(1)*RDc(IR)      
      ENDIF
      ENDDO ! IR KR

      USEA(J,I,L)=sumUSEA/sumRD     
      
      ENDDO ! J 2

      ENDDO ! I KC
	
	ENDDO ! L NWL            
      
	RETURN 
	END SUBROUTINE USU_WL
!c
!c***************************************
!c
!c ** calculate and save UO.. WL original kernels using U...
!c    original kernels 
!c
      SUBROUTINE USUO_WL(  keyEL,                                     &
                           KN1,grid1,KM,ANGLE,WAVEL,KRE,KIM,ARE,AIM,  &
                           KN,grid,NWL,WL,KR,KC,SD,                   &
                           xnmin,xnmax,xkmin,xkmax,                   &
                           KERNELS1,KERNELS2  )   
      !use alloc1 
      !use alloc2 
      use mod_alloc_kernels
      use mod_intrpl_linear
	    use mod_intrpl_spline
      USE mod_par_DLS
      USE mod_par_DLS_bin
	  
      type(kernels_triangle_bin),  intent(in)     ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2
      
      real grid(KNpar),grid1(KN1par),ANGLE(KMpar) 
!c      real*4 R(KRpar), RD(KRpar)
	real   SD(KCpar,KNpar)
	integer KR,KC,NWL,KN1,KN,KM
	real*4 xnmin,xnmax,xkmin,xkmax  &
            ,RABRE(KREpar), RABIM(KIMpar)
	integer nn31,nn32

	real*4       UR11(KMpar,KNpar,KIMpar,KREpar)  &
              , UR12(KMpar,KNpar,KIMpar,KREpar)  &
              , UR22(KMpar,KNpar,KIMpar,KREpar)  &
              , UR33(KMpar,KNpar,KIMpar,KREpar)  &
              , UR34(KMpar,KNpar,KIMpar,KREpar)  &
              , UR44(KMpar,KNpar,KIMpar,KREpar)  &
              , UREA(2    ,KNpar,KIMpar,KREpar)

   real Y(2),ARE(KREpar), AIM(KIMpar)
	real WL(KWLpar)
	integer keyEL
!c **  for SPLINE subroutine
      DOUBLE PRECISION :: XARG, YFIT
      DOUBLE PRECISION, DIMENSION(KN1par) :: XXS1, YYS1, b, c, d

!cl      real*4 T_UO, T_U=0.
!cl      real*4 tarray(2)
!cl      real*4 external dtime, etime
!cl                                 if(keySUB.eq.0) then
!cl								    T_UO=dtime(tarray) !+++
!cl									T_UO=tarray(1)+tarray(2)
!cl                                 endif
      PI=ACOS(-1.)
      PI2=2.*PI
!c
!c Define N and K index (redefine KRE & KIM)
!c
!c  N
      ind1=0
	do i=1,KRE-1
!c	write(*,*) 'i=',i,' xnmin=',xnmin,' ARE(i+1)=',ARE(i+1)
      if(xnmin.le.ARE(i+1)) then
	nn11=i
	ind1=ind1+1
	endif
      if(ind1.eq.1) EXIT
	enddo ! i
!c	 
      ind1=0
	do i=KRE,2,-1
!c	write(*,*) 'i=',i,' xnmax=',xnmax,' ARE(i-1)=',ARE(i-1)
      if(xnmax.ge.ARE(i-1)) then
	nn21=i
	ind1=ind1+1
	endif
      if(ind1.eq.1) EXIT
	enddo ! i
!c
      if(xnmin.eq.ARE(1)) then
	nn11=1
	endif ! xnmin 
      if(xnmax.eq.ARE(KRE)) then
	nn21=KRE
	endif ! xnmax
      nn31=nn21-nn11+1
!c	write(*,*) 'nn11=',nn11,' nn21=',nn21,' nn31=',nn31
!c  K
      ind2=0
	do i=1,KIM-1
!c	write(*,*) 'i=',i,' xkmin=',xkmin,' AIM(i+1)=',AIM(i+1)
      if(xkmin.le.AIM(i+1)) then
	nn12=i
	ind2=ind2+1
	endif
      if(ind2.eq.1) EXIT
	enddo ! i
	 
      ind2=0
	do i=KIM,2,-1
!c	write(*,*) 'i=',i,' xkmax=',xkmax,' AIM(i-1)=',AIM(i-1)
      if(xkmax.ge.AIM(i-1)) then
	nn22=i
	ind2=ind2+1
	endif
      if(ind2.eq.1) EXIT
	enddo ! i

      if(xkmin.eq.AIM(1)) then
	nn12=1
	endif ! xkmin 
      if(xkmax.eq.AIM(KIM)) then
	nn22=KIM
	endif ! xkmax
      nn32=nn22-nn12+1
!c	write(*,*) 'nn12=',nn12,' nn22=',nn22,' nn32=',nn32

      DO IR=1,KR
      DO L=1,NWL
!cl       cinv=WAVEL/WL(L)
      DO I=1,KN

      IP0=0
      IP1=0
      PO=grid(I)*PI2/WL(L)            
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
!c
!c **  SCATTERING MATRIX ELEMENTS
!c

	  XXS1(1:KN1)=grid1(1:KN1)*PI2/WAVEL
	  XARG=PO

       ind1=nn11-1      
	DO N=1,nn31
	 ind1=ind1+1
	 ind2=nn12-1
!c      write(*,*) 'TEST_N: N=',N,' ind1=',ind1
	DO K=1,nn32
       ind2=ind2+1
!c      write(*,*) 'TEST_K: N=',N,' ind2=',ind2
      DO J=1,KM
!c
!c ** U11     
!c
	    YYS1(1:KN1)=KERNELS1%U11(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=EXP(YFIT)
      UR11(J,I,K,N)=Y(1)/WL(L)
!c
!c ** U12   
!c
      IF(keyEL.GT.1) THEN
	    YYS1(1:KN1)=KERNELS1%U12(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
      UR12(J,I,K,N)=Y(1)/WL(L)
      ENDIF
!c
!c ** U22    
!c
      IF(keyEL.GT.2) THEN      
	    YYS1(1:KN1)=KERNELS1%U22(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=EXP(YFIT)
      UR22(J,I,K,N)=Y(1)/WL(L)
      ENDIF
!c
!c ** U33     
!c
      IF(keyEL.GT.3) THEN 
	    YYS1(1:KN1)=KERNELS1%U33(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
      IF(ANGLE(J).LE.40.) THEN
      Y(1)=EXP(Y(1))  	  
      ENDIF ! ANGLE
      UR33(J,I,K,N)=Y(1)/WL(L)
      ENDIF
!c
!c ** U34     
!c
      IF(keyEL.GT.4) THEN
	    YYS1(1:KN1)=KERNELS1%U34(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
      UR34(J,I,K,N)=Y(1)/WL(L)
      ENDIF
!c
!c ** U44    
!c
      IF(keyEL.GT.5) THEN
	    YYS1(1:KN1)=KERNELS1%U44(IR,J,1:KN1,ind2,ind1)
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT
      IF(ANGLE(J).LE.50.) THEN
      Y(1)=EXP(Y(1))
      ENDIF ! ANGLE
      UR44(J,I,K,N)=Y(1)/WL(L)
      ENDIF
      
      ENDDO   ! J KM
!c
!c ** EXTINCTION & ABSORPTION
!c
!cl	XXS1(1:KN1)=grid1(1:KN1)*PI2/WAVEL
!cl	XARG=PO

      DO J=1,2   
	    YYS1(1:KN1)=KERNELS1%UEA(IR,J,1:KN1,ind2,ind1)*WAVEL
      CALL intrpl_spline(KN1, XXS1(1:KN1), YYS1(1:KN1), XARG, YFIT, 1, 1, b, c, d)
      Y(1)=YFIT	      
      UREA(J,I,K,N)=Y(1)/WL(L)

      ENDDO ! J 2

      ENDDO ! K
	ENDDO ! N

      ENDDO ! I KN
!c      write(*,*) 'TEST after KN DO loop'

!c Calculate KC original Kernels 

	DO N=1,nn31
	DO K=1,nn32
	DO M=1,KC

      KERNELS2%UOEA(IR,1,M,K,N,L)=DOT_PRODUCT(UREA(1,1:KN,K,N),SD(M,1:KN))
      KERNELS2%UOEA(IR,2,M,K,N,L)=DOT_PRODUCT(UREA(2,1:KN,K,N),SD(M,1:KN))
      do j=1,KM
      KERNELS2%UO11(IR,j,M,K,N,L)=DOT_PRODUCT(UR11(j,1:KN,K,N),SD(M,1:KN))
      enddo ! j
        
      if(keyEL.gt.1) then
      do j=1,KM
      KERNELS2%UO12(IR,j,M,K,N,L)=DOT_PRODUCT(UR12(j,1:KN,K,N),SD(M,1:KN))
      enddo ! j
      endif
      if(keyEL.gt.2) then
      do j=1,KM
      KERNELS2%UO22(IR,j,M,K,N,L)=DOT_PRODUCT(UR22(j,1:KN,K,N),SD(M,1:KN))
      enddo ! j
	  endif
      if(keyEL.gt.3) then
      do j=1,KM
      KERNELS2%UO33(IR,j,M,K,N,L)=DOT_PRODUCT(UR33(j,1:KN,K,N),SD(M,1:KN))
      enddo ! j
	  endif
	  if(keyEL.gt.4) then
      do j=1,KM
      KERNELS2%UO34(IR,j,M,K,N,L)=DOT_PRODUCT(UR34(j,1:KN,K,N),SD(M,1:KN))
	  enddo ! j
      endif 
	  if(keyEL.gt.5) then
      do j=1,KM
      KERNELS2%UO44(IR,j,M,K,N,L)=DOT_PRODUCT(UR44(j,1:KN,K,N),SD(M,1:KN))
	  enddo ! j
      endif 

	ENDDO ! M KC
	ENDDO ! K
	ENDDO ! N

	ENDDO ! L NWL
	ENDDO ! IR KR

!c redefine ARE,AIM arrays
       ind1=nn11-1      
	DO N=1,nn31
	 ind1=ind1+1
       RABRE(N)=ARE(ind1)
	ENDDO ! N
	ARE(:)=0.
	ARE(1:nn31)=RABRE(1:nn31)
      KRE=nn31
!c	write(*,*) '1; KRE=',KRE,' nn31=',nn31
!c	write(*,*) ARE(1),' <ARE< ',ARE(KRE) 

	 ind2=nn12-1
	DO K=1,nn32
	 ind2=ind2+1
	 RABIM(K)=AIM(ind2)
	ENDDO ! K
	AIM(1:nn32)=0.
	AIM(1:nn32)=RABIM(1:nn32)
	KIM=nn32
!c	write(*,*) '1: KIM=',KIM,' nn32=',nn32
!c	write(*,*) AIM(1),' <AIM< ',AIM(KIM) 
!cl							if(keySUB.eq.0) then
!cl								   T_UO=dtime(tarray)
!cl	                               T_UO=tarray(1)+tarray(2)
!cl								   T_U=T_U+T_UO 
!cl
!cl      WRITE(*,61) T_U/60.     
!cl                            endif
!   61 format('  Interpolation UO ........ ',f8.3,' min.')

      RETURN
      END SUBROUTINE USUO_WL  
	
!c **************************************************************** 

      SUBROUTINE USU_WL_RD(key_RD,keyEL,keySUB,KR,R,RD,  &
                           RATIO,                        &
                           KM,KRE,KIM,ARE,AIM,           &
                           RN,RK,KC,LB,LE,               &
                           KERNELS2) 
      !use alloc2
      use mod_alloc_kernels      
      use mod_intrpl_linear
      USE mod_par_DLS
      USE mod_par_DLS_bin
!  -------------------------------------------------------------------------	  
      type(kernels_lognormal_bin), intent(in)  ::  KERNELS2

      real*4 :: R(KRpar), RATIO(KR1par)
      real   :: RD(KRpar)
      real   :: US11,US12,US22,US33,US34,US44,USEA
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
      real :: RDc(KR),AA(6),BB(6),AB(6),sumUS(6)
      real :: X(2),Y(2),ARE(KREpar),AIM(KIMpar)
!      real LINEAR
!      real,external	:: LINEAR
	  integer :: KC,LB,LE,keyEL,keySUB
	  real    :: RN(KWLpar), RK(KWLpar)
!  -------------------------------------------------------------------------

      PI=ACOS(-1.)
      PI2=2.*PI

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
	  
      DO L=LB,LE ! WL loop start

      RL=RN(L)
      RI=RK(L)
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
        IF(RL.LT.ARE(1).and.keySUB.eq.0) THEN
        WRITE(*,*) 'n=',RN(L),' is out of the range:',  &
                        ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(1)
        ENDIF
        I0=1
        I1=2
        RN(L)=ARE(1)
        RL=ARE(1)
      ENDIF
      IF(RL.GE.ARE(KRE)) THEN
        IF(RL.GT.ARE(KRE).and.keySUB.eq.0) THEN
        WRITE(*,*) 'n=',RN(L),' is out of the range:',  &
                        ARE(1),'< n <',ARE(KRE)
        WRITE(*,*) 'n has been changed',RN,' => ',ARE(KRE)
        ENDIF
        I0=KRE-1
        I1=KRE
        RN(L)=ARE(KRE)
        RL=ARE(KRE)
      ENDIF
      IF(RI.LE.AIM(1)) THEN
        IF(RI.LT.AIM(1).and.keySUB.eq.0) THEN
        WRITE(*,*) 'k=',RK(L),' is out of the range:',  &
                        AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(1)
        ENDIF
        K0=1
        K1=2
        RK(L)=AIM(1)
        RI=AIM(1)
      ENDIF
      IF(RI.GE.AIM(KIM)) THEN
        IF(RI.GT.AIM(KIM).and.keySUB.eq.0) THEN
        WRITE(*,*) 'k=',RK(L),' is out of the range:',  &
                       AIM(1),'< k <',AIM(KIM)
        WRITE(*,*) 'k has been changed',RK,' => ',AIM(KIM)
        ENDIF
        K0=KIM-1
        K1=KIM
        RK(L)=AIM(KIM)
        RI=AIM(KIM)
      ENDIF
!C      write(*,*) I0,I1,K0,K1,' I0,I1,K0,K1'
!C      WRITE(*,*) ARE(I0),ARE(I1),' ARE'
!C      WRITE(*,*) AIM(K0),AIM(K1),' AIM'
      DO I=1,KC
!c
!c **  SCATTERING MATRIX ELEMENTS
!c
      DO J=1,KM
      sumUS(:6)=0.

      DO IR=1,KR ! axis ratio loop
	  
        IF(RATIO(IR).NE.R(IR)) THEN
        WRITE(*,*) 'R=',R(IR),' .NE. RATIO=',RATIO(IR)
        WRITE(*,*) 'STOP in DLS_intrpl_bin.f90'
		STOP
        ENDIF

      X(1)=ARE(I0)
      X(2)=ARE(I1)
!c
!c ** U11     AA(1)
!c
      Y(1)=KERNELS2%UO11(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO11(IR,J,I,K0,I1,L)
      AA(1)=LINEAR( X, Y, 2, RL )
!cl      if(J.eq.1.and.I.eq.1) then
!cl      write(*,*) 'AA(1)=',AA(1)
!cl	endif
!c
!c ** U12     AA(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO12(IR,J,I,K0,I1,L)
      AA(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     AA(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO22(IR,J,I,K0,I1,L)
      AA(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     AA(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO33(IR,J,I,K0,I1,L)
      AA(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     AA(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO34(IR,J,I,K0,I1,L)
      AA(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     AA(6)
!c
      if(keyEL.gt.5) then
	   Y(1)=KERNELS2%UO44(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UO44(IR,J,I,K0,I1,L)
      AA(6)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U11     BB(1)
!c
      Y(1)=KERNELS2%UO11(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO11(IR,J,I,K1,I1,L)
      BB(1)=LINEAR( X, Y, 2, RL )
!c
!c ** U12     BB(2)
!c
      if(keyEL.gt.1) then
      Y(1)=KERNELS2%UO12(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO12(IR,J,I,K1,I1,L)
      BB(2)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U22     BB(3)
!c
      if(keyEL.gt.2) then
      Y(1)=KERNELS2%UO22(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO22(IR,J,I,K1,I1,L)
      BB(3)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U33     BB(4)
!c
      if(keyEL.gt.3) then
      Y(1)=KERNELS2%UO33(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO33(IR,J,I,K1,I1,L)
      BB(4)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U34     BB(5)
!c
      if(keyEL.gt.4) then
      Y(1)=KERNELS2%UO34(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO34(IR,J,I,K1,I1,L)
      BB(5)=LINEAR( X, Y, 2, RL )
      endif
!c
!c ** U44     BB(6)
!c
      if(keyEL.gt.5) then
      Y(1)=KERNELS2%UO44(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UO44(IR,J,I,K1,I1,L)
      BB(6)=LINEAR( X, Y, 2, RL )
      endif

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))

      DO II=1,keyEL
      Y(1)=AA(II)
      Y(2)=BB(II)
      AB(II)=LINEAR( X, Y, 2, log(RI) )
      sumUS(II)=sumUS(II)+AB(II)*RDc(IR)
      ENDDO ! II

      ENDDO ! IR RD

      US11(J,I,L)=sumUS(1)/sumRD
      if(keyEL.gt.1) US12(J,I,L)=sumUS(2)/sumRD
      if(keyEL.gt.2) US22(J,I,L)=sumUS(3)/sumRD
      if(keyEL.gt.3) US33(J,I,L)=sumUS(4)/sumRD
      if(keyEL.gt.4) US34(J,I,L)=sumUS(5)/sumRD
      if(keyEL.gt.5) US44(J,I,L)=sumUS(6)/sumRD

      ENDDO   ! J KM
!c
!c ** EXTINCTION & ABSORPTION
!c
      DO J=1,2
      sumUSEA=0.
      DO IR=1,KR

      X(1)=ARE(I0)
      X(2)=ARE(I1)
      
      Y(1)=KERNELS2%UOEA(IR,J,I,K0,I0,L)
      Y(2)=KERNELS2%UOEA(IR,J,I,K0,I1,L)
      AA(1)=LINEAR( X, Y, 2, RL )

      Y(1)=KERNELS2%UOEA(IR,J,I,K1,I0,L)
      Y(2)=KERNELS2%UOEA(IR,J,I,K1,I1,L)
      BB(1)=LINEAR( X, Y, 2, RL )

      X(1)=log(AIM(K0))
      X(2)=log(AIM(K1))
      Y(1)=AA(1)
      Y(2)=BB(1)
      AB(1)=LINEAR( X, Y, 2, log(RI) )

      sumUSEA=sumUSEA+AB(1)*RDc(IR)      

      ENDDO ! IR KR

      USEA(J,I,L)=sumUSEA/sumRD     
      
      ENDDO ! J 2

      ENDDO ! I KC
	
	ENDDO ! L NWL            
      
	RETURN 
	END SUBROUTINE USU_WL_RD

!c***************************************
