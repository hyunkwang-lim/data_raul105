! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module phase_func
      implicit none 
      contains

!c **************************************************************** 	 

      subroutine SINT(X1, Y1, KM, xnorm)
      use mod_intrpl_linear
!c ** Simpson integration
!c ** xh - angle step (degree)
      integer KM, KSIMP, KSIMP1, IK, II	 
      real :: X1( KM ), Y1( KM ), Y(2)
      real xnorm, xh, pi
      real XA, F11
      KSIMP=721
      xh=180./float(KSIMP-1)
      KSIMP1=INT((180.-X1(1))/xh+1.)
!c	write(*,*) 'KSIMP1=',KSIMP1
      pi=ACOS(-1.)
      Y(1)=Y1(1)
      Y(2)=Y1(KM)
      Y1(1:KM)=LOG(Y1(1:KM))
      xnorm=0.
      XA=X1(1)
      do IK=1,KSIMP1
         II=IK/2
         if(IK .eq. 1) then
            F11=1./3.*Y(1)
         else if(IK .eq. KSIMP) then
            F11=1./3.*Y(2)
         else if(II*2 .lt. IK) then
            F11=2./3.*EXP(LINEAR(X1, Y1, KM, XA))
         else if(II*2 .eq. IK) then
            F11=4./3.*EXP(LINEAR(X1, Y1, KM, XA))
         endif
         xnorm=xnorm+F11*SIN(XA*pi/180.)
         XA=XA+xh
      enddo ! IK
      xnorm=0.5*xnorm*xh*pi/180.
!cl      write(*,'(''xnorm='',e13.5)') xnorm
      Y1(1:KM)=EXP(Y1(1:KM))/xnorm

      return
      end subroutine SINT

!c ****************************************************************
! under development
      subroutine SINT_SPLINE(X1, Y1, KM, xnorm)
      use mod_intrpl_spline
!c ** Simpson integration
!c ** xh - angle step (degree)
      integer KM, KSIMP, KSIMP1, IK, II 	 
      real :: X1( KM ), Y1( KM ), Y(2)
      real xnorm, xh, pi, F11
      INTEGER :: icurr
      DOUBLE PRECISION :: XA, YFIT
      DOUBLE PRECISION, DIMENSION(KM) :: XXS2, YYS2, b, c, d

      KSIMP=721
      xh=180./float(KSIMP-1)
      KSIMP1=INT((180.-X1(1))/xh+1.)
!c	write(*,*) 'KSIMP1=',KSIMP1
      pi=ACOS(-1.)
      Y(1)=Y1(1)
      Y(2)=Y1(KM)
      Y1(1:KM)=LOG(Y1(1:KM))

      XXS2(1:KM)=X1(1:KM)
      YYS2(1:KM)=Y1(1:KM)
      icurr = 1
	      
      xnorm=0.
      XA=X1(1)

      do IK=1,KSIMP1
        !IF(IK .EQ. 2) THEN
        !  key_spln = 1
        !ELSE
        !  key_spln = 2
        !ENDIF ! IK
        II=IK/2

        if(IK .eq. 1) then
          F11=1./3.*Y(1)
        GOTO 1
        endif
        if(IK .eq. KSIMP) then
          F11=1./3.*Y(2)
          GOTO 1
        endif

        IF(II*2 .LT. IK) THEN
          CALL intrpl_spline(KM, XXS2(1:KM), YYS2(1:KM), XA, YFIT, 1, icurr, b, c, d)
          F11=2./3.*EXP(YFIT)
        GOTO 1
        ENDIF
        IF(II*2 .EQ. IK) THEN
          CALL intrpl_spline(KM, XXS2(1:KM), YYS2(1:KM), XA, YFIT, 1, icurr, b, c, d)
          F11=4./3.*EXP(YFIT)
        GOTO 1
        ENDIF
1       CONTINUE
        xnorm=xnorm+F11*SIN(XA*pi/180.)
        XA=XA+xh
        icurr = icurr + 1
      enddo ! IK

      xnorm=0.5*xnorm*xh*pi/180.
!cl      write(*,'(''xnorm='',e13.5)') xnorm
      Y1(1:KM)=EXP(Y1(1:KM))/xnorm

      return
      end subroutine SINT_SPLINE
!c	
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine ASYPAR(X2, Y2, KM, g)

      use mod_intrpl_linear
! --------------------------------------------------
      integer           , intent(in)  :: KM
      real,dimension(KM), intent(in)  :: X2, Y2
      real              , intent(out) :: g
! --------------------------------------------------
!c ** Asymmetry parameter
!c ** xh - angle step (degree)
      real :: X1( KM ), Y1( KM ), Y(2)
!      real,external	:: LINEAR
      real  :: xh, pi	  
      real  :: XA, xpi, F11, XX
      integer :: ii,ik,ksimp
! --------------------------------------------------

      X1(1:KM)=X2(1:KM)
      Y1(1:KM)=LOG(Y2(1:KM))
      KSIMP=721
      xh=180./float(KSIMP-1)
      pi=ACOS(-1.)
      xpi=pi/180.
      Y(1)=Y1(1)
      Y(2)=Y1(KM)
      g=0.
      XA=0.
	  
      do IK=1,KSIMP
         II=IK/2
         if(IK .eq. 1) then
            F11=1./3.*Y(1)
         else if(IK .eq. KSIMP) then
            F11=1./3.*Y(2)
         else if(II*2 .lt. IK) then
            F11=2./3.*EXP(LINEAR(X1, Y1, KM, XA))
         else if(II*2 .eq. IK) then
            F11=4./3.*EXP(LINEAR(X1, Y1, KM, XA))
         endif
         XX=XA*xpi
         g=g+F11*SIN(XX)*COS(XX)
         XA=XA+xh
      enddo ! IK
      g=0.5*g*(xh*xpi)

      return
      end subroutine ASYPAR

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module phase_func


