! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! file contains :

! SUBROUTINE RDMG
! SUBROUTINE RDMU1
! SUBROUTINE RDMU2 
!
!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE RDMG(R0,N,EM,SGM,R)
!C RANDOM NUMBER GENERATOR
!C R = exp(-(x-m)**2/2/s**2)/sqrt(2*pi)/s
!C--- HISTORY
!C 88.10.14 CREATED
!C--- INPUT
!C R0    D    INITIAL CONDITION (0 TO 1)
!C N     I    NUMBER OF RANDOM NUMBERS
!C EM    R    Mean (m)
!C SGM   R    Sigma (s)
!C--- OUTPUT
!C R    R(N)  Normal distribution RANDOM NUMBERS (-inf, +inf)
!C$ENDI
!CDD      implicit real * 8 (a-h, o-z)
      !USE mo_par_inv, only : KMESS
      IMPLICIT NONE
	  
      SAVE   ! ??? ask Oleg
      INTEGER,PARAMETER :: NU=12
      REAL*8 R0
      !REAL,DIMENSION(KMESS) :: R
      REAL,DIMENSION(N)  :: R

      REAL,DIMENSION(NU) :: W
      INTEGER :: N,L,I
      REAL    :: S,SGM,EM

      DO  L=1,N
      CALL RDMU1(R0,NU,W)
      S=0.
      DO  I=1,NU
      S=S+W(I)
      ENDDO ! I
      R(L)=SGM*(S-NU*0.5)+EM
      ENDDO ! L
  
      RETURN
      END SUBROUTINE RDMG

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE RDMU1(R0,N,R)
!C RANDOM NUMBER GENERATOR
!C--- HISTORY
!C 88.10.14 CREATED
!C--- INPUT
!C R0    D    INITIAL CONDITION (0 TO 1)
!C N     I    NUMBER OF RANDOM NUMBERS
!C--- OUTPUT
!C R    R(N)  UNIFORM RANDOM NUMBERS ( 0 TO 1)
!C$ENDI
      IMPLICIT NONE
	  
      SAVE  ! ??? ask Oleg
!CDD      implicit real * 8 (a-h, o-z)
      INTEGER :: N,I
      REAL,DIMENSION(N) :: R
      REAL*8 R0,R1,PI
      PARAMETER (PI=3.141592653589793)
      DO  I=1,N
      R0=(PI+R0)**5
      R1=INT(R0)
      R0=R0-R1
      R(I)=R0
      ENDDO ! I
	  
      RETURN
      END SUBROUTINE RDMU1

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE RDMU2(N,R)
!		use IFPORT  ! ifort version
!C RANDOM NUMBER GENERATOR
!C--- HISTORY
!C 88.10.14 CREATED
!C--- INPUT
!C R0    D    INITIAL CONDITION (0 TO 1)
!C N     I    NUMBER OF RANDOM NUMBERS
!C--- OUTPUT
!C R    R(N)  UNIFORM RANDOM NUMBERS ( 0 TO 1)
!C$ENDI
!      SAVE ! disabled
!CDD      implicit real * 8 (a-h, o-z)
!	------------------------------------------------------------------------------
      integer,intent(in)	::	N
!	------------------------------------------------------------------------------

!		   integer		::	cur_time
!		   integer(4)		::	cur_time ! for ifort ???
      integer(8)		::	cur_time
      integer		::	I
      real(8)	::	R0_local,R1
!		   real(8)	::	R0
      real(8),dimension(N) :: R
! REAL*8 R0,R1,PI
      REAL*8 PI
      PARAMETER (PI=3.141592653589793)
!	------------------------------------------------------------------------------
      CALL SYSTEM_CLOCK(COUNT=cur_time)
!		   cur_time = time8()	!  gfortran version
      !write(*,*) 'cur_time=',cur_time,'  SYSTEM_CLOCK'
!      cur_time = time() 	! ifort version
      !write(*,*) 'cur_time=',cur_time,'  time()'

!	------------------------------------------------------------------------------
!		  R0_local = 0.0
		  R0_local = mod(dble(cur_time),PI)/PI

      DO I=1,N
		  R0_local=(PI+R0_local)**5
		  R1=INT(R0_local)
		  R0_local=R0_local-R1
		  R(I)=R0_local
	    ENDDO ! I

      RETURN
      END SUBROUTINE RDMU2

!sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
