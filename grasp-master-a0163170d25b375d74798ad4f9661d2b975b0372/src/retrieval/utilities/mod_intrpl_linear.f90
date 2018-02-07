! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

        !> @file mod_intrpl_linear.f90
        !> Module contains linear interpolation routines
        !>
        !> @author Oleg Dubovik
        !>

      module mod_intrpl_linear

      contains
      
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      FUNCTION  LINEAR( X, Y, M, X1 )
! ----------------------------------------------------------------------
!       Linear interpolation function
! ----------------------------------------------------------------------
      IMPLICIT  REAL   (A-H, O-Z)
      REAL, INTENT(IN) :: X( * ), Y( * )
      REAL, INTENT(IN) :: X1
      INTEGER, INTENT(IN) :: M
      REAL :: LINEAR
! ----------------------------------------------------------------------
      INTEGER :: N
! ----------------------------------------------------------------------
      LINEAR = 0.0

      IF(X(2) .GT. X(1)) THEN
        IF(X1 .LT. X(1)) THEN
          LINEAR = Y(1)   + ( X1-X(1) )*( Y(2)-Y(1) )/( X(2)-X(1) )
        ELSE IF ( X1 .GT. X(M) ) THEN
          !LINEAR = Y(M)   + ( X1-X(M) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
          LINEAR = Y(M-1) + ( X1-X(M-1) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
        ELSE
          DO N=2,M
            IF(X1 .GE. X(N-1) .AND. X1 .LE. X(N)) &
            LINEAR = Y(N-1) + ( X1-X(N-1) )*( Y(N)-Y(N-1) )/( X(N)-X(N-1) )
!C*****************************
            IF(X1 .EQ. X(N-1)) LINEAR = Y(N-1)
            IF(X1 .EQ. X(N))   LINEAR = Y(N)
!C*****************************
          ENDDO ! N
        END IF ! X1 .LT. X(1)
      ELSE
        IF(X1 .GT. X(1)) THEN
          LINEAR = Y(1)   + ( X1-X(1) )*( Y(2)-Y(1) )/( X(2)-X(1) )
        ELSE IF ( X1 .LT. X(M) ) THEN
          !LINEAR = Y(M)   + ( X1-X(M) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
          LINEAR = Y(M-1) + ( X1-X(M-1) )*( Y(M)-Y(M-1) )/( X(M)-X(M-1) )
        ELSE
          IF( M .EQ. 2 .AND. (X(2) .EQ. X(1)) ) THEN
            IF(Y(2) .NE. Y(1)) THEN
              WRITE(*,'(a,i0,4(2x,a,e12.5))') &
              'M = ',M,'X(2) =',X(2),'X(1) =',X(1),'Y(2) =',Y(2),'Y(1) =',Y(1)
              STOP 'STOP in FUNCTION LINEAR (interpolation)'
            ENDIF
            LINEAR = Y(1)
          ELSE
            DO N=2,M
              IF(X1 .LE. X(N-1) .AND. X1 .GE. X(N)) &
              LINEAR = Y(N-1) + ( X1-X(N-1) )*( Y(N)-Y(N-1) )/( X(N)-X(N-1) )
!C*****************************
              IF(X1 .EQ. X(N-1)) LINEAR = Y(N-1)
              IF(X1 .EQ. X(N))   LINEAR = Y(N)
!C*****************************
            ENDDO ! N
          ENDIF ! M .EQ. 2 .AND.
        ENDIF ! X1 .GT. X(1)
      ENDIF ! X(2) .GT. X(1)

      RETURN
      END FUNCTION LINEAR

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      FUNCTION  LINEAR_LN (X, Y, M, X1 )
! ----------------------------------------------------------------------
!       Linear interpolation function.                                
! ----------------------------------------------------------------------
      IMPLICIT  REAL   (A-H, O-Z)
      REAL, INTENT(IN) :: X( * ), Y( * )
      REAL, INTENT(IN) :: X1
      INTEGER, INTENT(IN) :: M
      REAL ::	LINEAR_LN
! ----------------------------------------------------------------------
      INTEGER :: N
! ----------------------------------------------------------------------
      LINEAR_LN = 0.0

      IF(X(2) .GT. X(1)) THEN
        IF(X1 .LT. X(1)) THEN
          LINEAR_LN = LOG(Y(1))   + ( X1-X(1) )*( LOG(Y(2))-LOG(Y(1)) )/( X(2)-X(1) )
          LINEAR_LN = EXP(LINEAR_LN)
        ELSE IF ( X1 .GT. X(M) )  THEN
          !LINEAR_LN = LOG(Y(M))   + ( X1-X(M) )*( LOG(Y(M))-LOG(Y(M-1)) )/( X(M)-X(M-1) )
          LINEAR_LN = LOG(Y(M-1)) + ( X1-X(M-1) )*( LOG(Y(M))-LOG(Y(M-1)) )/( X(M)-X(M-1) )
          LINEAR_LN = EXP(LINEAR_LN)
        ELSE
          DO N=2,M
            IF(X1 .GE. X(N-1) .AND. X1 .LE. X(N)) THEN
              LINEAR_LN = LOG(Y(N-1)) + ( X1-X(N-1) )*( LOG(Y(N))-LOG(Y(N-1)) )/( X(N)-X(N-1) )
              LINEAR_LN = EXP(LINEAR_LN)
            ENDIF ! X1 .GE. X(N-1) .AND. X1 .LE. X(N)
!C*****************************
            IF(X1.EQ.X(N-1)) LINEAR_LN = Y(N-1)
            IF(X1.EQ.X(N))   LINEAR_LN = Y(N)
!C*****************************
          ENDDO ! N
        ENDIF !  X1.LT.X(1)
      ELSE 
        IF(X1 .GT. X(1)) THEN
          LINEAR_LN = LOG(Y(1))   + ( X1-X(1) )*( LOG(Y(2))-LOG(Y(1)) )/( X(2)-X(1) )
          LINEAR_LN = EXP(LINEAR_LN)
        ELSE IF(X1 .LT. X(M)) THEN
          !LINEAR_LN = LOG(Y(M))   + ( X1-X(M) )*( LOG(Y(M))-LOG(Y(M-1)) )/( X(M)-X(M-1) )
          LINEAR_LN = LOG(Y(M-1)) + ( X1-X(M-1) )*( LOG(Y(M))-LOG(Y(M-1)) )/( X(M)-X(M-1) )
          LINEAR_LN = EXP(LINEAR_LN)
        ELSE
          IF( M .EQ. 2 .AND. (X(2) .EQ. X(1)) ) THEN
            IF(Y(2) .NE. Y(1)) THEN
              WRITE(*,'(a,i0,4(2x,a,e12.5))') &
              'M = ',M,'X(2) =',X(2),'X(1) =',X(1),'Y(2) =',Y(2),'Y(1) =',Y(1)
              STOP 'STOP in FUNCTION LINEAR_LN (interpolation)'
            ENDIF
            LINEAR_LN = Y(1)
          ELSE
            DO N = 2, M
              IF(X1 .LE. X(N-1) .AND. X1 .GE. X(N)) THEN
                LINEAR_LN = LOG(Y(N-1)) + ( X1-X(N-1) )*( LOG(Y(N))-LOG(Y(N-1)) )/( X(N)-X(N-1) )
                LINEAR_LN = EXP(LINEAR_LN)
              ENDIF ! N
!C*****************************
              IF(X1 .EQ. X(N-1)) LINEAR_LN = Y(N-1)
              IF(X1 .EQ. X(N))   LINEAR_LN = Y(N)
!C*****************************
            ENDDO ! N
          ENDIF ! M .EQ. 2 .AND.
        ENDIF ! X1 .GT. X(1)
      ENDIF ! X(2) .GT. X(1)
      
      RETURN
      END FUNCTION  LINEAR_LN 

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      end module mod_intrpl_linear
