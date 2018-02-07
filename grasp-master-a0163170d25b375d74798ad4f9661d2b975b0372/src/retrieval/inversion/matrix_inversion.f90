! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!***********************************************************************

      SUBROUTINE ITERQ( IM, KN, UFS, FFS, QS )
   
!C******************************************************
!C  THIS SUBROUTINE IMPLEMENTS q-th (general) ITERATIONS
!C******************************************************
!C  INPUT:
!C         IM    I   IM=1 - solving linear system by steepest decent method
!C                   IM=2 - linear system can be solved by SVD standard routine
!C         KN    I        - number of parameters
!C         UFS   R(KN,KN) - "Fisher matrix" normalized by
!C                          maximum diagonal element
!C         FFS   R(KM)    - "gradient" vector normalized by
!C                          maximum diagonal element    
!C  OUTPUT:
!C         QS    R(KN)    - solution
!C         RES   R        - residual produced by DAP
!C******************************************************
      IMPLICIT NONE
! ------------------------------------------------------------
! IN :
      INTEGER              ,INTENT(IN)  :: IM,KN
      REAL,DIMENSION(KN,KN),INTENT(IN)  :: UFS
      REAL,DIMENSION(KN)   ,INTENT(IN)  :: FFS

! ------------------------------------------------------------
! OUT :
      REAL,DIMENSION(KN)   ,INTENT(OUT) :: QS
! ------------------------------------------------------------
! LOCAL :
      INTEGER                 :: I,info       
      CHARACTER*64            :: ERR=' '
      REAL*8,DIMENSION(KN,KN) :: UFM
      REAL*8,DIMENSION(KN)    :: FFM,QM
      REAL*8                  :: DT,DIF,DRES
      REAL                    :: AAM,RES
! ------------------------------------------------------------
!write(*,*) 'in ITERQ'      
      QS(:) = 0.0
      AAM = 0.0
      AAM = 1.0
!      DO J=1,KN
!OD      IF(AAM .LT. UFS(J,J)) AAM = UFS(J,J)
!      ENDDO

      FFM(1:KN) = FFS(1:KN)/AAM
      UFM(1:KN,1:KN) = UFS(1:KN,1:KN)/AAM

      SELECT CASE(IM)
      CASE(1) 
         DO I=1,KN        
         QM(I) = FFM(I)/SUM(ABS(UFM(I,1:KN)))
         QS(I) = QM(I)
         ENDDO
!      CASE(2)
! system can be solved using the singular value decomposition (SVD)
! QM(1:KN) - solution
!         QS(1:KN)=QM(1:KN)
      CASE DEFAULT
         WRITE(*,*) 'IM =',IM,'  unknown value'
         STOP 'STOP in ITERQ'
      END SELECT

      DRES = 0.0
      DO I=1,KN
      DIF = SUM(UFM(I,1:KN)*QM(1:KN))-FFM(I)
      DRES = DRES+DIF*DIF
!OD      WRITE(*,*) DIF,FFM(I),I,' DIF,FFM(I),I'
      ENDDO ! I
      DRES = SQRT(DRES/KN)
      RES = DRES
!OD      WRITE(*,*) RES," residual of matrix inversion"
      RETURN
      END SUBROUTINE ITERQ
      
!***********************************************************************

