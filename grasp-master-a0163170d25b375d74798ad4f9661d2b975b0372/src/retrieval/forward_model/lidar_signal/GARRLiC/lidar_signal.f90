! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../../../constants_set/mod_globals.inc"
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE lidar_signal_elastic ( HOBS_km, NVERT, HVP_km,   &
                                        NSD, AVP_norm, EXTA, LRA, &
                                        MVP_norm, EXTM,           &
                                        LS )
      
      USE mod_molecular_scattering, only : LRM
      USE mod_par_inv, only  : KVERTM
      USE mod_par_os,  only  : KSD
      USE mod_stop_report

! forward lidar problem for Single pixel and Single wavelength 
! with Molecular correction (extinction and backscatter)
! HOBS_km is altitude of obserwation above sea level in km (not in actual use)
! NVERT - number of heights for vertical aerosol profile
! HVP_km(NVERT) – altitudes of layers of vertical profiles in descending order
! NSD - number of aerosol fractions 
! APV_norm(NSD,NVERT) - Aerosol Vertival Profile normalized (dimentionless)
! EXTA(NSD) aerosol extinction in 1/km
! LRA - Lidar Ratio of Aerosol in Sr
! MVP_norm molecular vertical profile normalized (dimentionless)
! EXTM - molecular extinction in 1/km
! LS(NVERT) - ouput vector containing lidar measurement calculaions

      IMPLICIT NONE

      INTEGER,                    INTENT(IN)     :: NSD, NVERT ! number of aerosol modes, number of AVP_norm
      REAL, DIMENSION(KVERTM,KSD),INTENT(IN)     :: AVP_norm   ! normalized Aerosol extinction Vertical Profile (AVP)
      REAL, DIMENSION(KVERTM),    INTENT(IN)     :: HVP_km     ! Vertical Profile Heights 
      REAL, DIMENSION(KSD),       INTENT(IN)     :: EXTA, LRA  ! aerosol extinction and lidar ratio
      REAL,                       INTENT(IN)     :: HOBS_km    ! meters obove sea level
      REAL,                       INTENT(IN)     :: EXTM       ! aerosol extinction (scattering)
      REAL, DIMENSION(KVERTM),    INTENT(IN)     :: MVP_norm   ! normalized Molecular extinction Vertical Profile (MVP)
      REAL, DIMENSION(KVERTM),    INTENT(OUT)    :: LS         ! lidar signal
!-- internal variables ------------------------------------------- 
      INTEGER                     :: ISD,IVERT  ! indexes
      REAL, DIMENSION(KSD)        :: AVP_INT    ! AVP_int integrated aerosol profile  
      REAL                        :: LS_INT     ! accumulator for lidar signal norm
      REAL                        :: MVP_INT    ! integrated molecular profile
      REAL                        :: BA,EA      ! aerosol backscatter and extinction 
      REAL                        :: BM,EM      ! molecular backscatter and extinction
      REAL                        :: DH         ! height step for profile interation
!------------------------------------------------------------------------------------------ 
! [EXTA] = 1/km
! [LRA]  = Sr
! [EXTM] = 1/km
! [LRM]  = Sr
!
! [APV_norm]  = 
! [MPV_norm]  = 
 
      IF(HOBS_km .GT. HVP_km(1)) THEN
        WRITE(tmp_message,'(2(a,es9.4),a)') &
        'HOBS_km =',HOBS_km ,' .GT. HVP_km(1) =', HVP_km(1),'  - invalid measurements height'
        G_ERROR(trim(tmp_message))
      ENDIF

!WRITE(*,*) 'NVERT=', NVERT
!WRITE(*,*) 'EXTA=', EXTA
!WRITE(*,*) 'HOBS_km=', HOBS_km

      LS(:)      = 0.0      
      AVP_INT(:) = 0.0
      MVP_INT    = 0.0
      LS_INT     = 0.0

      do IVERT=NVERT-1,1,-1
          DH = HVP_km(IVERT)-HVP_km(IVERT+1)
          AVP_INT(1:NSD) = AVP_INT(1:NSD) + 0.5*(AVP_norm(IVERT,1:NSD)+AVP_norm(IVERT+1,1:NSD))*DH
          MVP_INT        = MVP_INT + 0.5*(MVP_norm(IVERT)+MVP_norm(IVERT+1))*DH
          BA = SUM(EXTA(1:NSD)*AVP_norm(IVERT,1:NSD)/LRA(1:NSD))
!          BA = 0.332739*AVP_norm(IVERT,1)/69.5112
          EA = SUM(EXTA(1:NSD)*AVP_INT(1:NSD))
!          EA = 0.332739*AVP_INT(1)
!          write(*,*) 'ivert, EXTA', IVERT,0.332739*AVP_norm(IVERT,1)
          BM = EXTM*MVP_norm(IVERT)/LRM
          EM = EXTM*MVP_INT
          LS(IVERT) = (BA+BM)*EXP(-2.*(EA+EM))
      enddo ! IVERT

!write(*,*) 'EA=',EA,'  EXTA=',SUM(EXTA(1:NSD)),'  AVP_INT(1:NSD): ',AVP_INT(1:NSD)
!write(*,*) 'EM=',EM,'  EXTM=',EM,              '  MVP_INT       : ',MVP_INT

! Lidar Signal normalization
      DO IVERT=NVERT-2,1,-1
        DH = HVP_km(IVERT)-HVP_km(IVERT+1)
        LS_INT = LS_INT + 0.5*(LS(IVERT)+LS(IVERT+1))*DH 
      ENDDO ! IVERT
      LS(1:NVERT-1) = LS(1:NVERT-1)/LS_INT

      IF(ANY(LS(1:NVERT-1) .LE. 0. .OR. IsNaN(LS(1:NVERT-1)))) THEN
        GOTO 33
        WRITE(*,*) 'ERROR IN lidar_signal_elastic_v0:'
        WRITE(*,*) 'INVALID OUTPUT'
        WRITE(*,*) 'DEBUG INFORMATION:'
        WRITE(*,*) 'LS:'
        WRITE(*,'(10e14.5)') LS(1:NVERT-1)
        WRITE(*,*) 'AVP_INT: ', AVP_INT(1:NSD)
        WRITE(*,*) 'MVP_INT: ', MVP_INT
        WRITE(*,*) 'EA=',EA,'  BA=',BA
        WRITE(*,*) 'EM=',EM,'  BM=',BM  
        WRITE(*,*) 'LS_INT=',LS_INT 
        WRITE(*,*) 'SUBROUTINE INPUTS:'
        WRITE(*,*) 'LRA=',LRA,' LRM=',LRM
        DO ISD=1,NSD
          WRITE(*,*) 'AVP_norm: ISD=',ISD
          WRITE(*,'(10e14.5)') AVP_norm(1:NVERT,ISD) 
        ENDDO ! ISD   
        WRITE(*,*) 'HVP in km:'
        WRITE(*,'(10f14.5)') HVP_km(1:NVERT)     
        WRITE(*,*) 'MVP_norm:'
        WRITE(*,'(10e14.5)') MVP_norm(1:NVERT)
        WRITE(*,*) 'NVERT=', NVERT
        WRITE(*,*) 'EXTA=',SUM(EXTA(1:NSD)),'  EXTA(1:NSD)=', EXTA(1:NSD)
        WRITE(*,*) 'HOBS_km=', HOBS_km
        !WRITE(*,*) 'HVP_km AVP_norm(1:NSD) : '
        !DO IVERT=1,NVERT
          !WRITE(*,*) HVP_km(IVERT),AVP_norm(IVERT,1:NSD)
        !ENDDO
33      CONTINUE
        write(tmp_message,'(a)') 'invalid output, negative or NaN'
        G_ERROR(trim(tmp_message))
      ENDIF ! ANY(LS(1:NVERT) .LT. 0. .OR. ....
                        
      RETURN
      END SUBROUTINE lidar_signal_elastic

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
 
      subroutine profile_normalization (NH, H, PROF_in, PROF_out)

      implicit none
!     INPUT
      integer,intent(in)             :: NH                    ! number of heights and points in profile
      real,intent(in)                :: H(NH),PROF_in(NH)     ! height and profile that be normalized
!     LOCAL
      integer                        :: ih
      real                           :: norma
!     OUTPUT
      real,intent(out)               :: PROF_out(NH)          ! normalized profile
!     ----------------------------------------------------------

      PROF_out(:) = 0.0
     
         norma=0
         do ih=2,NH
            norma = norma + 0.5 * (PROF_in(ih-1)+PROF_in(ih))  &
                                                       * (H(ih-1)-H(ih))
         enddo

         PROF_out = PROF_in/norma
      

!!begin       checking the normalization 
!         norma=0
!         do ih=2,NH
!            norma = norma + 0.5*(PROF_out(ih-1)+PROF_out(ih))  &
!                                                        *(H(ih-1)-H(ih))
!         enddo
!         write(*,*) norma,'  - if equal 1 than it is good normalization'
!!end

      end subroutine profile_normalization

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
