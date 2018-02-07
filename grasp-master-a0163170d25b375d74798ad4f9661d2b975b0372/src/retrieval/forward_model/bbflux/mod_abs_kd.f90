! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

      MODULE mod_abs_kd
!XH   this module contains the constants required for k-distribution 
!XH   subroutines
      IMPLICIT NONE
!
      INTEGER, PARAMETER ::	NCORPS	= 7
      INTEGER, PARAMETER ::	NLEVEL	= 50

!XH   for k-distribution
      INTEGER, PARAMETER :: MAXKD = 15
      INTEGER, PARAMETER :: MAXP  = 23
!XH   maximum number of wavelengths affected by absorbing gases
      INTEGER, PARAMETER :: MAXNU = 208

      TYPE DATA_ABS_GS_WL
!XH      describes data of absorption for single wavelength

!XH      wavelength corresponding to the current dataset
         REAL    :: WVL
!XH      tau for continuum
         REAL, DIMENSION(NLEVEL-1)              :: COEFCN
!XH      number of elements in AIKD
         INTEGER, DIMENSION(NCORPS)             :: NEXP
!XH      weights for K-distribution
         REAL, DIMENSION(MAXKD,NCORPS)          :: AIKD
!XH      tau for vertical layers
         REAL, DIMENSION(NLEVEL-1,MAXKD,NCORPS) :: COEFKD
      END TYPE DATA_ABS_GS_WL

      TYPE DATA_ABS
         INTEGER :: NWL
!XH      number of levels
         INTEGER :: NLV
!XH      number of layers
         INTEGER :: NLY
!XH      solar zenith angle
         REAL    :: SZA
!XH      altitudes of the levels stored in CEOFKD and AIKD
         REAL,    DIMENSION(NLEVEL) :: HLV
!XH      altitudes of the layers stored in CEOFKD and AIKD
         REAL,    DIMENSION(NLEVEL-1) :: HLY
!XH      set of data for absorbing gases at all wavelengths
         TYPE(DATA_ABS_GS_WL), DIMENSION(MAXNU) :: ABS_GS_WL
      END TYPE DATA_ABS

!XH   for solar spectrum
      INTEGER, PARAMETER :: NBNU   = 208
!XH   default solar constant 1361.1 W/m2
      REAL,    PARAMETER :: CSLR   = 1361.1
      REAL,    PARAMETER, DIMENSION(NBNU) :: FSLR =                     &
        (/0.144706953E+1,0.154556000E+1,0.166134000E+1,0.177756953E+1,  &
          0.188813961E+1,0.200285006E+1,0.212966967E+1,0.224530005E+1,  &
          0.237714005E+1,0.249806929E+1,0.263052917E+1,0.275070024E+1,  &
          0.288664985E+1,0.300513029E+1,0.313730979E+1,0.326431060E+1,  &
          0.339966989E+1,0.353551888E+1,0.368415022E+1,0.383008051E+1,  &
          0.398967957E+1,0.409932947E+1,0.422120094E+1,0.438673067E+1,  &
          0.458120966E+1,0.472347927E+1,0.486300135E+1,0.503059196E+1,  &
          0.515899992E+1,0.530883837E+1,0.549096966E+1,0.562026978E+1,  &
          0.575858927E+1,0.593737030E+1,0.606163931E+1,0.618625021E+1,  &
          0.631413841E+1,0.637730074E+1,0.650247765E+1,0.655056047E+1,  &
          0.660543013E+1,0.664381886E+1,0.665074015E+1,0.669554043E+1,  &
          0.672023916E+1,0.675603008E+1,0.679007864E+1,0.680182838E+1,  &
          0.681483030E+1,0.686350060E+1,0.694996071E+1,0.701644993E+1,  &
          0.706062126E+1,0.717784834E+1,0.718845892E+1,0.719963169E+1,  &
          0.724680948E+1,0.725563002E+1,0.726125908E+1,0.722606230E+1,  &
          0.724480963E+1,0.725711107E+1,0.725141048E+1,0.725760984E+1,  &
          0.726437044E+1,0.728695869E+1,0.721799850E+1,0.727764940E+1,  &
          0.725099039E+1,0.726626015E+1,0.730128956E+1,0.733327818E+1,  &
          0.738649940E+1,0.740124130E+1,0.740665102E+1,0.746838951E+1,  &
          0.746516037E+1,0.745713091E+1,0.747367096E+1,0.746265125E+1,  &
          0.751808071E+1,0.747275019E+1,0.752733850E+1,0.743650961E+1,  &
          0.749948931E+1,0.746387911E+1,0.750858068E+1,0.747383022E+1,  &
          0.745444012E+1,0.745608997E+1,0.712881947E+1,0.733367062E+1,  &
          0.698421860E+1,0.734103918E+1,0.732415295E+1,0.734712982E+1,  &
          0.726125860E+1,0.735255003E+1,0.729576063E+1,0.727475023E+1,  &
          0.726430845E+1,0.727113914E+1,0.726604986E+1,0.725602102E+1,  &
          0.713171005E+1,0.714738798E+1,0.714892960E+1,0.717704105E+1,  &
          0.711783028E+1,0.704155874E+1,0.705714130E+1,0.710239029E+1,  &
          0.704902935E+1,0.702091885E+1,0.693199062E+1,0.701413012E+1,  &
          0.701186895E+1,0.689467001E+1,0.692974901E+1,0.276803722E+2,  &
          0.274372444E+2,0.264879265E+2,0.266201496E+2,0.258811131E+2,  &
          0.257061520E+2,0.249477959E+2,0.243669853E+2,0.231927319E+2,  &
          0.227000523E+2,0.215444107E+2,0.209714794E+2,0.194459324E+2,  &
          0.194949493E+2,0.192061272E+2,0.181654053E+2,0.184181347E+2,  &
          0.174161949E+2,0.171397915E+2,0.164732265E+2,0.144773045E+2,  &
          0.130930386E+2,0.112844963E+2,0.121915741E+2,0.119985199E+2,  &
          0.109732065E+2,0.978662777E+1,0.656143045E+1,0.619820356E+1,  &
          0.558374119E+1,0.633936501E+1,0.650426149E+1,0.637373066E+1,  &
          0.483241129E+1,0.536343050E+1,0.536082268E+1,0.459261894E+1,  &
          0.402703285E+1,0.387259078E+1,0.446363926E+1,0.423248148E+1,  &
          0.306956339E+1,0.312597942E+1,0.293557811E+1,0.279447389E+1,  &
          0.257196021E+1,0.217034912E+1,0.181004930E+1,0.190915537E+1,  &
          0.203536296E+1,0.159504116E+1,0.918867528E+0,0.975803971E+0,  &
          0.392557293E+0,0.703308761E+0,0.575877190E+0,0.785403073E+0,  &
          0.766851902E+0,0.676903844E+0,0.273765117E+0,0.319081813E+0,  &
          0.266469538E+0,0.135761753E+0,0.131270260E+0,0.127161711E+0,  &
          0.124792971E+0,0.157211468E+0,0.119562685E+0,0.945537537E-1,  &
          0.105255023E+0,0.948318914E-1,0.106592536E+0,0.979281366E-1,  &
          0.961348340E-1,0.850596726E-1,0.124257594E+0,0.879886970E-1,  &
          0.906372145E-1,0.772183314E-1,0.610966384E-1,0.727205649E-1,  &
          0.563923009E-1,0.553333387E-1,0.333273299E-1,0.221149456E-1,  &
          0.190008096E-1,0.170744043E-1,0.143835442E-1,0.133104045E-1/)

      CONTAINS

      SUBROUTINE GET_ABS_DATA(PATH,IABS,ICONT,IATM,PSRF,UH2O,UO3,       &
                              ABS_DATA)
!
      CHARACTER(*),                 INTENT(IN) :: PATH
      INTEGER,                      INTENT(IN) :: IATM
      REAL,                         INTENT(IN) :: PSRF, UH2O, UO3
      INTEGER,   DIMENSION(NCORPS), INTENT(IN) :: IABS, ICONT
      TYPE(DATA_ABS),               INTENT(OUT):: ABS_DATA
!
      INTEGER, DIMENSION(NCORPS,MAXNU)         :: NEXP
      REAL,  DIMENSION(NCORPS,MAXNU,MAXKD)     :: AL
      REAL,  DIMENSION(NCORPS,MAXNU,MAXP,MAXKD):: XKI, AA, BB
      REAL,  DIMENSION(NLEVEL,13)              :: DONUSER
      REAL,  DIMENSION(NLEVEL)                 :: DENS, T, P, ALT
      REAL,  DIMENSION(NLEVEL-1)               :: TAUCONT
      REAL,  DIMENSION(NCORPS,NLEVEL)          :: RO
!
      REAL    :: XK, AIK, PRS, TMP, WNM
      INTEGER :: INU, J, K, IK
!
!.......................................................................
!     absorbants
!     lines :
!        1-H2O  2-CO2  3-O3  4-N2O  5-CO  6-CH4  7-O2
!     continuum :
!        1-H2O  2-CO2  3-O3  4-N2   5-O2  6-NO2  7-SO2
!.......................................................................
!XH   Atmopsheric Profile
      CALL DATATM(PATH,IATM,PSRF,UH2O,UO3,                              &
                  DONUSER,ALT,P,T,RO,DENS)
!XH   read absorption coefficients for H2O, CO2(other gases), and O3
      DO K=1, 3
         CALL READ_ABSN(PATH,K,IABS(K),                                 &
                        NEXP,AL,XKI,AA,BB)
      END DO
!
      ABS_DATA%NWL=NBNU
      ABS_DATA%NLV=NLEVEL
      ABS_DATA%HLV(1:NLEVEL)=ALT(1:NLEVEL)
      ABS_DATA%NLY=NLEVEL-1
      ABS_DATA%HLY(1:NLEVEL-1)=0.5*(ALT(1:NLEVEL-1)+ALT(2:NLEVEL))
!
      DO INU = 1, NBNU
         IF (INU .LT. 120) THEN
            WNM= 2550.0+100.0*(INU-1  )
         ELSE
            WNM=14600.0+400.0*(INU-120)
         END IF
         CALL CONTINUUMNEW(REAL(WNM,8),NLEVEL,ICONT,DONUSER,            &
                           TAUCONT)
         ABS_DATA%ABS_GS_WL(INU)%WVL   =1.0E+04/WNM
         ABS_DATA%ABS_GS_WL(INU)%NEXP  =NEXP(:,INU)
         ABS_DATA%ABS_GS_WL(INU)%COEFCN=TAUCONT
         DO K = 1, 3
            DO IK = 1, NEXP(K,INU)
               ABS_DATA%ABS_GS_WL(INU)%AIKD(IK,K)=AL(K,INU,IK)
               DO J=1, NLEVEL-1
                  PRS = 0.5*(P(j)+P(j+1))
                  TMP = 0.5*(T(j)+T(j+1))
!                 k-distribution coefficients:
!                 al(k,inu,ik) : weight of the k-distribution for
!                                the absorber k, spectral interval
!                                inu, and bin ik.
!                 xk           : absorption coefficent for the layer
!                                j and the spectral interval inu
                  XK  = 0.0
                  IF (NEXP(K,INU) .GT. 1)                               &
                     CALL COEFABS(K,INU,IK,PRS,TMP,XKI,AA,BB,           &
                                  XK)
                  ABS_DATA%ABS_GS_WL(INU)%COEFKD(J,IK,K)=XK*RO(K,J)
               END DO
            END DO
         END DO
      END DO
!
      RETURN
      END SUBROUTINE GET_ABS_DATA

      SUBROUTINE ABS_GAS_PROFILE ( ABS_DATA,       &  ! IN
                                   IW,             &
                                   I1, I2, I3,     &
                                   !NGAS,            &
                                   NLV_MAX,        &
                                   NLV, HLV, PROF, &  ! OUT
                                   EXT             &
                                 )
      USE mod_par_inv,   ONLY : KVERTM
      IMPLICIT NONE
!
      INTEGER,                   INTENT(IN) :: IW!,NGAS
      INTEGER,                   INTENT(IN) :: I1,I2,I3
!XH   data of gas absorption
      TYPE(DATA_ABS),            INTENT(IN) :: ABS_DATA

      INTEGER,                  INTENT(IN)  :: NLV_MAX
      INTEGER,                  INTENT(INOUT) :: NLV
      REAL, DIMENSION(NLV_MAX), INTENT(INOUT) :: HLV, PROF
      REAL,                     INTENT(OUT) :: EXT
!
      INTEGER :: NLY, I
      REAL, DIMENSION(NLV_MAX) :: HLY, TAU
!
!XH   number of levels
      NLV = ABS_DATA%NLV
      HLV(NLV:1:-1) = ABS_DATA%HLV(1:NLV)
      IF (NLV_MAX .LT. NLV) THEN
         WRITE(*,*) 'ABS_GAS_DATA in mod_abs_kd : ',                    &
                    'number of levels', NLV, ' > ',                     &
                    'maximum of levels', NLV_MAX
         STOP
      END IF
!XH   absorbing gas optical thickness
      NLY = ABS_DATA%NLY
      HLY(NLY:1:-1) = ABS_DATA%HLY(1:NLY)
      TAU(NLY:1:-1) = ABS_DATA%ABS_GS_WL(IW)%COEFCN(1:NLY)              &
                     +ABS_DATA%ABS_GS_WL(IW)%COEFKD(1:NLY,I1,1)         &
                     +ABS_DATA%ABS_GS_WL(IW)%COEFKD(1:NLY,I2,2)         &
                     +ABS_DATA%ABS_GS_WL(IW)%COEFKD(1:NLY,I3,3)
!XH   calculate gas profile, assuming prof = 0.0 at TOA
      PROF(1) = 0.0
      DO I = 2, NLV
!XH      TAU(I-1) = 0.5*[PROF(I-1)+PROF(I)]*[HLV(I-1)-HLV(I)]
         PROF(I) = 2.0*TAU(I-1)/(HLV(I-1)-HLV(I))-PROF(I-1)
      END DO
!XH   total extinction of all absorbing gases for all layers
      EXT = SUM(TAU(1:NLY))
!
      RETURN
      END SUBROUTINE ABS_GAS_PROFILE

      SUBROUTINE ABS_GAS_INT(                                           &
                              IW,I1,I2,I3,                              &
                              ABS_DATA,                                 &
                              NBV,                                      &
                              SL,SQ,SU,SL_TMP,SQ_TMP,SU_TMP,            &
                              UFG,DFG,UFG_TMP,DFG_TMP,                  &
                              NLV,                                      &
                              UFX0,DFX0,UFX,DFX,                        &
                              UFX0_TMP,DFX0_TMP,UFX_TMP,DFX_TMP         &
                            )
      USE mod_par_OS,    ONLY : NBVM,KNT
      IMPLICIT NONE
!
      INTEGER,                   INTENT(IN)    :: IW
      INTEGER,                   INTENT(IN)    :: I1,I2,I3
!XH   data of gas absorption
      TYPE(DATA_ABS),            INTENT(IN)    :: ABS_DATA

      INTEGER,                   INTENT(IN),    OPTIONAL :: NLV,NBV
      REAL,                      INTENT(IN),    OPTIONAL :: UFG,DFG
      REAL, DIMENSION(KNT),      INTENT(IN),    OPTIONAL :: UFX,DFX
      REAL, DIMENSION(KNT),      INTENT(IN),    OPTIONAL :: UFX0,DFX0
      REAL, DIMENSION(2*NBVM),   INTENT(IN),    OPTIONAL :: SL,SQ,SU
!
      REAL,                      INTENT(INOUT), OPTIONAL :: UFG_TMP
      REAL,                      INTENT(INOUT), OPTIONAL :: DFG_TMP
      REAL, DIMENSION(KNT),      INTENT(INOUT), OPTIONAL :: UFX_TMP
      REAL, DIMENSION(KNT),      INTENT(INOUT), OPTIONAL :: DFX_TMP
      REAL, DIMENSION(KNT),      INTENT(INOUT), OPTIONAL :: UFX0_TMP
      REAL, DIMENSION(KNT),      INTENT(INOUT), OPTIONAL :: DFX0_TMP
      REAL, DIMENSION(2*NBVM),   INTENT(INOUT), OPTIONAL :: SL_TMP
      REAL, DIMENSION(2*NBVM),   INTENT(INOUT), OPTIONAL :: SQ_TMP
      REAL, DIMENSION(2*NBVM),   INTENT(INOUT), OPTIONAL :: SU_TMP
!
      REAL :: WT
!
!XH   weight from K-distribution for current wavelength
      WT = ABS_DATA%ABS_GS_WL(IW)%AIKD(I1,1)                            &
          *ABS_DATA%ABS_GS_WL(IW)%AIKD(I2,2)                            &
          *ABS_DATA%ABS_GS_WL(IW)%AIKD(I3,3)
      IF (PRESENT(NBV)) THEN
         UFG_TMP=UFG_TMP+UFG*WT
         DFG_TMP=DFG_TMP+DFG*WT
!XH      integration for Stokes vector
         SL_TMP(1:NBV)=SL_TMP(1:NBV)+SL(1:NBV)*WT
         IF (PRESENT(SQ)) THEN
            SQ_TMP(1:NBV)=SQ_TMP(1:NBV)+SQ(1:NBV)*WT
            SU_TMP(1:NBV)=SU_TMP(1:NBV)+SU(1:NBV)*WT
         END IF
      END IF
      IF (PRESENT(NLV)) THEN
!XH      integration for flux
         UFX0_TMP(1:NLV)=UFX0_TMP(1:NLV)+UFX0(1:NLV)*WT
         DFX0_TMP(1:NLV)=DFX0_TMP(1:NLV)+DFX0(1:NLV)*WT
         UFX_TMP(1:NLV) = UFX_TMP(1:NLV)+ UFX(1:NLV)*WT
         DFX_TMP(1:NLV) = DFX_TMP(1:NLV)+ DFX(1:NLV)*WT
      END IF
!
      RETURN
      END SUBROUTINE ABS_GAS_INT
!
      END MODULE mod_abs_kd
