! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!
! Module for BRDF/BPDF and BRM (Bidirectional Reflection Matrix) calculation
! (Contact P.Litvinov for further details 
! (Pavel.Litvinov@univ-lille1.fr;
! PVLitvinov@mail.ru)

! Physical description can be found in the papers:

! For land reflection models:

! 1.    Pavel Litvinov, Otto Hasekamp, Brian Cairns. 
!  Models for surface reflection of radiance and polarized radiance: 
! Comparison with airborne multi-angle photopolarimetric measurements and implications 
! for modeling top-of-atmosphere measurements. Remote Sensing of Environment (2011), 115, 781–792

! 2.	Pavel Litvinov, Otto Hasekamp, Brian Cairns, Michael Mishchenko. 
! Semi-empirical BRDF and BPDF models applied to the problem of aerosol retrievals over land: 
! testing an airborne data and implications for modeling of top-of-atmosphere measurements.
! In book: Polarimetric Detection, Characterization, and Remote Sensing. Springer, 2011.
!
! 3. Pavel Litvinov, Otto Hasekamp, Oleg Dubovik, Brian Cairns. 
! Model for land surface reflectance treatment: 
! Physical derivation, application for bare soil and evaluation on airborne 
! and satellite measurements. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.

! For ocean reflection models:

! 1. Mishchenko MI, Travis LD. Satellite retrieval of aerosol properties over the ocean using
! polarization as well as intensity of reflected sunlight. J. Geophys. Res., 1997;102:16989-
! 17013.
!
! 2. Tsang L, Kong JA, Ding K-H. Scattering of electromagnetic waves. Advanced topics.
! Wiley, New York, 2001. 
! 
! 3. 6.	Hasekamp, O., P. Litvinov, and A. Butz, 
! Aerosol properties over the ocean from PARASOL multi-angle photopolarimetric measurements, 
! J. Geophys. Res. (2011), doi:10.1029/2010JD015469, D14204.
!

! Ocean:
! 1). i_BRDF=0: isotripic case
! 2). i_BRDF=1: anisotropic case
!
! Land:
! Subroutine BRM chooses the BRM model according to the flags:
! 1). i_BRDF=0: RPV model in combination with BPDF model (see flag description for
!              i_BPDF). Function "BRDF_RPV"
! 2). i_BRDF=1: Ross-Li model  in combination with BPDF model (see flag description for
!              i_BPDF). Subroutine "LiSp_Ross"
! 3). i_BRDF=2: new physical model based on (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!           Here just intensity is calculated. BPDF model is chosen with flag i_BPDF.
!            Subroutine "Refl_LitvBRM_R11"
! 4). i_BRDF=9,iBPDF=9: new physical model (Litvinov et al. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!                   Subroutine  Refl_mat_4p(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,m_ind,Rm) for reflection matrix

! Surface polarization
! 1). i_BPDF=-1: no surface polarization
! 2). i_BPDF=0:  Maignan et al. BPDF model 
! 3). i_BPDF=1: BPDF model developted at SRON (Litinov et all.  Remote Sensing of Environment (2011), 115, 781–792;
!                                         Litinov et all.In book: Polarimetric Detection, Characterization, and Remote Sensing. Springer, 2011.)
! 4). i_BPDF=9: new physical model (Litvinov et all.) always with iBRDF=9. Subroutine  Refl_mat_4p


Module Mod_BRM
use mod_par_OS, only: NN0,NF0
implicit none
!      Integer, parameter:: NN0=40,NF0=40
      integer, parameter:: NN0_my=NN0, & ! NG0_my=21
                           M_max_max=NF0, nph=1, n_phi_max=nph*M_max_max !n_tet_max=100, M_max_max=100
      Integer,save:: i_gauleg=0
      Real(8),save:: Gw_p_old(1:n_phi_max),phi_old(1:n_phi_max)
!------------------------------------------------------
      Real(8), parameter :: dpi=3.141592653589793d0 !dpi=3.14159265358979323846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
      Real(8), parameter :: d180_pi=180.d0/dpi
      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
      Integer, parameter :: n_par_p=4, n_par_max=6
      Real(8) :: x_vect_p(1:n_par_p)
      Real(8) :: x_vect(1:n_par_max)
!------------------------------------------------------
!XH   Parameter related to the spectral dependent refractive index
!XH   IWAVE    : if IWAVE is true, then spectral dependent surface properties will be used; Otherwise, use refractive index 1.33 for ocean and 1.5 for land
      LOGICAL, PARAMETER :: IWAVE = .FALSE.
!------------------------------------------------------
      
Contains

REAL(8) FUNCTION BRDF_SNOW(MUS,MUV,PHI,D,M,WV)
!XH   Snow BRDF model - A. A. Kokhanovsky and F.-M. Breon, IEEE  Geosci. Remote Sens. Lett. 9, 5, 2012
!     MUS: cosine of solar zenith angle
!     MUV: cosine of viewing zenith angle
!     PHI: relative azimuth angle based on the POLDER convention
!     D  : average optical diameter of snow grains in micron
!     M  : proportional to mass concentration of pollutants in snow
!     WV : wavelength in micron
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: MUS,MUV,PHI,D,M,WV
      INTEGER             :: IWV
      REAL(8) :: THETA, PTHETA, ALPHA, KMUS, KMUV, R0, CHI
!XH   imaginary part of refractive index of ice, given by S. G. Warren and R. E. Brandt, J. Geophys. Res. 113, D14220, 2008
      INTEGER, PARAMETER :: NWV = 144
      REAL(8), PARAMETER, DIMENSION(NWV) :: WVS = (/1.990E-001,2.010E-001,2.019E-001,2.100E-001,2.500E-001,  &
                                                    3.000E-001,3.500E-001,3.900E-001,4.000E-001,4.100E-001,  &
                                                    4.200E-001,4.300E-001,4.400E-001,4.500E-001,4.600E-001,  &
                                                    4.700E-001,4.800E-001,4.900E-001,5.000E-001,5.100E-001,  &
                                                    5.200E-001,5.300E-001,5.400E-001,5.500E-001,5.600E-001,  &
                                                    5.700E-001,5.800E-001,5.900E-001,6.000E-001,6.100E-001,  &
                                                    6.200E-001,6.300E-001,6.400E-001,6.500E-001,6.600E-001,  &
                                                    6.700E-001,6.800E-001,6.900E-001,7.000E-001,7.100E-001,  &
                                                    7.200E-001,7.300E-001,7.400E-001,7.500E-001,7.600E-001,  &
                                                    7.700E-001,7.800E-001,7.900E-001,8.000E-001,8.100E-001,  &
                                                    8.200E-001,8.300E-001,8.400E-001,8.500E-001,8.600E-001,  &
                                                    8.700E-001,8.800E-001,8.900E-001,9.000E-001,9.100E-001,  &
                                                    9.200E-001,9.300E-001,9.400E-001,9.500E-001,9.600E-001,  &
                                                    9.700E-001,9.800E-001,9.900E-001,1.000E+000,1.010E+000,  &
                                                    1.020E+000,1.030E+000,1.040E+000,1.050E+000,1.060E+000,  &
                                                    1.070E+000,1.080E+000,1.090E+000,1.100E+000,1.110E+000,  &
                                                    1.120E+000,1.130E+000,1.140E+000,1.150E+000,1.160E+000,  &
                                                    1.170E+000,1.180E+000,1.190E+000,1.200E+000,1.210E+000,  &
                                                    1.220E+000,1.230E+000,1.240E+000,1.250E+000,1.260E+000,  &
                                                    1.270E+000,1.280E+000,1.290E+000,1.300E+000,1.310E+000,  &
                                                    1.320E+000,1.330E+000,1.340E+000,1.350E+000,1.360E+000,  &
                                                    1.370E+000,1.380E+000,1.390E+000,1.400E+000,1.410E+000,  &
                                                    1.420E+000,1.430E+000,1.440E+000,1.449E+000,1.460E+000,  &
                                                    1.471E+000,1.481E+000,1.493E+000,1.504E+000,1.515E+000,  &
                                                    1.527E+000,1.538E+000,1.563E+000,1.587E+000,1.613E+000,  &
                                                    1.650E+000,1.680E+000,1.700E+000,1.730E+000,1.760E+000,  &
                                                    1.800E+000,1.830E+000,1.840E+000,1.850E+000,1.855E+000,  &
                                                    1.860E+000,1.870E+000,1.890E+000,1.905E+000,1.923E+000,  &
                                                    1.942E+000,1.961E+000,1.980E+000,2.000E+00/)
      REAL(8), PARAMETER, DIMENSION(NWV) :: CHIS = (/9.565E-011,3.249E-011,2.000E-011,2.000E-011,2.000E-011,  &
                                                     2.000E-011,2.000E-011,2.000E-011,2.365E-011,2.669E-011,  &
                                                     3.135E-011,4.140E-011,6.268E-011,9.239E-011,1.325E-010,  &
                                                     1.956E-010,2.861E-010,4.172E-010,5.889E-010,8.036E-010,  &
                                                     1.076E-009,1.409E-009,1.813E-009,2.289E-009,2.839E-009,  &
                                                     3.461E-009,4.159E-009,4.930E-009,5.730E-009,6.890E-009,  &
                                                     8.580E-009,1.040E-008,1.220E-008,1.430E-008,1.660E-008,  &
                                                     1.890E-008,2.090E-008,2.400E-008,2.900E-008,3.440E-008,  &
                                                     4.030E-008,4.300E-008,4.920E-008,5.870E-008,7.080E-008,  &
                                                     8.580E-008,1.020E-007,1.180E-007,1.340E-007,1.400E-007,  &
                                                     1.430E-007,1.450E-007,1.510E-007,1.830E-007,2.150E-007,  &
                                                     2.650E-007,3.350E-007,3.920E-007,4.200E-007,4.440E-007,  &
                                                     4.740E-007,5.110E-007,5.530E-007,6.020E-007,7.550E-007,  &
                                                     9.260E-007,1.120E-006,1.330E-006,1.620E-006,2.000E-006,  &
                                                     2.250E-006,2.330E-006,2.330E-006,2.170E-006,1.960E-006,  &
                                                     1.810E-006,1.740E-006,1.730E-006,1.700E-006,1.760E-006,  &
                                                     1.820E-006,2.040E-006,2.250E-006,2.290E-006,3.040E-006,  &
                                                     3.840E-006,4.770E-006,5.760E-006,6.710E-006,8.660E-006,  &
                                                     1.020E-005,1.130E-005,1.220E-005,1.290E-005,1.320E-005,  &
                                                     1.350E-005,1.330E-005,1.320E-005,1.320E-005,1.310E-005,  &
                                                     1.320E-005,1.320E-005,1.340E-005,1.390E-005,1.420E-005,  &
                                                     1.480E-005,1.580E-005,1.740E-005,1.980E-005,3.442E-005,  &
                                                     5.959E-005,1.028E-004,1.516E-004,2.030E-004,2.942E-004,  &
                                                     3.987E-004,4.941E-004,5.532E-004,5.373E-004,5.143E-004,  &
                                                     4.908E-004,4.594E-004,3.858E-004,3.105E-004,2.659E-004,  &
                                                     2.361E-004,2.046E-004,1.875E-004,1.650E-004,1.522E-004,  &
                                                     1.411E-004,1.302E-004,1.310E-004,1.339E-004,1.377E-004,  &
                                                     1.432E-004,1.632E-004,2.566E-004,4.081E-004,7.060E-004,  &
                                                     1.108E-003,1.442E-003,1.614E-003,1.640E-003/)
!
!XH   interpolate log(CHI) linearly in log(WV)
      IF (WV .LT. 0.2D0) THEN
          STOP 'BRDF_SNOW in Mod_BRM.f90: wavelength less than the lower limit 0.2 micron'
      ELSE IF (WV .GT. 2.0D0) THEN
          STOP 'BRDF_SNOW in Mod_BRM.f90: wavelength greater than the upper limit 2.0 micron'
      ELSE
          DO IWV = 2, NWV
              IF (WVS(IWV-1) .LE. WV .AND. WV .LE. WVS(IWV)) THEN
!XH               LOG(CHI) = (LOG(CHIS(IWV))-LOG(CHIS(IWV-1)))/(LOG(WVS(IWV))-LOG(WVS(IWV-1)))*(LOG(WV)-LOG(WVS(IWV)))+LOG(CHIS(IWV))
                  CHI = CHIS(IWV)*DEXP(DLOG(CHIS(IWV)/CHIS(IWV-1))/DLOG(WVS(IWV)/WVS(IWV-1))*DLOG(WV/WVS(IWV)))
                  EXIT
              END IF
          END DO
      END IF
!
      ALPHA = SQRT(13.D0*D*4.D0*dpi*(CHI+M)/WV)
!XH   scattering angle THETA in degree
      THETA = d180_pi*DACOS(-MUS*MUV-SQRT((1.D0-MUS*MUS)*(1.D0-MUV*MUV))*DCOS(PHI*dpi_180))
!XH   approximate phase function of snow grain PTHETA
      PTHETA = 11.1D0*DEXP(-0.087D0*THETA)+1.1D0*DEXP(-0.014D0*THETA)
!XH   escape function K at MUS and MUV
      KMUS = 3.D0*(1.D0+2.D0*MUS)/7.D0
      KMUV = 3.D0*(1.D0+2.D0*MUV)/7.D0
!XH   reflectance without absorption R0
      R0 = 0.25D0*(1.247D0+1.186D0*(MUS+MUV)+5.157D0*MUS*MUV+PTHETA)/(MUS+MUV)
      BRDF_SNOW = R0*DEXP(-ALPHA*KMUS*KMUV/R0)
!
      RETURN
END FUNCTION BRDF_SNOW

COMPLEX(8) PURE FUNCTION REFR_OCEAN(WAVE,SALINITY,T_OCEAN)
!XH   Empirical equation for the refractive index of ocean (Applied Optics 34, 18, 1995)
!XH   SALINITY : salinity of seawater (g/kg)
!XH   T_OCEAN  : temperature of ocean surface (Celsius)
      REAL(8),    INTENT(IN) :: SALINITY
      REAL(8),    INTENT(IN) :: T_OCEAN
      REAL(8),    INTENT(IN) :: WAVE
!XH   CF0 - CF9: coefficients in empirical equation for the refractive index of ocean
      REAL(8), PARAMETER :: CF0 = 1.31405D0
      REAL(8), PARAMETER :: CF1 = 1.779D-4
      REAL(8), PARAMETER :: CF2 =-1.05D-6
      REAL(8), PARAMETER :: CF3 = 1.6D-8
      REAL(8), PARAMETER :: CF4 =-2.02D-6
      REAL(8), PARAMETER :: CF5 =15.868D0
      REAL(8), PARAMETER :: CF6 = 1.155D-2
      REAL(8), PARAMETER :: CF7 =-4.23D-3
      REAL(8), PARAMETER :: CF8 =-4.382D3
      REAL(8), PARAMETER :: CF9 = 1.1455D6
!
      REAL(8) :: WV, T, S
!
      IF (WAVE .GT. 0.7D0) THEN
          WV = 700.D0
      ELSE IF (WAVE .LT. 0.4D0) THEN
          WV = 400.D0
      ELSE
          WV = WAVE*1000.D0
      END IF

      IF (T_OCEAN .GT. 30.D0) THEN
          T = 30.D0
      ELSE IF (T_OCEAN .LT. 0.D0) THEN
          T = 0.D0
      ELSE
          T = T_OCEAN
      END IF

      IF (SALINITY .GT. 35.D0) THEN
          S = 35.D0
      ELSE IF (SALINITY .LT. 0.D0) THEN
          S = 0.D0
      ELSE
          S = SALINITY
      END IF

      REFR_OCEAN = CMPLX(CF0+(CF1+CF2*T+CF3*T*T)*S+CF4*T*T+(CF5+CF6*S+CF7*T)/WV+CF8/WV/WV+CF9/WV/WV/WV,0.D0)
!
      RETURN
END FUNCTION REFR_OCEAN


REAL(8) PURE FUNCTION FOAM_ALBEDO(WV)
!XH   Spectral dependent foam albedos are taken from R. Frouin Journal of Geophysical Research 101, C6, 1996
      REAL(8), INTENT(IN) :: WV
!
      INTEGER, PARAMETER :: NWV = 5
      REAL(8), PARAMETER, DIMENSION(NWV) :: WVS = (/0.44D0,0.85D0,1.02D0,1.65D0,2.25D0/)
      REAL(8), PARAMETER, DIMENSION(NWV) :: ALB = (/1.00D0,0.60D0,0.50D0,0.15D0,0.00D0/)
!
      INTEGER :: I
!
      IF (WV .LE. WVS(1)) THEN
          FOAM_ALBEDO = ALB(1)
      ELSE IF (WV .GE. WVS(NWV)) THEN
          FOAM_ALBEDO = ALB(NWV)
      ELSE
          DO I = 2, NWV
              IF (WV .LT. WVS(I)) THEN
                  FOAM_ALBEDO = ALB(I)+(ALB(I)-ALB(I-1))/(WVS(I)-WVS(I-1))*(WV-WVS(I))
                  EXIT
              END IF
          END DO
      END IF
!
      RETURN
END FUNCTION FOAM_ALBEDO

Subroutine BRM_ocean(n_el,i_BRDF,n_par,mu0,muv,phi,WAVE,x_vect,Rip)
!     i_BRDF,i_BPDF: flags for model
!     n_par        : number of model parameters
!     mu0,muv      : cosinus of solar an viewing zenith angles
!     phi          :  difference of azimuth angles of viewing and solar directions
!     WAVE         : wavelength in micron
!     x_vect       : vector of input parameters
!     Rip          : reflection matrix
      Implicit none
      Integer,    intent(in) :: n_el,i_BRDF
      Integer,    intent(in) :: n_par
      Real(8),    intent(in) :: mu0,muv,phi
      Real(8),    intent(in) :: x_vect(1:n_par)
      real(8),    intent(in) :: WAVE
      Real(8),    intent(out):: Rip(1:n_el,1:n_el)
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
      Real(8)    :: tet0,tetv,phi0,phiv
      complex(8) :: m
      real(8)    :: afoam

      if (IWAVE) then
!XH       Assume foam reflectance at 440nm is 22%
          afoam = FOAM_ALBEDO(WAVE)*0.22D0
!XH       Assume ocean salinity is 0.035 and temperature is 15 degree
          m     = REFR_OCEAN(WAVE,35.D0,15.D0)
      else
          afoam = 0.d0
          m     = cmplx(1.33d0,0.d0)
      end if

      phi0=0.d0
      phiv=phi*dpi_180
      
	  select case(i_BRDF)
      case(0)
!        simplest isotropic case
         call Refl_mat_ocean(n_el,.true.,n_par,mu0,muv,phi0,phiv,x_vect(1:n_par),m,afoam,Rip)
      case(1)
!        anisotropic case
         call Refl_mat_ocean(n_el,.false.,n_par,mu0,muv,phi0,phiv,x_vect(1:n_par),m,afoam,Rip)
	  case(9)
         call Refl_mat_3p_ocean(n_el,n_par,dacos(mu0),dacos(muv),  &
	                          phi0,phiv,x_vect(1:n_par),m,Rip)						  
      case default
         write(*,*) 'i_BRDF=',i_BRDF,' - unknown value (iBRM_water)'
         stop 'stop in Mod_BRM: BRM_ocean'
      
	  end select ! i_BRDF
!
End Subroutine BRM_ocean


!PL     i_BRDF,i_BPDF: flags for model
!PL     n_par is number of model parameters
!PL     n_par_i is number of parameters just for BRDF
!PL     mu0,muv: cosinus of solar an viewing zenith angles
!PL     phi:  difference of azimuth angles of viewing and solar directions
!PL     WAVE: wavelength in micron
!PL     x_vect: vector of input parameters
!PL     Rip: reflection matrix
Subroutine BRM(n_el,i_BRDF,i_BPDF,n_par,n_par_i,mu0,muv,phi,WAVE,x_vect,Rip)
      Implicit none
      Integer,  intent (in) :: n_el,i_BRDF,i_BPDF
      Integer,  intent (in) :: n_par,n_par_i
      Real(8),  intent (in) :: mu0,muv,phi
      Real(8),  intent (in) :: WAVE
      Real(8),  intent (in) :: x_vect(1:n_par)
      Real(8),  intent(out) :: Rip(1:n_el,1:n_el)
  
      Real(8) :: rho,vk,gteta,rhoH,Cv,R11
!      Real(8)    ::  x_vect_i(1:n_par)
      Real(8) :: tet0,tetv,phi0,phiv
      Complex(8) :: m
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

!XH    refractive index of land surface used in BPDF (may be changed below in case of snow surface)
       m = cmplx(1.5d0,0.d0)
!      Write(*,*) 'in  ', i_BRDF,i_BPDF,n_par,mu0,muv,phi,m
!      Write(*,*) x_vect
!      stop

!XH   test of snow surface
!      select case(3)

      select case(i_BRDF)
      case(0)
!        RPV model with HotSpot
         rho   = x_vect(1)
         vk    = x_vect(2)
         gteta = x_vect(3)
         if (n_par_i .le. 3) then
            rhoH  = rho
!           for Atakama test
!            rhoH=0.2775
         else
	        rhoH  = x_vect(4)
!			rhoH=0.2775
         end if
         R11=BRDF_RPV(mu0,muv,phi,rho,vk,gteta,rhoH)
      case(1)
!        Ross_Li Sparse model with HotSpot
         call LiSp_Ross(n_par_i,mu0,0.d0,muv,phi,x_vect(1:n_par_i),R11)
!PL       Write(*,'(4f18.10)') x_vect(1), x_vect(2), x_vect(3), R11
      case(2)
         tet0=dacos(mu0)
         tetv=dacos(muv)
         phi0=0.d0
         phiv=phi*dpi_180
!         call Refl_LitvBRM_fast(n_el,n_par,tet0,tetv,phi0,phiv,x_vect(1:n_par),m,Rip)
         call Refl_LitvBRM_R11(n_par_i,tet0,tetv,phi0,phiv,x_vect(1:n_par_i),R11)
      case(3)
!XH      snow surface
         R11=BRDF_SNOW(mu0,muv,phi,x_vect(2),x_vect(3),WAVE)
!         R11=x_vect(1)*BRDF_SNOW(mu0,muv,phi,x_vect(2),x_vect(3),WAVE)
!         R11=x_vect(1)+(1.d0-x_vect(1))*BRDF_SNOW(mu0,muv,phi,x_vect(2),x_vect(3),WAVE)
!XH      refractive index of snow surface used in BPDF
         m = cmplx(1.31d0,0.d0)
	  case(9)
         tet0=dacos(mu0)
         tetv=dacos(muv)
         phi0=0.d0
         phiv=phi*dpi_180
	     call Refl_mat_4p(n_el,n_par,tet0,tetv,phi0,phiv,x_vect(1:n_par),m,Rip)
!        call Refl_LitvBRM_fast(n_el,n_par,tet0,tetv,phi0,phiv,x_vect(1:n_par),m,Rip)
!		 call Refl_mat_4p_simp(n_el,n_par,tet0,tetv,phi0,phiv,x_vect(1:n_par),m,Rip)
      case default
         write(*,*) 'i_BRDF=',i_BRDF,' - unknown value'
         stop 'stop in Mod_BRM'
      end select ! i_BRDF

      select case(i_BPDF)
      case(0)
!        one parametric BPDF model of Maignan and Breon
         Cv=x_vect(n_par_i+1)
!XH      modified Maign and Breon model based on Fresnel scattering matrix R_Fren_My
         call Maign_Breon_XH(n_el,mu0,muv,phi,m,Cv,Rip)
!        BRDF + Maignan_Fr_matr
         Rip(1,1)= Rip(1,1)+R11
!XH      to prevent unphysical HDR when applying BRDF of snow surface
!         Rip(1,1)= R11
      case(1)
!        BRDF + My BPDF model developed in SRON
         x_vect_p(1) = 1.5d0
         x_vect_p(2) = x_vect(n_par_i+1)
         x_vect_p(3) = x_vect(n_par_i+2)
         x_vect_p(4) = 1.5d0
         phiv=phi*dpi_180
         call R_Fren_PL_SRON(n_el,n_par_p,mu0,dpi,muv,phiv,x_vect_p,Rip)
!         Rip(1,1)= R11
         Rip(1,1)= Rip(1,1)+R11
      case(-1)
!        BRDF only
         Rip=0.d0
         Rip(1,1)= R11
      end select ! i_BPDF
!
End Subroutine

!PL RPV model with coe. back. effect "CB_peak_new(tet)"	  
Real(8) FUNCTION BRDF_RPV(cs,cv,phi,rho,vk,gt,rhoH)
!     RPV model
!XH   observation geometry is based on the POLDER convention
      Implicit none
      Real(8), intent(in):: cs,cv,phi,rho,vk,gt,rhoH
      Real(8) :: xx,yy,zz,FF1,ww,aa,FF2,vv,G,FF3,tet
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
!
      xx=dabs(cs)**(vk-1.d0)
      yy=dabs(cv)**(vk-1.d0)
      zz=(dabs(cs)+dabs(cv))**(1.d0-vk)
      FF1=rho*xx*yy/zz

      xx=dsqrt(1.d0-cs*cs)
      yy=dsqrt(1.d0-cv*cv)
      ww=cs*cv+xx*yy*dcos(phi*dpi_180)
      aa=1.d0+gt*gt+2.d0*gt*ww
      FF2=(1.d0-gt*gt)/(aa**1.5d0)

      vv=xx/cs
      ww=yy/cv
      G=dsqrt(vv*vv+ww*ww-2.d0*vv*ww*dcos(phi*dpi_180))
!      FF3=1+(1-rho)/(1+G)
      FF3=1.d0+(1.d0-rhoH)/(1.d0+G)

!	  tet=-dacos(ww)
!	  FF3=1.d0+(CB_peak_new(tet)-1.d0)/(1.d0+G) 
!	  FF3=CB_peak_new(tet)+(1.d0-rhoH)/(1.d0+G)

      BRDF_RPV=FF1*FF2*FF3
!
	  return
end FUNCTION BRDF_RPV


!PL     Li Sparse- Ross model with azimuth function
!PL     x_vect_i(3): geometric kernel
!PL     x_vect_i(2): volumatric kernel
!PL     x_vect_i(1): isotropic kernel
Subroutine LiSp_Ross(n_par_i,mu0,phi0,mu,Phir,x_vect_i,LiD_R)
      Implicit none
!     Parameter dl_par1 is diferent for soil and vegetation surfaces
!     soil dl_par1=1.
!     veg  dl_par1=2.
      Real(8), parameter :: dl_par1=2.d0, alpha0=1.5d0
      Integer, intent(in):: n_par_i
      Real(8), intent(in):: phi0,Phir,x_vect_i(1:n_par_i)
      Real(8),intent(out):: LiD_R
      Real(8) :: cos_scat,sin_scat,tet_sc,Ross_K,LiD_K,mu0,mu
      Integer :: i
      Real(8) :: pr_t1,tan_i,tan_v,phi,D,sec_i_v,cos_t,sin_t,t,O,muv_lim,mu0_lim
      Real(8) :: ni_x,ni_y,ni_z, nv_x,nv_y,nv_z, phi_inc, phi_v, sin_v,sin_inc,del_phi1
      Real(8) :: Hot_Sp,alpha
      Real(8) :: Lid_R_limit
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

      tan_i=dsqrt(1.d0-mu0*mu0)/mu0
      tan_v=dsqrt(1.d0-mu*mu)/mu

      phi_v=Phir*dpi_180
      phi_inc=phi0*dpi_180

      sin_inc=dsqrt(1.d0-mu0*mu0)
      ni_x=sin_inc*dcos(phi_inc)
      ni_y=sin_inc*dsin(phi_inc)
      ni_z=mu0 !(cos of solar zenith angle)

      sin_v=dsqrt(1.d0-mu*mu)
      nv_x=sin_v*dcos(phi_v)
      nv_y=sin_v*dsin(phi_v)
      nv_z=mu !(cos of viewing zenith angle)

      del_phi1=phi_v-phi_inc
      cos_scat=-(ni_x*nv_x+ni_y*nv_y+ni_z*nv_z)
!XH   for potential computational error in cos_scat
      if (cos_scat .lt. -1.d0) then
         cos_scat=-1.d0
         sin_scat=0.d0
         tet_sc=dpi
      else if (cos_scat .gt. 1.d0) then
         cos_scat=1.d0
         sin_scat=0.d0
         tet_sc=0.d0
      else
         sin_scat=dsqrt(1.d0-cos_scat*cos_scat)
         tet_sc=dacos(cos_scat)
      end if

      alpha=(dpi-tet_sc)*d180_pi
      Hot_Sp=(1.d0+1.d0/(1.d0+alpha/alpha0))

      Ross_K=Hot_Sp*((dpi*0.5d0-tet_sc)*cos_scat+sin_scat)/(mu+mu0)-dpi*0.25d0
  
!     here positive angle corresponds to backscattering direction
      sec_i_v=(1.d0/mu+1.d0/mu0)
      D=dsqrt(tan_i*tan_i+tan_v*tan_v-2.d0*tan_i*tan_v*dcos(del_phi1))
  
      cos_t=dsqrt(D*D+(tan_i*tan_v*dsin(del_phi1))**2)/sec_i_v*dl_par1
      If (dabs(cos_t) .gt. 1.d0) cos_t=1.d0

      sin_t=dsqrt(1.d0-cos_t*cos_t)
      t=dacos(cos_t)
      O=1.d0/dpi*(t-sin_t*cos_t)*sec_i_v
 
      LiD_K=O-sec_i_v+(1.d0-cos_scat)/mu/mu0*0.5d0

      Lid_R_limit=(1.d0-x_vect_i(1))*0.5d0
!      Lid_R_limit=(1.d0-x_vect_i(1))/4.d0
      If (LiD_K*x_vect_i(3)*x_vect_i(1)*mu .gt. Lid_R_limit) LiD_K=Lid_R_limit/x_vect_i(3)/x_vect_i(1)/mu
!PL      If (LiD_K*x_vect_i(3)*mu .gt. Lid_R_limit) LiD_K=Lid_R_limit/x_vect_i(3)/mu

      LiD_R=(LiD_K*x_vect_i(3)+Ross_K*x_vect_i(2)+1.d0)*x_vect_i(1)
!PL      LiD_R=LiD_K*x_vect_i(3)+(Ross_K*x_vect_i(2)+1.d0)*x_vect_i(1)
!      LiD_R=LiD_K*x_vect_i(3)+Ross_K*x_vect_i(2)+x_vect_i(1)

      if (LiD_R .lt. 0.d0) LiD_R=0.d0
!
End Subroutine LiSp_Ross


!PL     one parametric BPDF model of Maignan and Breon
Subroutine Maign_Breon_XH(n_el,mu0,muv,phi,m,Cmgn,R_fr)
!XH   modified one parametric BPDF model of Maignan and Breon based on Fresnel reflection matrix
      Implicit none
      integer,    intent(in) :: n_el
      Complex(8), intent(in) :: m
      Real(8),    intent(in) :: mu0,muv,phi,Cmgn
      Real(8),    intent(out):: R_fr(1:n_el,1:n_el)
!
      Real(8) :: koef,alf
!
!XH   alf is the phase angle which is equal to PI minus scattering angle
      call R_Fren_My(n_el,mu0,muv,phi,m,R_fr,alf)
      koef=Cmgn*dexp(-dtan(alf*0.5d0))*0.25d0/(mu0+muv)
      R_fr=koef*R_fr
!
      return
end subroutine Maign_Breon_XH


!PL one parametric BPDF model of Maignan and Breon
Subroutine Maign_Breon(mu01,muv1,phi,m,Cmgn,R_fr)
  Implicit none
  Complex(8), intent(in):: m
  Real(8), intent(in):: mu01,muv1,phi,Cmgn
  Real(8), intent(out):: R_fr(1:4,1:4)
  
  Real(8), parameter:: st_min=0.001,sin_eps=0.00000001
  Real(8) zz,alf,uu,koef,mu_alf !,pi
  Complex(8) r1,r2,r_par1,r_par2,r_perp1,r_perp2,r_par,r_perp,&
             r_par_per1,r_par_per2
  Real(8) r_par_sq,r_perp_sq
  Real(8) del_phi,sin_ti,sin_tv,sin_tet,cos_tet,cos_etv,cos_eti, &
          sin_eti,sin_etv,cos_2etv,sin_2etv,cos_2eti,sin_2eti
  Real(8) L_v(1:4,1:4), L_inc(1:4,1:4), F_fr(1:4,1:4),mu0,muv,tet0v,sin_phi

!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

  mu0=mu01
  muv=muv1
  !pi=4.d0*datan(1.d0)

!  sin_eps=dsin(st_min*dpi_180)

  sin_tv=dsqrt(1.d0-muv*muv)
  sin_ti=dsqrt(1.d0-mu0*mu0)
  zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
  cos_tet=-zz
  sin_tet=dsqrt(1.d0-cos_tet*cos_tet)


  If (dabs(sin_tv) .lt. sin_eps) then
   tet0v=dacos(muv)	  
   tet0v=dabs(tet0v-st_min*dpi_180)
   muv=dcos(tet0v)
   sin_tv=dsqrt(1.d0-muv*muv)
   zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
   cos_tet=-zz
   sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
  endif
 
  If (dabs(sin_ti) .lt. sin_eps) then
   tet0v=dacos(mu0)	  
   tet0v=dabs(tet0v-st_min*dpi_180)
   mu0=dcos(tet0v)
   sin_ti=dsqrt(1.d0-mu0*mu0)
   zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
   cos_tet=-zz
   sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
  endif

  if (dabs(sin_tet) .lt. sin_eps) then
   tet0v=dacos(mu0)	  
   tet0v=dabs(tet0v-st_min*dpi_180)
   muv=dcos(tet0v)
   sin_tv=dsqrt(1.d0-muv*muv)
   zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
   cos_tet=-zz
   sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
  endif

  if (dabs(sin_tet) .lt. sin_eps) write(*,*) 'sin_tet',sin_tet
  if (dabs(sin_tv)  .lt. sin_eps) write(*,*) 'sin_tv',sin_tv
  if (dabs(sin_ti)  .lt. sin_eps) write(*,*) 'sin_ti',sin_ti

  alf=dacos(zz)
  uu=dtan(alf*0.5d0)
  mu_alf=dcos(alf*0.5d0)
  koef=Cmgn*dexp(-uu)*0.25d0/(mu0+muv)

  r1=m*m*mu_alf
  r2=cdsqrt(m*m-1.d0+mu_alf*mu_alf)
  r_par1=r1-r2
  r_par2=r1+r2

  r_perp1=mu_alf-r2
  r_perp2=mu_alf+r2

  r_par=r_par1/r_par2
  r_perp=r_perp1/r_perp2

  r_par_sq=dble(r_par*dconjg(r_par))
  r_perp_sq=dble(r_perp*dconjg(r_perp))
  r_par_per1=r_par*dconjg(r_perp)
  r_par_per2=dconjg(r_par_per1)

  F_fr=0.d0
  F_fr(1,1)=(r_par_sq+r_perp_sq)*koef*0.5d0
  F_fr(2,2)= F_fr(1,1)
  F_fr(2,1)=(r_par_sq-r_perp_sq)*koef*0.5d0
  F_fr(1,2)= F_fr(2,1)
  F_fr(3,3)=dble(r_par_per1+r_par_per2)*koef*0.5d0
  F_fr(4,4)=F_fr(3,3)
  F_fr(3,4)=dble((0,-1.d0)*(r_par_per1-r_par_per2))*koef*0.5d0
  F_fr(4,3)=-F_fr(3,4)

  del_phi=phi*dpi_180-dpi
  cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
  cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

  sin_etv=dsin(del_phi)*sin_ti/sin_tet
  sin_eti=dsin(del_phi)*sin_tv/sin_tet

  cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
  sin_2etv=2.d0*cos_etv*sin_etv

  cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
  sin_2eti=2.d0*cos_eti*sin_eti

  L_v=0.d0
  L_inc=0.d0

  sin_phi=dsin(phi*dpi_180)
  If (dabs(sin_phi) .lt. sin_eps/10.d0) then
   L_v(1,1)=1.d0
   L_v(2,2)=1.d0
   L_v(2,3)=0
   L_v(3,3)=1.d0
   L_v(3,2)=0
   L_v(4,4)=1.d0

   L_inc(1,1)=1.d0
   L_inc(2,2)=1.d0
   L_inc(2,3)=0
   L_inc(3,3)=1.d0
   L_inc(3,2)=0
   L_inc(4,4)=1.d0
  else
   L_v(1,1)=1.d0
   L_v(2,2)=cos_2etv
   L_v(2,3)=sin_2etv
   L_v(3,3)=cos_2etv
   L_v(3,2)=-sin_2etv
   L_v(4,4)=1.d0

   L_inc(1,1)=1.d0
   L_inc(2,2)=cos_2eti
   L_inc(2,3)=sin_2eti
   L_inc(3,3)=cos_2eti
   L_inc(3,2)=-sin_2eti
   L_inc(4,4)=1.d0

!!!!!!PL   
!    L_v(2,3)=- L_v(2,3)
!	L_v(3,2)=- L_v(3,2)

!	L_inc(2,3)=- L_inc(2,3)
!	L_inc(3,2)=- L_inc(3,2)

  endif

  R_fr=Matmul(L_v,Matmul(F_fr,L_inc))

end Subroutine Maign_Breon


!PL one parametric BPDF model of Maignan and Breon
Subroutine Fresnel_prov(mu0,muv,phi,xind,Cmgn,R_fr)
  Implicit none
  Real(8), intent(in):: xind
  Real(8), intent(in):: mu0,muv,phi,Cmgn
  Real(8), intent(out):: R_fr(1:4,1:4)

  Real(8) yy,zz,xx,rl,rr,mu_alf,uu,alf,koef !,pi
!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
  
  !pi=4.d0*datan(1.d0)
  zz=mu0*muv+dsqrt(1.d0-mu0*mu0)*dsqrt(1.d0-muv*muv)*dcos(phi*dpi_180)
  IF(zz .GT. 1.0) zz = 1.0
  IF(zz .LT.-1.0) zz =-1.0
  alf=dacos(zz)
  uu=dtan(alf*0.5d0)
  mu_alf=dcos(alf*0.5d0)
  koef=Cmgn*dexp(-uu)*0.25d0/(mu0+muv)

  xx=mu_alf
  yy=dsqrt(xind*xind-1.0+xx*xx)
  zz=xx*xind*xind
  rl=(zz-yy)/(zz+yy)
  rr=(xx-yy)/(xx+yy)
  R_fr(1,1)=0.5*(rl*rl+rr*rr)*koef
  R_fr(1,2)=0.5*(rl*rl-rr*rr)*koef
  R_fr(2,1)=R_fr(1,2)
  R_fr(2,2)=R_fr(1,1)
  R_fr(3,3)=rl*rr*koef
 
 end Subroutine


!PL     Fourier expansion
!PL     For land  surfaces: land_percent=100 !i_ocean_land=2
!PL     For water  surfaces: land_percent=0
!PL     For mixed water/land  surfaces: 0 < land_percent < 100
!PL     M_max is maximum number of expansion
!PL     n_par is number of model parameters
!PL     n_par_i is number of parameters just for BRDF
!PL     x_vect: vector of input parameters
!PL     WAVE: wavelength in micron
!PL     R_exp: matrix of Fourier expansion coeffcients
      Subroutine BRM_Fexp_OSH(n_el,land_percent,M_max,n_tet,mu,WAVE,     &
	                          i_BRDF,i_BPDF,n_par,n_par_i,x_vect,        &
							  i_BRM,n_par_water,x_vect_water,            &
							  R_exp)
!XH   surface properties may be wavelength dependent
      Implicit none
      integer,    intent(in) :: n_el
      Integer,    intent(in) :: n_par,n_par_i,M_max, n_tet, i_BRDF,i_BPDF
	  Real,       intent(in) :: land_percent
      Real(8),    intent(in) :: WAVE
      Real(8),    intent(in) :: mu(1:n_tet),x_vect(1:n_par)
	  Integer,    intent(in) :: n_par_water,i_BRM
      Real(8),    intent(in) :: x_vect_water(1:n_par_water)
      Real(8),    intent(out):: R_exp(1:n_el,1:n_el,1:n_tet,1:n_tet,0:M_max)
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
      Integer :: i1,i2,i3,j1,j2,i_m,n_phi
      Real(8) :: low_lp, up_lp, Rm(1:n_el,1:n_el),Rm_ocean(1:n_el,1:n_el)
      Real(8) :: phi(1:nph*M_max),Gw_p(1:nph*M_max)
!      Real(8) :: mu_mphi(0:M_max,1:nph*M_max),sin_mphi(0:M_max,1:nph*M_max), Int_phi, sin_cos
      Real(8) :: mu_mphi,sin_mphi !(1:nph*M_max),sin_mphi(1:nph*M_max)
!PL      Real(8) :: Rmat(1:n_tet,1:n_tet,1:nph*M_max,1:n_el,1:n_el)

      n_phi=nph*M_max
 
      If (M_max .gt. M_max_max) then
         Write(*,*) 'There are problems with number of ', &
                    'M_max_max (Mod_BRM.f90). Increase M_max_max=',M_max
	     stop
      end if

      If (n_tet .gt. NN0_my) then
	     Write(*,*) 'There are problems with number of', &
                    ' angles (Mod_BRM.f90). Increase NN0_my=',n_tet
	     stop
      end if

      If (n_phi .gt. n_phi_max) then
	     Write(*,*) 'There are problems with n_phi_max in Mod_BRM.f90.', &
	                'Increase phi_max ',n_phi
         stop
      end if

      low_lp=0.d0
      up_lp=dpi
!      up_lp=2.d0*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
!      If (i_num .eq. 0) then
!      If (i_fd .le. 1 .or. i_fd .eq. 2) then
      If (i_gauleg .eq. 1) then
!      If ((i_gauleg .eq. 1) .and. (i_fd.le.1)) then
         Do i3=1, n_phi
            Gw_p(i3)=Gw_p_old(i3)
            phi(i3) = phi_old(i3)
         end do
      else
         call gauss_l(n_phi,n_phi,low_lp,up_lp,phi,Gw_p)
         Do i3=1, n_phi
            Gw_p_old(i3)=Gw_p(i3)
	        phi_old(i3) = phi(i3)
         end do !i3
         i_gauleg=1
      end if !i_gauleg
!     end if !i_fd

      R_exp=0.d0

      Do i1=1, n_tet
!        here we use reciprocity principle
         Do i2=i1, n_tet
            Do i3=1, n_phi
	           select case(nint(land_percent))
               case(0)
!PL               The input geometry for BRM and BRM_ocean is similar to POLDER geometry, but Fourier expansion inside of subroutine
!PL               is reculclated for the geometry related to incident direction. Therefore it is
!PL               necessary to take '180.0-phi(i3)'
	              call BRM_ocean(n_el,i_BRM,n_par_water,mu(i1),mu(i2),(180.d0-phi(i3)*d180_pi), &
				                 WAVE,x_vect_water(1:n_par_water),Rm)
               case(100)
                  call BRM(n_el,i_BRDF,i_BPDF,n_par,n_par_i,mu(i1),mu(i2),(180.d0-phi(i3)*d180_pi),&
				                 WAVE,x_vect(1:n_par),Rm)
	           case default

                  If ((land_percent .lt. 0) .or. (land_percent .gt. 100)) then
					write(*,*) 'land_percent=',land_percent,' - unknown value'
					stop 'stop in BRM_Fexp_OSH. "land_percent" is not correct'
				  endif 

                  call BRM_ocean(n_el,i_BRM,n_par_water,mu(i1),mu(i2),(180.d0-phi(i3)*d180_pi), &
				                 WAVE,x_vect_water(1:n_par_water),Rm_ocean)

                  call BRM(n_el,i_BRDF,i_BPDF,n_par,n_par_i,mu(i1),mu(i2),(180.d0-phi(i3)*d180_pi),&
				                 WAVE,x_vect(1:n_par),Rm)

				  Rm=Rm_ocean*(1. - 0.01*land_percent)+0.01*land_percent*Rm				 

	           end select

			   do i_m=0,M_max
			     mu_mphi=Gw_p(i3)*dcos(i_m*phi(i3))/dpi
				 sin_mphi=Gw_p(i3)*dsin(i_m*phi(i3))/dpi
			     R_exp(1,1,i1,i2,i_m)=R_exp(1,1,i1,i2,i_m)+mu_mphi*Rm(1,1)
                 if (n_el .ge. 3) then
                  R_exp(2,1,i1,i2,i_m)=R_exp(2,1,i1,i2,i_m)+mu_mphi*Rm(2,1)
                  R_exp(1,2,i1,i2,i_m)=R_exp(1,2,i1,i2,i_m)+mu_mphi*Rm(1,2)
                  R_exp(2,2,i1,i2,i_m)=R_exp(2,2,i1,i2,i_m)+mu_mphi*Rm(2,2)
                  R_exp(3,3,i1,i2,i_m)=R_exp(3,3,i1,i2,i_m)+mu_mphi*Rm(3,3)
                  R_exp(3,1,i1,i2,i_m)=R_exp(3,1,i1,i2,i_m)+sin_mphi*Rm(3,1)
                  R_exp(3,2,i1,i2,i_m)=R_exp(3,2,i1,i2,i_m)+sin_mphi*Rm(3,2)
                  R_exp(1,3,i1,i2,i_m)=R_exp(1,3,i1,i2,i_m)+sin_mphi*Rm(1,3)
                  R_exp(2,3,i1,i2,i_m)=R_exp(2,3,i1,i2,i_m)+sin_mphi*Rm(2,3)
				 endif
			    enddo !i_m

            End do !i3
         End do !i2
      End do !i1

      do i_m=0,M_max
         do i1=1,n_tet
            do i2=i1,n_tet
!XH            applying reciprocity principle
               if (i1 .ne. i2) then
                  R_exp(1,1,i2,i1,i_m)=R_exp(1,1,i1,i2,i_m)
                  if (n_el .ge. 3) then
                     R_exp(2,1,i2,i1,i_m)=R_exp(1,2,i1,i2,i_m)
                     R_exp(1,2,i2,i1,i_m)=R_exp(2,1,i1,i2,i_m)
                     R_exp(2,2,i2,i1,i_m)=R_exp(2,2,i1,i2,i_m)
                     R_exp(3,3,i2,i1,i_m)=R_exp(3,3,i1,i2,i_m)
                     R_exp(3,1,i2,i1,i_m)=R_exp(1,3,i1,i2,i_m)
                     R_exp(3,2,i2,i1,i_m)=R_exp(2,3,i1,i2,i_m)
                     R_exp(1,3,i2,i1,i_m)=R_exp(3,1,i1,i2,i_m)
                     R_exp(2,3,i2,i1,i_m)=R_exp(3,2,i1,i2,i_m)
!                     if (n_el .eq. 4) then
!                        R_exp(4,3,i2,i1,i_m)=-R_exp(3,4,i1,i2,i_m)
!                        R_exp(3,4,i2,i1,i_m)=-R_exp(4,3,i1,i2,i_m)
!                        R_exp(4,4,i2,i1,i_m)=R_exp(4,4,i1,i2,i_m)
!                        R_exp(4,1,i2,i1,i_m)=-R_exp(1,4,i1,i2,i_m)
!                        R_exp(4,2,i2,i1,i_m)=-R_exp(2,4,i1,i2,i_m)
!                        R_exp(1,4,i2,i1,i_m)=-R_exp(4,1,i1,i2,i_m)
!                        R_exp(2,4,i2,i1,i_m)=-R_exp(4,2,i1,i2,i_m)
!                     end if
                  end if
               end if 
            end do !i2
         end do !i1
      end do !i_m

!     original code by Pavel
!      Do i_m=0,M_max
!         Do i3=1, n_phi
!            mu_mphi(i_m,i3) =dcos(i_m*phi(i3))
!            sin_mphi(i_m,i3)=dsin(i_m*phi(i3))
!         end do
!      end do
!
!      Do i1=1,n_tet
!         Do i2=i1,n_tet
!            Do j1=1, n_el
!               Do j2=1, n_el
!	              Do i_m=0,M_max
!                     Int_phi=0.d0
!                     Do i3=1, n_phi
!!                       here we use the mirror symmetry of the scattering object
!	                    If (((j1 .eq. 3) .and. (j2 .eq. 1)) .or. &
!	                        ((j1 .eq. 1) .and. (j2 .eq. 3)) .or. &
!		                    ((j1 .eq. 2) .and. (j2 .eq. 3)) .or. &
!		                    ((j1 .eq. 3) .and. (j2 .eq. 2)) .or. &
!		                    ((j1 .eq. 2) .and. (j2 .eq. 4)) .or. &
!		                    ((j1 .eq. 4) .and. (j2 .eq. 2))) then
!		                   sin_cos=sin_mphi(i_m,i3)
!	                    else
!		                   sin_cos= mu_mphi(i_m,i3)
!                        end if
!                        Int_phi=Int_phi+Gw_p(i3)*Rmat(i1,i2,i3,j1,j2)*sin_cos
!	                 end do !i3
!                     R_exp(j1,j2,i1,i2,i_m)=Int_phi/dpi
!!                    here we use reciprocity principle
!                     If (i1 .ne. i2) then
!                        R_exp(j2,j1,i2,i1,i_m)=Int_phi/dpi
!!PL                     recip. princ. for 4,j1 elements
!!PL                      If (j2 .eq. 4) R_exp(j2,j1,i2,i1,i_m)=-R_exp(j2,j1,i2,i1,i_m)
!		             end if
!	              end do!i_m
!	           end do!j2
!	        end do!j1
!         end do!i2
!      end do!i1
!
      End Subroutine BRM_Fexp_OSH



!PL  Coherent Backscattering for infinite medium
 Real(8) Function CB_peak_new(tet)
 Implicit none
 Real(8),parameter:: eps=0.00001,kl_ext=30.d0,Coh_peak=0.5
 Real(8), intent(in):: tet
 Real(8) q!,pi

  q=2.*kl_ext*dcos(tet*0.5d0)

  If (dabs(tet-dpi) .gt. eps)  then
   CB_peak_new=1.+Coh_peak*datan(q)/q
  else 
   CB_peak_new=1.+Coh_peak
  endif

 end function


!PL CB for infinite medium (test)
 Real(8) Function CB_peak_new_pr(tet,kl_ext)
 Implicit none
 Real(8),parameter:: eps=0.00001
 Real(8), intent(in):: tet
 Real(8) kl_ext,q,pi,CB_ampl
 
  CB_ampl=2.
!  CB_ampl=1./1.!0.8
 ! kl_ext=50

  pi=4.d0*datan(1.d0)
  q=2.*kl_ext*dcos(tet/2.)

  If (dabs(tet-pi) .gt. eps)  then
   CB_peak_new_pr=1.+datan(q)/q/CB_ampl
  else 
   CB_peak_new_pr=1.+1./CB_ampl
  endif

!  CB_peak_new=CB_peak_new*CB_ampl

end function



!PL reflection matrix for ocean based on
!PL isotropic and anisotropic Cox-Munk model
Subroutine Refl_mat_ocean(n_el,iso,n_par,mu0,muv,phi0,phiv,x_vect,m,A_foam,Rm)
      Implicit none
      Logical,    intent(in) :: iso
      Integer,    intent(in) :: n_el,n_par
      Real(8),    intent(in) :: mu0,phi0,muv,phiv,x_vect(1:n_par)
      Real(8),    intent(in) :: A_foam
      Complex(8), intent(in) :: m
      Real(8),    intent(out):: Rm(1:n_el,1:n_el)
!
      Real(8), parameter :: dpi_180=dpi/180.d0
      !Real(8), parameter :: A_foam=0.2d0
      Real(8) :: phii,mu01,muv1,phiv1
      Real(8) :: sig1,alb,Fren_fr,mu_n1
      Real(8) :: f_pdf1,del_phi,Shad,cos_tet,tet_sc
      Real(8) :: phiw, sigmau2, sigmac2
      Real(8) :: zx,zy,zu,zc,sin0,sinv,cphi0,cphiv,sphi0,sphiv

      alb    = x_vect(1)
      Fren_fr= x_vect(2)
      sig1   = x_vect(3)
      If (.not. iso) then
!         phiw    = x_vect(4)*dpi_180
!PL      Here phiw in written in radian !!!!!!!!!!
         phiw    = x_vect(4)
         sigmac2 = (2.d0*sig1-0.003d0)/0.00512d0*0.00192d0+0.003d0
         sigmau2 = (2.d0*sig1-0.003d0)/0.00512d0*0.00316d0
!         sigmac2 = sig1
!         sigmau2 = sig1
      end if

      del_phi = phiv-phi0
      sin0    = dsqrt(1.d0-mu0*mu0)
      sinv    = dsqrt(1.d0-muv*muv)
      cos_tet =-mu0*muv-sin0*sinv*dcos(del_phi)
      tet_sc  = dacos(cos_tet)
      mu_n1   = (mu0+muv)/(2.d0*dsin(0.5d0*tet_sc))

      If (iso) then
         f_pdf1=dexp(-(1.d0-mu_n1*mu_n1)/(2.d0*sig1*mu_n1*mu_n1))/(2.d0*dpi*sig1*mu_n1**3)

      else
         cphi0 =-dcos(phi0)
         cphiv = dcos(phiv)
         sphi0 =-dsin(phi0)
         sphiv = dsin(phiv)
         zx    = (sinv*cphiv-sin0*cphi0)/(mu0+muv)
         zy    = (sinv*sphiv-sin0*sphi0)/(mu0+muv)
         zu    = dcos(phiw)*zx+dsin(phiw)*zy
         zc    =-dsin(phiw)*zx+dcos(phiw)*zy
         f_pdf1= dexp(-0.5d0*(zc*zc/sigmac2+zu*zu/sigmau2))
         f_pdf1= f_pdf1/(2.d0*dpi*sqrt(sigmac2*sigmau2)*mu_n1**3)
      end if

      call R_Fren_My(n_el,mu0,muv,del_phi*d180_pi,m,Rm)

!      phii=phi0+dpi
!      mu01=mu0
!      muv1=muv
!      phiv1=phiv
!      call R_Fren(mu01,phii,muv1,phiv1,Rm)

      call Shad_func(mu0,muv,sig1,Shad)

      Rm=Fren_fr*Shad*f_pdf1*dpi/(4.d0*mu0*muv*mu_n1)*Rm
!      Rm=Fren_fr*f_pdf1*dpi/(4.d0*mu0*muv*mu_n1)*Rm
      Rm(1,1)=Fren_fr*alb+Rm(1,1)+(1.d0-Fren_fr)*A_foam
!      Rm(1,1)=alb*Shad+Rm(1,1)!+(1.d0-Fren_fr)*A_foam
!      Rm(1,1)=(Fren_fr*alb+Rm(1,1)+(1.d0-Fren_fr)*A_foam)*Shad
!
      return
End Subroutine Refl_mat_ocean




!PL Subroutine for Ocen surface based on 
!PL Pavel Litvinov, Otto Hasekamp, Brian Cairns. Remote Sensing of Environment (2011), 115, 781–792
!PL (Experimental version)
Subroutine  Refl_mat_3p_ocean_dep(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,Rm)
!XH   add input variable n_el
!     New physical model (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!     n_par - the number of model parameters
!     n_tet - the number of viewing zenith angle
!     tet0  - solar zenith angle
!     tet_v - array of viewing zenith angles
!     phi0  - azimuzth angles of solar direction
!     phiv  - azimuzth angles of viewing direction
!     x_vect- vector of input parameters
!     In given version the meaning of parameters is follows:
!        x_vect(1) is ocean abedo
!        x_vect(2) is rms of surface slope (in scale 1)
!        x_vect(3) is Fraction of surface wich can provide Fresnel reflection
!     Rm: calculated reflection matrix
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sin_eps=0.00000001d0

      Integer :: i1,i2
      Real(8) :: del_phi
      Real(8) :: Randf1(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,Int_bet1_Fr,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_tet,sin_ti,sin_tv, cos_etv,cos_eti, sin_etv,sin_eti
      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi,Depol_22

      Real(8) Depol

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

	  Real(8) mu_n0,f_pdf0
      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

	  Fren_Fr=1.

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

!      Do i1=1,n_par
!         boolean=(boolean .and. (dabs(x_vect_old(i1)-x_vect(i1)) .le. eps))
!      end do
!XH   take advantage of new feature in Fortran95
      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s     =x_vect(1)
         Depol     =x_vect(2)
         sig1      =x_vect(3)
		 k_shad    =sig1

		 mu0=dcos(tet0)
         muv=dcos(tetv)
 
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
         If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
         end if
 
         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

		 call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,dcmplx(1.33d0,0.d0),Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
        
         mu_n0 = (mu0+muv)/(2.d0*dsin(0.5d0*tet))
         f_pdf0=dexp(-(1.d0-mu_n0*mu_n0)/(2.d0*sig1*mu_n0*mu_n0))/(2.d0*dpi*sig1*mu_n0**3)

         Rm=Fren_fr*Shad*f_pdf0*Rm*dpi/(4.d0*mu0*muv*mu_n0)

!XH      new computation for Randf1 to reduce redundant calculations
         Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
                  Randf1(i1,i2)   = mu1v*mu2i*f_pdf1(i2)
               else
                  Randf1(i1,i2)   =0.d0
	           end if
            end do !i1
         end do !i2

!XH      new two dimensional integral taking advantage of the fact that n_alf>n_bet and n_alf is number of lines of Randf1
         Int_alfbet1   =0.d0
         Int_alfbet1_Fr=0.d0
         do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
         end do ! i2

		 Rm(1,1)=alb_s*Int_alfbet1/mu0/muv*Shad+Rm(1,1)

         if (n_el .ge. 3) then
!XH         for vector case only
            sin_tet=dsqrt(1.d0-cos_tet*cos_tet)

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

            P_Rel=-(1.d0-cos_tet*cos_tet)

			!Depol_22=0.8d0*dcos(tet/4.)*dcos(tet/4.)
			Depol_22=0. !0.8d0*(1.-dsin(tet/6.d0))

!			Rm(2,2)= Depol_22*alb_s*Ross_K*CB_peak_new(tet,kl_ext)*Int_alfbet1/mu0/muv*Shad*cos_2etv*cos_2eti+Rm(2,2)
!			Rm(2,2)= Rm(2,2)
            Rm(2,1)= Depol*P_Rel*cos_2etv*Int_alfbet1/mu0/muv*Shad+Rm(2,1)
            Rm(1,2)= Depol*P_Rel*cos_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,2)
            Rm(3,1)=-Depol*P_Rel*sin_2etv*Int_alfbet1/mu0/muv*Shad+Rm(3,1)
            Rm(1,3)= Depol*P_Rel*sin_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,3)  !!! Checked "-"
         end if ! n_el .ge. 3
         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_mat_3p_ocean_dep




!PL Subroutine for Ocen surface based on 
!PL Pavel Litvinov, Otto Hasekamp, Brian Cairns. Remote Sensing of Environment (2011), 115, 781–792
!PL (Experimental version)
Subroutine  Refl_mat_3p_ocean(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,m_ind,Rm)
!PL   Subroutine for Ocean surface based on
!     New physical model (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!     n_el: number of variables in Refl. Matr.
!     n_par - the number of model parameters
!     n_tet - the number of viewing zenith angle
!     tet0  - solar zenith angle
!     tet_v - array of viewing zenith angles
!     phi0  - azimuzth angles of solar direction
!     phiv  - azimuzth angles of viewing direction
!     x_vect- vector of input parameters
!     In given version the meaning of parameters is follows:
!        x_vect(1) is ocean abedo
!        x_vect(2) is Fraction of surface wich can provide Fresnel reflection
!        x_vect(3) is rms of surface slope (in scale 1)
!     Rm: calculated reflection matrix
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
	  Complex(8), intent(in) :: m_ind
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sig_min=0.02d0 !, sin_eps=0.00000001d0, kl_ext=90.d0

      Integer :: i1,i2
      Real(8) :: del_phi
      Real(8) :: Randf1(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,mu0,muv!,Int_bet1_Fr
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_ti,sin_tv !sin_tet,cos_etv, cos_eti, , sin_etv,sin_eti
!      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1 !,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi!,Depol_22

!      Real(8), parameter:: Depol=0.00d0

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

	  Real(8) mu_n0,f_pdf0
      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s     =x_vect(1)
         Fren_fr   =x_vect(2)
         sig1      =x_vect(3)
		 k_shad=1.

		 mu0=dcos(tet0)
         muv=dcos(tetv)

         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

         call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,m_ind,Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
        
         mu_n0 = (mu0+muv)/(2.d0*dsin(0.5d0*tet))
         f_pdf0=dexp(-(1.d0-mu_n0*mu_n0)/(2.d0*sig1*mu_n0*mu_n0))/(2.d0*dpi*sig1*mu_n0**3)

         Rm=Fren_fr*Shad*f_pdf0*Rm*dpi/(4.d0*mu0*muv*mu_n0)		 

         If (sig1 .ge. sig_min) then		 
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
          If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
          end if
 
          Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
                  Randf1(i1,i2)   = mu1v*mu2i*f_pdf1(i2)
               else
                  Randf1(i1,i2)   =0.d0
	           end if
            end do !i1
          end do !i2

          Int_alfbet1   =0.d0
          do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
          end do ! i2

		  Rm(1,1)=alb_s*Int_alfbet1/mu0/muv*Shad+Rm(1,1)
		 else
          Rm(1,1)=alb_s+Rm(1,1)
		 endif

         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_mat_3p_ocean



Subroutine  Refl_mat_4p_ocean(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,Rm)
!PL     New physical model (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!PL 	That is experimental version for future developments
!PL     n_par - the number of model parameters
!PL     n_tet - the number of viewing zenith angle
!PL     tet0  - solar zenith angle
!PL     tet_v - array of viewing zenith angles
!PL     phi0  - azimuzth angles of solar direction
!PL     phiv  - azimuzth angles of viewing direction
!PL     x_vect- vector of input parameters
!PL     In given version the meaning of parameters is follows:
!PL        x_vect(1) is abedo of a surface element
!PL        x_vect(2) is rms of surface slope (in scale 1)
!PL        x_vect(3) is Fraction of surface wich can provide Fresnel reflection
!PL        x_vect(4) is a parameter of shadowing function
!PL     An example of retrieved values of the parameters for soil (Mongu):
!PL        x_vect(1)=0.046 (440 nm)
!PL    	          =0.072 (490 nm)
!PL    	          =0.115 (550 nm)
!PL    	          =0.155 (670 nm)
!PL    	          =0.351 (865 nm)
!PL    	          =0.434 (1020 nm)
!PL        x_vect(2)=0.23
!PL        x_vect(3)=0.9
!PL        x_vect(4)=0.7
!PL     Rm: calculated reflection matrix
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sin_eps=0.00000001d0
      Real(8), parameter:: par_free=4.d0/3.d0,  alb_m=0.2d0!alb_m
!      Real(8), parameter:: Depol=0.005d0, sig2=0.005d0  !0.9!x_vect(5) ! for ocean
      Real(8), parameter:: Depol=0.00d0
	  Real(8) sig2  !0.9!x_vect(5)   ! for land

!      Real(8) alb_m

      Integer :: i1,i2
      Real(8) :: del_phi,Ross_K
      Real(8) :: Randf1(1:n_alf,1:n_bet), Randf1_fr(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,Int_bet1_Fr,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_tet,sin_ti,sin_tv, cos_etv,cos_eti, sin_etv,sin_eti
      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi,Randf0,Depol_22

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

!      Do i1=1,n_par
!         boolean=(boolean .and. (dabs(x_vect_old(i1)-x_vect(i1)) .le. eps))
!      end do
!XH   take advantage of new feature in Fortran95
      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s  =x_vect(1)
         sig1   =x_vect(2)
		 sig2   =x_vect(3)
         Fren_fr=1.!x_vect(3)
         select case(n_par)
         case(3)
            k_shad=sig1
         case(4)
            k_shad=x_vect(4)
         case default
            write(*,*) 'n_par=',n_par,' - unknown value'
            stop 'stop in Refl_mat_4p'
         end select

         mu0=dcos(tet0)
         muv=dcos(tetv)
 
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
         If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
         end if
 
         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

!XH      original computation for Randf1
!         Do i1=1,n_alf
!            Do i2=1,n_bet
!               mu1v=mu_beta(i2)*muv+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-muv*muv)*dcos(alfa(i1)-phiv)
!               mu2i=mu_beta(i2)*mu0+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-mu0*mu0)*dcos(alfa(i1)-phi0)
!               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
!                  call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
!                  Randf1(i1,i2)   =mu1v*mu2i*f_pdf1(i2)
!                  Randf1_fr(i1,i2)=mu1v*mu2i*f_pdf1(i2)*Randf1_fr(i1,i2)
!               else
!                  Randf1(i1,i2)   =0.d0
!                  Randf1_fr(i1,i2)=0.d0
!               end if
!            end do !i2
!         end do !i1

!XH      new computation for Randf1 to reduce redundant calculations
         Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
		          call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
                  Randf0          =mu1v*mu2i*f_pdf1(i2)
                  Randf1(i1,i2)   =Randf0
                  Randf1_fr(i1,i2)=Randf0*Randf1_fr(i1,i2)
               else
                  Randf1(i1,i2)   =0.d0
                  Randf1_fr(i1,i2)=0.d0
	           end if
            end do !i1
         end do !i2

!XH      original two dimensional integral
!         Int_alfbet1=0.d0
!         Int_alfbet1_Fr=0.d0
!
!         Do i1=1,n_alf
!            Int_bet1=0.d0
!            Int_bet1_Fr=0.d0
!            Do i2=1,n_bet
!               Int_bet1=Int_bet1+Gw_t(i2)*Randf1(i1,i2)
!               Int_bet1_Fr=Int_bet1_Fr+Gw_t(i2)*Randf1_Fr(i1,i2)
!            end Do !i2
!            Int_alfbet1=Int_alfbet1+Gw_p(i1)*Int_bet1
!            Int_alfbet1_Fr=Int_alfbet1_Fr+Gw_p(i1)*Int_bet1_Fr
!         end do !i1

!XH      new two dimensional integral taking advantage of the fact that n_alf>n_bet and n_alf is number of lines of Randf1
         Int_alfbet1   =0.d0
         Int_alfbet1_Fr=0.d0
         do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
            Int_bet1_Fr=sum(Gw_p*Randf1_fr(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
            Int_alfbet1_Fr=Int_alfbet1_Fr+Int_bet1_Fr*Gw_t(i2)
         end do ! i2

         call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,dcmplx(1.33d0,0.d0),Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
         Rm=Fren_fr*Rm*Int_alfbet1_Fr/mu0/muv*Shad
!         Ross_K=(1.-alb_m*dsin(tet/par_free)*dsin(tet/par_free))
!         Ross_K=1

!         Rm(1,1)=alb_s*Ross_K*CB_peak_new(tet,kl_ext)*Int_alfbet1/mu0/muv*Shad+Rm(1,1)
		 Rm(1,1)=alb_s*Int_alfbet1/mu0/muv*Shad+Rm(1,1)

         if (n_el .ge. 3) then
!XH         for vector case only
            sin_tet=dsqrt(1.d0-cos_tet*cos_tet)

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

            P_Rel=-(1.d0-cos_tet*cos_tet)

			!Depol_22=0.8d0*dcos(tet/4.)*dcos(tet/4.)
			Depol_22=0. !0.8d0*(1.-dsin(tet/6.d0))

!			Rm(2,2)= Depol_22*alb_s*Ross_K*CB_peak_new(tet,kl_ext)*Int_alfbet1/mu0/muv*Shad*cos_2etv*cos_2eti+Rm(2,2)
!			Rm(2,2)= Rm(2,2)
            Rm(2,1)= Depol*P_Rel*cos_2etv*Int_alfbet1/mu0/muv*Shad+Rm(2,1)
            Rm(1,2)= Depol*P_Rel*cos_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,2)
            Rm(3,1)=-Depol*P_Rel*sin_2etv*Int_alfbet1/mu0/muv*Shad+Rm(3,1)
            Rm(1,3)= Depol*P_Rel*sin_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,3)  !!! Checked "-"
         end if ! n_el .ge. 3
         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_mat_4p_ocean


!PL     New physical model (element R11) (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!PL     n_par - the number of model parameters
!PL     n_tet - the number of viewing zenith angle
!PL     tet0  - solar zenith angle
!PL     tet_v - array of viewing zenith angles
!PL     phi0  - azimuzth angles of solar direction
!PL     phiv  - azimuzth angles of viewing direction
!PL     x_vect- vector of input parameters
!PL     In given version the meaning of parameters is follows:
!PL        x_vect(1) is abedo of a surface element
!PL        x_vect(2) is rms of surface slope (in scale 1)
!PL        x_vect(3) is a parameter of shadowing function
!PL     An example of retrieved values of the parameters for soil (Mongu):
!PL        x_vect(1)=0.042 (440 nm)
!PL    	          =0.072 (490 nm)
!PL    	          =0.115 (550 nm)
!PL    	          =0.155 (670 nm)
!PL    	          =0.351 (865 nm)
!PL    	          =0.434 (1020 nm)
!PL        x_vect(2)=0.203
!PL        x_vect(3)=0.17
!PL     R11: calculated element R11 of the reflection matrix  
Subroutine  Refl_LitvBRM_R11(n_par,tet0,tetv,phi0,phiv,x_vect,R11)
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=5,n_alf=15! !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
      Real(8), intent (out):: R11

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sig_min=0.02d0
      Real(8), parameter:: par_free=4.d0/3.d0, alb_m=0.2d0!

      Integer :: i1,i2
      Real(8) :: R11_elem,R11_elem1,R11_elem2,R11_elem1_180,alb_m1,R11_elem2_180, &
	             R11_elem1_frac, R11_elem1_0, R11_elem2_0,R11_elem_0	
      Real(8) :: Randf1(1:n_alf,1:n_bet) !, Randf1_fr(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,k_shad,mu_n1
      Real(8) :: mu1v,mu2i,Shad
      Real(8) :: sin_ti,sin_tv 
      Real(8) :: Int_alfbet1
      Real(8) :: mubv,mubi,sinbv,sinbi

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old,f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)
      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         R11=Rm_old
      else
         alb_s  =x_vect(1)
         sig1   =x_vect(2)
!		 R11_elem1_frac=x_vect(3)
         select case(n_par)
         case(2)
          k_shad=1.
         case(3)
          k_shad=x_vect(3)
         case default
            write(*,*) 'n_par=',n_par,' - unknown value'
            stop 'stop in Refl_mat_4p'
         end select

         mu0=dcos(tet0)
         muv=dcos(tetv)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)
         R11_elem=(1.-alb_m*dsin(tet/par_free)*dsin(tet/par_free))

         If (sig1 .ge. sig_min) then
 !        If sig1 is the same then f_pdf1 is not calculated and taken from previous
 !        values since f_pdf1 is saved
          If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
          end if
 

          Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
                  Randf1(i1,i2) = mu1v*mu2i*f_pdf1(i2)
               else
                  Randf1(i1,i2)   =0.d0
	           end if
            end do !i1
          end do !i2

          Int_alfbet1   =0.d0
          do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
          end do ! i2

          R11=alb_s*R11_elem*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad
         else
		  R11=alb_s*R11_elem*CB_peak_new(tet)*Shad
         endif ! sig1

         Rm_old=R11
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean 
!
      return
     End Subroutine Refl_LitvBRM_R11



!PL     New physical model (element R11) (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!PL		Similar to Refl_LitvBRM_R11 with some experimental features
!PL     n_par - the number of model parameters
!PL     n_tet - the number of viewing zenith angle
!PL     tet0  - solar zenith angle
!PL     tet_v - array of viewing zenith angles
!PL     phi0  - azimuzth angles of solar direction
!PL     phiv  - azimuzth angles of viewing direction
!PL     x_vect- vector of input parameters
!PL     In given version the meaning of parameters is follows:
!PL        x_vect(1) is abedo of a surface element
!PL        x_vect(2) is rms of surface slope (in scale 1)
!PL        x_vect(3) is a parameter of shadowing function
!PL     An example of retrieved values of the parameters for soil (Mongu):
!PL        x_vect(1)=0.042 (440 nm)
!PL    	          =0.072 (490 nm)
!PL    	          =0.115 (550 nm)
!PL    	          =0.155 (670 nm)
!PL    	          =0.351 (865 nm)
!PL    	          =0.434 (1020 nm)
!PL        x_vect(2)=0.203
!PL        x_vect(3)=0.17
!PL     R11: calculated element R11 of the reflection matrix
Subroutine  Refl_LitvBRM_fast(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,m_ind,Rm)
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=5,n_alf=15! !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
	  Complex(8),intent(in) :: m_ind
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sin_eps=0.00000001d0
      Real(8), parameter:: par_free=4.d0/3.d0,  alb_m=0.2d0!alb_m
!      Real(8), parameter:: Depol=0.005d0, sig2=0.005d0  !0.9!x_vect(5) ! for ocean
      Real(8), parameter:: Depol=0.00d0, sig2=0.5 !0.9d0  !0.9!x_vect(5)   ! for land

!      Real(8) alb_m

      Integer :: i1,i2
      Real(8) :: del_phi,Ross_K
      Real(8) :: Randf1(1:n_alf,1:n_bet) !, Randf1_fr(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,Int_bet1_Fr,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_tet,sin_ti,sin_tv, cos_etv,cos_eti, sin_etv,sin_eti
      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi,Randf0,Depol_22

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

 	  Real(8) mu_n0,f_pdf0,sig3
      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

!      Do i1=1,n_par
!         boolean=(boolean .and. (dabs(x_vect_old(i1)-x_vect(i1)) .le. eps))
!      end do
!XH   take advantage of new feature in Fortran95
      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s  =x_vect(1)
         sig1   =x_vect(2)
         Fren_fr=x_vect(3)
         select case(n_par)
         case(3)
            k_shad=sig1
         case(4)
            k_shad=x_vect(4)
         case default
            write(*,*) 'n_par=',n_par,' - unknown value'
            stop 'stop in Refl_mat_4p'
         end select

         mu0=dcos(tet0)
         muv=dcos(tetv)
 
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
         If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
         end if
 
         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

!		 call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,m_ind,Rm)
!		 sig3=sig2+sig1
!         call Shad_func(mu0,muv,k_shad*sig3,Shad)
        
!         mu_n0 = (mu0+muv)/(2.d0*dsin(0.5d0*tet))
!         f_pdf0=dexp(-(1.d0-mu_n0*mu_n0)/(2.d0*sig3*mu_n0*mu_n0))/(2.d0*dpi*sig3*mu_n0**3)

!         Rm=Fren_fr*Shad*f_pdf0*Rm*dpi/(4.d0*mu0*muv*mu_n0)

!		 Cv=x_vect(n_par_i+1)
!      modified Maign and Breon model based on R_Fren_My
         call Maign_Breon_XH(n_el,mu0,muv,phiv/dpi_180,m_ind,Fren_fr,Rm)

		 

!XH      original computation for Randf1
!         Do i1=1,n_alf
!            Do i2=1,n_bet
!               mu1v=mu_beta(i2)*muv+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-muv*muv)*dcos(alfa(i1)-phiv)
!               mu2i=mu_beta(i2)*mu0+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-mu0*mu0)*dcos(alfa(i1)-phi0)
!               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
!                  call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
!                  Randf1(i1,i2)   =mu1v*mu2i*f_pdf1(i2)
!                  Randf1_fr(i1,i2)=mu1v*mu2i*f_pdf1(i2)*Randf1_fr(i1,i2)
!               else
!                  Randf1(i1,i2)   =0.d0
!                  Randf1_fr(i1,i2)=0.d0
!               end if
!            end do !i2
!         end do !i1

!XH      new computation for Randf1 to reduce redundant calculations
         Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
!		          call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
                  Randf0          =mu1v*mu2i*f_pdf1(i2)
                  Randf1(i1,i2)   =Randf0
!                  Randf1_fr(i1,i2)=Randf0*Randf1_fr(i1,i2)
               else
                  Randf1(i1,i2)   =0.d0
!                  Randf1_fr(i1,i2)=0.d0
	           end if
            end do !i1
         end do !i2

!XH      original two dimensional integral
!         Int_alfbet1=0.d0
!         Int_alfbet1_Fr=0.d0
!
!         Do i1=1,n_alf
!            Int_bet1=0.d0
!            Int_bet1_Fr=0.d0
!            Do i2=1,n_bet
!               Int_bet1=Int_bet1+Gw_t(i2)*Randf1(i1,i2)
!               Int_bet1_Fr=Int_bet1_Fr+Gw_t(i2)*Randf1_Fr(i1,i2)
!            end Do !i2
!            Int_alfbet1=Int_alfbet1+Gw_p(i1)*Int_bet1
!            Int_alfbet1_Fr=Int_alfbet1_Fr+Gw_p(i1)*Int_bet1_Fr
!         end do !i1

!XH      new two dimensional integral taking advantage of the fact that n_alf>n_bet and n_alf is number of lines of Randf1
         Int_alfbet1   =0.d0
!         Int_alfbet1_Fr=0.d0
         do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
!            Int_bet1_Fr=sum(Gw_p*Randf1_Fr(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
!            Int_alfbet1_Fr=Int_alfbet1_Fr+Int_bet1_Fr*Gw_t(i2)
         end do ! i2

!         call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,dcmplx(1.5d0,0.d0),Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
!         Rm=Fren_fr*Rm*Int_alfbet1_Fr/mu0/muv*Shad

        ! Rm=Rm*Shad
         Ross_K=(1.-alb_m*dsin(tet/par_free)*dsin(tet/par_free))
!         Ross_K=1

         Rm(1,1)=alb_s*Ross_K*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad+Rm(1,1)
         if (n_el .ge. 3) then
!XH         for vector case only
            sin_tet=dsqrt(1.d0-cos_tet*cos_tet)

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

!PL Resent  P_Rel=-(1.d0-cos_tet*cos_tet)

			!Depol_22=0.8d0*dcos(tet/4.)*dcos(tet/4.)
			Depol_22=0.8d0*(1.-dsin(tet/6.d0))

			Rm(2,2)= Depol_22*alb_s*Ross_K*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad*cos_2etv*cos_2eti+Rm(2,2)
!PL Resent             Rm(2,1)= Depol*P_Rel*cos_2etv*Int_alfbet1/mu0/muv*Shad+Rm(2,1)
!PL Resent             Rm(1,2)= Depol*P_Rel*cos_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,2)
!PL Resent             Rm(3,1)=-Depol*P_Rel*sin_2etv*Int_alfbet1/mu0/muv*Shad+Rm(3,1)
!PL Resent             Rm(1,3)= Depol*P_Rel*sin_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,3)  !!! Checked "-"
         end if ! n_el .ge. 3
         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_LitvBRM_fast



!PL     New physical model (reflection matrix  with some experimental features) 
!PL		(Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!PL     n_par - the number of model parameters
!PL     n_tet - the number of viewing zenith angle
!PL     tet0  - solar zenith angle
!PL     tet_v - array of viewing zenith angles
!PL     phi0  - azimuzth angles of solar direction
!PL     phiv  - azimuzth angles of viewing direction
!PL     x_vect- vector of input parameters
!PL     In given version the meaning of parameters is follows:
!PL        x_vect(1) is abedo of a surface element
!PL        x_vect(2) is rms of surface slope (in scale 1)
!PL        x_vect(3) is a parameter of shadowing function
!PL     An example of retrieved values of the parameters for soil (Mongu):
!PL        x_vect(1)=0.042 (440 nm)
!PL    	          =0.072 (490 nm)
!PL    	          =0.115 (550 nm)
!PL    	          =0.155 (670 nm)
!PL    	          =0.351 (865 nm)
!PL    	          =0.434 (1020 nm)
!PL        x_vect(2)=0.203
!PL        x_vect(3)=0.17
!PL     R11: calculated element R11 of the reflection matrix 
Subroutine  Refl_mat_4p_simp(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,m_ind,Rm)
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
	  Complex(8),intent(in) :: m_ind
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sin_eps=0.00000001d0
      Real(8), parameter:: par_free=4.d0/3.d0,  alb_m=0.2d0!alb_m
!      Real(8), parameter:: Depol=0.005d0, sig2=0.005d0  !0.9!x_vect(5) ! for ocean
      Real(8), parameter:: Depol=0.005d0, sig2=0.9 !0.9d0  !0.9!x_vect(5)   ! for land

!      Real(8) alb_m

      Integer :: i1,i2
      Real(8) :: del_phi,Ross_K
      Real(8) :: Randf1(1:n_alf,1:n_bet), Randf1_fr(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,Int_bet1_Fr,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_tet,sin_ti,sin_tv, cos_etv,cos_eti, sin_etv,sin_eti
      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi,Randf0,Depol_22

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

 	  Real(8) mu_n0,f_pdf0,sig3
      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

!      Do i1=1,n_par
!         boolean=(boolean .and. (dabs(x_vect_old(i1)-x_vect(i1)) .le. eps))
!      end do
!XH   take advantage of new feature in Fortran95
      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s  =x_vect(1)
         sig1   =x_vect(2)
         Fren_fr=x_vect(3)
         select case(n_par)
         case(3)
            k_shad=sig1
         case(4)
            k_shad=x_vect(4)
         case default
            write(*,*) 'n_par=',n_par,' - unknown value'
            stop 'stop in Refl_mat_4p'
         end select

         mu0=dcos(tet0)
         muv=dcos(tetv)
 
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
         If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
         end if
 
         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

		 call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,m_ind,Rm)
		 sig3=sig2 !+sig1
         call Shad_func(mu0,muv,k_shad*sig3,Shad)
        
         mu_n0 = (mu0+muv)/(2.d0*dsin(0.5d0*tet))
         f_pdf0=dexp(-(1.d0-mu_n0*mu_n0)/(2.d0*sig3*mu_n0*mu_n0))/(2.d0*dpi*sig3*mu_n0**3)

         Rm=Fren_fr*Shad*f_pdf0*Rm*dpi/(4.d0*mu0*muv*mu_n0)

!		 Cv=x_vect(n_par_i+1)
!      modified Maign and Breon model based on R_Fren_My
!         call Maign_Breon_XH(n_el,mu0,muv,phiv/dpi_180,m_ind,Fren_fr,Rm)

		 

!XH      original computation for Randf1
!         Do i1=1,n_alf
!            Do i2=1,n_bet
!               mu1v=mu_beta(i2)*muv+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-muv*muv)*dcos(alfa(i1)-phiv)
!               mu2i=mu_beta(i2)*mu0+dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*dsqrt(1.d0-mu0*mu0)*dcos(alfa(i1)-phi0)
!               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
!                  call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
!                  Randf1(i1,i2)   =mu1v*mu2i*f_pdf1(i2)
!                  Randf1_fr(i1,i2)=mu1v*mu2i*f_pdf1(i2)*Randf1_fr(i1,i2)
!               else
!                  Randf1(i1,i2)   =0.d0
!                  Randf1_fr(i1,i2)=0.d0
!               end if
!            end do !i2
!         end do !i1

!XH      new computation for Randf1 to reduce redundant calculations
         Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
!		          call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
                  Randf0          =mu1v*mu2i*f_pdf1(i2)
                  Randf1(i1,i2)   =Randf0
!                  Randf1_fr(i1,i2)=Randf0*Randf1_fr(i1,i2)
               else
                  Randf1(i1,i2)   =0.d0
!                  Randf1_fr(i1,i2)=0.d0
	           end if
            end do !i1
         end do !i2

!XH      original two dimensional integral
!         Int_alfbet1=0.d0
!         Int_alfbet1_Fr=0.d0
!
!         Do i1=1,n_alf
!            Int_bet1=0.d0
!            Int_bet1_Fr=0.d0
!            Do i2=1,n_bet
!               Int_bet1=Int_bet1+Gw_t(i2)*Randf1(i1,i2)
!               Int_bet1_Fr=Int_bet1_Fr+Gw_t(i2)*Randf1_Fr(i1,i2)
!            end Do !i2
!            Int_alfbet1=Int_alfbet1+Gw_p(i1)*Int_bet1
!            Int_alfbet1_Fr=Int_alfbet1_Fr+Gw_p(i1)*Int_bet1_Fr
!         end do !i1

!XH      new two dimensional integral taking advantage of the fact that n_alf>n_bet and n_alf is number of lines of Randf1
         Int_alfbet1   =0.d0
!         Int_alfbet1_Fr=0.d0
         do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
!            Int_bet1_Fr=sum(Gw_p*Randf1_Fr(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
!            Int_alfbet1_Fr=Int_alfbet1_Fr+Int_bet1_Fr*Gw_t(i2)
         end do ! i2

!         call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,dcmplx(1.5d0,0.d0),Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
         Rm=Rm*Int_alfbet1_Fr*Shad

        ! Rm=Rm*Shad
         Ross_K=(1.-alb_m*dsin(tet/par_free)*dsin(tet/par_free))
!         Ross_K=1

         Rm(1,1)=alb_s*Ross_K*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad+Rm(1,1)
         if (n_el .ge. 3) then
!XH         for vector case only
            sin_tet=dsqrt(1.d0-cos_tet*cos_tet)

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

            P_Rel=-(1.d0-cos_tet*cos_tet)

			!Depol_22=0.8d0*dcos(tet/4.)*dcos(tet/4.)
			Depol_22=0.8d0*(1.-dsin(tet/6.d0))

			Rm(2,2)= Depol_22*alb_s*Ross_K*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad*cos_2etv*cos_2eti+Rm(2,2)
            Rm(2,1)= Depol*P_Rel*cos_2etv*Int_alfbet1/mu0/muv*Shad+Rm(2,1)
            Rm(1,2)= Depol*P_Rel*cos_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,2)
            Rm(3,1)=-Depol*P_Rel*sin_2etv*Int_alfbet1/mu0/muv*Shad+Rm(3,1)
            Rm(1,3)= Depol*P_Rel*sin_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,3)  !!! Checked "-"
         end if ! n_el .ge. 3
         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_mat_4p_simp



!PL     New physical model (Litvinov et all. JQSRT (2012), 113, 2023-2039. doi:10.1016/j.jqsrt.2012.06.027.)
!PL     n_el: number of elements of the Refl. Matrix
!PL     n_par - the number of model parameters
!PL     n_tet - the number of viewing zenith angle
!PL     tet0  - solar zenith angle
!PL     tet_v - array of viewing zenith angles
!PL     phi0  - azimuzth angles of solar direction
!PL     phiv  - azimuzth angles of viewing direction
!PL     x_vect- vector of input parameters
!PL     In given version the meaning of parameters is follows:
!PL        x_vect(1) is abedo of a surface element
!PL        x_vect(2) is rms of surface slope (in scale 1)
!PL        x_vect(3) is Fraction of surface wich can provide Fresnel reflection
!PL        x_vect(4) is a parameter of shadowing function
!PL     An example of retrieved values of the parameters for soil (Mongu):
!PL        x_vect(1)=0.042 (440 nm)
!PL    	          =0.072 (490 nm)
!PL    	          =0.115 (550 nm)
!PL    	          =0.155 (670 nm)
!PL    	          =0.351 (865 nm)
!PL    	          =0.434 (1020 nm)
!PL        x_vect(2)=0.203
!PL        x_vect(3)=0.635
!PL        x_vect(4)=0.17
!PL     An example of retrieved values of the parameters for Banizoumbou:
!PL        x_vect(1)=0.072 (440 nm)
!PL    	          =0.126 (490 nm)
!PL    	          =0.229 (550 nm)
!PL    	          =0.371(670 nm)
!PL    	          =0.510 (865 nm)
!PL    	          =0.590 (1020 nm)
!PL        x_vect(2)=0.115
!PL        x_vect(3)=0.725
!PL        x_vect(4)=0.1
!PL     Rm: calculated reflection matrix

Subroutine  Refl_mat_4p(n_el,n_par,tet0,tetv,phi0,phiv,x_vect,m_ind,Rm)
      Implicit none
      Integer, parameter   :: n_bet=9,n_alf=22 !n_bet=20,n_alf=40!n_bet=8,n_alf=12 !n_bet=9,n_alf=25
      Integer, intent (in) :: n_el,n_par
      Real(8), intent (in) :: tetv,tet0,phi0,phiv,x_vect(1:n_par_max)
	  Complex(8),intent(in) :: m_ind
      Real(8), intent (out):: Rm(1:n_el,1:n_el)

!      Real(8),parameter:: sin_eps=0.00000001!, k_shad=1.6

      real(8), parameter:: eps=0.00000000001d0, sin_eps=0.00000001d0, sig_min=0.02d0
      Real(8), parameter:: par_free=4.d0/3.d0,  alb_m=0.2d0!alb_m
!      Real(8), parameter:: Depol=0.005d0, sig2=0.005d0  !0.9!x_vect(5) ! for ocean
      Real(8), parameter:: Depol=0.005d0, sig2=0.9 !0.9d0  !0.9!x_vect(5)   ! for land

!      Real(8) alb_m

      Integer :: i1,i2
      Real(8) :: del_phi,R11_elem
      Real(8) :: Randf1(1:n_alf,1:n_bet), Randf1_fr(1:n_alf,1:n_bet)
      Real(8) :: low_lp,up_lp
      Real(8) :: Int_bet1,Int_bet1_Fr,mu0,muv
      Real(8) :: cos_tet,tet
      Real(8) :: sig1,alb_s,Fren_fr,k_shad,mu_n1
      Real(8) :: phii,mu1v,mu2i,Shad, P_Rel
      Real(8) :: sin_tet,sin_ti,sin_tv, cos_etv,cos_eti, sin_etv,sin_eti
      Real(8) :: cos_2etv,sin_2etv, cos_2eti,sin_2eti
      Real(8) :: Int_alfbet1,Int_alfbet1_Fr
      Real(8) :: mubv,mubi,sinbv,sinbi,Randf0,Depol_22,mu_n0,f_pdf0

      Integer, save:: i_num
!     These parameter must be calculated just once
      Real(8), save:: tet0_old,tetv_old,phi0_old,phiv_old,x_vect_old(1:n_par_max),&
                      Rm_old(1:4,1:4),f_pdf1(1:n_bet)
      Real(8), save:: alfa(1:n_alf),mu_beta(1:n_bet),beta(1:n_bet),Gw_p(1:n_alf),Gw_t(1:n_bet)

      logical boolean
!
      low_lp=0.
      up_lp=2.*dpi

!     Here we use the fact that variable i_num is saved in the memory
!     and we need to calculate Gaussian quadratures just one time.
      If (i_num .eq. 0) then
         call gauss_l(n_alf,n_alf,low_lp,up_lp,alfa,Gw_p)
         call gauss_l(n_bet,n_bet,0.d0,1.d0,mu_beta,Gw_t)
         beta=dacos(mu_beta)
         x_vect_old=0.d0
         tet0_old  =0.d0
         tetv_old  =0.d0
         phi0_old  =0.d0
         phiv_old  =0.d0
         i_num=1
      end if
      If (i_num .eq. 0) then
         write(*,*) 'something wrong'
         stop
      end if

      boolean=((dabs(tet0_old-tet0) .le. eps) .and. (dabs(tetv_old-tetv) .le. eps) &
         .and. (dabs(phi0_old-phi0) .le. eps) .and. (dabs(phiv_old-phiv) .le. eps))

!      Do i1=1,n_par
!         boolean=(boolean .and. (dabs(x_vect_old(i1)-x_vect(i1)) .le. eps))
!      end do
!XH   take advantage of new feature in Fortran95
      boolean = (boolean .and. all(dabs(x_vect_old-x_vect) .le. eps))
 	  
      If (boolean) then
         Rm=Rm_old
      else
         alb_s  =x_vect(1)
         sig1   =x_vect(2)
         Fren_fr=x_vect(3)
         select case(n_par)
         case(3)
            k_shad=0.01!sig1
         case(4)
            k_shad=x_vect(4)
         case default
            write(*,*) 'n_par=',n_par,' - unknown value'
            stop 'stop in Refl_mat_4p'
         end select

         mu0=dcos(tet0)
         muv=dcos(tetv)

         phii=phi0+dpi
         del_phi=phiv-phii
         sin_ti=dsqrt(1.d0-mu0*mu0)
         sin_tv=dsqrt(1.d0-muv*muv)
         cos_tet=-(muv*mu0+sin_tv*sin_ti*dcos(phiv-phi0))
         tet=dacos(cos_tet)

         call R_Fren_My(n_el,mu0,muv,phiv/dpi_180,m_ind,Rm)
         call Shad_func(mu0,muv,k_shad*sig1,Shad)
         R11_elem=(1.-alb_m*dsin(tet/par_free)*dsin(tet/par_free))
		 Rm=Fren_fr*Rm*Shad/mu0/muv

         P_Rel=-(1.d0-cos_tet*cos_tet)
         sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
		!Depol_22=0.8d0*dcos(tet/4.)*dcos(tet/4.)
         Depol_22=0.8d0*(1.-dsin(tet/6.d0))


		 !If the surface is not flat in sclae1
         If (sig1 .ge. sig_min) then
!        If sig1 is the same then f_pdf1 is not calculated and taken from previous
!        values since f_pdf1 is saved
          If (dabs(x_vect_old(2)-x_vect(2)) .gt. eps) then
            Do i2=1,n_bet
               mu_n1=mu_beta(i2)
               f_pdf1(i2)=dexp(-(1.-mu_n1*mu_n1)/(2.d0*sig1)/mu_n1/mu_n1)/ &
	                               dpi/(2.*sig1)/mu_n1/mu_n1/mu_n1/mu_n1
            end do
          end if

          Do i2=1,n_bet
            mubv=mu_beta(i2)*muv
            mubi=mu_beta(i2)*mu0
            sinbv=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_tv
            sinbi=dsqrt(1.d0-mu_beta(i2)*mu_beta(i2))*sin_ti
            Do i1=1,n_alf
               mu1v=mubv+sinbv*dcos(alfa(i1)-phiv)
               mu2i=mubi+sinbi*dcos(alfa(i1)-phi0)
               if ((mu1v .ge. 0.) .and. (mu2i .ge. 0.)) then
		          call Pdf_Fren(alfa(i1),beta(i2),sig2,mu0,phii,muv,phiv,Randf1_fr(i1,i2),mu1v,mu2i)
                  Randf0          =mu1v*mu2i*f_pdf1(i2)
                  Randf1(i1,i2)   =Randf0
                  Randf1_fr(i1,i2)=Randf0*Randf1_fr(i1,i2)
               else
                  Randf1(i1,i2)   =0.d0
                  Randf1_fr(i1,i2)=0.d0
	           end if
            end do !i1
          end do !i2

          Int_alfbet1   =0.d0
          Int_alfbet1_Fr=0.d0
          do i2=1,n_bet
            Int_bet1   =sum(Gw_p*Randf1(:,i2))
            Int_bet1_Fr=sum(Gw_p*Randf1_Fr(:,i2))
            Int_alfbet1   =Int_alfbet1   +Int_bet1   *Gw_t(i2)
            Int_alfbet1_Fr=Int_alfbet1_Fr+Int_bet1_Fr*Gw_t(i2)
          end do ! i2

          Rm=Rm*Int_alfbet1_Fr
          Rm(1,1)=alb_s*R11_elem*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad+Rm(1,1)
          if (n_el .ge. 3) then 

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

			Rm(2,2)= Depol_22*alb_s*R11_elem*CB_peak_new(tet)*Int_alfbet1/mu0/muv*Shad*cos_2etv*cos_2eti+Rm(2,2)
            Rm(2,1)= Depol*P_Rel*cos_2etv*Int_alfbet1/mu0/muv*Shad+Rm(2,1)
            Rm(1,2)= Depol*P_Rel*cos_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,2)
            Rm(3,1)=-Depol*P_Rel*sin_2etv*Int_alfbet1/mu0/muv*Shad+Rm(3,1)
            Rm(1,3)= Depol*P_Rel*sin_2eti*Int_alfbet1/mu0/muv*Shad+Rm(1,3)  !!! Checked "-"
          end if ! n_el .ge. 3
	     else
           
          Rm=f_pdf0*Rm*dpi/4.d0/mu_n0
          Rm(1,1)=alb_s*R11_elem*CB_peak_new(tet)*Shad+Rm(1,1)
          if (n_el .ge. 3) then 

            if (dabs(sin_tet) .lt. sin_eps) sin_tet=sin_eps
            if (dabs(sin_tv)  .lt. sin_eps) sin_tv= sin_eps
            if (dabs(sin_ti)  .lt. sin_eps) sin_ti= sin_eps
 
            cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
            cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

            sin_etv=dsin(del_phi)*sin_ti/sin_tet
            sin_eti=dsin(del_phi)*sin_tv/sin_tet

            cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
            sin_2etv=2.d0*cos_etv*sin_etv

            cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
            sin_2eti=2.d0*cos_eti*sin_eti

			Rm(2,2)= Depol_22*alb_s*R11_elem*CB_peak_new(tet)*Shad*cos_2etv*cos_2eti+Rm(2,2)
            Rm(2,1)= Depol*P_Rel*cos_2etv*Shad+Rm(2,1)
            Rm(1,2)= Depol*P_Rel*cos_2eti*Shad+Rm(1,2)
            Rm(3,1)=-Depol*P_Rel*sin_2etv*Shad+Rm(3,1)
            Rm(1,3)= Depol*P_Rel*sin_2eti*Shad+Rm(1,3)  !!! Checked "-"
          end if ! n_el .ge. 3
		 endif !(sig1)

         Rm_old(1:n_el,1:n_el)=Rm
         x_vect_old=x_vect
         tet0_old  =tet0
         tetv_old  =tetv
         phi0_old  =phi0
         phiv_old  =phiv
      end if !boolean
!
      return
End Subroutine Refl_mat_4p



!PL     BPDF model developted at SRON (Litinov et all.  Remote Sensing of Environment (2011), 115, 781–792;
!PL     Litinov et all.In book: Polarimetric Detection, Characterization, and Remote Sensing. Springer, 2011.)
!PL     Calculation of the Fresnel Reflection matrix is performed using M.Mishchenko code for
!PL     Gaussian random rough surface
!PL        M. I. Mishchenko and L. D. Travis, Satellite retrieval
!PL        of aerosol properties over the ocean using polarization as well as
!PL        intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-
!PL        17013 (1997).
!PL     The new shadowing function is introduced and the dependence on dmui and dmur is modified
!PL     (1/(dmui+dmur) is used instead of 1/dmui/dmur)
!PL        x_vect_p(1) = Re(m)
!PL        x_vect_p(2) = sigma^2
!PL        x_vect_p(3) = alpha
!PL        x_vect_p(4) = k_tet
!PL     Usually parameters x_vect_p(1) and x_vect_p(4) can be fixed.
!PL     Thus just 2 parameters can be retrieved.
!PL     An example of the retrieved parameters from PARASOL/POLDER3:
!PL        x_vect_p(2)=0.146
!PL        x_vect_p(3)=0.8
Subroutine R_Fren_PL_SRON(n_el,n_par_p,dmui,phii,dmur,Phir,x_vect_p,R)
      IMPLICIT NONE
      Integer, intent(in) :: n_el,n_par_p
      Real(8), intent(in) :: x_vect_p(1:n_par_p)
      Real(8), intent(out):: R(n_el,n_el)
      Real(8), intent(in) :: dmui,phii,dmur,Phir
!
      Real(8)    :: sigma2,shad,darg,dcos_scat,dtet_sc
      Complex(8) :: cn1,cn2

      cn1=dcmplx(1.d0,0.d0)
      cn2=dcmplx(x_vect_p(1),0.d0)
      Sigma2=x_vect_p(2)

!XH   calculate n_el x n_el Fresnel reflection matrix
      call R_Fren(n_el,dmui,phii,dmur,phir,cn1,cn2,R,sigma2)
 
!     New Shadowing function
      dcos_scat=dmui*Dmur-dsqrt(1.d0-dmui**2)*dsqrt(1.D0-dmur**2)*dcos(phii-phir)
!XH   prevent computational error
      if (dcos_scat .gt. 1.d0) then
         dcos_scat=1.d0
         dtet_sc  =0.d0
      else if (dcos_scat .lt. -1.d0) then
         dcos_scat=-1.d0
         dtet_sc  = dpi
      else
         dtet_sc=dacos(dcos_scat)
      end if
      darg=x_vect_p(4)!*dtet_sc
      Shad=((1.d0+dcos(darg))/2.d0)**3
      if (abs(darg) .gt. dpi) Shad=0.d0
!     end New Fr with shad

      R=Shad*x_vect_p(3)*dmui*dmur/(dmui+dmur)*R
!      R=Shad*x_vect_p(3)*R
!
      return
End Subroutine R_Fren_PL_SRON


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PL  Probability Density Function (PDF) for Gaussian rough surfaces
Subroutine Pdf_Fren(alf,bet,sigma2,dmui1,phii,dmur1,Phir,Pdf_fr,mu1r,mu2i)
  implicit none
 ! IMPLICIT DOUBLE PRECISION (A-B,D-H,O-Z), COMPLEX*16 (C)
  !Integer, intent(in):: n_par_p
  !Double precision, intent(in):: x_vect_p(1:n_par_p)
  Double precision, intent(out):: PDF_fr
  Double precision, intent(in)::alf,bet,sigma2,phii,Phir,mu1r,mu2i
  Double precision, intent(inout):: dmui1,dmur1
  Double precision dmui,dmur
  Double precision S1,S2,S3,mu_ir,mu_sq,cot,T1,T2,P,ShadI,ShadR,Shad, &
                   fact,fact1
  Double precision cos_phi_inc,sin_phi_inc,cos_phi_r,sin_phi_r,sin_inc_ang,sin_ref_ang,vi1o,vi2o,vi3o,vr1o, &
                   vr2o,vr3o,vi1,vi2,vi3,vr1,vr2,vr3, &
                   unit1,unit2,unit3,mu_n2,mu_n4,coeff,d_exp,r_rot(3,3)!,pi
  
!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

   !pi = 2.d0*Dasin(1.d0)
   

   IF (DABS(dmui1-1D0).LT.1D-10) dmui1=0.999999999999D0
   IF (DABS(dmur1-1D0).LT.1D-10) dmur1=0.999999999999D0
   
   cos_phi_inc=DCOS(PHII)
   sin_phi_inc=DSIN(PHII)
   cos_phi_r=DCOS(PHIR)
   sin_phi_r=DSIN(PHIR)

   sin_inc_ang=DSQRT(1D0-dmui1*dmui1)
   sin_ref_ang=DSQRT(1D0-dmur1*dmur1)
 
   !coordinates of incident normal
   VI1o=sin_inc_ang*cos_phi_inc
   VI2o=sin_inc_ang*sin_phi_inc
   VI3o=-dmui1

   !coordinates of viewing normal
   VR1o=sin_ref_ang*cos_phi_r
   VR2o=sin_ref_ang*sin_phi_r
   VR3o=dmur1


!rotation matrix
   r_rot(1,1)=dcos(alf)*dcos(bet)
   r_rot(1,2)=dsin(alf)*dcos(bet)
   r_rot(1,3)=-dsin(bet)

   r_rot(2,1)=-dsin(alf)
   r_rot(2,2)=dcos(alf)
   r_rot(2,3)=0.d0

   r_rot(3,1)=dcos(alf)*dsin(bet)
   r_rot(3,2)=dsin(alf)*dsin(bet)
   r_rot(3,3)=dcos(bet)

!transformation on rotation
   VI1=r_rot(1,1)*VI1o+r_rot(1,2)*VI2o+r_rot(1,3)*VI3o
   VI2=r_rot(2,1)*VI1o+r_rot(2,2)*VI2o+r_rot(2,3)*VI3o
   VI3=r_rot(3,1)*VI1o+r_rot(3,2)*VI2o+r_rot(3,3)*VI3o

!   Write(*,'(2f10.5,3f20.16)') dacos(dmui1)*d180_pi,phii*d180_pi,VI1o,VI2o,VI3o
   !local solar zenith angle
   dmui=-VI3

!!!!!!!!!!!!!!!!!1  
 !  dmui=mu2i
   if (dabs(dmui-mu2i) .gt. 0.00001) write(*,*) 'Error dmui 1', dmui,mu2i


   if (dmui .lt. 0) Write(*,*) 'Error dmui 2'
   IF (DABS(dmui-1D0).LT.1D-10) dmui=0.999999999999D0
   sin_inc_ang=DSQRT(1D0-dmui*dmui)


   !local solar azimuth angle
   cos_phi_inc=VI1/sin_inc_ang
   sin_phi_inc=VI2/sin_inc_ang
!!!!!!!!!!!!!!!!!

   if (dabs(dabs(sin_phi_inc)-dsqrt(1-cos_phi_inc**2)) .gt. 0.00001) then
    Write(*,*) 'Error sin_phi_inc', sin_phi_inc,dsqrt(1-cos_phi_inc**2)
   end if

   VR1=r_rot(1,1)*VR1o+r_rot(1,2)*VR2o+r_rot(1,3)*VR3o
   VR2=r_rot(2,1)*VR1o+r_rot(2,2)*VR2o+r_rot(2,3)*VR3o
   VR3=r_rot(3,1)*VR1o+r_rot(3,2)*VR2o+r_rot(3,3)*VR3o

   !local viewing zenith angle
   dmur=VR3


!!!!!!!!!!!!!!!!!
   if (dabs(DMUr-mu1r) .gt. 0.00001) write(*,*) 'Error dmur', DMUr,mu1r

   if (dmur .lt. 0) Write(*,*) 'Error dmur'
   IF (DABS(dmur-1D0).LT.1D-10) dmur=0.999999999999D0
   sin_ref_ang=DSQRT(1.D0-dmur*dmur)
   !local viewing azimuth angle
   cos_phi_r=VR1/sin_ref_ang
   sin_phi_r=VR2/sin_ref_ang

   if (dabs(dabs(sin_phi_r)-dsqrt(1.-cos_phi_r**2)) .gt. 0.0001) then
    Write(*,*) 'Error sin_phi_r',dacos(dmur1)*d180_pi, sin_phi_r,dsqrt(1-cos_phi_r**2)
    Write(*,'(4d16.5)') VR2,VR1o,VR2o,VR3o
   end if

   UNIT1=VI1-VR1
   UNIT2=VI2-VR2
   UNIT3=VI3-VR3
   fact1=UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
   fact=DSQRT(1D0/fact1)
    
   mu_n2=UNIT3*UNIT3
   mu_n4=mu_n2*mu_n2
 
   coeff=1D0/(4D0*dmui*dmur*mu_n4*2D0*sigma2)  !original
   d_exp= -(UNIT1*UNIT1 + UNIT2*UNIT2)/(2D0*sigma2*mu_n2)
   d_exp=DEXP(d_exp)
   coeff=coeff*fact1*fact1*d_exp   !!!! mu_n2^4  !original

 
!  SHADOWING
   P=DACOS(-1D0)
   S1=DSQRT(2D0*sigma2/P)
   S3=1D0/(DSQRT(2D0*sigma2))
   S2=S3*S3
   mu_ir=dmui
   mu_sq=mu_ir*mu_ir
   cot=mu_ir/DSQRT(1D0-mu_sq)
   T1=DEXP(-cot*cot*S2)

!   T2=erfc_func(cot*S3)

   T2=DERFC(cot*S3)


!	T2=ERFC(cot*S3)
!	T2=1-ERF(cot*S3)
   ShadI=0.5D0*(S1*T1/cot-T2)
   mu_ir=dmur
   mu_sq=mu_ir*mu_ir
   cot=mu_ir/DSQRT(1D0-mu_sq)
   T1=DEXP(-cot*cot*S2)

!   T2=erfc_func(cot*S3)
   T2=DERFC(cot*S3)
!      T2=ERFC(cot*S3)
!	T2=1-ERF(cot*S3)
   ShadR=0.5D0*(S1*T1/cot-T2)
   Shad=1D0/(1D0+ShadI+ShadR)

   Pdf_fr=coeff*Shad

End Subroutine Pdf_Fren



!PL Fresnel reflection from Kirchhoff approximation for Guassian surface
!PL cn1 and cn2 are complex refractive index for two media
Subroutine R_Fren(n_el,dmui1,phii,dmur1,Phir,cn1,cn2,R,sigma2)
!XH   provide sigma2 for roughened surface; otherwise for flat surface
      IMPLICIT NONE
      integer,    intent(in) :: n_el
      Real(8),    intent(in) :: phii,Phir
      Real(8),    intent(in) :: dmui1,dmur1
      Real(8),    intent(in), optional :: sigma2
      Complex(8), intent(in) :: cn1,cn2
      Real(8),    intent(out):: R(n_el,n_el)
!
      Complex(8) :: fr_c2,c1,c2,perp,paral,CF11,CF12,CF21,CF22,ci, &
                    c21,c22,a1,a2,a3,a4,a5,a6
      Real(8)    :: cos_phi_inc,sin_phi_inc,cos_phi_r,sin_phi_r,sin_inc_ang,sin_ref_ang,vi1,vi2,vi3,vr1,vr2,vr3, &
                    unit1,unit2,unit3,fact1,fact,fr_c1,t_inc1,t_inc2,t_inc3,t_ref1,t_ref2,t_ref3,pi1,pi2,pi3,pr1,pr2,pr3, &
                    p_1,p_2,p_3,p_4,E1,E2,E3,E4,VP1,VP2,VP3,coeff,d_exp,F,F11,F12,F21,F22,module2,dmui,dmur,mu_n2,mu_n4
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------
      dmui=dmui1
      dmur=dmur1

      IF (DABS(dmui1-1.D0).LT.1.D-10) dmui=0.999999999999D0
      IF (DABS(dmur1-1.D0).LT.1.D-10) dmur=0.999999999999D0

      cos_phi_inc=DCOS(PHII)
      sin_phi_inc=DSIN(PHII)
      cos_phi_r=DCOS(PHIR)
      sin_phi_r=DSIN(PHIR)

      sin_inc_ang=DSQRT(1.D0-dmui*dmui)
      sin_ref_ang=DSQRT(1.D0-dmur*dmur)
      VI1=sin_inc_ang*cos_phi_inc
      VI2=sin_inc_ang*sin_phi_inc
      VI3=-dmui
      VR1=sin_ref_ang*cos_phi_r
      VR2=sin_ref_ang*sin_phi_r
      VR3=dmur
 
!     Local normal to surface (for spec. refl.)
      UNIT1=VI1-VR1
      UNIT2=VI2-VR2
      UNIT3=VI3-VR3
      fact1=UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
      fact=DSQRT(1.D0/fact1)

!     Fresnel coefficients
      fr_c1=fact*(UNIT1*VI1+UNIT2*VI2+UNIT3*VI3)
      fr_c2=1.D0 - (1.D0-fr_c1*fr_c1)*cn1*cn1/(cn2*cn2)
      fr_c2=CDSQRT(fr_c2)
      C1=cn1*fr_c1
      C2=cn2*fr_c2
      perp=(C1-C2)/(C1+C2)
      C1=cn2*fr_c1
      C2=cn1*fr_c2
      paral=(C1-C2)/(C1+C2)
 
!     Amplitude scattering matrix
      t_inc1=-dmui*cos_phi_inc
      t_inc2=-dmui*sin_phi_inc
      t_inc3=-sin_inc_ang
 
      t_ref1= dmur*cos_phi_r
      t_ref2= dmur*sin_phi_r
      t_ref3=-sin_ref_ang
 
      PI1=-sin_phi_inc
      PI2= cos_phi_inc
      PI3= 0.D0

      PR1=-sin_phi_r
      PR2= cos_phi_r
      PR3= 0.D0
 
      p_1=PI1*VR1+PI2*VR2+PI3*VR3
      p_2=PR1*VI1+PR2*VI2+PR3*VI3
      p_3=t_inc1*VR1+t_inc2*VR2+t_inc3*VR3
      p_4=t_ref1*VI1+t_ref2*VI2+t_ref3*VI3
 
      E1=p_1*p_2
      E2=p_3*p_4
      E3=p_3*p_2
      E4=p_1*p_4
 
      CF11= E1*perp+E2*paral
      CF12=-E3*perp+E4*paral
      CF21=-E4*perp+E3*paral
      CF22= E2*perp+E1*paral
 
      VP1=VI2*VR3-VI3*VR2
      VP2=VI3*VR1-VI1*VR3
      VP3=VI1*VR2-VI2*VR1
      module2=VP1*VP1+VP2*VP2+VP3*VP3
      module2=module2*module2

      if (present(sigma2)) then
!XH      roughened surface with RMS slope sigma2
         mu_n2=UNIT3*UNIT3
         mu_n4=mu_n2*mu_n2
         coeff=1.D0/(4.D0*dmui*dmur*module2*mu_n4*2.D0*sigma2)
!         coeff=1.D0/(4.D0*(dmui+dmur)*module2*mu_n4*2D0*sigma2)
         d_exp= -(UNIT1*UNIT1 + UNIT2*UNIT2)/(2.D0*sigma2*mu_n2)
         d_exp=DEXP(d_exp)
         coeff=coeff*fact1*fact1*d_exp
      else
!XH      flat surface
         coeff=1.D0/module2
      end if

      F=0.5D0*coeff

      F11=CDABS(CF11)
      F12=CDABS(CF12)
      F21=CDABS(CF21)
      F22=CDABS(CF22)
      F11=F11*F11
      F12=F12*F12
      F21=F21*F21
      F22=F22*F22

!XH   1 x 1
      R(1,1)=(F11+F12+F21+F22)*F
      if (n_el .ge. 3) then
!XH      3 x 3
         R(1,2)=(F11-F12+F21-F22)*F
         R(2,1)=(F11-F22+F12-F21)*F
         R(2,2)=(F11-F12-F21+F22)*F
 
         CI=(0.D0, -1.D0)
         C21=DCONJG(CF21)
         C22=DCONJG(CF22)
         a1=CF11*DCONJG(CF12)
         a2=CF11*C21
         a3=CF11*C22
         a4=CF12*C21
         a5=CF12*C22
         a6=CF21*C22
 
         R(1,3)=    (-a1-a6)*coeff
         R(2,3)=    (-a1+a6)*coeff
         R(3,1)=    (-a2-a5)*coeff
         R(3,2)=    (-a2+a5)*coeff
         R(3,3)=    ( a3+a4)*coeff
!         if (n_el .eq. 4) then
!!XH         4 x 4
!            R(1,4)=-CI*( a1+a6)*coeff
!            R(2,4)=-CI*( a1-a6)*coeff
!            R(3,4)= CI*( a3-a4)*coeff
!            R(4,1)= CI*( a2+a5)*coeff
!            R(4,2)= CI*( a2-a5)*coeff
!            R(4,3)=-CI*( a3+a4)*coeff
!            R(4,4)=    ( a3-a4)*coeff
!         end if
      end if
 !
      return
End Subroutine R_Fren


!PL Fresnel reflection for Guassian surface (the same as R_Fren for flat surface)
!PL m is relative complex refractive index of the medium
Subroutine R_Fren_My(n_el,mu01,muv1,phi,m,R_fr,alf1)
!XH   add variable n_el, output alf if needed
      Implicit none
      integer,    intent(in) :: n_el
      Complex(8), intent(in) :: m
      Real(8),    intent(in) :: mu01,muv1,phi
      Real(8),    intent(out):: R_fr(1:n_el,1:n_el)
      Real(8),    intent(out), optional :: alf1
  
      Real(8), parameter :: st_min=0.001,sin_eps=0.00000001
      Real(8)    :: zz,alf,mu_alf
      Complex(8) :: r1,r2,r_par1,r_par2,r_perp1,r_perp2,r_par,r_perp,&
                    r_par_per1,r_par_per2
      Real(8)    :: r_par_sq,r_perp_sq
      Real(8)    :: del_phi,sin_ti,sin_tv,sin_tet,cos_tet,cos_etv,cos_eti, &
                    sin_eti,sin_etv,cos_2etv,sin_2etv,cos_2eti,sin_2eti
      Real(8)    :: L_v(1:n_el,1:n_el), L_inc(1:n_el,1:n_el), F_fr(1:n_el,1:n_el),mu0,muv,tet0v,sin_phi
!------------------------------------------------------
!      Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!      Real(8), parameter :: dpi2=dpi*dpi
!      Real(8), parameter :: d180_pi=180.d0/dpi
!      Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

      mu0=mu01
      muv=muv1

!      sin_eps=dsin(st_min*dpi_180)

      sin_tv=dsqrt(1.d0-muv*muv)
      sin_ti=dsqrt(1.d0-mu0*mu0)
      zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
      cos_tet=-zz
      sin_tet=dsqrt(1.d0-cos_tet*cos_tet)

      If (dabs(sin_tv) .lt. sin_eps) then
         tet0v=dacos(muv)
         tet0v=dabs(tet0v-st_min*dpi_180)
         muv=dcos(tet0v)
         sin_tv=dsqrt(1.d0-muv*muv)
         zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
         cos_tet=-zz
         sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
      end if
 
      If (dabs(sin_ti) .lt. sin_eps) then
         tet0v=dacos(mu0)
         tet0v=dabs(tet0v-st_min*dpi_180)
         mu0=dcos(tet0v)
         sin_ti=dsqrt(1.d0-mu0*mu0)
         zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
         cos_tet=-zz
         sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
      end if

      if (dabs(sin_tet) .lt. sin_eps) then
         tet0v=dacos(mu0)
         tet0v=dabs(tet0v-st_min*dpi_180)
         muv=dcos(tet0v)
         sin_tv=dsqrt(1.d0-muv*muv)
         zz=muv*mu0+sin_ti*sin_tv*dcos(phi*dpi_180)
         cos_tet=-zz
         sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
      end if

      if (dabs(sin_tet) .lt. sin_eps) write(*,*) 'sin_tet',sin_tet
      if (dabs(sin_tv)  .lt. sin_eps) write(*,*) 'sin_tv',sin_tv
      if (dabs(sin_ti)  .lt. sin_eps) write(*,*) 'sin_ti',sin_ti

!XH   phase angle alf = pi - scattering angle
      alf=dacos(zz)
      if (present(alf1)) alf1=alf
      mu_alf=dcos(alf*0.5d0)

      r1=m*m*mu_alf
      r2=cdsqrt(m*m-1.d0+mu_alf*mu_alf)
      r_par1=r1-r2
      r_par2=r1+r2

      r_perp1=mu_alf-r2
      r_perp2=mu_alf+r2

      r_par=r_par1/r_par2
      r_perp=r_perp1/r_perp2

      r_par_sq=dble(r_par*dconjg(r_par))
      r_perp_sq=dble(r_perp*dconjg(r_perp))
      r_par_per1=r_par*dconjg(r_perp)
      r_par_per2=dconjg(r_par_per1)

      F_fr=0.d0
      F_fr(1,1)=(r_par_sq+r_perp_sq)*0.5d0
!XH   different treatments for n_el=1,3,4
      if (n_el .ge. 3) then
         F_fr(2,2)= F_fr(1,1)
         F_fr(2,1)=(r_par_sq-r_perp_sq)*0.5d0
         F_fr(1,2)= F_fr(2,1)
         F_fr(3,3)=dble(r_par_per1+r_par_per2)*0.5d0
!         if (n_el .eq. 4) then
!            F_fr(4,4)=F_fr(3,3)
!            F_fr(3,4)=dble((0,-1.d0)*(r_par_per1-r_par_per2))*0.5d0
!            F_fr(4,3)=-F_fr(3,4)
!         end if
         del_phi=phi*dpi_180-dpi
         cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
         cos_eti= (muv+mu0*cos_tet)/sin_tet/sin_ti !!! checked "-"

         sin_etv=dsin(del_phi)*sin_ti/sin_tet
         sin_eti=dsin(del_phi)*sin_tv/sin_tet

         cos_2etv=cos_etv*cos_etv-sin_etv*sin_etv
         sin_2etv=2.d0*cos_etv*sin_etv

         cos_2eti=cos_eti*cos_eti-sin_eti*sin_eti
         sin_2eti=2.d0*cos_eti*sin_eti

         L_v=0.d0
         L_inc=0.d0

         sin_phi=dsin(phi*dpi_180)
         If (dabs(sin_phi) .lt. sin_eps/10.d0) then
            L_v(1,1)=1.d0
            L_v(2,2)=1.d0
            L_v(2,3)=0
            L_v(3,3)=1.d0
            L_v(3,2)=0

            L_inc(1,1)=1.d0
            L_inc(2,2)=1.d0
            L_inc(2,3)=0
            L_inc(3,3)=1.d0
            L_inc(3,2)=0

!            if (n_el .eq. 4) then
!               L_v(4,4)=1.d0
!               L_inc(4,4)=1.d0
!            end if
         else
            L_v(1,1)=1.d0
            L_v(2,2)=cos_2etv
            L_v(2,3)=sin_2etv
            L_v(3,3)=cos_2etv
            L_v(3,2)=-sin_2etv

            L_inc(1,1)=1.d0
            L_inc(2,2)=cos_2eti
            L_inc(2,3)=sin_2eti
            L_inc(3,3)=cos_2eti
            L_inc(3,2)=-sin_2eti

!            if (n_el .eq. 4) then
!               L_v(4,4)=1.d0
!               L_inc(4,4)=1.d0
!            end if
         end if
         R_fr=Matmul(L_v,Matmul(F_fr,L_inc))
      else
         R_fr=F_fr
      end if ! n_el .ge. 3
!
      return
end Subroutine R_Fren_My


!!! Shadowing function
 Subroutine Shad_func(mu0,muv,sigma2,Shadow)
  Real(8), intent(in):: mu0,muv,sigma2
  Real(8), intent(out):: Shadow
  Real(8) S1,s3,s2,mu_ir,mu_sq,cot,t1,t2,ShadI,ShadR!,pi
!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

   !pi=dacos(-1.d0)
 !!!!!!!!another shadowing
   S1=DSQRT(2D0*sigma2/dpi)
   S3=1D0/(DSQRT(2D0*sigma2))
   S2=S3*S3
   mu_ir=mu0
   mu_sq=mu_ir*mu_ir
   cot=mu_ir/DSQRT(1D0-mu_sq)
   T1=DEXP(-cot*cot*S2)
   !CVF doesn't have DERFC. It was replaced by erfc_func
!   T2=erfc_func(cot*S3)
   T2=DERFC(cot*S3)
!	T2=ERFC(cot*S3)
!	T2=1-ERF(cot*S3)
   ShadI=0.5D0*(S1*T1/cot-T2)
   mu_ir=muv
   mu_sq=mu_ir*mu_ir
   cot=mu_ir/DSQRT(1D0-mu_sq)
   T1=DEXP(-cot*cot*S2)
   !CVF doesn't have DERFC. It was replaced by erfc_func
!   T2=erfc_func(cot*S3)
   T2=DERFC(cot*S3)
!      T2=ERFC(cot*S3)
!	T2=1-ERF(cot*S3)
   ShadR=0.5D0*(S1*T1/cot-T2)
   Shadow=1D0/(1D0+ShadI+ShadR)
 ! Shadow=1D0/(1D0+ShadI)*1D0/(1D0+ShadR)

 End Subroutine


!PL Gaussian quadrature 
SUBROUTINE gauss_l(ndim,ngauss,a_lim,b_lim,x_abscissa,w_gauss)

      Real(8),parameter:: eps= 1.d-13
      Integer m,i,j,ndim,ngauss,m_gauss
      Real(8) x_abscissa(ndim),w_gauss(ndim)
	  Real(8) a_lim,b_lim,x_aver,x_diff
	  Real(8) var1,var2,var3,var4,var_old,var_cur,var_new !,pi
!**********************************************************************
!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

      !pi=4.D0*datan(1.D0)
      m_gauss=(ngauss+1)/2
      x_aver=0.5D0*(a_lim+b_lim)
      x_diff=0.5D0*(b_lim-a_lim)

      DO i=1,m_gauss
        var1= dcos(dpi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
	    var3=1.d0
        do while(var3 .gt. eps) ! tl
		    var_old=0.D0
			var_cur=1.d0
            DO j=1,ngauss-1
	           var_new=(((dble(2*j)-1.d0)*var1*var_cur-(dble(j)-1.d0)*var_old)/dble(j))
   			   var_old=var_cur
			   var_cur=var_new
            ENDDO
			var_new=(((dble(2*j)-1.d0)*var1*var_cur-(dble(j)-1.d0)*var_old)/dble(j))
            var4=ngauss*(var1*var_new-var_cur)/(var1*var1-1.d0)
            var2=var1 
            var1= var2-var_new/var4
            var3=dabs(var1-var2)
        enddo ! while(var3 .gt. eps)

        x_abscissa(i)= x_aver-x_diff*var1
        x_abscissa(ngauss+1-i)= x_aver+x_diff*var1
        w_gauss(i)=2.D0*x_diff/((1.D0-var1*var1)*var4*var4)
        w_gauss(ngauss+1-i)= w_gauss(i)
	  end do
!**********************************************************************

 END Subroutine

!PL Error function
 Real(8) Function erfc_func(z)
 Real(8), parameter:: eps=0.000000000001d0
 Real(8) z, erf, erf_old, erf_n, rel_er!, pi
 Integer n
!------------------------------------------------------
!Real(8), parameter :: dpi=3.141592653589793 !23846264338327950
!!Real(8), parameter :: dpi2=dpi*dpi
!Real(8), parameter :: d180_pi=180.d0/dpi
!Real(8), parameter :: dpi_180=dpi/180.d0
!------------------------------------------------------

  If (dabs(z) .le. 1.d-50) then
   erf=0.d0
  else 
   If (dabs(z) .gt. 6) then
    if (z .lt. 0) erf=-1.d0
    if (z .gt. 0) erf=1.d0
   else
    !pi=2.d0*Dasin(1.d0)
    rel_er=1.d0
    erf_n=z
    erf=z
    n=0
    Do While (rel_er .gt. eps)
     n=n+1
     erf_old=erf
     erf_n=erf_n*(1.d0-2.d0*n)*z*z/n/(2.d0*n+1.d0)
     erf=erf+erf_n
     rel_er=dabs((erf-erf_old)/erf)
    End do
	erf=2.d0*erf/dsqrt(dpi)
   end if
  end if
  erfc_func=1.d0-erf
 End Function


 !!! linear interpolation 
 Subroutine Lin_interp(n_tet0,tet0,tet,M_tet,M_el)
  Implicit none
  integer, intent(in):: n_tet0
  Real(8), intent(in):: tet0(1:n_tet0),tet,M_tet(1:n_tet0)
  Real(8), intent(out):: M_el
  integer i_tet
  Real(8) x1,x2,y1,y2,a_c,b_c

    Do i_tet=2, n_tet0
     If (((tet .ge. tet0(i_tet-1)) .and. (tet .le. tet0(i_tet))) .or. &
         ((tet .le. tet0(i_tet-1)) .and. (tet .gt. tet0(i_tet)))) then
      x1=tet0(i_tet-1)
      x2=tet0(i_tet)
      y1=M_tet(i_tet-1)
      y2=M_tet(i_tet)
      a_c=(y2-y1)/(x2-x1)
      b_c=y1-a_c*x1
      M_el=a_c*tet+b_c
	  exit
     end if
    end do !i_tet
   
 End Subroutine Lin_interp
 
End Module Mod_BRM 
