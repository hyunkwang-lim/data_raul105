! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **


!   file contains:
!
!     SUBROUTINE radiative_transfer_SOS
!     subroutine surface1
!     subroutine aerosol1
!     SUBROUTINE phase_matrix_intrpl_NQDR
!     SUBROUTINE MODIF_DE_TRONC
!     SUBROUTINE rad_transf_MS
!     subroutine gauss
!     subroutine developpe_ocean_land
!     subroutine noyaux
!     subroutine fresnel
!     subroutine betal
!     subroutine legendre
!     subroutine discr_profile
!     subroutine root
!     subroutine profiles_new
!     subroutine get_WD
!     subroutine grid_heights
!     subroutine grid_heights_LUT
!     FUNCTION   ANGTURN
!     FUNCTION   LINEAR_LN


!   mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
!XH   truncation of the phase matrix
!XH      ITRONC=F : calculate strictly without truncation
!XH      ITRONC=T : calculate approximately with truncation

!XH   number of scattering angles of phase matrix NANG <= KANG
!   mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
#include "../../../constants_set/mod_globals.inc"
      module mod_os
      use mod_stop_report
!PL   Module for different BRDF/BPDF
      use Mod_BRM, only: x_vect,n_par_max,         &
                         BRM_Fexp_OSH,BRM,BRM_ocean !BRM_Fexp_opt
!      use mod_molecular_scattering
!------------------------------------------------------
      real, parameter :: rpi=3.141592653589793 !23846264338327950 !=3.141592653
!      real, parameter :: rpi2=rpi*rpi
      real, parameter :: r180_pi=180./rpi
      real, parameter :: rpi_180=rpi/180.
!------------------------------------------------------
      contains
	 
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE radiative_transfer_SOS ( RIN,                                          & ! IN
                                          IW, IMSC, NG, NN, NF,                         &
                                          iBRDF, iBPDF, iBRM_water, land_percent,       &
                                          NQDR, ITRONC, iop, NT1,                       &
                                          tetas, NBV, vis, fiv,                         &
                                          WAVE,                                         &
                                          n_par_land, n_par_land_i, surf_land_par_vect, &
                                          n_par_water, surf_water_par_vect,             &
                                          EXT, SSA,                                     &
                                          KMpar, NANG, ANGL,                            &
                                          HOBS_km, HGR_km, HMAX_atm_km, DISCRVD,        &
                                          laerosol, lsurface,                           &
                                          PF11_I, PF12_I, PF22_I, PF33_I,               &

                                          thd, SLout, SQout, SUout, SLPout,             & ! OUT
                                          UFG, DFG, NLV, HLV, UFX, DFX                  &
                                        )

!	------------------------------------------------------------------------------------
      use mod_retr_settings_derived_type
      use mod_par_OS
      use mod_intrpl_linear
!XH   2*KBF is the maximum length of the array containing surface parameters
      use mod_par_inv,   only : KBF
      use mod_vertical_distr_derived_type

      implicit none
      integer, parameter  ::  KANG=2*NG0T+1
!      integer, parameter  ::  KNG=NG0
!	------------------------------------------------------------------------------------
!   IN :
      type(retr_input_settings),         intent(in) ::  RIN
      integer,                           intent(in) ::  NQDR,KMpar,NANG,ITRONC ! ,NM
      real,dimension(KMpar),             intent(in) ::  ANGL
      real,dimension(KMpar,KSD),         intent(in) ::  PF11_I,PF12_I,PF22_I,PF33_I
      real,dimension(NMM+NMG),           intent(in) ::  EXT,SSA
      integer,                           intent(in) ::  IW,IMSC,NG,NN,NF,iBRDF,iBPDF,iBRM_water
      integer,                           intent(in) ::  NBV
      real,                              intent(in) ::  WAVE
      integer,                           intent(in) ::  iop !,ngas
      integer,                           intent(in) ::  n_par_land,n_par_land_i,n_par_water
      real,                              intent(in) ::  land_percent
      real, dimension(2*KBF),            intent(in) ::  surf_land_par_vect
      real, dimension(KBF),              intent(in) ::  surf_water_par_vect
      real,                              intent(in) ::  tetas
      real,dimension(2*NBVM),            intent(in) ::  vis,fiv
      real,                              intent(in) ::  HOBS_km,HGR_km,HMAX_atm_km
      INTEGER,                           intent(in) ::  NT1
      !REAL, DIMENSION(NMM+NMG),          intent(in) ::  ECH,sigma
      !INTEGER, DIMENSION(NMM+NMG),       intent(in) ::  JP

      !integer,                           intent(in) ::  max_NH0
      !integer,dimension(NMM+NMG),        intent(in) ::  NH0
      !real, dimension(max_NH0,NMM+NMG),  intent(in) ::  H0
      !real, dimension(max_NH0,NMM+NMG),  intent(in) ::  PROF0
      logical,                           intent(in) ::  laerosol,lsurface
      !integer,                           intent(in) ::  mol_prof_type, aer_prof_type !AL
      integer,               optional,   intent(in) ::  NLV
      real,DIMENSION(KNT),   optional,   intent(in) ::  HLV
      type(discret_vertical_distribution), intent(in) :: DISCRVD
!	--------------------------------------------------------------------------------------
!   OUT :
      real,dimension(2*NBVM),            intent(out) ::  thd
      real,dimension(2*NBVM),            intent(out) ::  SLout,SQout,SUout,SLPout
      real,                              intent(out) ::  UFG,DFG
      real,DIMENSION(KNT),   optional,   intent(out) ::  UFX,DFX
!	--------------------------------------------------------------------------------------
!   LOCAL :
      integer                               ::  keyEl, n_el, NM, ngas, NM1
      integer                               ::  NANG1
      REAL,DIMENSION(KANG)                  ::  ANGL1
      REAL,DIMENSION(KANG,KSD)              ::  PF11,PF12,PF22,PF33
!PL   save parameters
      LOGICAL,                        SAVE  ::  LFIRST_F = .true.
      LOGICAL,                        SAVE  ::  LFIRST_D = .true.
      REAL*8,DIMENSION(-NG0:NG0),     SAVE  ::  xmu_F,xg_F,xmu_D,xg_D   						 
      REAL*8,DIMENSION(-NN0:NN0),     SAVE  ::  rmu_F,rg_F,rmu_D,rg_D                                           
      REAL*8,DIMENSION(-NG0T:NG0T),   SAVE  ::  zmu_F,zg_F,zmu_D,zg_D	
      REAL,DIMENSION(-NN0:NN0),       SAVE  ::  thv_F,thv_D
      REAL,DIMENSION(0:2*NG0,NMM),    SAVE  ::  alp,bet,gam,zet
      INTEGER,                        SAVE  ::  NAV,NT
      REAL, DIMENSION(KNT),           SAVE  ::  EXTH
      REAL, DIMENSION(KNT-1,NMM+NMG), SAVE  ::  WD
!PL   Because of reciprocity principle (Rij(j,k,IS)=Rji(k,i,IS)),
!PL   it is possible to reduce the number of elements Rij
!PL   But for anizotropy surfaces this can not be done
!XH   need to double think about solutions to memory problems of SQ1,RQ1,R12,R13,R21,R22,R23,R31,R32,R33
      REAL,DIMENSION(2*NBVM),         SAVE  ::  SL1,RL1,SQ1,RQ1
      REAL,DIMENSION(NN0,NN0,0:NF0),  SAVE  ::  R11,R12,R13,     &
                                                R21,R22,R23,     &
                                                R31,R32,R33
!>>>>>>>
      !REAL,DIMENSION(-1:2*NG0,-NN0:NN0,0:NF0),TARGET, SAVE  ::  PSL_F,RSL_F,TSL_F
      !REAL,DIMENSION(-1:2*NG0,-NN0:NN0,0:NF0),TARGET, SAVE  ::  PSL_D,RSL_D,TSL_D
!>>>>>>> Tested. Works good

      REAL*8,DIMENSION(-NG0T:NG0T)          ::  zmu,zg
      REAL*8,DIMENSION(-NG0:NG0)            ::  xmu,xg
      REAL*8,DIMENSION(-NN0:NN0)            ::  rmu,rg
      REAL,  DIMENSION(-NN0:NN0)            ::  thv
      INTEGER,DIMENSION(2*NBVM)             ::  idir
      REAL,   DIMENSION(2*NBVM)             ::  coff,chi
      REAL,   DIMENSION(2*NBVM)             ::  cos_2etv,sin_2etv
      REAL,   DIMENSION(2*NBVM)             ::  AT
      real,dimension(-NG0:NG0,NMM)          ::  pha11,pha12,pha22,pha33
!      REAL, POINTER                         ::  PTR_PSL(:,:,:),  &
!                                                PTR_RSL(:,:,:),  &
!                                                PTR_TSL(:,:,:)
      integer                               ::  i,iv,ISD,j,IS,m,k
      integer                               ::  itest
      real                                  ::  aaa,alpm,betm,gamm, &
                                                rmus,ron,tang,      &
                                                xx,yy,zz
!	------------------------------------------------------------------------------------------------------
      real                                  ::  tetot_os,tevtot
      real,dimension(NMM+NMG)               ::  EXT_os,SSA_os
      REAL,DIMENSION(2*NBVM)                ::  DDL1, SL
      REAL,DIMENSION(2*NBVM)                ::  DDQ1, SQos,SUos,SLP,SQ,SU,SLPsig
      real                                  ::  cc,ss,DOMEGA
!XH   list of levels used internally in the radiative transfer calculation
      real,DIMENSION(KNT)                   ::  HLV0,UFX0,DFX0
      real                                  ::  FAC
      real(8)                               ::  x_vect_water(1:n_par_water), x_vect_land(1:n_par_land)
      logical                               ::  IPRI_additional_info
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
      NM    = DISCRVD%naer
      ngas  = DISCRVD%ngas
      keyEl = RIN%DLSF%keyEl
      IPRI_additional_info=RIN%IPRI_additional_info


      SLout( :) = 0.0
      SQout( :) = 0.0
      SUout( :) = 0.0
      SLPout(:) = 0.0
      thd(:)   = 0.0
      UFG      = 0.0
      DFG      = 0.0
      if (present(UFX)) then
         UFX0=0.0
         DFX0=0.0
      end if

!     n_par_land   - total number of surface parameters (BRDF+BPDF) for land
!     n_par_land_i - number of BRDF surface parameters for land
!     n_par_water  - number of parameters for water

      x_vect_land(1:n_par_land)=surf_land_par_vect(1:n_par_land)
      x_vect_water(1:n_par_water)=surf_water_par_vect(1:n_par_water)
	   
      if(NANG .gt. KANG) then
        write(tmp_message,'(2(a,i0))') 'NANG = ',NANG,' .gt. KANG = ',KANG
        G_ERROR(trim(tmp_message))
      endif

!PL   Calculating Gaussian weights (zg) and quadratures (zmu)
!PL   for phase matrix expansion (IMSC =0, 1 - forward calculations)
!PL   (IMSC=2 -Jacobean matrix calculations).
!PL   Gaussian quadratures for phase matrix expansion "zmu"
!PL   are calculated in the range -NQDR:NQDR
!PL   Therefore NG0>=NQDR !!!

	
!PL   Calculating Gaussian weights and quadratures
!PL   "xg" and "xmu" are used  for expansion of trancated Phase matrix
!PL   "rmu" and "rg" are used for intagrals calculations in RT.
!PL   (IMSC =0, 1 - forward calculations)
!PL   (IMSC=2 -Jacobean matrix calculations).
!PL   Array xmu,xg can be defined in the range  -NG:NG !!!
!PL   Array rmu,rg can be defined in the range  -NN:NN !!!
!PL   Therefore NN0>=NN

      select case (IMSC)      
!      case(0,1) ! observation modeling (multiple scattering) 
      case(0)
      IF (LFIRST_F) THEN
            zmu_F(-NG0T:NG0T) = 0.0
            zg_F (-NG0T:NG0T) = 0.0
            call GAUSS(NQDR,zmu_F,zg_F,NG0T)
            xmu_F(-NG0:NG0)=0.0
            xg_F (-NG0:NG0)=0.0
            call gauss(NG,xmu_F,xg_F,NG0)
            rmu_F(-NN0:NN0)=0.0
            rg_F (-NN0:NN0)=0.0			
            call gauss(NN,rmu_F,rg_F,NN0)
            thv_F(-NN0:NN0)=0.0			
            thv_F(-NN:NN)=acos(rmu_F(-NN:NN))*r180_pi
!            IF (IMSC.EQ.0) THEN
!               PSL_F(:,:,:)=0.0
!               TSL_F(:,:,:)=0.0
!               RSL_F(:,:,:)=0.0
!               Do IS=0,NF
!                  call legendre(IS,rmu_F,NG,NN,PSL_F(:,:,IS),TSL_F(:,:,IS),RSL_F(:,:,IS))
!               End do ! IS
!            END IF
            LFIRST_F = .false.
         END IF ! LFIRST_F
         xmu(-NG0:NG0) = xmu_F(-NG0:NG0)
         xg (-NG0:NG0) = xg_F (-NG0:NG0)
         rmu(-NN0:NN0) = rmu_F(-NN0:NN0)
         rg (-NN0:NN0) = rg_F (-NN0:NN0)
         thv(-NN0:NN0) = thv_F(-NN0:NN0)
         zmu(-NG0T:NG0T) = zmu_F(-NG0T:NG0T)
         zg (-NG0T:NG0T) = zg_F (-NG0T:NG0T)
!         IF (IMSC.EQ.0) THEN
!            PTR_PSL(-1:,-NN0:,0:) => PSL_F
!            PTR_TSL(-1:,-NN0:,0:) => TSL_F
!            PTR_RSL(-1:,-NN0:,0:) => RSL_F
!         END IF
      case(2) ! derivatives (multiple scattering)
         IF (LFIRST_D) THEN
            zmu_D(-NG0T:NG0T) = 0.0
            zg_D (-NG0T:NG0T) = 0.0
            call GAUSS(NQDR,zmu_D,zg_D,NG0T)
            xmu_D(-NG0:NG0)=0.0
            xg_D (-NG0:NG0)=0.0
            call gauss(NG,xmu_D,xg_D,NG0)
            rmu_D(-NN0:NN0)=0.0
            rg_D (-NN0:NN0)=0.0
            call gauss(NN,rmu_D,rg_D,NN0)		
            thv_D(-NN0:NN0)=0.0			
            thv_D(-NN:NN)=acos(rmu_D(-NN:NN))*r180_pi
!            PSL_D(:,:,:)=0.0
!            TSL_D(:,:,:)=0.0
!            RSL_D(:,:,:)=0.0
!            Do IS=0,NF
!               call legendre(IS,rmu_D,NG,NN,PSL_D(:,:,IS),TSL_D(:,:,IS),RSL_D(:,:,IS))
!            End do ! IS
            LFIRST_D = .false.
         END IF ! LFIRST_D
         xmu(-NG0:NG0) = xmu_D(-NG0:NG0)
         xg (-NG0:NG0) = xg_D (-NG0:NG0)
         rmu(-NN0:NN0) = rmu_D(-NN0:NN0)
         rg (-NN0:NN0) = rg_D (-NN0:NN0)
         thv(-NN0:NN0) = thv_D(-NN0:NN0)
         zmu(-NG0T:NG0T) = zmu_D(-NG0T:NG0T)
         zg (-NG0T:NG0T) = zg_D (-NG0T:NG0T)
!         PTR_PSL(-1:,-NN0:,0:) => PSL_D
!         PTR_TSL(-1:,-NN0:,0:) => TSL_D
!         PTR_RSL(-1:,-NN0:,0:) => RSL_D
      end select

!      if(IMSC .ne. 1) then
!      WRITE(*,*) 'xmu:'
!      WRITE(*,'(10e14.4)') xmu(-NG0:NG0)
!      WRITE(*,*) 'xg:'
!      WRITE(*,'(10e14.4)') xg(-NG0:NG0)
!      WRITE(*,*) 'rmu:'
!      WRITE(*,'(10e14.4)') rmu(-NN0:NN0)
!      WRITE(*,*) 'rg:'
!      WRITE(*,'(10e14.4)') rg(-NN0:NN0)
!      WRITE(*,*) 'IMSC=',IMSC,'  thv:'
!      WRITE(*,'(10e14.4)') thv(-NN:NN)

!      WRITE(*,'(10e14.4)') thv(-NN0:NN0)
!      WRITE(*,*) 'zmu:'
!      WRITE(*,'(10e14.4)') zmu(-NG0T:NG0T)
!      WRITE(*,*) 'zg:'
!      WRITE(*,'(10e14.4)') zg(-NG0T:NG0T)
!      endif
!	------------------------------------------------------------------------------------------------------

!XH   cosine of solar zenith angle - rmus
      rmus=cos(tetas*rpi_180)
!      do j=-NN,NN
!         thv(j)=acos(rmu(j))*r180_pi
!         WRITE(*,*) thv(j),j
!      end do ! j

!XH   observation geometries
!XH   itest = +1 : Satellite observations
!XH   itest = -1 : Ground based observations

      itest  =+1  ! Satellite observations
      idir   = 0
      coff   = 0.0

!      write(*,*) 'in OS_H NBV=',NBV
      do iv=1,NBV
         chi(iv)=ANGTURN(tetas,vis(iv),fiv(iv))
!PL       cos_2etv=-cc !!! Because of positive SQ1,SL1
!PL       sin_2etv=-ss !!! The reason is not clear yet
!PL       call ANGTURN_PL(tetas,vis(iv),fiv(iv),cos_2etv(iv),sin_2etv(iv))
         xx=cos(vis(iv)*rpi_180)
         if (itest .eq. +1 .and. xx .lt. 0.0) itest=-1  ! Ground based observations
!PL       zz=-rmus*xx-sqrt(1.0-rmus*rmus)*sqrt(1.0-xx*xx)*cos(fiv(iv)*rpi_180)
         zz=-rmus*xx-real(dsqrt(1.d0-rmus*rmus)*dsqrt(1.d0-xx*xx))*cos(fiv(iv)*rpi_180)
         if (zz .GT. 1.0) then
            write(*,*) 'Incorrect cosine of scattering angle',zz
            zz = 1.0
         else if (zz .LT. -1.0) then
            write(*,*) 'Incorrect cosine of scattering angle',zz
            zz =-1.0
         end if ! zz .GT. 1.0

         thd(iv)=acos(zz)*r180_pi
         if (IMSC .ne. 1) then
         if (thv(NN) .lt. vis(iv)) then
            if (abs(thv(-NN)-vis(iv)) .le. 1e-4) then !pl&tl 23-02-2015
               k = -NN
               idir(iv) = -NN
               coff(iv) = 0.0
            else
               k = NN-1
               do while (thv(k) .lt. vis(iv))
                  k = k-1
               end do
               idir(iv) = k
               coff(iv) = (vis(iv)-thv(k+1))/(thv(k)-thv(k+1))
            endif
         else
!XH         for the exact nadir direction       
            !k = NN-1
            idir(iv) = NN-1
            coff(iv) = 0.0
         end if
!         write(*,*) 'vis:'
!         write(*,'(10f16.4)') vis(1:NBV)
!         write(*,*) 'in OS_H k=',k,' NN=',NN,'  thv(k)=',thv(k),'  thv(k+1)=',thv(k+1),'  iv=',iv,'  vis(iv)=',vis(iv)
         endif ! IMSC .ne. 1
      end do !iv
!	------------------------------------------------------------------------------------------------------

!PL   Expansion coeff. for Rayleigh scattering
!PL   ron: molecular depolarization factor
!XH   gamm=alpm=0.0 for scalar case
      ron =0.0279  ! ron=0.014  ron=0.03  ron=0
      aaa =2.0*(1.0-ron)/(2.0+ron)
      betm=0.5*aaa
      if (keyEl .eq. 4) then
         gamm=-sqrt(1.5)*aaa
         alpm=       3.0*aaa
      else
         gamm=0.0
         alpm=0.0
      end if

!     Phase matrix interpolation
!PL   Subroutine "phase_matrix_intrpl_NQDR" interpolates Phase Matrix "PF11_I"
!PL   defined at angles "ANGL" into the Phase Matrix  "PFij" which is defined at
!PL   angles "ANGL1"=acos(zmu)

      If (laerosol) then
!         write(*,*) NM,NQDR,KMpar,NANG,'  - NM,NQDR,KMpar,NANG'
!         do i=1,NANG
!            write(*,'(f14.3,10e14.5)') ANGL(i),PF11_I(i,1:NM)
!         end do
!         write(*,*) 'zmu:'
!         write(*,'(10e14.5)') zmu(1:NQDR)
!         write(*,*) 'zg:'
!         write(*,'(10e14.5)') zg(1:NQDR)
         select case(IMSC)
         case (0,2)
            call phase_matrix_intrpl_NQDR (                                 &
                                             keyEl,                         & ! IN
                                             NM,NQDR,                       &
                                             zmu,zg,                        &
                                             KMpar,NANG,ANGL,               &
                                             PF11_I,PF12_I,PF22_I,PF33_I,   &
                                             NANG1,ANGL1,                   & ! OUT
                                             PF11,PF12,PF22,PF33            &
                                          )
            if ( stop_report%status ) return
         case (1) ! single scattering
            if (KANG .lt. KMpar) then
               write(tmp_message,'(2(a,i0))') 'KANG = ',KANG,' lt. KMpar = ',KMpar
               G_ERROR(trim(tmp_message))
            endif
            NANG1 = NANG
            ANGL1(1:NANG) = ANGL(1:NANG)
            PF11(1:NANG,1:NM) = PF11_I(1:NANG,1:NM)
            if (keyEl .eq. 4) then
               PF12(1:NANG,1:NM) = PF12_I(1:NANG,1:NM)
               PF22(1:NANG,1:NM) = PF22_I(1:NANG,1:NM)
               PF33(1:NANG,1:NM) = PF33_I(1:NANG,1:NM)
            endif ! keyEl .eq. 4
         end select ! IMSC
        
         NM1=NM+1
         EXT_os(1:NM1+ngas)=EXT(1:NM1+ngas) ! aerosol, molecular, gas 
         SSA_os(1:NM1+ngas)=SSA(1:NM1+ngas)
         tetot_os=SUM(EXT_os(1:NM1+ngas))
         tevtot=SUM(EXT(1:NM1+ngas))

!PL      Subroutine "MODIF_DE_TRONC" substitutes the Phase Matrix  "PFij" with the truncated
!PL      matrix, and modifies the AOD (tetot_os,EXT_os(1:NM))
!PL      and single scattering albedo "SSA_os" according to the truncation

!         write(*,*) ITRONC,NM,NQDR,itest
!         write(*,*) 'before MODIF_DE_TRONC: tetot_os',tetot_os
!         write(*,*) 'before MODIF_DE_TRONC: EXT_os(1:NM)',EXT_os(1:NM)
!         write(*,*) 'before MODIF_DE_TRONC: SSA_os(1:NM)',SSA_os(1:NM)

!         do i=1,NANG1
!            write(*,'(f14.3,10e14.5)') ANGL1(i),PF11(i,1:NM)
!         end do
         if (ITRONC .EQ. 1) then
            CALL MODIF_DE_TRONC_new (                                    &
                                     keyEl,NM,NQDR,itest,                & ! IN
                                     zmu,zg,NANG1,ANGL1,                 &
                                     NBV,thd,rmus,                       &
                                     tetot_os,EXT_os(1:NM),SSA_os(1:NM), & ! INOUT
                                     PF11,PF12,PF22,PF33,                &
                                     DDL1,DDQ1                           & ! OUT
                                    )
         end if

!PL      Subroutine "profiles_new" calculates aerosol-Reyleigh altitude
!PL      AOD profile "WD"
         if (present(UFX)) then
            call profiles_new (  IPRI_additional_info,     &
                                 tetot_os, EXT_os, SSA_os, &
                                 HGR_km, HOBS_km, HMAX_atm_km, NT1, &
                                 DISCRVD,   &
                                 NT, NAV, EXTH, WD, HLV0   & ! OUT
                              )
            if ( stop_report%status ) return
         else
            call profiles_new (  IPRI_additional_info,     &
                                 tetot_os, EXT_os, SSA_os, &
                                 HGR_km, HOBS_km, HMAX_atm_km, NT1, &
                                 DISCRVD,   &
                                 NT, NAV, EXTH, WD         & ! OUT
                              )
            if ( stop_report%status ) return
         end if

!PL      Reyleigh + Aeroslol phase matrix for incident and viewing geometries
!         write(*,*) NM,NT,NBV,NAV,'  - NM,NT,NBV,NAV'
!         write(*,*) 'ANGL1:'
!         do i=1,NANG1
!            write(*,'(f14.3,10e14.5)') ANGL1(i),PF11(i,1:NM),PF12(i,1:NM)
!         end do
!XH      the first order scattering is given by SL1 and SQ1
         call aerosol1 (keyEl,itest,     & ! IN
                        NM,NT,NBV,NAV,   &
                        NANG1,ANGL1,     &
                        rmus,thd,vis,    &
                        EXTH,WD,         &
                        betm,gamm,       &
                        PF11,PF12,       &

                        SL1,SQ1          & ! OUT
                       )

!PL      Diffuse caclulation (begin)
!        Interpolation of phase matrix at the required angles:
         IF (IMSC .NE. 1) THEN
!PL         PF11 is defined in the range (1:2*NQDR+1) (that is (-NQDR,NQDR))
!PL         pha11 is defined in the range (-NG,NG)
            IF (NQDR .EQ. NG) THEN
               pha11(-NG:NG,1:NM)=PF11(1:NG+NG+1,1:NM)
               if (keyEl .eq. 4) then
!XH               considering linear polarization
                  pha12(-NG:NG,1:NM)=PF12(1:NG+NG+1,1:NM)
                  pha22(-NG:NG,1:NM)=PF22(1:NG+NG+1,1:NM)
                  pha33(-NG:NG,1:NM)=PF33(1:NG+NG+1,1:NM)
               end if
            ELSE
               DO j=-NG,NG
                  tang=acos(xmu(j))*r180_pi
                  DO ISD=1,NM
                     pha11(j,ISD)=LINEAR_LN(ANGL1(1:NANG1),PF11(1:NANG1,ISD),NANG1,tang)
                  END DO
                  if (keyEl .eq. 4) then
!XH                  considering linear polarization
                     DO ISD=1,NM
                        pha12(j,ISD)=LINEAR(ANGL1(1:NANG1),PF12(1:NANG1,ISD),NANG1,tang)
                        pha22(j,ISD)=LINEAR_LN(ANGL1(1:NANG1),PF22(1:NANG1,ISD),NANG1,tang)
                        pha33(j,ISD)=LINEAR(ANGL1(1:NANG1),PF33(1:NANG1,ISD),NANG1,tang)
                     END DO
                  end if
               END DO
            END IF ! NQDR.EQ.NG

            do m=1,NM
               call betal(keyEl,NG,NG0,xmu,xg,pha11(:,m),pha12(:,m),pha22(:,m),pha33(:,m),  &
                                                    alp(0:NG0+NG0-2,m),bet(0:NG0+NG0-2,m),  &
                                                    gam(0:NG0+NG0-2,m),zet(0:NG0+NG0-2,m)  )
            end do ! m
         END IF  !PL IMSC .NE. 1
	
      end if  !PL If (laerosol)

!	  ------------------------------------------------------------------------------------------------------

!      isaut=int(land_percent/100)+1
!	  select case(isaut)
!PL   Ocean
!      case(1)
!         xind=Rem_water 
!PL   Land Surfaces with different BRDFs
!      case(2)
!         xind=Rem_land 
!      case default
!         write(*,*) 'isaut=',isaut,'  - unknown value'
!         stop 'stop in mod_os: radiative_transfer_SOS'
!	  end select ! isaut
	   
!      write(*,*) 'in radiative_transfer_SOS: NF,iBRDF,iBPDF,iBRM_water,isaut,IW : ',NF,iBRDF,iBPDF,iBRM_water,isaut,IW

! ------------------------------------------------------------------------------------------------------

!PL   surface (begin)
      if (lsurface) then
!PL      Surface reflection for incident and viewing geometries
         RL1(:)=0.0
         if (keyEl .eq. 1) then
!XH         scalar case
            n_el=1
         else if (keyEl .eq. 4) then
!XH         consider linear polarization
            RQ1(:)=0.0
            n_el=3
         end if
         if (itest .eq. +1) then   ! ask OD and PL
!XH         surface only constributes to upward radiance at TOA for the first order of scattering
            call surface1 (n_el,NBV,iBRDF,iBPDF,iBRM_water,land_percent,WAVE,              & ! IN
                           rmus,vis,fiv,n_par_land,n_par_land_i,x_vect_land(1:n_par_land), &
                           n_par_water,x_vect_water(1:n_par_water),                        &
                           RL1,RQ1                                                         & ! OUT
                          )
         end if ! itest

         if (IMSC .ne. 1) then
!PL            Subroutine "Developpe_ocean_land" calculates Fourier expansion coefficients for
!PL            surfaces
!PL            Because of reciprocity principle (Rij(j,k,IS)=Rji(k,i,IS)),
!PL            potentially, it is possible to reduce aamount of variables
!PL            R11,R21,R31,R22,R33,R32
!PL            Subroutine for surface expansion calculation (Different BRDF/BPDF)

               call developpe_ocean_land (n_el,land_percent,NN,NF,rmu,WAVE,                    &
                                          iBRDF,iBPDF,n_par_land,n_par_land_i,x_vect_land,     &
                                          iBRM_water,n_par_water,x_vect_water(1:n_par_water),  &
                                          R11,R12,R13,R21,R22,R23,R31,R32,R33)
           end if !IMSC
      end if !PL lsurface

!     TEST
!      write(*,*) ' IN OS_H  IMSC=',IMSC,' laerosol=',laerosol,' lsurface=',lsurface
!      write(*,*)
!      write(*,*) 'SL1:'
!      write(*,'(10e16.6)') SL1
!      write(*,*) 'SQ1:'
!      write(*,'(10e16.6)') SQ1
!      write(*,*)

!     TEST
!      write(*,*) 'RL1:'
!      write(*,'(10e16.6)') RL1
!      write(*,*) 'RQ1:'
!      write(*,'(10e16.6)') RQ1
!      write(*,*)

!     TEST
!      write(*,*) 'R11:'
!      write(*,'(10e16.6)') R11
!      write(*,*) 'R12:'
!      write(*,'(10e16.6)') R12
!      write(*,*) 'R13:'
!      write(*,'(10e16.6)') R13

!      write(*,*) 'R21:'
!      write(*,'(10e16.6)') R21
!      write(*,*) 'R22:'
!      write(*,'(10e16.6)') R22
!      write(*,*) 'R23:'
!      write(*,'(10e16.6)') R23

!      write(*,*) 'R31:'
!      write(*,'(10e16.6)') R31
!      write(*,*) 'R32:'
!      write(*,'(10e16.6)') R32
!      write(*,*) 'R22:'
!      write(*,'(10e16.6)') R33
!      write(*,*)


!      open(100,file=trim(RIN%DLSF%external_file_path)//'test_OSH.txt')
																																			
!      write(100,*) 'R11:'
!      write(100,'(10e16.6)') R11
!      write(100,*) 'R12:'
!      write(100,'(10e16.6)') R12
!      write(100,*) 'R13:'
!      write(100,'(10e16.6)') R13
									
!      write(100,*) 'R21:'
!      write(100,'(10e16.6)') R21
!      write(100,*) 'R22:'
!      write(100,'(10e16.6)') R22
!      write(100,*) 'R23:'
!      write(100,'(10e16.6)') R23
							        
!      write(100,*) 'R31:'
!      write(100,'(10e16.6)') R31
!      write(100,*) 'R32:'
!      write(100,'(10e16.6)') R32
!      write(100,*) 'R22:'
!      write(100,'(10e16.6)') R33
!      write(100,*)
									
!      close(100)

! ------------------------------------------------------------------------------------------------------
!PL   Single scattering + surface (end)
! ------------------------------------------------------------------------------------------------------

      SL(:) = 0.0
      if (keyEl .eq. 4) then
         SQ(:)  =0.0
         SU(:)  =0.0
         SQos(:)=0.0
         SUos(:)=0.0
      end if
	  
      IF (IMSC .NE. 1) THEN
         if (present(UFX)) then
!XH         calculate flux profile
            CALL rad_transf_MS (                                           &
                                 keyEl,NN,NG,NF,NBV,NM,                    & ! IN
                                 itest,NAV,NT,EXTH,WD,                     &
                                 fiv,rmu,rg,tetas,thv,                     &
                                 idir,coff,RIN%eps_err,                    &
                                 alpm,betm,gamm,                           &
                                 alp,bet,gam,zet,                          &
!                                 PTR_PSL,PTR_TSL,PTR_RSL,                  &
                                 R11,R12,R13,R21,R22,R23,R31,R32,R33,      &

                                 SL,SQos,SUos,                             & ! OUT
                                 UFG,DFG,UFX0,DFX0                         &
                               )
            if ( stop_report%status ) return
         else
!XH            do not calculate flux profile
               CALL rad_transf_MS (                                           &
                                    keyEl,NN,NG,NF,NBV,NM,                    & ! IN
                                    itest,NAV,NT,EXTH,WD,                     &
                                    fiv,rmu,rg,tetas,thv,                     &
                                    idir,coff,RIN%eps_err,                    &
                                    alpm,betm,gamm,                           &
                                    alp,bet,gam,zet,                          &
!                                    PTR_PSL,PTR_TSL,PTR_RSL,                  &
                                    R11,R12,R13,R21,R22,R23,R31,R32,R33,      &

                                    SL,SQos,SUos,                             & ! OUT
                                    UFG,DFG                                   &
                                  )
            if ( stop_report%status ) return
         end if
!            write(*,*) 'SL:'
!            write(*,'(10e16.6)') SL
!            write(*,*) 'SQos:'
!            write(*,'(10e16.6)') SQos
!            write(*,*) 'SUos:'
!            write(*,'(10e16.6)') SUos

!            NULLIFY (PTR_PSL,PTR_TSL,PTR_RSL)

         if (keyEl .eq. 4) then
!XH         rotation to the scattering plane

            do iv=1,NBV
!              Original
		       cc=cos(2.0*chi(iv))
!              Original
		       ss=sin(2.0*chi(iv))
!PL            Stokes vector in the scattering plane (SL, SQ)
!              Original
 		       SQ(iv)= cc*SQos(iv)+ss*SUos(iv)
!PL            U Stokes parameter in scattering plane
!              Original
		       SU(iv)=-ss*SQos(iv)+cc*SUos(iv)
!PL
!               cc=-cos(2.*chi(iv))
!               ss=-sin(2.*chi(iv))
!PL		
!               cc=cos_2etv(iv)
!PL		
!               ss=sin_2etv(iv)
!
!               SQ(iv)=-(cc*SQos(iv)-ss*SUos(iv)) !sign is changed
!               SU(iv)=ss*SQos(iv)+cc*SUos(iv)
            end do ! iv

         end if ! keyEl .eq. 4

         if (ITRONC .eq. 1 .and. itest .eq. -1) then
            SL(1:NBV)=SL(1:NBV)+DDL1(1:NBV)
            if (keyEl .eq. 4) SQ(1:NBV)=SQ(1:NBV)+DDQ1(1:NBV)
         end if

      END IF ! IMSC.NE.1

! ------------------------------------------------------------------------------------------------------

      AT(:)=0.0 
      if (itest .eq. +1) then ! Satellite ! ask OD or PL
         do j=1,NBV
            xx=cos(vis(j)*rpi_180)
!            if(xx. lt. 0.0) cycle ! do not delete
            AT(j)=exp(-EXTH(NT)/rmus)*exp(-(EXTH(NT)-EXTH(NAV))/xx)
         end do
      end if ! itest .eq. +1

      if (keyEl .eq. 4) then
!XH      considering lineaer polarization
!PL      Stokes vector in the scattering plane (SL, SQ)
         SL(1:NBV)=SL(1:NBV)+SL1(1:NBV)+RL1(1:NBV)*AT(1:NBV)
         SQ(1:NBV)=SQ(1:NBV)+SQ1(1:NBV)+RQ1(1:NBV)*AT(1:NBV)
	     do iv=1,NBV
!           Original
		    cc=cos(2.0*chi(iv))
!           Original
		    ss=sin(2.0*chi(iv))
!PL         Calculation of Q and U Stokes parameters in meridional plane
!           Original
		    SQos(iv)=cc*SQ(iv)-ss*SU(iv)
!           Original
		    SUos(iv)=ss*SQ(iv)+cc*SU(iv)

!PL         del
!		     SQos(iv)=cc*SQ(iv)-ss*SU(iv)
!		     SUos(iv)=ss*SQ(iv)+cc*SU(iv)
!PL	   
!            cc=-cos(2.*chi(iv))
!            ss=-sin(2.*chi(iv))
!PL    
!            cc=cos_2etv(iv)
!PL	   
!            ss=sin_2etv(iv)
!PL         Calculation of Q and U Stokes parameters in meridional plane
!            SQos(iv)=-cc*SQ(iv)+ss*SU(iv)
!            SUos(iv)=ss*SQ(iv)+cc*SU(iv)

            xx=SQ(iv)
            yy=SU(iv)
            SLP(iv)=sqrt(xx*xx+yy*yy)
            SLPsig(iv)=SLP(iv)
            if (SQ(iv) .lt. 0.0) then 
               SLPsig(iv)=-SLP(iv)
            endif
         end do ! iv

         select case(iop)
         case(0)
!           out is scattering plance
            do iv=1,NBV
               SQout(iv)=SQ(iv)
               SUout(iv)=SU(iv)
               SLPout(iv)=SLPsig(iv)
            end do
         case(1)
!           out is meridian plance
            do iv=1,NBV
               SQout(iv)=SQos(iv)
               SUout(iv)=SUos(iv)
               SLPout(iv)=SLPsig(iv)
            end do
         case default
            write(tmp_message,'(a,i0,a)') 'iop = ',iop,' - unknown value'
            G_ERROR(trim(tmp_message))
         end select
      else
!XH      for scalar case
         SL(1:NBV)=SL(1:NBV)+SL1(1:NBV)+RL1(1:NBV)*AT(1:NBV)
      end if

      do iv=1,NBV
         SLout(iv)=SL(iv)
!         IF (thd(iv).EQ.0) THEN
!!PL      Extinction in  the exact forward direction
!         if (ABS(thd(iv)) .lt. 0.05) then
!            DOMEGA=-2.0*rpi*(COS(0.6*rpi_180)-COS(0.0))
!!PL         This is an absolute flux  for F0=pi
!            SLout(iv)=pi*exp(-tevtot/rmus)+DOMEGA*SL(iv)
!!PL         This is TAU measured by sunphotometer
!            write(*,*) tevtot,'tevtot'
!!PL          SLout(iv)=-log(SLout(iv)/pi)*rmus
!         end if
      end do ! iv

!XH   interpolate flux to the user-defined levels
      if (present(UFX)) then
         do iv=1,NLV
            if (HLV(iv) .GE. HLV0(1)) then
               UFX(iv)=UFX0(1)
               DFX(iv)=DFX0(1)
            else if (HLV(iv) .LE. HLV0(NT)) then
               UFX(iv)=UFX0(NT)
               DFX(iv)=DFX0(NT)
            else
               do j=2,NT
                  if (HLV(iv) .EQ. HLV0(j-1)) then
                     UFX(iv)=UFX0(j-1)
                     DFX(iv)=DFX0(j-1)
                     exit
                  else if (HLV(iv) .EQ. HLV0(j)) then
                     UFX(iv)=UFX0(j)
                     DFX(iv)=DFX0(j)
                     exit
                  else if (HLV(iv) .LT. HLV0(j-1) .AND. HLV(iv) .GT. HLV0(j)) then
!XH                  interpolation in linear scale
                     FAC = (HLV(iv)-HLV0(j))/(HLV0(j-1)-HLV0(j))
                     UFX(iv)=UFX0(j)+FAC*(UFX0(j-1)-UFX0(j))
                     DFX(iv)=DFX0(j)+FAC*(DFX0(j-1)-DFX0(j))
!XH                  interpolation in logrithmic scale
!                     FAC = LOG((HLV(iv)+1.0)/(HLV0(j)+1.0))/LOG((HLV0(j-1)+1.0)/(HLV0(j)+1.0))
!                     UFX(iv)=UFX0(j)+FAC*(UFX0(j-1)-UFX0(j))
!                     DFX(iv)=DFX0(j)+FAC*(DFX0(j-1)-DFX0(j))
                     exit
                  end if
               end do
            end if
         end do
      end if
!
      
      RETURN
      END SUBROUTINE radiative_transfer_SOS

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
  
      subroutine surface1 (n_el,NBV,iBRDF,iBPDF,iBRM_water,land_percent,WAVE,        & ! IN
                           rmus,vis,fiv,n_par,n_par_i,x_vect_land,                   &
                           n_par_water,x_vect_water,                                 &
                           RL1,RQ1                                                   & ! OUT
                          ) 
!XH   optional parameters for vector case
      use mod_par_OS, only : NBVM
      implicit none

!	------------------------------------------------------------------------------------	
!   IN:
!XH   add variable n_el to determine dimension of Rip for scalar and vector cases
      integer,                          intent(in)    ::  NBV,iBRDF,iBPDF,iBRM_water,n_el
      real,dimension(2*NBVM),           intent(in)    ::  vis,fiv
      real,                             intent(in)    ::  land_percent,rmus,WAVE
      integer,                          intent(in)    ::  n_par,n_par_i,n_par_water
      Real(8),                          intent(in)    ::  x_vect_land(1:n_par),x_vect_water(1:n_par_water)
!	------------------------------------------------------------------------------------
!   OUT:
      real,dimension(2*NBVM),           intent(out)   ::  RL1,RQ1
!	------------------------------------------------------------------------------------
!   LOCAL:
      Real(8)  ::  Rip_land(1:n_el,1:n_el),Rip(1:n_el,1:n_el)
      real     ::  xxx!, AT
      integer  ::  j
!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------

      do j=1,NBV
         xxx=cos(vis(j)*rpi_180)
!         if (xxx .ge. 0.0) then ! do not delete
!            AT=exp(-EXTH(NT)/rmus)*exp(-(EXTH(NT)-EXTH(NAV))/xxx)
         select case(nint(land_percent))
         case(0)
 	        call BRM_ocean(n_el,iBRM_water,n_par_water,dble(rmus),dble(xxx),dble(fiv(j)),dble(WAVE), &
                         x_vect_water(1:n_par_water),Rip)
         case(100)
	        call BRM(n_el,iBRDF,iBPDF,n_par,n_par_i,dble(rmus),dble(xxx),dble(fiv(j)),dble(WAVE),    &
                         x_vect_land(1:n_par),Rip)
         case default
          call BRM_ocean(n_el,iBRM_water,n_par_water,dble(rmus),dble(xxx),dble(fiv(j)),dble(WAVE),   &
                         x_vect_water(1:n_par_water),Rip)
          call BRM(n_el,iBRDF,iBPDF,n_par,n_par_i,dble(rmus),dble(xxx),dble(fiv(j)),dble(WAVE),      &
                         x_vect_land(1:n_par),Rip_land)
          Rip=(1. - 0.01*land_percent)*Rip+0.01*land_percent*Rip_land
!            write(*,*) 'land_percent=',land_percent,'  - unknown'
!            stop 'stop in mod_os'
         end select ! land_percent
!         RQ1(j)=AT*rmus*dsqrt(Rip(2,1)*Rip(2,1)+Rip(3,1)*Rip(3,1))
!	      RL1(j)=AT*rmus*Rip(1,1)
         RL1(j)=rmus*Rip(1,1)
!XH      considering linear polarization
         if (n_el .eq. 3) RQ1(j)=rmus*dsqrt(Rip(2,1)*Rip(2,1)+Rip(3,1)*Rip(3,1))
!		  Write(*,*) 'surface 1: ',fiv(j), RL1(j),rmus*Rip(2,1),rmus*Rip(3,1)
!         end if ! xxx .ge. 0.0 ! do not delete
      end do ! j
!
      return
      end subroutine surface1

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine aerosol1 (keyEl,itest,     & ! IN
                           NM,NT,NBV,NAV,   &
                           NANG1,ANGL1,     &
                           rmus,thd,vis,    &
                           EXTH,WD,         &
                           betm,gamm,       &
                           PF11,PF12,       &
							   							   							   							 
                           SL1,SQ1          & ! OUT
                          ) 

!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------
      use mod_par_OS, only : NG0T,KSD,NBVM,KNT,NMM,NMG
      use mod_intrpl_linear
!      use mod_intrpl_spline
	  	  	  	    
      implicit none

      integer, parameter  ::  KANG=2*NG0T+1

!	------------------------------------------------------------------------------------
! IN:
      integer,                          intent(in)  ::  keyEl,NM,NT,NBV,NAV,itest
      integer,                          intent(in)  ::  NANG1
      real,dimension(KANG),             intent(in)  ::  ANGL1
      real,dimension(KNT-1,NMM+NMG),    intent(in)  ::  WD
      real,dimension(KNT),              intent(in)  ::  EXTH
      real,dimension(2*NBVM),           intent(in)  ::  thd,vis
      real,                             intent(in)  ::  rmus,betm,gamm
      real,dimension(KANG,KSD),         intent(in)  ::  PF11,PF12
!	------------------------------------------------------------------------------------
! OUT:
      real,dimension(2*NBVM),           intent(out) ::  SL1,SQ1
!	------------------------------------------------------------------------------------
! LOCAL:
      real,dimension(NM+1)       ::  f11,f12
      integer                    ::  j,m,n
      real                       ::  xx,xi,xp,xxx,YY,WW,VV
!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------
!write(*,*)						   
!write(*,*) 'NM,NT,NBV,NAV,NANG1:'
!write(*,*) NM,NT,NBV,NAV,NANG1
!write(*,*) 'ANGL1:'
!write(*,'(10e16.6)') ANGL1
!write(*,*) 'PF11:'
!write(*,'(10e16.6)') PF11
!write(*,*) 'PF12:'
!write(*,'(10e16.6)') PF12
!write(*,*) 'rmus:'
!write(*,'(10e16.6)') rmus
!write(*,*) 'thd:'
!write(*,'(10e16.6)') thd
!write(*,*) 'vis:'
!write(*,'(10e16.6)') vis
!write(*,*) 'EXTH:'
!write(*,'(10e16.6)') EXTH
!write(*,*) 'WD:'
!do j=1,NMM+NMG
!write(*,'(10e16.6)') WD(1:KNT-1,j)
!write(*,*)
!enddo 
!write(*,*) 'betm:'
!write(*,'(e16.6)') betm
!write(*,*) 'gamm:'
!write(*,'(e16.6)') gamm

!write(*,*)
	  	  	  
!PL   Single scattering + surface (begin)
      SL1(:)=0.
      if (keyEl .eq. 4) SQ1(:)=0.
      do j=1,NBV
         xxx=cos(vis(j)*rpi_180)
!PL      xx is recalculated several times!!!
         xx=cos(thd(j)*rpi_180)
!PL      f11(j,NM+1) - Reyleigh scattering for viewing geometry
         f11(NM+1)=1.0+betm*(3.*xx*xx-1.)*0.5
         do m=1,NM
!PL         Interpolation for viewing geometry
!PL         With option ITRONC.ne.0 single scattering is not calculated accurately
!PL         in forward direction,
!PL         since PF11 is already modified with trancation!!!
            f11(m)=LINEAR_LN(ANGL1(1:NANG1),PF11(1:NANG1,m),NANG1,thd(j))
         end do ! m

!XH      for vector case
         if (keyEl .eq. 4) then
            f12(NM+1)=sqrt(3./8.)*(1.-xx*xx)*gamm
            do m=1,NM
               f12(m)=LINEAR(ANGL1(1:NANG1),PF12(1:NANG1,m),NANG1,thd(j))
            end do ! m
         end if

         if (xxx .gt. 0.0 .and. NAV .lt. NT) then
!           if(itest .eq. +1) then
!PL         For airborne and satellite observation
            do n=NAV,NT-1
               WW=(EXTH(n+1)-EXTH(n))*(1./rmus+1./xxx)
               XX=(1.0-exp(-WW))/(1.0+xxx/rmus)
               YY=exp(-EXTH(n)/rmus)*exp(-(EXTH(n)-EXTH(NAV))/xxx)
!PL            Single aerosol+Reyleigh scattering
               xi=dot_product(f11,WD(n,1:NM+1))
               SL1(j)=SL1(j)+XX*YY*xi*0.25
!XH            for vector case
               if (keyEl .eq. 4) then
                  xp=dot_product(f12,WD(n,1:NM+1))
                  SQ1(j)=SQ1(j)-XX*YY*xp*0.25
               end if
            end do ! n
         end if ! xxx .gt. 0.0 .and. NAV .lt. NT

         if (xxx .lt. 0.0)then
!           if(itest .eq. -1)then
!PL         For ground based observation
            do n=1,NT-1
               WW=(EXTH(n+1)-EXTH(n))*(1./rmus+1./xxx)
               XX=(1.0-exp(-WW))/(1.0+xxx/rmus)
!AL            inserted correction xxx instead of xxx (typo?)
!	            W=(EXTH(n+1)-EXTH(n))/xxx
               VV=(EXTH(n+1)-EXTH(n))/xxx
!c             CORRECTION ALMUCANTAR
!c             Erreur 1 corrigee le 28 mai 2009
!c              if(rmus+xxx.lt.0.001) XX=(EXTH(n+1)-EXTH(n))/xxx
               if (abs(rmus+xxx) .lt. 0.001) XX=(EXTH(n+1)-EXTH(n))/xxx
!c             Erreur 2 corrigee le 28 mai 2009
!c              YY=-exp(-EXTH(n)/rmus)*exp((EXTH(NT)-EXTH(n))/xxx)
               YY=-exp(-EXTH(n)/rmus)*exp((EXTH(NT)-EXTH(n+1))/xxx)
!AL            inserted correction
               XX=XX*exp(VV)
!PL            Single aerosol+Reyleigh scattering
               xi=dot_product(f11,WD(n,1:NM+1))
               SL1(j)=SL1(j)+XX*YY*xi*0.25
!XH            for vector case
               if (keyEl .eq. 4) then
                  xp=dot_product(f12,WD(n,1:NM+1))
                  SQ1(j)=SQ1(j)-XX*YY*xp*0.25
               end if
!PL            The sign changes: -XX*YY*xp*0.25. So SQ1 for Reyleigh is > 0 !!!
!PL            Because of this cs (ANGTURN) is oposit !!!!!
            end do ! n
         end if ! xxx .lt. 0.0
      end do ! j
!
      return
      end subroutine aerosol1
	  
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE phase_matrix_intrpl_NQDR (                                 &
                                             keyEl,                         & ! IN
                                             NM,NQDR,                       &
                                             zmu,zg,                        &
                                             KMpar,NANG,ANGL,               &
                                             PF11_I,PF12_I,PF22_I,PF33_I,   &
											
                                             NANG1,ANGL1,                   & ! OUT
                                             PF11,PF12,PF22,PF33            &
                                          )

!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------
      use mod_par_OS, only : NG0T,KSD
      use mod_intrpl_spline
	  	  	    
      implicit none
	  
      integer, parameter  ::  KANG=2*NG0T+1
!	------------------------------------------------------------------------------------
! IN :
      integer,                           intent(in)  ::  keyEl,NM,NQDR,NANG,KMpar
      real,dimension(KMpar),             intent(in)  ::  ANGL
      real,dimension(KMpar,KSD),         intent(in)  ::  PF11_I,PF12_I,PF22_I,PF33_I
      REAL*8,DIMENSION(-NG0T:NG0T),      intent(in)  ::  zmu,zg

!	--------------------------------------------------------------------------------------
! OUT :
      integer,                          intent(out)  ::  NANG1
      REAL,DIMENSION(KANG),             intent(out)  ::  ANGL1
      REAL,DIMENSION(KANG,KSD),         intent(out)  ::  PF11,PF12,PF22,PF33

!	--------------------------------------------------------------------------------------
! LOCAL :
      DOUBLE PRECISION :: XARG, YFIT
      DOUBLE PRECISION, DIMENSION(NANG) :: XS, YS, b, c, d
      REAL ::  zz0
      INTEGER :: I,j,k,ISD
!	------------------------------------------------------------------------------------------------------
!	-------------------------------------------------------------------------------------------------------
      NANG1=0
      ANGL1(:) =0.0
      PF11(:,:)=0.0
      if (keyEl .eq. 4) then
         PF12(:,:)=0.0
         PF22(:,:)=0.0
         PF33(:,:)=0.0
      end if

      if (NANG .gt. KANG) then
         write(tmp_message,'(2(a,i0))') 'NANG = ',NANG,' .gt. KANG = ',KANG
         G_ERROR(trim(tmp_message))
      endif

!C CORRIGER D'ABORD LA MATRICE PFij?
!C Tres peu d'influence. Modifie de qqs 10-3 en relatif

!PL   ISD is aerosol component
      DO ISD=1,NM
!XH      for P11
         IF (ANGL(2).GT.ANGL(1)) THEN
            XS(1:NANG)=ANGL(1:NANG)
            YS(1:NANG)=log(PF11_I(1:NANG,ISD))
         ELSE
            XS(1:NANG)=ANGL(NANG:1:-1)
            YS(1:NANG)=log(PF11_I(NANG:1:-1,ISD))
         END IF
!PL      Calling for Phase matrix Interpolation
!PL      Pase matrix Interpolation. Result: YFIT(PFij) at
!PL      Gaussian quadratures zmu
         DO j=-NQDR,NQDR
            XARG=acos(zmu(j))*r180_pi
            IF(XARG .GT. 180.0) XARG=180.0
            ANGL1(NQDR+j+1) = XARG
            CALL intrpl_spline(NANG, XS, YS, XARG, YFIT, 1, NQDR+j+1, b, c, d)
            PF11(NQDR+j+1,ISD)=exp(YFIT)
         END DO !j

!XH      for other elements in vector case
         if (keyEl .eq. 4) then
!XH         for P12
            IF (ANGL(2).GT.ANGL(1)) THEN
               YS(1:NANG)=PF12_I(1:NANG,ISD)
            ELSE
               YS(1:NANG)=PF12_I(NANG:1:-1,ISD)
            END IF
!PL         Calling for Phase matrix Interpolation
!PL         Phase matrix Interpolation. Result: YFIT(PFij) at
!PL         Gaussian quadratures zmu
            DO j=-NQDR,NQDR
               XARG=acos(zmu(j))*r180_pi
               IF(XARG .GT. 180.0) XARG=180.0
               ANGL1(NQDR+j+1) = XARG
               CALL intrpl_spline(NANG, XS, YS, XARG, YFIT, 1, NQDR+j+1, b, c, d)
               PF12(NQDR+j+1,ISD)=YFIT
            END DO !j
!XH         for P22
            IF (ANGL(2).GT.ANGL(1)) THEN
               YS(1:NANG)=log(PF22_I(1:NANG,ISD))
            ELSE
               YS(1:NANG)=log(PF22_I(NANG:1:-1,ISD))
            END IF
!PL         Calling for Phase matrix Interpolation
!PL         Pase matrix Interpolation. Result: YFIT(PFij) at
!PL         Gaussian quadratures zmu
            DO j=-NQDR,NQDR
               XARG=acos(zmu(j))*r180_pi
               IF(XARG .GT. 180.0) XARG=180.0
               ANGL1(NQDR+j+1) = XARG
               CALL intrpl_spline(NANG, XS, YS, XARG, YFIT, 1, NQDR+j+1, b, c, d)
               PF22(NQDR+j+1,ISD)=exp(YFIT)
            END DO !j
!XH         for P33
            IF (ANGL(2).GT.ANGL(1)) THEN
               YS(1:NANG)=PF33_I(1:NANG,ISD)
            ELSE
               YS(1:NANG)=PF33_I(NANG:1:-1,ISD)
            END IF
!PL         Calling for Phase matrix Interpolation
!PL         Pase matrix Interpolation. Result: YFIT(PFij) at
!PL         Gaussian quadratures zmu
            DO j=-NQDR,NQDR
               XARG=acos(zmu(j))*r180_pi
               IF(XARG .GT. 180.0) XARG=180.0
               ANGL1(NQDR+j+1) = XARG
               CALL intrpl_spline(NANG, XS, YS, XARG, YFIT, 1, NQDR+j+1, b, c, d)
               PF33(NQDR+j+1,ISD)=YFIT
            END DO !j
         end if
      END DO !PL ISD

!PL   Check for phase matrix normalization.
!PL   In ideal case Phase Matrix must be normalazed is such way that zz0=1.
!PL   To insure this the renormalization is performed here
      NANG1=2*NQDR+1
      do ISD=1,NM
         zz0 = 0.5*dot_product(PF11(1:NANG1,ISD),zg(-NQDR:NQDR))
         PF11(1:NANG1,ISD)=PF11(1:NANG1,ISD)/zz0
         if (keyEl .eq. 4) then
            PF12(1:NANG1,ISD)=PF12(1:NANG1,ISD)/zz0
            PF22(1:NANG1,ISD)=PF22(1:NANG1,ISD)/zz0
            PF33(1:NANG1,ISD)=PF33(1:NANG1,ISD)/zz0
         end if
      end do

      END SUBROUTINE phase_matrix_intrpl_NQDR




!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE MODIF_DE_TRONC_new (                                &
                                      keyEl,                         & ! IN
                                      NM,NQDR,itest,                 &
                                      zmu,zg,NANG1,ANGL1,            &
                                      NBV,thd,rmus,                  &

                                      tetot,EXT,SSA,                 & ! INOUT
                                      PF11,PF12,PF22,PF33,           &
      								  
                                      DDL1,DDQ1                      & ! OUT
                                    )

!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------

      use mod_par_OS, only : NG0T,KSD,NBVM
	  	  	    
      implicit none
	 
       integer, parameter  ::  KANG=2*NG0T+1
!	------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  keyEl,NM,NQDR,itest
      real*8,dimension(-NG0T:NG0T), intent(in)  ::  zmu,zg
      integer,                      intent(in)  ::  NBV
      real,dimension(2*NBVM),       intent(in)  ::  thd
      real,                         intent(in)  ::  rmus
      integer,                      intent(in)  ::  NANG1
      real,dimension(KANG),         intent(in)  ::  ANGL1
!	--------------------------------------------------------------------------------------
! OUT :
      real,dimension(NM),               intent(inout) ::  EXT,SSA
      real,                             intent(inout) ::  tetot
      real,dimension(KANG,KSD),         intent(inout) ::  PF11,PF12,PF22,PF33
      real,dimension(2*NBVM),           intent(out)   ::  DDL1,DDQ1
!	--------------------------------------------------------------------------------------
!	--------------------------------------------------------------------------------------
! LOCAL :
	   REAL                         ::  tdav,tdav1,tevtot
	   REAL,DIMENSION(KANG)         ::  P11av,P12av,P22av,P33av,Vav,Qav,Uav,Wav
 
!CD MODIFICATION SI TRONCATURE 0
	   REAL,DIMENSION(-NG0T:NG0T)   ::  Q11,Q12,Q22,Q33

	   REAL                         ::  zz0,stron,textnew,a,albnew,  & ! pi,
	                                    ANGTRONC,b,pente,tpointe,TT,xx,yy  
       INTEGER                      ::  j,k,ISD,JMAX
!	------------------------------------------------------------------------------------------------------
       REAL,DIMENSION(0:2*NG0T)      ::  alp,bet,gam,zet

	   real       :: ck,dd,ee,                   & 
                     p0,pp0,pp1,p1,p2,pp2,       &
                     ppri,pav,                   &
                     qpri,qqav,                  &
                     yyy,zzz
!PL  
      Integer, parameter:: N_od_tr=4
	  real Yav(2*NBVM,N_od_tr),QYav(2*NBVM,N_od_tr),coeff(N_od_tr)
					 
!	-------------------------------------------------------------------------------------------------------
      !pi = 3.141592653
      alp(:)=0
      bet(:)=0
      gam(:)=0
      zet(:)=0

      tevtot = tetot
	  
!C    CORRIGER D'ABORD LA MATRICE PFij?
!C    Tres peu d'influence. Modifie de qqs 10-3 en relatif

!c    MODIFICATION DE TRONCATURE 3
!c    Initialiser la pointe avant

      P11av(1:NANG1)=0.0
      if (keyEl .eq. 4) then
         P12av(1:NANG1)=0.0
         P22av(1:NANG1)=0.0
         P33av(1:NANG1)=0.0
	  endif

!c    On tronque a 16‚àû (MODIFIABLE-PARAMETRER?)
!PL   ANGTRONC is the angle of phase matrix truncation
      ANGTRONC=16.0*rpi_180

      k=NQDR	  
      do while(acos(zmu(k)) .lt. ANGTRONC)
         k=k-1
      enddo
       	  
      JMAX=k+NQDR+1
!PL   JMAX defines the maximum number of angles before the truncation
      tpointe=0.0
      tdav=0.0

!c    Boucle sur les modeles
!PL   additional trancated elements are added
      do ISD=1,NM
!c       Allure avant du modele
         xx=PF11(JMAX,ISD)
         yy=PF11(JMAX-1,ISD)
         pente=(log(xx)-log(yy))/(ANGL1(JMAX)-ANGL1(JMAX-1))
         b=-pente*0.5/ANGL1(JMAX)
         a=log(PF11(JMAX,ISD))+b*ANGL1(JMAX)*ANGL1(JMAX)
!c       Matrice tronquee
         do j=1,NANG1
            if(j.lt.JMAX)   TT=PF11(j,ISD)
            if(j.gt.JMAX-1) TT=exp(a-b*ANGL1(j)*ANGL1(j))
            Vav(j)=PF11(j,ISD)-TT
            if (keyEl .eq. 4) then
               Qav(j)=Vav(j)*PF12(j,ISD)/PF11(j,ISD)
!PL            additional trancated elements F22a,F33a
               Uav(j)=Vav(j)*PF22(j,ISD)/PF11(j,ISD)
               Wav(j)=Vav(j)*PF33(j,ISD)/PF11(j,ISD)
!PL            Truncated phase matrix
               PF12(j,ISD)=PF12(j,ISD)*TT/PF11(j,ISD)
               PF22(j,ISD)=PF22(j,ISD)*TT/PF11(j,ISD)
               PF33(j,ISD)=PF33(j,ISD)*TT/PF11(j,ISD)
            end if			 
            PF11(j,ISD)=TT
         enddo ! j
         !do j=1,NANG
         !   write(6,*)ANGL(j),Vav(j),Qav(j),PF11(j,ISD)
         !enddo
!        Re-normalisation
		 zz0=0.5*dot_product(PF11(1:NQDR+NQDR+1,ISD),zg(-NQDR:NQDR))
         if (zz0 .ne. 1.0) then
		    stron=zz0
            !write(6,*)'stron',stron
		    do j=1,NANG1
               PF11(j,ISD)=PF11(j,ISD)/zz0
			   Vav(j)=Vav(j)/abs(1.0-zz0)
               if (keyEl .eq. 4) then
			      PF12(j,ISD)=PF12(j,ISD)/zz0
			      PF22(j,ISD)=PF22(j,ISD)/zz0
			      PF33(j,ISD)=PF33(j,ISD)/zz0
			
			      Qav(j)=QAV(j)/abs(1.0-zz0)
			      Uav(j)=Uav(j)/abs(1.0-zz0)
			      Wav(j)=WAV(j)/abs(1.0-zz0)
		       endif
		    enddo ! j
         else
            stron=1.0
         endif ! zz0
				
         !write(6,*) "zz0=",zz0
         !do j=1,NANG1
         !   write(6,*) ANGL1(j),PF11(j,ISD),PF12(j,ISD) !,PF22(j,ISD),PF33(j,ISD)
         !enddo
!c       Correction des epaisseurs optiques abs et diff
		 textnew=EXT(ISD)*(1.0-SSA(ISD)+stron*SSA(ISD))
		 tdav1=EXT(ISD)*SSA(ISD)*(1.0-stron)
		 albnew=stron*SSA(ISD)/(1.0-SSA(ISD)+stron*SSA(ISD))
		 EXT(ISD)=textnew
		 SSA(ISD)=albnew
		 tpointe=tpointe+tdav1

         P11av(1:NANG1)=P11av(1:NANG1)+tdav1*Vav(1:NANG1)
         if (keyEl .eq. 4) then
		    P12av(1:NANG1)=P12av(1:NANG1)+tdav1*Qav(1:NANG1)
		    P22av(1:NANG1)=P22av(1:NANG1)+tdav1*Uav(1:NANG1)
		    P33av(1:NANG1)=P33av(1:NANG1)+tdav1*Wav(1:NANG1)
		 endif

		 tdav=tdav+tdav1

	  enddo ! ISD

      tetot=tetot-tpointe

      if (itest .eq. +1) RETURN
!     Ground based observations

  	  Q11(:)=0.0
      Q11(-NQDR:NANG1-NQDR-1)=P11av(1:NANG1)/tpointe
	  if (keyEl .eq. 4) then
         Q12(:)=0.0
         Q22(:)=0.0
         Q33(:)=0.0
		 Q12(-NQDR:NANG1-NQDR-1)=P12av(1:NANG1)/tpointe
		 Q22(-NQDR:NANG1-NQDR-1)=P22av(1:NANG1)/tpointe
		 Q33(-NQDR:NANG1-NQDR-1)=P33av(1:NANG1)/tpointe
      end if

!c    Correction des mesures dans la pointe tronquee
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  
!PL   1). It is necessary to check wether formula (6) in Moris description is correct!!!!!!!
!         ck=(exp(-tetot/rmus)-exp(-tevtot/rmus))*0.25
!PL   Expansion of the trancated peak
!PL   This expansions are initialized just for ground measurements
!PL   when cos(vis(iv))<0 (for sunphotometer)
!PL   the number of expansion coefficients is 2*NQDR-2
!PL   Therefore the tha range of arrays xalp,xbet,xgam,xzet
!PL   can be (0:2*NQDR-2) !!! The corresponding changes should be
!PL   also done in subroutine betal
!PL   new Resulting forward scattering phase functions Yav(n), Qav(n), for n=2, 3..
!PL   new scatterings in the forward peak
	  
      DDL1(:) = 0.0
      DDQ1(:) = 0.0
      call betal( keyEl,NQDR,NG0T,zmu,zg,Q11,Q12,Q22,Q33,                                  &   ! in
                  alp(0:2*NG0T-2),bet(0:2*NG0T-2),gam(0:2*NG0T-2),zet(0:2*NG0T-2) )            ! out
      call CONVOL_PAV_PL( keyEl,ANGTRONC,N_od_tr,NG0T,2*NQDR-2,alp,bet,gam,                &
                          2*NBVM,NBV,thd*rpi_180,Yav,QYav)
 
      coeff(1)=1.0
	  do k=2,N_od_tr
         coeff(k)=coeff(k-1)*tdav/float(k)/rmus
	  enddo !k

!PL	  Write(*,*) 'coef', coeff,CTOT
!PL   old: ck=(exp(-tetot/rmus)-exp(-tevtot/rmus))*0.25/CTOT
!PL   new
      ck=tdav*0.25/rmus*exp(-tevtot/rmus)

      do j=1,NBV
         !Write(*,*),'Yav ',thd(j),Yav(j,1),Yav(j,2),Yav(j,3),Yav(j,4)!,QYav(j,2)
	  
!PL      Correction of I (scalar and vector case)
         pav=Yav(j,1)
	     do k=2,N_od_tr
            pav=pav+Yav(j,k)*coeff(k)
         enddo !k
         DDL1(j)=pav*ck

!PL      Correction of Q (vector case)
         if (keyEl .eq. 4) then
		    qqav=QYav(j,1)
	        do k=2,N_od_tr
               qqav=qqav+QYav(j,k)*coeff(k)
            enddo !k
            DDQ1(j)=-qqav*ck
         endif
         !Write(*,*) 'coef ', coeff(1),coeff(2),coeff(3)
         !Write(*,*),'QYav ',thd(j),QYav(j,1),QYav(j,2),QYav(j,3),QYav(j,4)!,QYav(j,2)
      enddo !j

      RETURN
      END SUBROUTINE MODIF_DE_TRONC_new


!PL   Convolution calculation
      Subroutine CONVOL_PAV_PL( keyEl, ANGTR,NDAV,NP,NBUAV,             &
                                alpav,betav,gamav,NV,NBV,thd,Yav,Qav )
!
!     Lumiere diffusee n fois dans la pointe; polar. paral. ou perp. au plan de diffusion:
!     L(n)=K(n)*p(n,thd); Lpol(n)=K(n)*q(n,thd) (K(n):cf. sub.AUREOLE);
!     p(n)=S bet(n,l)P(l,thd); q(n)=S gam(n,l)P2(l,thd)
!     bet(1,l)=betav(l); gam(1,l)=gamav(l)
!     bet(n+1,l)=bet(n,l)*betav/(2*l+1)
!     gam(n+1,l)=[bet(n,l)*gamav(l)+gam(n,l)*alpav(l)]/(2*l+1)
!
       implicit none
	   integer, intent(in) :: keyEl,NDAV,NP,NBUAV,NV,NBV
	   real,    intent(in) :: ANGTR,thd(NV)
	   real,dimension(0:2*NP),  intent(in) :: alpav,betav,gamav
       real,dimension(NV,NDAV), intent(out):: Yav,Qav

       real, dimension(NDAV,0:2*NP):: gam,bet 

	   integer j,k,l,m
	   real pa,pb,pc,xx,dd,ee

!PL    pi = 3.141592653
       bet(1,:)=betav
       do k=2,NDAV
          do l=0,NBUAV
             bet(k,l)=bet(k-1,l)*betav(l)/(2*l+1)
          enddo
       enddo
       Yav(:,:)=0.0

!PL    For vector case
       if (keyEl .eq. 4) then
          gam(1,:)=gamav(:)
          do k=2,NDAV
             do l=0,NBUAV
                gam(k,l)=(bet(k-1,l)*gamav(l)+gam(k-1,l)*alpav(l))/(2*l+1)
             enddo
          enddo
          Qav(:,:)=0.0
       endif

       do j=1,NBV
!PL       if(thd(j) .le. 2*ANGTR) then
          xx=cos(thd(j))
!PL       Correction of I
          pa=0.
          pb=1.
          do k=0,NBUAV
             pc=((2*k+1.)*xx*pb-(k)*pa)/(k+1)
             do m=1,NDAV
                Yav(j,m)=Yav(j,m)+bet(m,k)*pb
             enddo
             pa=pb
             pb=pc
          enddo !k

          if (keyEl .eq. 4) then
             pa=0.0
             pb=3.0*(1.0-xx*xx)/2.0/sqrt(6.0)
             do k=2,NBUAV
                dd=(2*k+1.)/sqrt((k+3.)*(k-1.))
                ee=sqrt((k+2.)*(k-2.))/(2*k+1.)
                pc=dd*(xx*pb-ee*pa)
                do m=1,NDAV
                   Qav(j,m)=Qav(j,m)+gam(m,k)*pb
                enddo
                pa=pb
                pb=pc
             enddo !k
		  endif ! keyEl .eq. 4
!         endif ! thd(j) .le. 2*ANGTR
       enddo !j
       return
       end Subroutine CONVOL_PAV_PL
	  

!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

       SUBROUTINE MODIF_DE_TRONC (                               &
                                  keyEl,NM,NQDR,itest,           & ! IN
                                  zmu,zg,NANG1,ANGL1,            &
                                  NBV,thd,rmus,                  &

                                  tetot,EXT,SSA,                 & ! INOUT
                                  PF11,PF12,PF22,PF33,           &
      								  
                                  DDL1,DDQ1                      & ! OUT
                                 )

!	------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------
!XH    optional parameters for vector case
       use mod_par_OS, only : NG0T,KSD,NBVM
	  	  	    
       implicit none
	  
       integer, parameter  ::  KANG=2*NG0T+1
!	------------------------------------------------------------------------------------
! IN :
      integer,                      intent(in)  ::  keyEl,NM,NQDR,itest
      real*8,dimension(-NG0T:NG0T), intent(in)  ::  zmu,zg
      integer,                      intent(in)  ::  NBV
      real,dimension(2*NBVM),       intent(in)  ::  thd
      real,                         intent(in)  ::  rmus
      integer,                      intent(in)  ::  NANG1
      real,dimension(KANG),         intent(in)  ::  ANGL1
!	--------------------------------------------------------------------------------------
! OUT :
      real,dimension(NM),               intent(inout) ::  EXT,SSA
      real,                             intent(inout) ::  tetot
      real,dimension(KANG,KSD),         intent(inout) ::  PF11,PF12,PF22,PF33
      real,dimension(2*NBVM),           intent(out)   ::  DDL1,DDQ1
!	--------------------------------------------------------------------------------------
! LOCAL :
	   REAL                         ::  tdav,tdav1,tevtot
	   REAL,DIMENSION(KANG)         ::  P11av,P12av,Vav,Qav
 
!CD MODIFICATION SI TRONCATURE 0
	   REAL,DIMENSION(-NG0T:NG0T) ::  Q11,Q12,Q22,Q33

	   REAL                         ::  zz0,stron,textnew,a,albnew,      &
	                                    ANGTRONC,b,pente,tpointe,TT,xx,yy  
       INTEGER                      ::  j,k,ISD,JMAX
!	------------------------------------------------------------------------------------------------------
!XH    modified upper bound of arrays from 2*NG0T to 2*NG0T-2
!       REAL,DIMENSION(0:2*NG0T)      ::  alp,bet,gam,zet
       REAL,DIMENSION(0:2*NG0T-2)    ::  alp,bet,gam,zet

       real       ::  ck,dd,ee,                   &
                      p0,pav,pp0,pp1,p1,p2,pp2,   &
                      ppri,psec,ptri,             &
                      qpri,qqav,qsec,             &
                      yyy,zzz

!	------------------------------------------------------------------------------------------------------
!	-------------------------------------------------------------------------------------------------------
      
!write(*,*) 'tetot,EXT,SSA: ',tetot,EXT,SSA	
      tevtot = tetot
	  
!C CORRIGER D'ABORD LA MATRICE PFij?
!C Tres peu d'influence. Modifie de qqs 10-3 en relatif

!c  MODIFICATION DE TRONCATURE 3
!c     Initialiser la pointe avant
      tpointe=0.0
!XH   initialize P12av only for vector case
      P11av(1:NANG1)=0.0
      if (keyEl .eq. 4) P12av(1:NANG1)=0.0
!c On tronque a 16‚àû (MODIFIABLE-PARAMETRER?)
!PL ANGTRONC is the angle of phase matrix truncation
      ANGTRONC=16.0*rpi_180

      k=NQDR	  
      do while(acos(zmu(k)) .lt. ANGTRONC)
      k=k-1
      enddo
       	  
      JMAX=k+NQDR+1
!PL JMAX defines the maximum number of angles before the truncation	  

      tdav=0.0

!c    Boucle sur les modeles
      do ISD=1,NM
!c       Allure avant du modele
         xx=PF11(JMAX,ISD)
         yy=PF11(JMAX-1,ISD)
         pente=(log(xx)-log(yy))/(ANGL1(JMAX)-ANGL1(JMAX-1))
         b=-pente*0.5/ANGL1(JMAX)
         a=log(PF11(JMAX,ISD))+b*ANGL1(JMAX)*ANGL1(JMAX)
!c       Matrice tronquee
!XH      only need to calculate PF11,Vav for scalar case
         do j=1,NANG1
            if(j .lt. JMAX)   TT=PF11(j,ISD)
            if(j .gt. JMAX-1) TT=exp(a-b*ANGL1(j)*ANGL1(j))
            Vav(j)=PF11(j,ISD)-TT
            if (keyEl .eq. 4) then
               Qav(j)=Vav(j)*PF12(j,ISD)/PF11(j,ISD)
               PF12(j,ISD)=PF12(j,ISD)*TT/PF11(j,ISD)
               PF22(j,ISD)=PF22(j,ISD)*TT/PF11(j,ISD)
               PF33(j,ISD)=PF33(j,ISD)*TT/PF11(j,ISD)
            end if
            PF11(j,ISD)=TT
         end do ! j
!c       do j=1,NANG
!c          write(6,*)ANGL(j),Vav(j),Qav(j),PF11(j,ISD)
!c       enddo
!c       Re-normalisation
!XH      perform renormalization in vector form
!XH      only need to renormalize PF11,Vav for scalar case
         zz0=0.5*dot_product(PF11(1:NQDR+NQDR+1,ISD),zg(-NQDR:NQDR))
!         zz0=0.0
!         do j=-NQDR,NQDR
!            k=j+NQDR+1
!            zz0=zz0+0.5*PF11(k,ISD)*zg(j)
!         enddo
         if (zz0.ne.1.0) then
            stron=zz0
            do j=1,NANG1
               PF11(j,ISD)=PF11(j,ISD)/zz0
               Vav(j)=Vav(j)/abs(1.0-zz0)
               if (keyEl .eq. 4) then
                  PF12(j,ISD)=PF12(j,ISD)/zz0
                  PF22(j,ISD)=PF22(j,ISD)/zz0
                  PF33(j,ISD)=PF33(j,ISD)/zz0
                  Qav(j)=QAV(j)/abs(1.0-zz0)
                end if
            end do ! j
         else
            stron=1.0
         endif ! zz0
				
         !write(6,*) "zz0=",zz0
         !do j=1,NANG1
         !write(6,*) ANGL1(j),PF11(j,ISD),PF12(j,ISD) !,PF22(j,ISD),PF33(j,ISD)
         !enddo
!c       Correction des epaisseurs optiques abs et diff
         textnew=EXT(ISD)*(1.0-SSA(ISD)+stron*SSA(ISD))
         tdav1=EXT(ISD)*SSA(ISD)*(1.0-stron)
         albnew=stron*SSA(ISD)/(1.0-SSA(ISD)+stron*SSA(ISD))
         EXT(ISD)=textnew
         SSA(ISD)=albnew
         tpointe=tpointe+tdav1
!XH      only need to calculate P11av for scalar case
         P11av(1:NANG1)=P11av(1:NANG1)+tdav1*Vav(1:NANG1)
         if (keyEl .eq. 4) P12av(1:NANG1)=P12av(1:NANG1)+tdav1*Qav(1:NANG1)
         tdav=tdav+tdav1
	  end do ! ISD

      tetot=tetot-tpointe

      if(itest .eq. +1) RETURN

!c    Ground based observations
!XH   only need to calculate Q11 for scalar case
      Q11(:)=0.0
      Q11(-NQDR:NANG1-NQDR-1)=P11av(1:NANG1)/tpointe
      if (keyEl .eq. 4) then
         Q12(:)=0.0
         Q22(:)=0.0
         Q33(:)=0.0
         Q12(-NQDR:NANG1-NQDR-1)=P12av(1:NANG1)/tpointe
      end if

!c    Correction des mesures dans la pointe tronquee
      ck=(exp(-tetot/rmus)-exp(-tevtot/rmus))*0.25
!PL   Expansion of the trancated peak
!PL   This expansions are initialized just for ground measurements
!PL   when cos(vis(iv))<0 (for sunphotometer)
!PL   the number of expansion coefficients is 2*NQDR-2
!PL   Therefore the tha range of arrays xalp,xbet,xgam,xzet
!PL   can be (0:2*NQDR-2) !!! The corresponding changes should be
!PL   also done in subroutine betal
!XH   only need to calculate DDL1 for scalar case
      DDL1(:) = 0.0
      if (keyEl .eq. 4) DDQ1(:) = 0.0
      call betal(keyEl,NQDR,NG0T,zmu,zg,Q11,Q12,Q22,Q33,  &   ! in
                 alp,bet,gam,zet                  )           ! out

      yyy=tdav*0.5/rmus
      zzz=tdav*tdav/(6.*rmus*rmus)
!XH   separate the calculation of DDL1 and DDQ1 for vector and scalar case
!XH   need to check the following code: qsec
      qsec=0.0  !tl
      do j=1,NBV
         xx=cos(thd(j)*rpi_180)
         ppri=0.0
         psec=0.0
         ptri=0.0
         p0  =0.0
         p1  =1.0

         do k=0,2*NQDR-2
            p2=((2*k+1)*xx*p1-k*p0)/(k+1)
            ppri=ppri+bet(k)*p1
            psec=psec+bet(k)*bet(k)*p1/(2*k+1)
            ptri=ptri+bet(k)*bet(k)*bet(k)*p1/((2*k+1)*(2*k+1))
            p0=p1
            p1=p2
         end do
         pav=ck*(ppri+yyy*psec+zzz*ptri)/(1.0+yyy+zzz)
         DDL1(j)=pav

         if (keyEl .eq. 4) then
            qpri=0.0
            pp0 =0.0
            pp1 =3.0*(1.0-xx*xx)*0.5/sqrt(6.0)
            do k=2,2*NQDR-2
               dd=(2*k+1.)/sqrt((k+3.)*(k-1.))
               ee=sqrt((k+2.)*(k-2.))/(2*k+1)
               pp2=dd*(xx*pp1-ee*pp0)
               qpri=qpri+gam(k)*pp1
               qsec=qsec+gam(k)*gam(k)*pp1/(2*k+1)
               pp0=pp1
               pp1=pp2
            end do
            qqav=-ck*(qpri+yyy*qsec)/(1.0+yyy)
            DDQ1(j)=qqav
         end if
      end do
!c    FIN DES MODIFS
      RETURN
      END SUBROUTINE MODIF_DE_TRONC
  
!	sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      SUBROUTINE rad_transf_MS (                                         &
                                 keyEl,NN,NG,NF,NBV,NM,                  & ! IN
                                 itest,NAV,NT,EXTH,WD,                   &
                                 fiv,rmu,rg,tetas,thv,                   &
                                 idir,coff,eps_err,                      &
                                 alpm,betm,gamm,                         & 
                                 alp,bet,gam,zet,                        &
!                                 PTR_PSL,PTR_TSL,PTR_RSL,                &
                                 R11,R12,R13,R21,R22,R23,R31,R32,R33,    &
																					
                                 SL,SQos,SUos,                           & ! OUT 
                                 UFG,DFG,UFX,DFX                         &
                               )
!XH   optional flux calculation
      use mod_par_OS
	  	  	    
      implicit none
	  
!	------------------------------------------------------------------------------------
!   IN :
      integer,                               intent(in)  ::  keyEl,NN,NG,NF
      integer,                               intent(in)  ::  itest,NBV,NAV,NT
      integer,                               intent(in)  ::  NM
      INTEGER,DIMENSION(2*NBVM),             intent(in)  ::  idir
      real,                                  intent(in)  ::  tetas,eps_err
      real,dimension(KNT-1,NMM+NMG),         intent(in)  ::  WD
      real,dimension(KNT),                   intent(in)  ::  EXTH
      real,dimension(2*NBVM),                intent(in)  ::  fiv
      REAL(8),DIMENSION(-NN0:NN0),           intent(in)  ::  rmu,rg
      REAL   ,DIMENSION(-NN0:NN0),           intent(in)  ::  thv
      REAL,   DIMENSION(2*NBVM),             intent(in)  ::  coff
!      REAL, POINTER,                         intent(in)  ::  PTR_PSL(:,:,:),  &
!                                                             PTR_RSL(:,:,:),  &
!                                                             PTR_TSL(:,:,:)
      REAL,DIMENSION(NN0,NN0,0:NF0),         intent(in)  ::  R11,R12,R13,R21,R22,R23,R31,R32,R33
      real,                                  intent(in)  ::  betm,alpm,gamm
      real,dimension(0:2*NG0,NMM),           intent(in)  ::  bet,alp,gam,zet
!	--------------------------------------------------------------------------------------
!   OUT :
      REAL,DIMENSION(2*NBVM),                intent(out) ::  SL,SQos,SUos
      real,                                  intent(out) ::  UFG,DFG
      real,dimension(KNT),          optional,intent(out) ::  UFX,DFX
!	--------------------------------------------------------------------------------------
!   LOCAL :
!XH   IM  : upward radiance at observation level (NAV) from 1 to NN; downward radiance at observation level (NAV) from -NN to -1
      REAL,DIMENSION(-NN0:NN0)           ::  IM,QM,UM,I1NAV
      REAL,DIMENSION(-NN0:NN0,NN0,NMM)   ::  P11,P21,P22,P31,P32,P33
      REAL,DIMENSION(-NN0:NN0,KNT-1)     ::  SFI,SFQ,SFU
      REAL,DIMENSION(-NN0:NN0,KNT)       ::  I1,Q1,U1
      REAL,DIMENSION(0:2*NG0,-NN0:NN0)   ::  PSL,RSL,TSL
!XH   IZX : sum of all order-of-scattering terms for 0 Fourier component at each level
!XH   IZG : sum of all order-of-scattering terms for 0 Fourier component at ground level
!      REAL,DIMENSION(-NN0:NN0,KNT)       ::  IZX
      REAL,DIMENSION(:,:),allocatable    ::  IZX
      REAL,DIMENSION(-NN0:NN0)           ::  IZG
!      REAL,DIMENSION(-NN0:NN0,-NN0:NN0,KNT)   ::   B11, B21, B22, B31, B32, B33

!PL   Because of reciprocity principle (Rij(j,k,IS)=Rji(k,i,IS)),
!PL   it is possible to reduce the number of elements Rij
!PL   But for anizotropy surfaces this can not be done
!	------------------------------------------------------------------------------------------------------

      integer    ::  k,iv,j,IS,m,l,n0,ND,NMAX
      real       ::  aj,ak,bj,bk,        &
                     PP,RP,TP,TT,RR,TR
      real       ::  b11,b12,b13,        &
                     b21,b22,b23,        &
                     b31,b32,b33
!      integer    ::  j2,j3,k2,k3
!	-------------------------------------------------------------------------------------------------------
      real       ::  ang,                &
                     C,cths,rmus,        &
                     RS,RT,U,            &
                     X,xi,XQ,XU,xx,      &
                     Y,YI,YQ,YU,z,       &
                     RS1
!PL
!tl      real, parameter  :: eps_err=0.0001 !PL original eps_err=0.0001
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------

!      do j=-NN,NN
!         thv(j)=acos(rmu(j))*r180_pi
!         WRITE(*,*) j,thv(j)
!      end do ! j
!      write(*,*) 'in rad_transf_MS: NF : ',NF

      if (present(UFX)) then
!XH      assign array dynamically for flux calculation
         allocate(IZX(-NN:NN,NT))
      end if

!PL   Fourier series
!PL   The number of Fourier expansion term NF is equal to NN
      Do IS=0,NF ! PL and OD changes !Do IS=0,50 !100
!        BEGIN of successive order
!        for IS loop, IS = 50 is never reached
!        Plm(j,k) j:outgoing, k:incident

         if (IS .le. 2) then
!XH         consider Rayleigh scattering for Fourier terms below the second order
            NMAX=NM+1
!PL         IS Fourier component for Reyleigh scattering
            do k=1,NN
               ak=rmu(k)
               bk=sqrt(1.0-ak*ak)
               do j=-NN,NN
                  aj=rmu(j)
                  bj=sqrt(1.0-aj*aj)
                  select case(IS)
                  case(0)
                     P11(j,k,NMAX)=1.0+0.25*betm*(3.0*aj*aj-1.0)*(3.0*ak*ak-1.0)
                     if (keyEl .eq. 4) then
                        P21(j,k,NMAX)=sqrt(0.09375)*gamm*bj*bj*(3.0*ak*ak-1.0)
                        P22(j,k,NMAX)=        0.375*alpm*bj*bj*bk*bk
                        P31(j,k,NMAX)=0.0
                        P32(j,k,NMAX)=0.0
                        P33(j,k,NMAX)=0.0
                     end if
                  case(1)
                     P11(j,k,NMAX)=1.5*betm*aj*bj*ak*bk
                     if (keyEl .eq. 4) then
                        P21(j,k,NMAX)=-sqrt(0.375)*gamm*aj*bj*ak*bk
                        P22(j,k,NMAX)=        0.25*alpm*aj*bj*ak*bk
                        P31(j,k,NMAX)= sqrt(0.375)*gamm   *bj*ak*bk
                        P32(j,k,NMAX)=       -0.25*alpm   *bj*ak*bk
                        P33(j,k,NMAX)=        0.25*alpm   *bj   *bk
                     end if
                  case(2)
                     P11(j,k,NMAX)=0.375*betm*(1.0-aj*aj)*(1.0-ak*ak)
                     if (keyEl .eq. 4) then
                        P21(j,k,NMAX)=0.25*sqrt(0.375)*gamm*(1.0+aj*aj)*(1.0-ak*ak)
                        P22(j,k,NMAX)=          0.0625*alpm*(1.0+aj*aj)*(1.0+ak*ak)
                        P31(j,k,NMAX)=-0.5*sqrt(0.375)*gamm*aj         *(1.0-ak*ak)
                        P32(j,k,NMAX)=          -0.125*alpm*aj         *(1.0+ak*ak)
                        P33(j,k,NMAX)=            0.25*alpm*aj*ak
                     end if
                  case default
                     write(tmp_message,'(a,i0,a)') 'IS = ',IS,'  - unknown value'
                     G_ERROR(trim(tmp_message))
                  end select ! IS
               end do ! j
            end do ! k
         else
!XH         don't have to consider Rayleigh scattering for Fourier terms higher than the second order
		    NMAX=NM
         end if ! IS .le. 2

         if (keyEl .eq. 4) then
!XH         for vector case
            PSL(:,:)=0.0
	        TSL(:,:)=0.0
	        RSL(:,:)=0.0
!PL         Legendre functions calculations
            call legendre(IS,rmu,NG,NN,PSL,TSL,RSL)
!PL         IS Fourier component of aerosol
            do k=1,NN
               do j=-NN,NN
!PL               aerosol components
                  P21(j,k,1:NM)=0.0
                  P22(j,k,1:NM)=0.0
                  P31(j,k,1:NM)=0.0
                  P32(j,k,1:NM)=0.0
                  P33(j,k,1:NM)=0.0
                  do l=IS,NG+NG-2
!                     PP=PTR_PSL(l,j,IS)*PTR_PSL(l,k,IS)
!                     RP=PTR_RSL(l,j,IS)*PTR_PSL(l,k,IS)
!                     TP=PTR_TSL(l,j,IS)*PTR_PSL(l,k,IS)
!                     TT=PTR_TSL(l,j,IS)*PTR_TSL(l,k,IS)
!                     RR=PTR_RSL(l,j,IS)*PTR_RSL(l,k,IS)
!                     TR=PTR_TSL(l,j,IS)*PTR_RSL(l,k,IS)
!                     RT=PTR_RSL(l,j,IS)*PTR_TSL(l,k,IS)
                     RP=RSL(l,j)*PSL(l,k)
                     TP=TSL(l,j)*PSL(l,k)
                     TT=TSL(l,j)*TSL(l,k)
                     RR=RSL(l,j)*RSL(l,k)
                     TR=TSL(l,j)*RSL(l,k)
                     RT=RSL(l,j)*TSL(l,k)

                     P21(j,k,1:NM)=P21(j,k,1:NM)+gam(l,1:NM)*RP
                     P31(j,k,1:NM)=P31(j,k,1:NM)-gam(l,1:NM)*TP
                     P22(j,k,1:NM)=P22(j,k,1:NM)+alp(l,1:NM)*RR+zet(l,1:NM)*TT
                     P33(j,k,1:NM)=P33(j,k,1:NM)+alp(l,1:NM)*TT+zet(l,1:NM)*RR
                     P32(j,k,1:NM)=P32(j,k,1:NM)-alp(l,1:NM)*TR-zet(l,1:NM)*RT
                  end do ! l
               end do ! j
            end do ! k
            QM(-NN:NN)=0.0
            UM(-NN:NN)=0.0
         else
!XH         for scalar case
            PSL(:,:)=0.0
            call legendre(IS,rmu,NG,NN,PSL)
         end if
!XH      for both scalar and vector cases
         IM(-NN:NN)=0.0
         do k=1,NN
            do j=-NN,NN
               P11(j,k,1:NM)=0.0
               do l=IS,NG+NG-2
                  PP=PSL(l,j)*PSL(l,k)
                  P11(j,k,1:NM)= P11(j,k,1:NM)+bet(l,1:NM)*PP
               end do
            end do
         end do
         IM(-NN:NN)=0.0

!PL      Interpolation of IS-Fouries component (Pij) of
!PL      aerosol and Reyligh phase matrices
!PL      for incident geometry (tetas) (PRij)
!PL      Heare the case of passive remote sensing is assumed
!PL      (unpolarized incident radiation). Therefore the inerpolation is
!PL      performed just for the elements P11, P21 and P31!!!
!PL      The sign of P31 and R31 is changed to satisfy angular convention defined 
!PL      for Polder data (Stokes parameter U).	 
!PL      for incident geometry (tetas)
         rmus=cos(tetas*rpi_180)
         if (thv(NN) .lt. tetas) then
            n0 = NN-1
            do while(thv(n0).lt.tetas)
	           n0 = n0-1
            end do
            if (n0 .gt. 0) then
               cths = (rmus-rmu(n0+1))/(rmu(n0)-rmu(n0+1))
            else
!XH            when solar zenith angle is close to 90 degrees
               cths = (rmus-rmu(1))/(-rmu(1))
            end if
         else
!XH         for the exact nadir direction
            n0  =NN-1
            cths=0.0
         end if
!PL      Mixing aerosols and Reyleigh components according to the weights WD
         do l=1,NT-1
            C= 0.125 * ( exp( - EXTH(l) / rmus ) + exp( - EXTH(l+1) / rmus ) )
            do j=-NN,NN
               if (n0 .gt. 0) then
                  SFI(-j,l)=C*dot_product((1.0-cths)*P11(j,n0+1,1:NMAX)+cths*P11(j,n0,1:NMAX),WD(l,1:NMAX))
                  if (keyEl .eq. 4) then
!XH                  for vector case
                     SFQ(-j,l)= C*dot_product((1.0-cths)*P21(j,n0+1,1:NMAX)+cths*P21(j,n0,1:NMAX),WD(l,1:NMAX))
                     SFU(-j,l)=-C*dot_product((1.0-cths)*P31(j,n0+1,1:NMAX)+cths*P31(j,n0,1:NMAX),WD(l,1:NMAX))
                  end if
               else
!XH               when solar zenith angle is close to 90 degrees

!XH               assuming plane-parallel model
                  SFI(-j,l)=C*dot_product((1.0-cths)*P11(j,1,1:NMAX),WD(l,1:NMAX))
                  if (keyEl .eq. 4) then
!XH                  for vector case
                     SFQ(-j,l)= C*dot_product((1.0-cths)*P21(j,1,1:NMAX),WD(l,1:NMAX))
                     SFU(-j,l)=-C*dot_product((1.0-cths)*P31(j,1,1:NMAX),WD(l,1:NMAX))
                  end if

!!XH               assuming a more realistic case
!                  SFI(-j,l)=C*dot_product(P11(j,1,1:NMAX),WD(l,1:NMAX))
!                  if (keyEl .eq. 4) then
!!XH                  for vector case
!                     SFQ(-j,l)= C*dot_product(P21(j,1,1:NMAX),WD(l,1:NMAX))
!                     SFU(-j,l)=-C*dot_product(P31(j,1,1:NMAX),WD(l,1:NMAX))
!                  end if
               end if
            end do ! j
         end do ! l

!XH      calculation for first order scattering

!PL      Interpolation of IS-Fouries component (Rij) of
!PL      the surface reflection matrix
!PL      for incident geometry (tetas).
!PL      Combining aerosol-Ryeligh with surface reflection (I1-U1)
!PL      and integration over TAU for upward direction
         I1(0,NT)=0.0
         I1(0, 1)=0.0
         if (keyEl .eq. 4) then
!XH         for vector case
            Q1(0,NT)=0.0
            Q1(0, 1)=0.0
            U1(0,NT)=0.0
            U1(0, 1)=0.0
         end if

!XH      integration over TAU for upward direction
         RT=EXP(-EXTH(NT)/rmus)
!PL      Eq.(31e) J.Lenoble, M.Herman et al, JQSRT, 2007
!XH      surface contribution to upward direction, interpolate into solar direction
         if (n0 .gt. 0) then
            I1(1:NN,NT)=((1.0-cths)*R11(1:NN,n0+1,IS)+cths*R11(1:NN,n0,IS))*RT
            if (keyEl .eq. 4) then
               Q1(1:NN,NT)=((1.0-cths)*R21(1:NN,n0+1,IS)+cths*R21(1:NN,n0,IS))*RT
               U1(1:NN,NT)=((1.0-cths)*R31(1:NN,n0+1,IS)+cths*R31(1:NN,n0,IS))*RT
            end if
         else
!XH         when solar zenith angle is close to 90 degrees

!XH         assuming plane-parallel model
            I1(1:NN,NT)=((1.0-cths)*R11(1:NN,1,IS))*RT
            if (keyEl .eq. 4) then
               Q1(1:NN,NT)=((1.0-cths)*R21(1:NN,1,IS))*RT
               U1(1:NN,NT)=((1.0-cths)*R31(1:NN,1,IS))*RT
            end if

!!XH         assuming a more realistic case
!            I1(1:NN,NT)=(R11(1:NN,1,IS))*RT
!            if (keyEl .eq. 4) then
!               Q1(1:NN,NT)=(R21(1:NN,1,IS))*RT
!               U1(1:NN,NT)=(R31(1:NN,1,IS))*RT
!            end if
         end if
!PL      Combination of single scattering from surface
!PL      and aerosol-Reyleigh (SFI-SFU) component for each layer.
!PL         Eq.31(b). (single scattering)

         do l=NT-1,1,-1
            do j=1,NN
               C=EXP( (EXTH(l) - EXTH(l+1) ) / RMU(j) )
               I1(j,l)=C*I1(j,l+1)+(1.0-C)*SFI(j,l)
               if (keyEl .eq. 4) then
                  Q1(j,l)=C*Q1(j,l+1)+(1.0-C)*SFQ(j,l)
                  U1(j,l)=C*U1(j,l+1)+(1.0-C)*SFU(j,l)
               end if
            end do ! j
         end do ! l

!PL      integration over TAU for downward direction
!PL      Eq.(31a) (single scattering)
         I1(-NN:-1,1)=0.0
         if (keyEl .eq. 4) then
            Q1(-NN:-1,1)=0.0
            U1(-NN:-1,1)=0.0
         end if
         do l=2,NT
            do j=-NN,-1
               C=EXP( (EXTH(l) - EXTH(l-1) ) / RMU(j))
               I1(j,l)=C*I1(j,l-1)+(1.0-C)*SFI(j,l-1)
               if (keyEl .eq. 4) then
                  Q1(j,l)=C*Q1(j,l-1)+(1.0-C)*SFQ(j,l-1)
                  U1(j,l)=C*U1(j,l-1)+(1.0-C)*SFU(j,l-1)
               end if
            end do ! j
         end do ! l

!XH      contribution of first order scattering to flux
!PL      Single scattering (atmosphere + surface)		 
         if (IS .eq. 0) then
            I1NAV(-NN:NN) = I1(-NN:NN,NAV)
            IZG(  -NN:NN) = I1(-NN:NN,NT )
            if (present(UFX)) then
!XH            IS=0 component is equivalent to integral of radiance over azimuth angle
               IZX( -NN:NN,1:NT)=I1(-NN:NN,1:NT)
            end if
         end if ! IS .eq. 0
!XH      the end of the first order scattering calculation

!c       ****** BOUCLE SUR LES DIFFUSIONS
!c       ****** NOYAUX DE DIFFUSION DES MELANGES  BIJ(j,k,M) j:diffuse, k:incident
!PL      Mixing aerosols and Reyleigh PIJ  components
!PL      according to the weights WD (Bij include single scattering)
!        do 700 l=1,NT-1
!           do 701 k=1,NN
!              do 702 j=-NN,NN
!                 if (j.eq.0) goto 702
!                 B11(j,k,l)=0.0
!                 B21(j,k,l)=0.0
!                 B22(j,k,l)=0.0
!                 B31(j,k,l)=0.0
!                 B32(j,k,l)=0.0
!                 B33(j,k,l)=0.0
!                 do 703 m=1,NMAX
!                    B11(j,k,l)=B11(j,k,l)+P11(j,k,m)*WD(l,m)
!                    B21(j,k,l)=B21(j,k,l)+P21(j,k,m)*WD(l,m)
!                    B22(j,k,l)=B22(j,k,l)+P22(j,k,m)*WD(l,m)
!                    B31(j,k,l)=B31(j,k,l)+P31(j,k,m)*WD(l,m)
!                    B32(j,k,l)=B32(j,k,l)+P32(j,k,m)*WD(l,m)
!                    B33(j,k,l)=B33(j,k,l)+P33(j,k,m)*WD(l,m)
!                 703 continue
!                 B11(-j,-k,l)=B11(j,k,l)
!                 B21(-j,-k,l)=B21(j,k,l)
!                 B22(-j,-k,l)=B22(j,k,l)
!                 B31(-j,-k,l)=-B31(j,k,l)
!                 B32(-j,-k,l)=-B32(j,k,l)
!                 B33(-j,-k,l)=B33(j,k,l)

!                 Because of reciprocity principle
!                 it is possible to remove B12,B13,B23
!PL               B12(k,j,l)=B21(j,k,l)
!PL               B13(k,j,l)=B31(j,k,l)
!PL               B23(k,j,l)=B32(j,k,l)
!PL               B12(-k,-j,l)=B21(-j,-k,l)
!PL               B13(-k,-j,l)=B31(-j,-k,l)
!PL               B23(-k,-j,l)=B32(-j,-k,l)

!              702 continue
!           701 continue
!        700 continue
!        --------------------------
!             j   k
!        B11(-NN,NN)  <=> P11(-NN,NN)
!        B11(-NN,-NN) <=> P11(NN,NN)

!        B11(NN,-NN) <=> P11(-NN,NN)
!        B11(NN,NN)  <=> P11(NN,NN)
 
!        --------------------------

!        B31(-NN,NN)  <=> P31(-NN,NN)
!        B31(-NN,-NN) <=> -P31(NN,NN)

!        B31(NN,-NN) <=> -P31(-NN,NN)
!        B31(NN,NN)  <=> P31(NN,NN)
!        --------------------------
!        B32(-NN,NN)  <=> P32(-NN,NN)
!        B32(-NN,-NN) <=> -P32(NN,NN)

!        B32(NN,-NN) <=> -P32(-NN,NN)
!        B32(NN,NN)  <=> P32(NN,NN)
!        --------------------------
!c       ICI RETOUR ND

!PL      Order of scattering cycle is starting here !!!
!         write(*,*) 'before ND loop in rad_transf_MS'
         ND=1
         LOOP_scat_order: Do
            ND=ND+1
!PL         1.Multiplication of matrix Bij and vector (I1-Q1)
!PL         2.Integration using Gaussian quadrature formulas (rg(k))
!PL           SFI-SFU contain result for different layers
!PL           X=X+U*(B11(j,k,l)*YI+B12(j,k,l)*YQ+B13(j,k,l)*YU)
!PL           Y=Y+U*(B21(j,k,l)*YI+B22(j,k,l)*YQ+B23(j,k,l)*YU)
            if (keyEl .eq. 4) then
!XH            for vector case
               SFI(:,:)=0.0
               SFQ(:,:)=0.0
               SFU(:,:)=0.0
               do l=1,NT-1
!PL Multiple scattering at level "l" for downward direction		    
                  do j=-NN,-1
                     X=0.0
                     Y=0.0
                     Z=0.0
                     do k=1,NN
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        YQ=(Q1(k,l+1)+Q1(k,l))
                        YU=(U1(k,l+1)+U1(k,l))
                        do m=1,NMAX
                           X=X+U*( P11(j,k,m)*YI + P21(-k,-j,m)*YQ - P31(-k,-j,m)*YU )*WD(l,m)
                           Y=Y+U*( P21(j,k,m)*YI + P22( j, k,m)*YQ - P32(-k,-j,m)*YU )*WD(l,m)
                           Z=Z+U*( P31(j,k,m)*YI + P32( j, k,m)*YQ + P33( j, k,m)*YU )*WD(l,m)
                        end do ! m
                     end do ! k
                     do k=-NN,-1
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        YQ=(Q1(k,l+1)+Q1(k,l))
                        YU=(U1(k,l+1)+U1(k,l))
                        do m=1,NMAX
                           X=X+U*( P11(-j,-k,m)*YI + P21(-k,-j,m)*YQ - P31(-k,-j,m)*YU )*WD(l,m)
                           Y=Y+U*( P21(-j,-k,m)*YI + P22(-j,-k,m)*YQ - P32(-k,-j,m)*YU )*WD(l,m)
                           Z=Z+U*(-P31(-j,-k,m)*YI - P32(-j,-k,m)*YQ + P33(-j,-k,m)*YU )*WD(l,m)
                        end do ! m
                     end do ! k
                     SFI(j,l)=X*0.25
                     SFQ(j,l)=Y*0.25
                     SFU(j,l)=Z*0.25
                  end do ! j

!                  write(*,*) '11 : IS,ND,l,I1(l,1),SFI(l,1): ',IS,ND,l,I1(l,1),SFI(l,1)
!PL Multiple scattering at level "l" for upward direction					  

                  do j=1,NN
                     X=0.0
                     Y=0.0
                     Z=0.0
                     do k=1,NN
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        YQ=(Q1(k,l+1)+Q1(k,l))
                        YU=(U1(k,l+1)+U1(k,l))
                        do m=1,NMAX
                           X=X+U*( P11(j,k,m)*YI + P21(k,j,m)*YQ + P31(k,j,m)*YU )*WD(l,m)
                           Y=Y+U*( P21(j,k,m)*YI + P22(j,k,m)*YQ + P32(k,j,m)*YU )*WD(l,m)
                           Z=Z+U*( P31(j,k,m)*YI + P32(j,k,m)*YQ + P33(j,k,m)*YU )*WD(l,m)
                        end do ! m
                     end do ! k
                     do k=-NN,-1
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        YQ=(Q1(k,l+1)+Q1(k,l))
                        YU=(U1(k,l+1)+U1(k,l))
                        do m=1,NMAX
                           X=X+U*( P11(-j,-k,m)*YI + P21( k, j,m)*YQ + P31( k, j,m)*YU )*WD(l,m)
                           Y=Y+U*( P21(-j,-k,m)*YI + P22(-j,-k,m)*YQ + P32( k, j,m)*YU )*WD(l,m)
                           Z=Z+U*(-P31(-j,-k,m)*YI - P32(-j,-k,m)*YQ + P33(-j,-k,m)*YU )*WD(l,m)
                        end do ! m
                     end do ! k
                     SFI(j,l)=X*0.25
                     SFQ(j,l)=Y*0.25
                     SFU(j,l)=Z*0.25
                  end do ! j
               end do ! l
            else
!XH            for scalar case
               SFI(:,:)=0.0
               do l=1,NT-1
                  do j=-NN,-1
                     X=0.0
                     do k=1,NN
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        X=X+U*YI*dot_product(P11(j,k,1:NMAX),WD(l,1:NMAX))
                     end do ! k
                     do k=-NN,-1
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        X=X+U*YI*dot_product(P11(-j,-k,1:NMAX),WD(l,1:NMAX))
                     end do ! k
                     SFI(j,l)=X*0.25
                  end do ! j
                  do j=1,NN
                     X=0.0
                     do k=1,NN
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        X=X+U*YI*dot_product(P11(j,k,1:NMAX),WD(l,1:NMAX))
                     end do ! k
                     do k=-NN,-1
                        U=rg(k)
                        YI=(I1(k,l+1)+I1(k,l))
                        X=X+U*YI*dot_product(P11(-j,-k,1:NMAX),WD(l,1:NMAX))
                     end do ! k
                     SFI(j,l)=X*0.25
                  end do ! j
               end do ! l
            end if
!           integration over tau for upward radiance
!PL         Next order of scattering (for upward direction)  Eq.(31b)
!PL         Eq.(31f) J.Lenoble, M.Herman et al, JQSRT, 2007
            if (keyEl .eq. 4) then
!XH            for vector case
               do j=1,NN
                  I1(j,NT)=2.0*dot_product(rg(1:NN),( I1(-1:-NN:-1,NT)*R11(j,1:NN,IS)  &
                                                     +Q1(-1:-NN:-1,NT)*R12(j,1:NN,IS)  &
                                                     +U1(-1:-NN:-1,NT)*R13(j,1:NN,IS)))
                  Q1(j,NT)=2.0*dot_product(rg(1:NN),( I1(-1:-NN:-1,NT)*R21(j,1:NN,IS)  &
                                                     +Q1(-1:-NN:-1,NT)*R22(j,1:NN,IS)  &
                                                     +U1(-1:-NN:-1,NT)*R23(j,1:NN,IS)))
                  U1(j,NT)=2.0*dot_product(rg(1:NN),( I1(-1:-NN:-1,NT)*R31(j,1:NN,IS)  &
                                                     +Q1(-1:-NN:-1,NT)*R32(j,1:NN,IS)  &
                                                     +U1(-1:-NN:-1,NT)*R33(j,1:NN,IS)))
               end do
            else
!XH            for scalar case
               do j=1,NN
                  I1(j,NT)=2.0*dot_product(rg(1:NN),I1(-1:-NN:-1,NT)*R11(j,1:NN,IS))
               end do
            end if
            do l=NT-1,1,-1
               do j=1,NN
                  C=EXP( ( EXTH(l) - EXTH(l+1) ) / RMU(j) )
                  I1(j,l)=C*I1(j,l+1)+(1.0-C)*SFI(j,l)
                  if (keyEl .eq. 4) then
!XH                  for vector case
                     Q1(j,l)=C*Q1(j,l+1)+(1.0-C)*SFQ(j,l)
                     U1(j,l)=C*U1(j,l+1)+(1.0-C)*SFU(j,l)
                  end if
               end do ! j
            end do ! l
!           integration over tau for downward radiance
!PL         Next order of scattering (for downward direction) Eq.(31a)
            I1(-NN:-1,1)=0.0
            if (keyEl .eq. 4) then
!XH            for vector case
               Q1(-NN:-1,1)=0.0
               U1(-NN:-1,1)=0.0
            end if
            do l=2,NT
               do j=-NN,-1
                  C=EXP( ( EXTH(l) - EXTH(l-1) ) / RMU(j) )
                  I1(j,l)=C*I1(j,l-1)+(1.0-C)*SFI(j,l-1)
                  if (keyEl .eq. 4) then
!XH                  for vector case
                     Q1(j,l)=C*Q1(j,l-1)+(1.0-C)*SFQ(j,l-1)
                     U1(j,l)=C*U1(j,l-1)+(1.0-C)*SFU(j,l-1)
                  end if
               end do ! j
            end do ! l

!PL         Sum of scatering orders (upward and downward)
            if (itest .eq. +1) then
!XH            upward radiance
               IM(1:NN)=IM(1:NN)+I1(1:NN,NAV)
               if (keyEl .eq. 4) then
!XH               for vector case
                  QM(1:NN)=QM(1:NN)+Q1(1:NN,NAV)
                  UM(1:NN)=UM(1:NN)+U1(1:NN,NAV)
               end if
            else
!XH            downward radiance
               IM(-NN:-1)=IM(-NN:-1)+I1(-NN:-1,NAV)
               if (keyEl .eq. 4) then
!XH               for vector case
                  QM(-NN:-1)=QM(-NN:-1)+Q1(-NN:-1,NAV)
                  UM(-NN:-1)=UM(-NN:-1)+U1(-NN:-1,NAV)
               end if
            end if ! itest.eq.+1

!XH         flux calculation only related to 0 Fourier component (IS=0)
            if (IS.eq.0)then
               RS=sum(I1(1:NN,NAV)/I1NAV(1:NN))/NN
               I1NAV(1:NN) = I1(1:NN,NAV)
               IZG(-NN:NN) = IZG(-NN:NN)+I1(-NN:NN,NT)
               if (present(UFX)) then
                  IZX(-NN:NN,1:NT)=IZX(-NN:NN,1:NT)+I1(-NN:NN,1:NT)
               end if
!PL            Aproximate calculation of multiple scattering after 50 orders of scattering
               if (ND.eq.50)then
                  RS1 = RS/(1.0-RS)
                  if (itest .eq. +1) then
                     IM(  1:NN) = IM(  1:NN)+RS1*I1(  1:NN,NAV)
                  else
                     IM(-NN:-1) = IM(-NN:-1)+RS1*I1(-NN:-1,NAV)
                  end if
                  IZG(-NN:NN)=IZG(-NN:NN)+RS1*I1(-NN:NN,NT)
                  if (present(UFX)) then
                     IZX(-NN:NN,1:NT)=IZX(-NN:NN,1:NT)+RS1*I1(-NN:NN,1:NT)
                  end if
                  exit LOOP_scat_order
               end if
            end if

!           if (IS.eq.0) then
!              write(*,*)
!              write(*,*) '1: I1:'
!              write(*,'(10e16.6)') I1
!              write(*,*) 'Q1:'
!              write(*,'(10e16.6)') Q1
!              write(*,*) 'Q1:'
!              write(*,'(10e16.6)') U1
!              write(*,*)
!           end if

!XH         check convergence of radiance calculation
!XH         only need to check I1 for scalar case
            if (itest .eq. +1) then
               z = max(0.0,maxval(abs(I1(  1:NN,NAV))))
               if (keyEl .eq. 4) z = max(z,maxval(abs(Q1(  1:NN,NAV))),maxval(abs(U1(  1:NN,NAV))))
            else
               z = max(0.0,maxval(abs(I1(-NN:-1,NAV))))
               if (keyEl .eq. 4) z = max(z,maxval(abs(Q1(-NN:-1,NAV))),maxval(abs(U1(-NN:-1,NAV))))
            end if ! itest.eq.+1
            if (z.le.eps_err) Exit LOOP_scat_order

         end do LOOP_scat_order
!PL      Here is the end of the Order of scattering cycle

!        if (IS.eq.0) then
!           write(*,*)
!           write(*,*) '2: I1:'
!           write(*,'(10e16.6)') I1
!           write(*,*) 'Q1:'
!           write(*,'(10e16.6)') Q1
!           write(*,*) 'Q1:'
!           write(*,'(10e16.6)') U1
!           write(*,*)
!        end if

!        calculation of total diffuse radiance and check convergence of Fourier series
         xx=0.0
         do iv=1,NBV
            k=idir(iv)
            XI=IM(k+1)+coff(iv)*(IM(k)-IM(k+1))
            xx=max(xx,abs(XI))
            if (IS .eq. 0)then
               SL(iv)=XI
            else
               ang=IS*(fiv(iv)-180.0)*rpi_180
               SL(iv)  =SL(iv)+2.0*cos(ang)*XI
            end if ! IS .eq. 0
         end do ! iv
         if (keyEl .eq. 4) then
!XH         for vector case
            do iv=1,NBV
               k=idir(iv)
               XQ=QM(k+1)+coff(iv)*(QM(k)-QM(k+1))
               XU=UM(k+1)+coff(iv)*(UM(k)-UM(k+1))
               xx=max(xx,abs(XQ),abs(XU))
!               write(6,*) iv,coff(iv),k,IM(k),IM(k+1),XI,  &
!                          '  iv,coff(iv),k,IM(k),IM(k+1),XI'
!               stop
               if (IS .eq. 0)then
                  SQos(iv)=XQ
                  SUos(iv)=0.0
               else
                  ang=IS*(fiv(iv)-180.0)*rpi_180
!PL               Fourier series
!PL               Here the Stoks vector (SL,SQos,SUos) for coupled atmosphere-surface
!PL               system is calculated
                  SQos(iv)=SQos(iv)+2.0*cos(ang)*XQ
                  SUos(iv)=SUos(iv)+2.0*sin(ang)*XU
               end if ! IS .eq. 0
            end do ! iv
         end if
!         write(6,*)'xx ',xx,'  IS ',IS,'  ND ',ND
!		  write(6,*)'XI ',XI
!		  write(6,*)'XQ ',XQ
!		  write(6,*)'XU ',XU
         if (xx .lt. eps_err) exit
      end do ! IS
!      write(6,*)'xx ',xx,'  IS ',IS,'  ND ',ND

      UFG=                         2.0*dot_product(rg(1:NN),rmu(1:NN)*IZG(1:NN))
      DFG=rmus*exp(-EXTH(NT)/rmus)+2.0*dot_product(rg(1:NN),rmu(1:NN)*IZG(-1:-NN:-1))
      if (present(UFX)) then
         do l=1,NT
!XH         integrate over zenith angle to get flux for each level
            UFX(l)=                        2.0*dot_product(rg(1:NN),rmu(1:NN)*IZX( 1: NN   ,l))
            DFX(l)=rmus*exp(-EXTH(l)/rmus)+2.0*dot_product(rg(1:NN),rmu(1:NN)*IZX(-1:-NN:-1,l))
         end do

!XH      deallocate array for flux calculation
         deallocate(IZX)
      end if
!
      return
      end SUBROUTINE rad_transf_MS


! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

 
      Subroutine ANGTURN_PL(tsol,avis,afiv,cos_2etv,sin_2etv)
	  !PL cos_2etv=-cs !!! Because of positive SQ1,SL1
	  !PL sin_2etv=-ss !!! The reason is not clear yet
	  Implicit none
	   Real(8), parameter:: st_min=0.001,sin_eps=0.00000001
	   real, intent(in):: tsol,avis,afiv
	   real, intent(out):: cos_2etv,sin_2etv
       real(8) mu0,muv,del_phi,cos_etv,sin_etv,sin_ti,sin_tv,zz,cos_tet, &
	           sin_tet,tet0v
	   
	   mu0=dcos(dble(tsol*rpi_180))
       muv=dcos(dble(avis*rpi_180))
	   del_phi=afiv*rpi_180-rpi

	   sin_ti=dsqrt(1.d0-mu0*mu0)
	   sin_tv=dsqrt(1.d0-muv*muv)
       
	   zz=muv*mu0+sin_ti*sin_tv*dcos(dble(afiv*rpi_180))
       cos_tet=-zz
       sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
     
       If (dabs(sin_tv) .lt. sin_eps) then
        tet0v=dacos(muv)	  
        tet0v=dabs(tet0v-st_min*rpi_180)
        muv=dcos(tet0v)
        sin_tv=dsqrt(1.d0-muv*muv)
        zz=muv*mu0+sin_ti*sin_tv*dcos(dble(afiv*rpi_180))
        cos_tet=-zz
        sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
       endif
 
       If (dabs(sin_ti) .lt. sin_eps) then
        tet0v=dacos(mu0)	  
        tet0v=dabs(tet0v-st_min*rpi_180)
        mu0=dcos(tet0v)
        sin_ti=dsqrt(1.d0-mu0*mu0)
        zz=muv*mu0+sin_ti*sin_tv*dcos(dble(afiv*rpi_180))
        cos_tet=-zz
        sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
       endif

       if (dabs(sin_tet) .lt. sin_eps) then
        tet0v=dacos(mu0)	  
        tet0v=dabs(tet0v-st_min*rpi_180)
        muv=dcos(tet0v)
        sin_tv=dsqrt(1.d0-muv*muv)
        zz=muv*mu0+sin_ti*sin_tv*dcos(dble(afiv*rpi_180))
        cos_tet=-zz
        sin_tet=dsqrt(1.d0-cos_tet*cos_tet)
       endif


 	   cos_etv=(-mu0-muv*cos_tet)/sin_tet/sin_tv
       sin_etv=dsin(del_phi)*sin_ti/sin_tet
       
	   cos_2etv=-real(cos_etv*cos_etv-sin_etv*sin_etv)
       sin_2etv=-real(2.d0*cos_etv*sin_etv)

	  End Subroutine ANGTURN_PL

  

      FUNCTION ANGTURN(tsol,avis,afiv)
      !pi = 3.141592653
      x0=cos(tsol*rpi_180)
      x1=cos(avis*rpi_180)
      z=cos(afiv*rpi_180)
      x2=x0*x1+z*sqrt(1-x1*x1)*sqrt(1-x0*x0)
      sbeta=(x0-x1*x2)/sqrt(1-x2*x2)/sqrt(1-x1*x1)
      if(sbeta.gt.1.0)sbeta=1.0
      if(sbeta.lt.-1.0)sbeta=-1.0
      ANGTURN=acos(sbeta)-rpi*0.5
!c CORRECTION
      if(afiv.gt.180.)ANGTURN=-ANGTURN
      if(afiv.lt.0.1.or.afiv.gt.359.9		&
      .or.abs(afiv-180.0).lt.0.1)ANGTURN=rpi*0.5
      if(avis.lt.0.1.or.avis.gt.179.9)ANGTURN=-rpi*afiv/180.0+rpi*0.5
      RETURN
      END FUNCTION ANGTURN

!	  FUNCTION ANGTURN(tsol,avis,afiv)

!      !pi = 3.141592653
!	  pi=3.141592653589793

!	  write(*,'(2f25.15)') pi/180.,pi
!	  write(*,'(2f25.15)')rpi_180,rpi

!      x0=cos(tsol*pi/180.)
!      x1=cos(avis*pi/180.)
!      z=cos(afiv*pi/180.)
!      x2=x0*x1+z*sqrt(1-x1*x1)*sqrt(1-x0*x0)
!      sbeta=(x0-x1*x2)/sqrt(1-x2*x2)/sqrt(1-x1*x1)
!      if(sbeta.gt.1.0)sbeta=1.0
!      if(sbeta.lt.-1.0)sbeta=-1.0
!      ANGTURN=acos(sbeta)-pi*0.5!/2.
!!c CORRECTION
!      if(afiv.gt.180.)ANGTURN=-ANGTURN
!      if(afiv.lt.0.1.or.afiv.gt.359.9		&
!      .or.abs(afiv-180.0).lt.0.1)ANGTURN=pi*0.5!/2.0
!      if(avis.lt.0.1.or.avis.gt.179.9)ANGTURN=-pi*afiv/180.0+pi*0.5!/2.
!      RETURN
!      END FUNCTION ANGTURN


! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine gauss(MM,AMU,PMU,NPAR)
!C     ORDRE DE LA QUADRATURE N=2*MM-2 
      !use mod_par_OS, only : IIX

!CD      PARAMETER (IX=600)
!CD      PARAMETER(NG0=91)
      !PARAMETER (IX=2*IIX)
      IMPLICIT Real(8) (A-H,O-Z)
!tl      DIMENSION Z(IX),PA(IX),W(IX),R(IX),AMU(-NG0:NG0),PMU(-NG0:NG0)
      DIMENSION Z(2*MM),PA(2*MM),W(2*MM),R(2*MM),AMU(-NPAR:NPAR),PMU(-NPAR:NPAR)
      TOL = 1.0D-15
      PI  = 3.141592653589793D+00
      PI2 = PI*PI
      N=2*MM-2
      XL=-1.
      XU=+1.
      AA = 2.0D+00/PI2
      AB = -62.0D+00/(3.0D+00*PI2*PI2)
      AC = 15116.0D+00/(15.0D+00*PI2*PI2*PI2)
      AD = -12554474.D+00/(105.0D+00*PI2*PI2*PI2*PI2)
      PA(1) = 1.D0
      EN = N
      NP1 = N+1
      U= 1.0D+00-(2.0D+00/PI)**2
      D = 1.0D+00/DSQRT((EN+0.5D+00)**2+U*0.25D+00)

      DO I= 1,N
      SM = I
      AZ = 4.0D+00*SM-1.0D+00
      AE = AA/AZ
      AF = AB/AZ**3
      AG = AC/AZ**5
      AH = AD/AZ**7
      Z(I) = 0.25D+00*PI*(AZ+AE+AF+AG+AH)
      ENDDO ! I
	  
      DO K = 1,N
      X = COS(Z(K)*D)
      XDD=1.d0
      do while(XDD .gt. 0.d0) ! tl
         PA(2) = X
         DO NN = 3,NP1
         ENN = NN-1
         PA(NN) =		&
         ((2.0D+00*ENN-1.0D+00)*X*PA(NN-1)-(ENN-1.0D+00)*PA(NN-2))/ENN
         ENDDO ! NN
         PNP = EN*(PA(N)-X*PA(NP1))/(1.0D+00-X*X)
         XI = X-PA(NP1)/PNP
         XD = ABS(XI-X)
         XDD = XD-TOL
         X = XI
      enddo ! while(XDD .gt. 0.d0)
      R(K) = X
      W(K) = 2.0D+00*(1.0D+00-X*X)/(EN*PA(N))**2
      ENDDO ! K
	  
      AP = (XU-XL)/2.D0
      BP = (XU+XL)/2.D0
      do I=1,MM-1
      K=MM-I
      AMU(K)=BP+AP*R(I)
      PMU(K)=AP*W(I)
      AMU(-K)=-AMU(K)
      PMU(-K)=PMU(K)
      enddo
      AMU(-MM)=-1.
      AMU(MM)=1.
      AMU(0)=0.
      PMU(0)=0.
      PMU(-MM)=0.
      PMU(MM)=0.
      RETURN
      END subroutine gauss

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

!PL   new Subroutine for BRDF Fourier expansion
!PL   Pij(NN0,NN0,0:2*NG0-2), but just 0:NN expansion are used !!!

      Subroutine developpe_ocean_land ( n_el,land_percent,NN,M_surf,rmu,WAVE,       &
                                        iBRDF,iBPDF,n_par,n_par_i,x_vect,           &
                                        iBRM_water,n_par_water,x_vect_water,        &
                                        P11,P12,P13,P21,P22,P23,P31,P32,P33)
      use mod_par_OS, only : NN0,NG0,NF0
      Integer,                               intent(in)  :: n_el,iBRDF,iBPDF,iBRM_water,NN,M_surf
      integer,                               intent(in)  :: n_par,n_par_i,n_par_water
      Real,                                  intent(in)  :: land_percent,WAVE
      Real(8),                               intent(in)  :: rmu(-NN0:NN0)
      Real(8),                               intent(in)  :: x_vect(1:n_par),x_vect_water(1:n_par_water)
      real,dimension(NN0,NN0,0:NF0),         intent(out) :: P11,P12,P13,P21,P22,P23,P31,P32,P33
!XH   R_exp requires lots of memory
!XH   introduce n_el to determine the dimension of R_exp we need
      Real(8) :: R_exp(n_el,n_el,NN,NN,0:M_surf)
	  Integer :: j,k,IS
	  Double complex Complex_land,Complex_water
      		 
      iBRM1=iBRM_water		 
      If (iBRM_water .eq. 1) iBRM1=0 !Switch to isotropic model
          call BRM_Fexp_OSH(n_el,land_percent,M_surf,NN,rmu(1:NN),dble(WAVE),       &
                            iBRDF,iBPDF,n_par,n_par_i,x_vect(1:n_par),              &
                            iBRM1,n_par_water,x_vect_water(1:n_par_water),          &
                            R_exp)
!PL   R_exp(mu0,muv) but P(muv,mu0), then jk -> kj
      do k=1,NN
         do j=1,NN
!PL         Since in PL BRDFs (R_exp) the second index corresponds to viewing
!PL         direction whereas the first one corresponds to incident (solar)
!PL         direction, then:
            P11(j,k,0:M_surf)=R_exp(1,1,k,j,0:M_surf)*rmu(k)
            if (n_el .ge. 3) then
!XH            for vector case only
               P21(j,k,0:M_surf)=R_exp(2,1,k,j,0:M_surf)*rmu(k)
               P31(j,k,0:M_surf)=R_exp(3,1,k,j,0:M_surf)*rmu(k)
               P12(j,k,0:M_surf)=R_exp(1,2,k,j,0:M_surf)*rmu(k)
               P22(j,k,0:M_surf)=R_exp(2,2,k,j,0:M_surf)*rmu(k)
               P32(j,k,0:M_surf)=R_exp(3,2,k,j,0:M_surf)*rmu(k)
               P13(j,k,0:M_surf)=R_exp(1,3,k,j,0:M_surf)*rmu(k)
               P23(j,k,0:M_surf)=R_exp(2,3,k,j,0:M_surf)*rmu(k)
               P33(j,k,0:M_surf)=R_exp(3,3,k,j,0:M_surf)*rmu(k)
            end if
         end do !j
      end do !k
!
      return
      end subroutine developpe_ocean_land

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      !!PL Fourier expansion coeficient calculation for the elelments of
	  !!PL phase matrix using generalized Legendre functions

      subroutine noyaux(NG,xmu,NN,rmu,alp,bet,gam,zet,P11,P12,P13,		&
      P21,P22,P23,P31,P32,P33)
      use mod_par_OS, only : NN0,NG0
      Real(8) rmu,xmu
!CD      PARAMETER(NN0=85,NG0=91)
      dimension rmu(-NN0:NN0),xmu(-NG0:NG0), 		& 
      alp(0:2*NG0-2),bet(0:2*NG0-2),gam(0:2*NG0-2),zet(0:2*NG0-2),		&
      P11(NN0,NN0,0:2*NG0-2),P12(NN0,NN0,0:2*NG0-2),		&
      P13(NN0,NN0,0:2*NG0-2),		&
      P21(NN0,NN0,0:2*NG0-2),P22(NN0,NN0,0:2*NG0-2),		&
      P23(NN0,NN0,0:2*NG0-2),		&
      P31(NN0,NN0,0:2*NG0-2),P32(NN0,NN0,0:2*NG0-2),		&
      P33(NN0,NN0,0:2*NG0-2),		&
      PSL(-1:2*NG0,-NN0:NN0),RSL(-1:2*NG0,-NN0:NN0),		&
      TSL(-1:2*NG0,-NN0:NN0)
      do 200 m=0,2*NG-2
       call legendre(m,rmu,NG,NN,PSL,TSL,RSL)
       do 201 j=1,NN
        do 202 k=1,NN
         P11(j,k,m)=0.0
         P12(j,k,m)=0.0
         P13(j,k,m)=0.0
         P21(j,k,m)=0.0
         P22(j,k,m)=0.0
         P23(j,k,m)=0.0
         P31(j,k,m)=0.0
         P32(j,k,m)=0.0
         P33(j,k,m)=0.0
         DO 203 L=m,2*NG-2
          TT=TSL(L,j)*TSL(L,-k)
          TR=TSL(L,j)*RSL(L,-k)
          RR=RSL(L,j)*RSL(L,-k)
          RT=RSL(L,j)*TSL(L,-k)
          P11(j,k,m)=P11(j,k,m)+bet(L)*PSL(L,j)*PSL(L,-k)
          P21(j,k,m)=P21(j,k,m)+gam(L)*RSL(L,j)*PSL(L,-k)
          P31(j,k,m)=P31(j,k,m)-gam(L)*TSL(L,j)*PSL(L,-k)
          P12(j,k,m)=P12(j,k,m)+gam(L)*PSL(L,j)*RSL(L,-k)
          P22(j,k,m)=P22(j,k,m)+alp(L)*RR+zet(L)*TT
          P32(j,k,m)=P32(j,k,m)-alp(L)*TR-zet(L)*RT
          P13(j,k,m)=P13(j,k,m)-gam(L)*PSL(L,j)*TSL(L,-k)
          P23(j,k,m)=P23(j,k,m)-alp(L)*RT-zet(L)*TR
          P33(j,k,m)=P33(j,k,m)+alp(L)*TT+zet(L)*RR
  203 continue
  202 continue
  201 continue
  200 continue
      return
      end subroutine noyaux

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	  !! PL Fresnel matrix calculation

      subroutine fresnel(NG,xmu,xind,t11,t12,t22,t33)
      use mod_par_OS, only : NG0
      Real(8) xmu
!CD      PARAMETER(NG0=91)
      dimension		&
      t11(-NG0:NG0),t12(-NG0:NG0),t22(-NG0:NG0),t33(-NG0:NG0),		&
      xmu(-NG0:NG0)
      !pi = 3.141592653

      do j=-NG,NG
       av=acos(xmu(j))*r180_pi
       ai=(180.-av)*0.5
       xx=cos(ai*rpi_180)
       yy=sqrt(xind*xind-1.0+xx*xx)
       zz=xx*xind*xind
       rl=(zz-yy)/(zz+yy)
       rr=(xx-yy)/(xx+yy)
       t11(j)=0.5*(rl*rl+rr*rr)
       t12(j)=0.5*(rl*rl-rr*rr)
       t22(j)=t11(j)
       t33(j)=rl*rr
      enddo
      return
      end subroutine fresnel

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!PL   Expansion coeficientcal calculation
      subroutine betal(keyEl,LL,NG0,zmu,zp,f11,f12,f22,f33,     &  ! in
                       alpa,beta,gama,zeta           )             ! out
!XH   optional parameters for vector case
      integer,                                intent(in)  :: keyEl, LL, NG0
      real*8, dimension(-NG0:NG0),            intent(in)  :: zmu,zp
      real,   dimension(-NG0:NG0),            intent(in)  :: f11,f12,f22,f33
!
      real,   dimension(0:2*NG0-2),           intent(out) :: beta,alpa,gama,zeta
!
      integer :: k,j
      real,   dimension(-1:2*NG0)  :: pl,pol,ppl,pml
      real,   dimension(0:2*NG0-2) :: betap,betam

!XH   it can be simplified to series expansion in Legendre polynomials for scalar case
      pl(-1)=0.                                                                 
      pl(0) =1.
      beta(0:LL+LL-2)=0.0
      do j=-LL,LL
         xx=zmu(j)
         do k=0,LL+LL-2
            pl(k+1)=((2*k+1.)*xx*pl(k)-k*pl(k-1))/(k+1.)
            beta(k)=beta(k)+zp(j)*pl(k)*f11(j)*(k+0.5)
         end do ! k
      end do ! j
      if (keyEl .eq. 4) then
         gama(0:LL+LL-2) =0.
         alpa(0:LL+LL-2) =0.
         zeta(0:LL+LL-2) =0.
         betap(0:LL+LL-2)=0.
         betam(0:LL+LL-2)=0.
         pol(0)=0.0
         pol(1)=0.0
         ppl(0)=0.0
         ppl(1)=0.0
         pml(0)=0.0
         pml(1)=0.0
         do j=-LL,LL
            xx=zmu(j)
            pol(2)=3.0*(1.0-xx*xx)*0.5/sqrt(6.0)
            ppl(2)=(1.+xx)*(1.+xx)*0.25
            pml(2)=(1.-xx)*(1.-xx)*0.25
            do k=0,2*LL-2
               if (k .gt. 1)then
                  dd=(2*k+1.)/sqrt((k+3.)*(k-1.))
                  ee=sqrt((k+2.)*(k-2.))/(2*k+1.)
                  cc=(k+1.)*(k+2.)*(k-2.)/k/(k+3.)/(k-1.)
                  bb=(2*k+1.)/k/(k+3.)/(k-1)
                  pol(k+1)=dd*(xx*pol(k)-ee*pol(k-1))
                  ppl(k+1)=bb*(k*(k+1.)*xx-4.)*ppl(k)-cc*ppl(k-1)
                  pml(k+1)=bb*(k*(k+1.)*xx+4.)*pml(k)-cc*pml(k-1)
               end if ! k .gt. 1
               gama(k)=gama(k)+zp(j)*pol(k)*f12(j)*(k+0.5)
               betap(k)=betap(k)+zp(j)*ppl(k)*(f22(j)+f33(j))*(k+0.5)
               betam(k)=betam(k)+zp(j)*pml(k)*(f22(j)-f33(j))*(k+0.5)
               zeta(k)=(betap(k)-betam(k))*0.5
               alpa(k)=(betap(k)+betam(k))*0.5
            end do ! k
         end do ! j
      end if
!
      return
      end subroutine betal

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

!PL   Legendre function calculation
      subroutine legendre(IS,rmu,JJ,LL,    &
                          PSL,TSL,RSL )
!XH   only need PSL for scalar case
!XH   JJ : number of terms in phase matrix expansion
!XH   LL : number of terms in spherical surface integration
      use mod_par_OS, only : NN0,NG0
      integer,                                    intent(in)  :: IS, JJ, LL
      real(8),dimension(-NN0:NN0),                intent(in)  :: rmu
      real,dimension(0:2*NG0,-NN0:NN0),          intent(out) :: PSL
      real,dimension(0:2*NG0,-NN0:NN0),optional, intent(out) :: RSL, TSL
!      real,dimension(-1:2*NG0,-NN0:NN0),          intent(out) :: PSL
!      real,dimension(-1:2*NG0,-NN0:NN0),optional, intent(out) :: RSL, TSL
!
      integer :: I, J, K, LP, LM
      real :: X,A,B,C,D,E,F,XX,YY
!
      select case(IS)
      case(0)
         do J=-LL,LL
            C=RMU(J)
            PSL(0,J)=1.
            PSL(1,J)=C
            PSL(2,J)=1.5*C*C-0.5
            if (present(TSL)) then
               TSL(1,J)=0.
               TSL(2,J)=0.
               RSL(1,J)=0.
               RSL(2,J)=1.5/sqrt(6.0)*(1.-C*C)
            end if
         end do ! J
      case(1)
         do J=-LL,LL
            C=RMU(J)
            X=1.-C*C
            PSL(0,J)=0.
            PSL(1,J)=SQRT(X*0.5)
            PSL(2,J)=C*sqrt(3.0)*PSL(1,J)
            if (present(TSL)) then
               TSL(1,J)=0.
               TSL(2,J)=-SQRT(X)*0.5
               RSL(1,J)=0.
               RSL(2,J)=C*TSL(2,J)
            end if
         end do ! J
      case default
         A=1.
         do I=1,IS  
            A=A*SQRT(1.0*(I+IS)/I)*0.5
         end do
         B=A*SQRT(IS/(IS+1.))*SQRT((IS-1.)/(IS+2.))
         do J=-LL,LL
            C=RMU(J)
            XX=1.-C*C
            YY=IS*0.5-1.
            PSL(IS-1,J)=0.
            PSL(IS,J)=A*XX**(IS*0.5)
            if (present(TSL)) then
               X = B*XX**YY
               TSL(IS-1,J)=0.
               TSL(IS,J)=2.*C*X
               RSL(IS-1,J)=0.
               RSL(IS,J)=(1.+C*C)*X
            end if
         end do
      end select ! is

      K=MAX(2,IS)
      if (K .ne. 2*JJ-2) then
         do L=K,2*JJ-1
            LP=L+1
            LM=L-1
            A=(2*L+1.)/SQRT((L+IS+1.)*(L-IS+1.))
            B=SQRT(FLOAT((L+IS)*(L-IS)))/(2.*L+1.)
            PSL(LP,-LL:LL)=A*(RMU(-LL:LL)*PSL(L,-LL:LL)-B*PSL(LM,-LL:LL))
            if (present(TSL)) then
               D=(L+1.)*(2*L+1.)/SQRT((L+3.)*(L-1.)*(L+IS+1.)*(L-IS+1.))
               E=SQRT((L+2.)*(L-2.)*(L+IS)*(L-IS))/(L*(2.*L+1))
               F=2.*IS/(L*(L+1.))
               RSL(LP,-LL:LL)=D*(RMU(-LL:LL)*RSL(L,-LL:LL)-F*TSL(L,-LL:LL)-E*RSL(LM,-LL:LL))
               TSL(LP,-LL:LL)=D*(RMU(-LL:LL)*TSL(L,-LL:LL)-F*RSL(L,-LL:LL)-E*TSL(LM,-LL:LL))
            end if
         end do
      end if ! K .ne. 2*JJ-2
	     
      return     
      end subroutine legendre

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine profiles_new (	 IPRI_additional_info,              &
                                 TEXT, EXT, SSA,                    &
                                 HGR_km, HOBS_km, HMAX_atm_km, NT1, &
                                 DISCRVD,                           &
                                 NT, NAV, EXTH, WD, HLV             & ! OUT
                              )

      use mod_par_OS, only : KNT, NMM, NMG, KVERT_WD
      use mod_intrpl_linear
      use mod_vertical_distr_derived_type

      implicit none
!	------------------------------------------------------------------------------------------------------
      logical,                          intent(in)  ::  IPRI_additional_info
      integer,                          intent(in)	::  NT1
      real,                             intent(in)  ::  TEXT, HGR_km
      real,                             intent(in)  ::  HOBS_km, HMAX_atm_km
      type(discret_vertical_distribution), intent(in) :: DISCRVD
      real,dimension(NMM+NMG),          intent(in)  ::  EXT,SSA
      integer,                          intent(out) ::  NT
      integer,                          intent(out) ::  NAV
      real,dimension(KNT),              intent(out) ::  EXTH
      real,dimension(KNT-1,NMM+NMG),    intent(out) ::  WD
      real,dimension(KNT),optional,     intent(out) ::  HLV
!	------------------------------------------------------------------------------------------------------
      real, dimension(KVERT_WD)          ::  H
      real, dimension(NMM+NMG,KVERT_WD)  ::  tauH
      real, dimension(KVERT_WD,NMM+NMG)  ::  PROF
      integer :: natm
!	------------------------------------------------------------------------------------------------------
      integer                     ::  i,i1
!	------------------------------------------------------------------------------------------------------  
      real, parameter             ::  tiny = 1e-6
      integer                     ::  NH, ip, ip1, NTEMP, NTEMP1
      real                        ::  dtau, dtaup, &
                                      TEXT_HP, TEXT_HP0
										
      real,dimension(KVERT_WD)        ::  H_temp
      real,dimension(NMM+NMG,KVERT_WD)::  tauH_temp
      real,dimension(KNT)             ::  EXTH_temp
      real,dimension(KNT-1,NMM+NMG)   ::  WD_temp
      real,dimension(KNT)             ::  HL, HL_temp
      real,dimension(NMM+NMG)         ::  xnorm, tauHP, EXT_HP, EXT_HP0
      integer                         ::  iplane
      real                            ::  xnorm0
      real                            ::  DH

! iplane = 1  - satellite data
!        = 2  - plane  meas 
!        = 3  - ground based meas 

!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
!  a1 and b1 - coefficients for straight y(h) between grid heghts
!  y(i-1)=a1*h(i-1)+b1                    
!  y(i)  =a1*h(i)+b1 
!  a1=(y(i)-y(i-1))/(h(i)-h(i-1))
!  b1=y(i-1)-a1*h(i-1) 
!
!  NT - number of layers ; L - layer index
!  dh=h(i)-h(i-1) < 0.0 
!  dtau=TEXT/NT 
!  diff=dtau-tauS(i-1) ; tauS(i-1)=(y(h(L-1))+y(h(i-1))*(h(L-1)-h(i-1))
!  
!  0.5*(y(i-1)+y(i))*(h(i-1)-h(i))=diff 
!  0.5*(a1*h(i-1)+b+a1*(h(i-1)+dh)+b))*(-dh)=diff
!  0.5*a1*dh**2+(a1*h(i-1)+b)+diff=0.0
!
!  A = 0.5*a1
!  B = a1*h(i-1)+b1 ; B = y(h(i-1))
!  C = diff
!
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
      natm = DISCRVD%natm ! number of aerosol+clouds+molecular+gas
      WD(:,:) = 0.0
      EXTH(:) = 0.0      
      HL(:)   = 0.0
         
      if(HOBS_km .gt. HMAX_atm_km) then
         iplane=1 ! satellite 
      else if(HOBS_km .gt. HGR_km .and. HOBS_km .lt. HMAX_atm_km) then
         iplane=2 ! plane     
      else if(HOBS_km .eq. HGR_km) then
         iplane=3 ! ground based
      else if(HGR_km .ge. HMAX_atm_km) then
         write(tmp_message,'(2(a,f8.4))') &
         'HGR_km =',HGR_km,' .ge.  HMAX_atm_km =',HMAX_atm_km
         G_ERROR(trim(tmp_message))
      endif
        
      NH=KVERT_WD
!	------------------------------------------------------------------------------------------------------
! calculate NT - number of layers with tau<=dtau
      dtau = 0.02
      NT = int(TEXT/dtau)
      if ((MOD(TEXT,dtau)) .gt. 0.) NT = NT+1
      if (NT .gt. NT1-1) NT=NT1-1
      dtau=TEXT/real(NT)
      !write(*,*) 'new: NT=',NT,'  NT1=',NT1,'  dtau =',dtau,'  tetot =',text

      NH = DISCRVD%nh
      H(1:NH) = DISCRVD%h_km(1:NH)
      do i1=1,natm
        tauH(i1,1:NH) = DISCRVD%val(1:NH,i1)*EXT(i1)/DISCRVD%norm(i1)
      end do ! i1

      IF (iplane .ne. 2) THEN
         TEXT_HP0 = 0.0
         EXT_HP0(:) = 0.0
         !HL_temp(:) = 0.0
         call get_WD (     IPRI_additional_info,&
                           natm,NT,KVERT_WD,      & ! IN
                           TEXT,EXT,SSA,dtau,   &
                           TEXT_HP0,EXT_HP0,    &
                           NH,H,tauH,           &
                           EXTH,WD,HL_temp      & ! OUT
                     )
         if ( stop_report%status ) return

         HL(1:NT+1)=HL_temp(1:NT+1)

      ELSE ! iplane=2, plane observations

         do i1=1,natm	   
         tauHP(i1) = linear(H,tauH(i1,:),NH,HOBS_km)
         enddo ! i1	
		 	
         EXT_HP(:)=0.0
         ip=1		

         do while(HOBS_km .le. H(ip))		 		 
         ip=ip+1
         EXT_HP(1:natm)=EXT_HP(1:natm)+0.5*(tauH(1:natm,ip-1)+tauH(1:natm,ip))*(H(ip-1)-H(ip))
         enddo ! while(HP_km .lt. H(ip))

         EXT_HP(1:natm)=EXT_HP(1:natm)-0.5*(tauH(1:natm,ip-1)+tauH(1:natm,ip))*(H(ip-1)-H(ip))+  &
		                             0.5*(tauH(1:natm,ip-1)+tauHP(1:natm))*(H(ip-1)-HOBS_km)									 
         TEXT_HP=SUM(EXT_HP(1:natm))
         ip1=ip
         !write(*,*) 'ip=',ip,'  HOBS_km=',HOBS_km,'  H(ip)=',H(ip)	   

! above plane layer

         NTEMP = 0
         NTEMP1 = 0
         NTEMP = int(TEXT_HP/dtau )
         !if( MOD(TEXT_HP,dtau) .gt. 0. ) NTEMP=NTEMP+1
         if(NTEMP .gt. NT1-1) NTEMP = NT1-1
         NTEMP1 = NTEMP
         dtaup = TEXT_HP/real(NTEMP)
         !write(*,*) 'dtaup=',dtaup,'  dtau=',dtau,'  NTEMP_above=',NTEMP
         !write(*,*) 'ip=',ip,'  HOBS_km=',HOBS_km,'  TEXT_HP=',TEXT_HP

         H_temp(1:ip-1) = H(1:ip-1)
         H_temp(ip) = HOBS_km
         tauH_temp(1:natm,1:ip-1) = tauH(1:natm,1:ip-1)
         tauH_temp(1:natm,ip) = tauHP(1:natm)
         !EXTH_temp(:) = 0.0
         !WD_temp(:,:) = 0.0
         !HL_temp(:) = 0.0
         TEXT_HP0 = 0.0
         EXT_HP0(1:natm) = 0.0
         call get_WD (     IPRI_additional_info,      &
                           natm,NTEMP,KVERT_WD,         & ! IN
                           TEXT_HP,EXT_HP,SSA,dtaup,  &    
                           TEXT_HP0,EXT_HP0,          &
                           ip,H_temp,tauH_temp,       &
                           EXTH_temp,WD_temp,HL_temp  & ! OUT
                     )
         if ( stop_report%status ) return

         EXTH(1:NTEMP+1) = EXTH_temp(1:NTEMP+1)
         WD(1:NTEMP,1:natm) = WD_temp(1:NTEMP,1:natm)
         HL(1:NTEMP+1) = HL_temp(1:NTEMP+1)
          
! under plane layer

         NTEMP = int((TEXT-TEXT_HP)/dtau)
         if((MOD((TEXT-TEXT_HP),dtau)) .gt. 0.) NTEMP = NTEMP+1
         if(NTEMP1+NTEMP .gt. NT1-1) NTEMP = (NT1-1)-NTEMP1
         !if(NTEMP .le. 0) then
         !write(*,*) 'NTEMP_under_plane=',NTEMP,' .le. 0'
         !stop 'stop in profiles_new'			   
         !endif ! NTEMP .le. 0

         dtaup = (TEXT-TEXT_HP)/real(NTEMP)
	   
         if(NTEMP1+NTEMP+1 .ne. NT+1) then
            write(tmp_message,'(2(a,i0))') 'NTEMP+1 = ',NTEMP1+NTEMP+1,' .NE. NT+1 = ',NT+1
            G_ERROR(trim(tmp_message))
         endif

         H_temp(1) = HOBS_km
         H_temp(2:NH-ip+2) = H(ip:NH)
         tauH_temp(1:natm,1) = tauHP(1:natm)
         tauH_temp(1:natm,2:NH-ip+2) = tauH(1:natm,ip:NH)
         !EXTH_temp(:)=0.0
         !WD_temp(:,:)=0.0
         !HL_temp(:)=0.0
         ip = NH-ip+2
         TEXT_HP0 = TEXT_HP
         EXT_HP0(1:natm) = EXT_HP(1:natm)
         call get_WD (     IPRI_additional_info,         &
                           natm,NTEMP,KVERT_WD,            & ! IN
                           TEXT,EXT,SSA,dtaup,           &
                           TEXT_HP0,EXT_HP0,             &
                           ip,H_temp,tauH_temp,          &
                           EXTH_temp,WD_temp,HL_temp     & ! OUT
                     )
         if ( stop_report%status ) return

         do i=1,NTEMP+1
         EXTH_temp(i) = EXTH_temp(i)+TEXT_HP
         enddo ! i	
		
         EXTH(NTEMP1+1:NT+1) = EXTH_temp(1:NTEMP+1)
         WD(NTEMP1+1:NT,1:natm) = WD_temp(1:NTEMP,1:natm)
         HL(NTEMP1+1:NT+1) = HL_temp(1:NTEMP+1)
	   	   
      ENDIF ! iplane.ne.2
		
      !do i=1,NT+1
		  !write(*,'(a,i5,f15.6,10e14.6)') 'i,HL(i),EXTH(i),TEXT: ',  &
      !								 i,HL(i),EXTH(i),TEXT									
		  !enddo ! i	

      !do i=1,NT
		  !write(*,*) i,WD(i,1:natm),'  - i,WD(i,1:natm)'
		  !enddo ! i		

      NT = NT+1

      select case(iplane)
      case(1)
         NAV = 1         ! satellite
      case(2)        
         NAV = NTEMP1+1  ! plane
      case(3)
         NAV = NT        ! ground based
      case default 
         write(tmp_message,'(a,i0,a)') 'iplane = ',iplane,'  - unknown value'
         G_ERROR(trim(tmp_message))
      end select ! iplane
!XH   optional output of altitude for each level
      if (present(HLV)) HLV = HL
!
      return
      end subroutine profiles_new

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
subroutine root(A,B,C,x1,x2)

!quadratic equation solutions 	  
 
      implicit none 
!  ------------------------------------------------------------------------------------------------------	  
      real, intent(in)    ::  A, B, C
      real, intent(out)   ::  x1, x2
!  ------------------------------------------------------------------------------------------------------	  	  	  
      real                ::  D
      real                ::  sgn
      real                ::  tiny
      logical             ::  lind
!  ------------------------------------------------------------------------------------------------------
! A careful floating point computer implementation combines several strategies to produce a robust result. 
! Assuming the discriminant, b2 − 4ac, is positive and b is nonzero, the computation would be as follows:
! x1 = (-b - sgn( b ) * sqrt( b^2 - 4ac) ) / (2*a), x2 = 2c / (-b -sgn( b ) * sqrt( b^2 - 4ac ) = c / (a*x1). 
! Here sgn denotes the sign function, where sgn(b) is 1 if b is positive and −1 if b is negative. 
! This avoids cancellation problems between b and the square root of the discriminant by ensuring that 
! only numbers of the same sign are added.
!  ------------------------------------------------------------------------------------------------------
      x1 = 0.0
      x2 = 0.0
      
      tiny = 1e-30
      lind=.false.
      
      D = B*B-4.0*A*C     

!      if(abs(D) .le. tiny) then
      if(D .eq. 0.0) then
         x1 = -0.5*B/A 
         x2 = -0.5*B/A   
         lind=.true.
      elseif(D .gt. 0.0) then      
        if(B .ne. 0.0) then
! Floating-point implementation
          if(B .lt. 0.0) then
            sgn = -1.0
          else
            sgn = +1.0
          endif ! B .lt. 0.0            
          x1 = 0.5*(-B-sgn*sqrt(D))/A
          x2 = C/(A*x1)				   
          lind=.true.
        else
          x1 = 0.5*(-B+sqrt(D))/A
          x2 = 0.5*(-B-sqrt(D))/A	
          lind=.true.
        endif ! B .ne. 0.0			   	  
!      elseif(D .lt. -tiny) then			 
      elseif(D .lt. 0.0) then			 
        write(tmp_message,'(a,es11.4,a)') 'D =',D,'  Determinant < 0.0'
        G_ERROR(trim(tmp_message))
      endif ! abs(D) .le. tiny

      if(.not. lind) then
        write(tmp_message,'(a)') 'roots are not defined'
        G_ERROR(trim(tmp_message))
      endif ! .not. lind
      
return
end subroutine root 

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss	

      subroutine get_WD (IPRI_additional_info,  &
                         natm,NT,max_NH,        & ! IN
                         TEXT,EXT,SSA,dtau,     &
                         TEXT_HP0,EXT_HP0,      &
                         NH,H,TH,               &
                         EXTH,WD,HL             & ! OUT
                        )


      use mod_par_OS, only : KNT,NMM,NMG
!	------------------------------------------------------------------------------------------------------
      implicit none 
!	------------------------------------------------------------------------------------------------------
      integer,                        intent(in)  ::  natm, NT, NH, max_NH
      logical,                        intent(in)  ::  IPRI_additional_info
      real,                           intent(in)  ::  TEXT, dtau
      real,dimension(NMM+NMG),        intent(in)  ::  EXT, SSA
      real, dimension(max_NH),        intent(in)  ::  H
      real, dimension(NMM+NMG,max_NH),intent(in)  ::  TH
      real, dimension(NMM+NMG),       intent(in)  ::  EXT_HP0 
      real,                           intent(in)  ::  TEXT_HP0

      real,dimension(KNT),            intent(out) ::  EXTH
      real,dimension(KNT-1,NMM+NMG),  intent(out) ::  WD
      real,dimension(KNT),            intent(out) ::  HL
!	------------------------------------------------------------------------------------------------------
!     natm         - number of component (aerosol+molecular+gas)
!     NH         - number of grid heights for profiles
!     H          - grid heights
!     TH         - PROF*EXT/norm - calculated at grid heights (see subroutine discr_profile)
!     TEXT       - total EXT
!     EXT        - EXT for each component
!     SSA        - single scattering albedo for each component
!     dtau       - total EXT of layer
!     TEXT_HP0   - total EXT at plane heght
!     EXT_HP0    - EXT for each component at plane heght
!
!     TEXT_HP0=0.0 and EXT_HP0=0.0 for satellite and ground observations
!
!     tau(0)     - tau for sum of aerosol component
!     tau(1:NMM) - tau for each aerosol component

      real, dimension(0:NMM+NMG,KNT)   ::  TL
      real, dimension(0:NMM+NMG)       ::  TH1, TH2, a1, b1, &
                                           TS1, TS2, TS  
      real, dimension(1:NMM+NMG,KNT)   ::  EXTH1		
      real                             ::  x1, x2, DH1, DH2,       &
                                           A, B, C, TTH,           &
                                           H1, H2, diff, diff_test  

      integer                          ::  i, i1, IL
      real                             ::  tiny, tiny1
      integer :: istop
      istop = 0
!	------------------------------------------------------------------------------------------------------  
!	------------------------------------------------------------------------------------------------------
!     a1 and b1 - coefficients for straight y(h) between grid heights
!     y(i-1)=  a1*h(i-1)+b1
!     y(i)  =  a1*h(i)+b1
!     a1    =  (y(i)-y(i-1))/(h(i)-h(i-1))
!     b1    =  y(i-1)-a1*h(i-1)
!
!     NT - number of layers ; L - layer index
!     dh     = h(i)-h(i-1) < 0.0
!     dtau   = TEXT/NT
!     diff   = dtau-TS(i-1)
!     TS(i-1)= (y(h(L-1))+y(h(i-1))*(h(L-1)-h(i-1))
!  
!     0.5*(y(i-1)+y(i))*(h(i-1)-h(i))=diff
!     0.5*(a1*h(i-1)+b+a1*(h(i-1)+dh)+b))*(-dh)=diff
!     0.5*a1*dh**2+(a1*h(i-1)+b)+diff=0.0
!
!     A = 0.5*a1
!     B = a1*h(i-1)+b1 ; B = y(h(i-1))
!     C = diff
!	------------------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------------------
      tiny  = 1e-4
      tiny1 = 1e-30
!	------------------------------------------------------------------------------------------------------
      EXTH(:)=0.0	
      WD(:,:)=0.0				
				
!     Calculate heights, EXTH (at levels) and WD (for layers) with tau<=dtau

      !write(*,*) natm, NT, KNT, max_NH, dtau,' - natm, NT, KNT, max_NH, dtau,'
      !write(*,*) 'TEXT    =',TEXT,    '  EXT:     ',EXT(1:natm),'  SSA:   ',SSA(1:natm)
      !write(*,*) 'TEXT_HP0=',TEXT_HP0,'  EXT_HP0: ',EXT_HP0(1:natm)

      !do i=1,NH
      !   write(*,*) 'i=',i,'  H(i)=',H(i),'  TH(1:natm,i)',TH(1:natm,i),sum(TH(1:natm,i))
      !end do

!XH   integration from the top
      !TTH=0.0
      !do i=2,NH
      !   TTH=TTH+0.5*(SUM(TH(1:natm,i-1)+TH(1:natm,i)))*(H(i-1)-H(i))
      !end do
      !write(*,*) 'TTH=',TTH,' TEXT=',TEXT-TEXT_HP0,' TEXT_HP0=',TEXT_HP0,' EXT_HP0=',EXT_HP0

!XH   integration from the bottom
      !TTH=0.0
      !do i=NH,2,-1
      !   TTH=TTH+0.5*(SUM(TH(1:natm,i-1)+TH(1:natm,i)))*(H(i-1)-H(i))
      !end do
      !write(*,*) 'TTH=',TTH,' TEXT=',TEXT-TEXT_HP0,' TEXT_HP0=',TEXT_HP0,' EXT_HP0=',EXT_HP0

      HL(:)    = 0.0		
      HL(1)    = H(1)
      HL(NT+1) = H(NH)
       
      TL(:,:)       = 0.0
      TL(1:natm,   1) = TH(1:natm, 1) 
      TL(1:natm,NT+1) = TH(1:natm,NH) 
      TL(0,   1) = SUM(TH(1:natm, 1)) 
      TL(0,NT+1) = SUM(TH(1:natm,NH)) 

      EXTH1(1:natm,1)    = 0.0
      EXTH1(1:natm,NT+1) = EXT(1:natm)-EXT_HP0(1:natm)

      IL=1
      TS(:) = 0.0
				
      LOOP_HEIGHT: do i=2,NH
         H1 = H(i-1)
         H2 = H(i)
         TH1(1:natm) = TH(1:natm,i-1)
         TH1(0)    = SUM(TH1(1:natm))
         TH2(1:natm) = TH(1:natm,i)
         TH2(0)    = SUM(TH2(1:natm))
		  
         TS1(0:natm) = TS(0:natm)
         TS2(0:natm) = 0.5*(TH1(0:natm)+TH2(0:natm))*(H1-H2)
         TS(0:natm)  = TS(0:natm) + TS2(0:natm)

!         write(*,*) '0: i=',i,'  IL=',IL, TS1(0), TS(0), dtau, &
!                                   abs(TS(0)-dtau)/dtau, tiny, &
!                ' - TS1, TS, dtau, abs(TS(0)-dtau)/dtau, tiny'
          
         if (IL .eq. NT .and. i .eq. NH) then
            IL=IL+1
            HL(IL) = H(NH)
            EXTH1(1:natm,IL)=EXT(1:natm)-EXT_HP0(1:natm)
            if(IPRI_additional_info) then
            if ((abs(TS(0)-dtau)/dtau) .ge. tiny) then
               write(*,'(a,4e18.6,i6,e18.6)') 'Warning in get_WD! dtau(NT),dtau,tiny,(dtau(NT)-dtau)/dtau),NT+1, TEXT:  ',  &
                                              TS(0),dtau,tiny,(abs(TS(0)-dtau)/dtau),NT+1, TEXT

!               istop = 1
            end if
            end if ! IPRI_additional_info
            cycle LOOP_HEIGHT
         end if

         a1(0:natm) = (TH2(0:natm)-TH1(0:natm))/(H2-H1)
         b1(0:natm) = TH1(0:natm)-a1(0:natm)*H1

!  --------------------------------------------------------------------------------

         if (abs(TS(0)-dtau)/dtau .le. tiny) then
            TS(0:natm) = 0.0
            IL=IL+1
            HL(IL) = H2
            TL(0:natm,IL)=TH2(0:natm)
            EXTH1(1:natm,IL)=EXTH1(1:natm,IL-1)+TS2(1:natm)
!            write(*,*) '1: i=',i,'  IL=',IL, TS1(0), TS(0), dtau, &
!            abs(TS(0)-dtau)/dtau, tiny,' - TS1, TS, dtau, abs(TS(0)-dtau)/dtau, tiny'
            cycle LOOP_HEIGHT
         else if (TS(0) .gt. dtau) then
            do while(TS(0) .gt. dtau)
               diff  = dtau-TS1(0)
               if (diff .lt. 0.0) then
                  write(tmp_message,'(a,es11.4,a)') 'diff =',diff,'  - diff can not be <0.'
                  G_ERROR(trim(tmp_message))
               end if
               A = 0.5*a1(0)
               B = TH1(0)
               C = diff
			 
               IL = IL+1
               if (IL .gt. NT+1) then
                  write(tmp_message,'(2(a,i0),a)') &
                  'IL = ',IL,'  NT+1 = ',NT+1,'  - Number of layers with tau<=dtau .GT. NT+1'
                  G_ERROR(trim(tmp_message))
               end if
!              ROOT
!               if (abs(A) .ge. tiny1) then
               if (A .ne. 0.0) then
                  call root(A,B,C,DH1,DH2)
                  x1=DH1+H1
                  x2=DH2+H1
                  if(x1 .le. H1 .and. x1 .ge. H2) HL(IL)=x1
                  if(x2 .le. H1 .and. x2 .ge. H2) HL(IL)=x2
               else ! A=0
                  HL(IL) = -C/B+H1
               end if ! A.ne.0.0
!               write(*,*) A,B,C,a1(0),b1(0),TH1(0),TH2(0),i,'  - A,B,C,a1,b1,TH1(0),TH2(0),i'
!               write(*,'(15x,7e14.6,a,i5)') DH1,DH2,H1,H2,x1,x2,HL(IL),'  - DH1,DH2,H1,H2,x1,x2,HL(IL),IL=',IL

               TL(0:natm,IL) = a1(0:natm)*HL(IL)+b1(0:natm)
!              test 1:
!               diff_test= -1.0*(A*(HL(IL)-H1)*(HL(IL)-H1)+B*(HL(IL)-H1))
!               if (abs(diff_test-diff)/diff .gt. tiny) then
!                  write(*,'(3(a,i4),4(a,e14.6))') 'CHECK1:  i-1=',i-1,'  i=',i,'  IL=',IL,  &
!                                                  '  HL(IL)=',HL(IL),'  diff_test=',diff_test,'  diff=',diff
!                  write(*,*) '1: diff_test=',diff_test,' .NE. diff=',diff
!                  stop 'stop in get_WD'
!               end if
!              test 2:
!               diff_test = 0.5*(TH1(0)+TL(0,IL))*(H1-HL(IL))
!               if (abs(diff_test-diff)/diff .gt. tiny) then
!                  write(*,'(3(a,i4),3(a,e14.6))') 'CHECK2:  i-1=',i-1,'  i=',i,'  IL=',IL,  &
!                                                  '  HL(IL)=',HL(IL),'  diff_test=',diff_test,'  diff=',diff
!                  write(*,*) '2: diff_test=',diff_test,' .NE. diff=',diff
!                  stop 'stop in get_WD'
!               end if

               EXTH1(1:natm,IL) = EXTH1(1:natm,IL-1)+TS1(1:natm)+0.5*(TH1(1:natm)+TL(1:natm,IL))*(H1-HL(IL))
               TS2(0:natm) = 0.5*(TL(0:natm,IL)+TH2(0:natm))*(HL(IL)-H2)
!               write(*,'(3e12.4,2f10.4,a)') TS2(0),TL(0,IL),TH2(0),HL(IL),H2,  &
!                                        ' - TS2(0),TL(0,IL),TH2(0),HL(IL),H2'

               if (abs((TS2(0)-dtau)/dtau) .le. tiny) then
                  TS(0:natm) = 0.0
                  IL=IL+1
                  HL(IL) = H2
                  TL(0:natm,IL)=TH2(0:natm)
                  EXTH1(1:natm,IL)=EXTH1(1:natm,IL-1)+TS2(1:natm)
!                  write(*,*) '2: i=',i,'  IL=',IL, TS1(0), TS(0), dtau, &
!                                            abs(TS(0)-dtau)/dtau, tiny, &
!                         ' - TS1, TS, dtau, abs(TS(0)-dtau)/dtau, tiny'
                  if (IL .eq. NT .and. i .eq. NH) then
                     IL=IL+1
                     HL(IL) = H(NH)
                     EXTH1(1:natm,IL)=EXT(1:natm)-EXT_HP0(1:natm)
                  end if
                  cycle LOOP_HEIGHT
               else if(TS2(0) .lt. dtau) then
                  TS(0:natm) = TS2(0:natm)
                  if (IL .eq. NT .and. i .eq. NH) then
                     IL=IL+1
                     HL(IL) = H(NH)
                     EXTH1(1:natm,IL)=EXT(1:natm)-EXT_HP0(1:natm)
                  end if
!                  write(*,*) '3: i=',i,'  IL=',IL, TS1(0), TS(0), dtau, &
!                                            abs(TS(0)-dtau)/dtau, tiny, &
!                         ' - TS1, TS, dtau, abs(TS(0)-dtau)/dtau, tiny'
                  cycle LOOP_HEIGHT
               else if(TS2(0) .gt. dtau) then
                  TS1(0:natm) = 0.0
                  TS(0:natm)  = TS2(0:natm)
                  H1        = HL(IL)
                  TH1(0:natm) = TL(0:natm,IL)
!                  write(*,*) '4: i=',i,'  IL=',IL, TS1(0), TS(0), dtau, &
!                                            abs(TS(0)-dtau)/dtau, tiny, &
!                         ' - TS1, TS, dtau, abs(TS(0)-dtau)/dtau, tiny'
               end if ! TS2(0) .eq. dtau

!               if (IL .eq. NT+1) then
!                  EXTH1(1:natm,IL)=EXT(1:natm)-EXT_HP0(1:natm)
!                  cycle LOOP_HEIGHT
!               end if
            end do ! while(TS(0) .gt. dtau)
         end if ! abs(TS(0)-dtau)/dtau .le. tiny
!  --------------------------------------------------------------------------------

      end do LOOP_HEIGHT

      if (IL .ne. NT+1) then
         write(tmp_message,'(2(a,i0))') &
         'Execution has to be terminated. IL = ',IL,' .ne. NT+1 = ',NT+1
         G_ERROR(trim(tmp_message))
      end if ! IL

      EXTH1(1:natm,NT+1)=EXT(1:natm)-EXT_HP0(1:natm)
		
!      if (ABS(HL(NT+1)-H(NH))/H(NH) .gt. tiny) then
!         write(*,*) '!!! WARNING in get_WD: HL(NT+1) .ne. H(NH)'
!         write(*,*) 'NT+1=',NT+1,'  HL(NT+1)=',HL(NT+1),'  H(NH)=',H(NH)
!         write(*,*) '!!! WARNING in get_WD: H(NH) => HL(NT+1)'
!      end if

!     Calculate output: EXTH(1:NT+1) and WD(1:NT,1:natm)

      do i=1,NT+1
         EXTH(i)=SUM(EXTH1(1:natm,i))
!         write(*,*) i,EXTH(i),EXTH1(1:natm,i),EXT(1:natm),'  - i,EXTH(i),EXTH(1:natm,i),EXT(1:natm)'
      end do ! i
		
!     if (abs(EXTH(NT+1)-(TEXT-TEXT_HP0))/(TEXT-TEXT_HP0) .gt. tiny) then
!        write(*,*)
!        write(*,'(a,i4,2(a,f12.8))') 'NT+1=',NT+1,'  TEXT=',TEXT-TEXT_HP0,'  EXTH(NT+1)=',EXTH(NT+1)
!        write(*,*) 'Total extinction .NE. EXTH(NT+1)'
!        stop 'stop in get_WD'
!     end if
						 
      do i=1,NT
         TTH=SUM((EXTH1(1:natm,i+1)-EXTH1(1:natm,i)))
!         write(*,*) i,TTH,dtau,'  - i,TTH,dtau'
      
         do i1=1,natm
            WD(i,i1)=WD(i,i1)+(EXTH1(i1,i+1)-EXTH1(i1,i))*SSA(i1)/TTH
         end do ! i1
      end do ! i
!AL for testing purposes only
!      do i1=1,natm
!        write(*,*)'Component:', i1
!        write(*,*)'WD='
!        do i=1,NT
!           write(*,*) WD(i,i1)
!        enddo
!      enddo
!      stop 'AL:test stop in get_WD'
!AL for testing purposes only

      if(istop .eq. 1) stop 'test stop in get_WD'
!	------------------------------------------------------------------------------------------------------
      return
      end subroutine get_WD
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss	 

end module mod_os
