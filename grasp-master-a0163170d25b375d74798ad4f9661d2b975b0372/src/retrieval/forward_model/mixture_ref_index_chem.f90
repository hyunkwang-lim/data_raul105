! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

SUBROUTINE MIXING(	instrument, num_wvln, wvln_um, rh_mxtr, &
                    fract_inslbl, fract_carbon, fract_iron, &
                    fract_slbl, fract_wtr, RREAL, RIMAG)

  implicit none 
 

  ! Global variables 

  ! ... input ... 
  character(*)                     instrument
  integer                          num_wvln
  real (selected_real_kind(p=15))  rh_mxtr      ! Valid range: 0--1
  real (selected_real_kind(p=15))  fract_inslbl ! Volume fraction of insoluble inclusions 
  real (selected_real_kind(p=15))  fract_carbon   ! Volume fraction of mixture that is composed of soot+BrC .
  real (selected_real_kind(p=15))  fract_iron   ! Volume fraction of mixture that is composed of iron (aka Hematite).

  ! ... output ...
  real   (selected_real_kind(p=15)), dimension(num_wvln) :: RREAL, RIMAG, wvln_um      
  real   (selected_real_kind(p=15))                         fract_slbl        ! Volume fraction of soluble inclusions 
  real   (selected_real_kind(p=15))                         fract_wtr         ! Volume fraction of mixture that is composed of water.



  ! Local variables 
  integer                           IW
!tl  integer                           i_rh
  character(12)                     slbl_aersl_mtrl !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                    !           from Dinar, Faraday disc, 2008___________|      |        |
                                                    !           from Irshad, ACPD, 2008_________________________|        |
                                                    !           could be Volz sea salt or OC w/o absorption______________| 

  character(12)                     aersl_mtrl      ! iron or quartz

  real (selected_real_kind(p=15)), dimension(num_wvln) ::       refr_iron,       refi_iron,  refr_carbon,       refi_carbon
  real (selected_real_kind(p=15)), dimension(num_wvln) ::     refr_quartz,     refi_quartz
!tl  real (selected_real_kind(p=15)), dimension(num_wvln) ::       refr_soot,       refi_soot


  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_inslbl
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_carbon
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr 
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_iron          

!   WRITE(*,*) rh_mxtr,fract_inslbl, fract_soot, fract_iron, 'rh_mxtr,fract_inslbl, fract_soot, fract_iron'

  ! Iron and quartz refractive indices:
  if     (trim(instrument) == 'aeronet'     ) then 

     call get_refindx_aeronet_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'parasol'     ) then

     call get_refindx_parasol_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'eumetsat'    ) then

     call get_refindx_eumetsat_wvlns (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'user_defined') then

     aersl_mtrl = 'iron'
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_iron  , refi_iron)

!     aersl_mtrl = 'quartz'
     aersl_mtrl = 'mix_dust'   !  modify by lei on 20160301 to use new Ref. index for mix_dust
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_quartz, refi_quartz)

! modify by lei on 20160407 use carbonaceous to replace bc
     aersl_mtrl= 'carbon'
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_carbon  , refi_carbon)
! modify by lei on 20160407 use carbonaceous to replace bc

  endif


!  refindx_soot   = cmplx(       1.68,        0.2)  ! Bond et al., 2006, Aerosol Sci. and Tech.
! modify by lei on 20160407 use carbonaceous to replace bc
  refindx_carbon   = cmplx(refr_carbon, refi_carbon)
  refindx_inslbl = cmplx(refr_quartz, refi_quartz)  ! From Krekov, 1992
  refindx_iron   = cmplx(refr_iron  , refi_iron  )  ! from Dubovik and Derimian via email, Fall 2009

!stop
  slbl_aersl_mtrl = 'ammnm_slft'  !  'ammnm_ntrt'  'ammnm_slft'
 
     call frwrd_mdl (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_inslbl, refindx_carbon, refindx_iron, &
                                                 rh_mxtr,   fract_inslbl,   fract_carbon,   fract_iron, &
                                                 fract_slbl,   fract_wtr   , refindx_mxtr                )  

  do IW = 1,num_wvln
   !  write(6,'(4(f6.4,2x),f8.6)') fract_slbl, fract_wtr,  wvln_um(IW), refindx_mxtr(IW)
     RREAL(IW) =  REAL(refindx_mxtr(IW))
     RIMAG(IW) = AIMAG(refindx_mxtr(IW))

  end do

!  write(*,*) rh_mxtr, fract_slbl, fract_wtr, refindx_mxtr(IW)



end subroutine MIXING




!******************************14/11/2016
SUBROUTINE MIXING_f (instrument, num_wvln, wvln_um, ISD, rh_mxtr_f, 		&
                    fract_quartz_f, fract_bc, fract_brc, 			&
                    fract_slbl_f, fract_wtr_f, RREAL, RIMAG)

  implicit none 
 

  ! Global variables 

  ! ... input ... 
  character(*)                     instrument
  integer                          num_wvln
  real (selected_real_kind(p=15))  rh_mxtr_f      ! Valid range: 0--1
  real (selected_real_kind(p=15))  fract_quartz_f ! Volume fraction of insoluble inclusions
  real (selected_real_kind(p=15))  fract_bc, fract_brc   ! Volume fraction of soot,BrC .

  integer                          ISD

  ! ... output ...
  real   (selected_real_kind(p=15)), dimension(num_wvln) :: RREAL, RIMAG, wvln_um      
  real   (selected_real_kind(p=15))                         fract_slbl_f        ! Volume fraction of soluble inclusions
  real   (selected_real_kind(p=15))                         fract_wtr_f         ! Volume fraction of mixture that is composed of water.


  ! Local variables 
  integer                           IW
  !tl  integer                           i_rh
  character(12)                     slbl_aersl_mtrl !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                    !           from Dinar, Faraday disc, 2008___________|      |        |
                                                    !           from Irshad, ACPD, 2008_________________________|        |
                                                    !           could be Volz sea salt or OC w/o absorption______________| 

  character(12)                     aersl_mtrl      ! iron or quartz

  real (selected_real_kind(p=15)), dimension(num_wvln) ::    refr_quartz,  refi_quartz,  refr_brc,  refi_brc, refr_iron  , refi_iron
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_quartz
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_brc
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_bc

  ! Iron and quartz refractive indices:
  if     (trim(instrument) == 'aeronet'     ) then 

     call get_refindx_aeronet_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'parasol'     ) then

     call get_refindx_parasol_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'eumetsat'    ) then

     call get_refindx_eumetsat_wvlns (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'user_defined') then

     aersl_mtrl = 'quartz'
!     aersl_mtrl = 'mix_dust'   !  modify by lei on 20160301 to use new Ref. index for mix_dust
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_quartz, refi_quartz)

    aersl_mtrl = 'brc'
    call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_brc  , refi_brc)

  endif

  refindx_bc   = cmplx(       1.7,        0.3)
  refindx_quartz = cmplx(refr_quartz, refi_quartz)  ! From Krekov, 1992
  refindx_brc   = cmplx(refr_brc, refi_brc)

  slbl_aersl_mtrl = 'ammnm_slft'  !  'ammnm_ntrt'  'ammnm_slft'
 

     call frwrd_mdl_f (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_quartz, refindx_bc, refindx_brc, &
                                                 rh_mxtr_f,   fract_quartz_f,   fract_bc,   fract_brc, &
                                                 fract_slbl_f,   fract_wtr_f   , refindx_mxtr                )
  do IW = 1,num_wvln
     RREAL(IW) =  REAL(refindx_mxtr(IW))
     RIMAG(IW) = AIMAG(refindx_mxtr(IW))


  end do

end subroutine MIXING_f




SUBROUTINE MIXING_c (	instrument, num_wvln, wvln_um, ISD, rh_mxtr_c, 		&
                    fract_quartz_c, fract_iron, 	&
                    fract_slbl_c, fract_wtr_c, RREAL, RIMAG)

  implicit none 
 

  ! Global variables 

  ! ... input ... 
  character(*)                     instrument
  integer                          num_wvln
  real (selected_real_kind(p=15))  rh_mxtr_c      ! Valid range: 0--1
  real (selected_real_kind(p=15))  fract_quartz_c ! Volume fraction of insoluble inclusions
  real (selected_real_kind(p=15))  fract_iron   ! Volume fraction of mixture that is composed of iron (aka Hematite).

  integer                          ISD

  ! ... output ...
  real   (selected_real_kind(p=15)), dimension(num_wvln) :: RREAL, RIMAG, wvln_um      
  real   (selected_real_kind(p=15))                         fract_slbl_c        ! Volume fraction of soluble inclusions
  real   (selected_real_kind(p=15))                         fract_wtr_c         ! Volume fraction of mixture that is composed of water.



  ! Local variables 
  integer                           IW

  character(12)                     slbl_aersl_mtrl !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                    !           from Dinar, Faraday disc, 2008___________|      |        |
                                                    !           from Irshad, ACPD, 2008_________________________|        |
                                                    !           could be Volz sea salt or OC w/o absorption______________| 

  character(12)                     aersl_mtrl      ! iron or quartz


  real (selected_real_kind(p=15)), dimension(num_wvln) ::    refr_quartz,  refi_quartz, refr_iron, refi_iron
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_quartz
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr 
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_iron          

  ! Iron and quartz refractive indices:
  if     (trim(instrument) == 'aeronet'     ) then 

     call get_refindx_aeronet_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'parasol'     ) then

     call get_refindx_parasol_wvlns  (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'eumetsat'    ) then

     call get_refindx_eumetsat_wvlns (num_wvln, refr_iron  , refi_iron  , &
                                                refr_quartz, refi_quartz  )
  elseif (trim(instrument) == 'user_defined') then

     aersl_mtrl = 'iron'
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_iron  , refi_iron)

     aersl_mtrl = 'quartz'
!     aersl_mtrl = 'mix_dust'   !  modify by lei on 20160301 to use new Ref. index for mix_dust
     call get_refindx_user_wvlns(aersl_mtrl, num_wvln, wvln_um, refr_quartz, refi_quartz)

  endif

  refindx_quartz = cmplx(refr_quartz, refi_quartz)  ! From Krekov, 1992
  refindx_iron   = cmplx(refr_iron  , refi_iron  )  ! from Dubovik and Derimian via email, Fall 2009

  slbl_aersl_mtrl = 'ammnm_slft'  !  'ammnm_ntrt'  'ammnm_slft'  'sea_salt'

     call frwrd_mdl_c (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_quartz, refindx_iron, &
                                                 rh_mxtr_c,   fract_quartz_c,  fract_iron, &
                                                 fract_slbl_c,   fract_wtr_c   , refindx_mxtr    )
  do IW = 1,num_wvln
     RREAL(IW) =  REAL(refindx_mxtr(IW))
     RIMAG(IW) = AIMAG(refindx_mxtr(IW))


  end do

end subroutine MIXING_c
!*******************************14/11/2016






subroutine frwrd_mdl (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_inslbl, refindx_carbon, refindx_iron,  &
                                              rh_mxtr,       fract_inslbl,   fract_carbon,   fract_iron,  &
                                           fract_slbl,       fract_wtr   , refindx_mxtr                 )  


  ! Written by Greg Schuster; initially distributed on 10/15/2008.
  ! This subroutine computes the refractive index of a 4-component aerosol mixture. 

  ! Ambient RH (rh_mxtr), fraction of insoluble aerosols (fract_inslbl), and BC fraction (fract_soot) 
  ! can be tuned to obtain a desired refractive index. 

  ! RH is used to compute the refractive index of the host solution (water and ammonium nitrate, for example)
  ! using partial molar theory (Tang papers). Once the refractive index of the host has been computed, the refractive
  ! index of the 4-component mixture is computed using the Maxwell-Garnett equations. 


  implicit none

  ! Global Variables 
  ! ... input ...
  integer            num_wvln
  real (selected_real_kind(p=15)),    dimension(num_wvln  ) ::  wvln_um       
  character(*)                                slbl_aersl_mtrl   !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                                !           from Dinar, Faraday disc, 2008___________|      |        |
                                                                !           from Irshad, ACPD, 2008_________________________|        |
                                                                !           could be Volz sea salt or OC w/o absorption______________| 

  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_inslbl    ! Refractive index of insoluble inclusions. See test_frwrd_mdl.f90 for examples.
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_carbon
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_iron          

  real (selected_real_kind(p=15))     rh_mxtr   ! Valid range: 0--0.999 or so.
  real (selected_real_kind(p=15))  fract_inslbl ! Volume fraction of insoluble inclusions 
  real (selected_real_kind(p=15))  fract_carbon   ! Volume fraction of mixture that is composed of soot (aka BC).
  real (selected_real_kind(p=15))  fract_iron   ! Volume fraction of mixture that is composed of iron (aka Hematite).


  ! ... output ...
  real     (selected_real_kind(p=15))                         fract_slbl       ! Volume fraction of soluble inclusions 
  real     (selected_real_kind(p=15))                         fract_wtr        ! Volume fraction of mixture that is composed of water.
  complex  (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr



  ! Local Variables 
  integer                                                      i_wave,            i_rhmin
  real    (selected_real_kind(p=15)), dimension(num_wvln) :: refr_host    , refi_host
  real    (selected_real_kind(p=15))                         host_slbl_volfrac                ! Volume Fraction of host solution composed of soluble inclusions 
  real    (selected_real_kind(p=15))                         host_slbl_masfrac                ! Mass   Fraction of host solution composed of soluble inclusions
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_host

  ! Note:
  ! fract_slbl + fract_wtr + frac_soot + fract_inslbl = 1
  ! host_slbl_volfrac                                 = fract_slbl / (fract_slbl + fract_wtr)



  ! Require RH above efflorescence:
  !  ammnm_slft: rh >= 0.37
  !  ammnm_ntrt: rh >= 0.28
  !    sea_salt: rh >= 0.47

  call  get_rhmin (slbl_aersl_mtrl, i_rhmin)

  if (rh_mxtr < i_rhmin/1000.) then
     write(6,*) 'FATAL: rh_mxtr < rh_min... continue w/o processing this RH'
     write(6,*) '       rh_mxtr: ', rh_mxtr
     write(6,*) '  i_rhmin/1000: ', i_rhmin/1000.
     stop
  end if



  do i_wave = 1, num_wvln

     call solution_refindx (rh_mxtr, wvln_um(i_wave), slbl_aersl_mtrl, refr_host(i_wave), refi_host(i_wave), &
                                                                           host_slbl_masfrac    , host_slbl_volfrac      )

     refindx_host(i_wave) = cmplx(refr_host(i_wave), refi_host(i_wave))

     call avg_refindx_3inclsns (  fract_carbon        ,   fract_inslbl            ,  fract_iron          ,                       &
                                refindx_carbon(i_wave), refindx_inslbl(i_wave),  refindx_iron(i_wave), refindx_host(i_wave), &
                                                                                                     refindx_mxtr(i_wave)  )
  enddo

  !fract_slbl =        host_slbl_volfrac  * (1.0 - fract_soot - fract_inslbl)
  !fract_wtr  = (1.0 - host_slbl_volfrac) * (1.0 - fract_soot - fract_inslbl)

  fract_slbl  =           host_slbl_volfrac  * (1.0 - fract_carbon - fract_inslbl - fract_iron)
  fract_wtr =  (1.0 - host_slbl_volfrac) * (1.0 - fract_carbon - fract_inslbl - fract_iron)


 
end subroutine frwrd_mdl




!*********************************14/11/2016
subroutine frwrd_mdl_f (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_quartz, refindx_bc, refindx_brc,  &
                                              rh_mxtr,       fract_quartz,   fract_bc,   fract_brc,  &
                                           fract_slbl,       fract_wtr   , refindx_mxtr                 )  

  implicit none

  ! Global Variables 
  ! ... input ...
  integer            num_wvln
  real (selected_real_kind(p=15)),    dimension(num_wvln  ) ::  wvln_um       
  character(*)                                slbl_aersl_mtrl   !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                                !           from Dinar, Faraday disc, 2008___________|      |        |
                                                                !           from Irshad, ACPD, 2008_________________________|        |
                                                                !           could be Volz sea salt or OC w/o absorption______________| 

  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_quartz    ! Refractive index of insoluble inclusions. See test_frwrd_mdl.f90 for examples.
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_bc
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_brc

  real (selected_real_kind(p=15))     rh_mxtr   ! Valid range: 0--0.999 or so.
  real (selected_real_kind(p=15))  fract_quartz ! Volume fraction of insoluble inclusions
  real (selected_real_kind(p=15))  fract_bc   ! Volume fraction of mixture that is composed of soot (aka BC).
  real (selected_real_kind(p=15))  fract_brc   ! Volume fraction of mixture that is composed of iron (aka Hematite).


  ! ... output ...
  real     (selected_real_kind(p=15))                         fract_slbl       ! Volume fraction of soluble inclusions 
  real     (selected_real_kind(p=15))                         fract_wtr        ! Volume fraction of mixture that is composed of water.
  complex  (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr



  ! Local Variables 
  integer                                                      i_wave,            i_rhmin
  real    (selected_real_kind(p=15)), dimension(num_wvln) :: refr_host    , refi_host
  real    (selected_real_kind(p=15))                         host_slbl_volfrac                ! Volume Fraction of host solution composed of soluble inclusions 
  real    (selected_real_kind(p=15))                         host_slbl_masfrac                ! Mass   Fraction of host solution composed of soluble inclusions
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_host


  ! Require RH above efflorescence:
  !  ammnm_slft: rh >= 0.37
  !  ammnm_ntrt: rh >= 0.28
  !    sea_salt: rh >= 0.47

  call  get_rhmin (slbl_aersl_mtrl, i_rhmin)

  if (rh_mxtr < i_rhmin/1000.) then
     write(6,*) 'FATAL: rh_mxtr < rh_min... continue w/o processing this RH'
     write(6,*) '       rh_mxtr: ', rh_mxtr
     write(6,*) '  i_rhmin/1000: ', i_rhmin/1000.
     stop
  end if

  do i_wave = 1, num_wvln

     call solution_refindx (rh_mxtr, wvln_um(i_wave), slbl_aersl_mtrl, refr_host(i_wave), refi_host(i_wave), &
                                                                           host_slbl_masfrac    , host_slbl_volfrac      )

     refindx_host(i_wave) = cmplx(refr_host(i_wave), refi_host(i_wave))

     call avg_refindx_3inclsns (  fract_quartz,   fract_bc,  fract_brc,                       &
                                refindx_quartz(i_wave), refindx_bc(i_wave),  refindx_brc(i_wave), refindx_host(i_wave), &
                                                                                                     refindx_mxtr(i_wave)  )
  enddo


  fract_slbl  =           host_slbl_volfrac  * (1.0 - fract_quartz - fract_bc - fract_brc)
  fract_wtr =  (1.0 - host_slbl_volfrac) * (1.0 - fract_quartz - fract_bc - fract_brc)


 
end subroutine frwrd_mdl_f



subroutine frwrd_mdl_c (num_wvln,  wvln_um, slbl_aersl_mtrl, refindx_quartz,  refindx_iron,  &
                                              rh_mxtr,       fract_quartz,     fract_iron,  &
                                           fract_slbl,       fract_wtr   , refindx_mxtr                 )

  implicit none

  ! Global Variables 
  ! ... input ...
  integer            num_wvln
  real (selected_real_kind(p=15)),    dimension(num_wvln  ) ::  wvln_um       
  character(*)                                slbl_aersl_mtrl   !  ammnm_slft, ammnm_ntrt, ammnm_bislft, sea_salt, hulis, new_salt, mr_145
                                                                !           from Dinar, Faraday disc, 2008___________|      |        |
                                                                !           from Irshad, ACPD, 2008_________________________|        |
                                                                !           could be Volz sea salt or OC w/o absorption______________| 

  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_quartz    ! Refractive index of insoluble inclusions. See test_frwrd_mdl.f90 for examples.
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_iron

  real (selected_real_kind(p=15))     rh_mxtr   ! Valid range: 0--0.999 or so.
  real (selected_real_kind(p=15))  fract_quartz ! Volume fraction of insoluble inclusions
  real (selected_real_kind(p=15))  fract_iron   ! Volume fraction of mixture that is composed of iron (aka Hematite).


  ! ... output ...
  real     (selected_real_kind(p=15))                         fract_slbl       ! Volume fraction of soluble inclusions 
  real     (selected_real_kind(p=15))                         fract_wtr        ! Volume fraction of mixture that is composed of water.
  complex  (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_mxtr



  ! Local Variables 
  integer                                                      i_wave,            i_rhmin
  real    (selected_real_kind(p=15)), dimension(num_wvln) :: refr_host    , refi_host
  real    (selected_real_kind(p=15))                         host_slbl_volfrac                ! Volume Fraction of host solution composed of soluble inclusions 
  real    (selected_real_kind(p=15))                         host_slbl_masfrac                ! Mass   Fraction of host solution composed of soluble inclusions
  complex (selected_real_kind(p=15)), dimension(num_wvln) :: refindx_host


  ! Require RH above efflorescence:
  !  ammnm_slft: rh >= 0.37
  !  ammnm_ntrt: rh >= 0.28
  !    sea_salt: rh >= 0.47

  call  get_rhmin (slbl_aersl_mtrl, i_rhmin)

  if (rh_mxtr < i_rhmin/1000.) then
     write(6,*) 'FATAL: rh_mxtr < rh_min... continue w/o processing this RH'
     write(6,*) '       rh_mxtr: ', rh_mxtr
     write(6,*) '  i_rhmin/1000: ', i_rhmin/1000.
     stop
  end if



  do i_wave = 1, num_wvln

     call solution_refindx (rh_mxtr, wvln_um(i_wave), slbl_aersl_mtrl, refr_host(i_wave), refi_host(i_wave), &
                                host_slbl_masfrac    , host_slbl_volfrac      )

     refindx_host(i_wave) = cmplx(refr_host(i_wave), refi_host(i_wave))

     call avg_refindx_2inclsns (    fract_quartz,  fract_iron,     &
                                 refindx_quartz(i_wave),  refindx_iron(i_wave), refindx_host(i_wave), &
                                    refindx_mxtr(i_wave)  )
  enddo

  fract_slbl  =           host_slbl_volfrac  * (1.0 - fract_quartz - fract_iron)
  fract_wtr =  (1.0 - host_slbl_volfrac) * (1.0 - fract_quartz - fract_iron)

 
end subroutine frwrd_mdl_c
!*********************************14/11/2016





subroutine get_rhmin (incl_aersl_mtrl, i_rhmin)

  ! This subroutine provides the minimum rh (x10) above crystallization for several aerosol species.

  implicit none

  ! Global variables 
  ! ... input ...
  character(*) incl_aersl_mtrl
  ! ... output ...
  integer i_rhmin
  


  if     (trim(incl_aersl_mtrl) == 'ammnm_slft') then
     i_rhmin = 370  
  elseif (trim(incl_aersl_mtrl) == 'ammnm_ntrt') then
     i_rhmin = 280
  elseif (trim(incl_aersl_mtrl) == 'sea_salt') then
     i_rhmin = 470
  else
     write(6,*) 'WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING'
     write(6,*) 
     write(6,*) 'incl_aersl_mtrl = '//incl_aersl_mtrl//' unrecognized by subroutines mixture_refindx_tm_host and get_rh_min'
     write(6,*) '                   ...assuming that this aerosol species is insoluble, and using i_rhmin=0'
     i_rhmin = 0
  endif


end subroutine get_rhmin




subroutine solution_refindx (rh, wvln_um, slbl_aersl_mtrl, refr_solution, refi_solution, solute_masfrac, solute_volfrac)

  implicit none 
  ! Global variables 
  ! ... input ...
  real (selected_real_kind(p=15)) rh
  real (selected_real_kind(p=15)) wvln_um
  character(*) slbl_aersl_mtrl
  ! ... output ...
  real (selected_real_kind(p=15)) refr_solution, refi_solution
  real (selected_real_kind(p=15)) solute_masfrac, solute_volfrac 

  ! Local Variables 
  real (selected_real_kind(p=15)) dnsty_solute ! Density of soluble aerosol
  real (selected_real_kind(p=15)) dnsty_sltn   ! Density of solution

  refr_solution  = 0.0
  refi_solution  = 0.0
  solute_masfrac = 0.0
  solute_volfrac = 0.0

  if     (trim(slbl_aersl_mtrl) == 'ammnm_slft') then          
     call ammso4 (rh, wvln_um, refr_solution, refi_solution, solute_masfrac)
     dnsty_solute = 0.9971 + (5.92e-3)                - (5.036e-5)                   + (1.024e-8)
     dnsty_sltn   = 0.9971 + (5.92e-3)*solute_masfrac - (5.036e-5)*solute_masfrac**2 + (1.024e-8)*solute_masfrac**3

!  else
elseif (trim(slbl_aersl_mtrl) == 'ammnm_ntrt') then
     call ammno3 (rh, wvln_um, refr_solution, refi_solution, solute_masfrac)
     dnsty_solute = 0.9971 + (4.05e-3)                - (9.000e-6)    
     dnsty_sltn   = 0.9971 + (4.05e-3)*solute_masfrac - (9.000e-6)*solute_masfrac**2 
!
!  elseif (trim(slbl_aersl_mtrl) == 'sea_salt') then
!     call seasal (rh, wvln_um, refr_solution, refi_solution, solute_masfrac)
!     dnsty_solute = 0.9971 + (7.4100e-3) - (3.741e-5) + (2.2252e-6) - (2.060e-8)
!     dnsty_sltn   = 0.9971 + (7.4100e-3)*solute_masfrac     - (3.741e-5)*solute_masfrac**2 &
!                           + (2.2252e-6)*solute_masfrac**3  - (2.060e-8)*solute_masfrac**4

  else 
     write(6,*) 'FATAL in subroutine solution_refindx; invalid slbl_aersl_mtrl: ',slbl_aersl_mtrl 
     stop
  endif

  solute_volfrac = solute_masfrac*dnsty_sltn/dnsty_solute



end subroutine solution_refindx




subroutine avg_refindx_3inclsns(fract_1, fract_2, fract_3, refindx_1, refindx_2, refindx_3, refindx_host, refindx_avg)


  ! Modified April 2008 by GLS.

  ! This subroutine computes the average complex refractive index 
  ! for a mixture of two "inclusions" embedded into an otherwise 
  ! homogeneous matrix using Maxwell-Garnett. 
  ! For calculations with a single inclusion, use
  !     fract_2 = 0.0
  !   refindx_2 = (0.0, 0.0)  
  
  ! Ref Bohren and Huffman, Sec 8.5, p213, and Eq. (8.50).


  implicit none 
  ! Global variables
  ! ... input ...
  real    (selected_real_kind(p=15))   fract_1,     fract_2,   fract_3 ! volume fraction of inclusions
  complex (selected_real_kind(p=15)) refindx_1,   refindx_2, refindx_3 ! refractive index of inclusions
  complex (selected_real_kind(p=15)) refindx_host                      ! refractive index of host matrix
  ! ... output ... 
  complex (selected_real_kind(p=15)) refindx_avg                       ! avg refractive index of composite
  
  ! Local variables 
  complex (selected_real_kind(p=15)) diel_1,   diel_2, diel_3 ! dielectric constant of inclusions
  complex (selected_real_kind(p=15)) diel_host                ! dielectric constant of host matrix
  complex (selected_real_kind(p=15)) diel_avg                 ! average dielectric constant



  call refindx_to_dielectric (refindx_1,    diel_1   )
  call refindx_to_dielectric (refindx_2,    diel_2   )
  call refindx_to_dielectric (refindx_3,    diel_3   )
  call refindx_to_dielectric (refindx_host, diel_host)

  call avg_dielectric3 (fract_1, fract_2, fract_3, diel_1, diel_2, diel_3, diel_host, diel_avg  )

  call dielectric_to_refindx (diel_avg, refindx_avg  )


end subroutine avg_refindx_3inclsns




!**********************************14/11/2016
subroutine avg_refindx_2inclsns(fract_1, fract_2, refindx_1, refindx_2, refindx_host, refindx_avg)


  ! Modified April 2008 by GLS.

  ! This subroutine computes the average complex refractive index 
  ! for a mixture of two "inclusions" embedded into an otherwise 
  ! homogeneous matrix using Maxwell-Garnett. 
  ! For calculations with a single inclusion, use
  !     fract_2 = 0.0
  !   refindx_2 = (0.0, 0.0)  
  
  ! Ref Bohren and Huffman, Sec 8.5, p213, and Eq. (8.50).


  implicit none 
  ! Global variables
  ! ... input ...
  real    (selected_real_kind(p=15))   fract_1,     fract_2   ! volume fraction of inclusions
  complex (selected_real_kind(p=15)) refindx_1,   refindx_2   ! refractive index of inclusions
  complex (selected_real_kind(p=15)) refindx_host                      ! refractive index of host matrix
  ! ... output ... 
  complex (selected_real_kind(p=15)) refindx_avg                       ! avg refractive index of composite
  
  ! Local variables 
  complex (selected_real_kind(p=15)) diel_1,   diel_2         ! dielectric constant of inclusions
  complex (selected_real_kind(p=15)) diel_host                ! dielectric constant of host matrix
  complex (selected_real_kind(p=15)) diel_avg                 ! average dielectric constant



  call refindx_to_dielectric (refindx_1,    diel_1   )
  call refindx_to_dielectric (refindx_2,    diel_2   )
  call refindx_to_dielectric (refindx_host, diel_host)

  call avg_dielectric2 (fract_1, fract_2, diel_1, diel_2, diel_host, diel_avg  )

  call dielectric_to_refindx (diel_avg, refindx_avg  )


end subroutine avg_refindx_2inclsns
!**********************************14/11/2016




subroutine refindx_to_dielectric (refindx, dielectric)

  ! Written 3/29/01 by GLS
  ! Given the complex refractive index, this subroutine computes
  ! the complex dielectric fcn. Ref Bohren & Huffman, p 227.

  implicit none 
  ! Global variables
  ! ... input ... 
  complex (selected_real_kind(p=15)) refindx    ! Refractive index
  ! ... output ... 
  complex (selected_real_kind(p=15)) dielectric ! Dielectric constant

  ! Local variables 
  real (selected_real_kind(p=15)) refr,  refi  ! Real and imaginary refractive index
  real (selected_real_kind(p=15)) dielr, dieli ! Real and imaginary dielectric constants

  ! Get real and imaginary components of refractive index
  refr =  real(refindx) 
  refi = aimag(refindx) 
  
  ! Calculate the real and imaginary components of dielectric constants
  dielr = refr**2 - refi**2
  dieli = 2.*refr*refi

  dielectric = cmplx(dielr, dieli)



end subroutine refindx_to_dielectric





subroutine dielectric_to_refindx (dielectric, refindx) 


  ! Written 3/29/01 by GLS
  ! This program converts complex dielectric constants to complex
  ! refractive index. Ref Bohren and Huffman, p 227.

  implicit none
  ! Global variables 
  ! ... input ...
  complex (selected_real_kind(p=15)) dielectric ! Complex dielectric function
  ! ... output
  complex (selected_real_kind(p=15)) refindx    ! Complex refractive index

  ! Local variables 
  real (selected_real_kind(p=15)) dielr, dieli  ! Complex dielectric function
  real (selected_real_kind(p=15)) refr,  refi   ! Complex refractive index


  ! Get real and imaginary components of dielectric function
  dielr =  real(dielectric) 
  dieli = aimag(dielectric) 

  ! Calculate real and imaginary components of refractive index
  refr = (((dielr**2 + dieli**2)**0.5 + dielr) / 2.)**0.5
  refi = (((dielr**2 + dieli**2)**0.5 - dielr) / 2.)**0.5

  refindx = cmplx(refr, refi) 

end subroutine dielectric_to_refindx





subroutine avg_dielectric3(fract_1, fract_2, fract_3, diel_1, diel_2, diel_3, diel_host, diel_avg)


  ! Written 3/29/01 by GLS
  ! This subroutine computes the average complex dielectric 
  ! function for a 2-component mixture of "inclusions" embedded
  ! in an otherwise homogeneous matrix using Maxwell-Garnett. 
  ! Ref Bohren and Huffman, Sec 8.5, p213, and Eq. (8.50).

  implicit none 

  ! Global variables:
  ! .... input ... 
  real    (selected_real_kind(p=15)) fract_1,  fract_2, fract_3  ! Volume fraction of embedded inclusions
  complex (selected_real_kind(p=15)) diel_1,   diel_2,  diel_3  ! Dielectric constant of embedded inclusions
  complex (selected_real_kind(p=15)) diel_host                   ! Dielectric constant of homogeneous matrix
  ! ... output 
  complex (selected_real_kind(p=15)) diel_avg           ! Average dielectric constant of the mixture

  ! Local variables:
  complex (selected_real_kind(p=15)) ratio_1, ratio_2, ratio_3   ! 


  ratio_1    = (diel_1 - diel_host) / (diel_1 + 2.*diel_host)
  ratio_2    = (diel_2 - diel_host) / (diel_2 + 2.*diel_host)
  ratio_3    = (diel_3 - diel_host) / (diel_3 + 2.*diel_host)


  diel_avg = diel_host * (1. + 3.*(     fract_1*ratio_1 + fract_2*ratio_2 + fract_3*ratio_3) /  &
                                  (1. - fract_1*ratio_1 - fract_2*ratio_2 - fract_3*ratio_3)    )

end subroutine avg_dielectric3



!****************************14/11/2016
subroutine avg_dielectric2(fract_1, fract_2, diel_1, diel_2, diel_host, diel_avg)


  ! Written 3/29/01 by GLS
  ! This subroutine computes the average complex dielectric 
  ! function for a 2-component mixture of "inclusions" embedded
  ! in an otherwise homogeneous matrix using Maxwell-Garnett. 
  ! Ref Bohren and Huffman, Sec 8.5, p213, and Eq. (8.50).

  implicit none 

  ! Global variables:
  ! .... input ... 
  real    (selected_real_kind(p=15)) fract_1,  fract_2 ! Volume fraction of embedded inclusions
  complex (selected_real_kind(p=15)) diel_1,   diel_2  ! Dielectric constant of embedded inclusions
  complex (selected_real_kind(p=15)) diel_host                   ! Dielectric constant of homogeneous matrix
  ! ... output 
  complex (selected_real_kind(p=15)) diel_avg           ! Average dielectric constant of the mixture

  ! Local variables:
  complex (selected_real_kind(p=15)) ratio_1, ratio_2


  ratio_1    = (diel_1 - diel_host) / (diel_1 + 2.*diel_host)
  ratio_2    = (diel_2 - diel_host) / (diel_2 + 2.*diel_host)


  diel_avg = diel_host * (1. + 3.*(     fract_1*ratio_1 + fract_2*ratio_2 ) /  &
                                  (1. - fract_1*ratio_1 - fract_2*ratio_2 )    )

end subroutine avg_dielectric2
!****************************14/11/2016




subroutine get_refindx_aeronet_wvlns(num_wvln, refr_iron      , refi_iron      , &
                                               refr_quartz    , refi_quartz      )


  implicit none
  
  ! Global variables 
  ! ... input ...
  integer num_wvln
  ! ... output ...
  real (selected_real_kind(p=15)), dimension(num_wvln) ::           refr_iron,           refi_iron
  real (selected_real_kind(p=15)), dimension(num_wvln) ::         refr_quartz,         refi_quartz


!  wvln_um = (/  .440,  .675,  .870, 1.020 /)
  
  refr_iron = (/ 2.8556, 2.9480, 2.7409, 2.6587 /)
  refi_iron = (/ 3.8076E-01, 4.2060E-03, 3.1190E-03, 2.8734E-05 /)

  refr_quartz = (/ 1.5400, 1.5350, 1.5230, 1.5200 /)
  refi_quartz = (/ 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04 /)


end subroutine get_refindx_aeronet_wvlns





subroutine get_refindx_parasol_wvlns(num_wvln, refr_iron      , refi_iron      , &
                                               refr_quartz    , refi_quartz      )

  implicit none

  ! Global variables 
  ! ... input ...
  integer num_wvln
  ! ... output ...
  real (selected_real_kind(p=15)), dimension(num_wvln) ::           refr_iron,           refi_iron
  real (selected_real_kind(p=15)), dimension(num_wvln) ::         refr_quartz,         refi_quartz
 

!  wvln_um = (/  .443,  .490,  .565,  .675,  .870, 1.020 /)
 
  refr_iron = (/ 2.8692, 3.0498, 3.0849, 2.9480, 2.7409, 2.6587 /)
  refi_iron = (/ 3.7009E-01, 2.1856E-01, 7.7542E-02, 4.2060E-03, 3.1190E-03, 2.8734E-05 /)
    
  refr_quartz = (/ 1.5400, 1.5400, 1.5400, 1.5350, 1.5230, 1.5200 /)
  refi_quartz = (/ 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04 /)
  
end subroutine get_refindx_parasol_wvlns





subroutine get_refindx_eumetsat_wvlns(num_wvln, refr_iron      , refi_iron      , &
                                                refr_quartz    , refi_quartz      )
  
  implicit none
  
  ! Global variables 
  ! ... input ...
  integer num_wvln
  ! ... output ...
  real (selected_real_kind(p=15)), dimension(num_wvln) ::           refr_iron,           refi_iron
  real (selected_real_kind(p=15)), dimension(num_wvln) ::         refr_quartz,         refi_quartz
 

!  wvln_um = (/  .443,  .510,  .640,  .865, 1.380, 1.640, 2.130 /)
  
  refr_iron = (/ 2.8692, 3.0900, 2.9954, 2.7450, 2.6378, 2.6264, 2.6100 /)
  refi_iron = (/ 3.7009E-01, 1.6802E-01, 1.4371E-02, 3.1045E-03, 1.0018E-05, 5.1209E-06, 4.8900E-06 /)
  
  refr_quartz = (/ 1.5400, 1.5400, 1.5400, 1.5235, 1.5148, 1.5044, 1.4974 /)
  refi_quartz = (/ 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04, 5.0000E-04 /)
  
end subroutine get_refindx_eumetsat_wvlns






subroutine get_refindx_user_wvlns(aersl_mtrl, num_wvln,  wvln_um, refr_mtrl, refi_mtrl)


  implicit none 

  ! Global variables 
  ! ...input... 
  character(*) aersl_mtrl
  integer      num_wvln
  real (selected_real_kind(p=15)), dimension(num_wvln) ::    wvln_um
  ! ... output...
  real (selected_real_kind(p=15)), dimension(num_wvln) :: refr_mtrl, refi_mtrl
  
  ! Local variables 

  real    (selected_real_kind(p=15)), dimension(      61) ::    wvln_tbl , refr_tbl , refi_tbl
  complex (selected_real_kind(p=15)), dimension(      61) :: refindx_tbl
  real    (selected_real_kind(p=15)), dimension(      25) ::    wvln_tblq, refr_tblq, refi_tblq 
  complex (selected_real_kind(p=15)), dimension(      25) :: refindx_tblq
  real    (selected_real_kind(p=15)), dimension(      25) ::    wvln_tblb, refr_tblb, refi_tblb
  complex (selected_real_kind(p=15)), dimension(      25) :: refindx_tblb


  refr_tbl  = 0.0; refi_tbl  = 0.0; refindx_tbl  = 0.0
  refr_tblq = 0.0; refi_tblq = 0.0; refindx_tblq = 0.0
  refr_tblb = 0.0; refi_tblb = 0.0; refindx_tblb = 0.0

  if (trim(aersl_mtrl) == 'iron') then 

     call get_refindx_iron (61, wvln_tbl, refindx_tbl) 
     refr_tbl =  real(refindx_tbl)
     refi_tbl = aimag(refindx_tbl)
     call  linear_near_neighbor_fit_dp (      61 , wvln_tbl, refr_tbl, &
                                        num_wvln, wvln_um , refr_mtrl )
     call  linear_near_neighbor_fit_dp (      61, wvln_tbl, refi_tbl , &
                                        num_wvln, wvln_um , refi_mtrl  )

  elseif (trim(aersl_mtrl) == 'quartz') then
     if (maxval(wvln_um) > 4.0) then
        write(6,*) 'FATAL: Quartz refractive index unavailable for wavelengths greater than 4 um'; stop
     endif

     call get_refindx_quartz (25, wvln_tblq, refindx_tblq)
     refr_tblq(1:25) =  real(refindx_tblq)
     refi_tblq(1:25) = aimag(refindx_tblq)
     call  linear_near_neighbor_fit_dp (      25 , wvln_tblq, refr_tblq, &
                                        num_wvln, wvln_um , refr_mtrl )
     call  linear_near_neighbor_fit_dp (      25, wvln_tblq, refi_tblq , &
                                        num_wvln, wvln_um , refi_mtrl  )

  elseif (trim(aersl_mtrl) == 'brc') then
     if (maxval(wvln_um) > 4.0) then
        write(6, *) 'FATAL: BrC refractive index unavailable for wavelenths greater than 4 um'; stop
     endif

     call get_refindx_brc (25, wvln_tblb, refindx_tblb)
     refr_tblb(1:25) = real(refindx_tblb)
     refi_tblb(1:25) = aimag(refindx_tblb)
     call  linear_near_neighbor_fit_dp (      25 , wvln_tblb, refr_tblb, &
                                        num_wvln, wvln_um , refr_mtrl )
     call  linear_near_neighbor_fit_dp (      25, wvln_tblb, refi_tblb , &
                                        num_wvln, wvln_um , refi_mtrl  )
!!!******add by lei for carbonaceous refractive index
  elseif (trim(aersl_mtrl) == 'carbon') then
     if (maxval(wvln_um) > 4.0) then
        write(6, *) 'FATAL: carbonaceous refractive index unavailable for wavelenths greater than 4 um'; stop
     endif

     call get_refindx_carbon (25, wvln_tblb, refindx_tblb)
     refr_tblb(1:25) = real(refindx_tblb)
     refi_tblb(1:25) = aimag(refindx_tblb)
     call  linear_near_neighbor_fit_dp (      25 , wvln_tblb, refr_tblb, &
                                        num_wvln, wvln_um , refr_mtrl )
     call  linear_near_neighbor_fit_dp (      25, wvln_tblb, refi_tblb , &
                                        num_wvln, wvln_um , refi_mtrl  )


!!!******add by lei for mix_dust refractive index

  elseif (trim(aersl_mtrl) == 'mix_dust') then
     if (maxval(wvln_um) > 4.0) then
        write(6,*) 'FATAL: Mix_dust refractive index unavailable for wavelengths greater than 4 um'; stop
     endif

     call get_refindx_mix_dust(25, wvln_tblq, refindx_tblq)
     refr_tblq(1:25) =  real(refindx_tblq)
     refi_tblq(1:25) = aimag(refindx_tblq)
     call  linear_near_neighbor_fit_dp (      25 , wvln_tblq, refr_tblq, &
                                        num_wvln, wvln_um , refr_mtrl )
     call  linear_near_neighbor_fit_dp (      25, wvln_tblq, refi_tblq , &
                                        num_wvln, wvln_um , refi_mtrl  )

!!!********end by lei

  else
     write(6,*) 'FATAL in subroutine get_refindx_user_wvlns; invalid aersl_mtrl: ', trim(aersl_mtrl); stop
  endif






end subroutine get_refindx_user_wvlns





subroutine linear_near_neighbor_fit_dp (num_knowns, xknown, yknown, &
                                        num_fits,   xfit,   yfit    )

  ! Written 12/00 by GLS
  ! For a given number (num_knowns) of values (xknown, yknown), this program 
  ! interpolates between these values with a linear fit at a desired 
  ! number (num_fits) of abscissas (xfit). If unknown values are located
  ! outside of the range of known values, the slope of the nearest 2 points
  ! are extrapolated. Output is yfit.

  ! NOTE: This is not a linear regression, but rather a linear interpolation
  ! between nearest neighbors.

  implicit none

  ! Global variables
  ! ... Input
  integer                        num_knowns ! # known abscissas and ordinates
  integer                        num_fits   ! Number of desired interpolations
  real (selected_real_kind(p=15)), dimension(num_knowns) :: xknown     ! Known abscissas
  real (selected_real_kind(p=15)), dimension(num_knowns) :: yknown     ! Known ordinates
  real (selected_real_kind(p=15)), dimension(num_fits)   :: xfit       ! Abscissas of desired interps
  ! ... Output
  real (selected_real_kind(p=15)), dimension(num_fits)   :: yfit       ! Interpolated ordinates (at xfit)

  ! Local variables
  integer i
  integer ifit       ! Index of desired interpolations
  integer i_high     ! Index of 1st known abscissa > than desired abscissa
  real (selected_real_kind(p=15))  xdum       ! Scalar representation of a desired abscissa
  real (selected_real_kind(p=15))  ydum       ! Scalar representation of a desired ordinate
  real (selected_real_kind(p=15)), dimension(num_knowns) :: xknown_dum, yknown_dum ! 



  ! Check to be sure that abscissa values are in ascending order:
  if (xknown(num_knowns) < xknown(1)) then
     write(6,*) "SUBROUTINE LINEAR_FIT: "
     write(6,*) "xknown found to be in descending order."
     write(6,*) "This is OK.....Reassigning...."
     write(6,*) 
     do i = 1,num_knowns
        xknown_dum(i) = xknown(num_knowns - (i-1))
        yknown_dum(i) = yknown(num_knowns - (i-1))
     enddo
  else
     xknown_dum = xknown
     yknown_dum = yknown
  endif



  ! Linear interpolation/extrapolation:
  yfit = 0.0
  do ifit = 1, num_fits
     
     ! Determine the location of the desired abscissas within 
     ! the known abscissas
     i_high = 1
     do while ( (xknown_dum(i_high) < xfit(ifit)) .and. (i_high <= num_knowns) )
        i_high = i_high + 1
     enddo

     xdum = xfit(ifit) 
     ! If the abscissa is located outside of the known span of abscissas, 
     ! extrapolate slope of interior points. Otherwise, interpolate 
     ! between nearest neighbors.
     if ( (i_high == 1) .or. (i_high == num_knowns+1) ) then
        call extrapolate_linear_fit_dp &
             (xdum, num_knowns, i_high, xknown_dum, yknown_dum, ydum)
     else 
        call interpolate_linear_fit_dp (xdum, i_high, xknown_dum, yknown_dum, ydum)
     endif
     yfit(ifit) = ydum

  end do


end subroutine linear_near_neighbor_fit_dp





subroutine interpolate_linear_fit_dp (xdum, i_high, xknown, yknown, ydum)

  implicit none
  ! Global variables 
  ! ... Input
  real (selected_real_kind(p=15))                 xdum
  integer               i_high
  real (selected_real_kind(p=15)), dimension(*) :: xknown 
  real (selected_real_kind(p=15)), dimension(*) :: yknown 
  ! ... Output
  real (selected_real_kind(p=15))                 ydum
  
  ! Local variables
  real                  slope


  slope = (yknown(i_high) - yknown(i_high-1)) / & 
          (xknown(i_high) - xknown(i_high-1))

  ydum = yknown(i_high) + slope*(xdum - xknown(i_high))

end subroutine interpolate_linear_fit_dp





subroutine extrapolate_linear_fit_dp &
     (xdum, num_knowns, i_high, xknown, yknown, ydum)


  implicit none
  ! Global variables
  ! ... Input
  real (selected_real_kind(p=15))                           xdum
  integer                        num_knowns
  integer                        i_high
  real (selected_real_kind(p=15)), dimension(num_knowns) :: xknown
  real (selected_real_kind(p=15)), dimension(num_knowns) :: yknown
  ! .... Output 
  real (selected_real_kind(p=15))                           ydum
  
  ! Local variables: 
  real (selected_real_kind(p=15))                           slope 
  
  
  if (i_high == 1) then
     slope = (yknown(2) - yknown(1)) / & 
             (xknown(2) - xknown(1))
     ydum  = yknown(1)  - slope*(xknown(1) - xdum)
  else if (i_high == num_knowns+1) then
     slope = (yknown(num_knowns) - yknown(num_knowns-1)) / &
             (xknown(num_knowns) - xknown(num_knowns-1)) 
     ydum  = yknown(num_knowns)  + slope*(xdum - xknown(num_knowns))
  else
     write(6,*) "FATAL ERROR: subroutine linear_fit.f90:"
     write(6,*) "i_high /= 1 or num_knowns+1 in subroutine extrapolate"
     stop
  endif

end subroutine extrapolate_linear_fit_dp





subroutine get_refindx_iron (num_tbl, wavelen, refindx_iron)

! These refractive indices were obtained from Oleg via email on Nov 5, 2009. 
! Oleg got them from Yevgeny Derimian.


 implicit none

  ! Global variables 
  ! ... input ... 
  integer                        num_tbl
!tl  integer                        I
  ! ... output ... 
  real    (selected_real_kind(p=15)), dimension(num_tbl) :: wavelen
  complex (selected_real_kind(p=15)), dimension(num_tbl) :: refindx_iron

!  real    , dimension(num_tbl) :: wavelen
!  complex , dimension(num_tbl) :: refindx_iron

  ! Local variables 
  real (selected_real_kind(p=15)),    dimension(num_tbl) :: refr_iron_Shettle, refi_iron_Shettle
!  real ,    dimension(num_tbl) :: refr_iron_Shettle, refi_iron_Shettle




  if (num_tbl == 61) then

     wavelen = (/  0.250000,  0.300000,  0.350000,  0.400000,   0.450000, &
                   0.500000,  0.550000,  0.600000,  0.650000,   0.700000, &
                   0.750000,  0.800000,  0.900000,  1.000000,   1.250000, &
                   1.500000,  1.750000,  2.000000,  2.500000,   3.000000, &
                   3.200000,  3.390000,  3.500000,  3.750000,   4.000000, &
                   4.500000,  5.000000,  5.500000,  6.000000,   6.200000, &
                   6.500000,  7.200000,  7.900000,  8.200000,   8.500000, &
                   8.700000,  9.000000,  9.200000,  9.500000,   9.800000, &
                   10.000000, 10.600000, 11.000000, 11.500000, 12.500000, &
                   13.000000, 14.000000, 14.800000, 15.000000, 16.400000, &
                   17.200001, 18.000000, 18.500000, 20.000000, 21.299999, &
                   22.500000, 25.000000, 27.900000, 30.000000, 35.000000, &
                   40.000000                                              /)
     
     refr_iron_Shettle = (/	2.070,	2.320,	2.480,	2.674,	2.901, &
         			3.087,	3.102,	3.045,	2.983,	2.913, &
                                2.856,	2.799,	2.716,	2.660,	2.644, &
                                2.632,	2.622,	2.610,	2.610,	2.610, &
                                2.610,	2.610,	2.610,	2.610,	2.603, &
                                2.593,	2.580,	2.567,	2.547,	2.537, &
                                2.520,	2.473,	2.427,	2.390,	2.363, &
                                2.343,	2.317,	2.290,	2.253,	2.217, &
                                2.197,	2.102,	2.037,	1.945,	1.746, &
                                1.636,	1.246,	0.769,	0.592,	0.244, &
                                0.350,	0.711,	1.468,	2.812,	1.643, &
                                2.940,	2.649,	0.558,	0.583, 13.862, &
                                7.370	/)
     
     refi_iron_Shettle = (/	1.330E+00,	1.180E+00,	9.730E-01,	5.230E-01,	3.452E-01,&
                                1.869E-01,	9.250E-02,	4.264E-02,	7.304E-03,	1.108E-03,&
                                2.012E-03,	2.916E-03,	3.206E-03,	3.000E-05,	1.417E-05,&
                                6.186E-06,	4.284E-06,	4.500E-06,	6.000E-06,	1.000E-05,&
                                1.100E-05,	1.100E-05,	1.200E-05,	1.200E-05,	4.900E-05,&
                                4.630E-05,	4.370E-05,	4.230E-04,	8.120E-04,	1.010E-03,&
                                1.140E-03,	1.220E-03,	1.230E-03,	2.170E-03,	3.100E-03,&
                                4.070E-03,	5.000E-03,	5.930E-03,	7.330E-03,	8.730E-03,&
                                9.670E-03,	1.396E-02,	1.670E-02,	2.130E-02,	3.660E-02,&
                                4.750E-02,	8.010E-02,	2.170E-01,	3.070E-01,	1.640E+00,&
                                2.380E+00,	3.350E+00,	4.140E+00,	1.470E+00,	2.570E+00,&
                                4.290E+00,	1.020E+00,	2.680E+00,	4.860E+00,	1.010E+01,&
                                3.870E-01	/)
!     Do I=1,25
!     refi_iron_Shettle(I)=refi_iron_Shettle(I)*3.0
!     write(*,*) refi_iron_Shettle(I),I
!     ENDDO
!    refi_iron_Shettle=refi_iron_Shettle*3.0
!     refr_iron_Shettle = refr_iron_Shettle * 0.7
  else
     write(6,*) 'FATAL in subroutine get_refindx_iron.f90; invalid value for num_tbl:', num_tbl
     stop
  end if
  

  refindx_iron = cmplx(refr_iron_Shettle, refi_iron_Shettle)

end subroutine get_refindx_iron





subroutine get_refindx_quartz (num_tbl25, wvln_tbl25, refindx_tbl25) 

  implicit none 

  ! Global variables 
  ! ... input ... 
  integer                                                       num_tbl25
  ! ... output ... 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) ::    wvln_tbl25
  complex (selected_real_kind(p=15)), dimension(num_tbl25) :: refindx_tbl25
  
  ! Local variables 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) :: refr_tbl25, refi_tbl25

  if (num_tbl25 == 25) then 
     wvln_tbl25 = (/  0.250000,  0.300000,  0.350000,  0.400000,  0.450000, &
                      0.500000,  0.550000,  0.600000,  0.650000,  0.700000, &
                      0.750000,  0.800000,  0.900000,  1.000000,  1.250000, &
                      1.500000,  1.750000,  2.000000,  2.500000,  3.000000, &
                      3.200000,  3.390000,  3.500000,  3.750000,  4.000000  /)
     
     refr_tbl25 = (/ 1.54, 1.54, 1.54, 1.54, 1.54, &
                     1.54, 1.54, 1.54, 1.54, 1.53, &
                     1.53, 1.53, 1.52, 1.52, 1.52, &
                     1.51, 1.50, 1.50, 1.49, 1.49, &
                     1.48, 1.48, 1.48, 1.48, 1.47  /)       
     
     !refr_tbl25( 1:25) = 1.45
     
     refi_tbl25( 1:25) = 5.000D-04
     refi_tbl25(   23) = 1.000D-04

     refindx_tbl25 = cmplx(refr_tbl25, refi_tbl25)
  else 
    write(6,*) 'FATAL in subroutine get_refindx_quartz.f90; invalid value for num_tbl:', num_tbl25
    stop
  end if

end subroutine get_refindx_quartz



subroutine get_refindx_brc (num_tbl25, wvln_tbl25, refindx_tbl25)

implicit none


  ! Global variables 
  ! ... input ... 
  integer                                                       num_tbl25
  ! ... output ... 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) ::    wvln_tbl25
  complex (selected_real_kind(p=15)), dimension(num_tbl25) :: refindx_tbl25
  
  ! Local variables 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) :: refr_tbl25, refi_tbl25

  if (num_tbl25 == 25) then
     wvln_tbl25 = (/  0.250000,  0.300000,  0.350000,  0.400000,  0.450000, &
                      0.500000,  0.550000,  0.600000,  0.650000,  0.700000, &
                      0.750000,  0.800000,  0.900000,  1.000000,  1.250000, &
                      1.500000,  1.750000,  2.000000,  2.500000,  3.000000, &
                      3.200000,  3.390000,  3.500000,  3.750000,  4.000000  /)
     
     refr_tbl25 = (/ 1.54, 1.54, 1.54, 1.54, 1.54, &
                     1.54, 1.54, 1.54, 1.54, 1.54, &
                     1.54, 1.54, 1.54, 1.52, 1.51, &
                     1.50, 1.49, 1.47, 1.47, 1.47, &
                     1.47, 1.47, 1.47, 1.47, 1.47  /)

!!!*****20160608 test by lei*********
!     refr_tbl25 = (/ 1.56, 1.56, 1.56, 1.56, 1.56, &
!                     1.56, 1.56, 1.56, 1.56, 1.56, &
!                     1.56, 1.56, 1.56, 1.54, 1.53, &
!                     1.52, 1.51, 1.49, 1.49, 1.49, &
!                     1.49, 1.49, 1.49, 1.49, 1.49  /)
!!!*****20160608 test by lei*********
!
     refi_tbl25 = (/ 0.0698, 0.0698, 0.0698, 0.0698, 0.0667, &
                     0.0453, 0.0295, 0.0177, 0.0081, 0.0046, &
                     0.0039, 0.0032, 0.0022, 0.0020, 0.0015, &
                     0.0010, 0.0006, 0.0002, 0.0002, 0.0002, &
                     0.0002, 0.0002, 0.0002, 0.0002, 0.0002  /)
!
!! improving BrC imaginary refractive index from Chakrabarty R. K. 2010 by lei on 20160203
!     refi_tbl25 = (/ 0.015, 0.015, 0.013, 0.0103, 0.0043, &
!                     0.0024, 0.0015, 0.0010, 0.0008, 0.0007, &
!                     0.0007, 0.0005, 0.0005, 0.0005, 0.0005, &
!                     0.0005, 0.0005, 0.0002, 0.0002, 0.0002, &
!                     0.0002, 0.0002, 0.0002, 0.0002, 0.0002  /)
!!

! on 20160613 by lei two times
!     refi_tbl25 = (/ 0.030, 0.030, 0.0260, 0.0206, 0.0086, &
!                     0.0048, 0.0030, 0.0020, 0.0016, 0.0014, &
!                     0.0014, 0.0010, 0.0010, 0.0010, 0.0010, &
!                     0.0008, 0.0006, 0.0002, 0.0002, 0.0002, &
!                     0.0002, 0.0002, 0.0002, 0.0002, 0.0002  /)

!! three times
!
!     refi_tbl25 = (/ 0.0450, 0.0450, 0.0390, 0.0309, 0.0129, &
!                     0.0072, 0.0045, 0.0030, 0.0024, 0.0021, &
!                     0.0021, 0.0015, 0.0015, 0.0015, 0.0015, &
!                     0.0010, 0.0006, 0.0002, 0.0002, 0.0002, &
!                     0.0002, 0.0002, 0.0002, 0.0002, 0.0002  /)
!



!! improving BrC imaginary refractive index from Chakrabarty R. K. 2010 by lei on 20160203

     refindx_tbl25 = cmplx(refr_tbl25, refi_tbl25)
  else 
    write(6,*) 'FATAL in subroutine get_refindx_brc.f90; invalid value for num_tbl:', num_tbl25
    stop
  end if

end subroutine get_refindx_brc

!!!!!!************** add by lei create carbonaceous refractive index table

subroutine get_refindx_carbon (num_tbl25, wvln_tbl25, refindx_tbl25)

implicit none


  ! Global variables 
  ! ... input ... 
  integer                                                       num_tbl25
  ! ... output ... 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) ::    wvln_tbl25
  complex (selected_real_kind(p=15)), dimension(num_tbl25) :: refindx_tbl25
  
  ! Local variables 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) :: refr_tbl25, refi_tbl25

  if (num_tbl25 == 25) then
     wvln_tbl25 = (/  0.250000,  0.300000,  0.350000,  0.400000,  0.450000, &
                      0.500000,  0.550000,  0.600000,  0.650000,  0.700000, &
                      0.750000,  0.800000,  0.900000,  1.000000,  1.250000, &
                      1.500000,  1.750000,  2.000000,  2.500000,  3.000000, &
                      3.200000,  3.390000,  3.500000,  3.750000,  4.000000  /)
     
     refr_tbl25 = (/ 1.62, 1.62, 1.62, 1.62, 1.62, &
                     1.62, 1.62, 1.62, 1.62, 1.62, &
                     1.62, 1.62, 1.62, 1.61, 1.60, &
                     1.59, 1.58, 1.57, 1.57, 1.57, &
                     1.57, 1.57, 1.57, 1.57, 1.57  /)


!     refi_tbl25 = (/ 0.0925, 0.0925, 0.0907, 0.0883, 0.0829, &
!                     0.0812, 0.0804, 0.0799, 0.0797, 0.0796, &
!                     0.0796, 0.0795, 0.0795, 0.0795, 0.0795, &
!                     0.0795, 0.0795, 0.0795, 0.0792, 0.0792, &
!                     0.0792, 0.0792, 0.0792, 0.0792, 0.0792  /)
!
!! improved by lei on 20160405
!
!     refr_tbl25 = (/ 1.663, 1.663, 1.663, 1.663, 1.663, &
!                     1.663, 1.663, 1.663, 1.663, 1.663, &
!                     1.663, 1.663, 1.663, 1.649, 1.642, &
!                     1.635, 1.628, 1.614, 1.614, 1.614, &
!                     1.614, 1.614, 1.614, 1.614, 1.614  /)
!
     refi_tbl25 = (/ 0.0635, 0.0635, 0.0617, 0.0593, 0.0539, &
                     0.0522, 0.0514, 0.0509, 0.0507, 0.0506, &
                     0.0506, 0.0505, 0.0505, 0.0505, 0.0505, &
                     0.0505, 0.0505, 0.0502, 0.0502, 0.0502, &
                     0.0502, 0.0502, 0.0502, 0.0502, 0.0502  /)
!






!! improved by lei on 20160405



     refindx_tbl25 = cmplx(refr_tbl25, refi_tbl25)
  else 
    write(6,*) 'FATAL in subroutine get_refindx_brc.f90; invalid value for num_tbl:', num_tbl25
    stop
  end if

end subroutine get_refindx_carbon



!!!!!!************** add by lei create mix_dust refractive index table

subroutine get_refindx_mix_dust (num_tbl25, wvln_tbl25, refindx_tbl25)

implicit none


  ! Global variables 
  ! ... input ... 
  integer                                                       num_tbl25
  ! ... output ... 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) ::    wvln_tbl25
  complex (selected_real_kind(p=15)), dimension(num_tbl25) :: refindx_tbl25
  
  ! Local variables 
  real    (selected_real_kind(p=15)), dimension(num_tbl25) :: refr_tbl25, refi_tbl25

  if (num_tbl25 == 25) then
     wvln_tbl25 = (/  0.250000,  0.300000,  0.350000,  0.400000,  0.450000, &
                      0.500000,  0.550000,  0.600000,  0.650000,  0.700000, &
                      0.750000,  0.800000,  0.900000,  1.000000,  1.250000, &
                      1.500000,  1.750000,  2.000000,  2.500000,  3.000000, &
                      3.200000,  3.390000,  3.500000,  3.750000,  4.000000  /)
     
     refr_tbl25 = (/ 1.522, 1.522, 1.522, 1.522, 1.520, &
                     1.520, 1.517, 1.518, 1.520, 1.515, &
                     1.511, 1.507, 1.496, 1.491, 1.489, &
                     1.484, 1.479, 1.479, 1.475, 1.471, &
                     1.467, 1.467, 1.468, 1.468, 1.463  /)
     

     
     refi_tbl25 = (/ 0.00133, 0.00115, 0.00089, 0.00067, 0.00057, &
                     0.00052, 0.00046, 0.00047, 0.00053, 0.00057, &
                     0.00058, 0.00058, 0.00060, 0.00060, 0.00063, &
                     0.00068, 0.00073, 0.00076, 0.00528, 0.02000, &
                     0.02156, 0.02286, 0.02319, 0.02312, 0.02312  /)



     refindx_tbl25 = cmplx(refr_tbl25, refi_tbl25)
  else 
    write(6,*) 'FATAL in subroutine get_refindx_mix_dust.f90; invalid value for num_tbl:', num_tbl25
    stop
  end if

end subroutine get_refindx_mix_dust





! ********************************************************************************************

      subroutine mixture_chemical_components ( NSD, ISD, nfract, vfract_ap, GOUT_chem_pixel, vfract, chem_name )

      use mod_retr_general_output_derived_type

      implicit none
! --------------------------------------------------------------------------------------------
      integer, intent(in) :: NSD, ISD, nfract
      real, dimension(nfract), intent(in) :: vfract_ap ! not normalized volume fraction of chem component
      type(output_pixel_chemistry), intent(inout) :: GOUT_chem_pixel
      real, dimension(nfract), intent(out) :: vfract ! normalized volume fraction of chem component
      character(12), dimension(nfract), intent(out) :: chem_name
! --------------------------------------------------------------------------------------------
      integer :: ifract
      real :: xnorm
! --------------------------------------------------------------------------------------------
      chem_name(:) = ''
      vfract(:) = 0.0

       xnorm = sum(vfract_ap(1:nfract))
      vfract(1:nfract) = vfract_ap(1:nfract)/xnorm ! normalized volume fractions of chem components


      if ( NSD .eq. 1 ) then ! one component aerosol

        do ifract=1,nfract
          select case(ifract)
          case(1)
            chem_name(ifract) = 'carbon'
            GOUT_chem_pixel%fsoot(ISD) = vfract(ifract)
! remove brc element in volume mixture model by lei on 20160330
!          case(2)
!            chem_name(ifract) = 'brc'
!            GOUT_chem_pixel%fslbl(ISD) = vfract(ifract)
          case(2)
            chem_name(ifract) = 'mix_dust'  ! change 'quartz' to 'mix_dust' by lei using mix_dust refractive index
            GOUT_chem_pixel%fquartz(ISD) = vfract(ifract)
          case(3)
            chem_name(ifract) = 'iron'
            GOUT_chem_pixel%firon(ISD) = vfract(ifract)
          case(4)
            chem_name(ifract) = 'water'
            GOUT_chem_pixel%fwtr(ISD) = vfract(ifract)
          end select
        enddo ! ifract

!!*********improve two modes by lei on 20160608**start*****

      elseif ( NSD .eq. 2 ) then ! 2 component aerosol
        if( ISD .eq. 1 ) then ! fine mode

          do ifract=1,nfract
            select case(ifract)
          case(1)
            chem_name(ifract) = 'bc'
            GOUT_chem_pixel%fsoot(ISD) = vfract(ifract)
          case(2)
            chem_name(ifract) = 'brc'
            GOUT_chem_pixel%fbrc(ISD) = vfract(ifract)
          case(3)
            chem_name(ifract) = 'quartz'
            GOUT_chem_pixel%fquartz(ISD) = vfract(ifract)
!          case(3)
!            chem_name(ifract) = 'iron'
!            GOUT_chem_pixel%firon(ISD) = vfract(ifract)
          case(4)
            chem_name(ifract) = 'water'
            GOUT_chem_pixel%fwtr(ISD) = vfract(ifract)
            end select
          enddo ! ifract
        elseif( ISD .eq. 2 ) then ! coarse mode

          do ifract=1,nfract
            select case(ifract)
!          case(1)
!            chem_name(ifract) = 'bc'
!            GOUT_chem_pixel%fsoot(ISD) = vfract(ifract)
!          case(2)
!            chem_name(ifract) = 'brc'
!            GOUT_chem_pixel%fslbl(ISD) = vfract(ifract)
!          case(1)
!            chem_name(ifract) = 'carbon'
!            GOUT_chem_pixel%fsoot(ISD) = vfract(ifract)
          case(1)
            chem_name(ifract) = 'quartz'
            GOUT_chem_pixel%fquartz(ISD) = vfract(ifract)
          case(2)
            chem_name(ifract) = 'iron'
            GOUT_chem_pixel%firon(ISD) = vfract(ifract)
          case(3)
            chem_name(ifract) = 'water'
            GOUT_chem_pixel%fwtr(ISD) = vfract(ifract)
            end select
          enddo ! ifract
        endif ! ISD .eq. 1
      endif ! NSD .eq. 1
!!*********improve two modes by lei on 20160608**end*****

      end subroutine mixture_chemical_components

! ********************************************************************************************

      subroutine complex_refr_index_practical_mixture( nfract,vfract,chem_name, &
                                                       num_wvln,wvln_um,RREAL,RIMAG )

      implicit none
! --------------------------------------------------------------------------------------------
      integer, intent(in) :: nfract, num_wvln
      real (selected_real_kind(p=15)), dimension(num_wvln), intent(in) :: wvln_um
      real, dimension(nfract), intent(in) :: vfract
      character(12), dimension(nfract), intent(in) :: chem_name
      real (selected_real_kind(p=15)), dimension(num_wvln), intent(out) :: RREAL, RIMAG
! --------------------------------------------------------------------------------------------
      real (selected_real_kind(p=15)), dimension(num_wvln) :: refr, refi
      integer :: ifract
! --------------------------------------------------------------------------------------------
      RREAL(:) = 0.0
      RIMAG(:) = 0.0

!!*********improve two modes by lei on 20160608**start*****
      do ifract=1,nfract
        select case(trim(chem_name(ifract)))

        case('bc')
            refr(1:num_wvln) = 1.7 ! Bond et al., 2006, Aerosol Sci. and Tech.
            refi(1:num_wvln) = 0.3 ! Bond et al., 2006, Aerosol Sci. and Tech.

        case('carbon')
        call get_refindx_user_wvlns(chem_name(ifract), num_wvln,wvln_um, refr,refi)

        case('brc')
            call get_refindx_user_wvlns(chem_name(ifract), num_wvln, wvln_um, refr, refi)

        case('mix_dust')   ! change 'quartz' to 'mix_dust' by lei using mix_dust refractive index
            call get_refindx_user_wvlns(chem_name(ifract), num_wvln, wvln_um, refr, refi)

        case('quartz')     ! from Krekov, 1992
            call get_refindx_user_wvlns(chem_name(ifract), num_wvln, wvln_um, refr, refi)

        case('iron')       ! from Dubovik and Derimian via email, Fall 2009
            call get_refindx_user_wvlns(chem_name(ifract), num_wvln, wvln_um, refr, refi)

        case('water')
            refr(1:num_wvln) = 1.33
            refi(1:num_wvln) = 0.0000
        end select
!!*********improve two modes by lei on 20160608**end*****

        RREAL(1:num_wvln) = RREAL(1:num_wvln) + vfract(ifract)*refr(1:num_wvln)
        RIMAG(1:num_wvln) = RIMAG(1:num_wvln) + vfract(ifract)*refi(1:num_wvln)
!write(*,'(a10,i4,f9.4,2e12.4,a)') trim(chem_name(ifract)),ifract,refr,refi,vfract(ifract), &
!'  - ifract,refr,refi,vfract  in complex_refr_index_practical_mixture'
      enddo ! ifract

      end subroutine complex_refr_index_practical_mixture

! ********************************************************************************************



