! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm

module mod_par_type_surface

  implicit none

  integer,parameter ::	par_type_surface_beg = 20000
! Land
      integer,parameter ::	par_type_SURF1_land_beg = 20100
        integer,parameter ::	par_type_SURF1_land_Ross_Li_BRDF  = 20101   
        integer,parameter ::	par_type_SURF1_land_RPV_BRDF      = 20102   
        integer,parameter ::	par_type_SURF1_land_Litvinov      = 20103   
        integer,parameter ::	par_type_SURF1_land_Litvinov_fast = 20104
      integer,parameter ::	par_type_SURF1_land_end = 20200
      integer,parameter ::	par_type_SURF2_land_beg = 20200
        integer,parameter ::	par_type_SURF2_land_Maignan_Breon = 20201   
        integer,parameter ::	par_type_SURF2_land_Litvinov      = 20202 
      integer,parameter ::	par_type_SURF2_land_end = 20300
! Water
      integer,parameter ::	par_type_SURF_water_beg = 20300
        integer,parameter ::	par_type_SURF_water_Cox_Munk_iso  = 20301 
        integer,parameter ::	par_type_SURF_water_Cox_Munk_ani  = 20302 
        integer,parameter ::	par_type_SURF_water_Litvinov      = 20303   
      integer,parameter ::	par_type_SURF_water_end = 20400
  integer,parameter ::	par_type_surface_end = 30000

      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module                

end module mod_par_type_surface

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
