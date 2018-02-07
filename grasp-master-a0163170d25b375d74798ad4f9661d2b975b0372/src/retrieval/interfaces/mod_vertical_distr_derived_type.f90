! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
module mod_vertical_distr_derived_type

      !use mod_time_utils
      !use iso_c_binding
      use mod_par_OS,  only : KSD, KVERT_WD, NMM, NMG
      !use mod_par_inv, only : NMM, NMG
      ! NMM=KSD+1 particles + molecular

      implicit none 
!	----------------------------------------------------------------------
	type discret_vertical_distribution

      integer :: nh ! number of heights for discret vertical distribution
      real, dimension(KVERT_WD) :: h_km  ! height values
      integer :: natm, naer, ncld, nmol, ngas
      real, dimension(NMM+NMG) :: norm  
      !character(len=3) :: atmcomp(4) ! aer / cld / mol / gas
      real, dimension(KVERT_WD,NMM+NMG) :: val  ! vertical distribution values

	end type discret_vertical_distribution

!	----------------------------------------------------------------------
      contains
      
      subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module()
      end subroutine dummy_subroutine_for_avoiding_warnings_in_an_empty_module                

end module mod_vertical_distr_derived_type

! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
