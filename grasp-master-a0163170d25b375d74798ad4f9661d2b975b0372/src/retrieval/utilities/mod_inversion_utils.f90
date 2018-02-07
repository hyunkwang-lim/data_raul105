! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

module mod_inversion_utils

contains

!	ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

function check_nan(num_el,input_array)
	implicit none
!	----------------------------------------------------------------------
	integer,intent(in)							::	num_el
	real,dimension(:),intent(in)		::	input_array
	logical							::	check_nan
!	----------------------------------------------------------------------
	integer							::	ii
!	----------------------------------------------------------------------
	do ii = 1,num_el	
!tl		if(isnan(input_array(ii))) then
		if(input_array(ii) .ne. input_array(ii)) then
			write(*,*) "Found an error in check_nan in module mod_ip_utility"
			write(*,*) "Array element number",ii,"is NaN"
			write(*,*) 
			check_nan = .false.
#ifdef DEBUG
      call abort
#endif
			return
		endif
	enddo
!	----------------------------------------------------------------------
	check_nan = .true.
!	----------------------------------------------------------------------
end function check_nan

!	ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

	function geom2scat_angl( tetas,		  &
                           vis,		    &
                           fis,		    &
                           scat_angl	&
                          )
! function returns scattering angle from geometry configuration
!	 vis   - zenith observation angles
!	 fis   - azimuth observation angles
!	 tetas - value of solar Zenith Angle
		implicit none
!  ---------------------------------------------------------------------
		logical		::	geom2scat_angl
!	----------------------------------------------------------------------
		real,intent(in)  :: tetas
		real,intent(in)  :: vis
		real,intent(in)  :: fis
		real,intent(out) :: scat_angl
!	----------------------------------------------------------------------
		real :: pi
		real :: zz,xx,rmus
!	----------------------------------------------------------------------
		if(tetas .lt. 0) then
			write(*,*) "Error: Solar Zenith Angle is less then zero"
			write(*,*) "In geom2scat_angl in module mod_inversion_utils"
			geom2scat_angl 	= .false. 
			return
		endif
		if(vis .lt. 0) then
			write(*,*) "Error: azimuth observation angle is less then zero"
			write(*,*) "In geom2scat_angl in module mod_inversion_utils"
			geom2scat_angl 	= .false. 
			return
		endif
! PL Azimuth angle can be negative.
		!if(fis .lt. 0) then
		!	write(*,*) "Error: zenith observation angle is less then zero"
		!	write(*,*) "In geom2scat_angl in module mod_inversion_utils"
		!	geom2scat_angl 	= .false.
		!	return
		!endif
!	----------------------------------------------------------------------
      pi = acos(-1.0)
      rmus = cos(pi*tetas/180.0)
      xx = cos(vis*pi/180.)
      zz = -rmus*xx-sqrt(1.-rmus*rmus)*sqrt(1.-xx*xx)*cos(fis*pi/180.)
      scat_angl = acos(zz)*180./pi
!	----------------------------------------------------------------------
		geom2scat_angl = .true.
!	----------------------------------------------------------------------
	end function geom2scat_angl
!	ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module mod_inversion_utils
