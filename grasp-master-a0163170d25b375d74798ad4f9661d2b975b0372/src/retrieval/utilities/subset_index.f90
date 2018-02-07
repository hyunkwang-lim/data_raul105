! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine returns indices in real array a for array b (subset if a) elements
        !> 
        !> @param[in]  tiny - tiny value
        !> @param[in]  na  - number of elements in array a
        !> @param[in]  a   - ascending order array
        !> @param[in]  nb  - number of elements in array b (b is a subset of a)
        !> @param[in]  b   - ascending order array
        !> @param[out] index - indices of array b elements in array a
        !> @param[out] status - status of return (false or true)

      subroutine R_subset_index ( tiny, na, a, nb, b, index, status )

      implicit none
!  ---------------------------------------------------------------------------------
      integer, intent(in)  :: na, nb
      real,    intent(in)  :: a(na), b(nb), tiny
      integer, intent(out) :: index(nb)
      logical, intent(out) :: status
!  ---------------------------------------------------------------------------------
      integer :: i, i1, icount
!  ---------------------------------------------------------------------------------
      index(1:nb) = -999
      status = .false.

! Validate size of arrays: size( a ) must be >= size( b )
      if ( nb .gt. na ) then
        write(*,'(/,2(a,i0,x),x,a)') 'nb = ',nb,'> na = ',nb,'in R_subset_index'
        return
      endif

! Search for indices
      icount = 0
      i1 = 1
      do while ( i1 .le. nb )
        do i=1,na
! AL          if ( abs((b(i1)-a(i))/b(i1)) .le. tiny ) then
          if ( abs(b(i1)-a(i)) .le. tiny ) then
            index(i1) = i
            icount = icount + 1
          exit
          endif
        enddo
      i1 = i1 + 1
      enddo 

! Validate index values
      do i=1,nb
        if(index(i) .eq. -999) then
         write(*,'(/,2(a,i0,2x),a)') 'i = ',i,'index = ',index,'in R_subset_index'
         write(*,'(a,/)') ' !!!  Bad index for array b element  !!!'
         return
        endif
      enddo

! Validate if indices found for all array b elements
      if ( icount .ne. nb ) then
        write(*,'(/,2(a,i0,2x),a)') 'nb = ',nb,'icount = ',icount,'in R_subset_index'
        write(*,'(a,/)') ' !!!  Can not find index for array b element  !!!'
        return
      endif

      status = .true.

      return
      end subroutine R_subset_index
  
! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

