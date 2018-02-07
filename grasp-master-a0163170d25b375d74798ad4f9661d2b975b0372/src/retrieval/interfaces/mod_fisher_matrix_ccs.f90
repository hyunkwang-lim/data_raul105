! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

#include "../constants_set/mod_globals.inc"
module mod_fisher_matrix_ccs

      use mod_par_inv, only : KPAR,nnz_par
      use mod_stop_report

      implicit none
! -----------------------------------------------------------------------
! Compressed Column Storage (CCS)
! The CCS format is specified by the arrays {val, row, nval}. 
! nnz - number of nonzeros 
! val - nonzeros 
! row stores the row indices of each nonzero
! nval stores the index of the elements in val which start 
! a column of UF matrix
! nval(ncol+1)=nnz+1 where ncol - number of columns in matrix
!
! in addition to CCS, col stores column index of each nonzero

      !type nonzero
         !integer, dimension(nnz_par)         :: row,col
         !real*8,  dimension(nnz_par)         :: val
         !integer, dimension(KPAR+1)          :: nval
      !end type nonzero    
      type nonzero
         integer, dimension(:), allocatable  :: row,col
         real*8,  dimension(:), allocatable  :: val
         integer, dimension(:), allocatable  :: nval
      end type nonzero    

! -----------------------------------------------------------------------

      contains

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine allocate_sparse_matrix_storage ( UFNZ_CCS )

      implicit none
! ----------------------------------------------------------------------- 
!     Compressed Column Storage (CCS)

      type(nonzero),    intent(inout)  :: UFNZ_CCS 
! ----------------------------------------------------------------------- 
      integer  ::  alloc_stat
! ----------------------------------------------------------------------- 
       
      if(.not. allocated(UFNZ_CCS%row)) then
        allocate(UFNZ_CCS%row(1:nnz_par),stat=alloc_stat)
        if (alloc_stat /= 0) then
        !write(*,*) 'size(UFNZ_CCS%row)*4',size(UFNZ_CCS%row)*4,'  nnz_par=',nnz_par
        write(tmp_message,'(2(a,i0))') &
        'error while trying to allocate UFNZ_CCS%row(1:nnz_par)'
        G_ERROR(trim(tmp_message))
        endif
      endif
      if(.not. allocated(UFNZ_CCS%col)) then
        allocate(UFNZ_CCS%col(1:nnz_par),stat=alloc_stat)
        if (alloc_stat /= 0) then
        !write(*,*) 'size(UFNZ_CCS%col)*4',size(UFNZ_CCS%col)*4,'  nnz_par=',nnz_par
        write(tmp_message,'(2(a,i0))') &
        'error while trying to allocate UFNZ_CCS%col(1:nnz_par)'
        G_ERROR(trim(tmp_message))
        endif
      endif
      if(.not. allocated(UFNZ_CCS%val)) then
        allocate(UFNZ_CCS%val(1:nnz_par),stat=alloc_stat)
        if (alloc_stat /= 0) then
        !write(*,*) 'size(UFNZ_CCS%val)*8',size(UFNZ_CCS%val)*8,'  nnz_par=',nnz_par
        write(tmp_message,'(2(a,i0))') &
        'error while trying to allocate UFNZ_CCS%val(1:nnz_par)'
        G_ERROR(trim(tmp_message))
        endif
      endif
      if(.not. allocated(UFNZ_CCS%nval)) then
        allocate(UFNZ_CCS%nval(1:KPAR+1),stat=alloc_stat)
        if (alloc_stat /= 0) then
        !write(*,*) 'size(UFNZ_CCS%val)*4',size(UFNZ_CCS%nval)*4,'  KPAR+1=',KPAR+1
        write(tmp_message,'(2(a,i0))') &
        'error while trying to allocate UFNZ_CCS%nval(1:KPAR+1)'
        G_ERROR(trim(tmp_message))
        endif
      endif
      
      return
      end subroutine allocate_sparse_matrix_storage

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine deallocate_sparse_matrix_storage ( UFNZ_CCS )

      implicit none
! ----------------------------------------------------------------------- 
!     Compressed Column Storage (CCS)

      type(nonzero),    intent(inout)  :: UFNZ_CCS 
! ----------------------------------------------------------------------- 
      integer  ::  alloc_stat
! ----------------------------------------------------------------------- 
       
         if(allocated(UFNZ_CCS%row)) then
            deallocate(UFNZ_CCS%row,stat=alloc_stat)
            if (alloc_stat /= 0) then
            write(tmp_message,'(2(a,i0))') &
            'error while trying to deallocate UFNZ_CCS%row(1:nnz_par)'
            G_ERROR(trim(tmp_message))
            endif
         endif
         if(allocated(UFNZ_CCS%col)) then 
            deallocate(UFNZ_CCS%col,stat=alloc_stat)
            if (alloc_stat /= 0) then
            write(tmp_message,'(2(a,i0))') &
            'error while trying to deallocate UFNZ_CCS%col(1:nnz_par)'
            G_ERROR(trim(tmp_message))
            endif
         endif
         if(allocated(UFNZ_CCS%val)) then 
            deallocate(UFNZ_CCS%val,stat=alloc_stat)
            if (alloc_stat /= 0) then
            write(tmp_message,'(2(a,i0))') &
            'error while trying to deallocate UFNZ_CCS%val(1:nnz_par)'
            G_ERROR(trim(tmp_message))
            endif
         endif
         if(allocated(UFNZ_CCS%nval)) then
            deallocate(UFNZ_CCS%nval,stat=alloc_stat)
            if (alloc_stat /= 0) then
            write(tmp_message,'(2(a,i0))') &
            'error while trying to deallocate UFNZ_CCS%nval(1:KPAR+1)'
            G_ERROR(trim(tmp_message))
            endif
         endif
         
      return
      end subroutine deallocate_sparse_matrix_storage

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine UFNZ_symmetry_check ( KN,nnz,UFNZ,  & ! IN
                                       UF            & ! OUT
                                     )
	  	  
      implicit none
 
! ----------------------------------------------------------------------- 
      integer,                    intent(in)  :: KN,nnz
      type(nonzero),              intent(in)  :: UFNZ
      real, dimension(KPAR,KPAR), intent(out) :: UF
! -----------------------------------------------------------------------      
      integer                       :: i,j  											   
! -----------------------------------------------------------------------

! CHECK out a symmetry of UFNZ
write(*,*) 'test in ufnz_symmetry_check; nnz=',nnz,'  KN=',KN
      UF(:,:) = 0.0
      do i=1, nnz
        UF(UFNZ%row(i),UFNZ%col(i))=UFNZ%val(i)
write(*,*) 'i=',i
      enddo ! i
      do j=1,KN
        do i=1,KN
          if(UF(i,j) .ne. UF(i,j)) then
            write(*,*) 'in UFNZ_symmetry_check j=',j,'  i=',i,'  UF(i,j)=',UF(i,j),'  UF(j,i)=',UF(j,i)
          endif 
        enddo ! i
      enddo  ! j

      return
      end subroutine UFNZ_symmetry_check

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine UF_nonzero_one_pixel (                 &
                                          KNSINGF,UFS,  & ! IN
                                          UFNZ,nnz      & ! OUT
                                      )

      use mod_par_inv, only : KPARS,KPAR
	  	  
      implicit none
 
! ----------------------------------------------------------------------- 
      integer,                      intent(in)    :: KNSINGF
      real, dimension(KPARS,KPARS), intent(in)    :: UFS

      type(nonzero),                intent(inout) :: UFNZ
      integer,                      intent(out)   :: nnz
! -----------------------------------------------------------------------      
      real                          :: XA
      integer                       :: I,J,irow_mtr,irow_vec  											   
      integer, dimension(KPAR)      :: UF_nrow_col
! -----------------------------------------------------------------------
   
      UF_nrow_col(:) = 0

      UFNZ%row(:) = 0
      UFNZ%col(:) = 0
      UFNZ%val(:) = 0.0	

      nnz=0
      do J=1,KNSINGF
         do I=1,KNSINGF
            XA=UFS(I,J)
		      if(XA .eq. 0.) cycle
            nnz=nnz+1
            irow_mtr=I
            irow_vec=nnz
            UFNZ%col(irow_vec) = J
            UFNZ%row(irow_vec) = irow_mtr				 
            UFNZ%val(irow_vec) = XA             
            UF_nrow_col(J)=UF_nrow_col(J)+1
            if(irow_mtr .eq. 0) then
               write(tmp_message,'(a,i0,2x,a)') &
               '1: irow_mtr = ',irow_mtr,'can not be 0.'
               G_ERROR(trim(tmp_message))
            endif
!                 UF((IIMAGEQ1-1)*KNSINGF+I,(IIMAGEQ1-1)*KNSINGF+I1)=  &
!                 UFS(IIMAGEQ1,I,I1)
         enddo ! I
      UFNZ%nval(J) = sum(UF_nrow_col(1:J))-UF_nrow_col(J)+1          	
      enddo ! J

      UFNZ%nval(KNSINGF+1) = nnz+1
	  
      return
      end subroutine UF_nonzero_one_pixel
      
      
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine UF1_nonzero (                  &
                              KNSINGF,npixels,  & ! IN
                              IP2,NIPPN,IPPN,   &
                              UFNZ,             & ! INOUT
                              UFNZ1,            &
                              nnz               & ! OUT
                             )

      use mod_par_inv, only : KPAR,nnz_par
      use mod_index_cloud
	  	  
      implicit none
 
! ----------------------------------------------------------------------- 
      integer,                               intent(in)  :: KNSINGF,npixels
      integer,                               intent(in)  :: IP2
      integer, dimension(KPAR),              intent(in)  :: NIPPN                         
      integer, dimension(KPAR,KPAR),         intent(in)  :: IPPN

      type(nonzero),                         intent(inout) :: UFNZ
      type(nonzero),                         intent(inout) :: UFNZ1
      integer,                               intent(out)   :: nnz
! -----------------------------------------------------------------------      
      integer                :: IIPP,IP3,IP4
      integer                :: I,I1,I2,II1,II2,J,KNF,KN1
      real                   :: AA
      integer                :: irow_vec,irow_vec1
								
      real, dimension(KPAR)       :: TEMP
      integer, dimension(2,KPAR)  :: UF1_nrow_col ! 1 - number of nz in column 
                                                  ! 2 - index of column first element in UFNZ1  
      logical, save               :: lfirst=.true.
! -----------------------------------------------------------------------
      KNF=npixels*KNSINGF
      UFNZ1%row(:)=0
      UFNZ1%col(:)=0
      UFNZ1%val(:)=0.0
      IIPP=0
      KN1=IIPP+IP2
      UF1_nrow_col(:,:)=0         

! ***  1
      nnz=0
      IF(NIPPN(IP2+1) .GT. 0) THEN
         IIPP=NIPPN(IP2+1)
         KN1=IIPP+IP2
!write(*,*) IP2,NIPPN(IP2+1),KN1,'  IP2,NIPPN(IP2+1),KN1 in nonzero'
         UF1_nrow_col(2,1) = 1
         DO J=1,IIPP	! do not change J (out) and I (inside) loop order	   
         DO I=1,IIPP
            irow_vec =ivec2 (				                       &
                              KNSINGF,npixels,             &
                              IPPN(IP2+1,I),IPPN(IP2+1,J)  &
                            )
			   
            IF(irow_vec .eq. -999) CYCLE
            AA=UFNZ%val(irow_vec)
            IF(AA .EQ. 0.) CYCLE  		   

            nnz=nnz+1
            irow_vec1 = UF1_nrow_col(2,J) + UF1_nrow_col(1,J)
            UF1_nrow_col(1,J) = UF1_nrow_col(1,J)+1         			   

            UFNZ1%row(irow_vec1) = I				 
            UFNZ1%col(irow_vec1) = J
            UFNZ1%val(irow_vec1) = AA        

!            UF1(I,J)=UF(IPPN(IP2+1,I),IPPN(IP2+1,J))
!if(j .eq. IIPP) write(*,*) 'in nonzero: j=',j,'  i=',i,'  UF1=',UFNZ1%val(irow_vec1)
         ENDDO ! I

!write(*,*) 'j=',j,'  ip2=',ip2,'  first=',UF1_nrow_col(2,J),' nz=',UF1_nrow_col(1,J)

         UF1_nrow_col(2,J+1) = UF1_nrow_col(2,J) + UF1_nrow_col(1,J) + IP2

         ENDDO ! J


!write(*,*) 'nnz=',nnz,' nnz1=',sum(UF1_nrow_col(1,1:IIPP))

!write(*,*) '1: j=',IIPP+1,'  first=',UF1_nrow_col(2,IIPP+1)

         if(IP2 .gt. 1) then
         do j=2,IP2
            UF1_nrow_col(2,IIPP+j) = UF1_nrow_col(2,IIPP+j-1) + KN1
!write(*,*) '2: j=',IIPP+j,'  first=',UF1_nrow_col(2,IIPP+j),'  KN1=',KN1,'  IIPP=',IIPP,'  KN=',KN

         enddo ! j
         endif ! IP2 .gt. 1      

      ENDIF ! NIPPN(IP2+1) .GT. 0
      
! ***  2
      DO IP3=1,IP2     
         TEMP(:)=0.0
         DO I=1,KNF
         AA=0.0
         DO I1=1,NIPPN(IP3)
         irow_vec =ivec2 (				           &
                           KNSINGF,npixels,  &
                           I,IPPN(IP3,I1)    &
							    )

         IF(irow_vec .eq. -999) CYCLE
         AA=AA+UFNZ%val(irow_vec)

!         AA=AA+UF(I,IPPN(IP3,I1))
         ENDDO ! I1
         TEMP(I)=AA
         ENDDO ! I

         j=0
         IF(NIPPN(IP2+1) .GT. 0) THEN
         DO I=1,NIPPN(IP2+1)
         IF(TEMP(IPPN(IP2+1,I)) .EQ. 0.) CYCLE
         j=j+1
!-----------------------------------------------------------------------------------------
! add columns
         nnz = nnz+1
         irow_vec1 = UF1_nrow_col(2,NIPPN(IP2+1)+IP3) + UF1_nrow_col(1,NIPPN(IP2+1)+IP3) 
         UF1_nrow_col(1,NIPPN(IP2+1)+IP3) = UF1_nrow_col(1,NIPPN(IP2+1)+IP3) + 1         
         
         UFNZ1%row(irow_vec1) = I				 
! NIPPN(IP2+1)+IP3 - column number
         UFNZ1%col(irow_vec1) = NIPPN(IP2+1)+IP3
         UFNZ1%val(irow_vec1) = TEMP(IPPN(IP2+1,I))        
!write(*,*) 'j=',NIPPN(IP2+1)+IP3,'  i=',i,'  irow_vec1=',irow_vec1,  &
!'  col=',UFNZ1%col(irow_vec1),'  row=',UFNZ1%row(irow_vec1),'  val=',UFNZ1%val(irow_vec1)

!-----------------------------------------------------------------------------------------
! add rows to existing NIPPN(IP2+1) columns
         nnz = nnz+1         
         irow_vec1 = UF1_nrow_col(2,I) + UF1_nrow_col(1,I)          
         UF1_nrow_col(1,I) = UF1_nrow_col(1,I) + 1
! NIPPN(IP2+1)+IP3 - line number for already existing columns except NIPPN(IP2+1)+IP3 line
         UFNZ1%row(irow_vec1) = NIPPN(IP2+1)+IP3				 
         UFNZ1%col(irow_vec1) = I
         UFNZ1%val(irow_vec1) = TEMP(IPPN(IP2+1,I))   
!write(*,*) 'i=',NIPPN(IP2+1)+IP3,'  j=',i,'  irow_vec1=',irow_vec1,  &
!'  col=',UFNZ1%col(irow_vec1),'  row=',UFNZ1%row(irow_vec1),'  val=',UFNZ1%val(irow_vec1)

!write(*,*) 'irow_vec1=',irow_vec1
!write(*,*) '--------------------------'
!-----------------------------------------------------------------------------------------
!           UF1(I,NIPPN(IP2+1)+IP3)=AATEM(IPPN(IP2+1,I))
!           UF1(NIPPN(IP2+1)+IP3,I)=AATEM(IPPN(IP2+1,I))
!if(NIPPN(IP2+1)+IP3 .eq. 365) &
!write(*,*) '#=',j,'  j=',NIPPN(IP2+1)+IP3,'  i=',i,'  irow_vec1=',irow_vec1,  &
!'  col=',UFNZ1%col(irow_vec1),'  row=',UFNZ1%row(irow_vec1),'  val=',UFNZ1%val(irow_vec1)
         ENDDO ! I
         ENDIF ! NIPPN(IP2+1).GT.0 

!stop
goto 333
! ** CHECK
      i=count(UFNZ1%row(:) .ne. 0 )
      j=count(UFNZ1%val(:) .ne. 0.) 
      write(*,'(3(a,i8))') '2: i=',i,'  j=',j,'  nnz=',nnz 

		do j=1,NIPPN(IP2+1)+IP2
		  if(UF1_nrow_col(1,j) .eq. 0) then
			write(*,'(2(a,i8))') 'j=',j,'  UF1_nrow_col(1,j)=',UF1_nrow_col(1,j)
			stop 'stop in subroutine UF1_nonzero'     
		  endif
          I1=0
          do i=1,UF1_nrow_col(1,j)
         irow_vec1 = UF1_nrow_col(2,j) + i -1
         !write(*,*) j,i,irow_vec1,                       &
         !'  UFNZ%row(irow_vec1)=',UFNZ1%row(irow_vec1),  &
			!'  UFNZ%col(irow_vec1)=',UFNZ1%col(irow_vec1),  &
			!'  UFNZ%val(irow_vec1)=',UFNZ1%val(irow_vec1)


            if(UFNZ1%row(i) .eq. 0) cycle
			I2=UFNZ1%row(i)
            if(I1 .eq. 0) then
              I1=I2
			  cycle
			endif
            if(I1 .lt. I2) cycle
         write(*,*) j,i,'  - j,i'
			write(*,'(3(a,i8),a,e14.4)')      &
			  'i1=',I1,                       &
			'  UFNZ%row(i1)=',UFNZ1%row(I1),  &
			'  UFNZ%col(i1)=',UFNZ1%col(I1),  &
			'  UFNZ%val(i1)=',UFNZ1%val(I1)
			write(*,'(3(a,i8),a,e14.4)')      &
			  'i2=',I2,                       &
			'  UFNZ%row(i2)=',UFNZ1%row(I2),  &
			'  UFNZ%col(i2)=',UFNZ1%col(I2),  &
			'  UFNZ%val(i2)=',UFNZ1%val(I2)
         
			!stop 'stop in subroutine UF1_nonzero'
            I1=I2
		  enddo ! i
		enddo ! j
write(*,*) '2: after CHECK'
stop
333 continue
		   DO IP4=1,IP2
         AA=0.0
         DO I1=1,NIPPN(IP4)
         AA=AA+TEMP(IPPN(IP4,I1))
         ENDDO ! I1
         IF(AA .EQ. 0.) CYCLE
         
!-----------------------------------------------------------------------------------------   
         nnz=nnz+1
         irow_vec1 = UF1_nrow_col(2,IIPP+IP3) + UF1_nrow_col(1,IIPP+IP3) 
         UF1_nrow_col(1,IIPP+IP3) = UF1_nrow_col(1,IIPP+IP3) + 1 
         
         UFNZ1%row(irow_vec1) = IIPP+IP4				 
         UFNZ1%col(irow_vec1) = IIPP+IP3
         UFNZ1%val(irow_vec1) = AA        
!write(*,*) 'irow_vec1=',irow_vec1

!write(*,*) IP3,IIPP+IP4,IIPP+IP3,KN1,'  - IP3,IIPP+IP4,IIPP+IP3,KN1'
!-----------------------------------------------------------------------------------------

!           UF1(IIPP+IP4,IIPP+IP3)=AA

         ENDDO ! IP4
!
      ENDDO ! IP3=1,IP2

!do j=1,IIPP+IP2
!   write(*,*) 'j=',j,'  nz=',UF1_nrow_col(1,j)
!enddo 
!stop
      if(lfirst) then
! ** CHECK
goto 444

      i=count(UFNZ1%row(:) .ne. 0 )
      j=count(UFNZ1%val(:) .ne. 0.) 
      write(*,'(3(a,i8))') '3: i=',i,'  j=',j,'  nnz=',nnz 

		do j=1,NIPPN(IP2+1)+IP2
		  if(UF1_nrow_col(1,j) .eq. 0) then
			write(*,'(2(a,i8))') 'j=',j,'  UF1_nrow_col(1,j)=',UF1_nrow_col(1,j)
			stop 'stop in subroutine UF1_nonzero'     
		  endif
          I1=0
          do i=1,UF1_nrow_col(1,j)
            if(UFNZ1%row(i) .eq. 0) cycle
			I2=UFNZ1%row(i)
            if(I1 .eq. 0) then
              I1=I2
			  cycle
			endif
            if(I1 .lt. I2) cycle
         write(*,*) j,i,'  - j,i'
			write(*,'(3(a,i8),a,e14.4)')      &
			  'i1=',I1,                       &
			'  UFNZ%row(i1)=',UFNZ1%row(I1),  &
			'  UFNZ%col(i1)=',UFNZ1%col(I1),  &
			'  UFNZ%val(i1)=',UFNZ1%val(I1)
			write(*,'(3(a,i8),a,e14.4)')      &
			  'i2=',I2,                       &
			'  UFNZ%row(i2)=',UFNZ1%row(I2),  &
			'  UFNZ%col(i2)=',UFNZ1%col(I2),  &
			'  UFNZ%val(i2)=',UFNZ1%val(I2)
			stop 'stop in subroutine UF1_nonzero'
            I1=I2
		  enddo ! i
		enddo ! j
write(*,*) '3: after CHECK'
444 continue

      do i=1,NIPPN(IP2+1)+IP2
         if(UF1_nrow_col(1,i) .eq. 0) then
          write(tmp_message,'(2(a,i0,3x))') &
          'i = ',i,'UF1_nrow_col(1,i) = ',UF1_nrow_col(1,i)
          G_ERROR(trim(tmp_message))
			endif
      enddo ! i
      nnz=sum(UF1_nrow_col(1,1:NIPPN(IP2+1)+IP2))
!		  write(*,*) 'nnz=',nnz,'  in UF1_nonzero'   
      if(nnz.gt.nnz_par) then
         write(tmp_message,'(2(a,i0,3x),a)') &
         'nnz = ',nnz,'nnz_par = ',nnz_par,'nnz > nnz_par'
         G_ERROR(trim(tmp_message))
      endif
      i=count(UFNZ1%row(:) .ne. 0 )
      j=count(UFNZ1%val(:) .ne. 0.) 
      if(i .ne. j .or. i .ne. nnz)   then
         write(tmp_message,'(3(a,i0,3x),9a)') 'i = ',i,'j = ',j,'nnz = ',nnz,'i.ne.j.or.i.ne.nnz', &
         NEW_LINE('A'), &
         'Memory distribution in preparation of UF1 needs to be modified.', &
         NEW_LINE('A'), &
         'Temporary solution:', &
         NEW_LINE('A'), &
         'change TTEST=100.,XTEST=100.,YTEST=100. patameters in INPUT', &
         NEW_LINE('A'), &
         'in order not to create groups of similar parameters for different pixeles.'
         G_ERROR(trim(tmp_message))
      endif ! i .ne. j .or.
      lfirst=.false.
      endif ! lfirst			
			
! ** Compress zeros

      do j=IIPP+IP2,2,-1
         ii1=UF1_nrow_col(2,j)
         ii2=UF1_nrow_col(2,j)+sum(UF1_nrow_col(1,j:IIPP+IP2))-1  
         i1 =UF1_nrow_col(2,j-1)+UF1_nrow_col(1,j-1)
         i2 =i1+sum(UF1_nrow_col(1,j:IIPP+IP2))-1
            !write(*,*) 'j=',j,'  ii2-ii1=',ii2-ii1,'  i2-i1=',i2-i1
              		  
         UFNZ1%row(i1:i2) = UFNZ1%row(ii1:ii2)
         UFNZ1%col(i1:i2) = UFNZ1%col(ii1:ii2)				 
         UFNZ1%val(i1:i2) = UFNZ1%val(ii1:ii2)			
      enddo ! j

      ii1=nnz+1
		ii2=(KN1-1)*KN1+UF1_nrow_col(1,KN1)         
		UFNZ1%row(ii1:ii2) = 0				 
		UFNZ1%col(ii1:ii2) = 0
		UFNZ1%val(ii1:ii2) = 0.0

! ** CHECK

      if(lfirst) then
         i=count(UFNZ1%row(nnz+1:KN1*KN1).ne.0 ) 
         if(i .ne. 0) write(*,*) '++++  i=',i

         i=count(UFNZ1%row(1:nnz) .ne. 0 )
         j=count(UFNZ1%val(1:nnz) .ne. 0.)
         if(i .ne. j .or. i .ne. nnz)   then
         write(tmp_message,'(3(a,i0))') &
         'in compress_zero_rows:  i = ',i,  &
		                     'j = ',j,'nnz = ',nnz
         G_ERROR(trim(tmp_message))
         endif
         do j=1,NIPPN(IP2+1)+IP2
         do i=2,UF1_nrow_col(1,j)
            i1=sum(UF1_nrow_col(1,1:j))-UF1_nrow_col(1,j)+i
            if(UFNZ1%row(i1) .gt. UFNZ1%row(i1-1)) cycle
            write(*,'(9(a,i8))') 'j=',j,'  ivec-1=',i1-1,'  ivec=',i1,    &
                              '  UFNZ%row(i-1)=',UFNZ1%row(i1-1),         &
                              '  UFNZ%col(i-1)=',UFNZ1%col(i1-1),         &
                              '  UFNZ%row(i)=',UFNZ1%row(i1),             &
                              '  UFNZ%col(i)=',UFNZ1%col(i1),             &
                              '  nrow_col(j)=',UF1_nrow_col(1,j),         &
                              '  ncol_mtr=',KN1						  

            if(i1 .lt. UF1_nrow_col(1,j))  &
            write(*,'(a,i0)') '  UFNZ1%row(i+1)=',UFNZ1%row(i1+1)
            write(*,'(a)') 'compress_zero_rows: UFNZ1%row(i-1).ge.UFNZ1%row(i)'
            G_ERROR(trim(tmp_message))
         enddo ! i
         enddo ! j	
      lfirst=.false.
		endif ! lfirst						

	   
! **  Fill up nval - place (pointer?) of the first nonzero value in each column in vector UFNZ%val 

      UFNZ1%nval(:) = 0
      do j=1,NIPPN(IP2+1)+IP2				 
      UFNZ1%nval(j) = sum(UF1_nrow_col(1,1:j))-UF1_nrow_col(1,j)+1          	
      enddo ! i
      UFNZ1%nval(NIPPN(IP2+1)+IP2+1) = nnz+1

      return
      end subroutine UF1_nonzero

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

      subroutine UF_nonzero (                         &
                              KNSINGF,npixels,UFS,    & ! IN
                              TCCORG,                 &
                              smim,                   &
                              IP2,                    &
                              UFNZ,                   & ! INOUT
                              nnz                     & ! OUT
                            )

      use mod_par_inv, only : KIMAGE,KPARS,KPAR
      use mod_index_cloud
	  	  
      implicit none
 
! ----------------------------------------------------------------------- 
      integer,                             intent(in)  :: KNSINGF,npixels,IP2
      real, dimension(KPARS,KPARS,KIMAGE), intent(in)  :: UFS
      real,                                intent(in)  :: TCCORG
      real, dimension(KIMAGE,KPARS,KIMAGE),intent(in)  :: smim

      type(nonzero),                       intent(inout) :: UFNZ
      integer,                             intent(out)   :: nnz
! -----------------------------------------------------------------------      
      real                                 :: XA
      integer                              :: I,J,IS,IS1,IS2, &  
                                              ipix,jpix,KNF,KN1, &
                                              icol,irow_mtr,irow_vec  											   
      integer, dimension(KPAR)             :: UF_nrow_col
      logical                              :: lform
      logical, save                        :: lfirst=.true.
! -----------------------------------------------------------------------
   
      if(IP2 .EQ. 0) then
         lform=.true.
      else
		   lform=.false.
      endif ! IP2 .EQ. 0

      KN1=KNSINGF+npixels-1
      KNF =npixels*KNSINGF
	  	  
      UF_nrow_col(:)=0

      UFNZ%row(:)=0
      UFNZ%col(:)=0
      UFNZ%val(:)=0.0	

!  J  jpix   - UF matrix columns 
!  I  ipix   - UF matrix rows
      nnz=0
      do jpix=1,npixels
      do J=1,KNSINGF
      icol=(jpix-1)*KNSINGF+J
      do ipix=1,npixels
         if(ipix .eq. jpix) then
         do I=1,KNSINGF
            XA=UFS(I,J,ipix)
		      if(XA .eq. 0.) cycle
            nnz=nnz+1
            irow_mtr=(ipix-1)*KNSINGF+I
            if(lform) then
            irow_vec=nnz
            else
            irow_vec=ivec1 (                     &
                              KNSINGF,npixels,   &
                              ipix,ipix,I,J      &
                           ) 
            endif ! lform
            UFNZ%col(irow_vec) = icol
            UFNZ%row(irow_vec) = irow_mtr				 
            UFNZ%val(irow_vec) = XA             
            UF_nrow_col(icol)=UF_nrow_col(icol)+1
            if(irow_mtr .eq. 0) then
            write(tmp_message,'(a,i0,2x,a)') &
            '1: irow_mtr = ',irow_mtr,'can not be 0.'
            G_ERROR(trim(tmp_message))
            endif
!                 UF((IIMAGEQ1-1)*KNSINGF+I,(IIMAGEQ1-1)*KNSINGF+I1)=  &
!                 UFS(IIMAGEQ1,I,I1)
         enddo ! I
         else
            if(ipix .lt. jpix) then
            XA=TCCORG*smim(jpix,J,ipix)
            else 
            XA=TCCORG*smim(ipix,J,jpix)
            endif ! ipix.lt.jpix
		      if(XA .eq. 0.) cycle
            nnz=nnz+1          
            irow_mtr=(ipix-1)*KNSINGF+J
            if(lform) then
            irow_vec=nnz
            else
            irow_vec=ivec1 (                     &
                              KNSINGF,npixels,   &
                              ipix,jpix,J,J      &
                           ) 
            endif ! lform
            UFNZ%col(irow_vec) = icol
            UFNZ%row(irow_vec) = irow_mtr				 
            UFNZ%val(irow_vec) = XA
            UF_nrow_col(icol)=UF_nrow_col(icol)+1
            if(irow_mtr .eq. 0) then
            write(tmp_message,'(a,i0,2x,a)') &
            '2: irow_mtr=',irow_mtr,'can not be 0.'
            G_ERROR(trim(tmp_message))
            endif
!               UF((IMSM(IS,IS1)-1)*KNSINGF+I,(IMSM(IS,IS2)-1)*KNSINGF+I)=  & 
!               TCCORG*SMIM(IS,I,IS1,IS2)
		    endif ! ipix .eq. jpix
      enddo ! ipix
      if(lform) UFNZ%nval(icol) = sum(UF_nrow_col(1:icol))-UF_nrow_col(icol)+1          	
      enddo ! J
      enddo ! jpix  

      if(lform) UFNZ%nval(KNF+1) = nnz+1

! ** CHECK

      if(lfirst .and. (.not. lform)) then
		do j=1,KNF
         if(UF_nrow_col(j) .eq. 0) then
            write(tmp_message,'(2(a,i0,3x))') &
            'j = ',j,'  UF_nrow_col(j) = ',UF_nrow_col(j)
            G_ERROR(trim(tmp_message))
         endif ! UF_nrow_col(j) .eq. 0
         IS1=0
         do i=1,KN1
            if(UFNZ%row(i) .eq. 0) cycle
            IS2=UFNZ%row(i)
            if(IS1 .eq. 0) then
            IS1=IS2
            cycle
            endif ! IS1 .eq. 0
            if(IS1 .lt. IS2) cycle
            write(tmp_message,'(3(a,i0),a,es11.4,a,3(a,i0),a,es11.4)') &
            'i1=',IS1,                          &
            '  UFNZ%row(i1) = ',UFNZ%row(IS1),  &
            '  UFNZ%col(i1) = ',UFNZ%col(IS1),  &
            '  UFNZ%val(i1) = ',UFNZ%val(IS1),  &
            NEW_LINE('A'), &
            'i2=',IS2,                          &
            '  UFNZ%row(i1) = ',UFNZ%row(IS2),  &
            '  UFNZ%col(i1) = ',UFNZ%col(IS2),  &
            '  UFNZ%val(i1) = ',UFNZ%val(IS2)
            G_ERROR(trim(tmp_message))
            IS1=IS2
         enddo ! i
      enddo ! j

      nnz=sum(UF_nrow_col(1:KNF))
!      write(*,*) 'nnz=',nnz,'  in UF_nonzero'
!      write(*,*) '1: nrow_col :'
!      write(*,'(10i8)') (UF_nrow_col(j),j=1,KNF)

      if(nnz .gt. nnz_par) then
         write(tmp_message,'(2(a,i0,3x),a)') &
         'nnz = ',nnz,'nnz_par = ',nnz_par,'nnz > nnz_par'
         G_ERROR(trim(tmp_message))
      endif
		  i=count(UFNZ%row(1:KNF*KN1) .ne. 0 )
		  j=count(UFNZ%val(1:KNF*KN1) .ne. 0.)
      if(i .ne. j .or. i .ne. nnz)   then
         write(tmp_message,'(3(a,i0,3x),2a)') 'i = ',i,'j = ',j,'nnz =',nnz, &
         NEW_LINE('A'), &
         'i .ne. j .or. i .ne. nnz'
         G_ERROR(trim(tmp_message))
		  endif
      lfirst=.false.
      endif ! lfirst						
	  
      return
      end subroutine UF_nonzero

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! irow_vec calculation using irow_mtr,icol_mtr

!   			             icol=(jpix-1)*KNSINGF+J

!                        irow_mtr=(ipix-1)*KNSINGF+I

!  ipix.eq.jpix      irow_vec=(icol-1)*(KNSINGF+npixels-1)+ipix-1+I

!  ipix.LT.jpix
!					          irow_vec=(icol-1)*(KNSINGF+npixels-1)+ipix    

!  ipix.GT.jpix 
!					          irow_vec=(icol-1)*(KNSINGF+npixels-1)+ipix-1+KNSINGF

      integer function ivec2 (                    &
                              KNSINGF,npixels,     &
                              irow_mtr,icol_mtr   &
                             ) 
 
! icol_mtr - number of matrix column
! irow_mtr - number of matrix rows	  	  	  

      implicit none
! -----------------------------------------------------------------------
      integer,  intent(in)   :: KNSINGF,npixels,irow_mtr,icol_mtr	  
! -----------------------------------------------------------------------
      integer                :: ipix,jpix,I,J
! -----------------------------------------------------------------------

			
      jpix=icol_mtr/KNSINGF
      J=icol_mtr-jpix*KNSINGF 
      if(J .GT. 0) jpix=jpix+1  
      if(J .EQ. 0)      J=KNSINGF 

		  ipix=irow_mtr/KNSINGF
      I=irow_mtr-ipix*KNSINGF
		  if(I .gt. 0) ipix=ipix+1  
		  if(I .eq. 0) I=KNSINGF  
		   			
      if(ipix .ne. jpix .and. I .ne. J) then
         ivec2=-999
         return
		  endif

      ivec2=ivec1 (                     &
                     KNSINGF,npixels,    &
                     ipix,jpix,I,J      &
                  )

!        write(*,'(8(a,i3))') 'in function ivec2:  jpix=',jpix,'  J=',J, &
!						    '  ipix=',ipix,'  I=',I,                     &
!	                        '  icol=',icol,'  irow=',ivec2    


      end function ivec2

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

! calculate irow_vec using ipix,jpix,I,J

      integer function ivec1 (                    &
                              KNSINGF,npixels,     &
                              ipix,jpix,I,J       &
                             ) 

! matrix column
! icol - number of matrix column
! for vector row calculation
! ipix - number of pixel
! I      - number of retrieved parameter	  	  	  

      implicit none
! -----------------------------------------------------------------------
      integer,  intent(in)  :: KNSINGF,npixels,ipix,jpix,I,J	  
! -----------------------------------------------------------------------
      integer :: icol_mtr
! -----------------------------------------------------------------------

! VECTOR

!  number of columns      npixels*KNSINGF

!  number of rows         KNSINGF+(npixels-1)

      icol_mtr=(jpix-1)*KNSINGF+J

      if(ipix .eq. jpix)  &
      ivec1=(icol_mtr-1)*(KNSINGF+npixels-1)+(ipix-1)+I

      if(ipix .lt. jpix)  &
      ivec1=(icol_mtr-1)*(KNSINGF+npixels-1)+ipix
	      
      if(ipix .gt. jpix)  &
      ivec1=(icol_mtr-1)*(KNSINGF+npixels-1)+(ipix-1)+KNSINGF


!      write(*,'(8(a,i3))') 'in function ivec1:  jpix=',jpix,'  J=',J, &
!						  '  ipix=',ipix,'  I=',I,                     &
!	                      '  icol=',icol_mtrx,'  irow_vec=',ivec1    

		   			
      end function ivec1

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! subroutine is not in use
      subroutine compress_zeros (                     &
                                  KNSINGF,npixels,    & ! IN
                                  nnz,                &
                                  UFNZ                & ! INOUT
                                )

      use mod_par_inv, only : KPAR
      use mod_index_cloud

      implicit none
! -----------------------------------------------------------------------
      integer,                  intent(in)    :: npixels,KNSINGF,nnz
      type(nonzero),            intent(inout) :: UFNZ
! -----------------------------------------------------------------------
      integer, dimension(KPAR)  :: nrow_col                   

      integer                   :: ncol_mtr,nrow_mtr,nel_vec,  &	        
                                   i,j,i1,i2,ii1,ii2  
! -----------------------------------------------------------------------
      UFNZ%nval(:) = 0

! ***  1  Shift nonzero rows up

! UF
      ncol_mtr=npixels*KNSINGF
      nrow_mtr=KNSINGF+npixels-1
      nel_vec=nrow_mtr*ncol_mtr
write(*,*) 'ncol_mtr=',ncol_mtr,'  nrow_mtr=',nrow_mtr,'  npixels=',npixels,'  KNSINGF=',KNSINGF,'  nel_vec=',nel_vec
      
      ii2 = nel_vec
      i2  = nel_vec - 1
      do i=nel_vec,2,-1
!? can be improved; do not move '0' elements at the end of vector
      if(UFNZ%row(i) .ne. 0 .and. UFNZ%row(i-1) .eq. 0) then
         ii1 = i
         i1  = i-1
         UFNZ%col(i1:i2) = UFNZ%col(ii1:ii2)
         UFNZ%row(i1:i2) = UFNZ%row(ii1:ii2)				 
         UFNZ%val(i1:i2) = UFNZ%val(ii1:ii2)
         UFNZ%row(ii2)   = 0				 
         UFNZ%col(ii2)   = 0
         UFNZ%val(ii2)   = 0.0
         ii2 = ii2 - 1
         i2  = i2  - 1
      endif
      enddo ! i
      
      i=count(UFNZ%row(1:nel_vec-1) .ne. 0 )
      j=count(UFNZ%val(1:nel_vec-1) .ne. 0.) 
      if(i .ne. j)   then
         write(*,*) 'i=',i,'  j=',j
         stop 'stop 1 in compress_zeros: i .ne. j'
      endif
      if(i .ne. nnz) then
         write(*,*) 'i=',i,'  nnz=',nnz
         stop 'stop 1 in compress_zeros: i .ne. nnz'
      endif

	   
! ***  2  Fill out nval - place of first value in each column in vector UFNZ%val 
!      write(*,*) '2: nrow_col :'
!      write(*,'(10i8)') (count(UFNZ%col(1:nel_vec) .eq. j ),j=1,ncol_mtr)
      
      nrow_col(:) = 0
      do j=1,ncol_mtr
      nrow_col(j) = count(UFNZ%col(1:nel_vec) .eq. j )
      enddo ! j
      
      i1=0     
      UFNZ%nval(1)=1
      do j=2,ncol_mtr				 
         if(nrow_col(j-1) .eq. 0) cycle 
         i1=i1+1
         i=UFNZ%nval(j-1)+nrow_col(j-1)
         UFNZ%nval(j) = i		
      enddo ! i
      UFNZ%nval(i1+1) = nnz+1
	  
      return
      end subroutine compress_zeros

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! sparse_solver
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine sparse_solver ( iu_main_output, & ! IN
                                 KN,nnz,UFNZ,       & 
                                 b,                 & ! INOUT
                                 solver_timer,      & 
                                 npixels,KNSINGF    & ! OPTIONAL
                               )
         use mod_par_inv, only : KPAR

         implicit none
! -----------------------------------------------------------------------
         integer,                 intent(in)    :: iu_main_output,  &
                                                   KN,nnz
         type(nonzero),           intent(in)    :: UFNZ
         real*8, dimension(KPAR), intent(inout) :: b
         real,                    intent(inout) :: solver_timer
         integer,optional,        intent(in)    :: npixels,KNSINGF
! -----------------------------------------------------------------------
         real*4                        ::  time_start,time_finish
! -----------------------------------------------------------------------
      call CPU_TIME(time_start)

      if(present(npixels)) then
#if defined(SUPERLU_MT)
         !write(iu_main_output,*) 'solver_SuperLU_MP'
         write(*,*) 'Development needed for solver_SuperLU_MT:'
         write(*,*) 'error estimations after multipixel retrieval'
         stop 'stop in sparse_solver'
         call solver_SuperLU_MT(KN,nnz,UFNZ,b,npixels,KNSINGF)
#elif defined(SUPERLU_43)
         !write(iu_main_output,*) 'solver_SuperLU_43'
         call solver_SuperLU_43(KN,nnz,UFNZ,b,npixels,KNSINGF)
#elif defined(VIENNA_CL)
         !write(iu_main_output,*) 'solver_ViennaCL'
         write(*,*) 'Development needed for solver_ViennaCL'
         write(*,*) 'error estimations after multipixel retrieval'
         stop 'stop in sparse_solver'
         call solver_viennacl(KN,nnz,UFNZ,b,npixels,KNSINGF)
#elif defined(MUMPS)
         !write(iu_main_output,*) 'solver_MUMPS'
         write(*,*) 'Development needed for solver_MUMPS'
         write(*,*) 'error estimations after multipixel retrieval'
         stop 'stop in sparse_solver'
         call solver_MUMPS(KN,nnz,UFNZ,b,npixels,KNSINGF)
#else
#error No SPARSE_SOLVER configured!
#endif
      else
#if defined(SUPERLU_MT)
         !write(iu_main_output,*) 'solver_SuperLU_MP'
         call solver_SuperLU_MT(KN,nnz,UFNZ,b)
         if ( stop_report%status ) return
#elif defined(SUPERLU_43)
         !write(iu_main_output,*) 'solver_SuperLU_43'
         call solver_SuperLU_43(KN,nnz,UFNZ,b)
         if ( stop_report%status ) return
#elif defined(VIENNA_CL)
         !write(iu_main_output,*) 'solver_ViennaCL'
         call solver_viennacl(KN,nnz,UFNZ,b)
         if ( stop_report%status ) return
#elif defined(MUMPS)
         !write(iu_main_output,*) 'solver_MUMPS'
         call solver_MUMPS(KN,nnz,UFNZ,b,npixels)
         if ( stop_report%status ) return
#else
#error No SPARSE_SOLVER configured!
#endif      
      endif ! present(npixels) .and. present(KNSINGF)

      call CPU_TIME(time_finish) 
      solver_timer = solver_timer + (time_finish-time_start)

        return
    end subroutine sparse_solver
	
#ifdef SUPERLU_43
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! solver_SuperLU
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine solver_SuperLU_43 ( KN,nnz,UFNZ,       & ! IN
                                     b,                 & ! INOUT
                                     npixels,KNSINGF    & ! OPTIONAL
                                    )

      use mod_par_inv, only : KPAR
!      use mod_index_cloud

!     iso_c_binding is only needed to get the type integer(kind=c_intptr_t)
!     integer(kind=c_intptr_t) is the official integer type to hold pointers
!     in fortran 2003+.
!     When iso_c_binding is not available (especially with old compilers),
!     one can remove the 'use iso_c_binding' statement. One has then to replace
!     the type integer(kind=c_intptr_t) with;
!     * integer*4 on 32-bits configurations (32-bits pointers)
!     * integer*8 on 64-bits configurations (64-bits pointers)
      use iso_c_binding ! c_intptr_t
	  	  
      implicit none
 
!   ----------------------------------------------------------------------- 
      integer,                  intent(in)    :: KN,nnz
      type(nonzero),            intent(in)    :: UFNZ
      real*8, dimension(KPAR),  intent(inout) :: b
      integer,optional,         intent(in)    :: npixels,KNSINGF
!	------------------------------------------------------------------------------------------
!	------------------------------------------------------------------------------------------
! PhF
      integer                   ::  n, nrhs, ldb, info, iopt
!      integer*8  ::  factors
! FD
      integer(kind=c_intptr_t)  ::  factors
!	------------------------------------------------------------------------------------------
      real*8,dimension(:),allocatable  ::  b1
      integer                          ::  ipix,i,j,i1,ierr
!	------------------------------------------------------------------------------------------
      nrhs = 1
      n    = KN
      ldb  = KN
!*
!* First, factorize the matrix. The factors are stored in *factors* handle.     
      iopt = 1
      call c_fortran_dgssv (                                                       &
                              iopt, n, nnz, nrhs,                                  &
                              UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                              b(1:KN), ldb, factors, info                          &
                           )
!*
      if (info .eq. 0) then
!tl         write (*,*) 'Factorization succeeded'
      else
         write(tmp_message,'(a,i0)') 'INFO from factorization = ',info
         G_ERROR(trim(tmp_message))
      endif

      if(present(npixels)) then
      ! This part is used for error estimations
        allocate(b1(1:KPAR),stat=ierr)
        if (ierr /= 0) then
          write(tmp_message,'(a)') &
          'error while trying to allocate b1(1:KPAR) in solver_SuperLU_43'
          G_ERROR(trim(tmp_message))
        endif
        do ipix=1,npixels
        do i=1,KNSINGF
          j = KNSINGF*(ipix-1)+i
          b1(:) = 0.0d0
          b1(j) = 1.0d0
!*
!* Second, solve the system using the existing factors.
          iopt = 2
          call c_fortran_dgssv (                                                       &
                                  iopt, n, nnz, nrhs,                                  &
                                  UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                                  b1(1:KN), ldb, factors, info                         &
                               )
!*
          if (info .eq. 0) then
          !write (*,*) 'SuperLU:  Solve succeeded'
!tl         write (*,*) (b(i1), i1=1,10)
!           write (*,'(10e14.4)') (b(i1), i1=1,n )
          else
            write(tmp_message,'(a,i0)') 'INFO from triangular solve = ',info
            G_ERROR(trim(tmp_message))
          endif
          b(j) = b1(j)
        enddo ! i
        enddo ! ipix
        deallocate(b1,stat=ierr)
        if (ierr /= 0) then
          write(tmp_message,'(a)') &
          'error while trying to deallocate b1(1:KPAR) in solver_SuperLU_43'
          G_ERROR(trim(tmp_message))
        endif
      else
!*
!* Second, solve the system using the existing factors.
        iopt = 2
        call c_fortran_dgssv (                                                       &
                                iopt, n, nnz, nrhs,                                  &
                                UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                                b(1:KN), ldb, factors, info                          &
                             )
!*
        if (info .eq. 0) then
          !write (*,*) 'SuperLU:  Solve succeeded'
!tl         write (*,*) (b(i1), i1=1,10)
!           write (*,'(10e14.4)') (b(i1), i1=1,n )
        else
          write(tmp_message,'(a,i0)') 'INFO from triangular solve = ',info
          G_ERROR(trim(tmp_message))
        endif
      endif ! present(npixels) .and. present(KNSINGF)
!*
!* Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_fortran_dgssv (                                                       &
                              iopt, n, nnz, nrhs,                                  &
                              UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                              b(1:KN), ldb, factors, info                          &
						   )
!*      
      return
      end subroutine solver_SuperLU_43
#endif
#ifdef SUPERLU_MT
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! solver_SuperLU_MT
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine solver_SuperLU_MT (  KN,nnz,UFNZ,     & ! IN
                                      b,               & ! INOUT
                                      npixels,KNSINGF  & ! OPTIONAL
                                   )
      use mod_par_inv, only : KPAR
      use iso_c_binding ! c_intptr_t
	  	  
      implicit none
 !   ----------------------------------------------------------------------- 
      integer,                  intent(in)    :: KN,nnz
      type(nonzero),            intent(in)    :: UFNZ
      real*8, dimension(KPAR),  intent(inout) :: b
      integer,optional,         intent(in)    :: npixels,KNSINGF
!	------------------------------------------------------------------------------------------
      integer n, nrhs, ldb, info, iopt
      integer(kind=c_intptr_t) factors
!	------------------------------------------------------------------------------------------
!   ind <= 1 - factorize the matrix 
!        > 1 - use previous factorization  
!	------------------------------------------------------------------------------------------

      nrhs = 1
      n    = KN
      ldb  = KN
	  
! Unlike the regular SuperLU 4.3, the mt version does not separate factorization
! from solving. both is done at once.

!*
!* First, factorize & solve the matrix. The factors are stored in *factors* handle.
      iopt = 1
      call c_bridge_pdgssv ( 6,                                                   &
                             iopt, n, nnz, nrhs,                                  &
	                           UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                             b(1:KN), ldb, factors, info                          &
						               )
!*
      if (info .eq. 0) then
!tl         write (*,*) 'Factorization & solve succeeded'
      else
         write(tmp_message,'(a,i0)') 'INFO from factorization & solve = ',info
         G_ERROR(trim(tmp_message))
      endif

!* Last, free the storage allocated inside SuperLU
      iopt = 3
      call c_bridge_pdgssv ( 6,                                                   &
	                           iopt, n, nnz, nrhs,                                  &
	                           UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                             b(1:KN), ldb, factors, info                          &
						               )
!*
      return
      end subroutine solver_SuperLU_MT
#endif
#ifdef VIENNA_CL
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! solver_viennacl
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine solver_viennacl (  KN,nnz,UFNZ,     & ! IN
                                    b,               & ! INOUT
                                    npixels,KNSINGF  & ! OPTIONAL
                                 )

        use mod_par_inv, only : KPAR
        implicit none
					  		
		interface 
			subroutine viennacl_solve(m, nnz, vals, rowidx, colptr, b, info) bind (C, name="viennacl_solve")
			use, intrinsic :: ISO_C_BINDING
			implicit none
			
			integer(C_INT), value			:: m, nnz
			real(C_DOUBLE), dimension(nnz)	:: vals
			integer(C_INT), dimension(nnz)	:: rowidx
			integer(C_INT), dimension(m+1)	:: colptr
			real(C_DOUBLE), dimension(m)	:: b
			integer(C_INT) 					:: info
			
			end subroutine viennacl_solve
		end interface

!   ----------------------------------------------------------------------- 
        integer,                  intent(in)    :: KN,nnz
        type(nonzero),            intent(in)    :: UFNZ
        real*8, dimension(KPAR),  intent(inout) :: b
        integer,optional,         intent(in)    :: npixels,KNSINGF
!	------------------------------------------------------------------------------------------
        integer	:: info = 0
!	------------------------------------------------------------------------------------------
        call viennacl_solve (                                                      &
                              KN, nnz,                                             &
                              UFNZ%val(1:nnz), UFNZ%row(1:nnz), UFNZ%nval(1:KN+1), &
                              b, info                                              &
                            )
        if (info .ne. 0) then
            write (tmp_message,'(a,i0)') "cusparse failed: info = ", info
            G_ERROR(trim(tmp_message))
        endif
        
        return
    end subroutine solver_viennacl
#endif
#ifdef MUMPS
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
! solver_MUMPS
! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
      subroutine solver_MUMPS ( KN,nnz,UFNZ,     & ! IN
                                b,               & ! INOUT
                                npixels,KNSINGF  & ! OPTIONAL
                              )

      use mod_par_inv, only : KPAR
      use mod_index_cloud

      implicit none

      INCLUDE 'mpif.h'
      INCLUDE 'dmumps_struc.h'	  	  

!	------------------------------------------------------------------------------------------
      TYPE (DMUMPS_STRUC) :: mumps_par
      INTEGER             :: IERR
!	------------------------------------------------------------------------------------------ 
!   ----------------------------------------------------------------------- 
      integer,                  intent(in)    :: KN,nnz
      type(nonzero),            intent(in)    :: UFNZ
      real*8, dimension(KPAR),  intent(inout) :: b
      integer,optional,         intent(in)    :: npixels,KNSINGF

!	------------------------------------------------------------------------------------------
      integer            :: i,i1,nnz_diag
      logical            :: lsym = .true.
!	------------------------------------------------------------------------------------------

!*

!C Define a communicator for the package.                                                                                                                                                                                
      mumps_par%COMM = MPI_COMM_WORLD
!C  Initialize an instance of the package                                                                                                                                                                                
!C  for L U factorization (sym = 0, with working host)                                                                                                                                                                   
      mumps_par%JOB = -1
      mumps_par%PAR =  1
      if(lsym) then
      mumps_par%SYM =  2  ! 2 - A is general symmetric                          
      else
      mumps_par%SYM =  0  ! 0 - A is unsymmetric
      endif ! lsym
	  
      CALL DMUMPS(mumps_par)

!C  Define problem on the host (processor 0)                                                                                                                                                                             
      IF ( mumps_par%MYID .eq. 0 ) THEN
         mumps_par%N    = KN		
         mumps_par%NRHS = 1 !mumps_par%N
         if(lsym) then
            nnz_diag=0
            do i=1,nnz
            if(UFNZ%col(i).eq.UFNZ%row(i)) nnz_diag=nnz_diag+1
            enddo ! i
            mumps_par%NZ   = (nnz-nnz_diag)/2+nnz_diag  
            !write(*,*) '** symmetric matrix option used'
         else
            mumps_par%NZ   = nnz  
            !write(*,*) '** unsymmetric matrix option used'
         endif ! lsym
		
         mumps_par%ICNTL(1)=0 ! 0 NO ERRORS
         mumps_par%ICNTL(2)=0 ! 0 NO WARNINGS
         mumps_par%ICNTL(3)=0 ! 0 REMOVE MUMPS MAIN VERBOSE

         ALLOCATE( mumps_par%IRN ( mumps_par%NZ ) )
         ALLOCATE( mumps_par%JCN ( mumps_par%NZ ) )
         ALLOCATE( mumps_par%A   ( mumps_par%NZ ) )

         if(lsym) then
! Lower triangle of symmetric matrix; mumps_par%SYM =  2 - A is general symmetric
            i1=0
            do i=1,nnz
            if(UFNZ%row(i).le.UFNZ%col(i)) then
            i1=i1+1
            mumps_par%IRN(i1) = UFNZ%col(i) ! CCS to MUMPS storage conversion (matrix is symmetric)  
            mumps_par%JCN(i1) = UFNZ%row(i) ! CCS to MUMPS storage conversion (matrix is symmetric)
            mumps_par%A  (i1) = UFNZ%val(i)
            endif
            enddo ! i
            !WRITE( 6, * ) ' i1=', i1 ,' mumps_par%NZ=',mumps_par%NZ,'  in solver_MUMPS'
            if(i1.ne.mumps_par%NZ) then
            write(tmp_message,'(2(a,i0,3x),2a)') &
            'i1 = ', i1 ,'mumps_par%NZ = ',mumps_par%NZ, &
            NEW_LINE('A'), &
            'i1 .ne. mumps_par%NZ'
            G_ERROR(trim(tmp_message))
            endif
         else
!  Unsymmetric matrix; mumps_par%SYM =  0 - A is unsymmetric
            mumps_par%IRN(1:nnz) = UFNZ%col(1:nnz) ! CCS to MUMPS storage conversion (matrix is symmetric)  
            mumps_par%JCN(1:nnz) = UFNZ%row(1:nnz) ! CCS to MUMPS storage conversion (matrix is symmetric)
            mumps_par%A  (1:nnz) = UFNZ%val(1:nnz)
         endif ! lsym

         ALLOCATE( mumps_par%RHS(mumps_par%N*mumps_par%NRHS))
         mumps_par%RHS(1:mumps_par%N)  = b(1:mumps_par%N)                                                                                                                                                                                                        

      END IF ! mumps_par%MYID .eq. 0

!C  Call package for solution                                                                                                                                                                                            
      mumps_par%JOB = 6
      CALL DMUMPS(mumps_par)
!C  Solution has been assembled on the host                                                                                                                                                                              
      IF ( mumps_par%MYID .eq. 0 ) THEN
         !WRITE( 6, * ) ' MUMPS:  Solve succeeded'
         !WRITE( 6, * ) ' Solution is ',(mumps_par%RHS(I),I=1,mumps_par%N*mumps_par%NRHS)
         b(1:mumps_par%N) = mumps_par%RHS(1:mumps_par%N)
      END IF ! mumps_par%MYID .eq. 0
	  
!tl	  IF ( mumps_par%MYID .eq. 0 ) THEN
!tl       b(1:mumps_par%N) = mumps_par%RHS(1:mumps_par%N)
!tl	  ENDIF

!tl	  CALL MPI_BCAST(b,KN,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)  ! not sure if it is necessary ???
	  
!C  Deallocate user data                                                                                                                                                                                                 
      IF ( mumps_par%MYID .eq. 0 )THEN
         DEALLOCATE( mumps_par%IRN )
         DEALLOCATE( mumps_par%JCN )
         DEALLOCATE( mumps_par%A   )
         DEALLOCATE( mumps_par%RHS )
      END IF ! mumps_par%MYID .eq. 0
!C  Destroy the instance (deallocate internal data structures)                                                                                                                                                           
      mumps_par%JOB = -2
      CALL DMUMPS(mumps_par)
!* 
	  
      return
      end subroutine solver_MUMPS

#endif

! sssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

end module mod_fisher_matrix_ccs





