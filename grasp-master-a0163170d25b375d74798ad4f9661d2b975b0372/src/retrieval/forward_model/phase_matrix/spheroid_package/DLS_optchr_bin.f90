! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

!c ****************************************************************

      SUBROUTINE OPTCHAR_bin(key,keyEL,keySUB,keyLS,key_RD1      &
                        ,key_RD                                  &
                        ,NWL,WL,RN,RK,KN,grid,SD                 &
                        ,KR,R,RD,KM,ANGLE,xext,xabs,xsca,albedo  &
                        ,f11,f12,f22,f33,f34,f44                 &
                        ,KC,RC2,LB,LE                            &
                        ,pomin,pomax,xnmin,xnmax,xkmin,xkmax     &                
                        ,KERNEL,XGRID1                           &
                        ,KERNELS1,KERNELS2 )                          
! RD version
!c **************************************************************** c
!c **   08/26/02                                                 ** c
!c **   Subroutine calculates optical characteristics for given  ** c
!c **   size distribution, refractive index, axis ratio          ** c
!c **   distribution and wavelength.                             ** c
!c **                                                            ** c
!c **   In case RN or RK or R() is out of correspoding ranges    ** c 
!c **   subroutine changes RN or RK or R() for edge value and    ** c
!c **   gives optical characteristics for new values             ** c 
!c **************************************************************** c
!c **                                                            ** c
!c ** INPUT:                                                     ** c
!c **                                                            ** c 
!c **   key  = 0 - calculate opt.characteristics from original   ** c
!c **              kernels 'Rke...'                              ** c
!c **          1 - read original kernels, create WL kernels      ** c
!c **              and save them in dirname_N  and calculate     ** c
!c **              opt.characteristics                           ** c
!c **          2 - read WL kernels and calculate                 ** c
!c **              opt.characteristics                           ** c
!c **   key_RD =1 - volume mixture of spheroids                  ** c
!c **           2 - surface area  mixture of spheroids           ** c
!c **   keyEL=  0 - calculate EXT                                ** c 
!c **           1 - EXT,F11                                      ** c
!c **           2 - EXT,F11,F12                                  ** c
!c **           3 - EXT,F11,F12,F22                              ** c
!c **           4 - EXT,F11,F12,F22,F33                          ** c
!c **           5 - EXT,F11,F12,F22,F33,F34                      ** c
!c **           6 - EXT,F11,F12,F22,F33,F34,F44                  ** c
!c **   key_grid1 read grid radii and scat.angles which were used** c
!c **             for kernel look up tables or fixed kernels     ** c
!c **          =0 - 'grid1.dat'                                  ** c
!c **           1 - 'grid1.dat.new'                              ** c
!c **   WL   - wavelength                                        ** c
!c **   RN   - real part of the refractive index                 ** c
!c **   RK   - imaginary part of the refractive index            ** c
!c **   KN   - number of grid radii                              ** c
!c **   grid(KN) - grid radii                                    ** c
!c **   SD(KN)   - size distribution for grid radii              ** c
!c **   distname_O - original kernel directory name              ** c
!c **   distname_F - .fix kernel directory name                  ** c
!c **   distname_N - new original kernel directory name          ** c
!c **   KR  - number of axis ratios                              ** c
!c **   R(KR)  - grid axis ratios                                ** c
!c **   RD(KR) - axis ratio distribution for grid axis ratios    ** c
!c **   KM   - number of scattering angles                       ** c
!c **   ANGLE(KM) - scattering angles                            ** c
!c **                                                            ** c
!c ** OUTPUT:                                                    ** c
!c **                                                            ** c
!c **   ext     - extinction                                     ** c
!c **   albedo  - absorption                                     ** c
!c **   f... - scattering matrix elements                        ** c
!c **************************************************************** c
      USE mod_par_DLS
      USE mod_par_DLS_bin
      use mod_intrpl_linear
      use phase_func
      use mod_type_DLS
      use mod_alloc_kernels
!  ---------------------------------------------------------------------	  	  	  
      type(kernels_triangle_bin),  intent(inout)  ::  KERNELS1
      type(kernels_lognormal_bin), intent(inout)  ::  KERNELS2
      
      integer,               intent(in)  ::  KM,KN,KC,KR,   &
	                                          LB,LE 
      real,dimension(KMpar), intent(in)  ::  ANGLE
	    integer,               intent(in)  ::  key,keyEL,keySUB,    &
	                                          keyLS,key_RD1,key_RD 
      integer,               intent(in)  ::  NWL
      real,dimension(KWLpar),intent(in)  ::  WL
      real,dimension(KWLpar),intent(in)  ::  RN,RK
      real,dimension(KNpar), intent(in)  ::  grid
      real,dimension(KCpar), intent(in)  ::  RC2
      real*4,dimension(KRpar),intent(in) ::  R
      real,dimension(KRpar), intent(in)  ::  RD      
	    real*4,                intent(in)  ::  pomin,pomax,  &
	                                           xnmin,xnmax,  & 
											                       xkmin,xkmax

      type(DLS_kernel_location),  intent(in)     ::  KERNEL
      type(DLS_kernel_grid1),     intent(in)     ::  XGRID1
      !type(DLS_kernel_grid1),     intent(inout)     ::  XGRID1
	    real,dimension(KCpar,KNpar),intent(inout)  ::  SD

!  ---------------------------------------------------------------------	  	  	  

      real,dimension(KMpar,KWLpar),intent(out) :: f11,f12,f22,f33,f34,f44      
      real,dimension(KWLpar),      intent(out) :: xext,xabs,xsca,albedo 
!  ---------------------------------------------------------------------	  	  	  
! local
      real                    ::  xnorm
      real*4                  ::  dlnr,dlnr1,coeff,coeff1
      real,dimension(KWLpar)  ::  g 
      real,dimension(2)       ::  X,Y 

      integer,          save  ::  NDP=0

!  ---------------------------------------------------------------------	  	  	  
      COMMON /US1/ US11(KMpar,KNpar,KWLpar) 
      COMMON /US2/ US12(KMpar,KNpar,KWLpar)
      COMMON /US3/ US22(KMpar,KNpar,KWLpar)
      COMMON /US4/ US33(KMpar,KNpar,KWLpar)
      COMMON /US5/ US34(KMpar,KNpar,KWLpar)
      COMMON /US6/ US44(KMpar,KNpar,KWLpar)
      COMMON /US0/ USEA(2,KNpar,KWLpar)
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------
!  ---------------------------------------------------------------------	  	  	  

!c
!c ** GET MATRICES
!c
!c      write(*,*) 'before USMATRIX_bin'

!      do itest=1,2 ! delete after test
!write(*,*) key,keyEL,keySUB,keyLS,key_RD1,key_RD,'  - key,keyEL,keySUB,keyLS,key_RD1,key_RD'
!write(*,*) KC,LB,LE,NWL,WL,RN,RK,KN,'  - KC,LB,LE,NWL,WL,RN,RK,KN'
!write(*,*) grid,'  - grid'
!write(*,*) KR,R,RD,'  - KR,R,RD'                       
!write(*,*) KM,ANGLE,'  - KM,ANGLE'
!write(*,*) dlnr,dlnr1,'  - dlnr,dlnr1'
!write(*,*) pomin,pomax,xnmin,xnmax,xkmin,xkmax,'  - pomin,pomax,xnmin,xnmax,xkmin,xkmax'
      
      CALL USMATRIX_bin(key,keyEL,keySUB,keyLS,key_RD1   &
                       ,key_RD                           &
                       ,NWL,WL,RN,RK,KN,grid             &
                       ,KR,R,RD,KM,ANGLE,dlnr,dlnr1      &
                       ,pomin,pomax                      &
                       ,KERNEL,XGRID1                    &
                       ,KC,SD,LB,LE                      &
                       ,xnmin,xnmax,xkmin,xkmax          &
                       ,KERNELS1,KERNELS2)

!         write(*,*) 'after USMATRIX_bin'

!         write(*,*) 'in OPTCHAR11: US11=',US11(1,1),
!     &    ' US11=',US11(KM,KN)
!c
!c ** CALCULATE APPROXIMATED OPTICAL CHARACTERISTICS      
!c
!cl      f11(:)=0.
!cl      f12(:)=0.
!cl      f22(:)=0.
!cl      f33(:)=0.
!cl      f34(:)=0.
!cl      f44(:)=0.
!cl      xext=0.
!cl      xabs=0.

!c	write(*,*) SD(1:KN)

      coeff1 = dlnr/dlnr1
!
! coeff=1e+3 because of at the time of kernel calculation cross section (um**2) was multiplied by coeffitient 1e-3 (???)
!            
! Kernels have been multiplied by coefficient 1e+3 (kernel units = 1/um), that is why coeff changed to 1.
!      coeff = 1e+3
      coeff = 1.
      if(key.eq.0) then
         SD(1,1:KN)=SD(1,1:KN)*coeff1 !*coeff
!c	else if(key.eq.2) then
      !else 
! SD should be used with same grid radii as one was used for Kernel_bin precalculation 
         !do i1=1,KC
         !SD(i1,1:KN)=SD(i1,1:KN) !*coeff*coeff1
         !enddo ! i1
      endif ! key

!c	write(*,*) 'in OPTCHAR - SD2:'
!c	write(*,*) SD(1:KN)
!      write(*,*) 'test DLS_optchr_bin.f US22(1,1:KN,1) :'
!	  write(*,'(5e14.5)') US22(1,1:KN,1)
!      write(*,*) 'test DLS_optchr_bin.f SD(1,1:KN) :'
!	  write(*,'(5e14.5)') SD(1,1:KN)	  
!	  write(*,*) 'LB=',LB,' LE=',LE
      	  	  	       	  
      IF(key.eq.0) THEN
      DO L=LB,LE  ! Wavelength loop
      xext(L) = coeff * DOT_PRODUCT(USEA(1,1:KN,L),SD(1,1:KN)) 
      xabs(L) = coeff * DOT_PRODUCT(USEA(2,1:KN,L),SD(1,1:KN))
      xsca(L) = xext(L)-xabs(L)
      albedo(L)=xsca(L)/xext(L)
!      write(*,*) 'dlnr=',dlnr
!      write(*,*) 'in optchar_bin SD: ',SD(1,:)
	  do j=1,KM
	  f11(j,L)=coeff*DOT_PRODUCT(US11(j,1:KN,L),SD(1,1:KN))                            
	  enddo ! j 
	    if(keySUB.eq.0) f11(1:KM,L)=f11(1:KM,L)/xsca(L)
	  if(keyEL.gt.1) then
	  do j=1,KM
	  f12(j,L)=coeff*DOT_PRODUCT(US12(j,1:KN,L),SD(1,1:KN))
	  enddo ! j
	    if(keySUB.eq.0) f12(1:KM,L)=-f12(1:KM,L)  &
                                    /xsca(L)/f11(1:KM,L)                             
	  endif   
      if(keyEL.gt.2) then
       do j=1,KM
       f22(j,L)=coeff*DOT_PRODUCT(US22(j,1:KN,L),SD(1,1:KN))                         
       enddo ! j
	    if(keySUB.eq.0) f22(1:KM,L)=f22(1:KM,L)   &
                                   /xsca(L)/f11(1:KM,L) 
      endif   
      if(keyEL.gt.3) then
       do j=1,KM
       f33(j,L)=coeff*DOT_PRODUCT(US33(j,1:KN,L),SD(1,1:KN))
       enddo ! j
	    if(keySUB.eq.0) f33(1:KM,L)=f33(1:KM,L)   &
                                   /xsca(L)/f11(1:KM,L) 
      endif   
      if(keyEL.gt.4) then
       do j=1,KM
       f34(j,L)=coeff*DOT_PRODUCT(US34(j,1:KN,L),SD(1,1:KN))
!c       write(*,*) 'US34:',US34(j,1,L)
       enddo ! j
	    if(keySUB.eq.0) f34(1:KM,L)=f34(1:KM,L)   &
                                   /xsca(L)/f11(1:KM,L) 
      endif 
      if(keyEL.gt.5) then
       do j=1,KM
       f44(j,L)=coeff*DOT_PRODUCT(US44(j,1:KN,L),SD(1,1:KN))
       enddo ! j
	    if(keySUB.eq.0) f44(1:KM,L)=f44(1:KM,L)   &
                                   /xsca(L)/f11(1:KM,L) 
      endif 

!      write(*,*)	  
!      write(*,*) 'test DLS_optchr_bin.f f22(1,1), f11(1,1), xsca(1) :'
!	  write(*,'(5e14.5)') f22(1,1),f11(1,1),xsca(1)
      ENDDO ! L

      ELSE ! key =1 or 2

      DO L=LB,LE  ! Wavelength loop
       xext(L) = coeff * DOT_PRODUCT(RC2(1:KC),USEA(1,1:KC,L))
       xabs(L) = coeff * DOT_PRODUCT(RC2(1:KC),USEA(2,1:KC,L))
       xsca(L) = xext(L)-xabs(L)
       albedo(L)=xsca(L)/xext(L)

       do j=1,KM
       f11(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US11(j,1:KC,L))
       enddo ! j
 	    if(keySUB.eq.0) f11(1:KM,L)=f11(1:KM,L)  &
                                   /xsca(L)        
      if(keyEL.gt.1) then
       do j=1,KM
       f12(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US12(j,1:KC,L))
       enddo ! j
	    if(keySUB.eq.0) f12(1:KM,L)=-f12(1:KM,L)  &
                                   /xsca(L)/f11(1:KM,L) 

      endif
      if(keyEL.gt.2) then
       do j=1,KM
       f22(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US22(j,1:KC,L))
       enddo ! j
	    if(keySUB.eq.0) f22(1:KM,L)=f22(1:KM,L)  &
                                   /xsca(L)/f11(1:KM,L) 

      endif
      if(keyEL.gt.3) then
       do j=1,KM
       f33(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US33(j,1:KC,L))
       enddo ! j
	    if(keySUB.eq.0) f33(1:KM,L)=f33(1:KM,L)  &
                                   /xsca(L)/f11(1:KM,L) 
      endif
      if(keyEL.gt.4) then
       do j=1,KM
       f34(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US34(j,1:KC,L))
       enddo ! j
	    if(keySUB.eq.0) f34(1:KM,L)=f34(1:KM,L)  &
                                   /xsca(L)/f11(1:KM,L) 

      endif 
      if(keyEL.gt.5) then
       do j=1,KM
       f44(j,L) = coeff * DOT_PRODUCT(RC2(1:KC),US44(j,1:KC,L))
       enddo ! j
	    if(keySUB.eq.0) f44(1:KM,L)=f44(1:KM,L)  &
                                   /xsca(L)/f11(1:KM,L) 	   
      endif 
      ENDDO ! L
	  
	ENDIF ! key
!c
!c ** SMOOTH f33 and f44 
!c
	do L=LB,LE
      if(keyEL.gt.3) then
       j=1
       do while(ANGLE(j).lt.(40.).or.j.gt.KM)
	 j1=j
	 j=j+1
	 enddo
       j=1
       do while(ANGLE(j).lt.(50.).or.j.gt.KM)
	 j2=j
	 j=j+1
	 enddo
       if(j2.gt.j1) then
	 X(1)=ANGLE(j1-1)
	 X(2)=ANGLE(j2+1)
	 Y(1)=f33(j1-1,L)
	 Y(2)=f33(j2+1,L)
	 Y(1)=LOG(Y(1))
	 Y(2)=LOG(Y(2))    
       do j=j1,j2
       f33(j,L)=EXP(LINEAR(X, Y, 2, ANGLE(j)))
	 enddo ! j
       endif ! j2&j1
      endif ! keyEL>3
!c
      if(keyEL.gt.5) then
      j=1
      do while(ANGLE(j).lt.(50.).or.j.gt.KM)
	j1=j
	j=j+1
	enddo
      j=1
      do while(ANGLE(j).lt.(60.).or.j.gt.KM)
	j2=j
	j=j+1
	enddo

      if(j2.gt.j1) then
	X(1)=ANGLE(j1-1)
	X(2)=ANGLE(j2+1)
	Y(1)=f44(j1-1,L)
	Y(2)=f44(j2+1,L)
	Y(1)=LOG(Y(1))
	Y(2)=LOG(Y(2))
      do j=j1,j2
      f44(j,L)=EXP(LINEAR(X, Y, 2, ANGLE(j)))
	  enddo ! j
      endif ! j2&j1
	  endif ! keyEL>5
      enddo ! L

!c ** Check f11 norma
!c
      if((KM .eq. 180 .or. KM .eq. 181) .and.   &
             (ANGLE(1) .eq. (0.) .and. ANGLE(KM) .eq. (180.))) then
       do L=LB,LE
       if(keySUB .eq. 0) then
          call SINT  ( ANGLE, f11(:,L), KM, xnorm )
!          call SINT  ( ANGLE, f11(:,L), KM, xnorm )
          call ASYPAR( ANGLE, f11(:,L), KM, g(L) )
!write(*,*) 'test test xnorm=',xnorm,'  xsca=',xsca(L)
!stop 'test stop in DLS_optchr_bin'
       else ! keySUB=1
          call SINT  ( ANGLE, f11(:,L)/xsca(L), KM, xnorm )
!         call SINT  ( ANGLE, f11(:,L)/xsca(L), KM, xnorm )
          call ASYPAR( ANGLE, f11(:,L)/xsca(L), KM, g(L) )
!write(*,*) 'test test xnorm=',xnorm,'  xsca=',xsca(L)
!stop 'test stop in DLS_optchr_bin'
        endif ! keySUB .eq. 0
	     enddo ! L
       endif ! KM 
      if(key.eq.0) then
       SD(1,1:KN)=SD(1,1:KN)/coeff1 !/coeff
	    !else 
       !do i1=1,KC
       !SD(i1,1:KN)=SD(i1,1:KN) !/coeff/coeff1
	     !enddo ! i1
	    endif ! key

!c ** WRITE APPROXIMATED OPTICAL CHARACHTERISTICS
!c
      IF(keySUB.eq.0) THEN      
        if(NDP.eq.0) then
         open(10,file='sca_mtrx.dat',status='unknown')
         NDP=1
	    else
         open(10,file='sca_mtrx.dat',status='unknown', &
         position='append')
	    endif

        write(10,*)
!c        write(10,*)'<<<     WELCOME TO NONSPHERICAL',
!c     &        ' AEROSOL WORLD     >>>'
!c        write(10,*)
        if(key.ne.0) then
          write(10,*)'Bin concentration:'
          do i=1,KC 
          write(10,'(i5,e17.4)') i,RC2(i)
          enddo ! i
          numKC=KC		
        else ! key=0
          numKC=1
        endif ! key.ne.0 

        write(10,*)'Size distribution:'
        write(10,*) ' r(mkm)             Sd(r):'
        do i=1,KN
        write(10,'(e11.4,22e12.4)') grid(i),(SD(i1,i),i1=1,numKC)
        enddo ! i
        if(key.eq.2) then
          write(10,*)'Size distribution (SUM of bins):'	  
          do i=1,KN
          write(10,'(e11.4,e12.4)') grid(i),  &
          !SUM(RC2(1:numKC)*SD(1:numKC,i))*(dlnr1/coeff)
          SUM(RC2(1:numKC)*SD(1:numKC,i)) !*dlnr1
          enddo      
        endif ! key=2
        write(10,*)
        write(10,*)'Axis ratio distribution:'
        write(10,*) '         R                Rd(R)'
        do i=1,KR
        write(10,'(2e17.5)') R(i),RD(i)
        enddo ! i

        do L=LB,LE    ! wavelength loop start
        write(10,*) 
        write(10,10) WL(L), RN(L), RK(L)
        write(10,*)
        if(key_RD.eq.1) then
         write(10,*) '        APPROXIMATED OPTICAL CHARACTERISTICS'
         write(10,*) '                 (volume mixture)'
        endif ! key_RD
        if(key_RD.eq.2) then 
         write(10,*) '        APPROXIMATED OPTICAL CHARACTERISTICS'
         write(10,*) '               (surface area mixture)'
        endif ! key_RD
		write(10,*) 
        write(10,11) xext(L), xabs(L), xsca(L), albedo(L)
 
        if(keyEL.eq.1) then
		 write(10,'(2x,''ANGLE'',7x,''F11'')') 
		 do j=1,KM       
		 write(10,12) ANGLE(j),f11(j,L)
		 enddo ! j
        endif ! keyEL=1

        if(keyEL.eq.2) then
         write(10,'(2x,''ANGLE'',7x,''F11'',7x,''-F12/F11'')') 
         do j=1,KM       
         write(10,12) ANGLE(j),f11(j,L),f12(j,L)
         enddo ! j
        endif ! keyEL=2

        if(keyEL.eq.3) then
         write(10,'(2x,''ANGLE'',7x,''F11'',7x,''-F12/F11'',7x,''F22/F11'')') 
         do j=1,KM       
         write(10,12) ANGLE(j),f11(j,L),f12(j,L),f22(j,L)
         enddo ! j
        endif !keyEL=3

        if(keyEL.eq.4) then
         write(10,'(2x,''ANGLE'',7x,''F11'',7x,''-F12/F11'',7x,''F22/F11'',7x,''F33/F11'')')
         do j=1,KM       
         write(10,12) ANGLE(j),f11(j,L),f12(j,L),f22(j,L)  &
                        ,f33(j,L)
         enddo ! j
        endif ! keyEL=4

        if(keyEL.eq.5) then
         write(10,'(2x,''ANGLE'',7x,''F11'',7x,''-F12/F11'',7x,''F22/F11'',7x,''F33/F11'',7x,''F34/F11'')') 
         do j=1,KM       
         write(10,12) ANGLE(j),f11(j,L),f12(j,L),f22(j,L)  &
                        ,f33(j,L),f34(j,L)
         enddo ! j
        endif ! keyEL=5

        if(keyEL.eq.6) then
         write(10,'(2x,''ANGLE'',7x,''F11'',7x,''-F12/F11'',7x,''F22/F11'',7x,''F33/F11'',7x,''F34/F11'',7x,''F44/F11'')') 
         do j=1,KM       
         write(10,12) ANGLE(j),f11(j,L),f12(j,L),f22(j,L)  &
                        ,f33(j,L),f34(j,L),f44(j,L)
         enddo ! j
        endif ! keyEL=6
!c
        if((KM.eq.180.or.KM.eq.181).and.  &
             (ANGLE(1).eq.(0.).and.ANGLE(KM).eq.(180.)))  &
         write(10,'(''asymmetry parameter='',e13.5)') g(L)

      ENDDO ! L       ! wavelength loop end

      close(10)

      ENDIF ! key_SUB=0
	  
!13    format('check: f11 norma=',f7.3,'   one=',f7.3) 
10    FORMAT('wavelength=',e11.4,'  n=',e11.4,'  k=',e11.4) 
11    FORMAT('ext=',e13.5,'  abs=',e13.5,  &
             '  sca=',e13.5,'  albedo=',e13.5) 
12    FORMAT(F7.2,6E14.5)

!      write(*,*) 'in OPTCHAR : '
!      write(*,*) 'NWL=',NWL,' WL=',WL(1),' RN=',RN(1),' RK=',RK(1)
!      write(*,*) 'f11=',f11(1,1),' f22=',f22(1,1),' f33=',f33(1,1)
	  
!      enddo ! itest  - delete after test

      RETURN 
      END SUBROUTINE OPTCHAR_bin

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
