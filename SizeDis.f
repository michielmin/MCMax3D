c-----------------------------------------------------------------------
c This subroutine calculates the grainsize distribution according to a
c MRN-like powerlaw.
c It reads the grainsize from rgrain (Because Grain()%rv is not allocated),
c and sets the abundance in abun() (which is warg() in Init.f) 
c
c Input keywords:
c (mrn_rmin)    lower bound of the grain size grid.
c (mrn_rmax)    upper bound of the grain size grid.
c (mrn_index)   slope of the size distribution (def: 3.5)
c (mrn_ngrains) number of grains to include in calculation
c
c  Abundances of particles outside the sizegrid (mrn_rmin,mrn_rmax) are 
c  set to zero.
c
c  Helper routines: gsd_edges
c
c-----------------------------------------------------------------------


	subroutine SetSizeDis(abun,ipart,izone)
	use GlobalSetup
	IMPLICIT NONE
	real*8 abun(1:maxns),rgrain(1:maxns),r_int(1:maxns+1)
	real*8 amin,amax,apow
	integer ipart,izone,i,nsize
	integer            :: ii

	nsize=Part(ipart)%nsize
	rgrain(1:nsize)=Part(ipart)%rv(1:nsize)

	amin=Zone(izone)%amin*1d-4
	amax=Zone(izone)%amax*1d-4
	apow=Zone(izone)%apow

      !  Cell edges
	call gsd_edges(rgrain(1:nsize),r_int(1:nsize+1),nsize)
      
      !  Calc distribution
	call gsd_calc_mrn(abun(1:nsize),rgrain(1:nsize),r_int(1:nsize+1),nsize,amin,amax,apow)

      !  Made it!
      !
	return
	end

c-----------------------------------------------------------------------
c This subroutine calculates the cell edges for the grain size bins.
c It reads the grainsize from rgrain (Because Grain()%rv is not allocated),
c and sets the minimum and maximum size in rgrain_min and rgrain_max
c
c Input keywords:
c (rgrain)      grain sizes
c (mrn_ngrains) number of grains to include in calculation
c
c  The grid edges of particles outside the sizegrid (mrn_rmin,mrn_rmax) are 
c  not set.
c
c-----------------------------------------------------------------------

	subroutine gsd_edges(rgrain,r_int,ngrains)
	implicit none
	integer ngrains
	real*8 rgrain(1:ngrains)
	real*8 logr_int(1:ngrains+1),r_int(1:ngrains+1)
	real*8 rmin,rmax
	integer i,j,ii

      !  Find cell interfaces for grain size grid
      !
	do ii=2,ngrains  ! inbetween 
		logr_int(ii)=( log(rgrain(ii-1))+log(rgrain(ii)) )/2d0
	enddo
	logr_int(1)=         2*log(rgrain(1))       -logr_int(2)    ! lower bound
	logr_int(ngrains+1)= 2*log(rgrain(ngrains)) -logr_int(ngrains) ! upper bound

      ! the final grid
	r_int(1:ngrains+1)= exp(logr_int(1:ngrains+1))

      !  Check if the grid is ok
	do ii=1,ngrains
		if (r_int(ii).ge.rgrain(ii).or.r_int(ii+1).le.rgrain(ii)) then 
    		write(*,'("grain",i02," doesnt fit in the grid")') ii
			write(*,'("check that rgrain",i02," satisfies ",
     1                f10.3," < ",f10.3," < ",f10.3)')
     2           ii,r_int(ii),rgrain(ii),r_int(ii+1)
			write(*,*) "grid doesn't match grain size -> stop 67876"
			stop 67876
		endif
	enddo

	return
	end

c-----------------------------------------------------------------------
c This subroutine calculates the mass per grain binsize, i.e. abundance
c-----------------------------------------------------------------------

	subroutine gsd_calc_mrn(abun,rgrain,r_int,ngrains,mrn_rmin,mrn_rmax,mrn_index)
	implicit none
	integer ngrains
	real*8 mrn_rmin,mrn_rmax,mrn_index
	real*8 abun(1:ngrains),rgrain(1:ngrains),r_int(1:ngrains+1)
	real*8 rmin,rmax,r_index,r_factor,abun_norm
	integer i,j,ii

      !
      !  NOTE: MRN distribution is: f(a)=a^-3.5  , general a^-n
      !        in units of mass:    f(m)=m^-11/6 , general m^-(n+2)/3
      !
      !        For the abundance (or mass per size bin), integrate m f(m) dm
      !
      !        abun = 3/(4-n) m^( (4-n)/3 )  , between m(r_in) and m(r_out)
      !          OR   (check this)
      !        abun = 3/(4-n) r^(4-n)
      !
      !        note that it changes sign at n=4, where abun=log(r)  (ln, not 10log)
      !
	r_index=4d0-mrn_index
	r_factor=3d0/r_index  ! needed for sign change
c      if (diagnose) write(*,'("n= ",f4.1," , rindex= ",f4.1)') mrn_index,r_index
      !
	do ii=1,ngrains
		rmin=r_int(ii)
		rmax=r_int(ii+1)

         !  check if grainsize is between (mrn_rmin,mrn_rmax)
		if (rmin.ge.mrn_rmax .or. rmax.le.mrn_rmin) then
			abun(ii)=0d0  ! min or max outside grid cell
		else
			if (rmin.le.mrn_rmin) rmin=mrn_rmin ! min in gridcell
			if (rmax.ge.mrn_rmax) rmax=mrn_rmax ! max in gridcell
			if (rmax.le.rmin) stop 88354  ! shouldn't be
                                
            !  calculate abundance
			if (abs(mrn_index-4d0).le.1d-6) then 
				abun(ii)=log(rmax)-log(rmin)
			else
				abun(ii)=r_factor*(rmax**r_index - rmin**r_index) ! not normalized
			endif
		endif

	enddo
      
      !  Calculate total mass & normalize
      !
	rmin=r_int(1)
	rmax=r_int(ngrains+1)
	if (rmin.ge.mrn_rmax .or. rmax.le.mrn_rmin) stop 16661 ! no grains between mrn_rmin,mrn_rmax
	if (rmin.le.mrn_rmin) rmin=mrn_rmin ! min in gridcell
	if (rmax.ge.mrn_rmax) rmax=mrn_rmax ! max in gridcell
                                
      !  calculate abundance
	if (abs(mrn_index-4d0).le.1d-6) then 
		abun_norm=log(rmax)-log(rmin)
	else
		abun_norm=r_factor*(rmax**r_index - rmin**r_index) ! not normalized
	endif
      
      !  Normalize abundances
	abun(1:ngrains)=abun(1:ngrains) / abun_norm
      
      !  MCMax crashes at low, non-zero abundances
	do ii=1,ngrains
		if (abun(ii) .le.1d-100) abun(ii)=0d0
	enddo

      !  Check normalization
      !
	if (abs(sum(abun(1:ngrains))-1d0).ge.1d-6) then
		write(*,'("Something wrong with abundance normalization")')
		do ii=1,ngrains
			write(*,'(" abun",i02,"= ",f10.8)') ii,abun(ii)
		enddo
		write(*,'("Sum=",f10.8," =!= 1 -> stop 63636")') sum(abun(1:ngrains))
		stop 63636
	endif


	end subroutine 

