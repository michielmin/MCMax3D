	subroutine integrate(spec,L)
	use GlobalSetup
	IMPLICIT NONE
	real*8 spec(nlam),L
	integer i
	real*8 nu1,nu2,Iv1,Iv2
	
	L=0d0
	do i=1,nlam-1
		nu1 = nu(i)
		nu2 = nu(i+1)
		Iv1 = spec(i)
		Iv2 = spec(i+1)
		L = L + ABS(nu1-nu2)*0.5*(Iv1+Iv2)
	enddo

	return
	end



	subroutine emit(phot,spec,Lttot)
	use GlobalSetup
	IMPLICIT NONE
	real*8 spec(nlam),Lr,Lt,random,Ltold,Lttot,x,y,z,r
	integer i,ii,iopac
	type(Photon) phot
	
	specemit=spec/Lttot
	column=0d0

	call randomdirection(phot%vx,phot%vy,phot%vz)
	phot%scatt=.false.
	phot%pol=.false.

	x=-phot%vy
	y=phot%vx
	z=0d0
	r=sqrt(x**2+y**2+z**2)
	phot%Sx=x/r
	phot%Sy=y/r
	phot%Sz=z/r

	Lr=random(idum)*Lttot

	Ltold=0d0
	Lt=0d0
	do i=1,nlam-1
		Lt=Lt+ABS(nu(i)-nu(i+1))*0.5*(spec(i)+spec(i+1))
		if(Lt.ge.Lr) then
			phot%wl1=(Lr-Ltold)/(Lt-Ltold)
			phot%wl2=1d0-phot%wl1
			phot%nu=nu(i)*phot%wl1+nu(i+1)*phot%wl2
			phot%lam=2.9979d14/phot%nu
			phot%ilam1=i
			phot%ilam2=i+1
			goto 1
		endif
		Ltold=Lt
	enddo
	phot%wl1=0d0
	phot%wl2=1d0
	phot%nu=nu(nlam)
	phot%lam=2.9979d14/phot%nu
	phot%ilam1=nlam-1
	phot%ilam2=nlam

	phot%x0=phot%x
	phot%y0=phot%y
	phot%z0=phot%z

1	continue

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine randomdirection(x,y,z)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x,y,z,r,random
	
1	continue
	x=2d0*random(idum)-1d0
	y=2d0*random(idum)-1d0
	z=2d0*random(idum)-1d0
	r=x**2+y**2+z**2
	if(r.gt.1d0) goto 1
	r=sqrt(r)
	x=x/r
	y=y/r
	z=z/r
	
	return
	end

