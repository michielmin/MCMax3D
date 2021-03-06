	module RandomWalkModule
	IMPLICIT NONE
	integer NY
	parameter(NY=1000)
	real*8 phi(NY),y(NY)
	end module RandomWalkModule
	
	
	real*8 function MinDist(phot)
	use Parameters
	IMPLICIT NONE
	type(photon) phot
	real*8 r,R1,R2,T1,T2,rho,vR1,vR2,vT1,vT2,ST1,ST2,z
	integer izone

	MinDist=1d200
	do izone=1,nzones
		if(phot%inzone(izone)) then
			rho=phot%xzone(izone)**2+phot%yzone(izone)**2
			r=sqrt(rho+phot%zzone(izone)**2)
			rho=sqrt(rho)
			z=abs(phot%zzone(izone))
			R1=Zone(izone)%R(phot%i1(izone))
			R2=Zone(izone)%R(phot%i1(izone)+1)
			T1=cos(Zone(izone)%theta(phot%i2(izone)))
			T2=cos(Zone(izone)%theta(phot%i2(izone)+1))
			ST1=cos(Zone(izone)%theta(phot%i2(izone)))
			ST2=cos(Zone(izone)%theta(phot%i2(izone)+1))
	
			vR1=abs(r-R1)
			vR2=abs(R2-r)
			vT1=abs(T1*rho-ST1*z)
			vT2=abs(T2*rho-ST2*z)
	
			if(vR1.lt.MinDist) MinDist=vR1
			if(vR2.lt.MinDist) MinDist=vR2
			if(vT1.lt.MinDist) MinDist=vT1
			if(vT2.lt.MinDist) MinDist=vT2
		endif
	enddo

	return
	end
	
	
	
	logical function RandomWalk(phot,dmin,nRWinteract)
	use Parameters
	use RandomWalkModule
	IMPLICIT NONE
	real*8 dmin,v,ran2,x,lr,increaseT,increaseTP,Kabs,Kext
	real*8 EJv,T,T0,spec(nlam),nRWinteract
	type(photon) phot
	integer i,iT,iy,iT0,iter,iT1,l,ii,iopac

	RandomWalk=.false.
	
	iT=int(C(phot%i,phot%j)%T/dT)
	if(iT.eq.0) iT=1
	if(iT.gt.TMAX) iT=TMAX
	if(iT.ne.C(phot%i,phot%j)%iTD) call DiffCoeffCell(phot%i,phot%j,iT)	

	Kext=C(phot%i,phot%j)%KDext
	
	lr=1d0/(AU*Kext*C(phot%i,phot%j)%dens)

	if(dmin.le.factRW*lr) return

	C(phot%i,phot%j)%randomwalk=.true.

	RandomWalk=.true.

	x=ran2(idum)
	iy=1
	call hunt(phi,NY,x,iy)

	iT0=iT
	T0=C(phot%i,phot%j)%T

	iter=0
1	continue

	v=-3d0*log(y(iy))*dmin**2/(lr*pi**2)

	totaldist=totaldist+v

	Kabs=C(phot%i,phot%j)%KDabs

	call randomdirection(phot%vx,phot%vy,phot%vz)
	phot%x=phot%x+dmin*phot%vx
	phot%y=phot%y+dmin*phot%vy
	phot%z=phot%z+dmin*phot%vz

	EJv=phot%E*v*AU*C(phot%i,phot%j)%dens

	C(phot%i,phot%j)%Eabs=C(phot%i,phot%j)%Eabs+EJv*Kabs!*C(phot%i,phot%j)%dens*C(phot%i,phot%j)%V
	if(.not.tcontact.or.tdes_iter) then
		do i=1,ngrains
			C(phot%i,phot%j)%EabsP(i)=C(phot%i,phot%j)%EabsP(i)+EJv*Kabs!*C(phot%i,phot%j)%dens*C(phot%i,phot%j)%V
		enddo
	endif

c==============================================================================
c Change 15-02-2012, do not use the random walk module to compute the EJv's
c
	C(phot%i,phot%j)%EJv=C(phot%i,phot%j)%EJv+EJv*Kabs/C(phot%i,phot%j)%dens
	nEJv=nEJv+EJv*Kabs/C(phot%i,phot%j)%dens
c	if(use_qhp) then
c		do ii=1,ngrains
c			if(Grain(ii)%qhp) then
c				C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)=C(phot%i,phot%j)%EJvQHP(Grain(ii)%qhpnr)
c     &							+EJv*C(phot%i,phot%j)%KDQHP*C(phot%i,phot%j)%w(ii)
c			endif
c		enddo
c	endif
	if(.not.tcontact.or.tdes_iter) then
		do i=1,ngrains
			C(phot%i,phot%j)%EJvP(i)=C(phot%i,phot%j)%EJvP(i)+EJv*Kabs/C(phot%i,phot%j)%dens
		enddo
	endif
c==============================================================================

	C(phot%i,phot%j)%T=increaseT(phot)

	if(.not.tcontact) then
		do i=1,ngrains
			C(phot%i,phot%j)%TP(i)=increaseTP(phot,i)
		enddo
	endif

	nRWinteract=v/dmin
	tautot=tautot+v/dmin

	if(computeLRF) then
		call addRW_LRF(phot%i,phot%j,iT,EJv)
	endif

	return
	end
	



	subroutine InitRandomWalk()
	use RandomWalkModule
	use Parameters
	IMPLICIT NONE
	integer i,j,n,nmax,iT,ii,iT1,iT2,diT
	real*8 lr,Kext
	nmax=1000
	
	write(*,'("Initializing Diffusion/Random Walk")')
	write(9,'("Initializing Diffusion/Random Walk")')
	y(1)=1d0
	phi(1)=1d0
	do i=2,NY
		y(i)=real(NY-i+1)/real(NY-1)*0.75d0
		phi(i)=0d0
		do n=1,nmax
			phi(i)=phi(i)+(-1d0)**(n+1)*y(i)**(n**2)
		enddo
		phi(i)=phi(i)*2d0
	enddo
	
	iT1=100d0/dT
	iT2=1600d0/dT
	diT=(iT2-iT1)/3
	do i=1,D%nR-1
	call tellertje(i,D%nR-1)
	do j=1,D%nTheta-1
		C(i,j)%thick=.false.
		do iT=iT1,iT2,diT
			call DiffCoeffCell(i,j,iT)
			C(i,j)%randomwalk=.false.
			Kext=C(i,j)%KDext
			lr=1d0/(AU*Kext*C(i,j)%dens)
			if((D%R(i+1)-D%R(i)).gt.factRW*lr) then
				C(i,j)%thick=.true.
			endif
		enddo
	enddo
	enddo


	return
	end
	



