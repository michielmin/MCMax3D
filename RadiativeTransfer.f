	subroutine RadiativeTransfer
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	character*500 fitsfile

	allocate(phot%i1(nzones))
	allocate(phot%i2(nzones))
	allocate(phot%i3(nzones))
	allocate(phot%inzone(nzones))
	allocate(phot%edgeNr(nzones))

	call InitRadiativeTransfer
	
	call output("==================================================================")
	call output("Emitting " // trim(int2string(Nphot,'(i10)')) // " photon packages")
	do i=1,Nphot
		call tellertje(i,Nphot)
		call EmitPhoton(phot)
		call InWhichZones(phot)
		call TravelPhoton(phot)
		call MCoutput(phot)
	enddo

	call output("==================================================================")
	call output("Writing Monte Carlo observables")
	do i=1,nMCobs
		fitsfile="MCout" // trim(int2string(i,'(i0.4)')) // ".fits"
		call writefitsfile(fitsfile,MCobs(i)%image,nlam,MCobs(i)%npix)
	enddo
c	call DetermineTemperatures
	
	return
	end

	subroutine InWhichZones(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 r
	integer i
	
	do i=1,nzones
		r=sqrt((phot%x-Zone(i)%x0)**2+(phot%y-Zone(i)%y0)**2+(phot%z-Zone(i)%z0)**2)
		phot%inzone(i)=.false.
		if(r.lt.Zone(i)%Rout.and.r.gt.Zone(i)%Rin) then
			phot%inzone(i)=.true.
			call determine_i1(phot,i)
			call determine_i2(phot,i)
			call determine_i3(phot,i)
		endif
	enddo
	return
	end
	

	subroutine TravelPhoton(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,imin
	logical leave
	real*8 minv,Kext,tau0,tau,GetKext,random
	type(Travel) Trac(nzones)
	type(Cell),pointer :: C
	type(Photon) phot
	integer iter

	tau0=-log(random(idum))
	iter=0
1	continue
	iter=iter+1
	
	call checkphot(phot)
	
	Kext=0d0
	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone))
			end select
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			Kext=Kext+GetKext(phot%ilam1,C)*phot%wl1+GetKext(phot%ilam2,C)*phot%wl2
		else
			select case(Zone(izone)%shape)
				case("SPH")
					call HitSph(phot,izone,Trac(izone))
			end select
		endif
	enddo
	minv=20d0*maxR
	leave=.true.
	do izone=1,nzones
		if(Trac(izone)%v.gt.0d0.and.Trac(izone)%v.lt.minv) then
			minv=Trac(izone)%v
			imin=izone
			leave=.false.
		endif
	enddo

	if(leave) goto 2

	if((tau0-Kext*minv).lt.0d0) then
		minv=tau0/Kext
		phot%edgeNr=0
		phot%x=phot%x+minv*phot%vx
		phot%y=phot%y+minv*phot%vy
		phot%z=phot%z+minv*phot%vz
		call Interact(phot)
		tau0=-log(random(idum))
		goto 1
	endif

	phot%x=phot%x+minv*phot%vx
	phot%y=phot%y+minv*phot%vy
	phot%z=phot%z+minv*phot%vz

	do izone=1,nzones
		if(Trac(izone)%v.le.minv.or.izone.eq.imin) then
			phot%i1(izone)=Trac(izone)%i1next
			phot%i2(izone)=Trac(izone)%i2next
			phot%i3(izone)=Trac(izone)%i3next

			if(phot%i1(izone).eq.-1) call determine_i1(phot,izone)
			if(phot%i2(izone).eq.-1) call determine_i2(phot,izone)
			if(phot%i3(izone).eq.-1) call determine_i3(phot,izone)

			phot%edgeNr(izone)=Trac(izone)%edgenext
			phot%inzone(izone)=.true.
			if(phot%i1(izone).gt.Zone(izone)%n1.or.
     &			phot%i2(izone).gt.Zone(izone)%n2.or.
     &			phot%i3(izone).gt.Zone(izone)%n3.or.
     &			phot%i1(izone).lt.1.or.
     &			phot%i2(izone).lt.1.or.
     &			phot%i3(izone).lt.1) then
				phot%inzone(izone)=.false.
			endif
		else
			phot%edgeNr(izone)=0
		endif
	enddo
	tau0=tau0-Kext*minv

	goto 1
2	continue
	
	return
	end
	

	subroutine checkphot(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	integer i
	real*8 r
	logical oops

	oops=.false.	
	do i=1,nzones
		r=sqrt((phot%x-Zone(i)%x0)**2+(phot%y-Zone(i)%y0)**2+(phot%z-Zone(i)%z0)**2)
		if(phot%inzone(i)) then
			if(r.lt.Zone(i)%R(phot%i1(i))*0.99.or.
     &		   r.gt.Zone(i)%R(phot%i1(i)+1)*1.01) then
     			print*,'zone',i,r/AU
     			print*,'oeps',phot%x/AU,phot%y/AU,phot%z/AU
				print*,Zone(i)%R(phot%i1(i))/AU
				print*,Zone(i)%R(phot%i1(i)+1)/AU
				print*,phot%i1
				print*,phot%inzone
     			oops=.true.
     		endif
		else
			if(r.gt.Zone(i)%R(1)*1.01.and.
     &		   r.lt.Zone(i)%R(Zone(i)%nr+1)*0.99) then
     			print*,'zone',i,r/AU
c				print*,'oeps',phot%x/AU,phot%y/AU,phot%z/AU
c				print*,phot%i1
				print*,phot%inzone
c				oops=.true.
     		endif
		endif
	enddo
c	if(oops) then
c		stop
c	endif
	return
	end
	

	subroutine Interact(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot

	call randomdirection(phot%vx,phot%vy,phot%vz)

	phot%x0=phot%x
	phot%y0=phot%y
	phot%z0=phot%z

	return
	end
	

	subroutine MCoutput(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 inp,tot,x,y,z
	integer i,j,ix,iy
	type(Photon) phot
	
	do i=1,nMCobs
		inp=MCobs(i)%x*phot%vx+MCobs(i)%y*phot%vy+MCobs(i)%z*phot%vz
		if(inp.gt.MCobs(i)%opening) then
			do j=1,nlam
				if(column(j).lt.1000d0) then
					specemit(j)=specemit(j)*exp(-column(j))
				else
					specemit(j)=0d0
				endif
			enddo
			call integrate(specemit,tot)
			specemit=specemit/tot
			specemit=1d23*phot%sI*specemit
			MCobs(i)%spec(1:nlam)=MCobs(i)%spec(1:nlam)+specemit(1:nlam)
			x=phot%x0
			y=phot%y0
			z=phot%z0
			call rotateZ(x,y,z,cos(MCobs(i)%phi),sin(MCobs(i)%phi))
			call rotateY(x,y,z,cos(MCobs(i)%theta),sin(MCobs(i)%theta))
			ix=real(MCobs(i)%npix)*(x+maxR)/(2d0*maxR)
			iy=real(MCobs(i)%npix)*(y+maxR)/(2d0*maxR)
			if(ix.lt.MCobs(i)%npix.and.iy.lt.MCobs(i)%npix.and.ix.gt.0.and.iy.gt.0) then
				MCobs(i)%image(ix,iy,1:nlam)=MCobs(i)%image(ix,iy,1:nlam)+specemit(1:nlam)
			endif
		endif
	enddo
	
	return
	end
	

	subroutine EmitPhoton(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	real*8 Ltot,r,random
	
	Ltot=0d0
	do i=1,nstars
		Ltot=Ltot+Star(i)%L
	enddo
	r=Ltot*random(idum)
	do i=1,nstars
		r=r-Star(i)%L
		if(r.lt.0d0) exit
	enddo

	call randomdirection(phot%x,phot%y,phot%z)

	phot%x=Star(i)%R*phot%x
	phot%y=Star(i)%R*phot%y
	phot%z=Star(i)%R*phot%z

	call emit(phot,Star(i)%F,Star(i)%L)

	phot%sI=Ltot/real(Nphot)
	phot%sQ=0d0
	phot%sU=0d0
	phot%sV=0d0

	call randomdirection(phot%vx,phot%vy,phot%vz)
	if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
		phot%vx=-phot%vx
		phot%vy=-phot%vy
		phot%vz=-phot%vz
	endif
	
	phot%x=phot%x+Star(i)%x
	phot%y=phot%y+Star(i)%y
	phot%z=phot%z+Star(i)%z
	
	phot%edgeNr=0
	phot%inzone=.false.
	
	return
	end
	


	subroutine InitRadiativetransfer
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,ilam,i
	type(Cell),pointer :: C
	real*8 GetKext

	call output("Initializing radiative transfer")
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
			do i2=1,Zone(izone)%n2
				do i3=1,Zone(izone)%n3
					C => Zone(izone)%C(i1,i2,i3)
					C%E=0d0
					C%Ni=0d0
					C%T=0d0
				enddo
			enddo
		enddo
	enddo
	
	do i=1,nMCobs
		MCobs(i)%spec=0d0
		MCobs(i)%image=0d0
	enddo
	
	return
	end


c This function gives the Planck mean opacity of the cell
	real*8 function GetKp(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C
	
	GetKp=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKp=GetKp+C%densP(ipart,isize,iT)*Part(ipart)%Kp(isize,iT,ilam)
			enddo
		enddo
	enddo
c multiply by the volume to get the total mass
	GetKp=GetKp*C%V
	
	return
	end
	

c This function gives the absorption cross section per cm
	real*8 function GetKabs(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C
	
	GetKabs=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKabs=GetKabs+C%densP(ipart,isize,iT)*Part(ipart)%Kabs(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	

c This function gives the extinction cross section per cm
	real*8 function GetKext(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C

	GetKext=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKext=GetKext+C%densP(ipart,isize,iT)*Part(ipart)%Kext(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	
	
c This function gives the scattering cross section per cm
	real*8 function GetKsca(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,ilam,ipart,iT,isize
	type(Cell) C

	GetKsca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKsca=GetKsca+C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	
	
	
c This subroutine gets the scattering matrix in a certain cell
	subroutine function SetFmatrix(ilam,C,F,Ksca)
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,ilam,ipart,iT,isize
	real*8 Ksca
	type(Mueller) F
	type(Cell) C
	
	Ksca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				Ksca=Ksca+C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
			enddo
		enddo
	enddo
	F%F11(1:180)=0d0
	F%F12(1:180)=0d0
	F%F22(1:180)=0d0
	F%F33(1:180)=0d0
	F%F34(1:180)=0d0
	F%F44(1:180)=0d0
	F%IF11=0d0
	F%IF12=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				F%F11(1:180)=F%F11(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F11(1:180)
				F%F12(1:180)=F%F12(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F12(1:180)
				F%F22(1:180)=F%F22(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F22(1:180)
				F%F33(1:180)=F%F33(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F33(1:180)
				F%F34(1:180)=F%F34(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F34(1:180)
				F%F44(1:180)=F%F44(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F44(1:180)
				F%IF11=F%IF11+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF11
				F%IF12=F%IF12+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF12
			enddo
		enddo
	enddo
	F%F11=F%F11/Ksca
	F%F12=F%F12/Ksca
	F%F22=F%F22/Ksca
	F%F33=F%F33/Ksca
	F%F34=F%F34/Ksca
	F%F44=F%F44/Ksca
	F%IF11=F%IF11/Ksca
	F%IF12=F%IF12/Ksca

	Ksca=Ksca/C%dens
	
	return
	end


	subroutine determine_i1(phot,izone)
	use GlobalSetup
	IMPLICIT NONE
	integer i,izone
	real*8 r2
	type(Photon) phot
	
	select case(Zone(izone)%shape)
		case("SPH")
			r2=(phot%x-Zone(izone)%x0)**2+(phot%y-Zone(izone)%y0)**2+(phot%z-Zone(izone)%z0)**2
		case("CYL")
			r2=(phot%x-Zone(izone)%x0)**2+(phot%y-Zone(izone)%y0)**2
	end select

	do i=1,Zone(izone)%nr
		if(r2.ge.Zone(izone)%R2(i).and.r2.le.Zone(izone)%R2(i+1)) exit
	enddo
	phot%i1(izone)=i
	
	return
	end
	


	subroutine determine_i2(phot,izone)
	use GlobalSetup
	IMPLICIT NONE
	integer i,izone
	real*8 theta,r
	type(Photon) phot
	
	r=sqrt((phot%x-Zone(izone)%x0)**2+(phot%y-Zone(izone)%y0)**2+(phot%z-Zone(izone)%z0)**2)
	theta=acos((phot%z-Zone(izone)%z0)/r)
	
	do i=1,Zone(izone)%nt
		if(theta.ge.Zone(izone)%theta(i).and.theta.le.Zone(izone)%theta(i+1)) exit
	enddo
	phot%i2(izone)=i
	
	return
	end
	
	subroutine determine_i3(phot,izone)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,izone
	real*8 phi,r
	type(Photon) phot
	
	r=sqrt((phot%x-Zone(izone)%x0)**2+(phot%y-Zone(izone)%y0)**2)
	phi=acos((phot%x-Zone(izone)%x0)/r)
	if((phot%y-Zone(izone)%y0).lt.0d0) phi=2d0*pi-phi
	
	do i=1,Zone(izone)%nr
		if(phi.ge.Zone(izone)%phi(i).and.phi.le.Zone(izone)%phi(i+1)) exit
	enddo
	phot%i3(izone)=i
	
	return
	end
	
	
	
	