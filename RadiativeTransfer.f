	subroutine RadiativeTransfer
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon),allocatable :: phot(:)
	integer i,j,maxnopenmp,omp_get_max_threads,iopenmp,omp_get_thread_num
	real*8 starttime,stoptime

	maxnopenmp=omp_get_max_threads()+1
	allocate(phot(maxnopenmp))	
	do i=1,maxnopenmp
		allocate(phot(i)%i1(nzones))
		allocate(phot(i)%i2(nzones))
		allocate(phot(i)%i3(nzones))
		allocate(phot(i)%inzone(nzones))
		allocate(phot(i)%edgeNr(nzones))
		allocate(phot(i)%KabsZ(nzones))
		allocate(phot(i)%xzone(nzones))
		allocate(phot(i)%yzone(nzones))
		allocate(phot(i)%zzone(nzones))
		allocate(phot(i)%vxzone(nzones))
		allocate(phot(i)%vyzone(nzones))
		allocate(phot(i)%vzzone(nzones))
	enddo

	call InitRadiativeTransfer
	
	call output("==================================================================")
	call output("Emitting " // trim(int2string(Nphot,'(i10)')) // " photon packages")

	call cpu_time(starttime)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i,iopenmp)
!$OMP& SHARED(Nphot,nzones,phot)
!$OMP DO SCHEDULE(STATIC,1)
	do i=1,Nphot
		call tellertje(i,Nphot)
		iopenmp=omp_get_thread_num()+1
		phot(iopenmp)%nr=i
		call EmitPhoton(phot(iopenmp))
		call InWhichZones(phot(iopenmp))
		call TravelPhoton(phot(iopenmp))
		call MCoutput(phot(iopenmp))
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call cpu_time(stoptime)
	call output("Radiative transfer time: " // trim(dbl2string(stoptime-starttime,'(f10.3)')) // " s")

	call DetermineTemperatures

	return
	end

	subroutine InWhichZones(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 r
	integer i
	
	do i=1,nzones
		r=sqrt(phot%xzone(i)**2+phot%yzone(i)**2+phot%zzone(i)**2)
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
	integer izone,imin,istar,status
	logical leave,hitstar0
	real*8 minv,tau0,tau,GetKext,random,GetKabs
	type(Travel) Trac(nzones)
	type(Travel) TracStar(nstars)
	type(Cell),pointer :: C
	type(Photon) phot

	tau0=-log(random(idum))
	do izone=1,nzones
		Trac(izone)%recompute=.true.
	enddo
	do istar=1,nstars
		TracStar(istar)%recompute=.true.
	enddo
	minv=0d0
	
1	continue

	phot%Kext=0d0
	phot%Kabs=0d0
	do izone=1,nzones
		phot%KabsZ(izone)=0d0
		if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			phot%Kext=phot%Kext+GetKext(phot%ilam1,C)*phot%wl1+GetKext(phot%ilam2,C)*phot%wl2
			phot%KabsZ(izone)=GetKabs(phot%ilam1,C)*phot%wl1+GetKabs(phot%ilam2,C)*phot%wl2
			phot%Kabs=phot%Kabs+phot%KabsZ(izone)
		endif
	enddo

	status=0
	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					if(Trac(izone)%recompute) then
						call TravelSph(phot,izone,Trac(izone),status)
					else
						Trac(izone)%v=Trac(izone)%v-minv
					endif
			end select
			Trac(izone)%recompute=.false.
		else
			select case(Zone(izone)%shape)
				case("SPH")
					if(Trac(izone)%recompute) then
						call HitSph(phot,izone,Trac(izone))
					else
						Trac(izone)%v=Trac(izone)%v-minv
					endif
			end select
			Trac(izone)%recompute=.false.
		endif
	enddo
	do istar=1,nstars
		if(TracStar(istar)%recompute) then
			call HitStar(phot,istar,TracStar(istar))
			TracStar(istar)%recompute=.false.
		else
			TracStar(istar)%v=TracStar(istar)%v-minv
		endif
	enddo

	if(status.gt.0) then
		call output("Something is wrong... Don't worry I'll try to fix it.")
		do izone=1,nzones
			Trac(izone)%recompute=.true.
		enddo
		do istar=1,nstars
			TracStar(istar)%recompute=.true.
		enddo
		minv=0d0
		call InWhichZones(phot)
		goto 1
	endif

	minv=20d0*maxR
	leave=.true.
	do izone=1,nzones
		if(Trac(izone)%v.gt.0d0.and.Trac(izone)%v.lt.minv) then
			minv=Trac(izone)%v
			imin=izone
			leave=.false.
		endif
	enddo
	hitstar0=.false.
	do istar=1,nstars
		if(TracStar(istar)%v.gt.0d0.and.TracStar(istar)%v.lt.minv) then
			minv=TracStar(istar)%v
			imin=0
			hitstar0=.true.
		endif
	enddo

	if(leave) goto 3
	
	if((tau0-phot%Kext*minv).lt.0d0) then
		minv=tau0/phot%Kext
		phot%edgeNr=0
		call increaseColumn(phot,minv)
		call TravelPhotonX(phot,minv)
		call AddEtrace(phot,minv)
		call Interact(phot)
		do izone=1,nzones
			Trac(izone)%recompute=.true.
		enddo
		do istar=1,nstars
			TracStar(istar)%recompute=.true.
		enddo
		tau0=-log(random(idum))
		if(random(idum).lt.fstop) then
			phot%sI=0d0
			goto 3
		endif
		phot%sI=phot%sI/(1d0-fstop)
		phot%sQ=phot%sQ/(1d0-fstop)
		phot%sU=phot%sU/(1d0-fstop)
		phot%sV=phot%sV/(1d0-fstop)
		goto 1
	endif

	call increaseColumn(phot,minv)
	call TravelPhotonX(phot,minv)
	call AddEtrace(phot,minv)
	tau0=tau0-phot%Kext*minv

	if(hitstar0) goto 3

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
			Trac(izone)%recompute=.true.
		else
			phot%edgeNr(izone)=0
		endif
	enddo

	goto 1
3	continue
	
	return
	end
	

	subroutine increaseColumn(phot,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 v
	integer izone,ipart,isize,iT
	
	do izone=1,nzones
		if(phot%inzone(izone)) then
		do ipart=1,npart
			do isize=1,Part(ipart)%nsize
				do iT=1,Part(ipart)%nT
					column(ipart,isize,iT)=column(ipart,isize,iT)
     &		+Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))%densP(ipart,isize,iT)*v
				enddo
			enddo
		enddo
		endif
	enddo
	
	return
	end


	subroutine AddEtrace(phot,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer izone
	real*8 v
	
	do izone=1,nzones
		if(phot%inzone(izone)) then
			Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))%Etrace=
     & 			Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))%Etrace+
     & 			v*phot%KabsZ(izone)*phot%sI
		endif
	enddo
  
	return
	end
	

	subroutine MCoutput(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 inp,tot,x,y,z,ext(nlam),fact
	integer i,j,ix,iy,ipart,isize,iT
	logical doMCobs(nMCobs),doany
	type(Photon) phot

	doany=.false.
	doMCobs=.false.
	do i=1,nMCobs
		inp=MCobs(i)%x*phot%vx+MCobs(i)%y*phot%vy+MCobs(i)%z*phot%vz
		if(inp.gt.MCobs(i)%opening) then
			doMCobs(i)=.true.
			doany=.true.
		endif
	enddo

	if(doany) then
		ext=0d0
		do ipart=1,npart
			do isize=1,Part(ipart)%nsize
				do iT=1,Part(ipart)%nT
					do j=1,nlam
						ext(j)=ext(j)+column(ipart,isize,iT)*Part(ipart)%Kext(isize,iT,j)
					enddo
				enddo
			enddo
		enddo
		do j=1,nlam
			if(ext(j).lt.1000d0) then
				specemit(j)=specemit(j)*exp(-ext(j))
			else
				specemit(j)=0d0
			endif
		enddo
		call integrate(specemit,tot)
		if(tot.gt.1d-100) then
			specemit=specemit/tot
		else
			doMCobs=.false.
		endif
	endif

!$OMP CRITICAL
	do i=1,nMCobs
		if(doMCobs(i)) then
			fact=1d23*phot%sI/MCobs(i)%f
			MCobs(i)%spec(1:nlam)=MCobs(i)%spec(1:nlam)+specemit(1:nlam)*fact
			x=phot%x0
			y=phot%y0
			z=phot%z0
			call rotateZ(x,y,z,MCobs(i)%cosp,-MCobs(i)%sinp)
			call rotateY(x,y,z,MCobs(i)%cost,MCobs(i)%sint)
			ix=real(MCobs(i)%npix)*(y+MCobs(i)%maxR)/(2d0*MCobs(i)%maxR)
			iy=MCobs(i)%npix-real(MCobs(i)%npix)*(x+MCobs(i)%maxR)/(2d0*MCobs(i)%maxR)
			if(ix.le.MCobs(i)%npix.and.iy.le.MCobs(i)%npix.and.ix.gt.0.and.iy.gt.0) then
				MCobs(i)%image(ix,iy,1:nlam)=MCobs(i)%image(ix,iy,1:nlam)+specemit(1:nlam)*fact
			endif
		endif
	enddo
!$OMP END CRITICAL
	
	return
	end
	

	subroutine EmitPhoton(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	integer i,izone
	real*8 Ltot,r,random,x,y,z,theta,x0,y0,z0,fbeam,xn,yn,zn
	
	fbeam=0.9d0
	
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

	r=random(idum)
	do izone=1,nzones
		r=r-Zone(izone)%fbeamS(i)
		if(r.lt.0d0) exit
	enddo
	if(izone.le.nzones) then
c beaming
		x0=Zone(izone)%x0-Star(i)%x
		y0=Zone(izone)%y0-Star(i)%y
		z0=Zone(izone)%z0-Star(i)%z
		r=sqrt(x0**2+y0**2+z0**2)
		x0=x0/r
		y0=y0/r
		z0=z0/r
		x=random(idum)
		y=random(idum)
		z=random(idum)
		xn=y0*z-z0*y
		yn=z0*x-x0*z
		zn=x0*y-y0*x
		r=sqrt(xn**2+yn**2+zn**2)
		xn=xn/r
		yn=yn/r
		zn=zn/r
1		continue
		if(random(idum).lt.fbeam) then
			theta=acos(1d0-(1d0-Zone(izone)%ctbeamS(i))*random(idum))
			phot%sI=phot%sI*Zone(izone)%EfbeamS(i)/fbeam
		else
			theta=acos(-1d0-(-1d0-Zone(izone)%ctbeamS(i))*random(idum))
			phot%sI=phot%sI*(1d0-Zone(izone)%EfbeamS(i))/(1d0-fbeam)
		endif
		phot%vx=x0
		phot%vy=y0
		phot%vz=z0
		call rotate(phot%vx,phot%vy,phot%vz,xn,yn,zn,theta)
		theta=2d0*pi*random(idum)
		call rotate(phot%vx,phot%vy,phot%vz,x0,y0,z0,theta)
		if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
			call randomdirection(phot%x,phot%y,phot%z)
			phot%x=Star(i)%R*phot%x
			phot%y=Star(i)%R*phot%y
			phot%z=Star(i)%R*phot%z
			phot%sI=Ltot/real(Nphot)
			goto 1
		endif
	else
		if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
			phot%vx=-phot%vx
			phot%vy=-phot%vy
			phot%vz=-phot%vz
		endif
	endif

	phot%x=phot%x+Star(i)%x
	phot%y=phot%y+Star(i)%y
	phot%z=phot%z+Star(i)%z
	
	phot%edgeNr=0
	phot%inzone=.false.

	phot%x0=phot%x
	phot%y0=phot%y
	phot%z0=phot%z

	call TranslatePhotonX(phot)
	call TranslatePhotonV(phot)
	
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
					C%Etrace=0d0
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
	real*8 function GetKp(iBB,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ipart,iT,isize,iBB
	type(Cell) C
	
	GetKp=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKp=GetKp+C%densP(ipart,isize,iT)*Part(ipart)%Kp(isize,iT,iBB)
			enddo
		enddo
	enddo
	
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
	
	
	
c This function gives the scattering cross section per cm in a certain angle
	real*8 function GetF11(ilam,iscat,C)
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,ilam,ipart,iT,isize,iscat
	type(Mueller) F
	type(Cell) C
	
	GetF11=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetF11=GetF11+C%densP(ipart,isize,iT)*
     &				Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F11(iscat)
			enddo
		enddo
	enddo
	
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
			r2=phot%xzone(izone)**2+phot%yzone(izone)**2+phot%zzone(izone)**2
		case("CYL")
			r2=phot%xzone(izone)**2+phot%yzone(izone)**2
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
	
	r=sqrt(phot%xzone(izone)**2+phot%yzone(izone)**2+phot%zzone(izone)**2)
	theta=acos(phot%zzone(izone)/r)
	
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
	
	r=sqrt(phot%xzone(izone)**2+phot%yzone(izone)**2)
	phi=acos(phot%xzone(izone)/r)
	if(phot%yzone(izone).lt.0d0) phi=2d0*pi-phi
	
	do i=1,Zone(izone)%np
		if(phi.ge.Zone(izone)%phi(i).and.phi.le.Zone(izone)%phi(i+1)) exit
	enddo
	phot%i3(izone)=i
	
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
	
	subroutine DetermineTemperatures
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,izone
	real*8 determineT
	type(Cell),pointer :: C
	
	do izone=1,nzones
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i1,i2,i3,C)
!$OMP& SHARED(Zone,izone)
!$OMP DO
		do i1=1,Zone(izone)%n1
			do i2=1,Zone(izone)%n2
				do i3=1,Zone(izone)%n3
					C=>Zone(izone)%C(i1,i2,i3)
					C%T=determineT(C)
				enddo
			enddo
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	enddo
	
	return
	end
	
	