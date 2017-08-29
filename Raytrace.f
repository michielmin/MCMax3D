	subroutine Raytrace(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,i,ilam,ilam0
	real*8 dlmin,fluxZ(nzones+nstars),fov0,zscale,Reddening,compute_dlam
	real*8 starttime,stoptime
	real*8,allocatable :: im(:,:,:)
	logical simpleobs
	character*500 MCfile

	nerrors=0
	simpleobs=.false.
	do ilam=1,nlam
		if(lam(ilam).ge.MCobs(iobs)%lam1.and.lam(ilam).le.MCobs(iobs)%lam2) then
			simpleobs=.true.
			exit
		endif
	enddo
	ilam0=0
	if(.not.simpleobs) then
		dlmin=lam(nlam)-lam(1)
		do ilam=1,nlam
			if(abs(lam(ilam)-(MCobs(iobs)%lam1+MCobs(iobs)%lam2)/2d0).lt.dlmin) then
				dlmin=abs(lam(ilam)-(MCobs(iobs)%lam1+MCobs(iobs)%lam2)/2d0)
				ilam0=ilam
			endif
		enddo
	endif

	call output("==================================================================")
	call output("Raytracing observation " // int2string(iobs,'(i0.4)'))
	call output(" wavelength: " // trim(dbl2string(MCobs(iobs)%lam1,'(f13.3)')) // " micron")
	if(MCobs(iobs)%lam2.ne.MCobs(iobs)%lam1) then
	call output(" to:         " // trim(dbl2string(MCobs(iobs)%lam2,'(f13.3)')) // " micron")
	endif
	call output(" fov:        " // trim(dbl2string(MCobs(iobs)%fov,'(f13.3)')) //  " arcsec")
	call output("             " // trim(dbl2string(2d0*MCobs(iobs)%maxR/AU,'(f13.3)')) // " AU")
	call output("==================================================================")
	call cpu_time(starttime)

	fov0=(2d0*MCobs(iobs)%maxR/AU)/(Distance/parsec)
	zscale=1d3*(MCobs(iobs)%npix/fov0)**2

	do ilam=1,nlam
c		MCobs(iobs)%spec(ilam)=sum(MCobs(iobs)%image(:,:,ilam))*Reddening(lam(ilam),compute_dlam(lam(ilam)),Av)/distance**2
		MCobs(iobs)%spec(ilam)=MCobs(iobs)%spec(ilam)*Reddening(lam(ilam),compute_dlam(lam(ilam)),Av)/distance**2
	enddo

	call SetupPaths(iobs)
	MCfile=trim(outputdir) // "RTSpec" // trim(int2string(iobs,'(i0.4)')) // trim(MCobs(iobs)%flag) // ".dat"
	open(unit=20,file=MCfile,RECL=6000)
	do ilam=1,nlam
		if((lam(ilam).ge.MCobs(iobs)%lam1.and.lam(ilam).le.MCobs(iobs)%lam2)
     &			.or.ilam.eq.ilam0) then
			call TraceScattField(iobs,ilam,MCobs(iobs)%Nphot)
			call FormalSolution(iobs,ilam,fluxZ)
			MCfile=trim(outputdir) // "RTout" // trim(int2string(iobs,'(i0.4)')) // "_" // 
     &				trim(int2string(int(lam(ilam)),'(i0.6)')) // trim(dbl2string(lam(ilam)-int(lam(ilam)),'(f0.2)')) // 
     &				trim(MCobs(iobs)%flag) // ".fits.gz"
			MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,1:4)=
     &			MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,1:4)*Reddening(lam(ilam),compute_dlam(lam(ilam)),Av)/distance**2
			MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,4)=sqrt(
     &			MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,2)**2+
     &			MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,3)**2)
			allocate(im(MCobs(iobs)%npix,MCobs(iobs)%npix,4))
			im=zscale*MCobs(iobs)%image
			call writefitsfile(MCfile,im,4,MCobs(iobs)%npix)
			deallocate(im)
			if(MCobs(iobs)%telescope) then
				call Convolution(MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,1:3),MCobs(iobs)%npix,lam(ilam),
     &					MCobs(iobs)%D,MCobs(iobs)%D2,MCobs(iobs)%SpW,fov0,MCobs(iobs)%width,MCobs(iobs)%snoise)
				MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,4)=sqrt(
     &				MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,2)**2+
     &				MCobs(iobs)%image(1:MCobs(iobs)%npix,1:MCobs(iobs)%npix,3)**2)
				MCfile=trim(outputdir) // "RToutObs" // trim(int2string(iobs,'(i0.4)')) // "_" // 
     &				trim(int2string(int(lam(ilam)),'(i0.6)')) // trim(dbl2string(lam(ilam)-int(lam(ilam)),'(f0.2)')) // 
     &				trim(MCobs(iobs)%flag) // ".fits.gz"
				allocate(im(MCobs(iobs)%npix,MCobs(iobs)%npix,4))
				im=zscale*MCobs(iobs)%image
				call writefitsfile(MCfile,im,4,MCobs(iobs)%npix)
				deallocate(im)
			endif
			MCobs(iobs)%spec(ilam)=sum(MCobs(iobs)%image(:,:,1))
			fluxZ=fluxZ*Reddening(lam(ilam),compute_dlam(lam(ilam)),Av)/distance**2
			write(20,*) lam(ilam),MCobs(iobs)%spec(ilam),fluxZ(1:nzones+nstars)
			call flush(20)
		endif
	enddo
	close(unit=20)

	call deallocatePaths

	call cpu_time(stoptime)
	call output("Raytrace time: " // trim(dbl2string(stoptime-starttime,'(f10.3)')) // " s")
	call output("==================================================================")
	
	return
	end

	subroutine TraceScattField(iobs,ilam,NphotMono)
	use GlobalSetup
	IMPLICIT NONE
	integer iobs,ilam,izone,iT,i,i1,i2,i3,iphot,istar,NphotMono
	real*8 GetKabs,Etot,Erandom,random,x,y,z,r
	logical emitfromstar
	type(Cell),pointer :: C
	type(Photon),allocatable :: phot(:)
	integer ispat,nspat
	integer,allocatable :: zspat(:),i1spat(:),i2spat(:),i3spat(:)
	real*8,allocatable :: Espat(:)
	integer maxnopenmp,iopenmp,omp_get_thread_num,omp_get_max_threads
	real*8 starttime,starttime_w,omp_get_wtime

	maxnopenmp=omp_get_max_threads()+1
	if(.not.use_multi) maxnopenmp=1
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

	Etot=0d0
	do i=1,nstars
		Etot=Etot+Star(i)%F(ilam)
	enddo
	
	call output("Tracing wavelength: " // dbl2string(lam(ilam),'(f10.4)') // " micron")
	call output("  number " // int2string(ilam,'(i6)'))
	
	nspat=0
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
		do i2=1,Zone(izone)%n2
		do i3=1,Zone(izone)%n3
			nspat=nspat+1
		enddo
		enddo
		enddo
	enddo

	allocate(Espat(nspat+1))
	allocate(zspat(nspat))
	allocate(i1spat(nspat))
	allocate(i2spat(nspat))
	allocate(i3spat(nspat))
	nspat=0
	Espat(1)=0d0
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
		do i2=1,Zone(izone)%n2
		do i3=1,Zone(izone)%n3
			C => Zone(izone)%C(i1,i2,i3)
			C%Escatt=0d0
			C%Qscatt=0d0
			C%Uscatt=0d0
			C%Vscatt=0d0
			iT=(C%T+0.5d0)/dTBB
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
			if(BB(ilam,iT).gt.0d0) then
				C%Elam=GetKabs(ilam,C)*BB(ilam,iT)*C%V
				if(C%Elam.gt.0d0) then
					Etot=Etot+C%Elam
					nspat=nspat+1
					Espat(nspat+1)=Espat(nspat)+C%Elam
					zspat(nspat)=izone
					i1spat(nspat)=i1
					i2spat(nspat)=i2
					i3spat(nspat)=i3
				endif
			endif
		enddo
		enddo
		enddo
	enddo

	call cpu_time(starttime)
	starttime_w=omp_get_wtime()


	do i=1,MultiNphotMono
		if(MultiNphotMono.gt.1) call output("Number " // trim(int2string(i,'(i4)')) // 
     &										" of " // trim(int2string(MultiNphotMono,'(i4)')))
!$OMP PARALLEL IF(use_multi.and.rt_multi)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(iphot,iopenmp,Erandom,emitfromstar,ispat,i1,i2,i3,istar,izone,x,y,z,r)
!$OMP& SHARED(phot,iobs,Espat,Star,MCobs,nspat,i1spat,i2spat,i3spat,zspat,Etot,nstars,ilam,NphotMono,
!$OMP&			starttime,starttime_w,MultiNphotMono,i)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do iphot=1,NphotMono
c		call tellertje(iphot,NphotMono)
		call tellertje_time(iphot+NphotMono*(i-1),NphotMono*MultiNphotMono,starttime,starttime_w)
		iopenmp=omp_get_thread_num()+1
		phot(iopenmp)%sI=Etot/real(NphotMono)/real(MultiNphotMono)
		phot(iopenmp)%sQ=0d0
		phot(iopenmp)%sU=0d0
		phot(iopenmp)%sV=0d0
		phot(iopenmp)%ilam1=ilam
		phot(iopenmp)%pol=.false.
2		Erandom=Etot*random(idum)
		emitfromstar=.false.
		do istar=1,nstars
			Erandom=Erandom-Star(istar)%F(ilam)
			if(Erandom.le.0d0) then
				emitfromstar=.true.
				goto 1
			endif
		enddo

		call hunt(Espat,nspat+1,Erandom,ispat)
		if(ispat.ge.nspat+1.or.ispat.le.0) goto 2
		izone=zspat(ispat)
		i1=i1spat(ispat)
		i2=i2spat(ispat)
		i3=i3spat(ispat)
1		continue
		if(emitfromstar) then
			call EmitPhotonStar(phot(iopenmp),istar)
		else
			call EmitPhotonMatter(phot(iopenmp),izone,i1,i2,i3)
		endif
		if(iobs.gt.0) then
			x=MCobs(iobs)%y*phot(iopenmp)%vz-MCobs(iobs)%z*phot(iopenmp)%vy
			y=MCobs(iobs)%z*phot(iopenmp)%vx-MCobs(iobs)%x*phot(iopenmp)%vz
			z=MCobs(iobs)%x*phot(iopenmp)%vy-MCobs(iobs)%y*phot(iopenmp)%vx
		else
			x=-phot(iopenmp)%vy
			y=phot(iopenmp)%vx
			z=0d0
		endif
		r=sqrt(x**2+y**2+z**2)
		phot(iopenmp)%Sx=x/r
		phot(iopenmp)%Sy=y/r
		phot(iopenmp)%Sz=z/r
		call InWhichZones(phot(iopenmp))
		call TravelPhotonMono(phot(iopenmp),iobs)
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	enddo
	
	return
	end


	subroutine EmitPhotonStar(phot,istar)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer istar,izone
	real*8 sI_in,r,x0,y0,z0,xn,yn,zn,x,y,z,theta,fbeam,random
	type(Photon) phot
	
	sI_in=phot%sI

	fbeam=0.9
	
	call randomdirection(phot%x,phot%y,phot%z)

	phot%x=Star(istar)%R*phot%x
	phot%y=Star(istar)%R*phot%y
	phot%z=Star(istar)%R*phot%z

	call randomdirection(phot%vx,phot%vy,phot%vz)

	r=random(idum)
	do izone=1,nzones
		r=r-Zone(izone)%fbeamS(istar)
		if(r.lt.0d0) exit
	enddo
	if(izone.le.nzones) then
c beaming
		x0=Zone(izone)%x0-Star(istar)%x
		y0=Zone(izone)%y0-Star(istar)%y
		z0=Zone(izone)%z0-Star(istar)%z
		r=sqrt(x0**2+y0**2+z0**2)
		if(r.lt.Zone(izone)%Rout) then
			if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
				phot%vx=-phot%vx
				phot%vy=-phot%vy
				phot%vz=-phot%vz
			endif
		else
		x0=x0/r
		y0=y0/r
		z0=z0/r
2		x=random(idum)
		y=random(idum)
		z=random(idum)
		r=sqrt(x**2+y**2+z**2)
		if(r.gt.1d0) goto 2
		xn=y0*z-z0*y
		yn=z0*x-x0*z
		zn=x0*y-y0*x
		r=sqrt(xn**2+yn**2+zn**2)
		xn=xn/r
		yn=yn/r
		zn=zn/r
1		continue
		if(random(idum).lt.fbeam) then
			theta=acos(1d0-(1d0-Zone(izone)%ctbeamS(istar))*random(idum))
			phot%sI=phot%sI*Zone(izone)%EfbeamS(istar)/fbeam
		else
			theta=acos(-1d0-(-1d0-Zone(izone)%ctbeamS(istar))*random(idum))
			phot%sI=phot%sI*(1d0-Zone(izone)%EfbeamS(istar))/(1d0-fbeam)
		endif
		phot%vx=x0
		phot%vy=y0
		phot%vz=z0
		call rotate(phot%vx,phot%vy,phot%vz,xn,yn,zn,theta)
		theta=2d0*pi*random(idum)
		call rotate(phot%vx,phot%vy,phot%vz,x0,y0,z0,theta)
		if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
			call randomdirection(phot%x,phot%y,phot%z)
			phot%x=Star(istar)%R*phot%x
			phot%y=Star(istar)%R*phot%y
			phot%z=Star(istar)%R*phot%z
			phot%sI=sI_in
			goto 1
		endif
		endif
	else
		if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
			phot%vx=-phot%vx
			phot%vy=-phot%vy
			phot%vz=-phot%vz
		endif
	endif

	phot%x=phot%x+Star(istar)%x
	phot%y=phot%y+Star(istar)%y
	phot%z=phot%z+Star(istar)%z
	
	phot%edgeNr=0
	phot%inzone=.false.

	call TranslatePhotonX(phot)
	call TranslatePhotonV(phot)

	return
	end
	
	
	subroutine EmitPhotonMatter(phot,izone,i1,i2,i3)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3	
	real*8 R,theta,phi,random
	type(Photon) phot
	
	select case(Zone(izone)%shape)
		case("SPH")
			R=Zone(izone)%R(i1)+(Zone(izone)%R(i1+1)-Zone(izone)%R(i1))*random(idum)
			theta=Zone(izone)%theta(i2)+(Zone(izone)%theta(i2+1)-Zone(izone)%theta(i2))*random(idum)
			phi=Zone(izone)%phi(i3)+(Zone(izone)%phi(i3+1)-Zone(izone)%phi(i3))*random(idum)
			phot%x=R*cos(phi)*sin(theta)
			phot%y=R*sin(phi)*sin(theta)
			phot%z=R*cos(theta)
		case default
			call output("Raytracing on non-spherical zones not yet possible")
			stop
	end select

	call randomdirection(phot%vx,phot%vy,phot%vz)

	call TranslatePhotonXinverse(phot,izone)

	phot%edgeNr=0
	phot%inzone=.false.

	call TranslatePhotonX(phot)
	call TranslatePhotonV(phot)
	
	return
	end
	
	
	

	subroutine TravelPhotonMono(phot,iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,imin,iobs,istar,status
	logical leave,inany,hitstar0
	real*8 minv,tau0,tau,GetKext,random,GetKabs,fstopmono,albedo,theta,sin2t,cos2t
	type(Travel) Trac(nzones)
	type(Travel) TracStar(nstars)
	type(Cell),pointer :: C
	type(Photon) phot
	real*8 x,y,z,r

	tau0=-log(random(idum))
	
	if(iobs.gt.0) then
		theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
		phot%iscat=180d0*theta/pi
		if(phot%iscat.lt.1) phot%iscat=1
		if(phot%iscat.gt.180) phot%iscat=180
		call MakeRotateStokes(phot,MCobs(iobs)%xup,MCobs(iobs)%yup,MCobs(iobs)%zup,
     &		MCobs(iobs)%x,MCobs(iobs)%y,MCobs(iobs)%z,sin2t,cos2t)
	endif
	
1	continue

	phot%Kext=0d0
	phot%Kabs=0d0
	inany=.false.
	do izone=1,nzones
		phot%KabsZ(izone)=0d0
		if(phot%inzone(izone)) then
			inany=.true.
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			phot%Kext=phot%Kext+GetKext(phot%ilam1,C)
			phot%KabsZ(izone)=GetKabs(phot%ilam1,C)
			phot%Kabs=phot%Kabs+phot%KabsZ(izone)
		endif
	enddo
	if(inany) then
		albedo=(phot%Kext-phot%Kabs)/phot%Kext
		fstopmono=1d0-albedo**0.25
	else
		albedo=0d0
		fstopmono=1d0
	endif

	status=0
	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone),status)
			end select
		else
			select case(Zone(izone)%shape)
				case("SPH")
					call HitSph(phot,izone,Trac(izone))
			end select
		endif
	enddo
	do istar=1,nstars
		call HitStar(phot,istar,TracStar(istar))
	enddo
	if(status.gt.0) then
		call output("Something is wrong...")
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
		call TravelPhotonX(phot,minv)
		call AddEtraceMono(phot,minv,sin2t,cos2t,(iobs.gt.0))
		call InteractMono(phot)
		if(iobs.gt.0) then
			x=MCobs(iobs)%y*phot%vz-MCobs(iobs)%z*phot%vy
			y=MCobs(iobs)%z*phot%vx-MCobs(iobs)%x*phot%vz
			z=MCobs(iobs)%x*phot%vy-MCobs(iobs)%y*phot%vx
			r=sqrt(x**2+y**2+z**2)
			x=x/r
			y=y/r
			z=z/r
			call RotateStokes(phot,x,y,z)
			phot%Sx=x
			phot%Sy=y
			phot%Sz=z
			call MakeRotateStokes(phot,MCobs(iobs)%xup,MCobs(iobs)%yup,MCobs(iobs)%zup,
     &			MCobs(iobs)%x,MCobs(iobs)%y,MCobs(iobs)%z,sin2t,cos2t)
			theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
			phot%iscat=180d0*theta/pi
			if(phot%iscat.lt.1) phot%iscat=1
			if(phot%iscat.gt.180) phot%iscat=180
		endif
		tau0=-log(random(idum))
		if(random(idum).lt.fstopmono) then
			goto 3
		endif
		phot%sI=phot%sI*albedo/(1d0-fstopmono)
		phot%sQ=phot%sQ*albedo/(1d0-fstopmono)
		phot%sU=phot%sU*albedo/(1d0-fstopmono)
		phot%sV=phot%sV*albedo/(1d0-fstopmono)
		goto 1
	endif

	call TravelPhotonX(phot,minv)
	if(inany) call AddEtraceMono(phot,minv,sin2t,cos2t,(iobs.gt.0))
	tau0=tau0-phot%Kext*minv

	if(hitstar0) goto 3

	do izone=1,nzones
		if(Trac(izone)%v.le.minv.or.izone.eq.imin) then
			if(Zone(izone)%warped.and.phot%i1(izone).ne.Trac(izone)%i1next) then
				phot%i1(izone)=Trac(izone)%i1next
				call TranslatePhotonWarp(phot,izone)
				phot%i2(izone)=-1
				phot%i3(izone)=-1
			else
				phot%i1(izone)=Trac(izone)%i1next
				phot%i2(izone)=Trac(izone)%i2next
				phot%i3(izone)=Trac(izone)%i3next
			endif

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

	goto 1
3	continue
	
	return
	end



	subroutine AddEtraceMono(phot,v,sin2t,cos2t,doangle)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	integer izone,ipart,isize,iT,iscat,ilam
	real*8 v,GetF11,F11,F12,F22,F33,F34,F44,f,sI,sQ,sU,sV,Qt,Ut
	real*8 sin2t,cos2t
	type(Cell),pointer :: C
	logical doangle
	
	iscat=phot%iscat
	ilam=phot%ilam1
	do izone=1,nzones
		if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			if(doangle) then
			F11=0d0
			F12=0d0
			F22=0d0
			F33=0d0
			F34=0d0
			F44=0d0
			do ipart=1,npart
				do isize=1,Part(ipart)%nsize
					do iT=1,Part(ipart)%nT
						f=C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
						F11=F11+f*Part(ipart)%F(isize,iT,ilam)%F11(iscat)
						F12=F12+f*Part(ipart)%F(isize,iT,ilam)%F12(iscat)
						F22=F22+f*Part(ipart)%F(isize,iT,ilam)%F22(iscat)
						F33=F33+f*Part(ipart)%F(isize,iT,ilam)%F33(iscat)
						F34=F34+f*Part(ipart)%F(isize,iT,ilam)%F34(iscat)
						F44=F44+f*Part(ipart)%F(isize,iT,ilam)%F44(iscat)
					enddo
				enddo
			enddo
			sI=F11*phot%sI+F12*phot%sQ
			Qt=F12*phot%sI+F22*phot%sQ
			Ut=F34*phot%sV+F33*phot%sU
			sQ=Qt*cos2t+Ut*sin2t
			sU=-Qt*sin2t+Ut*cos2t
			sV=-F34*phot%sU+F44*phot%sV
			C%Escatt=C%Escatt+v*sI
			C%Qscatt=C%Qscatt+v*sQ
			C%Uscatt=C%Uscatt+v*sU
			C%Vscatt=C%Vscatt+v*sV
			else
			C%Escatt=C%Escatt+v*phot%sI
			endif
		endif
	enddo

	return
	end
	


	subroutine InteractMono(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 GetKabs,random,K,GetKp,spec(nlam),KscaP,GetF11,KscaR
	real*8 T0,T1,epsT0,epsT1,kp,tot,GetKsca,increaseT,Ksca(nlam)
	integer izone,iT0,iT1,l,ipart,isize,iT,iscat
	type(Mueller) M
	type(Cell),pointer :: C
	
2	KscaR=(phot%Kext-phot%Kabs)*random(idum)
	do izone=1,nzones
		if(phot%inzone(izone)) then
		C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
		do ipart=1,npart
			do isize=1,Part(ipart)%nsize
				do iT=1,Part(ipart)%nT
					KscaP=C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,phot%ilam1)
					KscaR=KscaR-KscaP
					if(KscaR.lt.0d0) then
						M=Part(ipart)%F(isize,iT,phot%ilam1)
						goto 1
					endif
				enddo
			enddo
		enddo
		endif
	enddo
	goto 2
1	continue

	call scatangle(phot,M,iscat)

	call TranslatePhotonV(phot)

	return
	end
	



	subroutine SetupPaths(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,izone,ir,ip,i,j,ix,iy,ii
	real*8 random,phi,x,y,z,flux,A,R1,R2,Rad,scale(3)
	real*8,allocatable :: Rtau1(:)
	type(ZoneType),pointer :: Zo
	type(StarType),pointer :: St
	type(Photon),allocatable :: phot(:),phot0(:)
	integer maxnopenmp,iopenmp,omp_get_thread_num,omp_get_max_threads,maxcount
	logical error
	integer,allocatable :: t_node(:,:),t_neighbor(:,:)
	integer matri,nnodes
	real*8,allocatable :: xy(:,:)
	real*8 matrix(2,2)

	maxnopenmp=omp_get_max_threads()+1
	if(.not.use_multi) maxnopenmp=1
	allocate(phot(maxnopenmp))	
	allocate(phot0(maxnopenmp))	
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
		allocate(phot0(i)%i1(nzones))
		allocate(phot0(i)%i2(nzones))
		allocate(phot0(i)%i3(nzones))
		allocate(phot0(i)%inzone(nzones))
		allocate(phot0(i)%edgeNr(nzones))
		allocate(phot0(i)%KabsZ(nzones))
		allocate(phot0(i)%xzone(nzones))
		allocate(phot0(i)%yzone(nzones))
		allocate(phot0(i)%zzone(nzones))
		allocate(phot0(i)%vxzone(nzones))
		allocate(phot0(i)%vyzone(nzones))
		allocate(phot0(i)%vzzone(nzones))
	enddo
	
	allocate(Pimage(nzones+nstars))

	maxcount=0
	do izone=1,nzones
		maxcount=maxcount+2d0*Zone(izone)%n1*Zone(izone)%n2*Zone(izone)%n3
	enddo
	
	nnodes=0
	do izone=1,nzones+nstars
		if(izone.le.nzones) then
			Zo => Zone(izone)
			allocate(Rtau1(Zo%nt))
			call ComputeRtau1(Rtau1,Zo%nt,izone)
			x=Zo%x0
			y=Zo%y0
			z=Zo%z0
			Pimage(izone)%nr=3*(2*Zo%nR+2*Zo%nt+200)
			Pimage(izone)%np=min(max(Zo%np*4,50),MCobs(iobs)%np)
			allocate(Pimage(izone)%R(Pimage(izone)%nr))
			allocate(Pimage(izone)%P(Pimage(izone)%nr,Pimage(izone)%np))
			ir=0
			scale(1)=Zo%xscale
			scale(2)=Zo%yscale
			scale(3)=Zo%zscale
			do i=1,150,MCobs(iobs)%nr
				ir=ir+1
				Pimage(izone)%R(ir)=Zo%R(1)*real(i)/151d0
			enddo
			do ii=1,3
			do j=1,Zo%nt,MCobs(iobs)%nr
				if(Rtau1(j).gt.0d0) then
					ir=ir+1
					Pimage(izone)%R(ir)=abs(scale(ii)*(Rtau1(j)*sin((Zo%theta(j)+Zo%theta(j+1))/2d0)))
					ir=ir+1
					Pimage(izone)%R(ir)=abs(scale(ii)*(Rtau1(j)*sin(MCobs(iobs)%theta-Zo%theta0+(Zo%theta(j)+Zo%theta(j+1))/2d0)))
				endif
			enddo
			do i=1,Zo%nR,MCobs(iobs)%nr
				ir=ir+1
				Pimage(izone)%R(ir)=abs(scale(ii)*sqrt(Zo%R(i)*Zo%R(i+1)))
				ir=ir+1
				Pimage(izone)%R(ir)=abs(scale(ii)*sqrt(Zo%R(i)*Zo%R(i+1))*sin(MCobs(iobs)%theta-Zo%theta0))
			enddo
			ir=ir+1
			Pimage(izone)%R(ir)=Zo%Rout
			enddo
			Pimage(izone)%nr=ir
			deallocate(Rtau1)
		else
			St => Star(izone-nzones)
			x=St%x
			y=St%y
			z=St%z
			Pimage(izone)%nr=150
			Pimage(izone)%np=min(150,MCobs(iobs)%np)
			allocate(Pimage(izone)%R(Pimage(izone)%nr))
			allocate(Pimage(izone)%P(Pimage(izone)%nr,Pimage(izone)%np))
			do i=1,Pimage(izone)%nr
				Pimage(izone)%R(i)=St%R*(real(i)-0.9)/(real(Pimage(izone)%nr)-0.9)
			enddo
		endif

1		call sort(Pimage(izone)%R,Pimage(izone)%nr)
		do ir=1,Pimage(izone)%nr-1
			if(Pimage(izone)%R(ir).eq.Pimage(izone)%R(ir+1)) then
				do i=ir+1,Pimage(izone)%nr-1
					Pimage(izone)%R(i)=Pimage(izone)%R(i+1)
				enddo
				Pimage(izone)%nr=Pimage(izone)%nr-1
				goto 1
			endif
		enddo

		call rotateZ(x,y,z,MCobs(iobs)%cosp,-MCobs(iobs)%sinp)
		call rotateY(x,y,z,MCobs(iobs)%cost,MCobs(iobs)%sint)
		Pimage(izone)%x=x
		Pimage(izone)%y=y
		if(izone.le.nzones) then
			call output("Paths for zone " // trim(int2string(izone,'(i4)')))
		else
			call output("Paths for star " // trim(int2string(izone-nzones,'(i4)')))
		endif

		call output("Number of rays: nr x nphi = " // trim(int2string(Pimage(izone)%nr,'(i4)')) // "x" 
     &				// trim(int2string(Pimage(izone)%np,'(i4)')) // " = "  
     &				// trim(int2string(Pimage(izone)%nr*Pimage(izone)%np,'(i8)')))
		call tellertje(1,100)
		do ir=1,Pimage(izone)%nr
			iopenmp=omp_get_thread_num()+1
			call tellertje(ir+1,Pimage(izone)%nr+2)
			do ip=1,Pimage(izone)%np
				nnodes=nnodes+1
				if(ir.eq.1) then
					R1=Pimage(izone)%R(ir)
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
					R2=R1
				else if(ir.eq.Pimage(izone)%nr) then
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=Pimage(izone)%R(ir)
					R1=R2
				else
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
				endif
				Rad=R1+random(idum)*(R2-R1)
				phi=2d0*pi*(real(ip-1)+random(idum))/real(Pimage(izone)%np)

				phot(iopenmp)%x=Rad*cos(phi)+x
				phot(iopenmp)%y=Rad*sin(phi)+y
				phot(iopenmp)%z=z

				Pimage(izone)%P(ir,ip)%x=phot(iopenmp)%x
				Pimage(izone)%P(ir,ip)%y=phot(iopenmp)%y
			enddo
		enddo
		call tellertje(100,100)
	enddo
	
	allocate(xy(2,nnodes))
	i=0
	do izone=1,nzones+nstars
		do ir=1,Pimage(izone)%nr
			do ip=1,Pimage(izone)%np
				i=i+1
				xy(1,i)=Pimage(izone)%P(ir,ip)%x
				xy(2,i)=Pimage(izone)%P(ir,ip)%y
			enddo
		enddo
	enddo
	matrix(1,1)=1d0
	matrix(1,2)=0d0
	matrix(2,1)=0d0
	matrix(2,2)=1d0
	matri=2*nnodes-3
	allocate(t_node(3,matri))
	allocate(t_neighbor(3,matri))
	call dtris2_lmap (nnodes, xy, matrix, nTracePaths, t_node, t_neighbor )

	call output("Number of rays: " // trim(int2string(nTracePaths,'(i8)')))

	allocate(TracePaths(nTracePaths))

!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,ip,R1,R2,i,Rad,phi,iopenmp,error)
!$OMP& SHARED(MCobs,iobs,phot,phot0,maxR,nzones,maxcount,TracePaths,xy,nTracePaths,t_node)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do i=1,nTracePaths
		iopenmp=omp_get_thread_num()+1
		call tellertje(i+1,nTracePaths+2)

		phot(iopenmp)%x=(xy(1,t_node(1,i))+xy(1,t_node(2,i))+xy(1,t_node(3,i)))/3d0
		phot(iopenmp)%y=(xy(2,t_node(1,i))+xy(2,t_node(2,i))+xy(2,t_node(3,i)))/3d0
		phot(iopenmp)%z=0d0

		TracePaths(i)%x=phot(iopenmp)%x
		TracePaths(i)%y=phot(iopenmp)%y
		do j=1,3
			TracePaths(i)%xedge(j)=xy(1,t_node(j,i))
			TracePaths(i)%yedge(j)=xy(2,t_node(j,i))
		enddo

		call rotateY(phot(iopenmp)%x,phot(iopenmp)%y,phot(iopenmp)%z,cos(MCobs(iobs)%theta),-sin(MCobs(iobs)%theta))
		call rotateZ(phot(iopenmp)%x,phot(iopenmp)%y,phot(iopenmp)%z,cos(MCobs(iobs)%phi),sin(MCobs(iobs)%phi))
		phot(iopenmp)%vx=0d0
		phot(iopenmp)%vy=0d0
		phot(iopenmp)%vz=-1d0
		call rotateY(phot(iopenmp)%vx,phot(iopenmp)%vy,phot(iopenmp)%vz,cos(MCobs(iobs)%theta),-sin(MCobs(iobs)%theta))
		call rotateZ(phot(iopenmp)%vx,phot(iopenmp)%vy,phot(iopenmp)%vz,cos(MCobs(iobs)%phi),sin(MCobs(iobs)%phi))

		call TranslatePhotonX(phot(iopenmp))
		call TranslatePhotonV(phot(iopenmp))

		call TravelPhotonX(phot(iopenmp),-3d0*maxR)

		phot0(iopenmp)=phot(iopenmp)
		call RaytracePath(phot(iopenmp),TracePaths(i),.true.,maxcount,error)

		allocate(TracePaths(i)%v(TracePaths(i)%n+1))
		allocate(TracePaths(i)%inzone(TracePaths(i)%n+1,nzones))
		allocate(TracePaths(i)%C(TracePaths(i)%n+1,nzones))

		phot(iopenmp)=phot0(iopenmp)
		call RaytracePath(phot(iopenmp),TracePaths(i),.false.,maxcount,error)
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)


	return
	end


	subroutine ComputeRtau1(Rtau1,nt,izone)
	use GlobalSetup
	IMPLICIT NONE
	integer nt,izone,it,ilam,ir,ip
	real*8 Rtau1(nt),tau,GetKext,tau_tot,d,R

	do ilam=1,nlam-1
		if(lam(ilam).lt.0.55.and.lam(ilam+1).ge.0.55) exit
	enddo
	
	do it=1,nt
		Rtau1(it)=0d0
		do ip=1,Zone(izone)%np
		R=-1d0
		tau_tot=0d0
		do ir=1,Zone(izone)%nr
			d=Zone(izone)%R(ir+1)-Zone(izone)%R(ir)
			tau=d*GetKext(ilam,Zone(izone)%C(ir,it,ip))
			if(tau_tot+tau.gt.1d0) then
				R=Zone(izone)%R(ir)+(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))*(1d0-tau_tot)/tau
				goto 1
			endif
			tau_tot=tau_tot+tau
		enddo
1		continue
		if(R.lt.Rtau1(it)) Rtau1(it)=R
		enddo
	enddo
		
	return
	end
	

	subroutine RaytraceFluxZone(P,flux,Q,U,V,ilam,izone0,traceall)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Path) :: P
	real*8 flux,Q,U,V
	integer ilam,izone0,istar,k
	type(Cell),pointer :: C

	integer izone,imin,iobs,iT
	logical leave,inany,hitstar0,traceall
	real*8 minv,tau0,tau,GetKext,random,GetKabs,theta,tau_e,fact,exptau,frac
	real*8 albedo,emis,scat(4),Kext,Kabs,KabsZ(nzones),Ksca

	tau0=0d0
	fact=1d0
	flux=0d0
	Q=0d0
	U=0d0
	V=0d0

	P%HitZone=.false.
	
	do k=1,P%n

	Kext=0d0
	Kabs=0d0
	inany=.false.
	do izone=1,nzones
		KabsZ(izone)=0d0
		if(P%inzone(k,izone)) then
			inany=.true.
			C => P%C(k,izone)%C
			Kext=Kext+GetKext(ilam,C)
			KabsZ(izone)=GetKabs(ilam,C)
			Kabs=Kabs+KabsZ(izone)
		endif
	enddo
	Ksca=(Kext-Kabs)
	albedo=(Kext-Kabs)/Kext

	if(inany) then
		tau_e=Kext*P%v(k)

		emis=0d0
		scat=0d0
		if(traceall) then
			do izone=1,nzones
			if(P%inzone(k,izone)) then
				if(izone.eq.izone0) P%HitZone=.true.
				C => P%C(k,izone)%C
				iT=(C%T+0.5d0)/dTBB
				if(iT.lt.1) iT=1
				if(iT.gt.nBB) iT=nBB
				emis=emis+BB(ilam,iT)*KabsZ(izone)
				scat(1)=scat(1)+C%Escatt/C%V
				scat(2)=scat(2)+C%Qscatt/C%V
				scat(3)=scat(3)+C%Uscatt/C%V
				scat(4)=scat(4)+C%Vscatt/C%V
			endif
			enddo
		else
			if(P%inzone(k,izone0)) then
				P%HitZone=.true.
				C => P%C(k,izone0)%C
				iT=(C%T+0.5d0)/dTBB
				if(iT.lt.1) iT=1
				if(iT.gt.nBB) iT=nBB
				emis=emis+BB(ilam,iT)*KabsZ(izone0)
				scat(1)=scat(1)+C%Escatt/C%V
				scat(2)=scat(2)+C%Qscatt/C%V
				scat(3)=scat(3)+C%Uscatt/C%V
				scat(4)=scat(4)+C%Vscatt/C%V
			endif
		endif
		emis=emis*(1d0-albedo)/Kabs
		scat=scat/Kext

		if(tau_e.lt.1d-6) then
			flux=flux+(scat(1)+emis)*tau_e*fact
			Q=Q+scat(2)*tau_e*fact
			U=U+scat(3)*tau_e*fact
			V=V+scat(4)*tau_e*fact
			fact=fact*(1d0-tau_e)
		else
			exptau=exp(-tau_e)
			frac=(1d0-exptau)
			flux=flux+(scat(1)+emis)*frac*fact
			Q=Q+scat(2)*frac*fact
			U=U+scat(3)*frac*fact
			V=V+scat(4)*frac*fact
			fact=fact*exptau
		endif
		tau0=tau0+tau_e
	endif

	enddo

	if(P%istar.gt.0.and.traceall) then
		flux=flux+(Star(P%istar)%F(ilam)/(pi*Star(P%istar)%R**2))*exp(-tau0)
		if((izone0-nzones).eq.P%istar) P%HitZone=.true.
	endif

	return
	end



	subroutine RaytraceFluxStar(P,flux,ilam,istar0)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Path) :: P
	real*8 flux
	integer ilam,izone0,istar0
	type(Cell),pointer :: C

	integer izone,imin,iobs,k
	logical leave,inany,hitstar0
	real*8 minv,tau0,tau,GetKext,random,theta,tau_e,exptau
	real*8 Kext

	tau0=0d0
	flux=0d0

	if(P%istar.ne.istar0) return
	P%HitZone=.true.
	
	do k=1,P%n
	
	Kext=0d0
	inany=.false.
	do izone=1,nzones
		if(P%inzone(k,izone)) then
			inany=.true.
			C => P%C(k,izone)%C
			Kext=Kext+GetKext(ilam,C)
		endif
	enddo

	if(inany) then
		tau_e=Kext*P%v(k)
		tau0=tau0+tau_e	
	endif

	enddo

	flux=(Star(istar0)%F(ilam)/(pi*Star(istar0)%R**2))*exp(-tau0)

	return
	end




	subroutine RaytracePath(phot,P,justcount,maxcount,error)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	type(Path) :: P
	integer istar,maxcount
	logical justcount,error

	integer izone,imin,iobs,k,status
	logical leave,inany,hitstar0
	real*8 minv,theta
	type(Travel) Trac(nzones),TracStar(nstars)
	type(Cell),pointer :: C

	error=.false.
	phot%inzone=.false.
	phot%edgeNr=0
	P%istar=0

	do izone=1,nzones
		Trac(izone)%recompute=.true.
	enddo
	do istar=1,nstars
		TracStar(istar)%recompute=.true.
	enddo
	minv=0d0

	if(justcount) P%n=0

	k=0
1	continue

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
		nerrors=nerrors+1
		if(nerrors.gt.100) stop
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

	minv=1d8*maxR
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
			P%istar=istar
		endif
	enddo

	k=k+1
	if(k.gt.maxcount) then
		error=.true.
		if(justcount) P%n=k		
		return
	endif

	if(.not.justcount) then
		do izone=1,nzones
			P%inzone(k,izone)=phot%inzone(izone)
			if(P%inzone(k,izone)) then
				P%C(k,izone)%C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			endif
		enddo
		P%v(k)=minv
	endif

	if(leave) goto 3
		
	call TravelPhotonX(phot,minv)

	if(hitstar0) goto 3

	do izone=1,nzones
		if(Trac(izone)%v.le.minv.or.izone.eq.imin) then
			if(Zone(izone)%warped.and.phot%i1(izone).ne.Trac(izone)%i1next) then
				phot%i1(izone)=Trac(izone)%i1next
				call TranslatePhotonWarp(phot,izone)
				phot%i2(izone)=-1
				phot%i3(izone)=-1
			else
				phot%i1(izone)=Trac(izone)%i1next
				phot%i2(izone)=Trac(izone)%i2next
				phot%i3(izone)=Trac(izone)%i3next
			endif

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

	if(justcount) P%n=k

	return
	end




	subroutine FormalSolution(iobs,ilam,fluxZ)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,ilam,izone,ir,ip,i,j,ix,iy,ii,istar,nint,nthreads,ithread,id
	integer omp_get_max_threads,omp_get_thread_num
	real*8 random,phi,x,y,flux,A,R1,R2,Rad,scale(3),fluxZ(nzones+nstars),Q,U,V
	real*8 l1,l2,l3,tot,w1,w2,w3,pl,vv2,vv3,v1v2,v1v3,v2v3,aaa,bbb,v2x,v2y,v3x,v3y
	real*8,allocatable :: imageopenmp(:,:,:,:,:)

	MCobs(iobs)%image(:,:,1:4)=0d0		! use the wavelength dimension for the Stokes vectors

	nthreads=omp_get_max_threads()
	allocate(imageopenmp(nthreads,nzones+nstars,MCobs(iobs)%npix,MCobs(iobs)%npix,4))
	imageopenmp=0d0
	
	do izone=1,nzones+nstars
		call tellertje(1,100)
!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,flux,A,i,j,x,y,ix,iy,nint,Q,U,V,ithread,w1,w2,w3,tot,pl,l1,l2,l3,istar,
!$OMP& 			vv2,vv3,v1v2,v1v3,v2v3,aaa,bbb,v2x,v2y,v3x,v3y)
!$OMP& SHARED(TracePaths,izone,ilam,MCobs,iobs,nTracePaths,transmissiontracing,Zone,imageopenmp,nzones)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
		do ir=1,nTracePaths
			ithread=omp_get_thread_num()+1
			call tellertje(ir+1,nTracePaths+2)
c			call RaytraceFluxZone(TracePaths(ir),flux,Q,U,V,ilam,izone,.true.)
			if(izone.le.nzones) then
				call RaytraceFluxZone(TracePaths(ir),flux,Q,U,V,ilam,izone,.false.)
			else
				istar=izone-nzones
				call RaytraceFluxStar(TracePaths(ir),flux,ilam,istar)
				flux=flux*MCobs(iobs)%fstar
				Q=0d0
				U=0d0
				V=0d0
			endif

			if(TracePaths(ir)%HitZone) then
			l1=sqrt((TracePaths(ir)%xedge(1)-TracePaths(ir)%xedge(2))**2+(TracePaths(ir)%yedge(1)-TracePaths(ir)%yedge(2))**2)
			l2=sqrt((TracePaths(ir)%xedge(1)-TracePaths(ir)%xedge(3))**2+(TracePaths(ir)%yedge(1)-TracePaths(ir)%yedge(3))**2)
			l3=sqrt((TracePaths(ir)%xedge(3)-TracePaths(ir)%xedge(2))**2+(TracePaths(ir)%yedge(3)-TracePaths(ir)%yedge(2))**2)
			pl=(l1+l2+l3)/2d0
			A=sqrt(pl*(pl-l1)*(pl-l2)*(pl-l3))

			nint=25d0*A*(real(MCobs(iobs)%npix)/(2d0*MCobs(iobs)%maxR))**2
			if(nint.lt.25) nint=25
			if(nint.gt.1000) nint=1000
			flux=1d23*flux*A/(real(nint)*4d0*pi)
			Q=-1d23*Q*A/(real(nint)*4d0*pi)
			U=-1d23*U*A/(real(nint)*4d0*pi)
			V=1d23*V*A/(real(nint)*4d0*pi)
			do i=1,nint
1				w1=random(idum)
				w2=random(idum)
				if((w1+w2).gt.1d0) goto 1
				v2x=TracePaths(ir)%xedge(2)-TracePaths(ir)%xedge(1)
				v3x=TracePaths(ir)%xedge(3)-TracePaths(ir)%xedge(1)
				v2y=TracePaths(ir)%yedge(2)-TracePaths(ir)%yedge(1)
				v3y=TracePaths(ir)%yedge(3)-TracePaths(ir)%yedge(1)
				x=TracePaths(ir)%xedge(1)+w1*v2x+w2*v3x
				y=TracePaths(ir)%yedge(1)+w1*v2y+w2*v3y

				x=real(MCobs(iobs)%npix)-real(MCobs(iobs)%npix)*(x+MCobs(iobs)%maxR)/(2d0*MCobs(iobs)%maxR)+1d0
				y=real(MCobs(iobs)%npix)*(y+MCobs(iobs)%maxR)/(2d0*MCobs(iobs)%maxR)+1d0
				ix=y
				iy=x
				if(ix.le.MCobs(iobs)%npix.and.iy.le.MCobs(iobs)%npix.and.ix.gt.0.and.iy.gt.0) then
					imageopenmp(ithread,izone,ix,iy,1)=imageopenmp(ithread,izone,ix,iy,1)+flux
					imageopenmp(ithread,izone,ix,iy,2)=imageopenmp(ithread,izone,ix,iy,2)+Q
					imageopenmp(ithread,izone,ix,iy,3)=imageopenmp(ithread,izone,ix,iy,3)+U
					imageopenmp(ithread,izone,ix,iy,4)=imageopenmp(ithread,izone,ix,iy,4)+V
				endif
			enddo
			endif
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		call tellertje(100,100)
	enddo
	fluxZ=0d0
	do ix=1,MCobs(iobs)%npix
		do iy=1,MCobs(iobs)%npix
			MCobs(iobs)%image(ix,iy,1:4)=0d0
			do izone=1,nzones+nstars
			do ithread=1,nthreads
				MCobs(iobs)%image(ix,iy,1:4)=MCobs(iobs)%image(ix,iy,1:4)+imageopenmp(ithread,izone,ix,iy,1:4)
				fluxZ(izone)=fluxZ(izone)+imageopenmp(ithread,izone,ix,iy,1)
			enddo
			enddo
		enddo
	enddo

	deallocate(imageopenmp)
	
	return
	end


	subroutine deallocatePaths
	use GlobalSetup
	IMPLICIT NONE
	integer i
	
	do i=1,nTracePaths
		deallocate(TracePaths(i)%inzone)
		deallocate(TracePaths(i)%C)
		deallocate(TracePaths(i)%v)
	enddo
	
2	deallocate(Pimage)
	deallocate(TracePaths)	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine MakeRotateStokes(phot,x,y,z,vx,vy,vz,sin2t,cos2t)
	use GlobalSetup
	IMPLICIT NONE
	type(photon) phot
	real*8 x,y,z,cost,sint,cos2t,sin2t,een,vx,vy,vz
	
	cost=phot%Sx*x+phot%Sy*y+phot%Sz*z
	sint=(x-
     &	(phot%Sx*(vy**2+vz**2)-vx*(vy*phot%Sy+vz*phot%Sz))*cost)/
     &	(vy*phot%Sz-vz*phot%Sy)
	
	cos2t=cost**2-sint**2
	sin2t=2d0*sint*cost

	een=sqrt(sin2t**2+cos2t**2)
	sin2t=sin2t/een
	cos2t=cos2t/een

	return
	end
