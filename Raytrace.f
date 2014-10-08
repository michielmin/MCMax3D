	subroutine Raytrace(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,i,ilam,ilam0
	real*8 dlmin
	logical simpleobs

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
		do ilam=1,ilam
			if(abs(lam(ilam)-(lam1+lam2)/2d0).lt.dlmin) then
				dlmin=abs(lam(ilam)-(lam1+lam2)/2d0)
				ilam0=ilam
			endif
		enddo
	endif

	call SetupPaths(iobs)
	do ilam=1,nlam
		if((lam(ilam).ge.MCobs(iobs)%lam1.and.lam(ilam).le.MCobs(iobs)%lam2)
     &			.or.ilam.eq.ilam0) then
			call TraceScattField(iobs,ilam)
			call FormalSolution(iobs,ilam)
		endif
	enddo

	call deallocatePaths
	
	return
	end
	
	subroutine TraceScattField(iobs,ilam)
	use GlobalSetup
	IMPLICIT NONE
	integer iobs,ilam,izone,iT,i,i1,i2,i3,iphot,istar
	real*8 GetKabs,Etot,Erandom,random,x,y,z,r
	logical emitfromstar
	type(Cell),pointer :: C
	type(Photon),allocatable :: phot(:)
	integer ispat,nspat
	integer,allocatable :: zspat(:),i1spat(:),i2spat(:),i3spat(:)
	real*8,allocatable :: Espat(:)
	integer maxnopenmp,iopenmp,omp_get_thread_num,omp_get_max_threads

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
			iT=(C%T+0.5d0)/dTBB
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
			C%KabsL=GetKabs(ilam,C)
			C%Elam=C%KabsL*BB(ilam,iT)*C%V
			Etot=Etot+C%Elam
			if(C%Elam.gt.0d0) then
				nspat=nspat+1
				Espat(nspat+1)=Espat(nspat)+C%Elam
				zspat(nspat)=izone
				i1spat(nspat)=i1
				i2spat(nspat)=i2
				i3spat(nspat)=i3
			endif
		enddo
		enddo
		enddo
	enddo

	call tellertje(1,100)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(iphot,iopenmp,Erandom,emitfromstar,ispat,i1,i2,i3,istar,izone,x,y,z,r)
!$OMP& SHARED(phot,iobs,Espat,Star,MCobs,nspat,i1spat,i2spat,i3spat,zspat,Etot,idum,nstars,ilam)
!$OMP DO
	do iphot=1,MCobs(iobs)%Nphot
!$OMP CRITICAL
		call tellertje(iphot+1,MCobs(iobs)%Nphot+2)
!$OMP END CRITICAL
		iopenmp=omp_get_thread_num()+1
		phot(iopenmp)%sI=Etot/real(MCobs(iobs)%Nphot)
		phot(iopenmp)%sQ=0d0
		phot(iopenmp)%sU=0d0
		phot(iopenmp)%sV=0d0
		phot(iopenmp)%ilam1=ilam
		phot(iopenmp)%pol=.false.
2		Erandom=Etot*random(idum)
		emitfromstar=.false.
		do istar=1,nstars
			Erandom=Erandom-Star(istar)%F(ilam)
			if(Erandom.lt.0d0) then
				emitfromstar=.true.
				goto 1
			endif
		enddo

		call hunt(Espat,nspat+1,Erandom,ispat)
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
		x=-phot(iopenmp)%vy
		y=phot(iopenmp)%vx
		z=0d0
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

	call tellertje(100,100)
	
	return
	end


	subroutine EmitPhotonStar(phot,istar)
	use GlobalSetup
	IMPLICIT NONE
	integer istar
	type(Photon) phot
	
	call randomdirection(phot%x,phot%y,phot%z)

	phot%x=Star(istar)%R*phot%x
	phot%y=Star(istar)%R*phot%y
	phot%z=Star(istar)%R*phot%z

	call randomdirection(phot%vx,phot%vy,phot%vz)

	if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
		phot%vx=-phot%vx
		phot%vy=-phot%vy
		phot%vz=-phot%vz
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
	integer izone,imin,iobs
	logical leave,inany
	real*8 minv,tau0,tau,GetKext,random,GetKabs,fstopmono,albedo,theta
	type(Travel) Trac(nzones)
	type(Cell),pointer :: C
	type(Photon) phot

	tau0=-log(random(idum))
	
	theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
	phot%iscat=180d0*theta/pi
	if(phot%iscat.lt.1) phot%iscat=1
	if(phot%iscat.gt.180) phot%iscat=180
	
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

	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone))
			end select
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

	if(leave) goto 3
	
	if((tau0-phot%Kext*minv).lt.0d0) then
		minv=tau0/phot%Kext
		phot%edgeNr=0
		call TravelPhotonX(phot,minv)
		call AddEtraceMono(phot,minv)
		call InteractMono(phot)
		theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
		phot%iscat=180d0*theta/pi
		if(phot%iscat.lt.1) phot%iscat=1
		if(phot%iscat.gt.180) phot%iscat=180
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
	if(inany) call AddEtraceMono(phot,minv)
	tau0=tau0-phot%Kext*minv

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

	goto 1
3	continue
	
	return
	end



	subroutine AddEtraceMono(phot,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer izone
	real*8 v,GetF11
	type(Cell),pointer :: C
	
	do izone=1,nzones
		if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			C%Escatt=C%Escatt+v*GetF11(phot%ilam1,phot%iscat,C)*phot%sI
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
	type(ZoneType),pointer :: Zo
	type(Photon),allocatable :: phot(:)
	integer maxnopenmp,iopenmp,omp_get_thread_num,omp_get_max_threads

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
	
	allocate(Pimage(nzones))
	
	do izone=1,nzones
		Zo => Zone(izone)
		Pimage(izone)%nr=3*(2*Zo%nR+2*Zo%nt+200)
		Pimage(izone)%np=Zo%np*2
		allocate(Pimage(izone)%R(Pimage(izone)%nr))
		allocate(Pimage(izone)%P(Pimage(izone)%nr,Pimage(izone)%np))
		ir=0
		scale(1)=Zo%xscale
		scale(2)=Zo%yscale
		scale(3)=Zo%zscale
		do ii=1,3
		do i=1,200
			ir=ir+1
			Pimage(izone)%R(ir)=scale(ii)*Zo%R(1)*real(i)/201d0
		enddo
		do j=1,Zo%nt
			ir=ir+1
			Pimage(izone)%R(ir)=abs(scale(ii)*(sqrt(Zo%R(1)*Zo%R(2))*sin((Zo%theta(j)+Zo%theta(j+1))/2d0)))
			ir=ir+1
			Pimage(izone)%R(ir)=abs(scale(ii)*(sqrt(Zo%R(1)*Zo%R(2))*sin(MCobs(iobs)%theta-Zo%theta0+(Zo%theta(j)+Zo%theta(j+1))/2d0)))
		enddo
		do i=1,Zo%nR
			ir=ir+1
			Pimage(izone)%R(ir)=abs(scale(ii)*sqrt(Zo%R(i)*Zo%R(i+1)))
			ir=ir+1
			Pimage(izone)%R(ir)=abs(scale(ii)*sqrt(Zo%R(i)*Zo%R(i+1))*sin(MCobs(iobs)%theta-Zo%theta0))
		enddo
		enddo
		Pimage(izone)%nr=ir
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
		x=Zo%x0
		y=Zo%y0
		z=Zo%z0
		call rotateZ(x,y,z,MCobs(iobs)%cosp,-MCobs(iobs)%sinp)
		call rotateY(x,y,z,MCobs(iobs)%cost,MCobs(iobs)%sint)
		Pimage(izone)%x=x
		Pimage(izone)%y=y
		call output("Paths for zone " // trim(int2string(izone,'(i4)')))
		call output("Number of rays: nr x nphi = " // trim(int2string(Pimage(izone)%nr,'(i4)')) // "x" 
     &				// trim(int2string(Pimage(izone)%np,'(i4)')) // " = "  
     &				// trim(int2string(Pimage(izone)%nr*Pimage(izone)%np,'(i8)')))
		call tellertje(1,100)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,ip,R1,R2,i,Rad,phi,iopenmp)
!$OMP& SHARED(Pimage,izone,MCobs,idum,iobs,phot,maxR,x,y,z)
!$OMP DO
		do ir=1,Pimage(izone)%nr
			iopenmp=omp_get_thread_num()+1
!$OMP CRITICAL
			call tellertje(ir+1,Pimage(izone)%nr+2)
!$OMP END CRITICAL
			do ip=1,Pimage(izone)%np
				if(ir.eq.1) then
					R1=Pimage(izone)%R(ir)
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
				else if(ir.eq.Pimage(izone)%nr) then
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=Pimage(izone)%R(ir)
				else
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
				endif
				Rad=R1+random(idum)*(R2-R1)
				phi=2d0*pi*(real(ip-1)+random(idum))/real(Pimage(izone)%np)

				phot(iopenmp)%x=Rad*cos(phi)+x
				phot(iopenmp)%y=Rad*sin(phi)+y
				phot(iopenmp)%z=0d0+z

				call rotateY(phot(iopenmp)%x,phot(iopenmp)%y,phot(iopenmp)%z,cos(MCobs(iobs)%theta),-sin(MCobs(iobs)%theta))
				call rotateZ(phot(iopenmp)%x,phot(iopenmp)%y,phot(iopenmp)%z,cos(MCobs(iobs)%phi),sin(MCobs(iobs)%phi))
				phot(iopenmp)%vx=0d0
				phot(iopenmp)%vy=0d0
				phot(iopenmp)%vz=-1d0
				call rotateY(phot(iopenmp)%vx,phot(iopenmp)%vy,phot(iopenmp)%vz,cos(MCobs(iobs)%theta),-sin(MCobs(iobs)%theta))
				call rotateZ(phot(iopenmp)%vx,phot(iopenmp)%vy,phot(iopenmp)%vz,cos(MCobs(iobs)%phi),sin(MCobs(iobs)%phi))

				call TravelPhotonX(phot(iopenmp),-maxR)

				call TranslatePhotonX(phot(iopenmp))
				call TranslatePhotonV(phot(iopenmp))
				
				call RaytracePath(phot(iopenmp),Pimage(izone)%P(ir,ip))
			enddo
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		call tellertje(100,100)
	enddo
		
	
	
	return
	end



	subroutine RaytraceFlux(P,flux,ilam,izone0)
	use GlobalSetup
	IMPLICIT NONE
	type(Path),target :: P
	type(Path),pointer :: Pnow
	real*8 flux
	integer ilam,izone0,istar
	type(Cell),pointer :: C


	integer izone,imin,iobs,iT
	logical leave,inany,hitstar0
	real*8 minv,tau0,tau,GetKext,random,GetKabs,theta,tau_e,fact,exptau,frac
	real*8 albedo,emis,scat,Kext,Kabs,KabsZ(nzones)

	Pnow => P

	tau0=0d0
	fact=1d0
	flux=0d0
	
1	continue

	Kext=0d0
	Kabs=0d0
	inany=.false.
	do izone=1,nzones
		KabsZ(izone)=0d0
		if(Pnow%inzone(izone)) then
			inany=.true.
			C => Pnow%C(izone)%C
			Kext=Kext+GetKext(ilam,C)
			KabsZ(izone)=GetKabs(ilam,C)
			Kabs=Kabs+KabsZ(izone)
		endif
	enddo
	albedo=(Kext-Kabs)/Kext

	if(inany) then
		tau_e=Kext*Pnow%v

		emis=0d0
		scat=0d0
		if(Pnow%inzone(izone0)) then
			C => Pnow%C(izone0)%C
			iT=(C%T+0.5d0)/dTBB
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
			emis=emis+BB(ilam,iT)*KabsZ(izone0)
			scat=scat+C%Escatt/C%V
		endif
		emis=emis*(1d0-albedo)/Kabs
		scat=scat*albedo/Kext

		if(tau_e.lt.1d-6) then
			flux=flux+(scat+emis)*tau_e*fact
			fact=fact*(1d0-tau_e)
		else
			exptau=exp(-tau_e)
			frac=(1d0-exptau)
			flux=flux+(scat+emis)*frac*fact
			fact=fact*exptau
		endif
		tau0=tau0+tau_e	
	endif

	if(.not.Pnow%last) then
		Pnow => Pnow%next
		goto 1
	endif	

	return
	end



	subroutine RaytraceFluxOld(phot,flux,ilam,izone0)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 flux
	integer ilam,izone0,istar


	integer izone,imin,iobs,iT
	logical leave,inany,hitstar0
	real*8 minv,tau0,tau,GetKext,random,GetKabs,theta,tau_e,fact,exptau,frac
	real*8 albedo,emis,scat
	type(Travel) Trac(nzones),TracStar(nstars)
	type(Cell),pointer :: C

	tau0=0d0
	phot%inzone=.false.
	phot%edgeNr=0
	fact=1d0
	flux=0d0
	
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
	albedo=(phot%Kext-phot%Kabs)/phot%Kext

	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone))
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

	call TravelPhotonX(phot,minv)

	if(inany) then
		tau_e=phot%Kext*minv

		emis=0d0
		scat=0d0
		if(phot%inzone(izone0)) then
			C => Zone(izone0)%C(phot%i1(izone0),phot%i2(izone0),phot%i3(izone0))
			iT=(C%T+0.5d0)/dTBB
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
			emis=emis+BB(ilam,iT)*phot%KabsZ(izone0)
			scat=scat+C%Escatt/C%V
		endif
		emis=emis*(1d0-albedo)/phot%Kabs
		scat=scat*albedo/phot%Kext

		if(tau_e.lt.1d-6) then
			flux=flux+(scat+emis)*tau_e*fact
			fact=fact*(1d0-tau_e)
		else
			exptau=exp(-tau_e)
			frac=(1d0-exptau)
			flux=flux+(scat+emis)*frac*fact
			fact=fact*exptau
		endif
		tau0=tau0+tau_e	
	endif

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
		else
			phot%edgeNr(izone)=0
		endif
	enddo

	goto 1
3	continue
	

	return
	end



	subroutine RaytracePath(phot,P)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	type(Path),target :: P
	type(Path),pointer :: Pnow
	integer istar


	integer izone,imin,iobs
	logical leave,inany,hitstar0
	real*8 minv,theta
	type(Travel) Trac(nzones),TracStar(nstars)
	type(Cell),pointer :: C

	phot%inzone=.false.
	phot%edgeNr=0

	Pnow => P
	
1	continue

	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone))
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

	allocate(Pnow%inzone(nzones))
	allocate(Pnow%C(nzones))
	allocate(Pnow%next)
	do izone=1,nzones
		Pnow%inzone(izone)=phot%inzone(izone)
		if(Pnow%inzone(izone)) then
			Pnow%C(izone)%C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
		endif
	enddo
	Pnow%v=minv
	Pnow%last=.false.

	if(leave) goto 3
		
	call TravelPhotonX(phot,minv)

	if(inany) then
		
	endif

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
		else
			phot%edgeNr(izone)=0
		endif
	enddo

	Pnow => Pnow%next

	goto 1
3	continue
	
	Pnow%last=.true.

	return
	end




	subroutine FormalSolution(iobs,ilam)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,ilam,izone,ir,ip,i,j,ix,iy,ii
	real*8 random,phi,x,y,flux,A,R1,R2,Rad,scale(3)
	type(ZoneType),pointer :: Zo

	MCobs(iobs)%image(:,:,ilam)=0d0
	
	do izone=1,nzones
		Zo => Zone(izone)
		call output("Formal solution zone " // trim(int2string(izone,'(i4)')))

		call tellertje(1,100)
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,ip,R1,R2,flux,A,i,Rad,phi,x,y,ix,iy)
!$OMP& SHARED(Pimage,izone,ilam,MCobs,idum,iobs)
!$OMP DO
		do ir=1,Pimage(izone)%nr
!$OMP CRITICAL
			call tellertje(ir+1,Pimage(izone)%nr+2)
!$OMP END CRITICAL
			do ip=1,Pimage(izone)%np
				if(ir.eq.1) then
					R1=Pimage(izone)%R(ir)
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
				else if(ir.eq.Pimage(izone)%nr) then
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=Pimage(izone)%R(ir)
				else
					R1=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))
					R2=sqrt(Pimage(izone)%R(ir)*Pimage(izone)%R(ir+1))
				endif
				
				call RaytraceFlux(Pimage(izone)%P(ir,ip),flux,ilam,izone)

				if(ir.eq.1) then
					A=pi*(Pimage(izone)%R(ir+1)*Pimage(izone)%R(ir)-Pimage(izone)%R(ir)**2)/real(Pimage(izone)%np)
				else if(ir.eq.Pimage(izone)%nr) then
					A=pi*(Pimage(izone)%R(ir)**2-Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))/real(Pimage(izone)%np)
				else
					A=pi*(Pimage(izone)%R(ir+1)*Pimage(izone)%R(ir)-Pimage(izone)%R(ir)*Pimage(izone)%R(ir-1))/real(Pimage(izone)%np)
				endif


				do i=1,1000
					Rad=R1+random(idum)*(R2-R1)
					phi=2d0*pi*(real(ip-1)+random(idum))/real(Pimage(izone)%np)
					x=Rad*cos(phi)+Pimage(izone)%x
					y=Rad*sin(phi)+Pimage(izone)%y

					ix=real(MCobs(iobs)%npix)*(y+MCobs(iobs)%maxR)/(2d0*MCobs(iobs)%maxR)
					iy=MCobs(iobs)%npix-real(MCobs(iobs)%npix)*(x+MCobs(iobs)%maxR)/(2d0*MCobs(iobs)%maxR)

					if(ix.lt.MCobs(iobs)%npix.and.iy.lt.MCobs(iobs)%npix.and.ix.gt.0.and.iy.gt.0) then
						MCobs(iobs)%image(ix,iy,ilam)=MCobs(iobs)%image(ix,iy,ilam)+1d23*flux*A/1000d0/(32d0*pi)
					endif
				enddo
			enddo
		enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
		call tellertje(100,100)
	enddo
	
	return
	end


	subroutine deallocatePaths
	use GlobalSetup
	IMPLICIT NONE
	integer izone,ir,ip
	type(Path),pointer :: Pnow
	
	do izone=1,nzones
		do ir=1,Pimage(izone)%nr
			do ip=1,Pimage(izone)%np
				deallocate(Pimage(izone)%P(ir,ip)%inzone)
				deallocate(Pimage(izone)%P(ir,ip)%C)
				Pnow => Pimage(izone)%P(ir,ip)%next
1				continue
				if(Pnow%last) goto 2
				deallocate(Pnow%inzone)
				deallocate(Pnow%C)
				if(.not.Pnow%last) then
					Pnow => Pnow%next
					goto 1
				endif
			enddo
		enddo
	enddo
	
2	deallocate(Pimage)
	
	return
	end
