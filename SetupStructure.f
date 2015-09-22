	subroutine SetupStructure
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j,nf
	real*8 T,Planck
	real*8 x,y,z,inp,r
	
	maxR=0d0
	
	call SetupLam()
	
c setup the blackbodies to use
	allocate(BB(nlam,0:nBB))
	do i=1,nlam
		BB(i,0)=0d0
		do j=1,nBB
			T=dTBB*real(j)
			BB(i,j)=Planck(T,lam(i))
		enddo
	enddo

c setup the observation direction
	do i=1,nMCobs
		MCobs(i)%theta=MCobs(i)%theta*pi/180d0
		MCobs(i)%cost=cos(MCobs(i)%theta)
		MCobs(i)%sint=sin(MCobs(i)%theta)
		MCobs(i)%phi=MCobs(i)%phi*pi/180d0
		MCobs(i)%cosp=cos(MCobs(i)%phi)
		MCobs(i)%sinp=sin(MCobs(i)%phi)
		MCobs(i)%opening=cos(MCobs(i)%opening*pi/180d0)
		MCobs(i)%x=MCobs(i)%cosp*MCobs(i)%sint
		MCobs(i)%y=MCobs(i)%sinp*MCobs(i)%sint
		MCobs(i)%z=MCobs(i)%cost
		MCobs(i)%xup=-MCobs(i)%cosp*MCobs(i)%cost
		MCobs(i)%yup=-MCobs(i)%sinp*MCobs(i)%cost
		MCobs(i)%zup=MCobs(i)%sint
		
		allocate(MCobs(i)%image(MCobs(i)%npix,MCobs(i)%npix,nlam))
		allocate(MCobs(i)%spec(nlam))
		MCobs(i)%f=2d0*pi*(1d0-MCobs(i)%opening)
	enddo

c spiral waves
	do i=1,nSpirals
		Spiral(i)%r=Spiral(i)%r*AU
		Spiral(i)%w=Spiral(i)%w*AU
		Spiral(i)%phi=Spiral(i)%phi*pi/180d0
	enddo
	
	call output("==================================================================")

	do i=1,npart
		call SetupPart(i)
	enddo

	call output("==================================================================")

	do i=1,nstars
		call SetupStar(i)
	enddo
	
	call output("==================================================================")

	do i=1,nzones
		call SetupZone(i)
	enddo
	
	do i=1,nMCobs
		if(MCobs(i)%fov.gt.0d0) then
			MCobs(i)%maxR=MCobs(i)%fov*distance*AU/2d0
		else
			if(MCobs(i)%maxR.lt.0d0) then
				MCobs(i)%maxR=maxR*1.5d0
			else
				MCobs(i)%maxR=MCobs(i)%maxR*AU
			endif
			MCobs(i)%fov=2d0*MCobs(i)%maxR/(distance*AU)
		endif
	enddo

	distance=distance*parsec
		
	call SetupBeaming()

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine SetupLam
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 i,nl,j
	
	allocate(lam(nlam))
	allocate(nu(nlam))
	allocate(dnu(nlam))
	
	if(nzlam.le.0) then
		do i=1,nlam
			lam(i)=10d0**(log10(lam1)+(log10(lam2)-log10(lam1))*real(i-1)/real(nlam-1))
		enddo
	else
c	does not seem to work!!! Have to fix this!!!
		call output("NZLAM OPTION NOT ALWAYS WORKING PROPERLY YET!!!")
		nl=(nlam-nzlam)
		do i=1,nl
			lam(i)=10d0**(log10(lam1)+(log10(lam2)-log10(lam1))*real(i-1)/real(nl-1))
		enddo
		j=nl
		do i=1,nzlam
			lam(i+j)=10d0**(log10(zlam1)+(log10(zlam2)-log10(zlam1))*real(i-1)/real(nzlam-1))
		enddo
		call sort(lam,nlam)
	endif

	do i=1,nlam
		nu(i)=1d4*clight/lam(i)
	enddo

	do i=2,nlam-1
		dnu(i)=(nu(i-1)-nu(i+1))/2d0
	enddo
	dnu(nlam)=(nu(nlam-1)-nu(nlam))
	dnu(1)=(nu(1)-nu(2))

	open(unit=20,file=trim(outputdir) // 'lamgrid.dat')
	do i=1,nlam
		write(20,*) lam(i)
	enddo
	close(unit=20)

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine SetupPart(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,iT,is,iBB,ilam,j
	real*8 phi,thet,tot,tot2,fact
	real*8,allocatable :: spec(:)
	
	allocate(Part(ii)%rv(Part(ii)%nsize))
	allocate(Part(ii)%rho(Part(ii)%nT))
	allocate(Part(ii)%Kabs(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Ksca(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Kext(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%F(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Kp(Part(ii)%nsize,Part(ii)%nT,0:nBB))
	allocate(spec(nlam))
	
	select case(Part(ii)%ptype)
		case("COMPUTE")
			do is=1,Part(ii)%nsize
				do iT=1,Part(ii)%nT
					call ComputePart(Part(ii),ii,is,iT)
				enddo
			enddo
c		case("PARTFILE")
c			call ReadParticle(Part(ii),ii)
		case default
			call output("I did not understand what I was trying to do. Sorry!")
	end select
	
	do iBB=1,nBB
		do is=1,Part(ii)%nsize
			do iT=1,Part(ii)%nT
				call integrate(BB(1:nlam,iBB)*Part(ii)%Kabs(is,iT,1:nlam),Part(ii)%Kp(is,iT,iBB))
			enddo
		enddo
	enddo
	do is=1,Part(ii)%nsize
		do iT=1,Part(ii)%nT
			Part(ii)%Kp(is,iT,0)=0d0
		enddo
	enddo

	if(nspike.gt.0.and.nspike.le.180) call output("making the first " // trim(int2string(nspike,'(i3)')) // " degrees isotropic")
	do ilam=1,nlam
		do is=1,Part(ii)%nsize
			do iT=1,Part(ii)%nT
				tot=0d0
				tot2=0d0
				do j=1,180
					tot=tot+Part(ii)%F(is,iT,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
					tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
				enddo
				do j=1,180
					Part(ii)%F(is,iT,ilam)%F11(j)=tot2*Part(ii)%F(is,iT,ilam)%F11(j)/tot
					Part(ii)%F(is,iT,ilam)%F12(j)=tot2*Part(ii)%F(is,iT,ilam)%F12(j)/tot
					Part(ii)%F(is,iT,ilam)%F22(j)=tot2*Part(ii)%F(is,iT,ilam)%F22(j)/tot
					Part(ii)%F(is,iT,ilam)%F33(j)=tot2*Part(ii)%F(is,iT,ilam)%F33(j)/tot
					Part(ii)%F(is,iT,ilam)%F34(j)=tot2*Part(ii)%F(is,iT,ilam)%F34(j)/tot
					Part(ii)%F(is,iT,ilam)%F44(j)=tot2*Part(ii)%F(is,iT,ilam)%F44(j)/tot
				enddo

				if(nspike.gt.0.and.nspike.le.180) then
c the nspike parameter removes the n degree spike in the forward direction.
					do j=1,nspike
						fact=Part(ii)%F(is,iT,ilam)%F11(nspike+1)/Part(ii)%F(is,iT,ilam)%F11(j)
						Part(ii)%F(is,iT,ilam)%F12(j)=Part(ii)%F(is,iT,ilam)%F12(j)*fact
						Part(ii)%F(is,iT,ilam)%F22(j)=Part(ii)%F(is,iT,ilam)%F22(j)*fact
						Part(ii)%F(is,iT,ilam)%F33(j)=Part(ii)%F(is,iT,ilam)%F33(j)*fact
						Part(ii)%F(is,iT,ilam)%F34(j)=Part(ii)%F(is,iT,ilam)%F34(j)*fact
						Part(ii)%F(is,iT,ilam)%F44(j)=Part(ii)%F(is,iT,ilam)%F44(j)*fact
						Part(ii)%F(is,iT,ilam)%F11(j)=Part(ii)%F(is,iT,ilam)%F11(nspike+1)
					enddo

					tot=0d0
					tot2=0d0
					do j=1,180
						tot=tot+Part(ii)%F(is,iT,ilam)%F11(j)*sin(pi*(real(j)-0.5)/180d0)
						tot2=tot2+sin(pi*(real(j)-0.5)/180d0)
					enddo
					Part(ii)%Ksca(is,iT,ilam)=Part(ii)%Ksca(is,iT,ilam)*tot/tot2
					Part(ii)%Kext(is,iT,ilam)=Part(ii)%Kabs(is,iT,ilam)+Part(ii)%Ksca(is,iT,ilam)
					do j=1,180
						Part(ii)%F(is,iT,ilam)%F11(j)=tot2*Part(ii)%F(is,iT,ilam)%F11(j)/tot
						Part(ii)%F(is,iT,ilam)%F12(j)=tot2*Part(ii)%F(is,iT,ilam)%F12(j)/tot
						Part(ii)%F(is,iT,ilam)%F22(j)=tot2*Part(ii)%F(is,iT,ilam)%F22(j)/tot
						Part(ii)%F(is,iT,ilam)%F33(j)=tot2*Part(ii)%F(is,iT,ilam)%F33(j)/tot
						Part(ii)%F(is,iT,ilam)%F34(j)=tot2*Part(ii)%F(is,iT,ilam)%F34(j)/tot
						Part(ii)%F(is,iT,ilam)%F44(j)=tot2*Part(ii)%F(is,iT,ilam)%F44(j)/tot
					enddo
				endif
			enddo
		enddo
	enddo

	do ilam=1,nlam
		do is=1,Part(ii)%nsize
			do iT=1,Part(ii)%nT
				Part(ii)%F(is,iT,ilam)%IF11=0d0
				Part(ii)%F(is,iT,ilam)%IF12=0d0
				do j=1,180
					thet=pi*(real(j)-0.5d0)/180d0
					Part(ii)%F(is,iT,ilam)%IF11=Part(ii)%F(is,iT,ilam)%IF11+pi*sin(thet)
     &			*Part(ii)%F(is,iT,ilam)%F11(j)/180d0
					Part(ii)%F(is,iT,ilam)%IF12=Part(ii)%F(is,iT,ilam)%IF12+pi*sin(thet)
     &			*Part(ii)%F(is,iT,ilam)%F12(j)/180d0
				enddo
			enddo
		enddo
	enddo
	do j=0,360
		phi=pi*real(j-1)/179.5d0
		cos2phi(j)=cos(2d0*phi)
		sin2phi(j)=sin(2d0*phi)
	enddo

	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine SetupStar(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i
	real*8 Luminosity,tot,Planck
	
	call output("Setting up star: " // trim(int2string(ii,'(i0.4)')))
	Star(ii)%L=Star(ii)%L*Lsun
	allocate(Star(ii)%F(nlam))

	select case(Star(ii)%startype)
		case("PLANCK")
			do i=1,nlam
				Star(ii)%F(i)=Planck(Star(ii)%T,lam(i))
			enddo
			call integrate(Star(ii)%F,tot)
			Star(ii)%F=Star(ii)%L*Star(ii)%F/tot
		case("KURUCZ")
			call ReadKurucz(Star(ii)%T,Star(ii)%logg,lam,Star(ii)%F,nlam)
			call integrate(Star(ii)%F,tot)
			Star(ii)%F=Star(ii)%L*Star(ii)%F/tot
		case("FILE")
			call readstar(Star(ii)%file,lam,Star(ii)%F,nlam)
			call integrate(Star(ii)%F,tot)
			Star(ii)%F=Star(ii)%L*Star(ii)%F/tot
		case default
			call output("Nope, such stars do not exist in our universe...")
	end select

	Star(ii)%R=Rsun*sqrt(Star(ii)%L/Luminosity(Star(ii)%T,Rsun))

	call output("Stellar temperature: " // trim(dbl2string(Star(ii)%T,'(f14.2)')) // " K")
	call output("Stellar luminosity:  " // trim(dbl2string(Star(ii)%L/Lsun,'(f14.2)')) // " Lsun")
	call output("Stellar radius:      " // trim(dbl2string(Star(ii)%R/Rsun,'(f14.2)')) // " Rsun")

	open(unit=20,file=trim(outputdir) // 'star' // trim(int2string(ii,'(i0.4)')) // '.dat')
	do i=1,nlam
		write(20,*) lam(i),Star(ii)%F(i)
	enddo
	close(unit=20)

	Star(ii)%x=Star(ii)%x*AU
	Star(ii)%y=Star(ii)%y*AU
	Star(ii)%z=Star(ii)%z*AU

	if((Star(ii)%x+Star(ii)%R).gt.maxR) maxR=(abs(Star(ii)%x)+Star(ii)%R)
	if((Star(ii)%y+Star(ii)%R).gt.maxR) maxR=(abs(Star(ii)%y)+Star(ii)%R)
	if((Star(ii)%z+Star(ii)%R).gt.maxR) maxR=(abs(Star(ii)%z)+Star(ii)%R)
	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	

	subroutine SetupZone(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i1,i2,i3,i
	real*8 random
	character*500 filename
	
	call output("Setting up zone nr.: "// trim(int2string(ii,'(i4)')))

	Zone(ii)%theta0=Zone(ii)%theta0*pi/180d0
	Zone(ii)%cost0=cos(Zone(ii)%theta0)
	Zone(ii)%sint0=sin(Zone(ii)%theta0)
	Zone(ii)%phi0=Zone(ii)%phi0*pi/180d0
	Zone(ii)%cosp0=cos(Zone(ii)%phi0)
	Zone(ii)%sinp0=sin(Zone(ii)%phi0)
	
	call output("allocating memory")
	select case(Zone(ii)%shape)
		case("SPH","CYL")
			allocate(Zone(ii)%C(Zone(ii)%nr,Zone(ii)%nt,Zone(ii)%np))
			allocate(Zone(ii)%R(Zone(ii)%nr+1))
			allocate(Zone(ii)%theta(Zone(ii)%nt+1))
			allocate(Zone(ii)%phi(Zone(ii)%np+1))
			allocate(Zone(ii)%R2(Zone(ii)%nr+1))
			allocate(Zone(ii)%cost2(Zone(ii)%nt+1))
			allocate(Zone(ii)%tanx(Zone(ii)%np+1))
			allocate(Zone(ii)%tany(Zone(ii)%np+1))
			Zone(ii)%n1=Zone(ii)%nr
			Zone(ii)%n2=Zone(ii)%nt
			Zone(ii)%n3=Zone(ii)%np
		case("CAR")
			allocate(Zone(ii)%C(Zone(ii)%nx,Zone(ii)%ny,Zone(ii)%nz))
			allocate(Zone(ii)%x(Zone(ii)%nx+1))
			allocate(Zone(ii)%y(Zone(ii)%ny+1))
			allocate(Zone(ii)%z(Zone(ii)%nz+1))
			Zone(ii)%n1=Zone(ii)%nx
			Zone(ii)%n2=Zone(ii)%ny
			Zone(ii)%n3=Zone(ii)%nz
		case default
			call output("Interesting shape option ("// Zone(ii)%shape //"). But I don't get it...")
			stop
	end select
	do i1=1,Zone(ii)%n1
		call tellertje(i1,Zone(ii)%nr)
		do i2=1,Zone(ii)%n2
			do i3=1,Zone(ii)%n3
				allocate(Zone(ii)%C(i1,i2,i3)%densP(npart,maxns,maxnT))
			enddo
		enddo
	enddo

	call ScaleZoneInput(ii)

c avoid zones with exactly the same inner or outer radii
1	continue
	do i=1,ii-1
		if(Zone(ii)%Rin.eq.Zone(i)%Rin.or.Zone(ii)%Rin.eq.Zone(i)%Rout) then
			Zone(ii)%Rin=Zone(ii)%Rin*(1d0+1d-6*random(idum))
			call output("Slightly adjusting Rin")
			goto 1
		endif
		if(Zone(ii)%Rout.eq.Zone(i)%Rin.or.Zone(ii)%Rout.eq.Zone(i)%Rout) then
			Zone(ii)%Rout=Zone(ii)%Rout*(1d0-1d-6*random(idum))
			call output("Slightly adjusting Rout")
			goto 1
		endif
	enddo
	
	if(Nphot.eq.0) then
		filename=trim(outputdir) // "Zone" // trim(int2string(ii,'(i0.4)')) // ".fits.gz"
		call output("Reading file: "// trim(filename))
		call readstruct_fits(filename,ZoneStructOutput,nZoneStructOutput,ii)
		Zone(ii)%imidplane=Zone(ii)%nt/2+1
		do i=1,Zone(ii)%nt+1
			Zone(ii)%cost2(i)=cos(Zone(ii)%theta(i))**2
		enddo
		if((Zone(ii)%x0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%x0)+Zone(ii)%Rout)
		if((Zone(ii)%y0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%y0)+Zone(ii)%Rout)
		if((Zone(ii)%z0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%z0)+Zone(ii)%Rout)
		do i=1,Zone(ii)%nr+1
			Zone(ii)%R2(i)=Zone(ii)%R(i)**2
		enddo
		do i=1,Zone(ii)%np+1
			if(Zone(ii)%phi(i).lt.(pi/4d0)) then
				Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
				Zone(ii)%tany(i)=1d0
			else if(Zone(ii)%phi(i).lt.(3d0*pi/4d0)) then
				Zone(ii)%tanx(i)=-1d0
				Zone(ii)%tany(i)=-cos(Zone(ii)%phi(i))/sin(Zone(ii)%phi(i))
			else if(Zone(ii)%phi(i).lt.(5d0*pi/4d0)) then
				Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
				Zone(ii)%tany(i)=1d0
			else if(Zone(ii)%phi(i).lt.(7d0*pi/4d0)) then
				Zone(ii)%tanx(i)=-1d0
				Zone(ii)%tany(i)=-cos(Zone(ii)%phi(i))/sin(Zone(ii)%phi(i))
			else
				Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
				Zone(ii)%tany(i)=1d0
			endif
		enddo
	else
		call output("setup density structure")
		select case(Zone(ii)%denstype)
			case("DISK","TAUFILE")
				call SetupDisk(ii)
			case("SHELL")
				call SetupShell(ii)
			case("SPHERE")
				call SetupSphere(ii)
			case default
				call SetupSpecialZone(ii)
		end select	
	endif

	return
	end
	
	subroutine ScaleZoneInput(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii
	
	select case(Zone(ii)%sscaletype)
		case("AU")
			Zone(ii)%sscale=AU
		case("CM","cm")
			Zone(ii)%sscale=1d0
		case("RJ","Rj","rj","RJUP","RJup","Rjup","rjup")
			Zone(ii)%sscale=Rjup
		case("RE","Re","re","REARTH","REarth","Rearth","rearth")
			Zone(ii)%sscale=Rearth
		case default
			call output("Unknown size scaling type")
	end select
	select case(Zone(ii)%mscaletype)
		case("gram","gr","GRAM","GR")
			Zone(ii)%mscale=1d0
		case("MSUN","Msun","MSun")
			Zone(ii)%mscale=Msun
		case("MJ","Mj","mj","MJUP","MJup","Mjup","mjup")
			Zone(ii)%mscale=Mjup
		case("ME","Me","me","MEARTH","MEarth","Mearth","mearth")
			Zone(ii)%mscale=Mearth
		case("MZODI","Mzodi","mzodi","MZodi")
			Zone(ii)%mscale=MZodi
		case default
			call output("Unknown size scaling type")
	end select

	Zone(ii)%x0=Zone(ii)%x0*AU		!Zone(ii)%sscale
	Zone(ii)%y0=Zone(ii)%y0*AU		!Zone(ii)%sscale
	Zone(ii)%z0=Zone(ii)%z0*AU		!Zone(ii)%sscale
	Zone(ii)%rvortex=Zone(ii)%rvortex*AU	!Zone(ii)%sscale

	Zone(ii)%Rin=Zone(ii)%Rin*Zone(ii)%sscale
	Zone(ii)%Rout=Zone(ii)%Rout*Zone(ii)%sscale
	Zone(ii)%Rexp=Zone(ii)%Rexp*Zone(ii)%sscale
	Zone(ii)%Rsh=Zone(ii)%Rsh*Zone(ii)%sscale
	Zone(ii)%sh=Zone(ii)%sh*Zone(ii)%sscale

	Zone(ii)%dx=Zone(ii)%dx*Zone(ii)%sscale
	Zone(ii)%dy=Zone(ii)%dy*Zone(ii)%sscale
	Zone(ii)%dz=Zone(ii)%dz*Zone(ii)%sscale

	Zone(ii)%Mdust=Zone(ii)%Mdust*Zone(ii)%mscale
	

	return
	end
	

	subroutine SetupDisk(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,ir,it,ip,jj,njj,i,ips,ipt
	real*8 r,z,hr,f1,f2,theta,Mtot,w(npart,maxns),ha,f2a,delta,phi,phi0
	real*8 RoundOff, roundfac
	real*8,allocatable :: Aspiral(:,:,:),H(:,:),SD(:,:),A(:,:),alpha(:,:),Avortex(:,:,:,:)
	delta=2d0

	if(Zone(ii)%shape.eq.'CAR') then
		call output("A disk on a rectangular grid... Let's don't and say we did.")
		stop
	endif

	call SetupRadGrid(ii)
	call SetupThetaGrid(ii)
	call SetupPhiGrid(ii)
	call SetupVolume(ii)

			
	do i=1,npart
		if(Part(i)%nsize.gt.1) then
			call SetSizeDis(w(i,1:Part(i)%nsize),i,ii)
		else
			w(i,1)=1d0
		endif
	enddo

	allocate(Aspiral(nSpirals,Zone(ii)%nr,Zone(ii)%np))
	allocate(Avortex(Zone(ii)%nr,Zone(ii)%np,npart,0:maxns))
	allocate(A(Zone(ii)%nr,Zone(ii)%np))
	allocate(H(Zone(ii)%nr,Zone(ii)%np))
	allocate(SD(Zone(ii)%nr,Zone(ii)%np))
	allocate(alpha(Zone(ii)%nr,Zone(ii)%np))
	
	Mtot=0d0
	Aspiral=0d0
	do i=1,nSpirals
		call ComputeAspiral(ii,i,Aspiral,Zone(ii)%nr,Zone(ii)%np)
	enddo
	
	
	do ir=1,Zone(ii)%nr
		r=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
		if(Zone(ii)%roundtype.ne.'NONE') then
			roundfac=RoundOff(r/Zone(ii)%sscale,Zone(ii)%roundradius,Zone(ii)%roundtype,Zone(ii)%roundindex,Zone(ii)%rounddepth)
		else
			roundfac= 1d0
		endif
		
		do ip=1,Zone(ii)%np
			H(ir,ip)=Zone(ii)%sh*(r/Zone(ii)%Rsh)**Zone(ii)%shpow
			SD(ir,ip)=r**(-Zone(ii)%denspow)*exp(-(r/Zone(ii)%Rexp)**(Zone(ii)%gamma_exp)) * roundfac
			A(ir,ip)=pi*(Zone(ii)%R(ir+1)**2-Zone(ii)%R(ir)**2)/real(Zone(ii)%np)
			alpha(ir,ip)=Zone(ii)%alpha
			Mtot=Mtot+A(ir,ip)*SD(ir,ip)
		enddo
	enddo
	SD=SD*Zone(ii)%Mdust/Mtot
	if(Zone(ii)%denstype.eq.'TAUFILE') then
		call readtaufile(ii,SD,Zone(ii)%nr,Zone(ii)%np)
	endif

	call ComputeAvortex(ii,Avortex,SD,H)

	open(unit=20,file=trim(outputdir) // "surfacedens" // trim(int2string(ii,'(i0.4)')) // ".dat")
	do ir=1,Zone(ii)%nr
		write(20,*) sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))/AU,SD(ir,1)
	enddo
	close(unit=20)

	do ir=1,Zone(ii)%nr
		do ip=1,Zone(ii)%np
			do i=1,nSpirals
				H(ir,ip)=H(ir,ip)*(1d0+Aspiral(i,ir,ip)*Spiral(i)%Aheight)
				SD(ir,ip)=SD(ir,ip)*(1d0+Aspiral(i,ir,ip)*Spiral(i)%Adens)
				alpha(ir,ip)=alpha(ir,ip)*(1d0+Aspiral(i,ir,ip)*Spiral(i)%Aalpha)
			enddo
		enddo
	enddo
	
	Mtot=0d0
	do ir=1,Zone(ii)%nr
		call tellertje(ir,Zone(ii)%nr)
		do it=1,Zone(ii)%nt
			njj=10
			Zone(ii)%C(ir,it,:)%gasdens=0d0
			do ip=1,Zone(ii)%np
				do i=1,npart
					do ips=1,Part(i)%nsize
						Zone(ii)%C(ir,it,ip)%densP(i,ips,:)=0d0
					enddo
				enddo
			enddo
			do jj=1,njj
				theta=Zone(ii)%theta(it)+(Zone(ii)%theta(it+1)-Zone(ii)%theta(it))*real(jj)/real(njj+1)
				r=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
				z=abs(sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))*cos(theta))
				do ip=1,Zone(ii)%np
					hr=H(ir,ip)
					f1=SD(ir,ip)
					f2=exp(-(z/hr)**2)
					Zone(ii)%C(ir,it,ip)%gasdens=Zone(ii)%C(ir,it,ip)%gasdens+f1*f2/hr/real(njj)
				enddo
			enddo
			do ip=1,Zone(ii)%np
				Zone(ii)%C(ir,it,ip)%gasdens=Zone(ii)%C(ir,it,ip)%gasdens*Avortex(ir,ip,1,0)
				Mtot=Mtot+Zone(ii)%C(ir,it,ip)%gasdens*Zone(ii)%C(ir,it,ip)%V
			enddo
		enddo
	enddo
	Mtot=Zone(ii)%Mdust/Mtot
	do ir=1,Zone(ii)%nr
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				Zone(ii)%C(ir,it,ip)%gasdens=Zone(ii)%gas2dust*Zone(ii)%C(ir,it,ip)%gasdens*Mtot
			enddo
		enddo
	enddo
	call tellertje(1,100)
!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ir,it,jj,njj,theta,r,z,hr,f1,f2a,ha)
!$OMP& SHARED(Zone,npart,Mtot,w,delta,Part,ii,SD,H,alpha,Avortex)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do ir=1,Zone(ii)%nr
!$OMP CRITICAL
		call tellertje(ir+1,Zone(ii)%nr+2)
!$OMP END CRITICAL
		do it=1,Zone(ii)%nt
			njj=10
			do jj=1,njj
				theta=Zone(ii)%theta(it)+(Zone(ii)%theta(it+1)-Zone(ii)%theta(it))*real(jj)/real(njj+1)
				r=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
				z=abs(sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))*cos(theta))
				do ip=1,Zone(ii)%np
					hr=H(ir,ip)
					f1=SD(ir,ip)
					do i=1,npart
						do ips=1,Part(i)%nsize
							ha=(1d0+delta)**(-0.25)*
     &			sqrt(alpha(ir,ip)*
     &			Zone(ii)%gas2dust*Mtot*f1
     &			/(Part(i)%rv(ips)*Part(i)%rho(1)))
							ha=ha*hr/sqrt(1d0+ha**2)
							f2a=exp(-(z/ha)**2)
							Zone(ii)%C(ir,it,ip)%densP(i,ips,1)=Zone(ii)%C(ir,it,ip)%densP(i,ips,1)+
     &		Mtot*w(i,ips)*Zone(ii)%abun(i)*f1*f2a/ha/real(njj)
						enddo
					enddo
				enddo
			enddo
			do ip=1,Zone(ii)%np
				do i=1,npart
					do ips=1,Part(i)%nsize
						Zone(ii)%C(ir,it,ip)%densP(i,ips,1)=Zone(ii)%C(ir,it,ip)%densP(i,ips,1)*Avortex(ir,ip,i,ips)
					enddo
				enddo
			enddo
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(100,100)

	Mtot=0d0
	do ir=1,Zone(ii)%nr
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				Zone(ii)%C(ir,it,ip)%dens=0d0
				do i=1,npart
					do ips=1,Part(i)%nsize
						do ipt=1,Part(i)%nT
							if(Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt).lt.1d-50) Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)=1d-50
						enddo
						Zone(ii)%C(ir,it,ip)%dens=Zone(ii)%C(ir,it,ip)%dens+Zone(ii)%C(ir,it,ip)%densP(i,ips,1)
					enddo
				enddo
				if(Zone(ii)%C(ir,it,ip)%dens.lt.1d-40) Zone(ii)%C(ir,it,ip)%dens=1d-40
			enddo
			do ip=1,Zone(ii)%np
				Mtot=Mtot+Zone(ii)%C(ir,it,ip)%dens*Zone(ii)%C(ir,it,ip)%V
			enddo
		enddo
	enddo
	Mtot=Zone(ii)%Mdust/Mtot
	do ir=1,Zone(ii)%nr
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				do i=1,npart
					do ips=1,Part(i)%nsize
						do ipt=1,Part(i)%nT
							Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)=Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)*Mtot
						enddo
					enddo
				enddo
				Zone(ii)%C(ir,it,ip)%dens=Zone(ii)%C(ir,it,ip)%dens*Mtot
			enddo
		enddo
	enddo
	if(Zone(ii)%denstype.eq.'TAUFILE') then
		call renormtaufile(ii,SD,Zone(ii)%nr,Zone(ii)%np)
	endif

	deallocate(Aspiral)
	deallocate(A)
	deallocate(H)
	deallocate(SD)
	deallocate(alpha)

	
	return
	end

	subroutine readtaufile(ii,SD,nr,np)
	use GlobalSetup
	use Constants
	integer ii,nr,np,ir,ip
	real*8 SD(nr,np),r,tau0,Mtot
	
	open(unit=20,file=Zone(ii)%taufile,RECL=6000)
	
	ir=1
	read(20,*)
1	read(20,*,end=2) r,tau0
	r=r*1d5*1d6
	if(tau0.lt.1d-5) tau0=1d-6
	do while(r.gt.Zone(ii)%R(ir+1))
		SD(ir,1:np)=tau0
		ir=ir+1
		if(ir.gt.nr) goto 2
	enddo
	goto 1
2	close(unit=20)

	Mtot=0d0
	do ir=1,nr
	do ip=1,np
		A=pi*(Zone(ii)%R(ir+1)**2-Zone(ii)%R(ir)**2)/real(Zone(ii)%np)
		Mtot=Mtot+A*SD(ir,ip)
	enddo
	enddo
	SD=SD*Zone(ii)%Mdust/Mtot

c	open(unit=20,file='surfdens.dat')
c	do ir=1,nr
c		write(20,*) Zone(ii)%R(ir)/1d5/1d6,Zone(ii)%R(ir+1)/1d5/1d6,SD(ir,1)
c	enddo	

	return
	end
	
	subroutine renormtaufile(ii,SD,nr,np)
	use GlobalSetup
	use Constants
	integer ii,nr,np,ir,ip,ilam,i
	real*8 SD(nr,np),d,lam0,r,GetKext
	type(Cell),pointer :: C
	
	open(unit=20,file=Zone(ii)%taufile,RECL=6000)
	
	ir=1
	read(20,*)
1	read(20,*,end=2) r,tau0
	r=r*1d5*1d6
	if(tau0.lt.1d-5) tau0=1d-6
	do while(r.gt.Zone(ii)%R(ir+1))
		SD(ir,1:np)=tau0
		ir=ir+1
		if(ir.gt.nr) goto 2
	enddo
	goto 1
2	close(unit=20)

	lam0=0.55
	d=lam(nlam)-lam(1)
	ilam=1
	do i=1,nlam
		if(abs(lam0-lam(i)).lt.d) then
			d=abs(lam0-lam(i))
			ilam=i
		endif
	enddo
		
	do ir=1,nr
		r=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
		do ip=1,np
			tau0=0d0
			do it=1,Zone(ii)%nt
				C => Zone(ii)%C(ir,it,ip)
				tau0=tau0+GetKext(ilam,C)*r*abs(Zone(ii)%theta(it+1)-Zone(ii)%theta(it))
			enddo
			do it=1,Zone(ii)%nt
				C => Zone(ii)%C(ir,it,ip)
				C%dens=C%dens*SD(ir,ip)/tau0
				C%M=C%M*SD(ir,ip)/tau0
				C%densP=C%densP*SD(ir,ip)/tau0
				C%gasdens=C%gasdens*SD(ir,ip)/tau0
			enddo
		enddo
	enddo
	
	open(unit=20,file='surfdens.dat')
	do ir=1,nr
		tot=0d0
		do it=1,Zone(ii)%nt
			C => Zone(ii)%C(ir,it,1)
			tot=tot+C%dens*C%V
		enddo
		tot=tot*real(Zone(ii)%np)/(pi*(Zone(ii)%R(ir+1)**2-Zone(ii)%R(ir)**2))
		write(20,*) Zone(ii)%R(ir)/1d5/1d6,tot
	enddo	
	close(unit=20)
	
	return
	end
	

	subroutine ComputeAvortex(ii,Avortex,SD,H)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,ipart,isize,ip,ir,it
	real*8 Avortex(Zone(ii)%nr,Zone(ii)%np,npart,0:maxns),scale(npart,0:maxns)
	real*8 Hv,a,fx,St,wv,SD(Zone(ii)%nr,Zone(ii)%np),H(Zone(ii)%nr,Zone(ii)%np)
	real*8 phi,RR,a1,a2
	
	if(Zone(ii)%rvortex.lt.0d0) then
		Avortex=1d0
		return
	endif

	scale=0d0
	do ir=1,Zone(ii)%nr
		do ip=1,Zone(ii)%np
			RR=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
			phi=360d0*(real(ip)-0.5d0)/real(Zone(ii)%np)
			do ipart=1,npart
				do isize=0,Part(ipart)%nsize
					if(isize.eq.0) then
						St=0d0
					else
						St=Part(ipart)%rv(isize)*Part(ipart)%rho(1)*pi/(2d0*Zone(ii)%gas2dust*SD(ir,ip))
					endif
					wv=3d0/(2d0*(Zone(ii)%avortex-1d0))		! Kida
c					wv=sqrt(3d0/(Zone(ii)%avortex**2-1d0))	! GNG
					fx=2d0*wv*Zone(ii)%avortex-(2d0*wv**2+3d0)/(1d0+1d0/Zone(ii)%avortex**2)
					fx=sqrt(fx)
					if(Zone(ii)%dvortex.lt.0d0) then
						delta_St=Zone(ii)%alpha/(1d0+St)
					else
						delta_St=Zone(ii)%dvortex
					endif
					Hv=H(ir,ip)*sqrt(delta_St/(St+delta_St))/fx
					a1=abs(RR-Zone(ii)%rvortex)
					a2=abs(phi-Zone(ii)%phivortex)
					if(a2.gt.180d0) a2=360d0-a2
					a2=a2*pi*RR/180d0/Zone(ii)%avortex
					a=sqrt(a1**2+a2**2)
					Avortex(ir,ip,ipart,isize)=exp(-a**2/(2d0*Hv**2))
					do it=1,Zone(ii)%nt
						scale(ipart,isize)=scale(ipart,isize)+Avortex(ir,ip,ipart,isize)*Zone(ii)%C(ir,it,ip)%V
					enddo
				enddo
			enddo
		enddo
	enddo

	do ir=1,Zone(ii)%nr
		do ip=1,Zone(ii)%np
			do ipart=1,npart
				do isize=0,Part(ipart)%nsize
					Avortex(ir,ip,ipart,isize)=Avortex(ir,ip,ipart,isize)/scale(ipart,isize)
				enddo
			enddo
		enddo
	enddo
			
	
	return
	end
	
	subroutine ComputeAspiral(ii,ispiral,Aspiral,nr,np)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ip,nir,ir,jj,ii,nr,np,nip,kk,i,n,ispiral
	real*8 r,hr,phi0,Aspiral(nSpirals,nr,np),f,phimin,phimax,beta
	real*8,allocatable :: phi(:),rp(:)

	nir=10
	nip=10

	call output("Setting up spiral wave " // int2string(ispiral,'(i3)'))

	allocate(phi(nir*nr))
	allocate(rp(nir*nr))

	i=0
	if(Spiral(ispiral)%hr.lt.0d0) then
		hr=(Zone(ii)%sh*(Spiral(ispiral)%r/Zone(ii)%Rsh)**Zone(ii)%shpow)/Spiral(ispiral)%r
	else
		hr=Spiral(ispiral)%hr
	endif

	if(Spiral(ispiral)%beta.lt.0d0) then
		beta=1.5d0-Zone(ii)%shpow
	else
		beta=Spiral(ispiral)%beta
	endif

	phimin=1d200
	phimax=-1d200
	do ir=1,nr
		do jj=1,nir
			i=i+1
			f=(real(jj)-0.5)/real(nir)
			r=sqrt(Zone(ii)%R(ir)**(2d0-2d0*f)*Zone(ii)%R(ir+1)**(2d0*f))
			rp(i)=r
			phi(i)=Spiral(ispiral)%phi-Spiral(ispiral)%sign*(sign(1d0,r-Spiral(ispiral)%r)/hr)*(
     &			(r/Spiral(ispiral)%r)**(1d0+beta)*
     &			(1d0/(1d0+beta)-(1d0/(1d0-Spiral(ispiral)%alpha+beta))*
     &			(r/Spiral(ispiral)%r)**(-Spiral(ispiral)%alpha))
     &			-(1d0/(1d0+beta)-1d0/(1d0-Spiral(ispiral)%alpha+beta)))
			if(phi(i).lt.phimin) phimin=phi(i)
			if(phi(i).gt.phimax) phimax=phi(i)
		enddo
	enddo
	n=i


	call tellertje(1,np)
!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(ip,ir,phi0,f,r,i)
!$OMP& SHARED(Zone,Spiral,ii,ispiral,phimin,np,nr,rp,phi,phimax,nir,nip,Aspiral,n)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
	do ip=1,np
!$OMP CRITICAL
		call tellertje(ip+1,np+2)
!$OMP END CRITICAL
		do ir=1,nr
			Aspiral(ispiral,ir,ip)=0d0
		enddo
		phi0=2d0*pi*(real(ip)-0.5)/real(np)
		do while(phi0.gt.phimin)
			do i=1,n-1
				if(sign(1d0,phi0-phi(i))*sign(1d0,phi0-phi(i+1)).lt.0d0) goto 1
			enddo
			goto 2
1			continue
			do ir=1,nr
				do jj=1,nir
					f=(real(jj)-0.5)/real(nir)
					r=sqrt(Zone(ii)%R(ir)**(2d0-2d0*f)*Zone(ii)%R(ir+1)**(2d0*f))
					Aspiral(ispiral,ir,ip)=Aspiral(ispiral,ir,ip)+
     &						((r/Spiral(ispiral)%r)**(sign(1d0,Spiral(ispiral)%r-r)*Spiral(ispiral)%q))*
     &						exp(-(r-rp(i))**2/Spiral(ispiral)%w**2)/real(nir)				
				enddo
			enddo
2			phi0=phi0-2d0*pi
		enddo
		phi0=2d0*pi*(real(ip)-0.5)/real(np)+2d0*pi
		do while(phi0.lt.phimax)
			do i=1,n-1
				if(sign(1d0,phi0-phi(i))*sign(1d0,phi0-phi(i+1)).lt.0d0) goto 3
			enddo
			goto 4
3			continue
			do ir=1,nr
				do jj=1,nir
					f=(real(jj)-0.5)/real(nir)
					r=sqrt(Zone(ii)%R(ir)**(2d0-2d0*f)*Zone(ii)%R(ir+1)**(2d0*f))
					Aspiral(ispiral,ir,ip)=Aspiral(ispiral,ir,ip)+
     &					((r/Spiral(ispiral)%r)**(sign(1d0,Spiral(ispiral)%r-r)*Spiral(ispiral)%q))*
     &					exp(-(r-rp(i))**2/Spiral(ispiral)%w**2)/real(nir)				
				enddo
			enddo
4			phi0=phi0+2d0*pi
		enddo
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
	call tellertje(np,np)

	return
	end

	subroutine SetupRadGrid(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i,nspan,nlev,j,ir

	nspan=Zone(ii)%nr/21
	nlev=7
	if(Zone(ii)%thin) then
		nspan=1
		nlev=1
	endif
	if(Zone(ii)%nr.lt.nspan*nlev) then
		call output("You need more than " // trim(int2string(nspan*nlev,'(i2)')) // " radial points")
		stop
	endif

c setup initial radial grid
	do i=1,Zone(ii)%nr+1-nspan*nlev
		Zone(ii)%R(i)=10d0**(log10(Zone(ii)%Rin)+real(i-1)*log10(Zone(ii)%Rout/Zone(ii)%Rin)/real(Zone(ii)%nr-nspan*nlev))
	enddo
	
	ir=Zone(ii)%nr+1-nspan*nlev
	do j=1,nlev
		do i=1,nspan
			ir=ir+1
			Zone(ii)%R(ir)=sqrt(Zone(ii)%R(i)*Zone(ii)%R(i+1))
		enddo
		call sort(Zone(ii)%R(1:ir),ir)
	enddo

	open(unit=20,file=trim(outputdir) // 'radgrid' // trim(int2string(ii,'(i0.4)')) // '.dat')
	write(20,'("# radial grid")')
	write(20,'("# units ",a)') trim(Zone(ii)%sscaletype)
	do i=1,Zone(ii)%nr+1
		write(20,*) Zone(ii)%R(i)/Zone(ii)%sscale
	enddo		
	close(unit=20)

	if((Zone(ii)%x0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%x0)+Zone(ii)%Rout)
	if((Zone(ii)%y0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%y0)+Zone(ii)%Rout)
	if((Zone(ii)%z0+Zone(ii)%Rout).gt.maxR) maxR=(abs(Zone(ii)%z0)+Zone(ii)%Rout)

	do i=1,Zone(ii)%nr+1
		Zone(ii)%R2(i)=Zone(ii)%R(i)**2
	enddo

	return
	end
	

	subroutine SetupThetaGrid(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i
	real*8 tmax,random
	
	if(Zone(ii)%shape.eq.'SPH') then
		tmax=pi/2d0
	else if(Zone(ii)%shape.eq.'CYL') then
		tmax=Zone(ii)%tmax*pi/180d0
	endif

	if((2*int(Zone(ii)%nt/2)).ne.Zone(ii)%nt) then
		Zone(ii)%nt=Zone(ii)%nt-1
		call output("I need a midplane boundary, so an even number of theta cells")
		call output("changing nt to: " // trim(int2string(Zone(ii)%nt,'(i4)')))
	endif

c setup initial theta grid
	do i=1,Zone(ii)%nt+1
		Zone(ii)%theta(i)=tmax*(2d0*real(i-1)/real(Zone(ii)%nt)-1d0)**3+pi/2d0
	enddo

	Zone(ii)%imidplane=Zone(ii)%nt/2+1

	call sort(Zone(ii)%theta,Zone(ii)%nt+1)

	open(unit=20,file=trim(outputdir) // 'thetagrid' // trim(int2string(ii,'(i0.4)')) // '.dat')
	do i=1,Zone(ii)%nt+1
		write(20,*) Zone(ii)%theta(i)
	enddo
	close(unit=20)

	do i=1,Zone(ii)%nt+1
		Zone(ii)%cost2(i)=cos(Zone(ii)%theta(i))**2
	enddo
	
	return
	end
	

	subroutine SetupThetaGridEqui(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i
	real*8 tmax,random
	
	if(Zone(ii)%shape.eq.'SPH') then
		tmax=pi/2d0
	else if(Zone(ii)%shape.eq.'CYL') then
		tmax=Zone(ii)%tmax*pi/180d0
	endif

	if((2*int(Zone(ii)%nt/2)).ne.Zone(ii)%nt) then
		Zone(ii)%nt=Zone(ii)%nt-1
		call output("I need a midplane boundary, so an even number of theta cells")
		call output("changing nt to: " // trim(int2string(Zone(ii)%nt,'(i4)')))
	endif

c setup initial theta grid
	do i=1,Zone(ii)%nt+1
		Zone(ii)%theta(i)=tmax*(2d0*real(i-1)/real(Zone(ii)%nt)-1d0)+pi/2d0
	enddo

	Zone(ii)%imidplane=Zone(ii)%nt/2+1

	call sort(Zone(ii)%theta,Zone(ii)%nt+1)

	open(unit=20,file=trim(outputdir) // 'thetagrid' // trim(int2string(ii,'(i0.4)')) // '.dat')
	do i=1,Zone(ii)%nt+1
		write(20,*) Zone(ii)%theta(i)
	enddo
	close(unit=20)

	do i=1,Zone(ii)%nt+1
		Zone(ii)%cost2(i)=cos(Zone(ii)%theta(i))**2
	enddo
	
	return
	end
	


	subroutine SetupPhiGrid(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i

c setup initial phi grid
	open(unit=20,file=trim(outputdir) // 'phigrid' // trim(int2string(ii,'(i0.4)')) // '.dat')
	do i=1,Zone(ii)%np+1
		Zone(ii)%phi(i)=2d0*pi*real(i-1)/real(Zone(ii)%np)
		write(20,*) Zone(ii)%phi(i)
	enddo		
	close(unit=20)

	do i=1,Zone(ii)%np+1
		if(Zone(ii)%phi(i).lt.(pi/4d0)) then
			Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
			Zone(ii)%tany(i)=1d0
		else if(Zone(ii)%phi(i).lt.(3d0*pi/4d0)) then
			Zone(ii)%tanx(i)=-1d0
			Zone(ii)%tany(i)=-cos(Zone(ii)%phi(i))/sin(Zone(ii)%phi(i))
		else if(Zone(ii)%phi(i).lt.(5d0*pi/4d0)) then
			Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
			Zone(ii)%tany(i)=1d0
		else if(Zone(ii)%phi(i).lt.(7d0*pi/4d0)) then
			Zone(ii)%tanx(i)=-1d0
			Zone(ii)%tany(i)=-cos(Zone(ii)%phi(i))/sin(Zone(ii)%phi(i))
		else
			Zone(ii)%tanx(i)=sin(Zone(ii)%phi(i))/cos(Zone(ii)%phi(i))
			Zone(ii)%tany(i)=1d0
		endif
	enddo
		
	return
	end


	subroutine SetupVolume(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ir,it,ip,ii
	
	do ir=1,Zone(ii)%nr
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				if(Zone(ii)%shape.eq.'SPH') then
					Zone(ii)%C(ir,it,ip)%V=(2d0*pi/3d0)*(Zone(ii)%R(ir+1)**3-Zone(ii)%R(ir)**3)*
     &					abs(cos(Zone(ii)%theta(it))-cos(Zone(ii)%theta(it+1)))/real(Zone(ii)%np)
     				Zone(ii)%C(ir,it,ip)%V=Zone(ii)%C(ir,it,ip)%V*Zone(ii)%xscale*Zone(ii)%yscale*Zone(ii)%zscale
				else if(Zone(ii)%shape.eq.'CYL') then
					call output("Still have to do this...")
					stop
				endif
    		enddo
		enddo
	enddo
	
	return
	end
	

	subroutine SetupSphere(ii)
	IMPLICIT NONE
	integer ii
	call SetupShell(ii)

	return
	end
	
	
	subroutine SetupShell(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,ir,it,ip,i,ips,ipt,ilam
	real*8 r,z,hr,f1,f2,Mtot,w(npart,maxns),tau,GetKext,lam0,d

	if(Zone(ii)%shape.eq.'CAR'.or.Zone(ii)%shape.eq.'CYL') then
		call output("Free advice: a spherical shell is better put on a spherical grid.")
		stop
	endif

	call SetupRadGrid(ii)
	call SetupThetaGrid(ii)
	call SetupPhiGrid(ii)
	call SetupVolume(ii)

	Mtot=0d0
	do i=1,npart
		if(Part(i)%nsize.gt.1) then
			call SetSizeDis(w(i,1:Part(i)%nsize),i,ii)
		else
			w(i,1)=1d0
		endif
	enddo
	
	do ir=1,Zone(ii)%nr
		call tellertje(ir,Zone(ii)%nr)
		r=sqrt(Zone(ii)%R(ir)*Zone(ii)%R(ir+1))
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				Zone(ii)%C(ir,it,ip)%dens=r**(-Zone(ii)%denspow)*exp(-(r/Zone(ii)%Rexp)**2)
				Mtot=Mtot+Zone(ii)%C(ir,it,ip)%dens*Zone(ii)%C(ir,it,ip)%V
    		enddo
		enddo
	enddo
	do ir=1,Zone(ii)%nr
		do it=1,Zone(ii)%nt
			do ip=1,Zone(ii)%np
				Zone(ii)%C(ir,it,ip)%dens=Zone(ii)%C(ir,it,ip)%dens*Zone(ii)%Mdust/Mtot
				do i=1,npart
					do ips=1,Part(i)%nsize
						Zone(ii)%C(ir,it,ip)%densP(i,ips,1)=w(i,ips)*Zone(ii)%abun(i)*Zone(ii)%C(ir,it,ip)%dens
						do ipt=2,Part(i)%nT
							Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)=0d0
						enddo
					enddo
				enddo
				Zone(ii)%C(ir,it,ip)%M=Zone(ii)%C(ir,it,ip)%dens*Zone(ii)%C(ir,it,ip)%V
			enddo
		enddo
	enddo

	if(Zone(ii)%tau_V.gt.0d0) then
		lam0=0.55
		d=lam(nlam)-lam(1)
		ilam=1
		do i=1,nlam
			if(abs(lam0-lam(i)).lt.d) then
				d=abs(lam0-lam(i))
				ilam=i
			endif
		enddo
	
		tau=0d0
		do ir=1,Zone(ii)%nr
			tau=tau+GetKext(ilam,Zone(ii)%C(ir,Zone(ii)%imidplane,1))*
     &					(Zone(ii)%R(ir+1)-Zone(ii)%R(ir))
		enddo
		do ir=1,Zone(ii)%nr
			do it=1,Zone(ii)%nt
				do ip=1,Zone(ii)%np
					Zone(ii)%C(ir,it,ip)%dens=Zone(ii)%C(ir,it,ip)%dens*Zone(ii)%tau_V/tau
					Zone(ii)%C(ir,it,ip)%M=Zone(ii)%C(ir,it,ip)%M*Zone(ii)%tau_V/tau
					do i=1,npart
						do ips=1,Part(i)%nsize
							do ipt=1,Part(i)%nT
								Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)=
     &								Zone(ii)%C(ir,it,ip)%densP(i,ips,ipt)*Zone(ii)%tau_V/tau
							enddo
						enddo
					enddo
				enddo
			enddo
		enddo
	endif

	
	return
	end

c-----------------------------------------------------------------------
c
c   This subroutine calculates the shape for a rounded-off gap edge
c   It returns only a factor which is the decrease in density
c
c   + type='POW': a powerlaw surface density between up to rmax,
c                 proportional to r^-index.
c
c   + type='HYDRO':  a gaussian like shape that fits hydro simulations
c                    of embedded planets (Mulders et al. 2013)
c 					index is the width parameter (>0)
c
c
c   OUPUT:  scaling factor at radius r
c
c-----------------------------------------------------------------------

	real*8 function RoundOff(r,rmax,roundtype,index,depth)
	IMPLICIT NONE
	character*10 roundtype
	real*8 scale,r,rmax,index,depth
	
	if(r.ge.rmax) then
		!  do not apply a scale factor if outside rmax
		scale= 1d0
	else if (roundtype.eq.'HYDRO') then
		!  use a shape from Hydrodynamic simulations
		scale=exp( -(((1-r/rmax)/index)**3d0))
	else if (roundtype.eq.'POW') then
		!  use a powerlaw
		scale=(rmax/r)**(-index)
	else
		call output("Gap shape not recognized")
		stop
	endif
	
	!  return gap depth
	RoundOff=max(scale,depth)

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine sort(x,n)
	IMPLICIT NONE
	integer n,i,j,imin
	real*8 x(n),min
	
	do j=1,n-1
	min=x(j)
	imin=j
	do i=j,n
		if(x(i).lt.min) then
			min=x(i)
			imin=i
		endif
	enddo
	min=x(j)
	x(j)=x(imin)
	x(imin)=min
	enddo
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Planck(T,lam)
	IMPLICIT NONE
	real*8 T,k,c,h,nu,lam,pi
	real*16 x

	k=1.3807d-16
	c=2.9979d10
	h=6.6261d-27
	nu=c/(lam*1d-4)
	x=h*nu/(k*T)
	pi=3.1415926536

	if(x.gt.1d3) then
		Planck=0d0
	else
		Planck=4d0*pi*((2d0*h*nu**3/c**2)/(exp(x)-1d0))
	endif

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	real*8 function Luminosity(T,R)
	IMPLICIT NONE
	real*8 T,R,pi,sigma

	pi=3.1415926536
	sigma=5.6704d-5

	Luminosity=4d0*pi*sigma*R**2*T**4
	
	return
	end



	subroutine SetupBeaming()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,istar
	real*8 f,E
	
	do izone=1,nzones
		Zone(izone)%fbeam=Zone(izone)%fbeam/real(nzones)
		call DetermineBeamZ(izone)
	enddo

	return
	end
	


	subroutine DetermineBeamZ(izone)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,istar,i
	real*8 f,x,y,z,vx,vy,vz,r,d
	
	allocate(Zone(izone)%fbeamS(nstars))
	allocate(Zone(izone)%tbeamS(nstars))
	allocate(Zone(izone)%ctbeamS(nstars))
	allocate(Zone(izone)%EfbeamS(nstars))

	do istar=1,nstars
		r=(Zone(izone)%xscale*(Star(istar)%x-Zone(izone)%x0))**2+
     &    (Zone(izone)%yscale*(Star(istar)%y-Zone(izone)%y0))**2+
     &    (Zone(izone)%zscale*(Star(istar)%z-Zone(izone)%z0))**2
		if(r.lt.Zone(izone)%Rout**2) then
			Zone(izone)%fbeamS(istar)=0d0
			Zone(izone)%tbeamS(istar)=0d0
			Zone(izone)%EfbeamS(istar)=0d0
		else
			Zone(izone)%fbeamS(istar)=Zone(izone)%fbeam
			d=sqrt((Zone(izone)%x0-Star(istar)%x)**2
     &	+(Zone(izone)%y0-Star(istar)%y)**2+(Zone(izone)%z0-Star(istar)%z)**2)
			Zone(izone)%tbeamS(istar)=atan((Zone(izone)%Rout+Star(istar)%R)/d)
			Zone(izone)%ctbeamS(istar)=cos(Zone(izone)%tbeamS(istar))
			f=(1d0-Zone(izone)%ctbeamS(istar))/2d0
			Zone(izone)%EfbeamS(istar)=f
		endif
	enddo
	
	return
	end
	
	
	