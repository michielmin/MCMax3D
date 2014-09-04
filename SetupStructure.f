	subroutine SetupStructure
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 T,Planck
	
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
	xobs=cos(phi_obs)*sin(theta_obs)
	yobs=sin(phi_obs)*sin(theta_obs)
	zobs=cos(theta_obs)
	
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
	
	if(nzlam.le.0) then
		do i=1,nlam
			lam(i)=10d0**(log10(lam1)+(log10(lam2)-log10(lam1))*real(i-1)/real(nlam-1))
		enddo
	else
c	does not seem to work!!! Have to fix this!!!
		call output("NZLAM OPTION NOT ALWAYS WORKING PROPERLY YET!!!")
		nl=(nlam-nzlam)*(log10(zlam1/lam1)/log10(lam2*zlam1/(zlam2*lam1)))
		do i=1,nl
			lam(i)=10d0**(log10(lam1)+(log10(zlam1)-log10(lam1))*real(i-1)/real(nl-1))
		enddo
		j=i-1
		nl=(nlam-nzlam)-nl
		do i=1,nl
			lam(i+j)=10d0**(log10(zlam2)+(log10(lam2)-log10(zlam2))*real(i-1)/real(nl-1))
		enddo
		j=j+i-1
		do i=1,nzlam
			lam(i+j)=10d0**(log10(zlam1)+(log10(zlam2)-log10(zlam1))*real(i)/real(nzlam+1))
		enddo
		call sort(lam,nlam)
	endif

	do i=1,nlam
		nu(i)=1d4*clight/lam(i)
	enddo

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
	integer ii,iT,is,iBB
	real*8 spec(nlam)
	
	allocate(Part(ii)%rv(Part(ii)%nsize,Part(ii)%nT))
	allocate(Part(ii)%rho(Part(ii)%nT))
	allocate(Part(ii)%Kabs(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Ksca(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Kext(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%F(Part(ii)%nsize,Part(ii)%nT,nlam))
	allocate(Part(ii)%Kp(Part(ii)%nsize,Part(ii)%nT,nBB))
	
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
				spec=BB(1:nlam,iT)*Part(ii)%Kabs(is,iT,1:nlam)
				call integrate(spec,Part(ii)%Kp(is,iT,iBB))
			enddo
		enddo
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


	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	

	subroutine SetupZone(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii
	
	select case(Zone(ii)%shape)
		case("SPH","CYL")
			allocate(Zone(ii)%C(Zone(ii)%nr,Zone(ii)%nt,Zone(ii)%np))
			allocate(Zone(ii)%R(Zone(ii)%nr+1))
			allocate(Zone(ii)%theta(Zone(ii)%nt+1))
			allocate(Zone(ii)%phi(Zone(ii)%np+1))
		case("CAR")
			allocate(Zone(ii)%C(Zone(ii)%nx,Zone(ii)%ny,Zone(ii)%nz))
			allocate(Zone(ii)%x(Zone(ii)%nx+1))
			allocate(Zone(ii)%y(Zone(ii)%ny+1))
			allocate(Zone(ii)%z(Zone(ii)%nz+1))
		case default
			call output("Interesting shape option ("// Zone(ii)%shape //"). But I don't get it...")
			stop
	end select
	
	select case(Zone(ii)%denstype)
		case("DISK")
			call SetupDisk(ii)
		case("SHELL")
			call SetupShell(ii)
		case default
			call output("Really? A " // trim(Zone(ii)%denstype) // "-zone? That is an awesome idea! (but not yet possible)")
			stop
	end select	
	
	return
	end
	

	subroutine SetupDisk(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,i

	if(Zone(ii)%shape.eq.'CAR') then
		call output("A disk on a rectangular grid... Let's don't and say we did.")
		stop
	endif

c setup initial radial grid
	do i=1,Zone(ii)%nr+1
		Zone(ii)%R(i)=10d0**(log10(Zone(ii)%Rin)+real(i-1)*log10(Zone(ii)%Rout/Zone(ii)%Rin)/real(Zone(ii)%nr))
	enddo		


	
	return
	end


	subroutine SetupShell(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii
	
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
		Planck=(2d0*h*nu**3/c**2)/(exp(x)-1d0)
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



