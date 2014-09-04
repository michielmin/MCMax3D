	subroutine SetupStructure
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i,j
	real*8 T,Planck
	
	call SetupLam()
	
	allocate(BB(nlam,0:nBB))
	do i=1,nlam
		BB(i,0)=0d0
		do j=1,nBB
			T=dTBB*real(j)
			BB(i,j)=Planck(T,lam(i))
		enddo
	enddo
	
	do i=1,npart
		call SetupPart(i)
	enddo

	do i=1,nstars
c		call SetupStar(i)
	enddo
	
	do i=1,nzones
c		call SetupZone(i)
	enddo
	
	return
	end


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
	real*8 T,k,c,h,nu,lam
	real*16 x

	k=1.3807d-16
	c=2.9979d10
	h=6.6261d-27
	nu=c/(lam*1d-4)
	x=h*nu/(k*T)
	if(x.gt.1d3) then
		Planck=0d0
	else
		Planck=(2d0*h*nu**3/c**2)/(exp(x)-1d0)
	endif
c	Planck=Planck*1e23

	return
	end



