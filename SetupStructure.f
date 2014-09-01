	subroutine SetupStructure
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer i
	
	call SetupLam()
	
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
	
	if(nzlam.le.0) then
		do i=1,nlam
			lam(i)=10d0**(log10(lam1)+(log10(lam2)-log10(lam1))*real(i-1)/real(nlam-1))
		enddo
	else
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

	return
	end


	subroutine SetupPart(ii)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer ii,iT
	
	select case(Part(ii)%ptype)
		case("COMPUTE")
			do iT=1,Part(ii)%nT
				call ComputePart(Part(ii),ii,iT)
			enddo
c		case("PARTFILE")
c			call ReadParticle(Part(ii),ii)
		case default
			call output("I did not understand what I was trying to do. Sorry!")
	end select
		
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


