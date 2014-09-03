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

