	subroutine OutputMCobs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	character*500 MCfile
	integer i,j,nf
	real*8 f,x,y,z,inp

	call output("==================================================================")
	call output("Writing Monte Carlo observables")
	do i=1,nMCobs
	
		f=0d0
		nf=100000
		do j=1,nf
			call randomdirection(x,y,z)
			inp=MCobs(i)%x*x+MCobs(i)%y*y+MCobs(i)%z*z
			if(inp.gt.MCobs(i)%opening) f=f+1d0
		enddo
		f=4d0*pi*f/real(nf)
	
		MCfile=trim(outputdir) // "MCout" // trim(int2string(i,'(i0.4)')) // ".fits"
		call writefitsfile(MCfile,MCobs(i)%image,nlam,MCobs(i)%npix)
		MCfile=trim(outputdir) // "MCSpec" // trim(int2string(i,'(i0.4)')) // ".dat"
		open(unit=20,file=MCfile)
		do j=1,nlam
			write(20,*) lam(j),MCobs(i)%spec(j)/distance**2/f
		enddo
		close(unit=20)
	enddo
	
	return
	end
	
