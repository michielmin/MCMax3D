	subroutine OutputMCobs()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	character*500 MCfile
	integer i,j,nf
	real*8 f,x,y,z,inp,Reddening,compute_dlam

	call output("==================================================================")
	call output("Writing Monte Carlo observables")
	do i=1,nMCobs	
		MCfile=trim(outputdir) // "MCout" // trim(int2string(i,'(i0.4)')) // ".fits.gz"
		if(MCobs(i)%writeimage) call writefitsfile(MCfile,MCobs(i)%image,nlam,MCobs(i)%npix)
		MCfile=trim(outputdir) // "MCSpec" // trim(int2string(i,'(i0.4)')) // ".dat"
		open(unit=20,file=MCfile)
		do j=1,nlam
			MCobs(i)%spec(j)=sum(MCobs(i)%image(:,:,j))
			write(20,*) lam(j),MCobs(i)%spec(j)*Reddening(lam(j),compute_dlam(lam(j)),Av)/distance**2
		enddo
		close(unit=20)
	enddo
	
	return
	end
	
