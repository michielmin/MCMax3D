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
		if(MCobs(i)%mcout) then
c			MCfile=trim(outputdir) // "MCout" // trim(int2string(i,'(i0.4)')) // trim(MCobs(i)%flag) // ".fits.gz"
c			if(MCobs(i)%writeimage) call writefitsfile(MCfile,MCobs(i)%image,nlam,MCobs(i)%npix)
			MCfile=trim(outputdir) // "MCSpec" // trim(int2string(i,'(i0.4)')) // trim(MCobs(i)%flag) // ".dat"
			open(unit=20,file=MCfile)
			do j=1,nlam
c				MCobs(i)%spec(j)=sum(MCobs(i)%image(:,:,j))
				write(20,*) lam(j),MCobs(i)%spec(j)*Reddening(lam(j),compute_dlam(lam(j)),Av)/distance**2
			enddo
			close(unit=20)
		endif
	enddo
	
	return
	end
	
