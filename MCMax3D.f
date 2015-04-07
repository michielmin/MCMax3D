	program MCMax3D
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical converged
	integer iter,i
	character*500 VersionGIT,filename

	criticalerror=.false.

	call GetOutputDir
	open(unit=9,file=trim(outputdir) // "log.dat",RECL=6000)

c terms of use
	call output("==================================================================")
	call output("By using MCMax3D you agree to the terms of use.")
	call output("It basically means you offer me co-author rights on any paper.")
	call output("that uses results computed with MCMax3D.")

	call output("==================================================================")
	call output("Let's get the show on the road!!")
	call output("MCMax3D "//trim(VersionGIT()))
	call output("==================================================================")

c let's initialize everything first
	call Initialize
c now setup the basic structure
	call SetupStructure
	
	call OutputStats
	
	converged=.false.
	iter=0

	if(Nphot.gt.0) then
c iterate when needed
		do while(.not.converged.and.iter.le.maxiter)
c compute the next iteration of the structure
			call EvolveStructure(iter)
c do the radiative transfer
			call RadiativeTransfer

c is everything converged? (or iterations not needed)
c			converged=DetermineConverged(iter)
	
c write the structure to the output files
c			call OutputStructure(iter,converged)
			iter=iter+1
			call OutputMCobs
		enddo

c		call createUV()

		do i=1,nzones
			filename=trim(outputdir) // "Zone" // trim(int2string(i,'(i0.4)')) // ".fits.gz"
			call output("Writing file: "// trim(filename))
			call outputstruct_fits(filename,ZoneStructOutput,nZoneStructOutput,i)
		enddo
	endif

c ok, structure is done, now let's see what this looks like
	do i=1,nMCobs
		if(MCobs(i)%raytrace) call Raytrace(i)
	enddo

c well that's it. we seem to be done!
c have a good day
	call output("==================================================================")
	call output("Success!!")
	call output("Have a nice day!")
	call output("==================================================================")

	call system("touch " // trim(outputdir) // "done")

	end



