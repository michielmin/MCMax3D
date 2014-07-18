	program MCMax3D
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	logical converged
	integer iter
	character*500 VersionGIT

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
c	call SetupStructure
	
	converged=.false.
	iter=0

c iterate when needed
	do while(.not.converged.and.iter.le.maxiter)
c compute the next iteration of the structure
c		call EvolveStructure(iter)
c do the radiative transfer
c		call RadiativeTransfer

c is everything converged? (or iterations not needed)
c		converged=DetermineConverged(iter)
	
c write the structure to the output files
c		call OutputStructure(iter,converged)
		iter=iter+1
	enddo

c ok, structure is done, now let's see what this looks like
c	call Observations
	
c well that's it. we seem to be done!
c have a good day
	call output("==================================================================")
	call output("Success!!")
	call output("Have a nice day!")
	call output("==================================================================")

	end

