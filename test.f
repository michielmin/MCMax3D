	program test
	IMPLICIT NONE
	character*100 c
	
	call getarg(1,c)
	
	print*,index(c,'=')
	
	end
