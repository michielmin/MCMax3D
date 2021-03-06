	subroutine outputform(string,form)
	IMPLICIT NONE
	character string*(*)
	character,intent(in),optional :: form*(*)

	if(form.ne.' ') then
		write(*,form) trim(string)
		write(9,form) trim(string)
	else
		write(*,'(a)') trim(string)
		write(9,'(a)') trim(string)
	endif
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine flushoutput()
	IMPLICIT NONE
	
	call flush(6)
	call flush(9)

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine output(string)
	IMPLICIT NONE
	character string*(*)

	write(*,'(a)') trim(string)
	write(9,'(a)') trim(string)
	
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine output_erase(string)
	IMPLICIT NONE
	character string*(*)

	write(*,'(1a1,a,$)') char(13),trim(string)
	write(9,'(1a1,a,$)') char(13),trim(string)
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine ignorestar(un)
	IMPLICIT NONE
	integer un
	character c
1	read(un,fmt=3,end=2) c
	if(c.eq.'*') goto 1
	backspace(unit=un)
2	continue
3	format(a1)
	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	character*20 function int2string(i,form)
	IMPLICIT NONE
	integer i
	character,intent(in),optional :: form*(*)
	
	if(form.ne.' ') then
		write(int2string,form) i
	else
		write(int2string,*) i
	endif
	
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	character*20 function dbl2string(x,form)
	IMPLICIT NONE
	real*8 x
	character,intent(in),optional :: form*(*)
	
	if(form.ne.' ') then
		write(dbl2string,form) x
	else
		write(dbl2string,*) x
	endif
	
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje(i,n)
	IMPLICIT NONE
	integer i,n,f
	
c GFORTRAN requires interface for this subroutine
	interface
		subroutine outputform(string,form)
			character string*(*)
			character,intent(in),optional :: form*(*)
		end subroutine
	end interface
	
	if(i.eq.1) call output("....................")

	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.20d0*dble(i+1)/dble(n).gt.dble(f)) then
		call outputform(".",'(a1,$)')
		call flushoutput
	endif
	
	if(i.eq.n) call output("")

	return
	end

	logical function checktellertje(i,n)
	IMPLICIT NONE
	integer i,n,f

	checktellertje=.false.
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.20d0*dble(i+1)/dble(n).gt.dble(f)) then
		checktellertje=.true.
	endif

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje_time(i,n,starttime,starttime_w)
	IMPLICIT NONE
	integer i,n,f
	real*8 starttime,stoptime,xx,starttime_w,stoptime_w,omp_get_wtime
	
c GFORTRAN requires interface for this function
	INTERFACE 
		character*20 function dbl2string(x,form)
			real*8 x
			character,intent(in),optional :: form*(*)
		end function dbl2string
	end INTERFACE
	
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*dble(i-1)/dble(n).lt.dble(f)
     &   .and.20d0*dble(i+1)/dble(n).gt.dble(f)) then
		call cpu_time(stoptime)
		stoptime_w=omp_get_wtime()
		xx=100d0*dble(i)/dble(n)
c		call output(trim(dbl2string(1000d0*(stoptime-starttime)/real(i),'(f8.3)'))
c     &			//" ms per photon package. Approx " // 
c     &			trim(dbl2string((stoptime-starttime)*(n-i)/real(i),'(f8.2)'))
c     &			//" s left. (" //
c     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
		if(i.eq.1) then
			call output(trim(dbl2string(1000d0*(stoptime_w-starttime_w)/dble(i),'(f8.3)'))
     &			//" ms per photon package. Approx " // 
     &			trim(dbl2string((stoptime_w-starttime_w)*(n-i)/dble(i),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
		else
			call output_erase(trim(dbl2string(1000d0*(stoptime_w-starttime_w)/dble(i),'(f8.3)'))
     &			//" ms per photon package. Approx " // 
     &			trim(dbl2string((stoptime_w-starttime_w)*(n-i)/dble(i),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
		endif
	endif
	if(i.eq.n) then
		call output("")
		call output(trim(dbl2string((stoptime_w-starttime_w),'(f8.2)'))
     &				//" s walltime. " // 
     &				trim(dbl2string((stoptime-starttime),'(f8.2)'))
     &				//" s cpu time. ")
	endif

	return
	end



