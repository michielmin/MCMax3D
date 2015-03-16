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
	
	if(i.eq.1) call output("....................")

	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
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
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
		checktellertje=.true.
	endif

	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine tellertje_time(i,n,starttime)
	IMPLICIT NONE
	integer i,n,f
	real*8 starttime,stoptime,xx
	character*20 dbl2string
	
	if(i.eq.1) then
		call cpu_time(stoptime)
		xx=100d0*real(i)/real(n)
		call output(trim(dbl2string(1000d0*(stoptime-starttime)/real(i),'(f8.3)'))
     &			//" ms per photon package. Approx " // 
     &			trim(dbl2string((stoptime-starttime)*(n-i)/real(i),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
	endif
	
	f=int(20d0*dble(i)/dble(n))
	
	if(20d0*real(i-1)/real(n).lt.real(f)
     &   .and.20d0*real(i+1)/real(n).gt.real(f)) then
		call cpu_time(stoptime)
		xx=100d0*real(i)/real(n)
		call output(trim(dbl2string(1000d0*(stoptime-starttime)/real(i),'(f8.3)'))
     &			//" ms per photon package. Approx " // 
     &			trim(dbl2string((stoptime-starttime)*(n-i)/real(i),'(f8.2)'))
     &			//" s left. (" //
     &			trim(dbl2string(xx,'(f5.1)')) // " %)")
	endif

	return
	end



