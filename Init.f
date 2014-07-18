c=========================================================================================
c This subroutine reads in the input file. In the input file all other input files are
c specified.
c=========================================================================================
	subroutine Initialize()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first

	call SetDefaults

	allocate(key)
	first => key
	
	call GetKeywords(key)
c Count the number of zones, particles, and stars
c allocate the arrays
	key => first%next
	call CountStuff(key)

	key => first%next

	do while(.not.key%last)

	select case(key%key1)
		case("zone")
			call ReadZone(key)
		case default
			call output("Unknown keyword: " // trim(key%key1))
	end select

	key => key%next
	
	enddo
	

	call output("==================================================================")
	
	return
	end

	subroutine ReadZone(key)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey) key

	select case(key%key2)
		case("x")
			read(key%value,*) Zone(key%nr1)%x0
		case("y")
			read(key%value,*) Zone(key%nr1)%y0
		case("z")
			read(key%value,*) Zone(key%nr1)%z0
		case("xn")
			read(key%value,*) Zone(key%nr1)%xn
		case("yn")
			read(key%value,*) Zone(key%nr1)%yn
		case("zn")
			read(key%value,*) Zone(key%nr1)%zn
		case("shape")
			Zone(key%nr1)%shape=key%value(1:3)
		case default
			call output("Unknown zone keyword: " // trim(key%key2))
	end select

	return
	end
	


	subroutine GetKeywords(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer ncla	! number of command line arguments
	character*1000 readline,inputfile

	call getarg(1,inputfile)

	open(unit=20,file=inputfile,RECL=1000)
	
	ncla=-1

	key => firstkey
	
20	ncla=ncla+1
10	continue
	if(ncla.eq.0) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) readline
	else
		call getarg(1+ncla,readline)
		if(readline(1:2).eq.'-s') then
			ncla=ncla+1
			call getarg(1+ncla,readline)
			call output("Command line argument: " // trim(readline))
			ncla=ncla+1
		else
			if(readline.ne.' ') then
c				try to read another command line argument
				ncla=ncla+1
				goto 10
			else
c				all arguments are read
				goto 30
			endif
		endif
	endif

	if(readline.eq.' ') goto 10

	allocate(key%next)
	key => key%next
	key%last=.false.
	call get_key_value(readline,key%key1,key%key2,key%value,key%nr1,key%nr2)

c read another command, so go back
	goto 10

30	continue
	close(unit=20)
	allocate(key%next)
	key => key%next
	key%last=.true.

	return
	end
	

c=========================================================================================
c This subroutine just seperates the key and value component of a string given
c key=value syntax. Key is transformed to lowercase.
c=========================================================================================
	subroutine get_key_value(line,key1,key2,value,nr1,nr2)
	IMPLICIT NONE
	character*1000 line
	character*100 key1,key2,value
	integer i,nr1,nr2,ikey1,ikey2
	
	ikey1=index(line,'=')
	ikey2=index(line,':')

	nr1=1
	nr2=1
	if(ikey2.gt.0) then
		key1=line(1:ikey2-1)
		key2=line(ikey2+1:ikey1-1)
		call checknr(key1,nr1)
		call checknr(key2,nr2)
	else
		key1=line(1:ikey1-1)
		key2=' '
		call checknr(key1,nr1)
	endif


	value=line(index(line,'=')+1:len_trim(line))
	if(value(1:1).eq.'"'.or.value(1:1).eq."'") then
		value=value(2:len_trim(value)-1)
	endif

	do i=1,len_trim(key1)
		if(iachar(key1(i:i)).ge.65.and.iachar(key1(i:i)).le.90) then
			key1(i:i)=achar(iachar(key1(i:i))+32)
		endif
	enddo
	do i=1,len_trim(key2)
		if(iachar(key2(i:i)).ge.65.and.iachar(key2(i:i)).le.90) then
			key2(i:i)=achar(iachar(key2(i:i))+32)
		endif
	enddo

	return
	end
	
	
	subroutine checknr(key,nr)
	IMPLICIT NONE
	character*100 key
	integer nr,i,n
	
	n=len_trim(key)
	i=n
1	read(key(i:n),*,err=2) nr
	i=i-1
	if(i.eq.0) goto 2
	goto 1
2	key=key(1:i)
	if(i.eq.n) nr=1
	
	return
	end
	
c=========================================================================================
c This subroutine sets the default values for the global variables
c=========================================================================================
	subroutine SetDefaults()
	use GlobalSetup
	IMPLICIT NONE
	
	maxiter=0
	
	return
	end
	

	subroutine CountStuff(firstkey)
	use GlobalSetup
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key

	key => firstkey

	nzones=1
	npart=1
	nstars=1
	do while(.not.key%last)
		select case(key%key1)
			case("zone")
				if(key%nr1.gt.nzones) nzones=key%nr1
			case("star")
				if(key%nr1.gt.nstars) nstars=key%nr1
			case("part","computepart","opac")
				if(key%nr1.gt.npart) npart=key%nr1
		end select
		key=>key%next
	enddo

	call output('Number of zones:     ' // int2string(nzones,'(i4)'))
	call output('Number of stars:     ' // int2string(nstars,'(i4)'))
	call output('Number of particles: ' // int2string(npart,'(i4)'))

	allocate(Zone(nzones))
	allocate(Star(nstars))
	allocate(Part(npart))

	return
	end
	