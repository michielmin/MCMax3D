c=========================================================================================
c This subroutine reads in the input file. In the input file all other input files are
c specified.
c=========================================================================================
	subroutine Initialize()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(SettingKey),pointer :: key,first
	character*1000 command
	integer i

	idum=-42
	
	allocate(key)
	first => key
	
	call GetKeywords(key)
c Count the number of zones, particles, and stars
c allocate the arrays
	key => first%next
	call CountStuff(key)

	call SetDefaults

	key => first%next

	do while(.not.key%last)

	select case(key%key1)
		case("zone")
			call ReadZone(key)
		case("star")
			call ReadStar(key)
		case("mcobs")
			call ReadMCobs(key)
		case("part")
			Part(key%nr1)%file=trim(key%value)
			Part(key%nr1)%ptype="PARTFILE"
		case("computepart")
			call ReadComputePart(key)
		case("nlam")
			read(key%value,*) nlam
		case("lam")
			select case(key%nr1)
				case(1)
					read(key%value,*) lam1
				case(2)
					read(key%value,*) lam2
				case default
					call output("Unknown lam value")
			end select
		case("nzlam")
			read(key%value,*) nzlam
		case("zlam")
			select case(key%nr1)
				case(1)
					read(key%value,*) zlam1
				case(2)
					read(key%value,*) zlam2
				case default
					call output("Unknown zlam value")
			end select
		case("maxiter")
			read(key%value,*) maxiter
		case("distance")
			read(key%value,*) distance
		case("dirparticle","particledir")
			particledir=trim(key%value)
		case("fstop")
			read(key%value,*) fstop
		case default
			call output("Unknown keyword: " // trim(key%key1))
			criticalerror=.true.
	end select

	key => key%next
	
	enddo
	

	call output("==================================================================")

	if(criticalerror) then
		call output("Critical error encountered")
		stop
	endif

	maxns=0d0
	maxnT=0d0
	do i=1,npart
		if(Part(i)%nsize.gt.maxns) maxns=Part(i)%nsize
		if(Part(i)%nT.gt.maxnT) maxnT=Part(i)%nT
	enddo

	if(particledir.eq.' ') particledir=outputdir
	if(particledir(len_trim(particledir)-1:len_trim(particledir)).ne.'/') then
		particledir=trim(particledir) // '/'
	endif
	call output("Particle dir: " // trim(particledir))
	write(command,'("mkdir -p ",a)') trim(particledir)
	call system(command)

	allocate(specemit(nlam))
	allocate(column(npart,maxns,maxnT))
	allocate(KabsTotal(nzones,nlam))
	allocate(KscaTotal(nzones,nlam))
	allocate(i1totalAbs(nzones))
	allocate(i2totalAbs(nzones))
	allocate(i3totalAbs(nzones))
	allocate(i1totalSca(nzones))
	allocate(i2totalSca(nzones))
	allocate(i3totalSca(nzones))
	i1totalAbs=0
	i2totalAbs=0
	i3totalAbs=0
	i1totalSca=0
	i2totalSca=0
	i3totalSca=0

	return
	end

	subroutine ReadZone(key)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(SettingKey) key

	select case(key%key2)
		case("x")
			read(key%value,*) Zone(key%nr1)%x0
		case("y")
			read(key%value,*) Zone(key%nr1)%y0
		case("z")
			read(key%value,*) Zone(key%nr1)%z0
		case("theta")
			read(key%value,*) Zone(key%nr1)%theta0
		case("phi")
			read(key%value,*) Zone(key%nr1)%phi0
		case("xscale")
			read(key%value,*) Zone(key%nr1)%xscale
		case("yscale")
			read(key%value,*) Zone(key%nr1)%yscale
		case("zscale")
			read(key%value,*) Zone(key%nr1)%zscale
		case("nx")
			read(key%value,*) Zone(key%nr1)%nx
		case("ny")
			read(key%value,*) Zone(key%nr1)%ny
		case("nz")
			read(key%value,*) Zone(key%nr1)%nz
		case("nr")
			read(key%value,*) Zone(key%nr1)%nr
		case("nt")
			read(key%value,*) Zone(key%nr1)%nt
		case("np")
			read(key%value,*) Zone(key%nr1)%np
		case("shape")
			Zone(key%nr1)%shape=key%value(1:3)
		case("rin")
			read(key%value,*) Zone(key%nr1)%Rin
		case("rout")
			read(key%value,*) Zone(key%nr1)%Rout
		case("tmax")
			read(key%value,*) Zone(key%nr1)%tmax
		case("rexp")
			read(key%value,*) Zone(key%nr1)%Rexp
		case("iter")
			read(key%value,*) Zone(key%nr1)%iter
		case("denstype")
			Zone(key%nr1)%denstype=trim(key%value)
		case("gamma_exp")
			read(key%value,*) Zone(key%nr1)%gamma_exp
		case("denspow","sigmapow")
			read(key%value,*) Zone(key%nr1)%denspow
		case("mdust")
			read(key%value,*) Zone(key%nr1)%Mdust
		case("alpha")
			read(key%value,*) Zone(key%nr1)%alpha
		case("gas2dust")
			read(key%value,*) Zone(key%nr1)%gas2dust
		case("sh")
			read(key%value,*) Zone(key%nr1)%sh
		case("rsh")
			read(key%value,*) Zone(key%nr1)%Rsh
		case("shpow")
			read(key%value,*) Zone(key%nr1)%shpow
		case("dx")
			read(key%value,*) Zone(key%nr1)%dx
		case("dy")
			read(key%value,*) Zone(key%nr1)%dy
		case("dz")
			read(key%value,*) Zone(key%nr1)%dz
		case("sscale")
			Zone(key%nr1)%sscaletype=trim(key%value)
		case("mscale")
			Zone(key%nr1)%mscaletype=trim(key%value)
		case("amin")
			read(key%value,*) Zone(key%nr1)%amin
		case("amax")
			read(key%value,*) Zone(key%nr1)%amax
		case("apow")
			read(key%value,*) Zone(key%nr1)%apow
		case("adens")
			read(key%value,*) Zone(key%nr1)%Adens
		case("aheight")
			read(key%value,*) Zone(key%nr1)%Aheight
		case("aalpha")
			read(key%value,*) Zone(key%nr1)%Aalpha
		case("rspiral")
			read(key%value,*) Zone(key%nr1)%r_spiral
		case("phispiral")
			read(key%value,*) Zone(key%nr1)%phi_spiral
		case("alphaspiral")
			read(key%value,*) Zone(key%nr1)%alpha_spiral
		case("betaspiral")
			read(key%value,*) Zone(key%nr1)%beta_spiral
		case("wspiral")
			read(key%value,*) Zone(key%nr1)%w_spiral
		case default
			call output("Unknown zone keyword: " // trim(key%key2))
			criticalerror=.true.
	end select


	return
	end
	
	subroutine ReadStar(key)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey) key

	select case(key%key2)
		case("file")
			Star(key%nr1)%file=trim(key%value)
			Star(key%nr1)%startype='FILE'
		case("x")
			read(key%value,*) Star(key%nr1)%x
		case("y")
			read(key%value,*) Star(key%nr1)%y
		case("z")
			read(key%value,*) Star(key%nr1)%z
		case("l")
			read(key%value,*) Star(key%nr1)%L
		case("r")
			read(key%value,*) Star(key%nr1)%R
		case("t")
			read(key%value,*) Star(key%nr1)%T
		case("m")
			read(key%value,*) Star(key%nr1)%M
		case("logg")
			read(key%value,*) Star(key%nr1)%logg
		case("type")
			Star(key%nr1)%startype=trim(key%value)
		case default
			call output("Unknown star keyword: " // trim(key%key2))
			criticalerror=.true.
	end select

	return
	end
	

	subroutine ReadMCobs(key)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey) key

	select case(key%key2)
		case("theta")
			read(key%value,*) MCobs(key%nr1)%theta
		case("phi")
			read(key%value,*) MCobs(key%nr1)%phi
		case("opening")
			read(key%value,*) MCobs(key%nr1)%opening
		case("npix")
			read(key%value,*) MCobs(key%nr1)%npix
		case("maxr")
			read(key%value,*) MCobs(key%nr1)%maxR
		case("raytrace")
			read(key%value,*) MCobs(key%nr1)%raytrace
		case("nphot")
			read(key%value,*) MCobs(key%nr1)%Nphot
		case("lam")
			select case(key%nr2)
				case(1)
					read(key%value,*) MCobs(key%nr1)%lam1
				case(2)
					read(key%value,*) MCobs(key%nr1)%lam2
				case default
					call output("Unknown MCobs:lam value")
			end select
		case default
			call output("Unknown MCobs keyword: " // trim(key%key2))
			criticalerror=.true.
	end select

	return
	end
	

	subroutine ReadComputePart(key)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey) key

	Part(key%nr1)%ptype="COMPUTE"

	select case(key%key2)
		case("file")
			Part(key%nr1)%file=trim(key%value)
		case("ngrains","nsize")
			read(key%value,*) Part(key%nr1)%nsize
		case("nsubgrains")
			read(key%value,*) Part(key%nr1)%nsubgrains
		case("amin")
			read(key%value,*) Part(key%nr1)%amin
		case("amax")
			read(key%value,*) Part(key%nr1)%amax
		case("apow")
			read(key%value,*) Part(key%nr1)%apow
		case("fmax")
			read(key%value,*) Part(key%nr1)%fmax
		case("blend")
			read(key%value,*) Part(key%nr1)%blend
		case("porosity")
			read(key%value,*) Part(key%nr1)%porosity
		case("standard")
			Part(key%nr1)%standard=trim(key%value)
		case("fcarbon")
			read(key%value,*) Part(key%nr1)%fcarbon
		case default
			call output("Unknown computepart keyword: " // trim(key%key2))
			criticalerror=.true.
	end select

	return
	end
	



	subroutine GetKeywords(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer ncla	! number of command line arguments
	character*1000 readline,inputfile,command

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

	call getarg(2,readline)
	read(readline,*) Nphot
	if(Nphot.gt.0) then
		call output("Number of photon packages for radiative transfer: " // int2string(Nphot,'(i10)'))
	else
		call output("No radiative transfer")
	endif

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
	use Constants
	IMPLICIT NONE
	integer i
	
	maxiter=0
	
	lam1=0.1
	lam2=10000
	nlam=200
	zlam1=5
	zlam2=35
	nzlam=0
	distance=150d0
	fstop=0d0
	
	particledir=' '
	
	do i=1,nstars
		Star(i)%x=0d0
		Star(i)%y=0d0
		Star(i)%z=0d0
		Star(i)%L=1d0
		Star(i)%R=1d0
		Star(i)%T=5777d0
		Star(i)%logg=4d0
		Star(i)%startype='KURUCZ'
	enddo
	
	do i=1,nzones
		Zone(i)%nr=150
		Zone(i)%nt=80
		Zone(i)%np=90
		Zone(i)%shape='SPH'
		Zone(i)%x0=0d0
		Zone(i)%y0=0d0
		Zone(i)%z0=0d0
		Zone(i)%theta0=0d0
		Zone(i)%phi0=0d0
		Zone(i)%xscale=1d0
		Zone(i)%yscale=1d0
		Zone(i)%zscale=1d0
		Zone(i)%Rin=1d0
		Zone(i)%Rout=500d0
		Zone(i)%Rexp=100d0
		Zone(i)%iter=.false.
		Zone(i)%denstype='DISK'
		Zone(i)%denspow=1d0
		Zone(i)%gamma_exp=1d0
		Zone(i)%Mdust=1d-4
		Zone(i)%alpha=1d-2
		Zone(i)%gas2dust=100d0
		Zone(i)%sh=0.1
		Zone(i)%Rsh=1.0
		Zone(i)%shpow=1.1
		Zone(i)%tmax=45d0
		Zone(i)%sscaletype='AU'
		Zone(i)%mscaletype='Msun'
		Zone(i)%abun(1:npart)=1d0/real(npart)
		Zone(i)%amin=0.05
		Zone(i)%amax=3000d0
		Zone(i)%apow=3.5

		Zone(i)%Adens=0d0			! Amplitude of wave in density
		Zone(i)%Aheight=0d0			! Amplitude of wave in scaleheight
		Zone(i)%Aalpha=0d0			! Amplitude of wave in alpha
		Zone(i)%r_spiral=5d0		! launching point
		Zone(i)%phi_spiral=0d0		! launching point
		Zone(i)%alpha_spiral=1.5d0	! Kepler rotation
		Zone(i)%beta_spiral=-0.4d0	! Value from Muto = 0.4, setting it to negative computes it
		Zone(i)%w_spiral=0.2d0		! in terms of orbits
	enddo

	do i=1,npart
		Part(i)%standard='DIANA'
		Part(i)%amin=0.05
		Part(i)%amax=3000.0
		Part(i)%apow=3.5
		Part(i)%fmax=0.8
		Part(i)%porosity=0.25
		Part(i)%fcarbon=0.15
		Part(i)%nsize=1
		Part(i)%nT=1
		Part(i)%nsubgrains=1
	enddo
	
	do i=1,nMCobs
		MCobs(i)%theta=35d0
		MCobs(i)%phi=0d0
		MCobs(i)%opening=5d0
		MCobs(i)%npix=500
		MCobs(i)%raytrace=.false.
		MCobs(i)%Nphot=100000
		MCobs(i)%lam1=0.1d0
		MCobs(i)%lam2=3000d0
		MCobs(i)%maxR=-1d0
	enddo
	
	return
	end
	

	subroutine CountStuff(firstkey)
	use GlobalSetup
	IMPLICIT NONE
	type(SettingKey),target :: firstkey
	type(SettingKey),pointer :: key
	integer i

	key => firstkey

	nzones=1
	npart=1
	nstars=1
	nMCobs=1
	do while(.not.key%last)
		select case(key%key1)
			case("zone")
				if(key%nr1.gt.nzones) nzones=key%nr1
			case("star")
				if(key%nr1.gt.nstars) nstars=key%nr1
			case("part","opac","computepart")
				if(key%nr1.gt.npart) npart=key%nr1
			case("mcobs")
				if(key%nr1.gt.nMCobs) nMCobs=key%nr1
		end select
		key=>key%next
	enddo

	call output('Number of zones:        ' // int2string(nzones,'(i4)'))
	call output('Number of stars:        ' // int2string(nstars,'(i4)'))
	call output('Number of particles:    ' // int2string(npart,'(i4)'))
	call output('Number of observations: ' // int2string(npart,'(i4)'))

	allocate(Zone(nzones))
	allocate(Star(nstars))
	allocate(Part(npart))
	allocate(MCobs(nMCobs))

	do i=1,npart
		Part(i)%nT=0
	enddo

	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("computepart")
				if(key%key2.eq."tfile") then
					Part(key%nr1)%nT=Part(key%nr1)%nT+1
				endif
		end select
		key=>key%next
	enddo
	do i=1,npart
		if(Part(i)%nT.eq.0) Part(i)%nT=1
		allocate(Part(i)%file(Part(i)%nT))
		allocate(Part(i)%Tmax(Part(i)%nT))
	enddo

	do i=1,nzones
		allocate(Zone(i)%abun(npart))
	enddo

	return
	end
	
	
	subroutine GetOutputDir
	use GlobalSetup
	IMPLICIT NONE
	integer ncla
	character*500 readline,command

	outputdir='./'

	ncla=3
1	continue
	call getarg(ncla,readline)
	if(readline(1:2).eq.'-o') then
		call getarg(1+ncla,outputdir)
		if(outputdir(len_trim(outputdir)-1:len_trim(outputdir)).ne.'/') then
			outputdir=trim(outputdir) // '/'
		endif
		ncla=ncla+1
		goto 1
	endif
	ncla=ncla+1
	if(readline.ne.' ') goto 1

	call output("Output dir: " // trim(outputdir))
	write(command,'("mkdir -p ",a)') trim(outputdir)
	call system(command)
	
	return
	end


	