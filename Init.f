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
	integer i,j,omp_get_max_threads,omp_get_thread_num

	call SetZoneStructOutput()

	j=omp_get_max_threads()
!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(j)
!$OMP DO SCHEDULE(STATIC,j)
	do i=1,j
		idum=-42-omp_get_thread_num()
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL

	
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
		case("spiral")
			call ReadSpiral(key)
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
		case("lam_ref","ref_lam","lamref","reflam")
			read(key%value,*) lam_ref
		case("maxiter")
			read(key%value,*) maxiter
		case("gammauvdes")
			read(key%value,*) gammaUVdes
		case("distance")
			read(key%value,*) distance
		case("dirparticle","particledir")
			particledir=trim(key%value)
		case("fstop")
			read(key%value,*) fstop
		case("nspike")
			read(key%value,*) nspike
		case("av")
			read(key%value,*) Av
		case("adjustav")
			read(key%value,*) adjustAv
		case("abun_in_name")
			read(key%value,*) abun_in_name
		case("multi","openmp")
			read(key%value,*) use_multi
		case("rt_multi","rt_openmp")
			read(key%value,*) rt_multi
		case("multinphot","multinphotmono")
			read(key%value,*) MultiNphotMono
		case("hardedge")
			read(key%value,*) hardedge
		case("opendisk")
			read(key%value,*) opendisk
			opendisk=opendisk*pi/180d0
		case("transtrace")
			read(key%value,*) transmissiontracing
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

	j=omp_get_max_threads()
	if(.not.use_multi) j=1
!$OMP PARALLEL IF(use_multi)
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(j,nlam,npart,maxns,maxnT,nzones)
!$OMP DO SCHEDULE(STATIC,1)
	do i=1,j
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
	enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL


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
			if(key%nr2.le.1) then
				read(key%value,*) Zone(key%nr1)%theta0
			else if(key%nr2.eq.2) then
				read(key%value,*) Zone(key%nr1)%theta1
				Zone(key%nr1)%warped=.true.
			endif
		case("phi")
			if(key%nr2.le.1) then
				read(key%value,*) Zone(key%nr1)%phi0
			else if(key%nr2.eq.2) then
				read(key%value,*) Zone(key%nr1)%phi1
				Zone(key%nr1)%warped=.true.
			endif
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
		case("densfile")
			Zone(key%nr1)%densfile=trim(key%value)
		case("gamma_exp")
			read(key%value,*) Zone(key%nr1)%gamma_exp
		case("denspow","sigmapow")
			read(key%value,*) Zone(key%nr1)%denspow
		case("bp_alpha")
			read(key%value,*) Zone(key%nr1)%bp_alpha
		case("bp_beta")
			read(key%value,*) Zone(key%nr1)%bp_beta
		case("bp_p")
			read(key%value,*) Zone(key%nr1)%bp_p
		case("bp_a")
			read(key%value,*) Zone(key%nr1)%bp_A
		case("bp_eta")
			read(key%value,*) Zone(key%nr1)%bp_eta
		case("mdust")
			read(key%value,*) Zone(key%nr1)%Mdust
		case("tau","maxtau")
			read(key%value,*) Zone(key%nr1)%tau_V
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
		case("fbeam")
			read(key%value,*) Zone(key%nr1)%fbeam
		case("abun")
			read(key%value,*) Zone(key%nr1)%abun(key%nr2)
		case("thin")
			read(key%value,*) Zone(key%nr1)%thin
		case("taufile")
			Zone(key%nr1)%taufile=trim(key%value)
		case("avortex","chivortex")
			read(key%value,*) Zone(key%nr1)%avortex
		case("rvortex")
			read(key%value,*) Zone(key%nr1)%rvortex
		case("phivortex")
			read(key%value,*) Zone(key%nr1)%phivortex
		case("dvortex","deltavortex")
			read(key%value,*) Zone(key%nr1)%dvortex
		case("roundtype")
			read(key%value,*) Zone(key%nr1)%roundtype
		case("roundradius")
			read(key%value,*) Zone(key%nr1)%roundradius
		case("roundindex")
			read(key%value,*) Zone(key%nr1)%roundindex
		case("rounddepth")
			read(key%value,*) Zone(key%nr1)%rounddepth
		case("warp_pow","warppow","powwarp")
			read(key%value,*) Zone(key%nr1)%warp_pow
		case default
			call output("Unknown zone keyword: " // trim(key%key2))
			criticalerror=.true.
	end select


	return
	end

	subroutine ReadSpiral(key)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(SettingKey) key

	select case(key%key2)
		case("adens")
			read(key%value,*) Spiral(key%nr1)%Adens
		case("aheight")
			read(key%value,*) Spiral(key%nr1)%Aheight
		case("aalpha")
			read(key%value,*) Spiral(key%nr1)%Aalpha
		case("r")
			read(key%value,*) Spiral(key%nr1)%r
		case("rin")
			read(key%value,*) Spiral(key%nr1)%Rin
		case("rout")
			read(key%value,*) Spiral(key%nr1)%Rout
		case("hr")
			read(key%value,*) Spiral(key%nr1)%hr
		case("phi")
			read(key%value,*) Spiral(key%nr1)%phi
		case("alpha")
			read(key%value,*) Spiral(key%nr1)%alpha
		case("beta")
			read(key%value,*) Spiral(key%nr1)%beta
		case("w")
			read(key%value,*) Spiral(key%nr1)%w
		case("q")
			read(key%value,*) Spiral(key%nr1)%q
		case("sign")
			read(key%value,*) Spiral(key%nr1)%sign
		case("n")
			read(key%value,*) Spiral(key%nr1)%n
		case("a")
			if(key%nr2.eq.1) read(key%value,*) Spiral(key%nr1)%A1
			if(key%nr2.eq.2) read(key%value,*) Spiral(key%nr1)%A2
		case("type")
			Spiral(key%nr1)%type=trim(key%value)
		case default
			call output("Unknown spiral keyword: " // trim(key%key2))
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
		case("fov")
			read(key%value,*) MCobs(key%nr1)%fov
		case("raytrace")
			read(key%value,*) MCobs(key%nr1)%raytrace
		case("mcout")
			read(key%value,*) MCobs(key%nr1)%mcout
		case("image","writeimage")
			read(key%value,*) MCobs(key%nr1)%writeimage
		case("nphot")
			read(key%value,*) MCobs(key%nr1)%Nphot
		case("nr")
			read(key%value,*) MCobs(key%nr1)%nr
		case("np")
			read(key%value,*) MCobs(key%nr1)%np
		case("lam")
			select case(key%nr2)
				case(1)
					read(key%value,*) MCobs(key%nr1)%lam1
				case(2)
					read(key%value,*) MCobs(key%nr1)%lam2
				case default
					call output("Unknown MCobs:lam value")
			end select
		case("d")
			read(key%value,*) MCobs(key%nr1)%D
		case("d2")
			read(key%value,*) MCobs(key%nr1)%D2
		case("spw")
			read(key%value,*) MCobs(key%nr1)%SpW
		case("width")
			read(key%value,*) MCobs(key%nr1)%width
		case("snoise")
			read(key%value,*) MCobs(key%nr1)%snoise
		case("fstar")
			read(key%value,*) MCobs(key%nr1)%fstar
		case("telescope")
			read(key%value,*) MCobs(key%nr1)%telescope
		case("flag")
			MCobs(key%nr1)%flag=key%value
		case("tracestar")
			read(key%value,*) MCobs(key%nr1)%trace(nzones+key%nr2)
		case("tracezone")
			read(key%value,*) MCobs(key%nr1)%trace(key%nr2)
		case("next")
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
	integer j

	Part(key%nr1)%ptype="COMPUTE"

	select case(key%key2)
		case("file")
			Part(key%nr1)%file(1)=trim(key%value)
		case("tfile")
			do j=1,Part(key%nr1)%nT
				if(key%nr2.eq.int(Part(key%nr1)%Tmax(j))) then
					Part(key%nr1)%file(j)=trim(key%value)
				endif
			enddo
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
		case("tdesa")
			read(key%value,*) Part(key%nr1)%TdesA
		case("tdesb")
			read(key%value,*) Part(key%nr1)%TdesB
		case("abun")
			if(key%nr2.gt.50) stop "number of species in file too large"
			read(key%value,*) Part(key%nr1)%inp_abun(key%nr2)
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
	logical readfile

	call getarg(1,inputfile)

	open(unit=20,file=inputfile,RECL=1000)

	call getarg(2,readline)
	read(readline,*) Nphot
	if(Nphot.gt.0) then
		call output("Number of photon packages for radiative transfer: " // int2string(Nphot,'(i10)'))
	else
		call output("No radiative transfer")
	endif

	call system("cp " // trim(inputfile) // " " // trim(outputdir) // "input.dat")
	if(Nphot.gt.0) then
		open(unit=21,file=trim(outputdir) // "input.dat",RECL=1000,access='APPEND')
		write(21,'("*** command line keywords ***")')
	endif
	
	ncla=0

	key => firstkey

	readfile=.true.
	
	goto 10
20	readfile=.false.
	close(unit=20)
10	continue
	if(readfile) then
		call ignorestar(20)
		read(20,'(a1000)',end=20,err=20) readline
	else
		call getarg(3+ncla,readline)
		if(readline(1:2).eq.'-s') then
			ncla=ncla+1
			call getarg(3+ncla,readline)
			call output("Command line argument: " // trim(readline))
			if(Nphot.gt.0) write(21,'(a)') trim(readline)
			ncla=ncla+1
		else if(readline(1:2).eq.'-o') then
			ncla=ncla+2
			goto 10
		else
			if(readline.ne.' ') then
c				try to read another command line argument
				open(unit=20,file=readline,RECL=1000)
				readfile=.true.
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
	close(unit=21)
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
	if(ikey1.eq.0) ikey1=len_trim(line)+1

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
2	continue
	if(i.eq.n) then
		nr=0
	else
		read(key(i+1:n),*,err=3) nr	
	endif
3	continue
	key=key(1:i)
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
	if(maxnT.gt.1) maxiter=1
	
	lam1=0.1
	lam2=10000
	lam_ref=0.55d0
	nlam=200
	zlam1=5
	zlam2=35
	nzlam=0
	distance=150d0
	fstop=0d0
	nspike=0
	Av=0d0
	adjustAv=.false.
	gammaUVdes=0d0
	use_multi=.true.
	rt_multi=.true.
	transmissiontracing=.false.

	delta_St=1d0
	
	particledir=' '
	
	MultiNphotMono=1
	hardedge=.false.
	opendisk=-1d0
	
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
		Zone(i)%theta1=-500d0
		Zone(i)%phi1=-500d0
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
		Zone(i)%tau_V=-1d0
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
		zone(i)%thin=.false.
		Zone(i)%fbeam=0d0
		Zone(i)%reflect=.false.
		Zone(i)%avortex=4d0
		Zone(i)%rvortex=-1d0
		Zone(i)%phivortex=0d0
		Zone(i)%dvortex=1d-2
		Zone(i)%roundtype='NONE'
		Zone(i)%roundradius=50d0
		Zone(i)%roundindex=0.2d0
		Zone(i)%rounddepth=1d-20
		Zone(i)%warped=.false.
		Zone(i)%warp_pow=1d0
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
		Part(i)%nsubgrains=1
		Part(i)%inp_abun=-1d0
!c Default values are those of water ice		!	olivine
		Part(i)%TdesA=27.303					!	4.4491649
		Part(i)%TdesB=3.536						!	0.35676055
	enddo
	
	do i=1,nMCobs
		MCobs(i)%theta=35d0
		MCobs(i)%phi=0d0
		MCobs(i)%opening=5d0
		MCobs(i)%npix=500
		MCobs(i)%raytrace=.false.
		MCobs(i)%mcout=.true.
		MCobs(i)%writeimage=.true.
		MCobs(i)%Nphot=100000
		MCobs(i)%lam1=0.1d0
		MCobs(i)%lam2=3000d0
		MCobs(i)%maxR=-1d0
		MCobs(i)%fov=-1d0
		MCobs(i)%nr=1
		MCobs(i)%np=360
		MCobs(i)%telescope=.false.
		MCobs(i)%D=0d0
		MCobs(i)%D2=0d0
		MCobs(i)%SpW=0d0
		MCobs(i)%width=0d0
		MCobs(i)%snoise=0d0
		MCobs(i)%fstar=1d0
		MCobs(i)%flag=' '
	enddo

	do i=1,nSpirals
		Spiral(i)%Adens=0d0			! Amplitude of wave in density
		Spiral(i)%Aheight=0d0		! Amplitude of wave in scaleheight
		Spiral(i)%Aalpha=0d0		! Amplitude of wave in alpha
		Spiral(i)%r=5d0				! launching point
		Spiral(i)%hr=-1d0			! H/R at the launching point, negative means it is computed
		Spiral(i)%phi=0d0			! launching point
		Spiral(i)%alpha=1.5d0		! Kepler rotation
		Spiral(i)%beta=0.4d0		! Value from Muto = 0.4, setting it to negative computes it
		Spiral(i)%w=5d0				! in AU
		Spiral(i)%sign=1d0			! which way it rotates
		Spiral(i)%q=1.7d0			! powerlaw dependence of the amplitude
		Spiral(i)%Rin=0d0			! inner radius of the spiral
		Spiral(i)%Rout=1d200		! outer radius of the spiral
c for an Archimedean spiral:
c R = Spiral(i)%A1 + Spiral(i)%A2 * (phi-Spiral(i)%phi)^Spiral(i)%n
		Spiral(i)%A1=5d0
		Spiral(i)%A2=1d0
		Spiral(i)%n=1d0				! for Archimedean spiral
		Spiral(i)%type='PLANET'		! 'PLANET','ARCH'
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
	nMCobs=0
	nSpirals=0
	do while(.not.key%last)
		select case(key%key1)
			case("zone")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nzones) nzones=key%nr1
			case("star")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nstars) nstars=key%nr1
			case("part","opac","computepart")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.npart) npart=key%nr1
			case("mcobs")
				if(key%key2.eq.'next') nMCobs=nMCobs+1
				if(key%nr1.eq.0) key%nr1=nMCobs
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nMCobs) nMCobs=key%nr1
			case("spiral")
				if(key%nr1.eq.0) key%nr1=1
				if(key%nr2.eq.0) key%nr2=1
				if(key%nr1.gt.nSpirals) nSpirals=key%nr1
		end select
		key=>key%next
	enddo

	if(nMCobs.eq.0) nMCobs=1

	call output('Number of zones:        ' // int2string(nzones,'(i4)'))
	call output('Number of stars:        ' // int2string(nstars,'(i4)'))
	call output('Number of particles:    ' // int2string(npart,'(i4)'))
	call output('Number of spirals:      ' // int2string(nSpirals,'(i4)'))
	call output('Number of observations: ' // int2string(nMCobs,'(i4)'))

	allocate(Zone(nzones))
	allocate(Star(nstars))
	allocate(Part(npart))
	allocate(Spiral(max(nSpirals,1)))
	allocate(MCobs(max(nMCobs,1)))

	do i=1,nMCobs
		allocate(MCobs(i)%trace(nzones+nstars))
		MCobs(i)%trace=.true.
	enddo

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

	do i=1,npart
		Part(i)%nT=0
	enddo
	key => firstkey
	do while(.not.key%last)
		select case(key%key1)
			case("computepart")
				if(key%key2.eq."tfile") then
					Part(key%nr1)%nT=Part(key%nr1)%nT+1
					Part(key%nr1)%Tmax(Part(key%nr1)%nT)=key%nr2
				endif
		end select
		key=>key%next
	enddo
	maxnT=1
	do i=1,npart
		if(Part(i)%nT.eq.0) then
			Part(i)%nT=1
		else
			call sort(Part(i)%Tmax(1:Part(i)%nT),Part(i)%nT)
		endif
		if(Part(i)%nT.gt.maxnT) maxnT=Part(i)%nT
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


	