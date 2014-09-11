c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec,Lsun,sigma
	real*8 Mearth,Rearth,Mjup,Rjup
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598d13)
	parameter(parsec=3.08568025d18)
	parameter(Rsun=6.955d10)
	parameter(Msun=1.98892d33)
	parameter(Lsun=3.827d33)
	parameter(kb=1.3806503d-16)
	parameter(sigma=5.6704d-5)
	parameter(mp=1.67262178d-24)	!proton mass
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068d-27) ! cm^2 g/s
	parameter(Mearth=5.97219d27)
	parameter(Mjup=1.89813d30)
	parameter(Rearth=6.3781d8)
	parameter(Rjup=7.1492d9)
	
	end module Constants

c=========================================================================================
c global setup for MCMax3D
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
	integer nzones,nstars,npart,maxns,maxnT,nMCobs
	integer maxiter,Nphot,idum
	logical criticalerror
	real*8 maxR,distance
	character*500 outputdir,particledir

c string converting functions
	character*20 int2string,dbl2string
	external int2string,dbl2string

c wavelength grid
	integer nlam,nzlam
	real*8 lam1,lam2,zlam1,zlam2
	real*8,allocatable :: lam(:),nu(:)
	
c Planck functions
	integer nBB
	real*8 dTBB
	parameter(nBB=5000,dTBB=1d0)
	real*8,allocatable :: BB(:,:)		! dimension nlam,nBB
	
c particle scattering
	integer nspike
	real*8 cos2phi(0:360),sin2phi(0:360)

c multiwav thingies
	real*8,allocatable :: specemit(:)				! nlam
	real*8,allocatable :: column(:,:,:)				! npart,nsize,iT
	
c storage speed options
	real*8,allocatable :: KabsTotal(:,:),KscaTotal(:,:)		! nzones,nlam
	integer,allocatable :: i1totalAbs(:),i2totalAbs(:),i3totalAbs(:)	! nzones
	integer,allocatable :: i1totalSca(:),i2totalSca(:),i3totalSca(:)	! nzones

	type Mueller
		real*8 F11(180),F12(180),F22(180)
		real*8 F33(180),F44(180),F34(180)
		real*8 IF11,IF12
	end type Mueller

	type StarType
		real*8 x,y,z
		real*8 L,R,T,logg,M
		character*10 startype
		character*500 file
		real*8,allocatable :: F(:)								! dimension nlam
	end type StarType

	type Cell
		real*8 T,M,V,E,dens		! Temperature, Mass, Volume, Energy absorbed, density
		real*8 Etrace,gasdens
		integer Ni				! statistics
		real*8,allocatable :: densP(:,:,:)						! dimension npart,nsize,nT
		logical diff,randomwalk
	end type Cell
	
	type Photon
		real*8 x,y,z,vx,vy,vz,sI,sQ,sU,sV
		real*8,allocatable :: xzone(:),yzone(:),zzone(:),vxzone(:),vyzone(:),vzzone(:)
		real*8 Sx,Sy,Sz,lam,nu,x0,y0,z0
		integer,allocatable :: i1(:),i2(:),i3(:),edgeNr(:)	! dimension nzones
		logical,allocatable :: inzone(:)					! dimension nzones
		integer ilam1,ilam2,nr
		real*8 wl1,wl2,Kext,Kabs
		real*8,allocatable :: KabsZ(:)						! dimension nzones
		logical scatt,pol
	end type Photon

	type Particle
		real*8,allocatable :: rv(:)									! dimension nsize
		real*8,allocatable :: rho(:)								! dimension nT
		real*8,allocatable :: Tmax(:)								! dimension nT
		real*8 amin,amax,apow,fmax,porosity,fcarbon
		logical blend
		real*8,allocatable :: Kabs(:,:,:),Ksca(:,:,:),Kext(:,:,:)	! dimension nsize,nT,nlam
		type(Mueller),allocatable :: F(:,:,:)						! dimension nsize,nT,nlam
		real*8,allocatable :: Kp(:,:,:)								! dimension nsize,nT,nBB
		character*500,allocatable :: file(:)						! dimension nT
		character*20 standard,ptype
		integer nsize,nT,nsubgrains
		real dust_moment1,dust_moment2,dust_moment3,rvmin,rvmax
	end type Particle
	
	type ZoneType
		type(Cell),allocatable :: C(:,:,:)					! dimension nx,ny,nz
		integer nx,ny,nz
		integer nr,nt,np,n1,n2,n3,imidplane
		real*8 x0,y0,z0,phi0,theta0
		real*8 Rin,Rout,dx,dy,dz,tmax
		character*3 shape			! CAR, SPH, CYL
		character*10 sscaletype,mscaletype
		real*8 sscale,mscale
		character*10 denstype		! DISK, SHELL
		logical iter
		real*8 denspow,Mdust,alpha,Rexp,sh,Rsh,shpow,gamma_exp,gas2dust
		real*8 amin,amax,apow
		real*8,allocatable :: R(:),theta(:),phi(:),x(:),y(:),z(:),abun(:)
		real*8,allocatable :: R2(:),cost2(:),tanx(:),tany(:)
	end type ZoneType

	type MCobsType
		integer npix
		real*8,allocatable :: image(:,:,:)						! dimension npix,npix,nlam
		real*8,allocatable :: spec(:)							! dimension nlam
		real*8 x,y,z,theta,phi,opening
	end type MCobsType
	
	type(ZoneType),allocatable,target :: Zone(:)						! dimension nzones
	type(StarType),allocatable :: Star(:)						! dimension nstars
	type(Particle),allocatable :: Part(:)						! dimension npart
	type(MCobsType),allocatable :: MCobs(:)						! dimension nMCobs

	type SettingKey
		character*100 key1,key2,value
		integer nr1,nr2
		logical last
		type(SettingKey),pointer :: next
	end type SettingKey

	type Travel
		real*8 v
		integer i1next,i2next,i3next,edgenext
	end type Travel
	
	end module GlobalSetup

