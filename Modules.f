c=========================================================================================
c module containing the physical constants in cgs
c=========================================================================================
	module Constants
	IMPLICIT NONE
	real*8 pi,G,Msun,AU,clight,Rsun,mp,kb,hplanck,parsec,Lsun,sigma
	parameter(pi=3.14159265358979323846264338328d0)
	parameter(clight=2.9979245800d10) !cm/s
	parameter(AU=1.49598e13)
	parameter(parsec=3.08568025e18)
	parameter(Rsun=6.955e10)
	parameter(Msun=1.98892e33)
	parameter(Lsun=3.827e33)
	parameter(kb=1.3806503d-16)
	parameter(sigma=5.6704d-5)
	parameter(mp=1.67262178d-24)	!proton mass
	parameter(G=6.67300d-8) ! in cm^3/g/s^2
	parameter(hplanck=6.626068e-27) ! cm^2 g/s
	
	end module Constants

c=========================================================================================
c global setup for MCMax3D
c=========================================================================================
	module GlobalSetup
	IMPLICIT NONE
	integer nzones,nstars,npart
	integer maxiter,Nphot
	logical criticalerror
	character*500 outputdir,particledir

c string converting functions
	character*20 int2string,dbl2string
	external int2string,dbl2string

c wavelength grid
	integer nlam,nzlam
	real*8 lam1,lam2,zlam1,zlam2
	real*8,allocatable :: lam(:)
	
	type StarType
		real*8 x,y,z
		real*8 L,R,T
		character*10 startype
		real*8,allocatable :: F(:)								! dimension nlam
	end type StarType

	type Cell
		real*8 T,M,V,E			! Temperature, Mass, Volume, Energy absorbed
		integer Ni				! statistics
		real*8,allocatable :: Kabs(:),Ksca(:),Kext(:)			! dimension nlam
		real*8,allocatable :: dens(:)							! dimension npart
		real*8,allocatable :: LRFI(:),LRFQ(:),LRFU(:),LRFV(:)	! dimension nlam
		logical diff,randomwalk
		real*8 x1,x2,y1,y2,z1,z2	! cell edges
		real*8 r1,r2,t1,t2,p1,p2	! cell edges
	end type Cell
	
	type Photon
		real*8 x,y,z,vx,vy,vz,sI,sQ,sU,sV
		real*8 Sx,Sy,Sz
		integer,allocatable :: i(:),j(:),k(:)				! dimension nzones
		logical,allocatable :: inzone(:)					! dimension nzones
		integer ilam1,ilam2,edgeNr
		real*8 wl1,wl2
		logical onEdge,scatt,pol
	end type Photon
	
	type Mueller
		real*8 F11(180),F12(180),F22(180)
		real*8 F33(180),F44(180),F34(180)
		real*8 IF11,IF12
	end type Mueller

	type Particle
		real*8,allocatable :: rv(:,:)								! dimension nsize,nT
		real*8,allocatable :: Tmax(:)								! dimension nT
		real*8 rho,amin,amax,apow,fmax,porosity,fcarbon
		logical blend
		real*8,allocatable :: Kabs(:,:,:),Ksca(:,:,:),Kext(:,:,:)	! dimension nsize,nT,nlam
		type(Mueller),allocatable :: F(:,:,:)						! dimension nsize,nT,nlam
		real*8,allocatable :: Kp(:,:,:)								! dimension nsize,nT,nBB
		character*500,allocatable :: file(:)						! dimension nT
		character*20 standard,ptype
		integer nsize,nT
	end type Particle
	
	type ZoneType
		type(Cell),allocatable :: C(:,:,:)					! dimension nx,ny,nz
		integer nx,ny,nz
		integer nr,nt,np
		real*8 x0,y0,z0,xn,yn,zn
		character*3 shape			! CAR, SPH, CYL
	end type ZoneType
	
	type(ZoneType),allocatable :: Zone(:)						! dimension nzones
	type(StarType),allocatable :: Star(:)						! dimension nstars
	type(Particle),allocatable :: Part(:)						! dimension npart

	type SettingKey
		character*100 key1,key2,value
		integer nr1,nr2
		logical last
		type(SettingKey),pointer :: next
	end type SettingKey

	end module GlobalSetup

