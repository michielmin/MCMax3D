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

	end module GlobalSetup

