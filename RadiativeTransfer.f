	subroutine RadiativeTransfer
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	character*500 fitsfile
	
	call InitRadiativeTransfer
	
	call output("==================================================================")
	call output("Emitting " // trim(int2string(Nphot,'(i10)')) // " photon packages")
	do i=1,Nphot
		call tellertje(i,Nphot)
		call EmitPhoton(phot)
c		call TravelPhoton(phot)
		call MCoutput(phot)
	enddo

	call output("==================================================================")
	call output("Writing Monte Carlo observables")
	do i=1,nMCobs
		fitsfile="MCout" // trim(int2string(i,'(i0.4)')) // ".fits"
		call writefitsfile(fitsfile,MCobs(i)%image,nlam,MCobs(i)%npix)
	enddo
c	call DetermineTemperatures
	
	return
	end


	subroutine MCoutput(phot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 inp,tot,x,y,z
	integer i,j,ix,iy
	type(Photon) phot
	
	do i=1,nMCobs
		inp=MCobs(i)%x*phot%vx+MCobs(i)%y*phot%vy+MCobs(i)%z*phot%vz
		if(inp.gt.MCobs(i)%opening) then
			do j=1,nlam
				if(column(j).lt.1000d0) then
					specemit(j)=specemit(j)*exp(-column(j))
				else
					specemit(j)=0d0
				endif
			enddo
			call integrate(specemit,tot)
			specemit=specemit/tot
			specemit=1d23*phot%sI*specemit
			MCobs(i)%spec(1:nlam)=MCobs(i)%spec(1:nlam)+specemit(1:nlam)
			x=phot%x
			y=phot%y
			z=phot%z
			call rotateZ(x,y,z,cos(MCobs(i)%phi),sin(MCobs(i)%phi))
			call rotateY(x,y,z,cos(MCobs(i)%theta),sin(MCobs(i)%theta))
			ix=real(MCobs(i)%npix)*(x+maxR)/(2d0*maxR)
			iy=real(MCobs(i)%npix)*(y+maxR)/(2d0*maxR)
			print*,ix,iy
			MCobs(i)%image(ix,iy,1:nlam)=MCobs(i)%image(ix,iy,1:nlam)+specemit(1:nlam)
		endif
	enddo
	
	return
	end
	

	subroutine EmitPhoton(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	real*8 Ltot,r,random
	
	Ltot=0d0
	do i=1,nstars
		Ltot=Ltot+Star(i)%L
	enddo
	r=Ltot*random(idum)
	do i=1,nstars
		r=r-Star(i)%L
		if(r.lt.0d0) exit
	enddo

	call randomdirection(phot%x,phot%y,phot%z)

	phot%x=Star(i)%R*phot%x
	phot%y=Star(i)%R*phot%y
	phot%z=Star(i)%R*phot%z

	call emit(phot,Star(i)%F,Star(i)%L)

	phot%sI=Ltot/real(Nphot)
	phot%sQ=0d0
	phot%sU=0d0
	phot%sV=0d0

	call randomdirection(phot%vx,phot%vy,phot%vz)
	if(phot%x*phot%vx+phot%y*phot%vy+phot%z*phot%vz.lt.0d0) then
		phot%vx=-phot%vx
		phot%vy=-phot%vy
		phot%vz=-phot%vz
	endif
	
	phot%x=phot%x+Star(i)%x
	phot%y=phot%y+Star(i)%y
	phot%z=phot%z+Star(i)%z
	
	phot%edgeNr=0
	
	return
	end
	


	subroutine InitRadiativetransfer
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,ilam,i
	type(Cell),pointer :: C
	real*8 GetKext

	call output("Initializing radiative transfer")
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
			do i2=1,Zone(izone)%n2
				do i3=1,Zone(izone)%n3
					C => Zone(izone)%C(i1,i2,i3)
					C%E=0d0
					C%Ni=0d0
					C%T=0d0
				enddo
			enddo
		enddo
	enddo
	
	do i=1,nMCobs
		MCobs(i)%spec=0d0
		MCobs(i)%image=0d0
	enddo
	
	return
	end


c This function gives the Planck mean opacity of the cell
	real*8 function GetKp(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C
	
	GetKp=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKp=GetKp+C%densP(ipart,isize,iT)*Part(ipart)%Kp(isize,iT,ilam)
			enddo
		enddo
	enddo
c multiply by the volume to get the total mass
	GetKp=GetKp*C%V
	
	return
	end
	

c This function gives the absorption cross section per cm
	real*8 function GetKabs(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C
	
	GetKabs=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKabs=GetKabs+C%densP(ipart,isize,iT)*Part(ipart)%Kabs(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	

c This function gives the extinction cross section per cm
	real*8 function GetKext(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer ilam,ipart,iT,isize
	type(Cell) C

	GetKext=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKext=GetKext+C%densP(ipart,isize,iT)*Part(ipart)%Kext(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	
	
c This function gives the scattering cross section per cm
	real*8 function GetKsca(ilam,C)
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,ilam,ipart,iT,isize
	type(Cell) C

	GetKsca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKsca=GetKsca+C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
			enddo
		enddo
	enddo
	
	return
	end
	
	
	
c This subroutine gets the scattering matrix in a certain cell
	subroutine function SetFmatrix(ilam,C,F,Ksca)
	use GlobalSetup
	IMPLICIT NONE
	integer i1,i2,i3,ilam,ipart,iT,isize
	real*8 Ksca
	type(Mueller) F
	type(Cell) C
	
	Ksca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				Ksca=Ksca+C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
			enddo
		enddo
	enddo
	F%F11(1:180)=0d0
	F%F12(1:180)=0d0
	F%F22(1:180)=0d0
	F%F33(1:180)=0d0
	F%F34(1:180)=0d0
	F%F44(1:180)=0d0
	F%IF11=0d0
	F%IF12=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				F%F11(1:180)=F%F11(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F11(1:180)
				F%F12(1:180)=F%F12(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F12(1:180)
				F%F22(1:180)=F%F22(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F22(1:180)
				F%F33(1:180)=F%F33(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F33(1:180)
				F%F34(1:180)=F%F34(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F34(1:180)
				F%F44(1:180)=F%F44(1:180)+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F44(1:180)
				F%IF11=F%IF11+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF11
				F%IF12=F%IF12+C%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF12
			enddo
		enddo
	enddo
	F%F11=F%F11/Ksca
	F%F12=F%F12/Ksca
	F%F22=F%F22/Ksca
	F%F33=F%F33/Ksca
	F%F34=F%F34/Ksca
	F%F44=F%F44/Ksca
	F%IF11=F%IF11/Ksca
	F%IF12=F%IF12/Ksca

	Ksca=Ksca/C%dens
	
	return
	end
	