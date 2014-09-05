	subroutine RadiativeTransfer
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	
c	call InitRadiativeTransfer
	
	do i=1,Nphot
c		call EmitPhoton(phot)
c		call TravelPhoton(phot)
c		call MCoutput(phot)
	enddo

c	call DetermineTemperatures
	
	return
	end
	
	
	real*8 function GetKabs(izone,ilam,i1,i2,i3)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,n1,n2,n3,ilam,ipart,iT,isize
	
	select case(Zone(izone)%shape)
		case("SPH","CYL")
			n1=Zone(izone)%nr
			n2=Zone(izone)%nt
			n3=Zone(izone)%np
		case("CAR")
			n1=Zone(izone)%nx
			n2=Zone(izone)%ny
			n3=Zone(izone)%nz
	end select
	GetKabs=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKabs=GetKabs+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Kabs(isize,iT,ilam)
			enddo
		enddo
	enddo
	GetKabs=GetKabs/Zone(izone)%C(i1,i2,i3)%dens
	
	return
	end
	

	
	real*8 function GetKext(izone,ilam,i1,i2,i3)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,n1,n2,n3,ilam,ipart,iT,isize
	
	select case(Zone(izone)%shape)
		case("SPH","CYL")
			n1=Zone(izone)%nr
			n2=Zone(izone)%nt
			n3=Zone(izone)%np
		case("CAR")
			n1=Zone(izone)%nx
			n2=Zone(izone)%ny
			n3=Zone(izone)%nz
	end select
	GetKext=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKext=GetKext+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Kext(isize,iT,ilam)
			enddo
		enddo
	enddo
	GetKext=GetKext/Zone(izone)%C(i1,i2,i3)%dens
	
	return
	end
	
	
	real*8 function GetKsca(izone,ilam,i1,i2,i3)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,n1,n2,n3,ilam,ipart,iT,isize
	
	select case(Zone(izone)%shape)
		case("SPH","CYL")
			n1=Zone(izone)%nr
			n2=Zone(izone)%nt
			n3=Zone(izone)%np
		case("CAR")
			n1=Zone(izone)%nx
			n2=Zone(izone)%ny
			n3=Zone(izone)%nz
	end select
	GetKsca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				GetKsca=GetKsca+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
			enddo
		enddo
	enddo
	GetKsca=GetKsca/Zone(izone)%C(i1,i2,i3)%dens
	
	return
	end
	
	
	
	subroutine function SetFmatrix(izone,ilam,i1,i2,i3,F,Ksca)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,n1,n2,n3,ilam,ipart,iT,isize
	real*8 Ksca
	type(Mueller) F
	
	select case(Zone(izone)%shape)
		case("SPH","CYL")
			n1=Zone(izone)%nr
			n2=Zone(izone)%nt
			n3=Zone(izone)%np
		case("CAR")
			n1=Zone(izone)%nx
			n2=Zone(izone)%ny
			n3=Zone(izone)%nz
	end select
	Ksca=0d0
	do ipart=1,npart
		do isize=1,Part(ipart)%nsize
			do iT=1,Part(ipart)%nT
				Ksca=Ksca+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)
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
				F%F11(1:180)=F%F11(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F11(1:180)
				F%F12(1:180)=F%F12(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F12(1:180)
				F%F22(1:180)=F%F22(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F22(1:180)
				F%F33(1:180)=F%F33(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F33(1:180)
				F%F34(1:180)=F%F34(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F34(1:180)
				F%F44(1:180)=F%F44(1:180)+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F44(1:180)
				F%IF11=F%IF11+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
     &						Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF11
				F%IF12=F%IF12+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*
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

	Ksca=Ksca/Zone(izone)%C(i1,i2,i3)%dens
	
	return
	end
	