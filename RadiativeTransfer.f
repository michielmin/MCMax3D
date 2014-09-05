	subroutine RadiativeTransfer
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	
	call InitRadiativeTransfer
	
	do i=1,Nphot
c		call EmitPhoton(phot)
c		call TravelPhoton(phot)
c		call MCoutput(phot)
	enddo

c	call DetermineTemperatures
	
	return
	end
	
	
	subroutine InitRadiativeTransfer
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3,n1,n2,n3,ilam,ipart,iT,isize
	
	do izone=1,nzones
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
		do i1=1,n1
			do i2=1,n2
				do i3=1,n3
					Zone(izone)%C(i1,i2,i3)%Kabs=0d0
					do ipart=1,npart
						do isize=1,Part(ipart)%nsize
							do iT=1,Part(ipart)%nT
								Zone(izone)%C(i1,i2,i3)%Kabs(1:nlam)=Zone(izone)%C(i1,i2,i3)%Kabs(1:nlam)
     &									  +Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Kabs(isize,iT,1:nlam)
								Zone(izone)%C(i1,i2,i3)%Ksca(1:nlam)=Zone(izone)%C(i1,i2,i3)%Ksca(1:nlam)
     &									  +Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,1:nlam)
								Zone(izone)%C(i1,i2,i3)%Kext(1:nlam)=Zone(izone)%C(i1,i2,i3)%Kext(1:nlam)
     &									  +Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Kext(isize,iT,1:nlam)
							enddo
						enddo
					enddo
					do ilam=1,nlam
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F12(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F22(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F33(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F34(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F44(1:180)=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%IF11=0d0
						Zone(izone)%C(i1,i2,i3)%F(ilam)%IF12=0d0
						do ipart=1,npart
							do isize=1,Part(ipart)%nsize
								do iT=1,Part(ipart)%nT
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F11(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F12(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F12(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F22(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F22(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F33(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F33(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F34(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F34(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%F11(1:180)=Zone(izone)%C(i1,i2,i3)%F(ilam)%F44(1:180)
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%F44(1:180)
									Zone(izone)%C(i1,i2,i3)%F(ilam)%IF11=Zone(izone)%C(i1,i2,i3)%F(ilam)%IF11
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF11
									Zone(izone)%C(i1,i2,i3)%F(ilam)%IF12=Zone(izone)%C(i1,i2,i3)%F(ilam)%IF12
     &				+Zone(izone)%C(i1,i2,i3)%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,ilam)*Part(ipart)%F(isize,iT,ilam)%IF12
								enddo
							enddo
						enddo
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F11=Zone(izone)%C(i1,i2,i3)%F(ilam)%F11/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F12=Zone(izone)%C(i1,i2,i3)%F(ilam)%F12/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F22=Zone(izone)%C(i1,i2,i3)%F(ilam)%F22/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F33=Zone(izone)%C(i1,i2,i3)%F(ilam)%F33/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F34=Zone(izone)%C(i1,i2,i3)%F(ilam)%F34/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%F44=Zone(izone)%C(i1,i2,i3)%F(ilam)%F44/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%IF11=Zone(izone)%C(i1,i2,i3)%F(ilam)%IF11/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
						Zone(izone)%C(i1,i2,i3)%F(ilam)%IF12=Zone(izone)%C(i1,i2,i3)%F(ilam)%IF12/Zone(izone)%C(i1,i2,i3)%Ksca(ilam)
					enddo
				enddo
			enddo
		enddo
	enddo	
	
	return
	end
	