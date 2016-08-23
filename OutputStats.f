	subroutine OutputStats
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 lam0,d,GetKext,radtau,Mtot,tau
	real*8,allocatable :: rad(:)
	integer ilam,i,izone,ir,ip,it,is

	lam0=0.55
	d=lam(nlam)-lam(1)
	ilam=1
	do i=1,nlam
		if(abs(lam0-lam(i)).lt.d) then
			d=abs(lam0-lam(i))
			ilam=i
		endif
	enddo

	do i=1,npart
		do is=1,Part(i)%nsize
		do iT=1,Part(i)%nT
			call output("Opacity of particle " // trim(int2string(i,"(i4)")) // ":" //
     &				trim(dbl2string(Part(i)%Kext(is,iT,ilam),"(e14.4)")) //
     &				" g/cm^2  (" // trim(dbl2string(Part(i)%rv(is),"(e9.4)")) // " cm)")
		enddo
		enddo
	enddo
	
	do izone=1,nzones
		Zone(izone)%tau_V=0d0
		do ir=1,Zone(izone)%nr
			Zone(izone)%tau_V=Zone(izone)%tau_V+GetKext(ilam,Zone(izone)%C(ir,Zone(izone)%imidplane,1))*
     &					(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
		enddo
		call output("Optical depth through zone " // trim(int2string(izone,"(i4)")) // 
     &						":" // trim(dbl2string(Zone(izone)%tau_V,"(e14.4)")))

		open(unit=20,file=trim(outputdir) // "heightR" // trim(int2string(izone,'(i0.4)')) // ".dat",RECL=6000)
		allocate(rad(Zone(izone)%np))
		do it=1,Zone(izone)%nt
			do ip=1,Zone(izone)%np
				radtau=0d0
				rad(ip)=0d0
				do ir=1,Zone(izone)%nr
					tau=GetKext(ilam,Zone(izone)%C(ir,it,ip))*
     &					(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
					if(radtau+tau.gt.1d0) then
						rad(ip)=Zone(izone)%R(ir)+(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))*(1d0-radtau)/tau
						exit
					endif
					radtau=radtau+tau
				enddo
			enddo
			write(20,*) (rad(ip)*sin(Zone(izone)%theta(it))/AU,rad(ip)*cos(Zone(izone)%theta(it))/AU,ip=1,min(Zone(izone)%np,125))
		enddo
		deallocate(rad)
		close(unit=20)
	enddo

	if(adjustAv) then
		radtau=0d0
		do izone=1,nzones
			do ir=1,Zone(izone)%nr
				radtau=radtau+GetKext(ilam,Zone(izone)%C(ir,1,1))*
     &					(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
			enddo
		enddo
		Av=Av+2.5*log10(exp(-radtau))
		if(Av.lt.0d0) Av=0d0
		call output("Av adjusted to: " // dbl2string(Av,'(f5.2)'))
	endif
		
	Mtot=0d0
	do izone=1,nzones
		do ir=1,Zone(izone)%nr
		do it=1,Zone(izone)%nt
		do ip=1,Zone(izone)%np
			Mtot=Mtot+Zone(izone)%C(ir,it,ip)%V*Zone(izone)%C(ir,it,ip)%dens
		enddo
		enddo
		enddo
	enddo
	call output("Total dust mass:" // dbl2string(Mtot/Msun,  '(e14.4)') // " Msun")
	call output("                " // dbl2string(Mtot/Mearth,'(e14.4)') // " Mearth")
	
	return
	end
	
	
	