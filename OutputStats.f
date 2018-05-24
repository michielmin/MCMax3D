	subroutine OutputStats
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 lam0,d,GetKext,radtau,Mtot,tau,wl1,wl2,Ke
	real*8,allocatable :: rad(:)
	integer ilam,i,izone,ir,ip,it,is

	lam0=lam_ref
	d=lam(nlam)-lam(1)
	ilam=1
	if(lam0.le.lam(1)) then
		wl1=1d0
		wl2=0d0
		ilam=1
	else if(lam0.gt.lam(nlam)) then
		wl1=0d0
		wl2=1d0
		ilam=nlam-1
	else
		do i=1,nlam-1
			if(lam0.gt.lam(i).and.lam0.le.lam(i+1)) then
				ilam=i
			endif
		enddo
		wl1=(lam(ilam+1)-lam0)/(lam(ilam+1)-lam(ilam))
		wl2=1d0-wl1
	endif

	do i=1,npart
		do is=1,Part(i)%nsize
		do iT=1,Part(i)%nT
			Ke=wl1*Part(i)%Kext(is,iT,ilam)+wl2*Part(i)%Kext(is,iT,ilam+1)
			call output("Opacity of particle " // trim(int2string(i,"(i4)")) // ":" //
     &				trim(dbl2string(Ke,"(e14.4)")) //
     &				" g/cm^2  (" // trim(dbl2string(Part(i)%rv(is),"(e9.4)")) // " cm)")
		enddo
		enddo
	enddo
	
	do izone=1,nzones
		Zone(izone)%tau_V=0d0
		do ir=1,Zone(izone)%nr
			Ke=wl1*GetKext(ilam,Zone(izone)%C(ir,Zone(izone)%imidplane,1))
			Ke=Ke+wl2*GetKext(ilam+1,Zone(izone)%C(ir,Zone(izone)%imidplane,1))
			Zone(izone)%tau_V=Zone(izone)%tau_V+Ke*(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
		enddo
		call output("Optical depth through zone " // trim(int2string(izone,"(i4)")) // 
     &						":" // trim(dbl2string(Zone(izone)%tau_V,"(e14.4)")))
c GFORTRAN increased the record length, otherwise gfortran complains
		open(unit=20,file=trim(outputdir) // "heightR" // trim(int2string(izone,'(i0.4)')) // ".dat",RECL=10000)
		allocate(rad(Zone(izone)%np))
		do it=1,Zone(izone)%nt
			do ip=1,Zone(izone)%np
				radtau=0d0
				rad(ip)=0d0
				do ir=1,Zone(izone)%nr
					Ke=wl1*GetKext(ilam,Zone(izone)%C(ir,it,ip))+wl2*GetKext(ilam+1,Zone(izone)%C(ir,it,ip))
					tau=Ke*(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
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
				Ke=wl1*GetKext(ilam,Zone(izone)%C(ir,1,1))+wl2*GetKext(ilam+1,Zone(izone)%C(ir,1,1))
				radtau=radtau+Ke*(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
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
	
	
	
