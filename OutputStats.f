	subroutine OutputStats
	use GlobalSetup
	real*8 lam0,d,GetKext,radtau
	integer ilam,i,izone,ir

	lam0=0.55
	d=lam(nlam)-lam(1)
	ilam=1
	do i=1,nlam
		if(abs(lam0-lam(i)).lt.d) then
			d=abs(lam0-lam(i))
			ilam=i
		endif
	enddo
	
	do izone=1,nzones
		Zone(izone)%tau_V=0d0
		do ir=1,Zone(izone)%nr
			Zone(izone)%tau_V=Zone(izone)%tau_V+GetKext(ilam,Zone(izone)%C(ir,Zone(izone)%imidplane,1))*
     &					(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
		enddo
		call output("Optical depth through zone " // trim(int2string(izone,"(i4)")) // 
     &						":" // trim(dbl2string(Zone(izone)%tau_V,"(e14.4)")))
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
		


	return
	end
	
	
	