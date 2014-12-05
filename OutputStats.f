	subroutine OutputStats
	use GlobalSetup
	real*8 tau,lam0,d,GetKext
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
		tau=0d0
		do ir=1,Zone(izone)%nr
			tau=tau+GetKext(ilam,Zone(izone)%C(ir,Zone(izone)%imidplane,1))*
     &					(Zone(izone)%R(ir+1)-Zone(izone)%R(ir))
		enddo
		call output("Optical depth through zone " // trim(int2string(izone,"(i4)")) // 
     &						":" // trim(dbl2string(tau,"(e14.4)")))
	enddo

	return
	end
	
	
	