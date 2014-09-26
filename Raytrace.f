	subroutine Raytrace(iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer iobs,i,ilam,ilam0
	real*8 dlmin
	logical simpleobs

	simpleobs=.false.
	do ilam=1,nlam
		if(lam(ilam).ge.MCobs(iobs)%lam1.and.lam(ilam).le.MCobs(iobs)%lam2) then
			simpleobs=.true.
			exit
		endif
	enddo
	ilam0=0
	if(.not.simpleobs) then
		dlmin=lam(nlam)-lam(1)
		do ilam=1,ilam
			if(abs(lam(ilam)-(lam1+lam2)/2d0).lt.dlmin) then
				dlmin=abs(lam(ilam)-(lam1+lam2)/2d0)
				ilam0=ilam
			endif
		enddo
	endif

	do ilam=1,nlam
		if((lam(ilam).ge.MCobs(iobs)%lam1.and.lam(ilam).le.MCobs(iobs)%lam2)
     &			.or.ilam.eq.ilam0) then
			call TraceScattField(iobs,ilam)
			call FormalSolution(iobs,ilam)
		endif
	enddo
	
	return
	end
	
	subroutine TraceScattField(iobs,ilam)
	use GlobalSetup
	IMPLICIT NONE
	integer iobs,ilam,izone,iT,i,i1,i2,i3,iphot,istar
	real*8 GetKabs,Etot,Erandom,random
	logical emitfromstar
	type(Cell),pointer :: C
	type(Photon) phot
	
	Etot=0d0
	do i=1,nstars
		Etot=Etot+Star(i)%F(ilam)
	enddo
	
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
		do i2=1,Zone(izone)%n2
		do i3=1,Zone(izone)%n3
			C => Zone(izone)%C(i1,i2,i3)
			C%Escatt=0d0
			iT=(C%T+0.5d0)/dTBB
			if(iT.lt.1) iT=1
			if(iT.gt.nBB) iT=nBB
			C%KabsL=GetKabs(ilam,C)
			C%Elam=C%KabsL*BB(ilam,iT)*C%V
			Etot=Etot+C%Elam
		enddo
		enddo
		enddo
	enddo

	do iphot=1,MCobs(iobs)%Nphot
		phot%sI=Etot/real(MCobs(iobs)%Nphot)
		phot%ilam1=ilam
2		Erandom=Etot*random(idum)
		emitfromstar=.false.
		do istar=1,nstars
			Erandom=Erandom-Star(i)%F(ilam)
			if(Erandom.lt.0d0) then
				emitfromstar=.true.
				goto 1
			endif
		enddo
		do izone=1,nzones
			do i1=1,Zone(izone)%n1
			do i2=1,Zone(izone)%n2
			do i3=1,Zone(izone)%n3
				Erandom=Erandom-Zone(izone)%C(i1,i2,i3)%Elam
				if(Erandom.lt.0d0) goto 1
			enddo
			enddo
			enddo
		enddo
		call output("something is wrong...")
		goto 2
1		continue
		if(emitfromstar) then
			call EmitPhotonStar(istar)
		else
			call EmitPhotonMatter(izone,i1,i2,i3)
		endif
		call TravelPhotonMono(phot,iobs)
	enddo
	
	return
	end


	subroutine EmitPhotonStar(istar)
	use GlobalSetup
	IMPLICIT NONE
	integer istar
	
	return
	end
	
	
	subroutine EmitPhotonMatter(izone,i1,i2,i3)
	use GlobalSetup
	IMPLICIT NONE
	integer izone,i1,i2,i3	
	
	return
	end
	
	
	

	subroutine TravelPhotonMono(phot,iobs)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,imin,iobs
	logical leave
	real*8 minv,tau0,tau,GetKext,random,GetKabs,fstopmono,albedo,theta
	type(Travel) Trac(nzones)
	type(Cell),pointer :: C
	type(Photon) phot

	tau0=-log(random(idum))
	
	theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
	phot%iscat=180d0*theta/pi
	if(phot%iscat.lt.1) phot%iscat=1
	if(phot%iscat.gt.180) phot%iscat=180
	
1	continue

	phot%Kext=0d0
	phot%Kabs=0d0
	do izone=1,nzones
		phot%KabsZ(izone)=0d0
		if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			phot%Kext=phot%Kext+GetKext(phot%ilam1,C)
			phot%KabsZ(izone)=GetKabs(phot%ilam1,C)
			phot%Kabs=phot%Kabs+phot%KabsZ(izone)
		endif
	enddo

	albedo=(phot%Kext-phot%Kabs)/phot%Kext
	fstopmono=1d0-albedo**0.25

	do izone=1,nzones
		if(phot%inzone(izone)) then
			select case(Zone(izone)%shape)
				case("SPH")
					call TravelSph(phot,izone,Trac(izone))
			end select
		else
			select case(Zone(izone)%shape)
				case("SPH")
					call HitSph(phot,izone,Trac(izone))
			end select
		endif
	enddo
	minv=20d0*maxR
	leave=.true.
	do izone=1,nzones
		if(Trac(izone)%v.gt.0d0.and.Trac(izone)%v.lt.minv) then
			minv=Trac(izone)%v
			imin=izone
			leave=.false.
		endif
	enddo

	if(leave) goto 3
	
	if((tau0-phot%Kext*minv).lt.0d0) then
		minv=tau0/phot%Kext
		phot%edgeNr=0
		call TravelPhotonX(phot,minv)
		call AddEtraceMono(phot,minv)
		call InteractMono(phot)
		theta=acos(MCobs(iobs)%x*phot%vx+MCobs(iobs)%y*phot%vy+MCobs(iobs)%z*phot%vz)
		phot%iscat=180d0*theta/pi
		if(phot%iscat.lt.1) phot%iscat=1
		if(phot%iscat.gt.180) phot%iscat=180
		tau0=-log(random(idum))
		if(random(idum).lt.fstopmono) then
			goto 3
		endif
		phot%sI=phot%sI*albedo/(1d0-fstopmono)
		phot%sQ=phot%sQ*albedo/(1d0-fstopmono)
		phot%sU=phot%sU*albedo/(1d0-fstopmono)
		phot%sV=phot%sV*albedo/(1d0-fstopmono)
		goto 1
	endif

	call TravelPhotonX(phot,minv)
	call AddEtraceMono(phot,minv)
	tau0=tau0-phot%Kext*minv

	do izone=1,nzones
		if(Trac(izone)%v.le.minv.or.izone.eq.imin) then
			phot%i1(izone)=Trac(izone)%i1next
			phot%i2(izone)=Trac(izone)%i2next
			phot%i3(izone)=Trac(izone)%i3next

			if(phot%i1(izone).eq.-1) call determine_i1(phot,izone)
			if(phot%i2(izone).eq.-1) call determine_i2(phot,izone)
			if(phot%i3(izone).eq.-1) call determine_i3(phot,izone)

			phot%edgeNr(izone)=Trac(izone)%edgenext
			phot%inzone(izone)=.true.
			if(phot%i1(izone).gt.Zone(izone)%n1.or.
     &			phot%i2(izone).gt.Zone(izone)%n2.or.
     &			phot%i3(izone).gt.Zone(izone)%n3.or.
     &			phot%i1(izone).lt.1.or.
     &			phot%i2(izone).lt.1.or.
     &			phot%i3(izone).lt.1) then
				phot%inzone(izone)=.false.
			endif
		else
			phot%edgeNr(izone)=0
		endif
	enddo

	goto 1
3	continue
	
	return
	end



	subroutine AddEtraceMono(phot,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer izone
	real*8 v,GetF11
	type(Cell),pointer :: C
	
	do izone=1,nzones
		if(phot%inzone(izone)) then
			C%Escatt=C%Escatt+v*GetF11(phot%ilam1,phot%iscat,C)*phot%sI
		endif
	enddo
  
	return
	end
	


	subroutine InteractMono(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 GetKabs,random,K,GetKp,spec(nlam),KscaP,GetF11,KscaR
	real*8 T0,T1,epsT0,epsT1,kp,tot,GetKsca,increaseT,Ksca(nlam)
	integer izone,iT0,iT1,l,ipart,isize,iT,iscat
	type(Mueller) M
	type(Cell),pointer :: C
	
	KscaR=(phot%Kext-phot%Kabs)*random(idum)
	do izone=1,nzones
		if(phot%inzone(izone)) then
		C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
		do ipart=1,npart
			do isize=1,Part(ipart)%nsize
				do iT=1,Part(ipart)%nT
					KscaP=C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,phot%ilam1)
					KscaR=KscaR-KscaP
					if(KscaR.lt.0d0) then
						M=Part(ipart)%F(isize,iT,phot%ilam1)
						goto 1
					endif
				enddo
			enddo
		enddo
		endif
	enddo

1	continue
	call scatangle(phot,M,iscat)

	call TranslatePhotonV(phot)

	return
	end
	



	subroutine FormalSolution(iobs,ilam)
	use GlobalSetup
	IMPLICIT NONE
	integer iobs,ilam
	
	return
	end
		
	