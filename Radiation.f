	subroutine integrate(spec,L)
	use GlobalSetup
	IMPLICIT NONE
	real*8 spec(nlam),L
	integer i
	real*8 nu1,nu2,Iv1,Iv2
	
	L=0d0
	do i=1,nlam-1
		nu1 = nu(i)
		nu2 = nu(i+1)
		Iv1 = spec(i)
		Iv2 = spec(i+1)
		L = L + ABS(nu1-nu2)*0.5*(Iv1+Iv2)
	enddo

	return
	end



	subroutine emit(phot,spec,Lttot)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 spec(nlam),Lr,Lt,random,Ltold,Lttot,x,y,z,r
	integer i,ii,iopac
	type(Photon) phot
	
	specemit=spec/Lttot
	column=0d0

	call randomdirection(phot%vx,phot%vy,phot%vz)
	phot%scatt=.false.
	phot%pol=.false.

	x=-phot%vy
	y=phot%vx
	z=0d0
	r=sqrt(x**2+y**2+z**2)
	phot%Sx=x/r
	phot%Sy=y/r
	phot%Sz=z/r

	Lr=random(idum)*Lttot

	Ltold=0d0
	Lt=0d0
	do i=1,nlam-1
		Lt=Lt+ABS(nu(i)-nu(i+1))*0.5*(spec(i)+spec(i+1))
		if(Lt.ge.Lr) then
			phot%wl1=(Lr-Ltold)/(Lt-Ltold)
			phot%wl2=1d0-phot%wl1
			phot%nu=nu(i)*phot%wl1+nu(i+1)*phot%wl2
			phot%lam=1d4*clight/phot%nu
			phot%ilam1=i
			phot%ilam2=i+1
			goto 1
		endif
		Ltold=Lt
	enddo
	phot%wl1=0d0
	phot%wl2=1d0
	phot%nu=nu(nlam)
	phot%lam=1d4*clight/phot%nu
	phot%ilam1=nlam-1
	phot%ilam2=nlam

1	continue

	phot%UV=.false.
	if(phot%lam.gt.0.0953.and.phot%lam.le.0.206) phot%UV=.true.

	return
	end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


	subroutine randomdirection(x,y,z)
	use GlobalSetup
	IMPLICIT NONE
	real*8 x,y,z,r,random
	
1	continue
	x=2d0*random(idum)-1d0
	y=2d0*random(idum)-1d0
	z=2d0*random(idum)-1d0
	r=x**2+y**2+z**2
	if(r.gt.1d0) goto 1
	r=sqrt(r)
	x=x/r
	y=y/r
	z=z/r
	
	return
	end




	real*8 function increaseT(C)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Cell) C
	real*8 E1,kp0,kp1,GetKp
	integer i,j

	E1=C%E/C%V
	j=int(C%T/dTBB)
	if(j.eq.nBB) then
		increaseT=C%T
		return
	endif
	
	kp0=GetKp(j,C)

	increaseT=0d0
	do i=j,nBB-1
		kp1=GetKp(i+1,C)
		if(kp1.ge.E1) then
			increaseT=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dTBB
			return
		endif
		kp0=kp1
	enddo
c not found, starting from 1 K
	do i=1,j
		kp1=GetKp(i+1,C)
		if(kp0.le.E1.and.kp1.ge.E1) then
			increaseT=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dTBB
			return
		endif
		kp0=kp1
	enddo
	increaseT=real(nBB)*dTBB

	return
	end



	real*8 function determineTslow(C)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Cell) C
	real*8 E1,kp0,kp1,GetKp
	integer i,j

	E1=C%Etrace/C%V
	
	kp0=0d0

	determineTslow=0d0
	do i=1,nBB-1
		kp1=GetKp(i+1,C)
		if(kp0.le.E1.and.kp1.ge.E1) then
			determineTslow=(real(i)**4+(real(i+1)**4-real(i)**4)*(E1-kp0)/(kp1-kp0))**(0.25d0)*dTBB
			return
		endif
		kp0=kp1
	enddo
	determineTslow=real(nBB)*dTBB

	return
	end



	real*8 function determineT(C)
	use GlobalSetup
	IMPLICIT NONE
	type(Cell) C
	real*8 E1,E,T,Emin,Emax,GetKp
	integer i,ii,iopac,iTmin,iTmax,iT0,iT,niter

	E1=C%Etrace/C%V
	iTmin=1
	iTmax=nBB

	iT=C%T/dTBB
	if(iT.lt.1) iT=1
	if(iT.gt.nBB-1) iT=nBB-1

	E=GetKp(iT,C)
	Emax=GetKp(iTmax,C)
	Emin=GetKp(1,C)

	if(IsNaN(E1)) then
		determineT=3d0
		return
	endif

	if(E1.gt.Emax) then
		iT=nBB-1
		determineT=real(iT)*dTBB
		return
	endif

	if(E1.lt.Emin) then
		iT=1
		determineT=real(iT)*dTBB
		return
	endif

	iT0=iT
	niter=0
	do while(abs(iTmax-iTmin).gt.1.and.niter.le.nBB)
		niter=niter+1
		iT=(E1/E)**(0.25)*iT
		if(iT.eq.iT0) then
			if(E1.lt.E) iT=iT0-1
			if(E1.gt.E) iT=iT0+1
		endif
1		continue
		if(iT.le.iTmin) then
			iT=iTmin+1
			goto 1
		endif
		if(iT.ge.iTmax) then
			iT=iTmax-1
			goto 1
		endif
		E=GetKp(iT,C)
		if(E.ge.E1) then
			iTmax=iT
			Emax=E
		endif
		if(E.le.E1) then
			iTmin=iT
			Emin=E
		endif
		iT0=iT
	enddo

	if(iTmin.eq.iTmax) then
		determineT=real(iTmin)*dTBB
	else
		determineT=(real(iTmin)**4+(real(iTmax)**4-real(iTmin)**4)*(E1-Emin)/(Emax-Emin))**(0.25d0)*dTBB
	endif
	
	return
	end




	subroutine Interact(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	real*8 GetKabs,random,K,GetKp,spec(nlam),KscaP,GetF11,KscaR
	real*8 T0,T1,epsT0,epsT1,kp,tot,GetKsca,increaseT,Ksca(nlam)
	integer izone,iT0,iT1,l,ipart,isize,iT,iscat
	type(Mueller) M
	type(Cell),pointer :: C
	
	K=random(idum)*phot%Kext
	if(K.lt.phot%Kabs) then
c absorption and reemission
		spec=0d0
		kp=0d0
		do izone=1,nzones
		if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			C%E=C%E+phot%sI*phot%KabsZ(izone)/phot%Kabs
			T0=C%T
			T1=increaseT(C)

			iT0=int(T0/dTBB)
			iT1=int(T1/dTBB)

			epsT0=T0-real(iT0)*dTBB
			epsT1=T1-real(iT1)*dTBB
			if(iT0.ge.nBB-1) iT0=nBB-2
			if(iT1.ge.nBB-1) iT1=nBB-1
			if(iT0.lt.1) iT0=1
			if(iT1.lt.2) iT1=2

			if(i1totalAbs(izone).ne.phot%i1(izone).or.
     &		   i2totalAbs(izone).ne.phot%i2(izone).or.
     &		   i3totalAbs(izone).ne.phot%i3(izone)) then
				do l=1,nlam
					KabsTotal(izone,l)=GetKabs(l,C)
				enddo
				i1totalAbs(izone)=phot%i1(izone)
				i2totalAbs(izone)=phot%i2(izone)
				i3totalAbs(izone)=phot%i3(izone)
			endif
			if(iT0.eq.iT1) then
				do l=1,nlam
					spec(l)=spec(l)+(BB(l,iT0+1)-BB(l,iT0))*KabsTotal(izone,l)
				enddo
				kp=kp+GetKp(iT0+1,C)-GetKp(iT0,C)
			else
				do l=1,nlam
					spec(l)=spec(l)+(epsT1*BB(l,iT1+1)+(1d0-epsT1)*BB(l,iT1)-epsT0*BB(l,iT0+1)-(1d0-epsT0)*BB(l,iT0))*
     &							KabsTotal(izone,l)
				enddo
				kp=kp+epsT1*GetKp(iT1+1,C)+(1d0-epsT1)*GetKp(iT1,C)-epsT0*GetKp(iT0+1,C)-(1d0-epsT0)*GetKp(iT0,C)
			endif
			C%T=T1
		endif
		enddo
		
		call emit(phot,spec,kp)
	else
c scattering
2		KscaR=(phot%Kext-phot%Kabs)*random(idum)
		do izone=1,nzones
			if(phot%inzone(izone)) then
			C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
			do ipart=1,npart
				do isize=1,Part(ipart)%nsize
					do iT=1,Part(ipart)%nT
						KscaP=phot%wl1*C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,phot%ilam1)
						KscaR=KscaR-KscaP
						if(KscaR.lt.0d0) then
							M=Part(ipart)%F(isize,iT,phot%ilam1)
							goto 1
						endif
						KscaP=phot%wl2*C%densP(ipart,isize,iT)*Part(ipart)%Ksca(isize,iT,phot%ilam2)
						KscaR=KscaR-KscaP
						if(KscaR.lt.0d0) then
							M=Part(ipart)%F(isize,iT,phot%ilam2)
							goto 1
						endif
					enddo
				enddo
			enddo
			endif
		enddo
		goto 2
1		continue
		call scatangle(phot,M,iscat)

c		Ksca=0d0
c		do izone=1,nzones
c			if(phot%inzone(izone)) then
c				if(i1totalSca(izone).ne.phot%i1(izone).or.
c     &			   i2totalSca(izone).ne.phot%i2(izone).or.
c     &			   i3totalSca(izone).ne.phot%i3(izone)) then
c					C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
c					do l=1,nlam
c						KscaTotal(izone,l)=GetF11(l,iscat,C)
c					enddo
c					i1totalSca(izone)=phot%i1(izone)
c					i2totalSca(izone)=phot%i2(izone)
c					i3totalSca(izone)=phot%i3(izone)
c				endif
c				Ksca(1:nlam)=Ksca(1:nlam)+KscaTotal(izone,1:nlam)
c			endif
c		enddo
		do l=1,nlam
			Ksca(l)=Part(ipart)%Ksca(isize,iT,l)*Part(ipart)%F(isize,iT,l)%F11(iscat)
		enddo
		specemit(1:nlam)=specemit(1:nlam)*Ksca(1:nlam)
		call integrate(specemit,tot)
		specemit=specemit/tot
	endif

	phot%x0=phot%x
	phot%y0=phot%y
	phot%z0=phot%z

	call TranslatePhotonV(phot)

	return
	end
	


c------------------------------------------------------------------------
c This subroutine computes a scattering event for photon phot.
c
c The photon should be fully initialized. Used properties are:
c phot%i 			: radial cell number
c phot%j 			: theta cell number
c phot%ilam1		: wavelength number.
c phot%ilam2		: wavelength number.
c phot%wl1			: weight of lam(phot%ilam1)
c phot%wl2			: weight of lam(phot%ilam2)
c			The wavelength of the photon is:
c				phot%lam=phot%wl1*lam(phot%ilam1)+phot%wl2*lam(phot%ilam2)
c phot%pol			: logical determining if the photon is polarized
c phot%vx,vx,vz		: direction of propagation (changed on output)
c phot%Sx,Sy,Sz		: direction of the reference Q vector (changed on output)
c phot%E,Q,U,V		: Stokes vector of the photon (changed on output)
c
c The grain properties used are:
c M%F(i)	: contains all info on the Mueller matrix integrated over grains at wavelength point i
c	F(i)%F11,F12,F22,F33,F34,F44	: elements of the Mueller matrix (normalized)
c	F(i)%IF11						: integral of the F11 and F12 elements weighted with sin(theta)
c 
c------------------------------------------------------------------------
	subroutine scatangle(phot,M,iscat)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	real*8 theta,Fr,random,Fi,FiOld,w1,w2,s1,s2,x,y,z
	real*8 E,Q,U,V,F,IF11,IF12,I0,I1,Itot,phi,wp1,P,P2,r,thet
	type(photon) phot
	integer i,ii,iopac,iscat
	type(Mueller) M

	IF11=M%IF11*180d0
	IF12=M%IF12*180d0
	
	if(phot%pol) then
c Integrate the total scattered intensity over all phi angles
c sin2phi and cos2phi are initialized already to be sin(2*phi) and cos(2*phi)
	Itot=0d0
	do i=1,180
		Itot=Itot+IF11*phot%sI+IF12*(phot%sQ*cos2phi(i)+phot%sU*sin2phi(i))
	enddo
	Itot=random(idum)*Itot
c Determine the phi angle for scattering
	I1=0d0
	I0=0d0
	do i=1,180
		I0=I0+IF11*phot%sI+IF12*(phot%sQ*cos2phi(i)+phot%sU*sin2phi(i))
		if(I0.gt.Itot) then
			wp1=(I0-Itot)/(I0-I1)
			phi=pi*(real(i)-wp1)/180d0
			goto 2
		endif
		I1=I0
	enddo
	i=180
	wp1=0d0
	phi=pi
2	continue
	else
c If incoming photon is unpolarized, phi angle is random.
	phi=random(idum)*pi
	endif

c Scattering plane is defined. Rotate the Stokes vector to the new plane.
	x=phot%Sx
	y=phot%Sy
	z=phot%Sz
c First rotate the axis
	call rotate(x,y,z,phot%vx,phot%vy,phot%vz,phi)
c Then rotate the Stokes vector to the new reference axis
	if(phot%pol) call RotateStokes(phot,x,y,z)
	phot%Sx=x
	phot%Sy=y
	phot%Sz=z

	Fr=180d0*(phot%sI*M%IF11+phot%sQ*M%IF12)

c Now determine the scattering angle theta
	Fr=random(idum)*Fr
	Fi=0d0
	FiOld=0d0
	do i=1,180
		thet=pi*(real(i)-0.5d0)/180d0
		Fi=Fi+pi*sin(thet)*(M%F11(i)*phot%sI+M%F12(i)*phot%sQ)
		if(Fi.gt.Fr) then
			theta=real(i)-(Fi-Fr)/(Fi-FiOld)
			theta=theta*pi/180d0
			goto 1
		endif
		FiOld=Fi
	enddo
	theta=pi
	i=180
1	continue
c Due to symmetry, could also have been -theta
	if(random(idum).lt.0.5d0) theta=-theta
	
c Now determine the scattered stokes vector.
	E=0d0
	Q=0d0
	U=0d0
	V=0d0
	E=E+M%F11(i)*phot%sI
	E=E+M%F12(i)*phot%sQ
	Q=Q+M%F12(i)*phot%sI
	Q=Q+M%F22(i)*phot%sQ
	U=U+M%F33(i)*phot%sU
c if rotation is the counterclockwise F=-F
	F=sign(M%F34(i),theta)
	U=U+F*phot%sV
	V=V-F*phot%sU
	V=V+M%F44(i)*phot%sV

c Now renormalize the photon package to conserve energy.
	if(E.ne.0d0) then
		phot%sQ=phot%sI*Q/E
		phot%sU=phot%sI*U/E
		phot%sV=phot%sI*V/E
	endif

c Rotate the propagation vector.
	call rotate(phot%vx,phot%vy,phot%vz,phot%Sx,phot%Sy,phot%Sz,theta)
	r=sqrt(phot%vx**2+phot%vy**2+phot%vz**2)
	phot%vx=phot%vx/r
	phot%vy=phot%vy/r
	phot%vz=phot%vz/r
	r=sqrt(phot%Sx**2+phot%Sy**2+phot%Sz**2)
	phot%Sx=phot%Sx/r
	phot%Sy=phot%Sy/r
	phot%Sz=phot%Sz/r

c Photon is now polarized !
	if(i.ne.1.and.i.ne.180) phot%pol=.true.
	iscat=i

	return
	end

c------------------------------------------------------------------------
c This subroutine rotates the Stokes vector from 
c one reference plane to another.
c
c The initial reference plane is stored in the phot%Sx,Sy,Sz
c The new reference plane is x,y,z
c The Stokes vector phot%Q,U is changed, the phot%Sx,Sy,Sz are unaltered.
c Note that the vector x,y,z should have unit length!
c------------------------------------------------------------------------
	subroutine RotateStokes(phot,x,y,z)
	use GlobalSetup
	IMPLICIT NONE
	type(photon) phot
	real*8 x,y,z,cost,sint,Q,U,P,P2,cos2t,sin2t,sPP2
	
	P=(phot%sU/phot%sI)**2+(phot%sQ/phot%sI)**2

	if(P.eq.0d0) return

	cost=phot%Sx*x+phot%Sy*y+phot%Sz*z
	sint=(x-
     &	(phot%Sx*(phot%vy**2+phot%vz**2)-phot%vx*(phot%vy*phot%Sy+phot%vz*phot%Sz))*cost)/
     &	(phot%vy*phot%Sz-phot%vz*phot%Sy)
	
	cos2t=cost**2-sint**2
	sin2t=2d0*sint*cost
	
	Q=phot%sQ*cos2t+phot%sU*sin2t
	U=-phot%sQ*sin2t+phot%sU*cos2t

	P2=(U/phot%sI)**2+(Q/phot%sI)**2

	sPP2=sqrt(P/P2)
	phot%sU=sPP2*U
	phot%sQ=sPP2*Q
	
	return
	end



