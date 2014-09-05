	subroutine TravelSph(phot,izone,Trac)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	type(Travel) Trac
	integer izone

	real*8 b,r,R1,R2,T1,T2,vR1,vR2,vT1,vT2,vP1,vP2
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame
	real*8 xt,yt,zt

	xt=phot%x-Zone(izone)%x0
	yt=phot%y-Zone(izone)%y0
	zt=phot%z-Zone(izone)%z0

	r=xt**2+yt**2+zt**2
	R1=Zone(izone)%R2(phot%i1(izone))
	R2=Zone(izone)%R2(phot%i1(izone)+1)
	T1=Zone(izone)%cost2(phot%i2(izone))
	T2=Zone(izone)%cost2(phot%i2(izone)+1)

	b=2d0*(xt*phot%vx+yt*phot%vy+zt*phot%vz)

	if(phot%edgeNr(izone).eq.0) then
		hitR1=hitR(R1,r,b,vR1)
		hitR2=hitR(R2,r,b,vR2)
		hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
		hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
	else
	if(phot%edgeNr(izone).eq.1) then
		hitR1=.false.
		vR1=1d200
		hitR2=hitR(R2,r,b,vR2)
		hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
		hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
	else if(phot%edgeNr(izone).eq.2) then
		hitR1=hitR(R1,r,b,vR1)
		hitR2=.true.
		vR2=-b
		hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
		hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
	else if(phot%edgeNr(izone).eq.3) then
		hitR1=hitR(R1,r,b,vR1)
		hitR2=hitR(R2,r,b,vR2)
		if(Zone(izone)%theta(phot%i2(izone)).lt.(pi/2d0)) then
			hitT1=.false.
			vT1=1d200
			hitT2=hitT(zt,phot%vz,phot,T2,r,b,vT2)
		else
			hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
			hitT2=hitTsame(zt,phot%vz,T2,r,b,vT2)
		endif
	else if(phot%edgeNr(izone).eq.4) then
		hitR1=hitR(R1,r,b,vR1)
		hitR2=hitR(R2,r,b,vR2)
		if(Zone(izone)%theta(phot%i2(izone)).gt.(pi/2d0)) then
			hitT1=hitT(zt,phot%vz,phot,T1,r,b,vT1)
			hitT2=.false.
			vT2=1d200
		else
			hitT1=hitTsame(zt,phot%vz,T1,r,b,vT1)
			hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
		endif
	endif
	endif

	if(.not.hitR2) then
		print*,'Cannot hit outer boundary...'
		print*,phot%x/AU,phot%y/AU,phot%z/AU
		print*,phot%inzone
		print*,phot%i1
		print*,phot%i2
		print*,phot%i3
		stop
	endif

	Trac%v=1d200
	if(hitR1.and.vR1.lt.Trac%v) then
		Trac%v=vR1
		Trac%i1next=phot%i1(izone)-1
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=2
	endif
	if(hitR2.and.vR2.lt.Trac%v) then
		Trac%v=vR2
		Trac%i1next=phot%i1(izone)+1
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=1
	endif
	if(hitT1.and.vT1.lt.Trac%v) then
		Trac%v=vT1
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)-1
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=4
	endif
	if(hitT2.and.vT2.lt.Trac%v) then
		Trac%v=vT2
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)+1
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=3
	endif

	
	return
	end
	

	subroutine HitSph(phot,izone,Trac)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	type(Travel) Trac
	integer izone
	real*8 R1,R2,vR1,vR2,r,b
	logical hitR1,hitR2,hitR
	real*8 xt,yt,zt

	xt=phot%x-Zone(izone)%x0
	yt=phot%y-Zone(izone)%y0
	zt=phot%z-Zone(izone)%z0

	r=xt**2+yt**2+zt**2
	R1=Zone(izone)%R2(1)
	R2=Zone(izone)%R2(Zone(izone)%nr+1)

	b=2d0*(xt*phot%vx+yt*phot%vy+zt*phot%vz)

	hitR1=hitR(R1,r,b,vR1)
	hitR2=hitR(R2,r,b,vR2)
	Trac%v=1d200
	if(hitR1.and.vR1.lt.Trac%v) then
		Trac%v=vR1
		Trac%i1next=1
		Trac%i2next=-1
		Trac%i3next=-1
		Trac%edgenext=2
	endif
	if(hitR2.and.vR2.lt.Trac%v) then
		Trac%v=vR2
		Trac%i1next=Zone(izone)%nr
		Trac%i2next=-1
		Trac%i3next=-1
		Trac%edgenext=1
	endif
	
	
	return
	end
	
	
	
	logical function hitR(Rad,r,b,v)
	use GlobalSetup
	IMPLICIT NONE
	real*8 Rad,r,b,cc,discr,vr1,vr2,v,q
	
	hitR=.false.
	v=1d200

	cc=r-Rad
	discr=(b**2-4d0*cc)
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(b.gt.0d0) then
			q=-0.5d0*(b+discr)
		else
			q=-0.5d0*(b-discr)
		endif
		vr1=q
		vr2=cc/q
		if(vr1.gt.0d0) then
			v=vr1
			hitR=.true.
		endif
		if(vr2.gt.0d0.and.vr2.lt.v) then
			v=vr2
			hitR=.true.
		endif
	endif
	return
	end

	logical function hitT(z,vz,Thet,r,b,v)
	use GlobalSetup
	IMPLICIT NONE
	real*8 Thet,r,b,at,bt,ct,discr,vt1,vt2,v,q,z,vz

	hitT=.false.
	v=1d200

	at=Thet-vz*vz
	bt=Thet*b-2d0*z*vz
	ct=Thet*r-z*z
	discr=bt*bt-4d0*at*ct
	if(discr.ge.0d0) then
		discr=sqrt(discr)
		if(bt.gt.0d0) then
			q=-0.5d0*(bt+discr)
		else
			q=-0.5d0*(bt-discr)
		endif
		vt1=q/at
		vt2=ct/q
		if(vt1.gt.0d0) then
			v=vt1
			hitT=.true.
		endif
		if(vt2.gt.0d0.and.vt2.lt.v) then
			v=vt2
			hitT=.true.
		endif
	endif
	return
	end

	logical function hitTsame(z,vz,Thet,r,b,v)
	use GlobalSetup
	IMPLICIT NONE
	real*8 Thet,r,b,at,bt,v,z,vz

	hitTsame=.true.
	v=1d200

	bt=Thet*b-2d0*z*vz
	at=Thet-vz**2
	v=-bt/at
	if(v.le.0d0) hitTsame=.false.

	return
	end

!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------

	
	