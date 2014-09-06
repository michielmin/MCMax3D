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

	select case(phot%edgeNr(izone))
		case(1)
			hitR1=.false.
			vR1=1d200
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
			hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
		case(2)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=.true.
			vR2=-b
			hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
			hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
		case(3)
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
		case(4)
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
		case default
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,phot%vz,T1,r,b,vT1)
			hitT2=hitT(zt,phot%vz,T2,r,b,vT2)
	end select

	if(.not.hitR2) then
		print*,'Cannot hit outer boundary...'
		print*,phot%x/AU,phot%y/AU,phot%z/AU
		print*,phot%vx,phot%vy,phot%vz
		print*,phot%inzone
		print*,phot%i1
		print*,sqrt(phot%x**2+phot%y**2+phot%z**2)/AU
		print*,Zone(1)%R(phot%i1(1))/AU,Zone(1)%R(phot%i1(1)+1)/AU
		stop
	endif

	Trac%v=1d200
	if(hitR1.and.vR1.lt.Trac%v.and.vR1.gt.0d0) then
		Trac%v=vR1
		Trac%i1next=phot%i1(izone)-1
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=2
	endif
	if(hitR2.and.vR2.lt.Trac%v.and.vR2.gt.0d0) then
		Trac%v=vR2
		Trac%i1next=phot%i1(izone)+1
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=1
	endif
C	if(hitT1.and.vT1.lt.Trac%v.and.vT1.gt.0d0) then
C		Trac%v=vT1
C		Trac%i1next=phot%i1(izone)
C		Trac%i2next=phot%i2(izone)-1
C		Trac%i3next=phot%i3(izone)
C		Trac%edgenext=4
C	endif
C	if(hitT2.and.vT2.lt.Trac%v.and.vT2.gt.0d0) then
C		Trac%v=vT2
C		Trac%i1next=phot%i1(izone)
C		Trac%i2next=phot%i2(izone)+1
C		Trac%i3next=phot%i3(izone)
C		Trac%edgenext=3
C	endif

	
	return
	end
	

	subroutine HitSph(phot,izone,Trac)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	type(Travel) Trac
	integer izone
	real*8 R1,R2,vR1,vR2,r,b
	logical hitR1,hitR2,hitRin,hitRout
	real*8 xt,yt,zt

	xt=phot%x-Zone(izone)%x0
	yt=phot%y-Zone(izone)%y0
	zt=phot%z-Zone(izone)%z0

	r=xt**2+yt**2+zt**2
	R1=Zone(izone)%R2(1)
	R2=Zone(izone)%R2(Zone(izone)%nr+1)

	b=2d0*(xt*phot%vx+yt*phot%vy+zt*phot%vz)

	hitR1=.false.
	hitR2=.false.
	hitR1=hitRin(R1,r,b,vR1)
	hitR2=hitRout(R2,r,b,vR2)

	Trac%v=1d200
	if(hitR1.and.vR1.lt.Trac%v) then
		Trac%v=vR1
		Trac%i1next=1
		Trac%i2next=-1
		Trac%i3next=-1
		Trac%edgenext=1
	endif
	if(hitR2.and.vR2.lt.Trac%v) then
		Trac%v=vR2
		Trac%i1next=Zone(izone)%nr
		Trac%i2next=-1
		Trac%i3next=-1
		Trac%edgenext=2
	endif
	
	
	return
	end
	
	
	
	logical function hitR(Rad,r,b,v)
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


	
	logical function hitRin(Rad,r,b,v)
	IMPLICIT NONE
	real*8 Rad,r,b,cc,discr,vr1,vr2,v,q
	
	hitRin=.false.
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
		v=vr1
		hitRin=.true.
		if(vr2.gt.v) then
			v=vr2
			hitRin=.true.
		endif
	endif
	if(v.lt.0d0) hitRin=.false.

	return
	end


	
	logical function hitRout(Rad,r,b,v)
	IMPLICIT NONE
	real*8 Rad,r,b,cc,discr,vr1,vr2,v,q
	
	hitRout=.false.
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
		v=vr1
		hitRout=.true.
		if(vr2.lt.v) then
			v=vr2
			hitRout=.true.
		endif
	endif
	if(v.lt.0d0) hitRout=.false.
	return
	end


	logical function hitT(z,vz,Thet,r,b,v)
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

	
	