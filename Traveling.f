	subroutine TravelSph(phot,izone,Trac)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	type(Travel) Trac
	integer izone

	real*8 b,r,R1,R2,T1,T2,vR1,vR2,vT1,vT2,vP1,vP2
	logical hitR1,hitR2,hitR,hitT1,hitT2,hitT,hitTsame
	logical hitP1,hitP2,hitP,i1midplane,i2midplane
	real*8 xt,yt,zt,tanx1,tanx2,tany1,tany2,vxt,vyt,vzt

	xt=phot%xzone(izone)
	yt=phot%yzone(izone)
	zt=phot%zzone(izone)
	vxt=phot%vxzone(izone)
	vyt=phot%vyzone(izone)
	vzt=phot%vzzone(izone)

	r=xt**2+yt**2+zt**2
	R1=Zone(izone)%R2(phot%i1(izone))
	R2=Zone(izone)%R2(phot%i1(izone)+1)
	T1=Zone(izone)%cost2(phot%i2(izone))
	T2=Zone(izone)%cost2(phot%i2(izone)+1)
	tanx1=Zone(izone)%tanx(phot%i3(izone))
	tanx2=Zone(izone)%tanx(phot%i3(izone)+1)
	tany1=Zone(izone)%tany(phot%i3(izone))
	tany2=Zone(izone)%tany(phot%i3(izone)+1)

	b=2d0*(xt*vxt+yt*vyt+zt*vzt)

	i1midplane=(phot%i2(izone).eq.Zone(izone)%imidplane)
	i2midplane=((phot%i2(izone)+1).eq.Zone(izone)%imidplane)

	select case(phot%edgeNr(izone))
		case(1)
			hitR1=.false.
			vR1=1d200
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
		case(2)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=.true.
			vR2=-b
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
		case(3)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			if(Zone(izone)%theta(phot%i2(izone)).le.(pi/2d0).or.i1midplane) then
				hitT1=.false.
				vT1=1d200
			else
				hitT1=hitTsame(zt,vzt,T1,r,b,vT1)
			endif
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
		case(4)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			if(Zone(izone)%theta(phot%i2(izone)+1).ge.(pi/2d0).or.i2midplane) then
				hitT2=.false.
				vT2=1d200
			else
				hitT2=hitTsame(zt,vzt,T2,r,b,vT2)
			endif
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
		case(5)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=.false.
			vP1=1d200
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
		case(6)
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=.false.
			vP2=1d200
		case default
			hitR1=hitR(R1,r,b,vR1)
			hitR2=hitR(R2,r,b,vR2)
			hitT1=hitT(zt,vzt,T1,r,b,vT1,i1midplane)
			hitT2=hitT(zt,vzt,T2,r,b,vT2,i2midplane)
			hitP1=hitP(tanx1,tany1,xt,vxt,yt,vyt,vP1)
			hitP2=hitP(tanx2,tany2,xt,vxt,yt,vyt,vP2)
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
	if(hitT1.and.vT1.lt.Trac%v.and.vT1.gt.0d0) then
		Trac%v=vT1
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)-1
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=4
	endif
	if(hitT2.and.vT2.lt.Trac%v.and.vT2.gt.0d0) then
		Trac%v=vT2
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)+1
		Trac%i3next=phot%i3(izone)
		Trac%edgenext=3
	endif
	if(hitP1.and.vP1.lt.Trac%v.and.vP1.gt.0d0) then
		Trac%v=vP1
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)-1
		if(Trac%i3next.lt.1) Trac%i3next=Zone(izone)%np
		Trac%edgenext=6
	endif
	if(hitP2.and.vP2.lt.Trac%v.and.vP2.gt.0d0) then
		Trac%v=vP2
		Trac%i1next=phot%i1(izone)
		Trac%i2next=phot%i2(izone)
		Trac%i3next=phot%i3(izone)+1
		if(Trac%i3next.gt.Zone(izone)%np) Trac%i3next=1
		Trac%edgenext=5
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
	logical hitR1,hitR2,hitRin,hitRout
	real*8 xt,yt,zt,vxt,vyt,vzt

	xt=phot%xzone(izone)
	yt=phot%yzone(izone)
	zt=phot%zzone(izone)
	vxt=phot%vxzone(izone)
	vyt=phot%vyzone(izone)
	vzt=phot%vzzone(izone)

	r=xt**2+yt**2+zt**2
	R1=Zone(izone)%R2(1)
	R2=Zone(izone)%R2(Zone(izone)%nr+1)

	b=2d0*(xt*vxt+yt*vyt+zt*vzt)

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


	logical function hitT(z,vz,Thet,r,b,v,midplane)
	IMPLICIT NONE
	real*8 Thet,r,b,at,bt,ct,discr,vt1,vt2,v,q,z,vz
	logical midplane

	hitT=.false.
	v=1d200

	if(midplane) then
		v=-z/vz
		if(v.gt.0d0) hitT=.true.
		return
	endif

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

	logical function hitP(tanx,tany,x0,vx,y0,vy,v)
	IMPLICIT NONE
	real*8 tanx,tany,x0,vx,y0,vy,v
	
	hitP=.true.
	v=(tany*y0-tanx*x0)/(tanx*vx-tany*vy)
	
	if(v.lt.0d0) hitP=.false.
	
	return
	end


!c-----------------------------------------------------------------------
!c-----------------------------------------------------------------------


	subroutine TranslatePhotonX(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i

	do i=1,nzones
		phot%xzone(i)=phot%x
		phot%yzone(i)=phot%y
		phot%zzone(i)=phot%z
		phot%xzone(i)=phot%xzone(i)-Zone(i)%x0
		phot%yzone(i)=phot%yzone(i)-Zone(i)%y0
		phot%zzone(i)=phot%zzone(i)-Zone(i)%z0
c		phot%xzone(i)=phot%xzone(i)/Zone(i)%xscale
c		phot%yzone(i)=phot%yzone(i)/Zone(i)%yscale
c		phot%zzone(i)=phot%zzone(i)/Zone(i)%zscale
		call rotateZ(phot%xzone(i),phot%yzone(i),phot%zzone(i),cos(-Zone(i)%phi0),sin(-Zone(i)%phi0))
		call rotateY(phot%xzone(i),phot%yzone(i),phot%zzone(i),cos(Zone(i)%theta0),sin(Zone(i)%theta0))
	enddo

	return
	end
	

	subroutine TranslatePhotonV(phot)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i

	do i=1,nzones
		phot%vxzone(i)=phot%vx
		phot%vyzone(i)=phot%vy
		phot%vzzone(i)=phot%vz
		call rotateZ(phot%vxzone(i),phot%vyzone(i),phot%vzzone(i),cos(-Zone(i)%phi0),sin(-Zone(i)%phi0))
		call rotateY(phot%vxzone(i),phot%vyzone(i),phot%vzzone(i),cos(Zone(i)%theta0),sin(Zone(i)%theta0))
	enddo

	return
	end
	
	

	subroutine TravelPhotonX(phot,v)
	use GlobalSetup
	IMPLICIT NONE
	type(Photon) phot
	integer i
	real*8 v

	phot%x=phot%x+v*phot%vx
	phot%y=phot%y+v*phot%vy
	phot%z=phot%z+v*phot%vz
	do i=1,nzones
		phot%xzone(i)=phot%xzone(i)+phot%vxzone(i)*v
		phot%yzone(i)=phot%yzone(i)+phot%vyzone(i)*v
		phot%zzone(i)=phot%zzone(i)+phot%vzzone(i)*v
	enddo

	return
	end
	
	
	
	