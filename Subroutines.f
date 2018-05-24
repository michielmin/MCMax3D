	real*8 function random(idum)
	IMPLICIT NONE
	real*8 ran1
	integer idum

	random=ran1(idum)

	return
	end
	
      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
!$OMP THREADPRIVATE(iv,iy)
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END



      FUNCTION gasdev(idum)
      INTEGER idum
      REAL*8 gasdev
CU    USES ran1
      INTEGER iset
      REAL*8 fac,gset,rsq,v1,v2,ran2
c GFORTRAN random needs to be defined to get the correct return type
      REAL*8 random
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*random(idum)-1.
        v2=2.*random(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END



c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotate(x,y,z,u,v,w,theta)
	IMPLICIT NONE
	real*8 x,y,z,u,v,w,yy(3),theta,inp
	real*8 cost,sint,u2,v2,w2
	cost=cos(theta)
	sint=sin(theta)
	u2=u*u
	v2=v*v
	w2=w*w

	inp=x*u+y*v+z*w
	yy(1)=u*inp
     & +(x*(v2+w2)-u*(v*y+w*z))*cost
     & +(v*z-w*y)*sint
	yy(2)=v*inp
     & +(y*(u2+w2)-v*(u*x+w*z))*cost
     & +(w*x-u*z)*sint
	yy(3)=w*inp
     & +(z*(u2+v2)-w*(u*x+v*y))*cost
     & +(u*y-v*x)*sint
	x=yy(1)
	y=yy(2)
	z=yy(3)
	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateZ(x,y,z,cost,sint)
	IMPLICIT NONE
	real*8 x,y,z,xx,yy,r
	real*8 cost,sint

	xx=x*cost-y*sint
	yy=y*cost+x*sint
	x=xx
	y=yy
	

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine rotateY(x,y,z,cost,sint)
	IMPLICIT NONE
	real*8 x,y,z,xx,zz,r
	real*8 cost,sint

	xx=x*cost-z*sint
	zz=z*cost+x*sint
	x=xx
	z=zz

	return
	end
	
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------



      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END


c-----------------------------------------------------------------------
c The new readstar subroutine uses a boxcar filtering to read in 
c high resolution spectra.
c-----------------------------------------------------------------------

	subroutine regridstar(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n,nls
	real*8 grid(n),y(n),x0,y0,xedge(n+1)
	real*8 grid2(n),y2(n),tot(n)
	real*8,allocatable :: ls(:),Fs(:),dls(:)
	character*500 input
	logical truefalse,done(n)

	do i=1,n-1
		xedge(i+1)=sqrt(grid(i)*grid(i+1))
	enddo
	xedge(1)=grid(1)**2/xedge(2)
	xedge(n)=grid(n)**2/xedge(n-1)

	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	j=1
1	read(20,*,end=2,err=1) x0,y0
	j=j+1
	goto 1
2	continue
	nls=j-1
	close(unit=20)
	allocate(ls(nls))
	allocate(dls(nls))
	allocate(Fs(nls))

	open(unit=20,file=input,RECL=100)
	do i=1,nls
3		read(20,*,end=4,err=3) ls(i),Fs(i)
	enddo
4	close(unit=20)
	do i=2,nls-1
		dls(i)=(1d0/ls(i-1)-1d0/ls(i+1))/2d0
	enddo
	dls(1)=dls(2)
	dls(nls)=dls(nls-1)
	
	j=1
	tot=0d0
	y=0d0
	done=.false.
	do i=1,nls
5		continue
		if(ls(i).lt.xedge(j)) goto 6
		if(ls(i).gt.xedge(j+1)) then
			j=j+1
			if(j.gt.n) goto 7
			goto 5
		endif
		y(j)=y(j)+Fs(i)*dls(i)
		tot(j)=tot(j)+dls(i)
		done(j)=.true.
6		continue
	enddo
7	continue
	j=0
	do i=1,n
		if(done(i)) then
			y(i)=y(i)/tot(i)
		else
			j=j+1
			grid2(j)=grid(i)
		endif
	enddo
	if(j.ne.0) then
		call readstar_interpol(input,grid2,y2,j)
		j=0
		do i=1,n
			if(.not.done(i)) then
				j=j+1
				y(i)=y2(j)
			endif
		enddo
	endif
	
	return
	end
	

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

	subroutine readstar_interpol(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		if(y1.gt.1d-60.and.y0.gt.1d-60) then
			y(i)=10d0**(log10(y1)+(log10(grid(i))-log10(x1))*(log10(y0)-log10(y1))/(log10(x0)-log10(x1)))
		else
			y(i)=y1+(grid(i)-x1)*(y0-y1)/(x0-x1)
		endif
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y1*x1**2/grid(j)**2
	enddo
	close(unit=20)
	return
	end


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
	
	subroutine regridlog(input,grid,y,n)
	IMPLICIT NONE
	integer i,j,n
	real*8 grid(n),y(n),x0,y0,x1,y1
	real*8 lx0,ly0,lx1,ly1,lx,ly
	character*500 input
	logical truefalse
	inquire(file=input,exist=truefalse)
	if(.not.truefalse) then
		write(*,200) input(1:len_trim(input))
200		format('File "',a,'" does not exist')
		stop
	endif
	open(unit=20,file=input,RECL=100)
	i=1
1	read(20,*,end=102,err=1) x0,y0
103	if(x0.ge.grid(i)) then
		y(i)=y0
		i=i+1
		goto 103
	endif
100	read(20,*,end=102) x1,y1
	if(y1.le.0d0) y1=y0*1d-50
101	if(grid(i).le.x1.and.grid(i).gt.x0) then
		y(i)=10d0**(log10(y1)+(log10(grid(i)/x1))*(log10(y0/y1))/(log10(x0/x1)))
		i=i+1
		if(i.gt.n) goto 102
		goto 101
	endif
	x0=x1
	y0=y1
	goto 100
102	continue
	do j=i,n
		y(j)=y(i-1)*grid(i-1)/grid(j)
	enddo
	close(unit=20)
	do i=1,n
		if(y(i).le.0d0) y(i)=y0*1d-60
	enddo
	return
	end

