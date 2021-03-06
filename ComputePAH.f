	subroutine ComputePAH(p,amin,amax,apow)
	use Parameters
	IMPLICIT NONE
	type(particle) p
	real*8 HC,x,lam1,lam2,rV,amin,amax,apow,theta,Mc,fn_i,fn_n,fn,tot,tot2
	integer i,j,ilam,ia
	logical ionized
	parameter(Mc=12d0*1.66d-24) !mass of a carbon atom in gram
	real*8,allocatable :: Ka_n(:),Ka_i(:),Ks_n(:),Ks_i(:)

	p%qhp=.true.
	p%gascoupled=.true.
	
	allocate(Ka_n(nlam))
	allocate(Ka_i(nlam))
	allocate(Ks_n(nlam))
	allocate(Ks_i(nlam))

c	rV=sqrt( (amax**(3d0-apow)-amin**(3d0-apow))*(1d0-apow)/((amax**(1d0-apow)-amin**(1d0-apow))*(3d0-apow)) )
	rV=sqrt(amax*amin)

	p%Nc=(rV*1d3)**3d0*468d0

	if(p%Nc.lt.25) then
		HC=0.5d0
	else if(p%Nc.lt.100) then
		HC=0.5/sqrt(p%Nc/25d0)
	else
		HC=0.25d0
	endif

	p%Mc=12d0
	p%Td_qhp=450d0

	p%rv=(p%Nc/468d0)**(1d0/3d0)*1d-7
	p%rho(1:p%nopac)=p%Nc*((12d0+HC)*Mc/12d0)/(4d0*pi*p%rv**3/3d0)

	write(*,'("--------------------------------------------------------")')
	write(*,'("Computing PAH opacities")')
	write(*,'("PAH size:              ",f12.6)') p%rV*1d4
	write(*,'("Number of carbon atoms:",f12.1)') p%Nc
	write(9,'("--------------------------------------------------------")')
	write(9,'("Computing PAH opacities")')
	write(9,'("PAH size:              ",f12.6)') p%rV*1d4
	write(9,'("Number of carbon atoms:",f12.1)') p%Nc

	ionized=.false.
	call MakePAH(lam,Ka_n,Ks_n,p%Nc,HC,nlam,ionized)
	fn_n=0.473692
	ionized=.true.
	call MakePAH(lam,Ka_i,Ks_i,p%Nc,HC,nlam,ionized)
	fn_i=0.0983061

c	fn=1d0
c	fn=fn_n
c	p%Kabs(1,1:nlam)=10d0**(((fn-fn_i)*log10(Ka_n)+(fn_n-fn)*log10(Ka_i))/(fn_n-fn_i))
c	p%Ksca(1,1:nlam)=10d0**(((fn-fn_i)*log10(Ks_n)+(fn_n-fn)*log10(Ks_i))/(fn_n-fn_i))

c	fn=0d0
c	fn=fn_i
c	p%Kabs(2,1:nlam)=10d0**(((fn-fn_i)*log10(Ka_n)+(fn_n-fn)*log10(Ka_i))/(fn_n-fn_i))
c	p%Ksca(2,1:nlam)=10d0**(((fn-fn_i)*log10(Ks_n)+(fn_n-fn)*log10(Ks_i))/(fn_n-fn_i))

	p%Kabs(1,1:nlam)=Ka_n(1:nlam)
	p%Kabs(2,1:nlam)=Ka_i(1:nlam)
	p%Ksca(1,1:nlam)=Ks_n(1:nlam)
	p%Ksca(2,1:nlam)=Ks_i(1:nlam)

c	open(unit=32,file='PAH.dat',RECL=6000)
c	do i=1,nlam
c		write(32,*) lam(i),p%Kabs(1,i),p%Kabs(2,i),p%Ksca(1,i),p%Ksca(2,i),Ka_n(i),Ka_i(i),Ks_n(i),Ks_i(i)
c	enddo
c	close(unit=32)

	tot=0d0
	tot2=0d0
	do ia=1,180
		theta=(real(ia)-0.5d0)*pi/180d0
		tot=tot+((1d0+cos(theta)**2)/2d0)*sin(theta)
		tot2=tot2+sin(theta)
	enddo

	do i=1,nlam
		do ia=1,180
			theta=(real(ia)-0.5d0)*pi/180d0
			p%F(1,i)%F11(ia)=(tot2/tot)*(1d0+cos(theta)**2)/2d0
			p%F(1,i)%F12(ia)=-(tot2/tot)*(1d0-cos(theta)**2)/2d0
			p%F(1,i)%F22(ia)=(tot2/tot)*(1d0+cos(theta)**2)/2d0
			p%F(1,i)%F33(ia)=(tot2/tot)*cos(theta)
			p%F(1,i)%F34(ia)=0d0
			p%F(1,i)%F44(ia)=(tot2/tot)*cos(theta)
			p%F(2,i)%F11(ia)=(tot2/tot)*(1d0+cos(theta)**2)/2d0
			p%F(2,i)%F12(ia)=-(tot2/tot)*(1d0-cos(theta)**2)/2d0
			p%F(2,i)%F22(ia)=(tot2/tot)*(1d0+cos(theta)**2)/2d0
			p%F(2,i)%F33(ia)=(tot2/tot)*cos(theta)
			p%F(2,i)%F34(ia)=0d0
			p%F(2,i)%F44(ia)=(tot2/tot)*cos(theta)
		enddo
		p%Kext(1,i)=p%Kabs(1,i)+p%Ksca(1,i)
		p%Kext(2,i)=p%Kabs(2,i)+p%Ksca(2,i)
	enddo



	deallocate(Ka_n)
	deallocate(Ka_i)
	deallocate(Ks_n)
	deallocate(Ks_i)

	return
	end

	
	subroutine MakePAH(lam,Cabs,Csca,Nc,HC,nlam,ionized)
	IMPLICIT NONE
	integer nlam,i,j,ilam,nj
	parameter(nj=30)
	real*8 lam(nlam),Cabs(nlam),Csca(nlam),Nc,HC,x,cutoffPAH,Mc
	real*8 e1x(nlam),e2x(nlam),e1y(nlam),e2y(nlam),CabsGra
	real*8 a,fpah,qgra
	real*8 r,QEX,QSC,QAB,G
	character*100 filename
	logical ionized
	real*8 lj(nj),gj(nj),sjn(nj),sji(nj),S(nj),pi
	parameter(pi=3.1415926536)
	parameter(Mc=2e-23) ! in grams
	data (lj(j),j=1,30) / 0.0722, 0.2175,1.050,1.260,1.905,3.300,5.270,5.700,
     &	6.220,6.690,7.417,7.598,7.850,8.330,8.610,10.68,11.23,11.33,11.99,12.62,
     &	12.69,13.48,14.19,15.90,16.45,17.04,17.375,17.87,18.92,15.0 /
	data (gj(j),j=1,30) / 0.195,0.217,0.055,0.11,0.09,0.012,0.034,0.035,0.030,
     &	0.070,0.126,0.044,0.053,0.052,0.039,0.020,0.012,0.032,0.045,0.042,0.013,
     &	0.040,0.025,0.020,0.014,0.065,0.012,0.016,0.10,0.8 / 
	data (sjn(j),j=1,30) / 7.97d7,1.23d7,0.0,0.0,0.0,394.0,2.5,4.0,29.4,7.35,
     &	20.8,18.1,21.9,6.94,27.8,0.3,18.9,52.0,24.2,35.0,1.3,8.0,0.45,0.04,0.5,
     &	2.22,0.11,0.067,0.10,50.0 /
	data (sji(j),j=1,30) / 7.97d7,1.23d7,2.0d4,0.078,-146.5,89.4,20.0,32.0,
     &	235.0,59.0,181.0,163.0,197.0,48.0,194.0,0.3,17.7,49.0,20.5,31.0,1.3,8.0,
     &	0.45,0.04,0.5,2.22,0.11,0.067,0.17,50.0 /
	external Graphite_x,Graphite_y

	r=1d-3*(Nc/468d0)**(1d0/3d0)

	call RegridDataLNK(Graphite_x,lam,e1x,e2x,nlam,.false.)
	call RegridDataLNK(Graphite_y,lam,e1y,e2y,nlam,.false.)

	do ilam=1,nlam
		call Q_MIE(e1x(ilam),e2x(ilam),lam(ilam),r,QEX,QSC,QAB,G)
		CabsGra=1d-8*(qab*pi*r**2)/3d0
		Csca(ilam)=1d-8*(qsc*pi*r**2)/3d0
		call Q_MIE(e1y(ilam),e2y(ilam),lam(ilam),r,QEX,QSC,QAB,G)
		CabsGra=CabsGra+2d0*1d-8*(qab*pi*r**2)/3d0
		Csca(ilam)=Csca(ilam)+2d0*1d-8*(qsc*pi*r**2)/3d0
		
		CabsGra=CabsGra/Nc
		Csca(ilam)=Csca(ilam)/NC
		
		x=1d0/lam(ilam)
		do j=1,nj
			S(j)=2d0*gj(j)*lj(j)*1d-4/(pi*((lam(ilam)/lj(j)-lj(j)/lam(ilam))**2+gj(j)**2))
			if(ionized) then
				if(j.eq.6) then
c use the 3.3 micron feature from Visser 2007
					S(j)=S(j)*sjn(j)*1d-20/(1d0+41d0/(Nc-14d0))
				else
					S(j)=S(j)*sji(j)*1d-20
				endif
			else
				S(j)=S(j)*sjn(j)*1d-20
			endif
		enddo
		S(6)=S(6)*HC
		do j=14,22
			S(j)=S(j)*HC
		enddo
		if(x.gt.17.25d0) then
			Cabs(ilam)=CabsGra
		else if(x.gt.15d0) then
			Cabs(ilam)=(126.0-6.4943*x)*1e-18
		else if(x.gt.10d0) then
			Cabs(ilam)=S(1)+(-3.0+1.35*x)*1e-18
		else if(x.gt.7.7d0) then
			Cabs(ilam)=(66.302-24.367*x+2.950*x**2-0.1057*x**3)*1e-18
		else if(x.gt.5.9d0) then
			Cabs(ilam)=S(2)+(1.8687+0.1905*x+0.4175*(x-5.9)**2+0.04370*(x-5.9)**3)*1e-18
		else if(x.gt.3.3d0) then
			Cabs(ilam)=S(2)+(1.8687+0.1905*x)*1e-18
		else
			Cabs(ilam)=34.58*10d0**(-18d0-3.431/x)*cutoffPAH(lam(ilam),Nc,ionized)
			do j=3,nj
				Cabs(ilam)=Cabs(ilam)+S(j)
			enddo
		endif
		a=(50d-4/r)**3
		if(a.gt.1d0) a=1d0
		qgra=0.01d0
		fpah=(1d0-qgra)*a
		Cabs(ilam)=(Cabs(ilam)*fpah+CabsGra*(1d0-fpah))/Mc
		if(Cabs(ilam).lt.0d0) Cabs(ilam)=0d0
		Csca(ilam)=Csca(ilam)/Mc
	enddo
	
	
	
			
	return
	end
	



	real*8 function cutoffPAH(lam,Nc,ionized)
	IMPLICIT NONE
	real*8 lam,Nc,y,M
	logical ionized
			
	if(Nc.gt.40d0) then
		M=0.4*Nc
	else
		M=0.3*Nc
	endif
	
	if(ionized) then
		y=1d0/(2.282*M**(-0.5)+0.889)
	else
		y=1d0/(3.804*M**(-0.5)+1.052)
	endif		
	y=y/lam
	
	cutoffPAH=atan(10d3*(y-1d0)**3/y)/3.1415926536+0.5d0
	
	return
	end
	

	
***********************************************************************
*       New Mie subroutine that approximates for big grains           *
*                                                                     *
***********************************************************************
      SUBROUTINE Q_MIE(E1,E2,LAM,RAD,QEX,QSC,QAB,G)
      IMPLICIT real*8 (A-H,O-Z)
      real*8 LAM,RAD,T,QEX,QSC,QAB,E1,E2,G,EV


C
C     MIE THEORY EFFICIENCY FACTORS FOR SPHERICAL PARTICLES OF
C     RADIUS 'RAD' AT WAVELENGTH 'LAM'.
C     E=E1 + I*E2 IS THE SQUARE OF THE COMPLEX REFRACTIVE INDEX.
C     THE REFRACTIVE INDEX IS GIVEN BY SUBROUTINE 'EPS'
C                         
      COMPLEX*16 E,RM,Y,ZN,ZN1,ZN2,C,A,B,AO,RRAT,A1,ANM1,BNM1
	complex*16,allocatable :: AN(:)

	T=0.0

      E=DCMPLX(E1,-E2)
      E=E**2.


      X=6.2831853*RAD/LAM
      IF(X.LT.0.001)THEN
C
C        USE SMALL PARTICLE FORMULAE.
C        CHANGED CRITERIION FROM X < 0.01 TO 0.001 BECAUSE SILICATE
C        SCATTERING WAS NOT CORRECT.
C	15-8-2001: Changed scattering formula from QSC=(X**4/.375)*DBLE(C**2)
C				into the correct formula QSC=(X**4/.375)*DABS(C)**2
C	Michiel Min
C
         C=(E-1.)/(E+2.)
         QSC=(X**4/.375)*cDABS(C)**2
         A=DIMAG(-4.*C)
         B=DIMAG(-C*(E*E+27.*E+38.)/(2.*E+3.)/3.75)
         QAB=X*(A+X*X*B)
         QEX=QAB+QSC
C
C        G THE ASYMMETRY PARAMETER IS ALWAYS NEGLIGIBLE FOR SMALL PARTICLES.
C
         G=0.0
         RETURN
      END IF
C
C     FULL MIE THEORY CALCULATION.
C     RM - COMPLEX REFRACTIVE INDEX
C
      RM=CDSQRT(E)
      EN1=DBLE(RM)
      EN2=DIMAG(RM)
      Y=X*RM
      ZN2=DCMPLX(DCOS(X),-DSIN(X))
      ZN1=DCMPLX(DSIN(X),DCOS(X))
      RIND=EN1**2+EN2**2     ! Rind = |rm|�
      NTIL=1.5*SQRT(RIND)*X+1

c	Number of iterations changed to improve for small |m| (Michiel Min)
	if(real(ntil).lt.(1.5*x))ntil=1.5*x

      NTOT=MAX0(20,NTIL)
c
      if (ntot.le.70000) then    ! go ahead with full Mie theory
	allocate(AN(NTOT))
c
      AN(NTOT)=DCMPLX(0,0)
      SUME=0.
      SUMS=0.
      SUMG1=0.
      SUMG2=0.
      PSG1=0.
      PSG2=0.
      NTOTA=NTOT
  100 P=DFLOAT(NTOTA)
      AN(NTOTA-1)=P/Y-(1./(P/Y+AN(NTOTA)))
      NTOTA=NTOTA-1
      IF(NTOTA.EQ.1) GOTO 101
      GOTO 100
  101 AO1=DSIN(EN1*X)*DCOS(EN1*X)
      EN2P=-EN2
c      IF(EN2P*X.GE.44.)WRITE(6,*)'EN2P,X,LAM,RAD,E1,E2',EN2P,X,LAM,
c     >RAD,E1,E2
      if(EN2P*X.GE.350.) then
         AO=dcmplx(0.0,1.0)
      else
        AO2=DSINH(EN2P*X)*DCOSH(EN2P*X)
        AO3=(DSIN(EN1*X))**2+(DSINH(EN2P*X))**2
        AO=DCMPLX(AO1,AO2)
        AO=AO/AO3
      endif
      A1=-1./Y+(1./(1./Y-AO))
      RRAT=A1/AN(1)
      f=2.0/(x*x)
      DO 4 N=1,NTOT
         AN(N)=AN(N)*RRAT
    4 CONTINUE 
      DO 2 N=1,NTOT
         P=DFLOAT(N)
         ZN=DFLOAT(2*N-1)*ZN1/X-ZN2
         C=AN(N)/RM+P/X
         A=C*DBLE(ZN)-DBLE(ZN1)
         A=A/(C*ZN-ZN1)
         C=RM*AN(N)+P/X
         B=C*DBLE(ZN)-DBLE(ZN1)
         B=B/(C*ZN-ZN1)
C
C        PP, PPG1, PPG2 ARE CONSTANTS CONTAINING THE N TERMS IN THE 
C        SUMMATIONS.
C
         PP=DFLOAT(2*N+1)
C
         PSS=PP*(A*dCONJG(A)+B*dCONJG(B))
         PSE=PP*DBLE(A+B)
         IF(N.GT.1)THEN
C
C           CALCULATE G USING FORMULA ON P.128 OF VAN DE HULST'S BOOK.
C           HAVE REPLACED N BY (N-1) IN THE FORMULA SO THAT WE CAN USE
C           PREVIOUS A(N) AND B(N) INSTEAD OF A(N+1) AND B(N+1)
C
            REN=DFLOAT(N)
            PPG1=(REN-1.)*(REN+1.)/REN
            PPG2=(2.*REN-1.)/((REN-1.)*REN)
            PSG1=PPG1*DBLE(ANM1*dCONJG(A)+BNM1*dCONJG(B))
            PSG2=PPG2*DBLE(ANM1*dCONJG(BNM1))
         END IF
         SUME=SUME+PSE
         SUMS=SUMS+PSS
         SUMG1=SUMG1+PSG1
         SUMG2=SUMG2+PSG2
         D1=ABS(PSE/SUME)
         D2=ABS(PSS/SUMS)
C        IF(D1.LT.1.E-7.AND.D2.LT.1.E-7) GO TO 5
         PT=ABS(PSS/PP)
         IF(PT.LE.1.E-20) GOTO 5
C
C        SAVE PREVIOUS A AND B FOR CALCULATION OF G THE ASYMMETRY PARAMETER
C
         ANM1=A
         BNM1=B
         ZN2=ZN1
         ZN1=ZN
    2 CONTINUE
    5 F=2.0/(X*X)
      QEX=F*SUME
      QSC=F*SUMS
      QAB=F*(SUME-SUMS)
      G=2.0*F*(SUMG1+SUMG2)/QSC
	deallocate(AN)
      RETURN
      else
c               Geometrical optics for big spheres
      call geopt(rm,ans)
      qex =2.0d0
      g=9.23d-01   !approx true for D&L silicate.......
      qsc=ans
      end if
      return
      END          
c******************************************************************************
      subroutine geopt(m,ans)
c      intgrates the reflection coefficient
c      trapezium rule integration from 0 to pi/2
      implicit real*8 (a-h,o-z)
      complex*16 m
      a=0.0d0
      b=1.570796327d0
      nstrip = 5000
      tot=0
      h=(b-a)/dfloat(nstrip)   !strip width
      tot=tot+0.5*ref(m,a)*h   !1st term
      do i=1,nstrip-1
       x=a+h*dfloat(i)
       tot=tot+ref(m,x)*h      !middle terms
      end do
      tot=tot+0.5*ref(m,b)*h   !last term
      ans=1.+2.*tot    !ans is Qsca
      return
      end
      
c******************************************************************************
                                                                               
      function ref(m,thetai)
c         Calculates Reflection coeffs
      implicit real*8 (a-h,o-z)
      complex*16 sinTHETAt,cosTHETAt ,m,rpll,rper          
      sinTHETAt=sin(THETAi)/m
      cosTHETAt=cdsqrt(1-(sinTHETAt*sinTHETAt))
c       r for E parallel to plane
      rpll = (cosTHETAt-m*cos(THETAi)) / (cosTHETAt+m*cos(THETAi))
c       r for E perp. to plane
      rper = (cos(THETAi)-m*cosTHETAt) / (cos(THETAi)+m*cosTHETAt)
C       R = �(|rpll|�+|rper|�)
      R= (abs(rpll)*abs(rpll) + abs(rper)*abs(rper))/2.0
      ref=r*sin(THETAi)*cos(THETAi)
      return                                            
      end                 

C
C
      SUBROUTINE INTERP(X,Y,NPTS,NTERMS,XIN,YOUT)
      REAL*8 DELTAX,PROD,SUM,X(3000),Y(3000),XIN,YOUT
      REAL*8 DELTA(10),A(10)
      REAL*8 DENOM
**************************************************   
*     SEARCH FOR AN APPROPRIATE VALUE OF X(1)    *
**************************************************
   11 DO 19 I=1,NPTS
      IF (XIN-X(I)) 13,17,19
   13 I1=I-NTERMS/2
      IF(I1) 15,15,21
   15 I1=1
      GOTO 21
   17 YOUT=Y(I)
   18 GOTO 61
   19 CONTINUE
      I1=NPTS-NTERMS+1
   21 I2=I1+NTERMS-1
      IF (NPTS-I2) 23,31,31
   23 I2=NPTS
      I1=I2-NTERMS+1
   25 IF (I1) 26,26,31
   26 I1=1
   27 NTERMS=I2-I1+1
C
C  EVALUATE DEVIATIONS DELTA
C
   31 DENOM=X(I1+1)-X(I1)
      DELTAX=(XIN-X(I1))/DENOM
      DO 35 I=1,NTERMS
         IX=I1+I-1
   35 DELTA(I)=(X(IX)-X(I1)) / DENOM
**********************************************
*           ACCUMULATE COEFFICIENTS A        *
**********************************************
   40 A(1)=Y(I1)
   41 DO 50 K=2,NTERMS
         PROD=1.
         SUM=0.
         IMAX=K-1
         IXMAX=I1+IMAX
         DO 49 I=1,IMAX
            J=K-I
            PROD=PROD*(DELTA(K)-DELTA(J))
   49    SUM=SUM-A(J)/PROD
   50 A(K)=SUM+Y(IXMAX)/PROD
***********************************************
*         ACCUMULATE SUM OF EXPANSION         *
***********************************************
   51 SUM=A(1)
      DO 57 J=2,NTERMS
         PROD=1.
         IMAX=J-1
         DO 56 I=1,IMAX
   56    PROD=PROD*(DELTAX-DELTA(I))
   57 SUM=SUM+A(J)* PROD
   60 YOUT=SUM
   61 RETURN
      END                

