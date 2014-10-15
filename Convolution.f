	subroutine Convolution(im,IMDIM,lam0,Diam,Diam2,SpW,fov0,width0,snoise)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer n,IMDIM,NMAX_CONVOLUTION
	parameter(NMAX_CONVOLUTION=4000)
	complex*16,allocatable :: psf(:,:),image(:,:),seeing(:,:)
	complex*16 ic
	parameter(ic=(0d0,1d0))
	real*8 bessj1,x,y,r,rscale,im(IMDIM,IMDIM,3),tot,tot2,tot3
	integer nn(2),ndim,i,j,ix,iy,nsplit,isplit,isign,ii,jj,nsplitpix,n2
	real*8 lam0,Diam,fov0,fov,width,width0,phi,rw,Diam2,SpW,dx,psf2
	integer*8 plan
	integer FFTW_ESTIMATE
	parameter(FFTW_ESTIMATE=64)
	real*8 gasdev,ran2,snoise
	character*500 filename
	complex*16 ctot1,ctot2
	complex*16,allocatable :: CimI(:,:),CimQ(:,:)

	nsplitpix=0
	fov=fov0

1	nsplitpix=nsplitpix+1
	n=IMDIM*nsplitpix
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	x=real(n/4)*rscale*2d0*pi/real(n)
	if(x.lt.Diam.and.IMDIM*(nsplitpix+1)*4.lt.NMAX_CONVOLUTION) goto 1
	if(rscale.lt.(Diam/pi).and.rscale.lt.((3600d0*180d0/pi**2)*lam0/1d6/width0).and.IMDIM*(nsplitpix+1)*4.lt.NMAX_CONVOLUTION) goto 1


	n=1
2	n=n*2
	fov=fov0*real(n)/real(IMDIM)/real(nsplitpix)
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	if(n.lt.IMDIM*nsplitpix*4.and.n.lt.NMAX_CONVOLUTION) goto 2
	if((Diam/(rscale*2d0*pi/real(n))).lt.20d0.and.n.lt.NMAX_CONVOLUTION.and.Diam.gt.0d0) goto 2
	if((Diam2/(rscale*2d0*pi/real(n))).lt.2d0.and.n.lt.NMAX_CONVOLUTION.and.Diam2.gt.0d0) goto 2

	fov=fov0*real(n)/real(IMDIM)/real(nsplitpix)
	rscale=real(n)*lam0/(pi*1d6)*(3600d0*180d0/pi)/fov
	x=real(n/4)*rscale*2d0*pi/real(n)

	write(*,'("--------------------------------------------------------")')
	write(9,'("--------------------------------------------------------")')
	write(*,'("Adjusting pixel scale factor ",i3)') nsplitpix
	write(9,'("Adjusting pixel scale factor ",i3)') nsplitpix
	if(n.ge.10000) then
		write(*,'("FFT on grid of ",i5,"x",i5)') n,n
		write(9,'("FFT on grid of ",i5,"x",i5)') n,n
	else if(n.ge.1000) then
		write(*,'("FFT on grid of ",i4,"x",i4)') n,n
		write(9,'("FFT on grid of ",i4,"x",i4)') n,n
	else if(n.ge.100) then
		write(*,'("FFT on grid of ",i3,"x",i3)') n,n
		write(9,'("FFT on grid of ",i3,"x",i3)') n,n
	else
		write(*,'("FFT on grid of ",i2,"x",i2)') n,n
		write(9,'("FFT on grid of ",i2,"x",i2)') n,n
	endif
	write(*,'("Telescope aperture sampled with ",i5," pixels")') int(Diam/(rscale*2d0*pi/real(n)))
	write(9,'("Telescope aperture sampled with ",i5," pixels")') int(Diam/(rscale*2d0*pi/real(n)))

	if(n.lt.IMDIM*nsplitpix*4) then
		write(*,'("Warning!! cannot do the convolution properly!!")')
		write(9,'("Warning!! cannot do the convolution properly!!")')
	endif

	allocate(seeing(n,n))
	allocate(psf(n,n))
	allocate(image(n,n))

	width=width0

	if(width.gt.0d0) then
	seeing=0d0
	tot2=0d0
	n2=IMDIM*nsplitpix
	if(n2.gt.(n/2)) n2=n/2
	do i=1,n2
	do j=1,n2
		if(i.eq.1.and.j.eq.1) then
			seeing(i,j)=1d0
		else
		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)/rscale
		nsplit=int((5d0/rscale))+10

		do ii=1,nsplit
		do jj=1,nsplit

		x=real(i-1)+real(ii-1)/real(nsplit)
		y=real(j-1)+real(jj-1)/real(nsplit)
		r=sqrt(x**2+y**2)
		rw=r*fov/real(n)

		seeing(i,j)=seeing(i,j)+((1d0+(rw/width)**2)**(-2.5))/real(nsplit*nsplit)
		
		enddo
		enddo
		endif

		tot2=tot2+seeing(i,j)
		seeing(n+1-i,j)=seeing(i,j)
		tot2=tot2+seeing(n+1-i,j)
		seeing(n+1-i,n+1-j)=seeing(i,j)
		tot2=tot2+seeing(n+1-i,n+1-j)
		seeing(i,n+1-j)=seeing(i,j)
		tot2=tot2+seeing(i,n+1-j)
	enddo
	enddo


	seeing=seeing/tot2

	if(snoise.gt.0d0) then
		do i=1,n
		do j=1,n
			psf(i,j)=gasdev(idum)
		enddo
		enddo

		isign=1
		call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
		call dfftw_execute_dft(plan, psf, psf)
		call dfftw_destroy_plan(plan)

		psf=psf/real(n*n)

		ctot1=0d0
		do i=1,n
		do j=1,n
			ctot1=ctot1+psf(i,j)
		enddo
		enddo

		ctot2=snoise
c		psf=ctot1+(psf-ctot1)*ctot2
		psf=psf**ctot2
		psf=psf/cdabs(psf)

		seeing=seeing*psf
	
		tot=0d0
		do i=1,n
		do j=1,n
			tot=tot+seeing(i,j)
		enddo
		enddo
		seeing=seeing/tot
	endif
	
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, seeing,seeing,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, seeing, seeing)
	call dfftw_destroy_plan(plan)
	
c	write(filename,'("seeing.fits")')
c	call writefitsfile(filename,seeing(1:IMDIM,1:IMDIM),IMDIM)

	else
	seeing=1d0
	endif

	if(Diam.gt.0d0) then
	psf=0d0
	tot=0d0
	do i=1,n/2
	do j=1,n/2
		if(i.eq.1.and.j.eq.1.and.Diam2.le.0d0) then
			psf(i,j)=1d0
		else
		x=real(i-1)
		y=real(j-1)
		r=sqrt(x**2+y**2)/rscale
		nsplit=int((5d0/rscale))+10

		dx=rscale*2d0*pi/real(n)/real(nsplit)
		psf2=0d0
		do ii=1,nsplit
		do jj=1,nsplit

		x=real(i-1)+real(ii-1)/real(nsplit)
		y=real(j-1)+real(jj-1)/real(nsplit)
		r=sqrt(x**2+y**2)
		rw=r*fov/real(n)
		r=r*rscale*2d0*pi/real(n)
		x=x*rscale*2d0*pi/real(n)
		y=y*rscale*2d0*pi/real(n)

		if(r.le.Diam.and.r.ge.Diam2) then
			if(abs(x).ge.(SpW).and.abs(y).ge.(SpW)) then
				psf(i,j)=psf(i,j)+1d0/real(nsplit**2)
			endif
			if(abs(x).lt.(SpW).and.(abs(x)+dx).ge.(SpW).and.abs(y).ge.(SpW)) then
				psf(i,j)=psf(i,j)+(1d0/real(nsplit**2))*(abs(x)+dx-SpW)/dx
			endif
			if(abs(y).lt.(SpW).and.(abs(y)+dx).ge.(SpW).and.abs(x).ge.(SpW)) then
				psf(i,j)=psf(i,j)+(1d0/real(nsplit**2))*(abs(y)+dx-SpW)/dx
			endif
			if(cdabs(psf(i,j)).gt.1d0) psf(i,j)=1d0
		endif
		
		enddo
		enddo
		endif

		tot=tot+psf(i,j)
		psf(n+1-i,j)=psf(i,j)
		tot=tot+psf(n+1-i,j)
		psf(n+1-i,n+1-j)=psf(i,j)
		tot=tot+psf(n+1-i,n+1-j)
		psf(i,n+1-j)=psf(i,j)
		tot=tot+psf(i,n+1-j)
	enddo
	enddo
	
c	write(filename,'("aperture.fits")')
c	call writefitsfile(filename,psf(1:IMDIM,1:IMDIM),IMDIM)

	psf=psf/tot
	ndim=2
	nn(1:2)=n

	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, psf, psf)
	call dfftw_destroy_plan(plan)

	tot=0d0
	do i=1,n
	do j=1,n
		psf(i,j)=cdabs(psf(i,j))**2
		if(i.gt.IMDIM*nsplitpix.and.i.lt.(n-IMDIM*nsplitpix)) psf(i,j)=0d0
		if(j.gt.IMDIM*nsplitpix.and.j.lt.(n-IMDIM*nsplitpix)) psf(i,j)=0d0
		tot=tot+psf(i,j)
	enddo
	enddo
	psf=psf/tot

	isign=1

	call dfftw_plan_dft_2d(plan, n,n, psf,psf,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, psf, psf)
	call dfftw_destroy_plan(plan)

	psf=psf/real(n*n)

	else
	psf=1d0/real(n*n)
	endif

	image=0d0
	do i=1,IMDIM
	do j=1,IMDIM
		x=real(i)-real(IMDIM+1)/2d0
		y=real(j)-real(IMDIM+1)/2d0
		if(IMDIM.eq.(2*(IMDIM/2))) then
			x=abs(x)-0.5d0
			y=abs(y)-0.5d0
		endif
		r=sqrt(x**2+y**2)*fov/real(n)
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=im(i,j,1)/real(nsplitpix*nsplitpix)
		enddo
		enddo
	enddo
	enddo
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	image=image*psf*seeing
	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)

	allocate(CimI(IMDIM,IMDIM))

	do i=1,IMDIM
	do j=1,IMDIM
		CimI(i,j)=0d0
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			CimI(i,j)=CimI(i,j)+(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
		enddo
		enddo
	enddo
	enddo

	im(:,:,1)=cdabs(CimI)

	allocate(CimQ(IMDIM,IMDIM))
	image=0d0
	do i=1,IMDIM
	do j=1,IMDIM
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=im(i,j,2)/real(nsplitpix*nsplitpix)
		enddo
		enddo
	enddo
	enddo
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	image=image*psf*seeing
	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	do i=1,IMDIM
	do j=1,IMDIM
		CimQ(i,j)=0d0
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			CimQ(i,j)=CimQ(i,j)+(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
		enddo
		enddo
	enddo
	enddo
		
	im(:,:,2)=(cdabs(CimI+CimQ)-cdabs(CimI-CimQ))/2d0

	image=0d0
	do i=1,IMDIM
	do j=1,IMDIM
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj)=im(i,j,3)/real(nsplitpix*nsplitpix)
		enddo
		enddo
	enddo
	enddo
	isign=1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	image=image*psf*seeing
	isign=-1
	call dfftw_plan_dft_2d(plan, n,n, image,image,
     &                       isign, FFTW_ESTIMATE)
	call dfftw_execute_dft(plan, image, image)
	call dfftw_destroy_plan(plan)
	do i=1,IMDIM
	do j=1,IMDIM
		CimQ(i,j)=0d0
		do ii=1,nsplitpix
		do jj=1,nsplitpix
			CimQ(i,j)=CimQ(i,j)+real(image((i-1)*nsplitpix+ii,(j-1)*nsplitpix+jj))
		enddo
		enddo
	enddo
	enddo

	im(:,:,3)=(cdabs(CimI+CimQ)-cdabs(CimI-CimQ))/2d0


	deallocate(psf)
	deallocate(image)
	deallocate(seeing)
	deallocate(CimI)
	deallocate(CimQ)
	
	return
	end

