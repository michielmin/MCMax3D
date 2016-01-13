	subroutine ReadParticleFits(input,p,isize,iT)
	use GlobalSetup
	use Constants
	implicit none
	character*500 input
	character*80 comment,errmessage
	character*30 errtext
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval,rho_gr,tmp
	logical*4 :: anynull
	integer*4, dimension(4) :: naxes

	real*8,allocatable :: array(:,:),matrix(:,:,:)

	type(particle) p,p0,p1
	integer i,j,ia,iT,iread,nl_read,isize
	logical truefalse,Rayleigh
	real*8 l0,l1,tot,tot2,theta,asym,Pmax,HG,asym2,wasym2
	real rho_av

	allocate(p0%Kabs(1,1,nlam))
	allocate(p0%Ksca(1,1,nlam))
	allocate(p0%Kext(1,1,nlam))
	allocate(p0%F(1,1,nlam))
	allocate(p1%Kabs(1,1,nlam))
	allocate(p1%Ksca(1,1,nlam))
	allocate(p1%Kext(1,1,nlam))
	allocate(p1%F(1,1,nlam))

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,input,readwrite,blocksize,status)
	if (status /= 0) then
		call output("Error reading particle file: " // trim(input))
		call output("==================================================================")
		stop
	endif
	group=1
	firstpix=1
	nullval=-999

	call output("Reading particle file: " // trim(input))

	!------------------------------------------------------------------------
	! HDU0 : opacities
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

	nl_read = naxes(1)

	npixels=naxes(1)*naxes(2)

	! Read model info

	call ftgkye(unit,'a1',p%dust_moment1,comment,status)
	call ftgkye(unit,'a2',p%dust_moment2,comment,status)
	p%rv(isize)=sqrt(p%dust_moment2)*1d-4
	call ftgkye(unit,'a3',p%dust_moment3,comment,status)
	call ftgkye(unit,'density',rho_av,comment,status)
	p%rho(iT)=rho_av

c	call ftgkyj(unit,'mcfost2prodimo',mcfost(1)%mcfost2ProDiMo,comment,stat4)
 
	! read_image
	allocate(array(nl_read,4))

	call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)


	!------------------------------------------------------------------------
	! HDU 1: matrix
	!------------------------------------------------------------------------
	!  move to next hdu
	call ftmrhd(unit,1,hdutype,status)	

	! Check dimensions
	call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)

	! read_image
	allocate(matrix(nl_read,6,180))

	call ftgpvd(unit,group,firstpix,npixels,nullval,matrix,anynull,status)
   
				 
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   call output('error in reading fits file' // int2string(status,'(i6)'))

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif


	i=1
	iread=1
	l0=array(i,1)
	p0%Kext(1,1,1)=array(iread,2)
	p0%Kabs(1,1,1)=array(iread,3)
	p0%Ksca(1,1,1)=array(iread,4)
	do j=1,180
		p0%F(1,1,1)%F11(j)=matrix(iread,1,j)
		p0%F(1,1,1)%F12(j)=matrix(iread,2,j)
		p0%F(1,1,1)%F22(j)=matrix(iread,3,j)
		p0%F(1,1,1)%F33(j)=matrix(iread,4,j)
		p0%F(1,1,1)%F34(j)=matrix(iread,5,j)
		p0%F(1,1,1)%F44(j)=matrix(iread,6,j)
	enddo
103	if(l0.ge.lam(i)) then
		p%Kext(isize,iT,i)=p0%Kext(1,1,1)
		p%Ksca(isize,iT,i)=p0%Ksca(1,1,1)
		p%Kabs(isize,iT,i)=p0%Kabs(1,1,1)
		p%F(isize,iT,i)=p0%F(1,1,1)
		call tellertje(i,nlam)
		i=i+1
		goto 103
	endif
100	iread=iread+1
	if(iread.gt.nl_read) goto 102
	l1=array(iread,1)
	p1%Kext(1,1,1)=array(iread,2)
	p1%Kabs(1,1,1)=array(iread,3)
	p1%Ksca(1,1,1)=array(iread,4)
	do j=1,180
		p1%F(1,1,1)%F11(j)=matrix(iread,1,j)
		p1%F(1,1,1)%F12(j)=matrix(iread,2,j)
		p1%F(1,1,1)%F22(j)=matrix(iread,3,j)
		p1%F(1,1,1)%F33(j)=matrix(iread,4,j)
		p1%F(1,1,1)%F34(j)=matrix(iread,5,j)
		p1%F(1,1,1)%F44(j)=matrix(iread,6,j)
	enddo
101	if(lam(i).le.l1.and.lam(i).ge.l0) then
		p%Kext(isize,iT,i)=exp(log(p1%Kext(1,1,1))+log(lam(i)/l1)*log(p0%Kext(1,1,1)/p1%Kext(1,1,1))/log(l0/l1))
		p%Ksca(isize,iT,i)=exp(log(p1%Ksca(1,1,1))+log(lam(i)/l1)*log(p0%Ksca(1,1,1)/p1%Ksca(1,1,1))/log(l0/l1))
		p%Kabs(isize,iT,i)=exp(log(p1%Kabs(1,1,1))+log(lam(i)/l1)*log(p0%Kabs(1,1,1)/p1%Kabs(1,1,1))/log(l0/l1))
		p%F(isize,iT,i)%F11(1:180)=p1%F(1,1,1)%F11(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F11(1:180)-p1%F(1,1,1)%F11(1:180))/(l0-l1)
		p%F(isize,iT,i)%F12(1:180)=p1%F(1,1,1)%F12(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F12(1:180)-p1%F(1,1,1)%F12(1:180))/(l0-l1)
		p%F(isize,iT,i)%F22(1:180)=p1%F(1,1,1)%F22(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F22(1:180)-p1%F(1,1,1)%F22(1:180))/(l0-l1)
		p%F(isize,iT,i)%F33(1:180)=p1%F(1,1,1)%F33(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F33(1:180)-p1%F(1,1,1)%F33(1:180))/(l0-l1)
		p%F(isize,iT,i)%F34(1:180)=p1%F(1,1,1)%F34(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F34(1:180)-p1%F(1,1,1)%F34(1:180))/(l0-l1)
		p%F(isize,iT,i)%F44(1:180)=p1%F(1,1,1)%F44(1:180)+(lam(i)-l1)*(p0%F(1,1,1)%F44(1:180)-p1%F(1,1,1)%F44(1:180))/(l0-l1)
		call tellertje(i,nlam)
		i=i+1
		if(i.gt.nlam) goto 102
		goto 101
	endif
	l0=l1
	p0%Kext(1,1,1)=p1%Kext(1,1,1)
	p0%Ksca(1,1,1)=p1%Ksca(1,1,1)
	p0%Kabs(1,1,1)=p1%Kabs(1,1,1)
	p0%F(1,1,1)=p1%F(1,1,1)
	goto 100
102	continue
	do j=i,nlam
		call tellertje(j,nlam)
		p%Ksca(isize,iT,j)=p%Ksca(isize,iT,i-1)*(lam(i-1)/lam(j))**4
		p%Kabs(isize,iT,j)=p%Kabs(isize,iT,i-1)*(lam(i-1)/lam(j))**2
		p%Kext(isize,iT,j)=p%Kabs(isize,iT,j)+p%Ksca(isize,iT,j)
		p%F(isize,iT,j)=p%F(isize,iT,i-1)
	enddo

	do j=1,nlam
		tot=0d0
		tot2=0d0
		do i=1,180
			tot=tot+p%F(isize,iT,j)%F11(i)*sin(pi*(real(i)-0.5)/180d0)
			tot2=tot2+sin(pi*(real(i)-0.5)/180d0)
		enddo
		do i=1,180
			p%F(isize,iT,j)%F11(i)=tot2*p%F(isize,iT,j)%F11(i)/tot
			p%F(isize,iT,j)%F12(i)=tot2*p%F(isize,iT,j)%F12(i)/tot
			p%F(isize,iT,j)%F22(i)=tot2*p%F(isize,iT,j)%F22(i)/tot
			p%F(isize,iT,j)%F33(i)=tot2*p%F(isize,iT,j)%F33(i)/tot
			p%F(isize,iT,j)%F34(i)=tot2*p%F(isize,iT,j)%F34(i)/tot
			p%F(isize,iT,j)%F44(i)=tot2*p%F(isize,iT,j)%F44(i)/tot
		enddo
	enddo

	deallocate(p0%Kabs)
	deallocate(p0%Ksca)
	deallocate(p0%Kext)
	deallocate(p0%F)
	deallocate(p1%Kabs)
	deallocate(p1%Ksca)
	deallocate(p1%Kext)
	deallocate(p1%F)

	deallocate(array)
	deallocate(matrix)

	return

	end


