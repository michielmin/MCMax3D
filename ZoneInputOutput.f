	subroutine outputstruct_fits(filename,vars,nvars,izone)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,k,ipart,is,iT,izone
	character*7 vars(nvars),hdu
	character*500 filename
	real*8,allocatable :: array(:,:,:,:,:,:)
	integer status,unit,blocksize,bitpix,naxis,naxes(6)
	integer group,fpixel,nelements
	logical simple,extend,truefalse
	type(ZoneType),pointer :: ZZ

	ZZ => Zone(izone)

	open(unit=32,file='midplanedens.dat',RECL=6000)
	do i=1,ZZ%nr
		write(32,*) sqrt(ZZ%R(i)*ZZ%R(i+1))/AU,ZZ%C(i,ZZ%nt/2,1)%gasdens
	enddo
	close(unit=32)
	
	if(filename(len_trim(filename)-4:len_trim(filename)).eq.'.fits') then
		filename=trim(filename)//'.gz'
	endif

	inquire(file=filename,exist=truefalse)
	if(truefalse) then
		call output("FITS file already exists, overwriting")
		open(unit=90,file=filename)
		close(unit=90,status='delete')
	endif

	status=0
C	 Get an unused Logical Unit Number to use to create the FITS file
	call ftgiou(unit,status)
C	 create the new empty FITS file
	blocksize=1
	call ftinit(unit,filename,blocksize,status)

	simple=.true.
	extend=.true.
	group=1
	fpixel=1

	bitpix=-64
	naxis=4
	naxes(1)=Zone(izone)%nr
	naxes(2)=Zone(izone)%nt
	naxes(3)=Zone(izone)%np
	naxes(4)=3
	naxes(5)=1
	naxes(6)=1
	nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)

	! Write the required header keywords.
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

	! Write optional keywords to the header

	call ftpkyd(unit,'Rin',ZZ%Rin/AU,8,'[AU]',status)
	call ftpkyd(unit,'Rout',ZZ%Rout/AU,8,'[AU]',status)
	call ftpkyd(unit,'x0',ZZ%x0/AU,8,'[AU]',status)
	call ftpkyd(unit,'y0',ZZ%y0/AU,8,'[AU]',status)
	call ftpkyd(unit,'z0',ZZ%z0/AU,8,'[AU]',status)

	call ftpkyj(unit,'nR',ZZ%nr,' ',status)
	call ftpkyj(unit,'nTheta',ZZ%nt,' ',status)
	call ftpkyj(unit,'nPhi',ZZ%np,' ',status)
	call ftpkyj(unit,'npart',npart,' ',status)
	call ftpkyj(unit,'nsize',maxns,' ',status)
	call ftpkyj(unit,'nTemp',maxnT,' ',status)

	call ftpkyj(unit,'nHDU',nvars,' ',status)	
	do i=1,nvars
		write(hdu,'("HDU",i3)') i
		call ftpkys(unit,hdu,trim(vars(i)),'',status)
	enddo

	!  Write the array to the FITS file.

	!------------------------------------------------------------------------------
	! HDU 0: spatial grid
	!------------------------------------------------------------------------------

	allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))

	do i=1,ZZ%nr
		do j=1,ZZ%nt
			do k=1,ZZ%np
				array(i,j,k,1,1,1)=sqrt(ZZ%R(i)*ZZ%R(i+1))/AU
				array(i,j,k,2,1,1)=(ZZ%theta(j)+ZZ%theta(j+1))/2d0
				array(i,j,k,3,1,1)=(ZZ%phi(k)+ZZ%phi(k+1))/2d0
			enddo
		enddo
	enddo

	call ftpprd(unit,group,fpixel,nelements,array,status)
	
	deallocate(array)


	do ivars=1,nvars
		! create new hdu
		call ftcrhd(unit, status)
		bitpix=-64
		naxes(1)=ZZ%nr
		naxes(2)=ZZ%nt
		naxes(3)=ZZ%np
		naxes(4)=1
		naxes(5)=1
		naxes(6)=1
		naxis=3
		select case (vars(ivars))
			case ('RGRID')
				naxis=1
				naxes(1)=ZZ%nr+1
				naxes(2)=1
				naxes(3)=1
				naxes(4)=1
				naxes(5)=1
				naxes(6)=1
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr+1
					array(i,1,1,1,1,1)=ZZ%R(i)
				enddo
			case ('TGRID')
				naxis=1
				naxes(1)=ZZ%nt+1
				naxes(2)=1
				naxes(3)=1
				naxes(4)=1
				naxes(5)=1
				naxes(6)=1
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nt+1
					array(i,1,1,1,1,1)=ZZ%theta(i)
				enddo
			case ('PGRID')
				naxis=1
				naxes(1)=ZZ%np+1
				naxes(2)=1
				naxes(3)=1
				naxes(4)=1
				naxes(5)=1
				naxes(6)=1
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%np+1
					array(i,1,1,1,1,1)=ZZ%phi(i)
				enddo
			case ('DENS')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%dens
						enddo
					enddo
				enddo
			case ('TEMP')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%T
						enddo
					enddo
				enddo
			case ('COMP')
				naxis=6
				naxes(4)=npart
				naxes(5)=maxns
				naxes(6)=maxnT
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				array=0d0
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							do ipart=1,npart
								do is=1,Part(ipart)%nsize
									do iT=1,Part(ipart)%nT
										array(i,j,k,ipart,is,iT)=ZZ%C(i,j,k)%densP(ipart,is,iT)
									enddo
								enddo
							enddo
						enddo
					enddo
				enddo
			case ('GASDENS')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%gasdens
						enddo
					enddo
				enddo
			case ('NPHOT')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%Ni
						enddo
					enddo
				enddo
			case ('VOLUME')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%V
						enddo
					enddo
				enddo
			case ('EJv')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%E
						enddo
					enddo
				enddo
			case ('G0')
				nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
				allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							array(i,j,k,1,1,1)=ZZ%C(i,j,k)%G0
						enddo
					enddo
				enddo
			case default
				call output("Error in output file specification")
				print*,vars(ivars)
				close(unit=20)
				stop
		end select


		!  Write the required header keywords.
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		!  Write the array to the FITS file.
		call ftpprd(unit,group,fpixel,nelements,array,status)

		deallocate(array)
	enddo
	
	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	if (status.gt.0) then
	   print*,'error in export to fits file',status
	end if


	return
	end
	



	subroutine readstruct_fits(filename,vars,nvars,izone)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer nvars,ivars,i,j,k,ipart,is,iT,naxis,nhdu,izone,iread
	character*7 vars(nvars)
	character*500 filename
	logical doalloc,truefalse
	real*8,allocatable :: array(:,:,:,:,:,:)
	integer*4 :: status,stat2,stat3,readwrite,unit,blocksize,nfound,group
	integer*4 :: firstpix,nbuffer,npixels,hdunum,hdutype,ix,iz,ilam
	integer*4 :: istat,stat4,tmp_int,stat5,stat6
	real*8  :: nullval
	logical*4 :: anynull
	integer*4, dimension(6) :: naxes
	character*80 comment,errmessage
	character*30 errtext
	type(ZoneType),pointer :: ZZ

	ZZ => Zone(izone)

	! Get an unused Logical Unit Number to use to open the FITS file.
	status=0

	call ftgiou (unit,status)
	! Open file
	readwrite=0
	call ftopen(unit,filename,readwrite,blocksize,status)
	if (status /= 0) then
		call output("Density file not found")
		print*,trim(filename),status
		call output("--------------------------------------------------------")
		stop
	endif
	group=1
	firstpix=1
	nullval=-999


	!------------------------------------------------------------------------
	! HDU0 : grid
	!------------------------------------------------------------------------
	! Check dimensions
	call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)

	npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

	! Read model info

	call ftgkyd(unit,'Rin',ZZ%Rin,comment,status)
	ZZ%Rin=ZZ%Rin*AU
	call ftgkyd(unit,'Rout',ZZ%Rout,comment,status)
	ZZ%Rout=ZZ%Rout*AU

	call ftgkyj(unit,'nR',ZZ%nr,comment,status)
	call ftgkyj(unit,'nTheta',ZZ%nt,comment,status)
	call ftgkyj(unit,'nPhi',ZZ%np,comment,status)
	call ftgkyj(unit,'npart',npart,comment,status)
	call ftgkyj(unit,'nsize',maxns,comment,status)
	call ftgkyj(unit,'nTemp',maxnT,comment,status)

	call ftgkyj(unit,'nHDU',nhdu,comment,status)
	if(status.ne.0) then
		nhdu=nvars
		status=0
	endif

	do ivars=1,minval((/nhdu,nvars/))
		!  move to next hdu
		call ftmrhd(unit,1,hdutype,status)
		if(status.ne.0) then
			nhdu=ivars
			status=0
			goto 1
		endif
		naxis=3
		select case (vars(ivars))
			case ('RGRID')
				naxis=1
			case ('TGRID')
				naxis=1
			case ('PGRID')
				naxis=1
			case ('COMP')
				naxis=6
		end select

		! Check dimensions
		call ftgknj(unit,'NAXIS',1,naxis,naxes,nfound,status)

		do i=naxis+1,6
			naxes(i)=1
		enddo
		npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)

		! read_image
		allocate(array(naxes(1),naxes(2),naxes(3),naxes(4),naxes(5),naxes(6)))

		call ftgpvd(unit,group,firstpix,npixels,nullval,array,anynull,status)
   
		select case (vars(ivars))
			case ('RGRID')
				do i=1,ZZ%nr+1
					ZZ%R(i)=array(i,1,1,1,1,1)
				enddo
				ZZ%Rin=ZZ%R(1)
				ZZ%Rout=ZZ%R(ZZ%nr+1)
			case ('TGRID')
				do i=1,ZZ%nt+1
					ZZ%theta(i)=array(i,1,1,1,1,1)
				enddo
			case ('PGRID')
				do i=1,ZZ%np+1
					ZZ%phi(i)=array(i,1,1,1,1,1)
				enddo
			case ('DENS')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%dens=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('TEMP')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%T=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('COMP')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%densP=0d0
							do ipart=1,npart
								do is=1,Part(ipart)%nsize
									do iT=1,Part(ipart)%nT
										ZZ%C(i,j,k)%densP(ipart,is,iT)=array(i,j,k,ipart,is,iT)
									enddo
								enddo
							enddo
						enddo
					enddo
				enddo
			case ('GASDENS')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%gasdens=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('NPHOT')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%Ni=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('VOLUME')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%V=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('EJv')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%E=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('G0')
				do i=1,ZZ%nr
					do j=1,ZZ%nt
						do k=1,ZZ%np
							ZZ%C(i,j,k)%G0=array(i,j,k,1,1,1)
						enddo
					enddo
				enddo
			case ('SKIP')
c	just skip this hdu
			case default
				call output("Error in input file specification")
				stop
		end select
		deallocate(array)
	enddo

1	continue

	if(ivars.le.nvars) then
		iread=ivars
		do ivars=iread,nvars
			print*,vars(ivars)
			select case (vars(ivars))
				case ('GASDENS')
					call output("setting gas density according to gas2dust ratio")
					do i=1,ZZ%nr
						do j=1,ZZ%nt
							do k=1,ZZ%np
								ZZ%C(i,j,k)%gasdens=ZZ%C(i,j,k)%dens*ZZ%gas2dust
							enddo
						enddo
					enddo
				case ('NPHOT')
					call output("setting photon statistics to 100")
					do i=1,ZZ%nr
						do j=1,ZZ%nt
							do k=1,ZZ%np
								ZZ%C(i,j,k)%Ni=100
							enddo
						enddo
					enddo
				case ('SKIP')
c	just skip this hdu
				case default
					call output("Error missing HDU")
					stop
			end select
		enddo
	endif

	!  Close the file and free the unit number.
	call ftclos(unit, status)
	call ftfiou(unit, status)

	!  Check for any error, and if so print out error messages
	!  Get the text string which describes the error
	if (status > 0) then
	   call ftgerr(status,errtext)
	   print *,'FITSIO Error Status =',status,': ',errtext

	   !  Read and print out all the error messages on the FITSIO stack
	   call ftgmsg(errmessage)
	   do while (errmessage .ne. ' ')
		  print *,errmessage
		  call ftgmsg(errmessage)
	   end do
	endif

	return
	end


	subroutine SetZoneStructOutput()
	use GlobalSetup
	IMPLICIT NONE
	
	ZoneStructOutput( 1)='RGRID'
	ZoneStructOutput( 2)='TGRID'
	ZoneStructOutput( 3)='PGRID'
	ZoneStructOutput( 4)='DENS'
	ZoneStructOutput( 5)='TEMP'
	ZoneStructOutput( 6)='COMP'
	ZoneStructOutput( 7)='GASDENS'
	ZoneStructOutput( 8)='NPHOT'
	ZoneStructOutput( 9)='VOLUME'
	ZoneStructOutput(10)='EJv'
	ZoneStructOutput(11)='G0'
	nZoneStructOutput=11

	return
	end
	
