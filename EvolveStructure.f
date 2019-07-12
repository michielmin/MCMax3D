	subroutine EvolveStructure(iter)
	use GlobalSEtup
	IMPLICIT NONE
	integer iter
	
	if(maxnT.gt.1) then
		if(iter.gt.0) call createUV()
		call Topac(iter)
	endif
	
	return
	end
	
	

c-----------------------------------------------------------------------
c  This routine calculates the temperature dependent opacities.
c  
c
c-----------------------------------------------------------------------


      subroutine Topac(niter)
      use GlobalSetup
      implicit none
      integer i,is,iT,niter,izone,i1,i2,i3
      real*8 dens,A,B,dens0,T,w1,w2
      type(Cell),pointer :: C

      ! loop over all grains
      do i=1,npart
         if (Part(i)%nT.gt.1) then
			call output("Calculating opacities for particle " // int2string(i,'(i4)'))

			! loop over all cells
			do is=1,Part(i)%nsize
			call tellertje(is,Part(i)%nsize)
			do izone=1,nzones
!$OMP PARALLEL IF(.true.)
!$OMP& DEFAULT(NONE)
!$OMP& PRIVATE(i1,i2,i3,C,dens,iT,T,A,B,dens0,w1,w2)
!$OMP& SHARED(Zone,izone,i,is,Part,gammaUVdes,niter)
!$OMP DO
!$OMP& SCHEDULE(DYNAMIC, 1)
			do i1=1,Zone(izone)%n1
				do i2=1,Zone(izone)%n2
				do i3=1,Zone(izone)%n3
					C => Zone(izone)%C(i1,i2,i3)
                  ! calc opacity using cell temperature
					dens=0d0
					do iT=1,Part(i)%nT
						dens=dens+C%densP(i,is,iT)
						C%densP(i,is,iT)=0d0
					enddo
					if(niter.eq.0) then
						T=2.7d0
					else
						T=C%T
					endif
					if(T.lt.2.7d0) T=2.7d0
					if(T.le.Part(i)%Tmax(1)) then
						C%densP(i,is,1)=dens
					else if(T.ge.Part(i)%Tmax(Part(i)%nT)) then
						C%densP(i,is,Part(i)%nT)=dens
					else
						do iT=1,Part(i)%nT-1
							if(T.gt.Part(i)%Tmax(iT).and.T.le.Part(i)%Tmax(iT+1)) then
								w1=(Part(i)%Tmax(iT+1)-T)/(Part(i)%Tmax(iT+1)-Part(i)%Tmax(iT))
								w2=1d0-w1
								C%densP(i,is,iT)=w1*dens
								C%densP(i,is,iT+1)=w2*dens
							endif
						enddo
					endif
					if(niter.ne.0) then
						A=Part(i)%TdesA
						B=Part(i)%TdesB
						dens0=10d0**(A/B-1d4/(T*B)-log10(T))+gammaUVdes*C%G0/sqrt(T)
						if(C%dens.lt.dens0) then
							C%densP(i,is,1:Part(i)%nT-1)=0d0
							C%densP(i,is,Part(i)%nT)=dens
						endif
					endif
				enddo
				enddo
			enddo
!$OMP END DO
!$OMP FLUSH
!$OMP END PARALLEL
			enddo
			enddo
		endif
	enddo


	return
	end


	subroutine createUV()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer NphotMono,izone,i1,i2,i3,iter,ilam
	type(Cell),pointer :: C
	
	NphotMono=250000
	
	do izone=1,nzones
		do i1=1,Zone(izone)%n1
		do i2=1,Zone(izone)%n2
		do i3=1,Zone(izone)%n3
			C => Zone(izone)%C(i1,i2,i3)
			C%G0=0d0
		enddo
		enddo
		enddo
	enddo
			
	do ilam=1,nlam
		if(lam(ilam).gt.0.0912.and.lam(ilam).le.0.205) then
			call TraceScattField(0,ilam,NphotMono)

			do izone=1,nzones
			do i1=1,Zone(izone)%n1
			do i2=1,Zone(izone)%n2
			do i3=1,Zone(izone)%n3
				C => Zone(izone)%C(i1,i2,i3)
				C%G0=C%G0+C%Escatt*dnu(ilam)
			enddo
			enddo
			enddo
			enddo
		endif
	enddo

	do izone=1,nzones
		do i1=1,Zone(izone)%n1
		do i2=1,Zone(izone)%n2
		do i3=1,Zone(izone)%n3
			C => Zone(izone)%C(i1,i2,i3)
			C%G0=C%G0/5.33d-14/3d10/C%V
		enddo
		enddo
		enddo
	enddo


	
	return
	end


	subroutine createAV()
c	Calculates the visual extinction in a simple way.
c     AV is usefull for doing postprocessing chemistry.
c     currently only the radial visual extinction (from the star outward)
c     is calculated
c
c     Author: Rab C.
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,ilam
	real*8 :: tau,dr,drm1
	real*8 :: GetKext
	real*8, parameter :: ref_lam=0.55d0
	real*8, parameter :: avfac=2.5*log10(exp(1.0))
	integer ir,it,ip


c	find the wavelength index for (optical used for AV)
c     do a simple nearest neighbour interpolation
	do ilam=1,nlam-1
		if(ref_lam.gt.lam(ilam).and.ref_lam.le.lam(ilam+1)) then
			exit
		endif
	enddo


	do izone=1,nzones
		do ir=1,Zone(izone)%nr
c                 cell center radius as defined in subroutine outputstruct_fits
			dr=sqrt(Zone(izone)%R(ir)*Zone(izone)%R(ir+1))-Zone(izone)%R(ir)
			if (ir.eq.1) then
				drm1=0
			else
				drm1=Zone(izone)%R(ir)-sqrt(Zone(izone)%R(ir-1)*Zone(izone)%R(ir))
			endif
			do it=1,Zone(izone)%nt
			do ip=1,Zone(izone)%np
			Zone(izone)%C(ir,it,ip)%AVrad=0.0
			tau=GetKext(ilam,Zone(izone)%C(ir,it,ip))*dr
			if (ir.eq.1) then
				Zone(izone)%C(ir,it,ip)%AVrad=avfac*tau
			else
				tau=tau+GetKext(ilam,Zone(izone)%C(ir-1,it,ip))*drm1
				Zone(izone)%C(ir,it,ip)%AVrad=Zone(izone)%C(ir-1,it,ip)%AVrad+tau*avfac
			endif
			!write(*,*) ir,it,ip,"tau",tau

			enddo
			enddo
		enddo
	enddo

	write(*,*) "tauRad",Zone(1)%C(Zone(1)%nr-1,Zone(1)%imidplane,1)%AVrad

	end subroutine
