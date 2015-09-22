	subroutine SetupSpecialZone(izone)
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	integer izone,ir,it,ip,i,ips,ipt,it1,it2,ir1,ir2
	real*8 x,y,z,r,random,w(npart,maxns),theta,SD,H
	real*8 Mdot,dphi,t0,dt,vesc,t,drad,dtheta,Mtot,lam0,d,tau
	integer nrings,jj,njj,ilam
	type(ZoneType),pointer :: ZZ

	ZZ => Zone(izone)

c setup:
c C(ir,it,ip)%gasdens
c C(ir,it,ip)%dens
c C(ir,it,ip)%densP(ipart,isize,iT)

	do i=1,npart
		if(Part(i)%nsize.gt.1) then
			call SetSizeDis(w(i,1:Part(i)%nsize),i,izone)
		else
			w(i,1)=1d0
		endif
	enddo

	select case(ZZ%denstype)
		case ("RINGS")
			if(ZZ%shape.eq.'CAR') then
				call output("Rings on a rectangular grid... Let's don't and say we did.")
				stop
			endif

			call SetupRadGrid(izone)
			call SetupThetaGridEqui(izone)
			call SetupPhiGrid(izone)
			call SetupVolume(izone)

			nrings=10
			Mdot=1d-3*Msun/year
			Mdot=ZZ%Mdust/year
			dphi=2d0*pi/(10d0*year)
			t0=5d0*year
			dt=2d0*pi/dphi/real(ZZ%np)
			vesc=10d5
			drad=0.001
			dtheta=pi*0.025
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
				ZZ%C(ir,it,ip)%M=0d0
			enddo
			enddo
			enddo

			do i=1,nrings
				call randomdirection(x,y,z)
				r=random(idum)*(ZZ%Rout-ZZ%Rin)+ZZ%Rin
				theta=acos(z/sqrt(x*x+y*y+z*z))
				do it1=1,ZZ%nt-1
					if((theta-dtheta).gt.ZZ%theta(it1).and.(theta-dtheta).lt.ZZ%theta(it1+1)) exit
				enddo
				if(it1.gt.ZZ%nt) it1=ZZ%nt
				do it2=1,ZZ%nt-1
					if((theta+dtheta).gt.ZZ%theta(it2).and.(theta+dtheta).lt.ZZ%theta(it2+1)) exit
				enddo
				if(it2.gt.ZZ%nt) it2=ZZ%nt
				ip=random(idum)*ZZ%np+1
				t=0d0
				do while(t.lt.t0.and.r.lt.ZZ%Rout)
					do ir1=1,ZZ%nr-1
						if((r-drad*r).gt.ZZ%R(ir1).and.(r-drad*r).lt.ZZ%R(ir1+1)) exit
					enddo
					if(ir1.gt.ZZ%nr) ir1=ZZ%nr
					do ir2=1,ZZ%nr-1
						if((r+drad*r).gt.ZZ%R(ir2).and.(r+drad*r).lt.ZZ%R(ir2+1)) exit
					enddo
					if(ir2.gt.ZZ%nr) ir2=ZZ%nr
					do it=it1,it2
					do ir=ir1,ir2
						ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%M+Mdot*dt/(real(abs(it2-it1)+1))/(real(abs(ir2-ir1)+1))
					enddo
					enddo
					t=t+dt
					ip=ip+1
					r=r+vesc*dt
					if(ip.gt.ZZ%np) ip=1
				enddo
			enddo
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
				ZZ%C(ir,it,ip)%dens=ZZ%C(ir,it,ip)%M/ZZ%C(ir,it,ip)%V
				if(ZZ%C(ir,it,ip)%dens.lt.1d-30) then
					ZZ%C(ir,it,ip)%dens=1d-30
					ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%dens*ZZ%C(ir,it,ip)%V
				endif
				do i=1,npart
					do ips=1,Part(i)%nsize
						ZZ%C(ir,it,ip)%densP(i,ips,1)=w(i,ips)*ZZ%abun(i)*ZZ%C(ir,it,ip)%dens
						do ipt=2,Part(i)%nT
							ZZ%C(ir,it,ip)%densP(i,ips,ipt)=0d0
						enddo
					enddo
				enddo
				ZZ%C(ir,it,ip)%gasdens=ZZ%C(ir,it,ip)%dens*ZZ%gas2dust
			enddo
			enddo
			enddo
		case ("BETAPIC")
			if(ZZ%shape.eq.'CAR') then
				call output("Beta Pic on a rectangular grid... Let's don't and say we did.")
				stop
			endif

			call SetupRadGrid(izone)
			call SetupThetaGrid(izone)
			call SetupPhiGrid(izone)
			call SetupVolume(izone)

			lam0=0.65
			d=lam(nlam)-lam(1)
			ilam=1
			do i=1,nlam
				if(abs(lam0-lam(i)).lt.d) then
					d=abs(lam0-lam(i))
					ilam=i
				endif
			enddo
			tau=0d0
			do i=1,npart
				do ips=1,Part(i)%nsize
					tau=tau+w(i,ips)*ZZ%abun(i)*Part(i)%Ksca(ips,1,ilam)
				enddo
			enddo

			Mtot=0d0
			njj=10
			do ir=1,ZZ%nr
			r=sqrt(ZZ%R(ir)*ZZ%R(ir+1))
			do ip=1,ZZ%np
			do it=1,ZZ%nt
				ZZ%C(ir,it,ip)%dens=0d0
				SD=(ZZ%bp_eta/tau)/sqrt((r/ZZ%Rexp)**(-2d0*ZZ%bp_alpha)+(r/ZZ%Rexp)**(-2d0*ZZ%bp_beta))
				H=ZZ%bp_A*ZZ%Rexp*(r/ZZ%Rexp)**ZZ%shpow
				do jj=1,njj
					theta=ZZ%theta(it)+(ZZ%theta(it+1)-ZZ%theta(it))*real(jj)/real(njj+1)
					z=abs(r*cos(theta))
					ZZ%C(ir,it,ip)%dens=ZZ%C(ir,it,ip)%dens+(SD/H)*exp(-(abs(z)/H)**ZZ%bp_p)
				enddo
				ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%dens*ZZ%C(ir,it,ip)%V
				Mtot=Mtot+ZZ%C(ir,it,ip)%M
			enddo
			enddo
			enddo
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
c				ZZ%C(ir,it,ip)%dens=ZZ%C(ir,it,ip)%dens*ZZ%Mdust/Mtot
c				ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%M*ZZ%Mdust/Mtot
				if(ZZ%C(ir,it,ip)%dens.lt.1d-30) then
					ZZ%C(ir,it,ip)%dens=1d-30
					ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%dens*ZZ%C(ir,it,ip)%V
				endif
			enddo
			enddo
			enddo
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
				do i=1,npart
					do ips=1,Part(i)%nsize
						ZZ%C(ir,it,ip)%densP(i,ips,1)=w(i,ips)*ZZ%abun(i)*ZZ%C(ir,it,ip)%dens
						do ipt=2,Part(i)%nT
							ZZ%C(ir,it,ip)%densP(i,ips,ipt)=0d0
						enddo
					enddo
				enddo
				ZZ%C(ir,it,ip)%gasdens=ZZ%C(ir,it,ip)%dens*ZZ%gas2dust
			enddo
			enddo
			enddo
		case ("PLANETRING")
			if(ZZ%shape.eq.'CAR') then
				call output("Rings on a rectangular grid... Let's don't and say we did.")
				stop
			endif

			call SetupRadGrid(izone)
			call SetupThetaGrid(izone)
			call SetupPhiGrid(izone)
			call SetupVolume(izone)

			tau=0d0
			do i=1,npart
				do ips=1,Part(i)%nsize
					tau=tau+w(i,ips)*ZZ%abun(i)*3d0/(4d0*Part(i)%rv(ips)*Part(i)%rho(1))
				enddo
			enddo

			Mtot=0d0
			do ir=1,ZZ%nr
			r=sqrt(ZZ%R(ir)*ZZ%R(ir+1))
			do ip=1,ZZ%np
				do it=1,ZZ%nt
					ZZ%C(ir,it,ip)%dens=0d0
				enddo
				it=ZZ%imidplane
				ZZ%C(ir,it,ip)%dens=ZZ%tau_V*(pi*(ZZ%R(ir+1)**2-ZZ%R(ir)**2)/real(ZZ%np))/(tau*ZZ%C(ir,it,ip)%V)
				ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%dens*ZZ%C(ir,it,ip)%V
				Mtot=Mtot+ZZ%C(ir,it,ip)%M
			enddo
			enddo
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
				if(ZZ%C(ir,it,ip)%dens.lt.1d-30) then
					ZZ%C(ir,it,ip)%dens=1d-30
					ZZ%C(ir,it,ip)%M=ZZ%C(ir,it,ip)%dens*ZZ%C(ir,it,ip)%V
				endif
			enddo
			enddo
			enddo
			do ir=1,ZZ%nr
			do it=1,ZZ%nt
			do ip=1,ZZ%np
				do i=1,npart
					do ips=1,Part(i)%nsize
						ZZ%C(ir,it,ip)%densP(i,ips,1)=w(i,ips)*ZZ%abun(i)*ZZ%C(ir,it,ip)%dens
						do ipt=2,Part(i)%nT
							ZZ%C(ir,it,ip)%densP(i,ips,ipt)=0d0
						enddo
					enddo
				enddo
				ZZ%C(ir,it,ip)%gasdens=ZZ%C(ir,it,ip)%dens*ZZ%gas2dust
			enddo
			enddo
			enddo
		case default
			call output("Really? A " // trim(ZZ%denstype) // "-zone? That is an awesome idea! (but not yet possible)")
			stop
	end select
	
	return
	end


