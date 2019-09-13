	subroutine TraceToOutside()
	use GlobalSetup
	use Constants
	IMPLICIT NONE
	type(Photon) phot
	type(Path) :: P
	integer istar,maxcount
	logical justcount,error

	integer izone,imin,iobs,k,status,i1,i2,i3,iz
	logical leave,inany,hitstar0
	real*8 minv,R,theta,phi,random
	type(Travel) Trac(nzones),TracStar(nstars)
	type(Cell),pointer :: C

	do iz=1,nzones
		do i1=1,Zone(iz)%n1
		do i2=1,Zone(iz)%n2
		do i3=1,Zone(iz)%n3

			error=.false.

c first position photon somewhere in the cell
			select case(Zone(iz)%shape)
				case("SPH")
					R=Zone(iz)%R(i1)+(Zone(iz)%R(i1+1)-Zone(iz)%R(i1))*random(idum)
					theta=Zone(iz)%theta(i2)+(Zone(iz)%theta(i2+1)-Zone(iz)%theta(i2))*random(idum)
					phi=Zone(iz)%phi(i3)+(Zone(iz)%phi(i3+1)-Zone(iz)%phi(i3))*random(idum)
					phot%x=R*cos(phi)*sin(theta)
					phot%y=R*sin(phi)*sin(theta)
					phot%z=R*cos(theta)
				case default
					call output("Raytracing on non-spherical zones not yet possible")
					stop
			end select

c next define the direction of the photon (you need to adjust this to your needs)
			phot%vx=1d0
			phot%vy=0d0
			phot%vz=0d0

c Translate the photon position from the zone system to the general system
			call TranslatePhotonXinverse(phot,iz)

			phot%edgeNr=0
			phot%inzone=.false.

c Compute the position of the photon in all zone coordinate systems
			call TranslatePhotonX(phot)
			call TranslatePhotonV(phot)

c Determine in which zones the photon is located
			call InWhichZones(phot)

			do izone=1,nzones
				Trac(izone)%recompute=.true.
			enddo
			do istar=1,nstars
				TracStar(istar)%recompute=.true.
			enddo
			minv=0d0
	
			k=0
c Start the loop over the path
1			continue

c First determine which cell boundaries to hit
			status=0
			do izone=1,nzones
				if(phot%inzone(izone)) then
					select case(Zone(izone)%shape)
						case("SPH")
c							if(Trac(izone)%recompute) then
								call TravelSph(phot,izone,Trac(izone),status)
c							else
c								Trac(izone)%v=Trac(izone)%v-minv
c							endif
					end select
					Trac(izone)%recompute=.false.
				else
					select case(Zone(izone)%shape)
						case("SPH")
c							if(Trac(izone)%recompute) then
								call HitSph(phot,izone,Trac(izone))
c							else
c								Trac(izone)%v=Trac(izone)%v-minv
c							endif
					end select
					Trac(izone)%recompute=.false.
				endif
			enddo
c Determine which stars are hit
			do istar=1,nstars
c				if(TracStar(istar)%recompute) then
					call HitStar(phot,istar,TracStar(istar))
					TracStar(istar)%recompute=.false.
c				else
c					TracStar(istar)%v=TracStar(istar)%v-minv
c				endif
			enddo


			if(status.gt.0) then
				call output("Something is wrong... Don't worry I'll try to fix it.")
				nerrors=nerrors+1
				if(nerrors.gt.100) stop
				do izone=1,nzones
					Trac(izone)%recompute=.true.
				enddo
				do istar=1,nstars
					TracStar(istar)%recompute=.true.
				enddo
				minv=0d0
				call InWhichZones(phot)
				goto 1
			endif

c Determine what is hit first.
			minv=1d8*maxR
			leave=.true.
			do izone=1,nzones
				if(Trac(izone)%v.gt.0d0.and.Trac(izone)%v.lt.minv) then
					minv=Trac(izone)%v
					imin=izone
					leave=.false.
				endif
			enddo
			hitstar0=.false.
			do istar=1,nstars
				if(TracStar(istar)%v.gt.0d0.and.TracStar(istar)%v.lt.minv) then
					minv=TracStar(istar)%v
					imin=0
					hitstar0=.true.
				endif
			enddo

c Now you can store the path traveled in this particular cell
			do izone=1,nzones
				if(P%inzone(k,izone)) then
					C => Zone(izone)%C(phot%i1(izone),phot%i2(izone),phot%i3(izone))
c					Now store whatever you want:
c					- Travel distance: minv
c					- density: C%dens
c					- zone number: izone
c					- cell number: phot%i1(izone),phot%i2(izone),phot%i3(izone)
c						-> column density contribution minv*C%dens
				endif
			enddo

c Move the photon to the next cell boundary
			if(leave) goto 3
		
			call TravelPhotonX(phot,minv)

			if(hitstar0) goto 3

			do izone=1,nzones
				if(Trac(izone)%v.le.minv.or.izone.eq.imin) then
					if(Zone(izone)%warped.and.phot%i1(izone).ne.Trac(izone)%i1next) then
						phot%i1(izone)=Trac(izone)%i1next
						call TranslatePhotonWarp(phot,izone)
						phot%i2(izone)=-1
						phot%i3(izone)=-1
					else
						phot%i1(izone)=Trac(izone)%i1next
						phot%i2(izone)=Trac(izone)%i2next
						phot%i3(izone)=Trac(izone)%i3next
					endif

					if(phot%i1(izone).eq.-1) call determine_i1(phot,izone)
					if(phot%i2(izone).eq.-1) call determine_i2(phot,izone)
					if(phot%i3(izone).eq.-1) call determine_i3(phot,izone)

					phot%edgeNr(izone)=Trac(izone)%edgenext
					phot%inzone(izone)=.true.
					if(phot%i1(izone).gt.Zone(izone)%n1.or.
     &			phot%i2(izone).gt.Zone(izone)%n2.or.
     &			phot%i3(izone).gt.Zone(izone)%n3.or.
     &			phot%i1(izone).lt.1.or.
     &			phot%i2(izone).lt.1.or.
     &			phot%i3(izone).lt.1) then
						phot%inzone(izone)=.false.
					endif
					Trac(izone)%recompute=.true.
				else
					phot%edgeNr(izone)=0
				endif
			enddo

c And continue with next cell when the photon is not yet escaped.
			goto 1
3			continue

		enddo
		enddo
		enddo
	enddo

	return
	end


