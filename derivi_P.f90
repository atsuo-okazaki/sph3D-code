      SUBROUTINE derivi (dt,itime,x,y,z,vx,vy,vz, &
                    u,h,dvx,dvy,dvz,dup,dhp,npart,ntot,ireal)
!************************************************************
!                                                           *
!  This subroutine drives the computation of the forces on  *
!     every particle on the list.                           *
!                                                           *
!************************************************************

      use idims

      use constants
      use force, only: fx, fy, fz, du, dh
      use table
      use tlist
      use btree
      use densi
      use numpa
      use gravi
      use ener
      use kerne
      use divve
      use eosq
      use cgas
!---- neighbor.mod <-- neighbor.f90 for non-parallel simulation
!---- neighbor.mod <-- neighbor_P.f90 for parallel simulation
      use neighbor
      use carac
      use integ
      use typef
      use timei
      use logun
      use debug
      use rbnd
      use phase
      use ptmass
      use current
      use hagain
      use curlist
      use avail
      use perform
      use zzhp

      implicit none

      INTEGER(I4B) :: ireal(idim)
      INTEGER(I4B) :: itime, npart, ntot, nlst_in, nlst_end, i, j, &
                jneigh, ipart, jpart
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: vx(idim), vy(idim), vz(idim), u(idim), h(idim)
      REAL(DP) :: dvx(idim), dvy(idim), dvz(idim), dup(idim), &
                 dhp(idim)
      REAL(DP) :: dumrho(idim),dumpr(idim),dumvs(idim)
      REAL(DP) :: dt, deltat, deltarho
      CHARACTER(len=7) :: where='derivi'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine derivi')
!
!--Set constants first time around
!
      nlst_in = 1
      nlst_end = nlst
      IF (itrace.EQ.'all') WRITE (iprint, 99002) nlst_in, nlst_end
99002 FORMAT(' derivi ',I8, I8)
!
!--Initialise
!
!$OMP PARALLEL DO SCHEDULE (runtime) default(none) &
!$OMP shared(npart,iavail,iphase) &
!$OMP private(i)
      DO i = 1, npart
         IF (iphase(i).NE.-1) iavail(i) = 0
      END DO
!$OMP END PARALLEL DO
!
!--Compute the neighbor indexes & gravitational forces of the distant
!     particles for all the particles in the list
!
      IF (igrape.EQ.0) THEN
         CALL insulate(3, ntot, npart, x, y, z, pmass, h)
      ELSEIF (igrape.EQ.1) THEN
         CALL insulate(4, ntot, npart, x, y, z, pmass, h)
      ELSE
         CALL error(where,1)
      ENDIF
!
!--Find neighbours and calculate gravity on and from point masses
!     The point masses are no longer in the TREE/GRAPE - done separately for
!     higher accuracy
!
      IF (nptmass.GT.0 .AND. (.NOT.(iptintree))) &
           CALL gforspt(x,y,z,vx,vy,vz,h,npart)
!
!--Compute the pressure, divv, etc. on list particles but
!     do not search twice the list particles
!
      jneigh = 0

      IF (itiming) CALL getused(tdens1)

      CALL densityi(jneigh,npart,x,y,z,vx,vy,vz,u,h, &
                     nlst_in,nlst_end,llist,itime)

      IF (itiming) THEN
         CALL getused(tdens2)
         tdens = tdens + (tdens2 - tdens1)
      ENDIF
!
!--Predict the pressure, divv, etc. on list-particle neighbors
!
      DO i = 1, nlst
         ipart = llist(i)
         DO j = 1, ilen(ipart)
            jpart = neighb(j,ipart)
            IF (iscurrent(jpart).EQ.0) THEN
               iavail(jpart) = 1
            ENDIF
         END DO
      END DO
!$OMP PARALLEL default(none) &
!$OMP shared(npart,iavail,it1,imax,itime,it0,isteps,dt,imaxstep) &
!$OMP shared(divv,rho,dumrho,u,dumpr,dumvs,pr,vsound,ntot,ireal) &
!$OMP shared(iphase) &
!$OMP private(i,j,deltat,deltarho)
!$OMP DO SCHEDULE (runtime)
      DO i = 1, npart
         IF (iphase(i).NE.-1) THEN
            IF (iavail(i).EQ.1) THEN
               IF (it1(i).EQ.imax) THEN
                  deltat = dt*(itime - it0(i) - isteps(i)/2)/imaxstep
               ELSE
                  deltat = dt*(itime - it0(i))/imaxstep
               ENDIF
!
!--Update the density value at neighbor's locations
!--Avoid, though, abrupt changes in density
!
               deltarho = -deltat*divv(i)
               IF (ABS(deltarho).GT.rho(i)/2.) THEN
                  deltarho = SIGN(1.0d0,deltarho)*rho(i)/2.0d0
               ENDIF
               dumrho(i) = rho(i) + deltarho

               CALL eospg(i,u,dumrho,dumpr,dumvs)
               iavail(i) = 0
            ELSE
               dumrho(i) = rho(i)
               dumpr(i) = pr(i)
               dumvs(i) = vsound(i)
            ENDIF
         ENDIF
      END DO
!$OMP END DO

!$OMP DO SCHEDULE (runtime)
      DO i = npart + 1, ntot
         j = ireal(i)
         dumrho(i) = dumrho(j)
         dumpr(i) = dumpr(j)
         dumvs(i) = dumvs(j)
         divv(i) = 0.
         u(i) = u(j)
      END DO
!$OMP END DO
!$OMP END PARALLEL
!
!--Compute forces on EACH particle
!
      jneigh = 0

      IF (itiming) CALL getused(tforce1)

      CALL forcei(jneigh,nlst_in,nlst_end,llist,dt,itime,npart, &
           x,y,z,vx,vy,vz,u,h,dvx,dvy,dvz,dhp,dup,dumrho, &
           dumpr,dumvs)

      IF (itiming) THEN
         CALL getused(tforce2)
         tforce = tforce + (tforce2 - tforce1)
      ENDIF

      END SUBROUTINE derivi
