      SUBROUTINE psplit
!  ************************************************************
!  *                                                          *
!  *   This subroutine splits particles near an accretor      *
!  *   in the case of a binary.                               *
!  *                                                          *
!  *   When an SPH particle enters the Roche radius of        *
!  *   the accretor, it is replacesd by 13 child particles    *
!  *   with masses pmass(child)=pmass(parent)/13 and          *
!  *   smoothing lengths h(child)=h(parent)/13^{1/3}.         *
!  *   One of these child particles is placed at the position *
!  *   of the parent particle, and the other 12 children are  *
!  *   positioned at the vertices of an hexagonal closed      *
!  *   packed array, at an equal distance disk=1.5*h(child)   *
!  *   from the central child particle. For details, see      *
!  *   Kitsionas, Whitworth 2002, MNRAS, 330, 129 and         *
!  *   Kitsionas, Whitworth 2007, MNRAS, 378, 507.            *
!  *                                                          *
!  *   9 June 2012: Created.                                  *
!  *                                                          *
!  ************************************************************

      use mpi_mod 
      use idims

      use constants 
      use units 
      use logun 
      use part 
      use densi 
      use typef 
      use carac 
      use cgas 
      use kerne 
      use ener 
      use divve 
      use eosq 
      use bodys 
      use force, only: f1vx, f1vy, f1vz, f1u, f1h, &
                       f2vx, f2vy, f2vz, f2u, f2h
      use btree 
      use fracg 
      use polyk2 
      use phase 
      use ptmass 
      use nextmpt 
      use tlist 
      use ptdump 
      use current
      use curlist 
      use timei 
      use init 
      use ghost 
      use table 
      use active 
      use accnum 
      use varet 
      use binary 
      use accurpt 
      use delay 
      use neighbor 
      use new 
      use sphcom 
      use crpart 
      use misali 
      use capt 
      use winds 
      use split 
                                                                        
      implicit none

      INTEGER(I4B), parameter :: nchild = 13
      INTEGER(I4B) :: lstsplit(idim)
      INTEGER(I4B) :: iremove(idim) 
      INTEGER(I4B) :: nsplit, i, ii, iii, ip, j, &
               it0p, it1p, istepsp, iphasep, &
               iantigrp, ibelongp
      INTEGER(I4B) :: i1, nkillcur, iadd, nadd
      INTEGER(I4B), allocatable :: listkilled(:)
      REAL(DP), save :: xi, yi, zi, rx1, ry1, rz1, r2, &
               d12, q13, q23, roche, roche2, hroche2, cr13, &
               sqrt23, sqrt3, sixth, third, half, dist, pmassc, hc, &
               xp, yp, zp, vxp, vyp, vzp, pmassp, hp, &
               f1vxp, f1vyp, f1vzp, f1hp, f1up, f2vxp, &
               f2vyp, f2vzp, f2hp, f2up, dum2vxp, dum2vyp, &
               dum2vzp, up, rhop, dphitp, dgravp, potenp, &
               prp, divvp
      REAL(DP), allocatable :: dummx(:)
      REAL(DP), allocatable :: dummy(:)
      REAL(DP), allocatable :: dummz(:)
      REAL(DP), allocatable :: dummvx(:)
      REAL(DP), allocatable :: dummvy(:)
      REAL(DP), allocatable :: dummvz(:)
      REAL(DP), allocatable :: dummrho(:)
      REAL(DP), allocatable :: dummu(:)
      REAL(DP), allocatable :: dummh(:)
      REAL(DP), allocatable :: dummpmass(:)
      CHARACTER(len=7), parameter ::  where='psplit'
                                                                        
!
!--Split particles near the accretor (i<>isphcom)
!--   Store gohost particles' data to temporary arrays
      IF (nghost /= 0) THEN
         allocate (dummx(nghost)) 
         allocate (dummy(nghost)) 
         allocate (dummz(nghost)) 
         allocate (dummvx(nghost)) 
         allocate (dummvy(nghost)) 
         allocate (dummvz(nghost)) 
         allocate (dummrho(nghost)) 
         allocate (dummu(nghost)) 
         allocate (dummh(nghost)) 
         allocate (dummpmass(nghost)) 
         DO i=npart+1,npart+nghost
            i1 = i-npart
            dummx(i1) = x(i) 
            dummy(i1) = y(i) 
            dummvx(i1) = vx(i) 
            dummvy(i1) = vy(i) 
            dummvz(i1) = vz(i) 
            dummrho(i1) = rho(i) 
            dummu(i1) = u(i) 
            dummh(i1) = h(i) 
            dummpmass(i1) = pmass(i) 
         ENDDO
      ENDIF

!--   List up killed particles including with other timesteps
      nkillcur = COUNT( iphase(1:npart) == -1 )
      IF (nkillcur > 0) THEN
         allocate (listkilled(nkillcur))

         ii = 0
         DO i=1,npart
            IF (iphase(i) == -1) THEN
               ii = ii+1
               listkilled(ii) = i
            ENDIF
         ENDDO
      ENDIF

!---- Only allow splitting of GAS particles avaluated
!     at the CURRENT time step
!                                                                       
      iremove = -1
      DO i = 1, nlst0
         j = llist(i)
         IF (iphase(j) == 0 .AND. j <= npart) iremove(j) = 0
      END DO
!
!--   Check individual particles at the current timestep
!     whether they newly entered the Roche radius (or half that)
!     of the accretor.
!
      nsplit = 0
      DO iii = 1, nptmass
         i = listpm(iii)
         IF (i == isphcom) CYCLE

         xi = x(i)
         yi = y(i)
         zi = z(i)

!----    Calculation of the Roche lobe radius to test whether           
!        particles enter the Roche radii of the companion.              
         q13 =(pmass(listpm(2))/pmass(listpm(1)))**(1.0/3.0) 
         q23 = q13*q13 
         d12 = SQRT((x(listpm(2))-x(listpm(1)))**2 &
                +(y(listpm(2))-y(listpm(1)))**2 &
                +(z(listpm(2))-z(listpm(1)))**2)

!----    Approximate formula for Roche radius (Eggleton 1983)
         roche   = d12*0.49*q23/(0.6*q23+LOG(1.0+q13)) 
         roche2  = roche*roche 
         hroche2 = roche2*0.25
                                                                        
         DO j = 1, npart
            IF (iremove(j) /= 0) CYCLE 
                                                                        
            rx1 = x(j)-xi 
            ry1 = y(j)-yi 
            rz1 = z(j)-zi 
            r2 = rx1*rx1+ry1*ry1+rz1*rz1 
                                                                        
!----       Is the point-mass neighbour inside theRoche radius?         
            IF (r2 >= roche2) CYCLE
!            IF (r2 >= hroche2) CYCLE

            IF (ibelong(j) >= 0) THEN 
               ibelong(j) = -i 
               nsplit = nsplit+1 
               lstsplit(nsplit) = j 
!            ELSE
!               IF (r2 < hroche2 .AND. ibelong(j) == -i) THEN 
!                  ibelong(j) = -i*11
!                  nsplit = nsplit+1 
!                  lstsplit(nsplit) = j 
!               ENDIF
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
      half = 0.5d0 
      third = 1.0d0/3.0d0 
      sixth = half*third 
      sqrt3 = SQRT(3.0d0) 
      sqrt23 = SQRT(2.0d0/3.0d0) 
      cr13 = DBLE(nchild)**(-third)
      iadd = 0
      DO ii = 1, nsplit 
!----    Parent's quantities                                            
         ip = lstsplit(ii) 
         xp = x (ip)
         yp = y (ip)
         zp = z (ip)
         hp = h (ip)
         pmassp = pmass (ip)
         vxp = vx (ip)
         vyp = vy (ip)
         vzp = vz (ip)
         f1vxp = f1vx (ip)
         f1vyp = f1vy (ip)
         f1vzp = f1vz (ip)
         f1hp = f1h (ip)
         f1up = f1u (ip)
         f2vxp = f2vx (ip)
         f2vyp = f2vy (ip)
         f2vzp = f2vz (ip)
         f2hp = f2h (ip)
         f2up = f2u (ip)
         dum2vxp = dum2vx (ip)
         dum2vyp = dum2vy (ip)
         dum2vzp = dum2vz (ip)
         up = u (ip)
         rhop = rho (ip)
         dphitp = dphit (ip)
         dgravp = dgrav (ip)
         potenp = poten (ip)
         prp = pr (ip)
         divvp = divv (ip)
         it0p = it0 (ip)
         it1p = it1 (ip)
         istepsp = isteps (ip)
         iphasep = iphase(ip) 
         iantigrp = iantigr(ip) 
         ibelongp = ibelong(ip) 
                                                                        
!----    Each child's quantities                                             
         hc = hp*cr13 
         pmassc = pmassp/nchild 
         dist = 1.5d0*hc 
                                                                        
         h (ip) = hc 
         pmass(ip) = pmassc 
                                                                        
         DO i = 2, nchild
            iadd = iadd+1
            IF (iadd <= nkillcur) THEN
               j = listkilled(iadd)
            ELSE
               j = npart+iadd-nkillcur
            ENDIF
                                                                        
            SELECTCASE(i) 
            CASE(2) 
               x(j) = xp+dist 
               y(j) = yp 
               z(j) = zp 
            CASE(3) 
               x(j) = xp-dist 
               y(j) = yp 
               z(j) = zp 
            CASE(4) 
               x(j) = xp+half*dist 
               y(j) = yp+half*sqrt3*dist 
               z(j) = zp 
            CASE(5) 
               x(j) = xp+half*dist 
               y(j) = yp-half*sqrt3*dist 
               z(j) = zp 
            CASE(6) 
               x(j) = xp-half*dist 
               y(j) = yp+half*sqrt3*dist 
               z(j) = zp 
            CASE(7) 
               x(j) = xp-half*dist 
               y(j) = yp-half*sqrt3*dist 
               z(j) = zp 
            CASE(8) 
               x(j) = xp 
               y(j) = yp+third*sqrt3*dist 
               z(j) = zp-sqrt23*dist 
            CASE(9) 
               x(j) = xp-half*dist 
               y(j) = yp-sixth*sqrt3*dist 
               z(j) = zp-sqrt23*dist 
            CASE(10) 
               x(j) = xp+half*dist 
               y(j) = yp-sixth*sqrt3*dist 
               z(j) = zp-sqrt23*dist 
            CASE(11) 
               x(j) = xp 
               y(j) = yp-third*sqrt3*dist 
               z(j) = zp+sqrt23*dist 
            CASE(12) 
               x(j) = xp-half*dist 
               y(j) = yp+sixth*sqrt3*dist 
               z(j) = zp+sqrt23*dist 
            CASE(13) 
               x(j) = xp+half*dist 
               y(j) = yp+sixth*sqrt3*dist 
               z(j) = zp+sqrt23*dist 
            ENDSELECT 
                                                                        
            h (j) = hc
            pmass (j) = pmassc

            vx(j) = vxp 
            vy(j) = vyp 
            vz(j) = vzp 
            f1vx(j) = f1vxp
            f1vy(j) = f1vyp
            f1vz(j) = f1vzp
            f1h(j) = f1hp
            f1u(j) = f1up
            f2vx(j) = f2vxp
            f2vy(j) = f2vyp
            f2vz(j) = f2vzp
            f2h(j) = f2hp
            f2u(j) = f2up
            dum2vx(j) = dum2vxp
            dum2vy(j) = dum2vyp
            dum2vz(j) = dum2vzp
            u(j) = up
            rho(j) = rhop
            dphit(j) = dphitp
            dgrav(j) = dgravp
            poten(j) = potenp
            pr(j) = prp
            divv(j) = divvp
            it0(j) = it0p
            it1(j) = it1p
            isteps(j) = istepsp
            iphase(j) = iphasep
            iantigr(j) = iantigrp
            ibelong(j) = ibelongp
         ENDDO 
      ENDDO

      IF (nkillcur > 0) THEN
         deallocate (listkilled)
      ENDIF
      nadd = MAX(iadd-nkillcur, 0)

      IF (nghost /= 0) THEN
         DO i=1,nghost
            i1 = npart+nadd+i
            x(i1) = dummx(i)
            y(i1) = dummx(i)
            z(i1) = dummy(i)
            vx(i1) = dummz(i)
            vy(i1) = dummvx(i)
            vz(i1) = dummvy(i)
            rho(i1) = dummrho(i)
            u(i1) = dummu(i)
            h(i1) = dummh(i)
            pmass(i1) = dummpmass(i)
         ENDDO
         deallocate (dummx) 
         deallocate (dummy) 
         deallocate (dummz) 
         deallocate (dummvx) 
         deallocate (dummvy) 
         deallocate (dummvz) 
         deallocate (dummrho) 
         deallocate (dummu) 
         deallocate (dummh) 
         deallocate (dummpmass)
      ENDIF
                                                                        
!---- Never count split particles as reassigned particles!
      numsplit = numsplit+nsplit 
      npart = npart+nadd
      nactive = nactive+iadd
                                                                        
      END SUBROUTINE psplit
