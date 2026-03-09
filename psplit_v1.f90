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
      use binfile 
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

      INTEGER(I4B), parameter :: nchild = 13, nspmax = idim/nchild
      INTEGER(I4B) :: lstsplit(nspmax) 
      INTEGER(I4B) :: nsplit, i, ii, iii, ip, j, jj, nchild1, ishift, &
               isp, ibp, iap, iphasep
      REAL(DP), save :: xi, yi, zi, rx1, ry1, rz1, r2, &
               d12, q13, q23, roche, roche2, cr13, sqrt23, &
               xp, yp, zp, vxp, vyp, vzp, sqrt3, sixth, third, &
               half, rhop, dist, pmassc, hc, pmassp, hp, up, &
               potenp, dgravp 
      CHARACTER(len=7), parameter ::  where='psplit'
!!    LOGICAL, parameter :: ifirst=.true.
                                                                        
      nchild1 = nchild-1 
                                                                        
!--Only allow accretion of GAS particles evaluated at the CURRENT timest
!                                                                       
!!cc      DO i = 1, npart + nghost                                      
!!      DO i = 1, npart                                                 
!!         iremove(i) = -1                                              
!!      END DO                                                          
!!      DO i = 1, nlst0                                                 
!!         j = llist(i)                                                 
!!cc         IF (iphase(j) == 0) iremove(j) = 0                         
!!         IF (iphase(j) == 0 .AND. j <= npart) iremove(j) = 0          
!!      END DO                                                          
!                                                                       
!--Split particles near the accretor (i<>isphcom)                       
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
         q13 =(pmass(listpm(2))/pmass(listpm(1))) **(1.0/3.0) 
         q23 = q13*q13 
         d12 = SQRT((x(listpm(2))-x(listpm(1))) *(x(listpm(2)) &
            -x(listpm(1)))+(y(listpm(2))-y(listpm(1)))      &
            *(y(listpm(2))-y(listpm(1)))+(z(listpm(2)) &
            -z(listpm(1))) *(z(listpm(2))-z(listpm(1))))
!----    Eggleton (1983)                                                
         roche = d12*0.49*q23 /(0.6*q23+LOG(1.0+q13)) 
         roche2 = roche*roche 
                                                                        
!!         if (myrank.eq.0) then                                        
!!            WRITE (iprint,*) 'nptlist(',iii,')=',nptlist(iii)         
!!         endif                                                        
         DO jj = 1, nptlist(iii) 
            j = nearpt(iii,jj) 
            IF (iphase(j) /= 0.OR.j > npart) CYCLE 
                                                                        
            rx1 = x(j)-xi 
            ry1 = y(j)-yi 
            rz1 = z(j)-zi 
            r2 = rx1*rx1+ry1*ry1+rz1*rz1 
                                                                        
!--         Is the point-mass neighbour inside theRoche radius?         
!!            if (myrank.eq.0) then                                     
!!               WRITE (iprint,*) 'r2=',r2,', roche2=',roche2, &
!!                     ', ibelong(',j,')=',ibelong(j)                   
!!            endif                                                     
            IF (r2 < roche2 .AND. ibelong(j) >= 0) THEN 
               ibelong(j) = -i 
               nsplit = nsplit+1 
               lstsplit(nsplit) = j 
            ENDIF 
         ENDDO 
      ENDDO 
                                                                        
      ishift = nchild1*nsplit 
      DO i=npart+1,npart+nghost 
         x(i+ishift) = x(i) 
         y(i+ishift) = y(i) 
         vx(i+ishift) = vx(i) 
         vy(i+ishift) = vy(i) 
         vz(i+ishift) = vz(i) 
         rho(i+ishift) = rho(i) 
         u(i+ishift) = u(i) 
         h(i+ishift) = h(i) 
         pmass(i+ishift) = pmass(i) 
      ENDDO 
                                                                        
      half = 0.5d0 
      third = 1.0d0/3.0d0 
      sixth = half*third 
      sqrt3 = SQRT(3.0d0) 
      sqrt23 = SQRT(2.0d0/3.0d0) 
      cr13 = DBLE(nchild)**(-third) 
      DO ii = 1, nsplit 
!----    Parent's quantities                                            
         ip = lstsplit(ii) 
         xp = x(ip) 
         yp = y(ip) 
         zp = z(ip) 
         vxp = vx(ip) 
         vyp = vy(ip) 
         vzp = vz(ip) 
         rhop = rho(ip) 
         up = u(ip) 
         hp = h(ip) 
         pmassp = pmass(ip) 
         dgravp = dgrav(ip) 
         potenp = poten(ip) 
         iphasep = iphase(ip) 
         iap = iantigr(ip) 
         ibp = ibelong(ip) 
         isp = isteps(ip) 
                                                                        
!----    Child's quantities                                             
         hc = hp*cr13 
         pmassc = pmassp/nchild 
         dist = 1.5d0*hc 
                                                                        
         h(ip) = hc 
         pmass(ip) = pmassc 
                                                                        
         DO i = 2, nchild 
            j = npart+nchild1*(ii-1)+i-1 
                                                                        
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
                                                                        
            vx(j) = vxp 
            vy(j) = vyp 
            vz(j) = vzp 
            rho(j) = rhop 
            u(j) = up 
            h(j) = hc 
            pmass(j) = pmassc 
            dgrav(j) = dgravp 
            poten(j) = potenp 
            iphase(j) = iphasep 
            iantigr(j) = iap 
            ibelong(j) = ibp 
            isteps(j) = isp 
         ENDDO 
      ENDDO 
                                                                        
      numsplit = numsplit+nsplit 
      npart = npart+nchild1*nsplit 
                                                                        
      END SUBROUTINE psplit
