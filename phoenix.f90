      SUBROUTINE phoenix(ipt, i, idtsyn, itime)
!***********************************************************************
!                                                                      *
!  Assigns an accreted particle or particle that goes outside the      *
!     dead boundary a new position, velocity etc to allow accretion    *
!     flow to stablise.                                                *
!                                                                      *
!***********************************************************************

      use mpi_mod
      use idims

      use constants
      use part
      use timei
      use init
      use dum
      use phase
      use carac
      use varet
      use binary
      use ghost
      use ener
      use force, only: f1vx, f1vy, f1vz, f1u, f1h
      use logun
      use typef
      use diskbd
      use rbnd
      use units
      use new
      use cgas
      use sphcom
      use crpart
      use ptmass
      use useles
      use neighbor
      use outneigh
      use misali
      use maslos
      use capt
      use winds
      use unidis
      use cooldata

      implicit none

      INTEGER(I4B) :: ipt, i, idtsyn, itime, nshellx, itry, ind
      REAL(DP) :: ran1
      REAL(DP) :: radinject, radinject2, rnd1, rnd2, rnd3, rnd4, &
           ang1, ang2, s1, s2, s3, rxy, rxy2, r2, fractanhere, &
           accelcent, accelpt, r, tanmag, radmag, vk0, hr0, &
           height00, gg2, rnd, &
           x1, y1, z1, height0, height5, den1, alpha, xtild, ytild, &
           vx1, vy1, vz1, xx, yy, zz, vxx, vyy, vzz
      REAL(DP) :: rshellx, vwindx, vrotx, cosavoid, vangle0, &
           prob1, tmp1, azimuth0, phi1, r1, &
           costheta, sintheta, cosphi, sinphi, hsinva, areafrac, &
           hi1, hi2, cosvang1, cosvang2, r1_xy
      REAL(DP) :: x1temp, y1temp, z1temp, vx1temp, vy1temp, &
           vz1temp, canginj1, canginj2, sanginj1, sanginj2
      REAL(DP) :: sx, sy, sz, svx, svy, svz, &
           prob2, prob3
      CHARACTER(len=7) :: inject

!---- Specify the method for particle injection
      inject = 'random'
!!      inject = 'uniform'

!!      WRITE (iprint,*) 'ibound =',ibound,' in phoenix'
      IF (ibound <= 91) THEN
         radinject = 0.98*deadbound
         radinject2 = radinject**2

         DO
            rnd1 = ran1(1)
            rnd2 = ran1(1)
            rnd3 = ran1(1)

            ang1 = 2*pi*rnd1
            ang2 = pi*rnd2
            s2 = SIN(ang2)
!
!--Make UNIFORM distribution over the surface of a SPHERE.
!
            IF (rnd3 > s2) CYCLE

            x(i) = radinject*COS(ang1)*s2
            y(i) = radinject*SIN(ang1)*s2
            z(i) = radinject*COS(ang2)

            rxy2 = x(i)*x(i) + y(i)*y(i)
            rxy = SQRT(rxy2)
            r2 = rxy2 + z(i)*z(i)
!
!--If ibound=90, then inject all particles with same specific angular momentum
!     ibound=91, then inject particles with uniform angular velocity at
!                    injection radius, radinject
!
            IF (ibound == 90) THEN
               fractanhere = fractan
            ELSEIF (ibound == 91) THEN
               fractanhere = fractan*(rxy2/radinject2)
            ENDIF

            accelcent = rxy*(fractanhere*specang/rxy2)**2
            accelpt = ptmassin/r2*SIN(ATAN2(rxy,ABS(z(i))))
            IF (accelcent <= accelpt) EXIT
         END DO

         r = SQRT(r2)
!
!--Velocity equals tangential + radial contributions.  The tangential
!     velocity is defined as a fraction (fractanhere) of the specific angular
!     momentum (specang) of a particle in circular orbit about the total
!     mass of the binary at the semi-major axis of the system (=1 in
!     code units).
!     The radial velocity is defined as a fraction (fracradial) of the
!     orbital velocity at the semi-major axis.
!
         tanmag = fractanhere*specang/rxy
!
!--Maximum radial velocity gives kinetic energy of particle to be that of a
!     particle which has fallen from infinity, accounting for tangential
!     velocity.  Note that specang=SQRT(Mtot).
!
         radmag = fracradial*SQRT(2.0*specang*specang/r &
                                  - tanmag*tanmag)
         vx(i) = - radmag*x(i)/r - tanmag*y(i)/rxy
         vy(i) = - radmag*y(i)/r + tanmag*x(i)/rxy
         vz(i) = - radmag*z(i)/r
         u(i) = u(i)
         h(i) = deadbound/(DBLE(npart)**(1.0/3.0))
         IF (igrape == 1) hmax(i) = h(i)
         pmass(i) = partm

      ELSEIF (ibound == 92 .OR. ibound == 93) THEN

!-- Radial scale-height of the boundary layer
!   = a half of the thickness of the boundary layer
!!         rscalebnd = 0.5*(rshell-rmind)
!!         rscalebnd2 = rscalebnd*rscalebnd
!-- Note: height used is that for the isothermal disk
         vk0 = SQRT(pmass(ipt)/rmind)*udist/utime
         hr0 = SQRT(gamma*thermal1/(3.0/2.0/uergg))/vk0
!-- Confine the disk to |z|<=2H
         height00 = 2. * hr0 * rmind
!!         IF (height0 > zmax) THEN
!!            zmax = height0
!!         ENDIF
         gg2 = 1.

         rnd2 = ran1(1)

         radinject = rmind + (rshell-rmind)*0.5
!!         radinject = rmind + (rshell-rmind)*rnd1
!!         radinject1 = (rshell-rmind)*rnd1
!!         radinject2 = radinject1*radinject1
         ang1 = 2*pi*rnd2
         x1 = radinject*COS(ang1)
         y1 = radinject*SIN(ang1)
!-- Change the following line if the disk is not isothermal!
         height0 = height00 * (radinject/rmind)**1.5
         height5 = 0.5*height0

         DO
            rnd3 = ran1(1)
            z1 = 2.*(height0*rnd3 - height5)

!!            den1 = EXP(-0.5*radinject2/rscalebnd2 &
!!                       -0.5*z1*z1/(height5*height5))
            den1 = EXP(-0.5*(z1/(hr0*rmind*(radinject/rmind)**1.5))**2)
            rnd = ran1(1)
            IF (rnd <= den1) EXIT
         ENDDO

         alpha = SQRT(gg2*pmass(ipt)/radinject)
         xtild = x1/radinject
         ytild = y1/radinject
         vx1 = -alpha * ytild
         vy1 =  alpha * xtild
         vz1 =  0.0

!--      Rotate SPH particles about x- and y-axes
         IF (rangle(1) == 0.0 .AND. rangle(2) == 0.0) THEN
            x(i) = x1
            y(i) = y1
            z(i) = z1
            vx(i) = vx1
            vy(i) = vy1
            vz(i) = vz1
         ELSE
            xx = x1
            yy = y1
            zz = z1
            x(i) = arot(1,1)*xx+arot(1,2)*yy+arot(1,3)*zz
            y(i) = arot(2,1)*xx+arot(2,2)*yy+arot(2,3)*zz
            z(i) = arot(3,1)*xx+arot(3,2)*yy+arot(3,3)*zz
            vxx = vx1
            vyy = vy1
            vzz = vz1
            vx(i) = arot(1,1)*vxx+arot(1,2)*vyy &
                    +arot(1,3)*vzz
            vy(i) = arot(2,1)*vxx+arot(2,2)*vyy &
                    +arot(2,3)*vzz
            vz(i) = arot(3,1)*vxx+arot(3,2)*vyy &
                    +arot(3,3)*vzz
         ENDIF

         IF (ipt /= 0) THEN
            x(i) = x(i) + x(ipt)
            y(i) = y(i) + y(ipt)
            z(i) = z(i) + z(ipt)
            vx(i) = vx(i) + vx(ipt)
            vy(i) = vy(i) + vy(ipt)
            vz(i) = vz(i) + vz(ipt)
         ENDIF

!!         u(i) = u(i)
         u(i) = thermal1
!-- Before 17 April 2001
!!         h(i) = ((rshell*rshell-rmind*rmind)*pi*2.0*height0 &
!!                 /DBLE(nshell))**(1.0/3.0)
!-- Since 17 April 2001
!!         h(i) = SQRT((rshell*rshell-rmind*rmind)*pi &
!!                 /DBLE(inshell))
!-- Since 14 September 2003 (in case inshell=0)
         h(i) = ((rshell*rshell-rmind*rmind)*pi*2.0*height0 &
                 /DBLE(nshell))**(1.0/3.0)

         IF (igrape == 1) hmax(i) = h(i)
         pmass(i) = partm
!!         write (iprint,*) 'gamma =',gamma,', thermal =',thermal
!!         write (iprint,*) 'vk0 =',vk0,', hr0 =',hr0, &
!!                          ', pi =',pi,', height0 =',height0
!!         write (iprint,*) 'radinject =',radinject,' in phoenix: ', &
!!                     'h(',i,') =',h(i)
!!         WRITE (iprint,*) i,': x1=',x1,', y1=',y1,', z1=',z1, &
!!                          ' (z1/H=',z1/height5,')'
      ELSEIF (ibound >= 94 .AND. ibound <= 97) THEN
!--      Escape velocity and critical velocity of each star
!!         vcritx = SQRT(pmass(ipt)/rptmas(ipt))
!!         vescx  = SQRT(2.0)*vcritx

!--      Distribute particles around each star
         SELECT CASE (ipt)
         CASE (1)
            nshellx = nshell1(ipt)
            rshellx = rshell1
            vwindx = vwind1
            vrotx = vrot1
            pmass(i) = partmass(ipt)
            u(i) = thermal1
            IF (inject == 'uniform') THEN
               IF(MOD(nUD1,idiminj) == 0) THEN
                  nUD1 = 1
               ELSE
                  nUD1 = nUD1+1
               ENDIF
               nUD = nUD1
            ENDIF

         CASE (2)
            IF (nshell1(2) /= 0) THEN
               nshellx = nshell1(ipt)
               rshellx = rshell2
               vwindx = vwind2
               vrotx = vrot2
!----          Note that emdotratio=emdot2/emdot1, not emdot1/emdot2
!              Also note that the mass of wind particles from star2
!              is given by partm*(emdot2/nshell1(2))/(emdot1/nshell1(1))
               pmass(i) = partmass(ipt)
!!               IF (varsta == 'intener') THEN
                  u(i) = thermal2
!!               ELSE
!!                  u(i) = thermal &
!!                     /(emdotratio*REAL(nshell1(1))/REAL(nshell1(2)) &
!!                     /pmass(i))**(gamma-1.0)
!!               ENDIF
               IF (inject == 'uniform') THEN
                  IF(MOD(nUD2,idiminj) == 0) THEN
                     nUD2 = 1
                     anginj1 = 2*pi*ran1(1)
                     anginj2 = 2*pi*ran1(1)
                  ELSE
                     nUD2 = nUD2+1
                  ENDIF
                  nUD = nUD2
               ENDIF
            ELSE
               WRITE (iprint,"('+++ skipping the call for ipt=',i2, &
                 & ' because nshell1(',i2,') is zero. +++')") ipt, ipt
               RETURN
            ENDIF

         CASE default
            pmass(i) = partmass(ipt)
            u(i) = thermal3

         END SELECT

!----    Set the minimum theta if there is a disk, below which
!        no wind particles should be injected (22 March 2010,
!        A. Okazaki).
         !!IF ((ipt == isphcom .OR. ipt == nptmass+1) .AND. &
         !1   (ibound == 95 .OR. ibound == 97)) THEN
         !IF ((ipt == isphcom .OR. ipt == 3) .AND. &
         !   (ibound == 95 .OR. ibound == 97)) THEN
         IF ((ipt == isphcom .OR. ipt == 3)) THEN
!----       Note: height used is that for the isothermal disk
            vk0 = SQRT(pmass(isphcom)/rptmas(isphcom)) &
                       *udist/utime
            IF (ibound == 95 .OR. ibound == 97) THEN
               hr0 = SQRT(gamma*thermal3/(3.0/2.0/uergg))/vk0
!----          Avoid the innermost disk part confined to |z|<=2H
               cosavoid = 2. * hr0
            ELSE
               cosavoid = 0.0
            END IF
         ELSE
            vk0 = SQRT(pmass(3-isphcom)/rptmas(3-isphcom)) &
                       *udist/utime
            cosavoid = 0.0
         ENDIF

         IF (ipt <= nptmass) THEN
            IF (igeom == 8 .OR. igeom == 9) THEN
               vangle0 = ATAN2(z(2)-z(1), &
                         SQRT((x(2)-x(1))**2+(y(2)-y(1))**2))
               cosvang1 = COS(0.5*pi-vangle0-vangle1(ipt))
               cosvang2 = COS(0.5*pi-vangle0-vangle2(ipt))
            ENDIF

            IF (inject == 'random') THEN
               DO
                  IF (igeom == 8 .OR. igeom == 9) THEN
                     costheta = ((cosvang1+cosvang2)*0.5 &
                                +(cosvang2-cosvang1)*(ran1(1)-0.5))
                  ELSE
                     costheta = 2.0*(ran1(1)-0.5)
                  ENDIF
                  IF (igeom == 9) THEN
                     den1 = 1.0_DP
                  ELSE
!----                The following random number genaration is
!                    for setting the particle injection rate
!                    proportional to 1-(rot param)*sin(theta)**2.
!                    This is to take into account the effect of
!                    gravity darkening due to rapid stellar rotation.
!                    (26 Oct. 2023)
                     ! The following line had a serious bug. 
                     ! SQRT(1.0_DP - vrotx*(1.0_DP-costheta**2)) was used
                     ! instead of 1.0_DP - vrotx*vrotx*(1.0_DP-costheta**2).
                     ! This was corrected on 29 May 2026.
                     !!den1 = SQRT(1.0_DP - vrotx*(1.0_DP-costheta**2))
                     den1 = 1.0_DP - vrotx*vrotx*(1.0_DP-costheta**2)
                  END IF
                  rnd = ran1(1)
                  IF (rnd <= den1 .AND. &
                        ABS(costheta) >= cosavoid) EXIT
               END DO

               IF (igeom == 8 .OR. igeom == 9) THEN
                  azimuth0 = ATAN2(y(2)-y(1),x(2)-x(1))
                  phi1 = azimuth0+azimuth1(ipt) &
                         +(azimuth2(ipt)-azimuth1(ipt))*ran1(1)
               ELSE
                  phi1 = 2.0*pi*ran1(1)
               ENDIF

               r1 = 0.5*(rptmas(ipt)+rshellx)
               sintheta = SQRT(1.0-costheta*costheta)
               cosphi = COS(phi1)
               sinphi = SIN(phi1)
               x1 = r1*sintheta*cosphi
               y1 = r1*sintheta*sinphi
               z1 = r1*costheta
!!               WRITE (iprint,*) i,': x1=',x1,', y1=',y1,', z1=',z1
               !-- Each velocity component is first set in cm/s.
               !   Here, vwind1 and vwind2 are in cm/s and
               !   vrotx is rotation parameter in the range 0-1.
               r1_xy = r1*sintheta
               IF (r1_xy < tiny) THEN
                  vx1 = 0.0_DP
                  vy1 = 0.0_DP
                  vz1 = vwindx
               ELSE
                  IF (igeom == 9) THEN
                     !-- Roche lobe overflow
                     vx1 = vwindx*x1/r1 - vrotx*vk0*sintheta*y1/r1_xy
                     vy1 = vwindx*y1/r1 + vrotx*vk0*sintheta*x1/r1_xy
                     vz1 = vwindx*z1/r1
                  ELSE
                     vx1 = vwindx*SQRT(1.0_DP-vrotx*vrotx*sintheta*sintheta)*x1/r1 &
                        - vrotx*vk0*sintheta*y1/r1_xy
                     vy1 = vwindx*SQRT(1.0_DP-vrotx*vrotx*sintheta*sintheta)*y1/r1 &
                        + vrotx*vk0*sintheta*x1/r1_xy
                     vz1 = vwindx*SQRT(1.0_DP-vrotx*vrotx*sintheta*sintheta)*z1/r1
                  END IF
               END iF
            ELSE IF (inject == 'uniform') THEN
               r1 = 0.5*(rptmas(ipt)+rshellx)
               ind = lUD(nUD)
!X               WRITE(*,*) i-nptmass,lUD(i-nptmass),ind
!X               WRITE(*,*) xUD(lUD(i-nptmass))
!X               WRITE(*,*) yUD(lUD(i-nptmass))
!X               WRITE(*,*) zUD(lUD(i-nptmass))
               x1 = r1*xUD(ind)
               y1 = r1*yUD(ind)
               z1 = r1*zUD(ind)
!X               WRITE(*,*) x1
!X               WRITE(*,*) y1
!X               WRIssssssssssssssssssssssssssssssssssssssssTE(*,*) z1
               !-- Each velocity component is first set in cm/s.
               !   Here, vwindx is in cm/s and
               !   vrotx is rotation parameter in the range 0-1.
               !   Note that currently, gravity darkening effect
               !   is available only for option 'random'.
               !   (22 Oct. 2023)
               r1_xy = SQRT(x1+x1+y1*y1)
               sintheta = r1_xy/r1
               IF (r1_xy < tiny) THEN
                  vx1 = 0.0_DP
                  vy1 = 0.0_DP
                  vz1 = vwindx
               ELSE
                  vx1 = vwindx*x1/r1 - vrotx*vk0*sintheta*y1/r1_xy
                  vy1 = vwindx*y1/r1 + vrotx*vk0*sintheta*x1/r1_xy
                  vz1 = vwindx*z1/r1
               END IF

!UNI---        Rotate anginj1 about z-axis and then anginj2 about x-axis
!               canginj1=COS(anginj1)
!               sanginj1=SIN(anginj1)
!               canginj2=COS(anginj2)
!               sanginj2=SIN(anginj2)
!               x1temp =  x1*canginj1 + y1*sanginj1
!               y1temp = -x1*sanginj1 + y1*canginj1
!               x1 = x1temp
!               y1 = y1temp
!               vx1temp =  vx1*canginj1 + vy1*sanginj1
!               vy1temp = -vx1*sanginj1 + vy1*canginj1
!               vx1 = vx1temp
!               vy1 = vy1temp
!               y1temp =  y1*canginj2 + z1*sanginj2
!               z1temp = -y1*sanginj2 + z1*canginj2
!               y1 = y1temp
!               z1 = z1temp
!               vy1temp =  vy1*canginj2 + vz1*sanginj2
!               vz1temp = -vy1*sanginj2 + vz1*canginj2
!               vy1 = vy1temp
!               vz1 = vz1temp
!!               WRITE (iprint,*) i,': vx1=',vx1,', vy1=',vy1, &
!!                                ', vz1=',vz1
            ELSE
               WRITE (*,*) '### inject=',inject,' is undefined! ###'
               CALL quit
            ENDIF

!----       vwind1, vwind2, vrot1, and vrot2 are given
!           in cm/s. (A. Okazaki, 09/01/2007)
            vx1 = vx1 * utime/udist
            vy1 = vy1 * utime/udist
            vz1 = vz1 * utime/udist
!!            WRITE (iprint,*) i,': vx1=',vx1,', vy1=',vy1,', vz1=',vz1

         ELSE
            gg2 = 1.
!-- Confine the disk to |z|<=2H
            height00 = 2.0 * hr0 * rptmas(isphcom)

            rnd2 = ran1(1)

            IF (isphcom == 1) THEN
               radinject = (rptmas(1)+rshell1)*0.5
               height0 = height00 * (radinject/rptmas(1))**1.5
            ELSE
               radinject = (rptmas(2)+rshell2)*0.5
               height0 = height00 * (radinject/rptmas(2))**1.5
            ENDIF
            ang1 = 2*pi*rnd2
            x1 = radinject*COS(ang1)
            y1 = radinject*SIN(ang1)
!----       Change the following line if the disk is not isothermal!
            height5 = 0.5*height0

            DO
               rnd3 = ran1(1)
               z1 = 2.*(height0*rnd3 - height5)

               den1 = EXP(-0.5*(z1 &
                          /(hr0*rptmas(isphcom) &
                           *(radinject/rptmas(isphcom))**1.5))**2)
               rnd = ran1(1)
               IF (rnd <= den1) EXIT
            ENDDO

            alpha = SQRT(gg2*pmass(isphcom)/radinject)
            xtild = x1/radinject
            ytild = y1/radinject
            vx1 = -alpha * ytild
            vy1 =  alpha * xtild
            vz1 =  0.0
         ENDIF

!--      Rotate SPH particles about x- and y-axes
         IF (rangle(1) == 0.0 .AND. rangle(2) == 0.0) THEN
            x(i) = x1
            y(i) = y1
            z(i) = z1
            vx(i) = vx1
            vy(i) = vy1
            vz(i) = vz1
         ELSE
            xx = x1
            yy = y1
            zz = z1
            x(i) = arot(1,1)*xx+arot(1,2)*yy+arot(1,3)*zz
            y(i) = arot(2,1)*xx+arot(2,2)*yy+arot(2,3)*zz
            z(i) = arot(3,1)*xx+arot(3,2)*yy+arot(3,3)*zz
            vxx = vx1
            vyy = vy1
            vzz = vz1
            vx(i) = arot(1,1)*vxx+arot(1,2)*vyy &
                    +arot(1,3)*vzz
            vy(i) = arot(2,1)*vxx+arot(2,2)*vyy &
                    +arot(2,3)*vzz
            vz(i) = arot(3,1)*vxx+arot(3,2)*vyy &
                    +arot(3,3)*vzz
         ENDIF

         IF (ipt /= 0) THEN
            IF (ipt <= nptmass) THEN
               x(i) = x(i) + x(ipt)
               y(i) = y(i) + y(ipt)
               z(i) = z(i) + z(ipt)
               vx(i) = vx(i) + vx(ipt)
               vy(i) = vy(i) + vy(ipt)
               vz(i) = vz(i) + vz(ipt)
            !ELSEIF (ipt == nptmass+1) THEN
            ELSEIF (ipt == 3) THEN
               x(i) = x(i) + x(isphcom)
               y(i) = y(i) + y(isphcom)
               z(i) = z(i) + z(isphcom)
               vx(i) = vx(i) + vx(isphcom)
               vy(i) = vy(i) + vy(isphcom)
               vz(i) = vz(i) + vz(isphcom)
            ENDIF
         ENDIF

         IF (igeom == 8 .OR. igeom == 9) THEN
!----       areafac: corrected (20 April 2023)
            hsinva = 0.5d0*(SIN(vangle2(ipt))-SIN(vangle1(ipt)))
            areafrac = hsinva &
                      *(azimuth2(ipt)-azimuth1(ipt))/(2.0d0*pi)
         ELSE
            areafrac = 1.0d0
         ENDIF

!----    ibelong is the tag that specifies which point-mass
!        group the individual particles belong to.
         IF (ipt <= nptmass) THEN
            hi1 = ((rshellx**3-rptmas(ipt)**3)*4.0*pi/3.0*areafrac &
                /REAL(nshellx))**(1.0/3.0)
            hi2 = 0.5*deadbound/(DBLE(npart)**(1.0/3.0))
            h(i) = MIN(hi1,hi2)
!!            h(i) = ((rshellx**3-rptmas(ipt)**3)*4.0*pi/3.0*areafrac &
!!               /REAL(nshellx))**(1.0/3.0)
            ibelong(i) = ipt
            IF (igeom == 9) THEN
               iantigr(i) = 0
            ELSE
               iantigr(i) = 1
            END IF
         ELSE
!----       The following part is executed only if ibound=95 or 97
            IF (isphcom == 1) THEN
               hi1 = (4.0*pi*rptmas(isphcom)**2 &
                     *(rshell1-rptmas(isphcom)) &
                     /nshell1(ipt))**(1.0/3.0)
            ELSE
               hi1 = (4.0*pi*rptmas(isphcom)**2 &
                     *(rshell2-rptmas(isphcom)) &
                     /nshell1(ipt))**(1.0/3.0)
            ENDIF
            hi2 = 0.5*deadbound/(DBLE(npart)**(1.0/3.0))
            h(i) = MIN(hi1,hi2)
            !!            h(i) = hi1
            

            ibelong(i) = isphcom
            iantigr(i) = 0
         ENDIF

         IF (igrape == 1) hmax(i) = h(i)

      ELSEIF (ibound == 99) THEN
!--      Inject particles whose distribution is statistcally
!        equivalent with that of the captured particles
!        in the corresponding BeX simulation within a sigma
         sx = SQRT(sxinfl(iphsi))
         sy = SQRT(syinfl(iphsi))
         sz = SQRT(szinfl(iphsi))
         svx = SQRT(svxinfl(iphsi))
         svy = SQRT(svyinfl(iphsi))
         svz = SQRT(svzinfl(iphsi))

         DO
            rnd1 = ran1(1)
            s1 = 2.*(rnd1 - 0.5)
            prob1 = EXP(-0.5*s1*s1)/SQRT(2.0*pi)

            rnd2 = ran1(1)
            s2 = 2.*(rnd2 - 0.5)
            prob2 = EXP(-0.5*s2*s2)/SQRT(2.0*pi)

            rnd3 = ran1(1)
            s3 = 2.*(rnd3 - 0.5)
            prob3 = EXP(-0.5*s3*s3)/SQRT(2.0*pi)

            rnd4 = ran1(1)
            IF (rnd4 <= prob1*prob2*prob3) EXIT
         ENDDO

         x(i) = xinfl(iphsi) + sx*s1 + x(ipt)
         y(i) = yinfl(iphsi) + sy*s2 + y(ipt)
         z(i) = zinfl(iphsi) + sz*s3 + z(ipt)

         DO
            rnd1 = ran1(1)
            s1 = 2.*(rnd1 - 0.5)
            prob1 = EXP(-0.5*s1*s1)/SQRT(2.0*pi)

            rnd2 = ran1(1)
            s2 = 2.*(rnd2 - 0.5)
            prob2 = EXP(-0.5*s2*s2)/SQRT(2.0*pi)

            rnd3 = ran1(1)
            s3 = 2.*(rnd3 - 0.5)
            prob3 = EXP(-0.5*s3*s3)/SQRT(2.0*pi)

            rnd4 = ran1(1)
            IF (rnd4 <= prob1*prob2*prob3) EXIT
         ENDDO

         vx(i) = vxinfl(iphsi) + svx*s1 + vx(ipt)
         vy(i) = vyinfl(iphsi) + svy*s2 + vy(ipt)
         vz(i) = vzinfl(iphsi) + svz*s3 + vz(ipt)

         u(i) = thermal1

!----    The following setting results in an infinite loop,
!        making h negative. Why? (10 May 2008, A. Okazaki)
!!         h(i) = hinfl(iphsi)/emdot0**(1.0d0/3.0d0)
         h(i) = SQRT(xinfl(iphsi)*xinfl(iphsi) &
                     + yinfl(iphsi)*yinfl(iphsi) &
                     + zinfl(iphsi)*zinfl(iphsi)) &
                /SQRT(REAL(npart))
!!         IF (hmaximum > 0.0) THEN
!!            IF (h(i) > hmaximum) THEN
!!               h(i) = hmaximum
!!               ioutmax = ioutmax + 1
!!            ENDIF
!!         ENDIF
!         IF (igrape == 1) hmax(i) = h(i)
         pmass(i) = partm

!!         IF (myrank.eq.0) THEN
!!            WRITE (iprint,*) 'i=',i,' (npart=',npart,'): h=',h(i)
!!            WRITE (iprint,*) 'iphsi=',iphsi,': x(i)=',x(i), &
!!                         ', y(i)=',y(i), &
!!                         ', z(i)=',z(i)
!!         ENDIF
      ELSE
         IF (myrank.eq.0) WRITE (iprint,*) 'ERROR in phoenix'
         IF (myrank.eq.0) WRITE (iprint,*) 'ibound =',ibound
         CALL quit
      ENDIF

      dumx(i) = x(i)
      dumy(i) = y(i)
      dumz(i) = z(i)
      dumvx(i) = vx(i)
      dumvy(i) = vy(i)
      dumvz(i) = vz(i)
      dumu(i) = u(i)
      dumh(i) = h(i)

      f1vx(i) = 0.0
      f1vy(i) = 0.0
      f1vz(i) = 0.0
      f1u(i) = 0.0
      f1h(i) = 0.0

      dgrav(i) = 0.0
      hasghost(i) = 0

      itry = istepmax/2
!!!      itry = istepmin
!!!      itry = imaxstep
      DO
         IF (MOD(idtsyn, itry) /= 0) THEN
            itry = itry/2
         ELSE
            isteps(i) = itry
            EXIT
         ENDIF
      ENDDO

      it0(i) = itime
      it1(i) = it0(i) + isteps(i)/2
      iphase(i) = 0
      iinold(i)  = 0

      nreassign = nreassign + 1

!---- Initialize arrays for radiative cooling info
      IF (encal == 'c') THEN
         DO i=1,npart
            tempini(i) = 0.0
            tempfin(i) = tempini(i)
            tcool(i) = 0.0
         END DO
      END IF

      END SUBROUTINE phoenix
