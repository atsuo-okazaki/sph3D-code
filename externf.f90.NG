      SUBROUTINE externf(ipart, x, y, z, iexf)
!************************************************************
!                                                           *
!  This subroutine computes the effect of an external       *
!     force.                                                *
!                                                           *
!        (1) Vertical gravitational field                   *
!        (2) -------                                        *
!        (3) Accretion disk                                 *
!        (4) Rotating cylinder                              *
!        (5) Central point mass                             *
!        (6) Distant point mass                             *
!        (7) Antigravity for stellar winds                  *
!        (8) Radiative inhibition of stellar winds          *
!                                                           *
!************************************************************

      use idims

      use kerne
      use ener
      use force, only: fx, fy, fz, du, dh
      use units
      use xforce
      use xtorq
      use carac
      use ptmass
      use gravi
      use sphcom
      use winds
      use logun

      implicit none

      INTEGER(I4B) :: ipart, iexf, j
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: aminrate, aminx, aminy, aminz, d, d2, d3, &
               di, dunix, duniy, dx, dy, dz, gravxi, gravyi, gravzi, &
               h1, h2, omeg,omega, pmassj, poteni, rij, rij1, rij2, &
               uang, grav, xi, yi, zi, runix, runiy, runiz, xantigr, &
               xc, yc, zc, zdist, xmasj, xtrqcoef
!
!--Unit angular momentum
!
      uang = udist**2/utime
!
!--Gravitational field
!
      IF (iexf.EQ.1) THEN
         grav = 1.0E4*utime**2/udist

         fz(ipart) = fz(ipart) - grav
!
!--Accretion disk
!
      ELSEIF (iexf.EQ.3) THEN
         IF (ipart.LE.nptmass) RETURN
!-- Original part
!c         t1 = 1.
!c         q = 2.
!c         h0 = 0.8
!c         r0 = h0**2
!c         r02 = r0**2
!c         omega = h0/r02
!c         d2 = x(ipart)**2 + y(ipart)**2
!c         d = SQRT(d2)
!c         d3 = d2*d
!c         omeg = omega*(r0/d)**q
!c         h1 = omeg*d2
!c         h2 = h1**2
!c         r2 = d2 + z(ipart)**2
!c         r = SQRT(r2)
!c
!c         runix = x(ipart)/r
!c         runiy = y(ipart)/r
!c         runiz = z(ipart)/r
!c         dunix = x(ipart)/d
!c         duniy = y(ipart)/d
!c
!c         fx(ipart) = fx(ipart) - 0.9999*runix/r2 + t1*h2*dunix/d3
!c         fy(ipart) = fy(ipart) - 0.9999*runiy/r2 + t1*h2*duniy/d3
!c         fz(ipart) = fz(ipart) - 0.9999*runiz/r2

!-- 27 November 2000
!-- Torque on an innermost part of the disk is introduced to
!   mimic the torque exerted by the Be star on the disk and
!   prevent particles from accretion.
!      The functional form is the same as that of Case B
!   in Pringle (1991, MNRAS 248, 754-759), except that parameters
!   eps (xeps here) and beta (xbeta here) are chosen such that:
!       - the torque is strong enough to prevent accretion onto
!         point mass 1,
!       - the torque is not so strong that it affects the outer
!         disk structure.
!   In reality, there is no parameter range that completely
!   satisfies these conditions. Therefore, you have to compromize,
!   i.e., you have to allow accretion of a small number of
!   particles and small (but significant) effect on the disk
!   structure. The best combination among those I've tested is
!   xeps=0.03 and xbeta=4 (for 4U0115+63).
!      Important Note: Avoid using too large xbeta (xbeta>=5),
!         which strongly disturbes the disk. Also avoid using the
!         torque in a form of exponential function, which also
!         disturbes the disk.

         xc = x(isphcom)
         yc = y(isphcom)
         zc = z(isphcom)
         xi = x(ipart) - xc
         yi = y(ipart) - yc
         zi = z(ipart) - zc

         di = SQRT(xi*xi+yi*yi)
         dunix = xi/di
         duniy = yi/di

!-- Power-law type torque
         xtrqcoef = xeps / (di/rptmas(isphcom))**xbeta
         aminrate = pmass(isphcom)/(di*di) * xtrqcoef

         aminx = -aminrate*duniy
         aminy = aminrate*dunix
         aminz = 0.0

         fx(ipart) = fx(ipart) + aminx
         fy(ipart) = fy(ipart) + aminy
         fz(ipart) = fz(ipart) + aminz
!
!--Rotating cylinders
!
      ELSEIF (iexf.EQ.4) THEN
         omega = 0.6
         d2 = x(ipart)**2 + y(ipart)**2
         d = SQRT(d2)
         d3 = d2*d
         omeg = omega
         IF (d.GT.1.) omeg = 0.
         h1 = d2*omeg
         h2 = h1**2
         dunix = x(ipart)/d
         duniy = y(ipart)/d

         fx(ipart) = fx(ipart) + h2*dunix/d3
         fy(ipart) = fy(ipart) + h2*duniy/d3
!
!--Central point mass
!
      ELSEIF (iexf.EQ.5) THEN
         xi = x(ipart)
         yi = y(ipart)
         zi = z(ipart)
         d2 = (xi*xi + yi*yi + zi*zi + tiny)
         d = SQRT(d2)
         runix = xi/d
         runiy = yi/d
         runiz = zi/d

         fx(ipart) = fx(ipart) - xmass*runix/d2
         fy(ipart) = fy(ipart) - xmass*runiy/d2
         fz(ipart) = fz(ipart) - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d
!
!--Distant point mass
!
      ELSEIF (iexf.EQ.6) THEN
         zdist = 100.
         xi = x(ipart)
         yi = y(ipart)
         zi = zdist - z(ipart)
         d2 = (xi*xi + yi*yi + zi*zi)
         d = SQRT(d2)
         runix = xi/d
         runiy = yi/d
         runiz = zi/d

         fx(ipart) = fx(ipart) - xmass*runix/d2
         fy(ipart) = fy(ipart) - xmass*runiy/d2
         fz(ipart) = fz(ipart) - xmass*runiz/d2
         poten(ipart) = poten(ipart) - xmass/d

!---- Antigravity
      ELSEIF (iexf.EQ.7 .OR. iexf.EQ.8) THEN
         IF (iantigr(ipart).NE.0) THEN
!--Gravity and potential energy
            xi = x(ipart)
            yi = y(ipart)
            zi = z(ipart)
            gravxi = 0.
            gravyi = 0.
            gravzi = 0.
            poteni = 0.
            DO j=1,nptmass
               dx = xi - x(j)
               dy = yi - y(j)
               dz = zi - z(j)
               rij2 = dx*dx + dy*dy + dz*dz + tiny
               rij = SQRT(rij2)
               rij1 = 1./rij
               pmassj = pmass(j)

!--Unit vectors
               runix = dx*rij1
               runiy = dy*rij1
               runiz = dz*rij1
!c               phi = -rij1

               xmasj = pmassj/rij2
!c               gravxi = gravxi - xmasj*runix
!c               gravyi = gravyi - xmasj*runiy
!c               gravzi = gravzi - xmasj*runiz
!c               poteni = poteni + phi*pmassj

!---- Antigravity
               IF (iexf.EQ.7) THEN
                  IF (ibelong(ipart).EQ.j) THEN
                     xantigr = xantgrav(j)*akappa(ipart)
                  ELSE
                     xantigr = 1.0d0
                  ENDIF
               ELSE
                  xantigr = xantgrav(j)*akappa(ipart)
               ENDIF
               fx(ipart) = fx(ipart) + xantigr*xmasj*runix
               fy(ipart) = fy(ipart) + xantigr*xmasj*runiy
               fz(ipart) = fz(ipart) + xantigr*xmasj*runiz
!c               poten(ipart) = poten(ipart) - xantigr*phi*pmassj
            ENDDO
         ENDIF
      ENDIF

      END SUBROUTINE externf
