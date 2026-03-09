      SUBROUTINE hdot(npart, ipart, dt, h)
!************************************************************
!                                                           *
!  This subroutine computes the derivative of the smoothing *
!     length.                                               *
!                                                           *
!************************************************************

      use idims

      use force, only: fx, fy, fz, du, dh
      use densi
      use divve
      use logun
      use rbnd
      use diskbd
      use useles
      use outneigh
      use tlist
      use btree
      use phase
      use ptmass
      use integ
      use isnpt
      use glrho
      use neighbor
      use g3monit
      use soft

      implicit none

      INTEGER(I4B) :: npart, ipart, numneigh
      REAL(DP) :: h(idim)
      REAL(DP) :: dt, dhdt, dhi, dhs, dhl, dnsup, dninf, wsex1, &
               wsex2, wiex1, wiex2, hnewg
!
!--Compute derivative of h, try to enforce finite range
!     of neighbors
!
      dhdt = 0.05 / dt

      IF (iphase(ipart).GE.1) THEN
         dh(ipart) = 0.
         GOTO 15
      ELSE
         numneigh = nneigh(ipart)
         IF (isnearpt(ipart).GT.0) &
              numneigh = numneigh + numneighadd(ipart)
      ENDIF

      IF (h(ipart).LT.hmin .AND. numneigh.GT.neimin) THEN
         dh(ipart) = 0.
         GOTO 15
      END IF

      dhi = h(ipart) * divv(ipart) / rho(ipart) / 3.0
      dhs = -dhdt * h(ipart)
      dhl =  dhdt * h(ipart)
      dnsup = MAX(neimax - numneigh, -100)
      dninf = MAX(numneigh - neimin, -100)
      IF (dnsup.LT.nrange) THEN
         wsex1 = EXP(dnsup / 3.5)
         wsex2 = 1. / wsex1
         dhi = (wsex1 * dhi + wsex2 * dhs) / (wsex1 + wsex2)
      ELSEIF (dninf.LT.nrange) THEN
         wiex1 = EXP(dninf / 3.5)
         wiex2 = 1. / wiex1
         dhi = (wiex1 * dhi + wiex2 * dhl) / (wiex1 + wiex2)
      END IF

      dh(ipart) = dhi

 200  hnewg = h(ipart) + dh(ipart)*dt
      IF (hnewg.LE.0) THEN
!      IF (ABS(hnewg-h(ipart)).GE.ABS(h(ipart)/2.0)) THEN
!cc         WRITE(iprint,*) ' h < 0 ', dh(ipart), dhi,
!cc     &        ipart, numneigh, h(ipart), hnewg, divv(ipart)
         IF (dh(ipart).LT.0.)  THEN
            dh(ipart) = dh(ipart)/2.0
!cc            WRITE(iprint,*) ipart, dh(ipart)
            GOTO 200
         ENDIF
      ENDIF

!cc        dh(ipart) = 0.0

 15   CONTINUE

      END SUBROUTINE hdot
