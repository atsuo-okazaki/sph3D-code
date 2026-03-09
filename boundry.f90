      SUBROUTINE boundry(npart,list,nlst0,x,y,z,vx,vy,vz,u,h, &
                         idonebound)
!************************************************************
!                                                           *
!   This subroutine checks to see if particles will surpass *
!   boundary and resets them to the boundary if they do     *
!                                                           *
!************************************************************

      use idims

      use constants
      use typef
      use logun
      use rbnd
      use diskbd
      use debug
      use phase
      use ptmass
      use carac
      use sphcom

      implicit none

      INTEGER(I4B) :: list(idim), npart, nlst0, idonebound, i, j, ichan, &
             iinner, iouter
      REAL(DP) :: ran1, x(idim), y(idim), z(idim)
      REAL(DP) :: vx(idim), vy(idim), vz(idim), u(idim), h(idim)
      REAL(DP) :: alpha, delta, r, r2, rcyl2, rcyld, rmax2, rmaxd, &
             rmind2, rmindd, vi, zmind, zmaxd, vr, xi, yi, zi, &
             xc, yc, zc, vxc, vyc, vzc, vxi, vyi, vzi, vi0, &
             xtild, ytild

!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine boundary')

      iinner = 0
      iouter = 0
      idonebound = 0
      IF ( ibound.EQ.1 ) THEN
         DO 250 j = 1, nlst0
            i = list(j)
            IF (iphase(i).NE.0) GOTO 250
            delta = (1.3-ran1(1))*0.5*h(i)
            ichan = 0
            IF (x(i).LE.xmin) THEN
               ichan = ichan + 1
               vx(i) = 0.
               x(i) = xmin + delta
            ENDIF
            IF (x(i).GE.xmax) THEN
               ichan = ichan + 1
               vx(i) = 0.
               x(i) = xmax - delta
            ENDIF
            IF (y(i).LE.ymin) THEN
               ichan = ichan + 1
               vy(i) = 0.
               y(i) = ymin + delta
            ENDIF
            IF (y(i).GE.ymax) THEN
               ichan = ichan + 1
               vy(i) = 0.
               y(i) = ymax - delta
            ENDIF
            IF (z(i).LE.zmin) THEN
               ichan = ichan + 1
               vz(i) = 0.
               z(i) = zmin + delta
            ENDIF
            IF (z(i).GE.zmax) THEN
               ichan = ichan + 1
               vz(i) = 0.
               z(i) = zmax - delta
            ENDIF
            IF (ichan.NE.0) THEN
               iouter = iouter + 1
               idonebound = 1
            ENDIF
 250     CONTINUE

      ELSE IF ( ibound.EQ.2 ) THEN
         rmind2 = rmind * rmind
         IF (isphcom.EQ.0) THEN
            rcyl2 = rcyl * rcyl
            DO 260 j = 1, nlst0
               i = list(j)
               IF (iphase(i).NE.0) GOTO 260
               IF (i.LE.nptmass) GOTO 260
               delta = (1.3-ran1(1))*0.5*h(i)
               rcyld = rcyl - delta
               rmindd = rmind + delta
               zmind = zmin + delta
               zmaxd = zmax - delta
               ichan = 0
               r2 = x(i)*x(i) + y(i)*y(i)
               IF ( r2.GT.rcyl2) THEN
                  ichan = ichan + 1
                  r = SQRT(r2)
                  vr = (vx(i)*x(i) + vy(i)*y(i))/r
                  xi = x(i)*rcyld/r
                  yi = y(i)*rcyld/r
                  zi = z(i)
!
!--Remove velocities perpendicular to boundary
!
                  IF (vr.GT.0) THEN
                     vx(i) = vx(i) - vr*x(i)/r
                     vy(i) = vy(i) - vr*y(i)/r
                  ENDIF
!
!--Adjust to conserve angular momentum
!
                  IF (ifcor.EQ.1) THEN
                     vx(i) = vx(i) * y(i) / yi
                     vy(i) = vy(i) * x(i) / xi
                  ELSEIF (ifcor.EQ.2) THEN
                     vy(i) = vy(i) * z(i) / zi
                     vz(i) = vz(i) * y(i) / yi
                  ENDIF
                  y(i) = yi
                  x(i) = xi
               ENDIF

               IF (z(i).LE.zmin) THEN
                  ichan = ichan + 1
                  vz(i) = 0.
                  z(i) = zmind
               ENDIF
               IF (z(i).GE.zmax) THEN
                  ichan = ichan + 1
                  vz(i) = 0.
                  z(i) = zmaxd
!
!--Remove velocities perpendicular to boundary
!
                  IF (vr.GT.0) THEN
                     vx(i) = vx(i) - vr*x(i)/r
                     vy(i) = vy(i) - vr*y(i)/r
                  ENDIF
!
!--Adjust to conserve angular momentum
!
                  IF (ifcor.EQ.1) THEN
                     vx(i) = vx(i) * y(i) / yi
                     vy(i) = vy(i) * x(i) / xi
                  ELSEIF (ifcor.EQ.2) THEN
                     vy(i) = vy(i) * z(i) / zi
                     vz(i) = vz(i) * y(i) / yi
                  ENDIF
                  y(i) = yi
                  x(i) = xi
               ENDIF

               IF (z(i).LE.zmin) THEN
                  ichan = ichan + 1
                  vz(i) = 0.
                  z(i) = zmind
               ENDIF
               IF (z(i).GE.zmax) THEN
                  ichan = ichan + 1
                  vz(i) = 0.
                  z(i) = zmaxd
               ENDIF
               IF (ichan.NE.0) THEN
                  iouter = iouter + 1
                  idonebound = 1
               ENDIF
 260        CONTINUE
         ELSE
!-- Don't use the outer boundary for disks in binary systems
            rcyl2 = 9.99e19
            xc = x(isphcom)
            yc = y(isphcom)
            zc = z(isphcom)
            vxc = vx(isphcom)
            vyc = vy(isphcom)
            vzc = vz(isphcom)
            DO 270 j = 1, nlst0
               i = list(j)
               IF (iphase(i).NE.0) GOTO 270
               IF (i.LE.nptmass) GOTO 270
               delta = (1.3-ran1(1))*0.5*h(i)
               rcyld = SQRT(rcyl2) - delta
               rmindd = rmind + delta
               zmind = zmin + delta
               zmaxd = zmax - delta
               ichan = 0
               r2 = (x(i)-xc)*(x(i)-xc) &
                    + (y(i)-yc)*(y(i)-yc)

               IF ( r2.LT.rmind2) THEN
                  iinner = iinner + 1
                  r = sqrt(r2)
                  xi = (x(i)-xc)*rmindd/r
                  yi = (y(i)-yc)*rmindd/r
                  zi = z(i)-zc
!
!--Remove relative velocities perpendicular to boundary
!
!!cc                  vr = ((vx(i)-vxc)*(x(i)-xc)
!!cc     &                 + ((vy(i)-vyc)*(y(i)-yc)))/r
!!                  vr = ((vx(i)-vxc)*xi
!!     &                 + ((vy(i)-vyc)*yi))/rmindd
!!                  IF (vr.LT.0.) THEN
!!cc                     vx(i) = vx(i) - vr*(x(i)-xc)/r
!!cc                     vy(i) = vy(i) - vr*(y(i)-yc)/r
!!                     vx(i) = vx(i) - vr*xi/rmindd
!!                     vy(i) = vy(i) - vr*yi/rmindd
!!                  ENDIF
!
!--Adjust to conserve angular momentum
!
!!                  IF (ifcor.EQ.1) THEN
!!                     vx(i) = (vx(i)-vxc)
!!     &                       * (y(i)-yc) / yi
!!     &                       + vxc*yc
!!                     vy(i) = (vy(i)-vyc)
!!     &                       * (x(i)-xc) / xi
!!     &                       + vyc*xc
!!                  ELSEIF (ifcor.EQ.2) THEN
!!                     vy(i) = (vy(i)-vyc)
!!     &                       * (z(i)-zc) / zi
!!     &                       + vyc*zc
!!                     vz(i) = (vz(i)-vzc)
!!     &                       * (y(i)-yc) / yi
!!     &                       + vzc*yc
!!                  ENDIF
                  x(i) = xi + xc
                  y(i) = yi + yc
                  z(i) = zi + zc
!--Have a particle rotate at the Keplerian speed
                  alpha = SQRT(pmass(isphcom)/rmindd)
                  xtild = xi/rmindd
                  ytild = yi/rmindd
                  vx(i) = -alpha * ytild + vxc
                  vy(i) =  alpha * xtild + vyc
                  vz(i) =  vzc
               ENDIF

               IF (ichan.NE.0) THEN
                  iouter = iouter + 1
                  idonebound = 1
               ENDIF
 270        CONTINUE

         ENDIF

      ELSE IF ( ibound.EQ.3 .OR. ibound.EQ.6) THEN
         IF (isphcom.EQ.0) THEN
            rmax2 = rmax * rmax
            DO 300 j = 1, nlst0
               i = list(j)
               IF (iphase(i).NE.0) GOTO 300
               r2 = (x(i)*x(i)+y(i)*y(i)+z(i)*z(i))
               IF ( r2.GT.rmax2 ) THEN
                  iouter = iouter + 1
                  idonebound = 1
                  r=sqrt(r2)
!
!--Put particle inside boundary
!
                  delta = (1.3-ran1(1))*0.5*h(i)
                  rmaxd = rmax - delta
                  zi = z(i)*rmaxd/r
                  yi = y(i)*rmaxd/r
                  xi = x(i)*rmaxd/r
!
!--Remove velocities perpendicular to boundary
!
                  vr = ( vx(i)*x(i) + vy(i)*y(i) &
                         + vz(i)*z(i) ) / r
                  IF ( vr.GT.0 ) THEN
                      vx(i) = vx(i) - vr * x(i) / r
                      vy(i) = vy(i) - vr * y(i) / r
                      vz(i) = vz(i) - vr * z(i) / r
                  ENDIF
!
!--Adjust to conserve angular momentum
!
                  IF (ifcor.EQ.1) THEN
                     vx(i) = vx(i) * y(i) / yi
                     vy(i) = vy(i) * x(i) / xi
                  ELSEIF (ifcor.EQ.2) THEN
                     vy(i) = vy(i) * z(i) / zi
                     vz(i) = vz(i) * y(i) / yi
                  ENDIF
                  z(i) = zi
                  y(i) = yi
                  x(i) = xi
               ENDIF
 300        CONTINUE
         ELSE
            xc = x(isphcom)
            yc = y(isphcom)
            zc = z(isphcom)
            vxc = vx(isphcom)
            vyc = vy(isphcom)
            vzc = vz(isphcom)
            rmind2 = rmind * rmind
            DO 310 j = 1, nlst0
               i = list(j)
               IF (iphase(i).NE.0) GOTO 310
               IF (i.LE.nptmass) GOTO 310

               xi = x(i) - xc
               yi = y(i) - yc
               zi = z(i) - zc
               vxi = vx(i) - vxc
               vyi = vy(i) - vyc
               vzi = vz(i) - vzc

               r2 = xi*xi+yi*yi+zi*zi
               IF ( r2.LT.rmind2 ) THEN
                  iinner = iinner + 1
                  idonebound = 1
                  r=SQRT(r2)
                  vi0 = SQRT(vxi*vxi+vyi*vyi+vzi*vzi)
!
!--Put particle outside boundary
!
!!                  delta = (1.3-ran1(1))*0.5*h(i)
                  delta = rmind * 0.01
                  rmindd = rmind + delta
                  xi = xi*rmindd/r
                  yi = yi*rmindd/r
                  zi = zi*rmindd/r
                  r = rmindd
!
!--Remove velocities perpendicular to boundary
!--Reflective boundary condition
                  vr = ( vxi*xi + vyi*yi &
                         + vzi*zi ) / r
                  IF ( vr.LT.0.0 ) THEN
                      vxi = vxi - vr * xi / r
                      vyi = vyi - vr * yi / r
                      vzi = vzi - vr * zi / r
!!                      vxi = vxi - 2.0 * vr * xi / r
!!                      vyi = vyi - 2.0 * vr * yi / r
!!                      vzi = vzi - 2.0 * vr * zi / r
                  ENDIF
!
!--Adjust to conserve angular momentum
                  IF (ifcor.EQ.1) THEN
                     vx(i) = vx(i) * y(i) / yi
                     vy(i) = vy(i) * x(i) / xi
                  ELSEIF (ifcor.EQ.2) THEN
                     vy(i) = vy(i) * z(i) / zi
                     vz(i) = vz(i) * y(i) / yi
                  ENDIF
!--Adjust to conserve energy
!!                  vi = SQRT(vxi*vxi+vyi*vyi+vzi*vzi)
!!                  vxi = vxi * vi0/vi
!!                  vyi = vyi * vi0/vi
!!                  vzi = vzi * vi0/vi
!--Have a particle rotate at the Keplerian speed
                  alpha = SQRT(pmass(isphcom)/rmindd)
                  vi = SQRT(vxi*vxi+vyi*vyi+vzi*vzi)
                  vxi = vxi * alpha/vi
                  vyi = vyi * alpha/vi
                  vzi = vzi * alpha/vi

                  x(i) = xi + xc
                  y(i) = yi + yc
                  z(i) = zi + zc
                  vx(i) = vxi + vxc
                  vy(i) = vyi + vyc
                  vz(i) = vzi + vzc
               ENDIF
 310        CONTINUE
         ENDIF
      ENDIF

      IF (ibound.EQ.2 .AND. iinner.NE.0) WRITE (iprint,99009) iinner
      IF (ibound.EQ.3 .AND. iinner.NE.0) WRITE (iprint,99009) iinner
      IF (iouter.NE.0) WRITE (iprint,99010) iouter

99009 FORMAT(' number of corrections inner boundary:',I6)
99010 FORMAT(' number of corrections outer boundary:',I6)

      END SUBROUTINE boundry
