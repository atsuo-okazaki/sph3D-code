      SUBROUTINE ghostp1(npart, x, y, z, vx, vy, vz, u, h)
!************************************************************
!                                                           *
!  This subroutine computes the list of ghost particles for *
!     treating the boundaries.                              *
!                                                           *
!************************************************************

      use idims

      use carac
      use ghost
      use densi
      use rbnd
      use diskbd
      use logun
      use debug
      use phase
      use kerne

      implicit none

      INTEGER(I4B) :: npart, i, nptot, ntot
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: vx(idim), vy(idim), vz(idim), h(idim), u(idim)
      REAL(DP) :: delta, delta2, dxmax, dxmax2, dxmin, dxmin2, &
               dymax, dymax2, dymin, dymin2, dzmax, dzmax2, dzmin, dzmin2, &
               hi, pmassi, radius2, radk2, rhoi, ui, vxi, vyi, vzi, &
               xi, yi, zi
      CHARACTER(len=7) :: where='ghostp1'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine ghostp1')

      nghost = 0
!
!--Find ghost particles (for all particles within radkernel*h of boundary)
!
      DO 200 i = 1, npart
         hasghost(i) = 0
         IF (iphase(i).NE.0) GOTO 200
         xi = x(i)
         yi = y(i)
         zi = z(i)

         vxi = vx(i)
         vyi = vy(i)
         vzi = vz(i)
         hi = h(i)
         pmassi = pmass(i)
         rhoi = rho(i)
         ui = u(i)
         delta = 0.1*hi
!
!--X axis
!
         dxmin = (xi - xmin)/hi
         IF (dxmin.GT.delta .AND. dxmin.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

         dxmax = (xmax - xi)/hi
         IF (dxmax.GT.delta .AND. dxmax.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
!
!--Y axis
!
         dymin = (yi - ymin)/hi
         IF (dymin.GT.delta .AND. dymin.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

         dymax = (ymax - yi)/hi
         IF (dymax.GT.delta .AND. dymax.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
!
!--Z axis
!
         dzmin = (zi - zmin)/hi
         IF (dzmin.GT.delta .AND. dzmin.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

         dzmax = (zmax - zi)/hi
         IF (dzmax.GT.delta .AND. dzmax.LT.radkernel) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
!
!--Edges
!
         dxmin2 = dxmin*dxmin
         dymin2 = dymin*dymin
         dzmin2 = dzmin*dzmin
         dxmax2 = dxmax*dxmax
         dymax2 = dymax*dymax
         dzmax2 = dzmax*dzmax
         delta2 = delta*delta

         radius2 = dxmin2 + dymin2
         radk2 = radkernel*radkernel
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                          dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                          dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                          dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                          dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

         radius2 = dxmin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                          dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                          dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                          dzmin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                          dzmax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

         radius2 = dzmin2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. &
                          dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmin2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmin.GT.delta .AND. &
                          dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmax2 + dymin2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. &
                          dymin.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dzmax2 + dymax2
         IF (radius2.GT.delta2 .AND. dzmax.GT.delta .AND. &
                          dymax.GT.delta .AND. radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
!
!--Corners
!
         radius2 = dxmin2 + dymin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                dymin.GT.delta .AND. dzmin.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                dymin.GT.delta .AND. dzmax.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                dymax.GT.delta .AND. dzmin.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmin2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmin.GT.delta .AND. &
                dymax.GT.delta .AND. dzmax.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmin - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                dymin.GT.delta .AND. dzmin.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymin2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                dymin.GT.delta .AND. dzmax.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymin - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmin2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                dymax.GT.delta .AND. dzmin.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmin - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF
         radius2 = dxmax2 + dymax2 + dzmax2
         IF (radius2.GT.delta2 .AND. dxmax.GT.delta .AND. &
                dymax.GT.delta .AND. dzmax.GT.delta .AND. &
                radius2.LT.radk2) THEN
            hasghost(i) = 1
            nghost = nghost + 1
            nptot = MIN0(npart + nghost, idim)
            ireal(nptot) = i
            x(nptot) = 2.0*xmax - xi
            y(nptot) = 2.0*ymax - yi
            z(nptot) = 2.0*zmax - zi
            vx(nptot) = -vxi
            vy(nptot) = -vyi
            vz(nptot) = -vzi
            h(nptot) = hi
            pmass(nptot) = pmassi
            rho(nptot) = rhoi
            u(nptot) = ui
            iphase(nptot) = 0
         ENDIF

 200  CONTINUE

      ntot = npart + nghost

      WRITE (iprint, *) 'npart, nghost', npart, nghost
      IF (ntot.GT.idim) CALL error(where, ntot)

      END SUBROUTINE ghostp1
