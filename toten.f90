      SUBROUTINE toten
!************************************************************
!                                                           *
!  This routine computes all energies per unit total mass   *
!     kinetic, rotational, potential and internal.          *
!     It also computes the escapor mass.                    *
!                                                           *
!************************************************************

      use idims

      use part
      use densi
      use carac
      use cgas
      use gtime
      use fracg
      use ener
      use bodys
      use eosq
      use kerne
      use varet
      use logun
      use debug
      use phase

      implicit none

      INTEGER(I4B) :: i
      REAL(DP) :: gama1, drdt, dlnrdt, poteni, r2xy, r2yz, &
                  rscale, rvx, rvz, tinout, total, tpot, ttherm, &
                  vtot2, xi, yi, zi, vxi, vyi, vzi
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine toten')
!
!--Compute first mechanical energies and escapors
!
      gama1 = gamma - 1.
      tkin = 0.
      trotz = 0.
      trotx = 0.
      tgrav = 0.
      tpot = 0.
      tterm = 0.
      escap = 0.
!
!--Scaling factors
!
      CALL scaling(gt, rscale, drdt, dlnrdt)

      IF (idebug(1:5).EQ.'toten') THEN
         WRITE (iprint, 99002) rscale
99002    FORMAT (1X, 'rscale =',1PE12.5)
      ENDIF
      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            xi = x(i)
            yi = y(i)
            zi = z(i)
            vxi = vx(i)
            vyi = vy(i)
            vzi = vz(i)
!
!--Kinetic energy
!
            vtot2 = vxi**2 + vyi**2 + vzi**2
            tkin = tkin + pmass(i)*vtot2
!
!--Rotational energy around z
!
            r2xy = xi*xi + yi*yi
            rvz = xi*vyi - yi*vxi
            IF (r2xy.NE.0) trotz = trotz + pmass(i)*rvz*rvz/r2xy
!
!--Rotational energy around x
!
            r2yz = yi*yi + zi*zi
            rvx = yi*vzi - zi*vyi
            IF (r2yz.NE.0) trotx = trotx + pmass(i)*rvx*rvx/r2yz
!
!--Potential energy
!
            poteni = pmass(i)*(poten(i) + dgrav(i))/rscale
            tgrav = tgrav + poteni
            IF (idebug(1:5).EQ.'toten') THEN
               IF (i.LE.10) THEN
                  WRITE (iprint, 99003) i,pmass(i),poten(i), &
                                        dgrav(i),poteni
99003             FORMAT (1X, 'i=',I3,': pmass=',1PE9.2, &
                          ', poten=',1PE9.2,', dgrav=',1PE9.2, &
                          ' --> poteni=',1PE9.2)
               ENDIF
            ENDIF
!
!--Escapors
!
            IF (x(i)*vxi + y(i)*vyi + z(i)*vzi.LT.0.) vtot2 = 0.
            tinout = 0.5*pmass(i)*vtot2
            tpot = poteni
            IF (varsta.EQ.'intener') THEN
               ttherm = u(i)*pmass(i)
            ELSE
               ttherm = pmass(i)*pr(i)/(gama1*rho(i))
            ENDIF
            total = tinout + tpot + ttherm
            IF (total.GT.0.) escap = escap + pmass(i)
         ENDIF
      END DO
!
!--Normalisations
!
      tkin = 0.5*tkin
      trotz = 0.5*trotz
      trotx = 0.5*trotx
      tgrav = 0.5*tgrav
      escap = escap
!
!--Thermal energy
!
      IF (varsta.NE.'entropy') THEN
!
!--Variable of state is specific internal energy
!
         DO i = 1, npart
            IF (iphase(i).EQ.0) tterm = tterm + pmass(i)*u(i)
         END DO
!
!--Variable of state is specific entropy
!
      ELSEIF (gama1.EQ.0.) THEN
         tterm = 1.5*u(1)
      ELSE
         DO i = 1, npart
            IF (iphase(i).EQ.0) tterm = &
                        tterm + pmass(i)*u(i)*rho(i)**gama1/gama1
         END DO
         tterm = tterm
      ENDIF

      IF (idebug(1:5).EQ.'toten') THEN
         WRITE (iprint, 99004) tkin, trotz, trotx, tgrav, tterm
99004    FORMAT (1X, 'tkin =',1PE12.5,', trotz =',1PE12.5,', trotx =', &
                 1PE12.5,', tgrav =',1PE12.5,', tterm =',1PE12.5)
      ENDIF

      END SUBROUTINE toten


