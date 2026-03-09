      SUBROUTINE chanref(idir)
!************************************************************
!                                                           *
!  This subroutine transforms the various quantities        *
!     between two frame of references.                      *
!                                                           *
!************************************************************

      use idims

      use constants
      use part
      use ener
      use densi
      use rotat
      use gtime
      use logun
      use debug
      use typef

      implicit none

      INTEGER(I4B) :: idir, i
      REAL(DP) :: drdt, dlnrdt, or, or1, r2, rscale, sdens, sergg
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine chanref')
!
!--Scaling factors
!
      CALL scaling(gt, rscale, drdt, dlnrdt)
      r2 = drdt/rscale**2
      sdens = 1./rscale**3
      sergg = 1./rscale
      or = omeg0*rscale
      or1 = omeg0/rscale
!
!--Transform into original unit system
!
      IF ( idir.NE.1 ) THEN
         IF ( ifcor.EQ.1 ) THEN
            DO 50 i = 1, npart
               vx(i) = vx(i)/rscale - x(i)*r2 + or1*y(i)
               vy(i) = vy(i)/rscale - y(i)*r2 - or1*x(i)
               vz(i) = vz(i)/rscale - z(i)*r2

               x(i) = x(i)/rscale
               y(i) = y(i)/rscale
               z(i) = z(i)/rscale
               h(i) = h(i)/rscale

               rho(i) = rho(i)/sdens
               u(i) = u(i)/sergg
 50         CONTINUE
         ELSE IF ( ifcor.EQ.2 ) THEN
            DO 80 i = 1, npart
               vx(i) = vx(i)/rscale - x(i)*r2
               vy(i) = vy(i)/rscale - y(i)*r2 + or1*z(i)
               vz(i) = vz(i)/rscale - z(i)*r2 - or1*y(i)

               x(i) = x(i)/rscale
               y(i) = y(i)/rscale
               z(i) = z(i)/rscale
               h(i) = h(i)/rscale

               rho(i) = rho(i)/sdens
               u(i) = u(i)/sergg
 80         CONTINUE
         ENDIF
         RETURN
!
!--Transform into expanding frame units
!
      ELSE
         IF ( ifcor.EQ.1 ) THEN
            DO 100 i = 1, npart
               vx(i) = rscale*vx(i) + x(i)*drdt - or*y(i)
               vy(i) = rscale*vy(i) + y(i)*drdt + or*x(i)
               vz(i) = rscale*vz(i) + z(i)*drdt

               x(i) = x(i)*rscale
               y(i) = y(i)*rscale
               z(i) = z(i)*rscale
               h(i) = h(i)*rscale

               rho(i) = rho(i)*sdens
               u(i) = u(i)*sergg
               poten(i) = poten(i)*sergg
 100        CONTINUE
         ELSE IF ( ifcor.EQ.2 ) THEN
            DO 120 i = 1, npart
               vx(i) = rscale*vx(i) + x(i)*drdt
               vy(i) = rscale*vy(i) + y(i)*drdt - or*z(i)
               vz(i) = rscale*vz(i) + z(i)*drdt + or*y(i)

               x(i) = x(i)*rscale
               y(i) = y(i)*rscale
               z(i) = z(i)*rscale
               h(i) = h(i)*rscale

               rho(i) = rho(i)*sdens
               u(i) = u(i)*sergg
               poten(i) = poten(i)*sergg
 120        CONTINUE
         ENDIF
         RETURN
      ENDIF

      END SUBROUTINE chanref
