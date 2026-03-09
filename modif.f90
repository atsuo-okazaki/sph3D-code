      SUBROUTINE modif
!************************************************************
!                                                           *
!  This subroutine allows modifications to be made during   *
!     the transfer of dumps.                                *
!                                                           *
!************************************************************

      use idims

      use part
      use carac
      use kerne
      use densi
      use typef
      use trans
      use expan
      use rotat
      use units
      use logun
      use debug

      implicit none

      INTEGER(I4B) :: i, ncount, need
      REAL(DP) :: d, d2, rsup, rsup2, rx, step, totm, totmin, &
               vrad
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine modif')
!
!--Find particles to change
!
      ncount = 0
      step = 0.005
      rsup = 0.000
      need = npart*frac + 1
!
!--Find radius containing particle fraction
!
 100  DO i = 1, npart
         d = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
         IF (d.LT.rsup) THEN
            IF (ncount + 1.GT.need) THEN
               rsup = rsup - step
               step = step/1.5
               IF (step.GE.1.E-4) GOTO 300
               GOTO 400
            ENDIF
            ncount = ncount + 1
         ENDIF
      END DO
      IF (ncount.EQ.need) GOTO 400
 300  rsup = rsup + step
      ncount = 0
      GOTO 100
!
!--Make the change
!
 400  rsup2 = rsup**2
      WRITE (*, *) 'rsup=', rsup
      totm = 0.
      totmin = 0.
      DO 500 i = 1, npart
         d2 = x(i)**2 + y(i)**2 + z(i)**2
         IF (ichang.EQ.1) THEN
            IF (d2.LE.rsup2) THEN
               u(i) = u(i) + energc
               totmin = totmin + pmass(i)
            ENDIF
            totm = totm + pmass(i)
            GOTO 500
         ENDIF
         IF (ichang.EQ.2) THEN
            IF (d2.LE.rsup2) THEN
               u(i) = u(i) + energc
               totmin = totmin + pmass(i)
            ENDIF
            totm = totm + pmass(i)
            rx = SQRT(x(i)**2 + y(i)**2 + z(i)**2)
            vrad = vexpan*( - 13.0*rx**10 - 7.0*rx**5)/20.0
            vx(i) = vrad*x(i)/rnorm - omeg0*y(i)
            vy(i) = vrad*y(i)/rnorm + omeg0*x(i)
            vz(i) = vrad*z(i)/rnorm
         ENDIF
 500  CONTINUE
      WRITE (*, *) 'mass ratio ', totmin/totm

      END SUBROUTINE modif
