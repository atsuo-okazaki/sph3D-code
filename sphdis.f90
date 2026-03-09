      SUBROUTINE sphdis(igeom, idist, np, h1, facx, facy, facz, &
                         delx, dely, nx, ny, nz)
!************************************************************
!                                                           *
!  This subroutine positions particles in a spherical       *
!     distribution                                          *
!                                                           *
!************************************************************

      use idims

      use constants
      use part
      use rbnd
      use diskbd
      use debug
      use logun
      use maspres
      use ptmass

      implicit none

      INTEGER(I4B) :: igeom, idist, np, nx, ny, nz, i, j, k, iprofr, &
               nspace
      REAL(DP) :: ran1, velweight, probav(idim)
      REAL(DP) :: vel(64,64,64)
      REAL(DP) :: h1, facx, facy, facz, delx, dely, denprof, probavn, &
               probr, profr, r1, r2, rc2, rcyl2, rm2n, rmax2, rmind2, &
               rxy2, theta1, theta2, third, xmax5, ymax5, zmax5, &
               rm20, velmax, weight, x1, y1, z1, xi, yi, zi
      CHARACTER(len=20) :: filename
      CHARACTER(len=1) :: prof, iok
!
!--Allow for tracing flow
!
      third = 1./3.
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine sphdis ')
!
!--Set Centrally Condensed Spherical Distribution
!
      IF ((idist.EQ.1) .OR. (idist.EQ.2)) THEN
         CALL unifdis(igeom, idist, np, h1, facx, facy, facz, &
                         delx, dely, nx, ny, nz, 0)
         rcyl2 = rcyl * rcyl
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax
         WRITE (*, 88001)
88001    FORMAT (' enter profile along r ',/, &
                 '                  r^-1 : 1',/, &
                 '                  r^-2 : 2',/, &
                 '           exponential : e',/, &
                 '             velocity  : v')
         READ (*, 99004) prof
         IF (prof.NE.'e') THEN
            IF (prof.EQ.'1') iprofr = 1
            IF (prof.EQ.'2') iprofr = 2
            profr = iprofr/2.
         ENDIF
!
!--Deform particle grid
!
         IF (iprofr.EQ.1) THEN
            DO i = 1, npart
               xi = x(i)
               yi = y(i)
               zi = z(i)
               rxy2 = xi*xi + yi*yi
               r2 = rxy2 + zi*zi
               theta1 = ATAN2(yi,xi)
               theta2 = ATAN2(zi,SQRT(rxy2))
               r1 = r2**0.75
               x(i) = r1*COS(theta1)*COS(theta2)
               y(i) = r1*SIN(theta1)*COS(theta2)
               z(i) = r1*SIN(theta2)
            END DO
         ELSEIF (iprofr.EQ.2 .OR. prof.EQ.'e' .OR. prof.EQ.'v') THEN
            WRITE (*,99002)
99002       FORMAT ('NOT IMPLEMENTED')
            CALL quit
         ENDIF
      ELSE
         npart = np + nptmass
         rcyl2 = rcyl * rcyl
         rmind2 = rmind * rmind
         rmax2 = rmax * rmax
         WRITE (*, 88001)
         READ (*, 99004) prof
         IF (prof.NE.'e' .AND. prof.NE.'v') THEN
            IF (prof.EQ.'1') iprofr = 1
            IF (prof.EQ.'2') iprofr = 2
            profr = iprofr/2.
         ENDIF

         IF (prof.EQ.'v') THEN
            WRITE (*, 88401)
88401       FORMAT ('Enter name of file')
            READ (*,99010) filename
99010       FORMAT(A20)

            WRITE (*,*) 'Enter size of velocity files (e.g.N=32^3)'
            READ (*,*) nspace

            OPEN (45,FILE=filename,FORM='unformatted')
            READ (45) (((vel(i,j,k), i=1,nspace),j=1,nspace), &
                 k=1,nspace)
            CLOSE (45)

            velmax = 0.
            DO k = 1, nspace
               DO j = 1, nspace
                  DO i = 1, nspace
                     velmax = MAX(velmax, vel(i,j,k))
                  END DO
              END DO
            END DO
            WRITE (*,*) 'velmax ',velmax
         ENDIF

         rm20 = rmax2 * 0.58 * 0.58
         xmax5 = xmax * 0.5
         ymax5 = ymax * 0.5
         zmax5 = zmax * 0.5
         probavn = 0.

         DO i = nptmass + 1, npart
            IF (MOD(i,1000).EQ.0) write (*,*) i
 100        x1 = 2.*(xmax*ran1(1) - xmax5)
            y1 = 2.*(ymax*ran1(1) - ymax5)
            z1 = 2.*(zmax*ran1(1) - zmax5)

            IF (prof.EQ.'v') THEN
               weight = velweight(nspace,x1,y1,z1,vel)/velmax
               IF (ran1(1).GT.weight) GOTO 100
            ENDIF

            rc2 = x1*x1 + y1*y1
            r2 = rc2 + z1*z1
            IF ((igeom.EQ.2) .AND. (rc2.GT.rcyl2)) GOTO 100
            IF ((igeom.EQ.3) .AND. (r2.GT.rmax2)) GOTO 100
            IF ((igeom.EQ.4) .AND. ((rc2.GT.rcyl2) .OR. &
                 (rc2.LT.rmind2))) GOTO 100
            IF ((igeom.EQ.7) .AND. ((r2.GT.rmax2) &
                 .OR. (r2.LT.rmind2))) GOTO 100

            IF (prof.NE.'e') THEN
               probr = 0.01 * (rmax2/r2) ** profr
            ELSE
               probr = 0.05 * 20.0 * EXP(-1.0*r2/rm20)
            ENDIF

            IF (ran1(1).GT.probr) GOTO 100

            probav(i) = probr
            probavn = probavn + probav(i)
            x(i) = x1
            y(i) = y1
            z(i) = z1
         END DO

         probavn = probavn / DBLE(npart - nptmass)
         WRITE (*,*) probavn
         DO i = nptmass + 1, npart
            denprof = (probavn / probav(i)) ** third
            h(i) = h1 * denprof
         END DO
      ENDIF

      DO i = nptmass + 1, npart
         disfrac(i) = 1.0
      END DO

      END  SUBROUTINE sphdis
