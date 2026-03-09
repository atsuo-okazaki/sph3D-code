      SUBROUTINE addump
!************************************************************
!                                                           *
!  This routine reads two existing dumps and create a third *
!  one out of it.                                           *
!                                                           *
!************************************************************

      use idims

      use units
      use constants
      use part
      use carac
      use cgas
      use densi
      use varet
      use typef
      use new
      use bodys
      use ener
      use kerne
      use logun
      use debug
      use polyk2
      use sphcom

      CHARACTER(len=7) ::  where='addump'
      CHARACTER(len=1) :: iok

!--allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine addump')

      fmas1 = 0.
      fmas2 = 0.
      WRITE (*, 99002)
99002 FORMAT (' PARTICLE SET UP', //)
      WRITE (*, 99003)
99003 FORMAT (' self-gravity included? (y/n)')
      READ (*, 99004) iok
99004 FORMAT (A1)
      igphi = 0
      IF ( iok.EQ.'y' ) igphi = 1
      WRITE (*, 99005)
99005 FORMAT (' give the position and the velocity of the cm of body 1')
      READ (*, *) xx1, yy1, zz1m, vvx1, vvy1, vvz1
      WRITE (*, 99006)
99006 FORMAT (' give the position and the velocity of the cm of body 2')
      READ (*, *) xx2, yy2, zz2, vvx2, vvy2, vvz2
      WRITE (*, 99007)
99007 FORMAT (' do you want to rotate body 2? (y/n)')
      READ (*, 99004) iok
      rotang = 0.
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99008)
99008    FORMAT (' give rotation angle (rad.)')
         READ (*, *) rotang
      ENDIF
!
!--read first object
!
      WRITE (*, 99009)
99009 FORMAT (' reading first file and storing data')
      READ (idisk1, END=300)  udist, umass, utime, &
                    np1, n11, n12, t1, gamma1, rhozero1, RK21, &
                             (h(i), i=1, np1), escap1, tkin1, &
                             tgrav1, tterm1, (x(i), i=1, np1), &
                             (y(i), i=1, np1), (z(i), i=1, np1), &
                             (vx(i), i=1, np1), &
                             (vy(i), i=1, np1), (vz(i), i=1, np1), &
                             (u(i), i=1, np1), (pmass(i), i=1, np1), &
                             (rho(i), i=1, np1)
!
!--put first object into place
!
      DO 100 j = 1, np1
         fmas1 = fmas1 + pmass(j)
         x(j) = x(j) + xx1
         y(j) = y(j) + yy1
         z(j) = z(j) + zz1
         vx(j) = vvx1
         vy(j) = vvy1
         vz(j) = vvz1
         dgrav(j) = 0.
 100  CONTINUE
!
!--read second dump
!
      WRITE (*, 99010)
99010 FORMAT (' reading second file and storing data')
      READ (idisk2, END=300)  udist, umass, utime, &
                             np2, n21, n22, t2, gamma2, rhozero2, RK22, &
                             (h(i+np1), i=1, np2), escap2, tkin2, &
                             tgrav2, tterm2, (x(i+np1), i=1, np2), &
                             (y(i+np1), i=1, np2), (z(i+np1), i=1, np2), &
                             (vx(i+np1), i=1, np2), &
                             (vy(i+np1), i=1, np2), &
                             (vz(i+np1), i=1, np2), (u(i+np1), i=1, np2) &
                             , (pmass(i+np1), i=1, np2), &
                             (rho(i+np1), i=1, np2)
!
!--rotate second object if needed
!
      crotang = cos(rotang)
      srotang = sin(rotang)
      DO 200 j = np1 + 1, np1 + np2
         fmas2 = fmas2 + pmass(j)
         x(j) = crotang*x(j) + srotang*y(j) + xx2
         y(j) = -srotang*x(j) + crotang*y(j) + yy2
         z(j) = z(j) + zz2
         vx(j) = vvx2
         vy(j) = vvy2
         vz(j) = vvz2
         dgrav(j) = 0.
 200  CONTINUE

      WRITE (*, 99011)
99011 FORMAT (' what is the equation of state variable:', /, &
              ' specific internal energy :  intener', /, &
              ' specific entropy         :  entropy')
      READ (*, 99012) varsta
99012 FORMAT (A7)
      WRITE (*, 99013)
99013 FORMAT (' is the equation of state a gamma-law? (y/n)')
      READ (*, 99004) iok
      IF ( iok.EQ.'y' ) THEN
         WRITE (*, 99014) gamma1, gamma2
99014    FORMAT (' what is the value of the adiabatic index gamma?', /, &
                 ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) gamma
         WRITE (*, 98014) RK21, RK22
98014    FORMAT (' what is the value of the adiabatic constant RK2?', /, &
                 ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) RK2
         WRITE (*, 98015) rhozero1, rhozero2
98015    FORMAT (' what is the value of the  constant rhozero?', /, &
                 ' object 1 had ', 1PE12.5, ' object 2 had ', 1PE12.5)
         READ (*, *) rhozero
      ENDIF
      npart = np1 + np2
      n1 = np1
      n2 = np2
      WRITE (*, 99015) n1, n2, npart
99015 FORMAT (' set up completed', /, &
              ' particles in object 1          :', I6, /, &
              ' particles in object 2          :', I6, /, &
              ' total number of particles used :', I6)
      WRITE (*, 99016)
99016 FORMAT (//, ' END SET UP')

      IF ( idebug(1:6).EQ.'addump' ) THEN
         WRITE (iprint, 99017) (x(i), i=1, npart)
         WRITE (iprint, 99017) (y(i), i=1, npart)
         WRITE (iprint, 99017) (z(i), i=1, npart)
         WRITE (iprint, 99017) (vx(i), i=1, npart)
         WRITE (iprint, 99017) (vy(i), i=1, npart)
         WRITE (iprint, 99017) (vz(i), i=1, npart)
         WRITE (iprint, 99017) (u(i), i=1, npart)
         WRITE (iprint, 99017) (pmass(i), i=1, npart)
         WRITE (iprint, 99017) (rho(i), i=1, npart)
99017    FORMAT (1X, 5(1PE12.5,1X))
      ENDIF
      GOTO 400

 300  CALL error(where, 1)

 400  RETURN
      END SUBROUTINE addump
