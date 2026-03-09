      SUBROUTINE header(where)
!************************************************************
!                                                           *
!  This routine writes on first page of listing the value   *
!     of all variable defined at the start of the run.      *
!                                                           *
!************************************************************

      use idims

      use constants
      use part
      use units
      use typef
      use gtime
      use dissi
      use numpa
      use bodys
      use rbnd
      use diskbd
      use varet
      use rotat
      use btree
      use integ
      use kerne
      use logun
      use files
      use polyk2
      use cgas
      use vargam
      use pres
      use soft
      use ptmass
      use nextmpt
      use active
      use ptsoftx
      use phase
      use binary
      use ptbin
      use physeos
      use useles
      use sphcom

      implicit none

      INTEGER(I4B) :: i
      REAL(DP) :: angmom, rscale, drdt, dlnrdt
      REAL(DP) :: velo, tcomp, tff, rhocrt, rhocrt2, rhocrt3, &
             rp1, rc1, rc2

      CHARACTER(len=35) :: var
      CHARACTER(len=7) :: where
      CHARACTER(len=7) :: where2='header'
!
!--Scaling factors
!
      CALL scaling(gt, rscale, drdt, dlnrdt)
!
!--Write units
!
      angmom = umass*dble(udist**2)/dble(utime)
      velo = udist/utime
      WRITE (iprint, 99001) umass, udist, udens, utime, velo, uergg, &
                            angmom
99001 FORMAT (//, ' The computations are done in the following units', &
              /, ' units of :  mass       :', 1PD12.4, &
              '   distance    :', ES12.4, /, &
              '             density    :', ES12.4, '   time        :', &
              ES12.4, /, '             velocity   :', ES12.4, &
              '   energy/mass :', ES12.4, /, &
              '             ang. mom.  :', ES12.4, //)
!
!--Write options
!
      IF (varsta.EQ.'entropy') THEN
         var = 'specific entropy'
      ELSE
         var = 'specific internal energy'
      ENDIF
      WRITE (iprint, 99002) var
99002 FORMAT (' Variable of state used : density and ', A35, /)
!
!--Mass fractions
!
      WRITE (iprint, 99003) npart, nactive, (fmas1 + fmas2), n1, &
              fmas1, n2, fmas2
99003 FORMAT (/, ' Total number of particles used : ', I8, /, &
              ' Number of active particles     : ', I8, /, &
              ' Total mass                     :', ES12.4, /, &
              ' Distribution : object 1 :  number of particles : ', I8, &
              /, '                            mass of object      : ', &
              ES12.4, /, &
              '                object 2 :  number of particles : ', I8, &
              /, '                            mass of object      : ', &
              ES12.4, /)

!
!--Calculate tff
!
      tcomp = SQRT((3 * pi) / (32 * rhozerox))
!c      tcomp = SQRT((3 * pi) / (32 * rhozero))
      tff = tcomp * utime
      WRITE (iprint, 99004) tff, tcomp
99004 FORMAT(/,' The free fall time is : ',ES14.7,/, &
                  ' in computational units: ',ES14.7)

      IF (where(1:6).NE.'newrun') THEN
!
!--Options of the code
!
         WRITE (iprint, 99005) encal, igrp, igphi, ifsvi, ifcor, &
                 ibound, iexf, iener, ichoc, iexpan, damp

99005    FORMAT (/, ' The following options were set :', /, &
                 ' energy calculation     : ', A1, /, &
                 ' pressure gradients     : ', I2, &
                 '   self-gravity     : ', I2, /, &
                 ' artificial viscosity   : ', I2, &
                 '   coriolis         : ', I2, /, &
                 ' boundary type          : ', I2, &
                 '   external force   : ', I2, /, &
                 ' energy conservation    : ', I2, &
                 '   shock heating    : ', I2, /, &
                 ' homologous exp.        : ', I2, &
                 '   general damping  : ', F4.1, /)
!
!--Gravity calculations
!
         IF (igphi.EQ.1) THEN
            IF (igrape.EQ.1) THEN
               WRITE (iprint, 98001)
98001          FORMAT (' Gravity calculated using: GRAPE board')
            ELSEIF (igrape.EQ.0) THEN
               WRITE (iprint, 98002)
98002          FORMAT (' Gravity calculated using: Binary Tree')
            ELSE
               CALL error(where2, 1)
            ENDIF

            IF (isoft.EQ.1) THEN
               WRITE (iprint, 98003) psoft
98003          FORMAT ('                           Softening = ', &
                    ES12.5)
            ELSEIF (isoft.EQ.0) THEN
               WRITE (iprint, 98004)
98004          FORMAT ('                           Kernel Softening')
            ELSE
               CALL error(where2, 2)
            ENDIF
         ENDIF
!
!--Rotation of expansion velocity
!
         IF (iexpan.NE.0) THEN
            WRITE (iprint, 99006) drdt*udist/utime
99006       FORMAT (/, &
                 ' Calculations done in a frame homologously expanding' &
                 , ' at ', ES12.5, ' cm/s', /)
         ENDIF
         IF (ifcor.GT.0) THEN
            WRITE (iprint, 99007) omeg0/utime
99007       FORMAT (/, &
                 ' Calculations done in a frame in uniform rotation', &
                 ' at ', ES12.5, ' /s ', /)
         ENDIF
!
!--Boundaries
!
         IF (ibound.EQ.1) THEN
            WRITE (iprint, 99008) xmin, xmax, ymin, ymax, zmin, zmax
99008       FORMAT (/,' Boundary type : reflective ', /, &
                    ' position :  cartesian    :  xmin : ', F7.3, &
                    '  xmax :', F7.3, /, &
                    '                             ymin : ', F7.3, &
                    '  ymax :', F7.3, /, &
                    '                             zmin : ', F7.3, &
                    '  zmax :', F7.3, /)
         ELSEIF (ibound.EQ.2) THEN
            WRITE (iprint, 99009) rcyl, zmin, zmax
99009       FORMAT(/,' Boundary type : cylindrical reflective ',/, &
                       '           rcyl : ', F7.3, /, &
                       '           zmin : ', F7.3, /, &
                       '           zmax : ', F7.3, /)
            IF (rmind.NE.0) WRITE (iprint, 99010) rmind
99010       FORMAT (/,' accretion disk inner boundary :',F7.3,/)
         ELSEIF (ibound.EQ.3) THEN
            IF (isphcom.EQ.0) THEN
               WRITE (iprint, 99011) 'rmax', rmax
99011          FORMAT (/,' Boundary type : spherical reflective ',/, &
                        '         ',A5,' : ', F7.3, /)
            ELSE
               WRITE (iprint, 99011) 'rmind', rmind
            ENDIF
         ELSEIF (ibound.EQ.7) THEN
            WRITE (iprint,99501) pext
99501       FORMAT (/,' Boundary type: constant pressure',/, &
                 '          pressure: ',ES12.4,/)
!c            exttemp = 2.0/3.0*extu/Rg*gmw*uergg
!c            WRITE (iprint,99501) extdens,extu,exttemp
!c99501       FORMAT (/,' Boundary type: constant pressure',/,
!c     &           '           density: ',ES12.4,/,
!c     &           '   internal energy: ',ES12.4,/,
!c     &           '       temperature: ',ES12.4,/)
         ELSEIF (ibound.EQ.8 .OR. ibound.GE.90) THEN
            WRITE (iprint,99502) deadbound, specang, fractan, fracradial
99502       FORMAT (/,' Boundary type: dead particle',/, &
                 '            radius: ',ES12.4,/, &
                 '   spec. ang. mom.: ',ES12.4,/, &
                 '  frac. tangential: ',ES12.4,/, &
                 '  frac. radial    : ',ES12.4,/)
         ENDIF
!
!--Critical densities for variable gamma
!
         WRITE (iprint, 99504) gamma
99504    FORMAT (' gamma                : ', ES12.3)
         WRITE (iprint, 99503) gmw
99503    FORMAT (/,' Mean molecular weight: ', ES12.3, /)

         IF (encal.EQ.'v') THEN
            rhocrt = rhocrit * udens
            rhocrt2 = rhocrit2 * udens
            rhocrt3 = rhocrit3 * udens
            WRITE (iprint, 99012) gam, rhocrt, gam, gamdh, rhocrt2, &
                   gamdh, gamah, rhocrt3
99012       FORMAT (/,' Critical densities for changing gamma ', /, &
            '           from 1.00    to ',F8.3,': ', ES12.3, /, &
            '           from ',F8.3,'to ',F8.3,': ', ES12.3, /, &
            '           from ',F8.3,'to ',F8.3,': ', ES12.3, //)
         ELSEIF (encal.EQ.'x') THEN
            rp1 = rhophys1 * udens
            rc1 = rhochange1 * udens
            rc2 = rhochange2 * udens
            WRITE (iprint, 99040) gamphys1, rp1, gamphys1, &
                   gamphys2, rc1, gamphys2, gamphys3, rc2
99040       FORMAT (/,' Physical equation of state ', /, &
            '           from 1.00    to ',F8.3,': ', ES12.3, /, &
            '           from ',F8.3,'to ',F8.3,': ', ES12.3, /, &
            '           from ',F8.3,'to ',F8.3,': ', ES12.3, //)
         ENDIF
!
!--Print out constants used for integration
!
         WRITE (iprint, 99013) alpha, beta, acc, tol, tolptm, tolh
99013    FORMAT (/,' Numerical constants used in this run :', /, &
                 ' artificial viscosity  alpha     : ', ES12.3, /, &
                 '                       beta      : ', ES12.3, /, &
                 ' binary tree accuracy param.     : ', ES12.3, /, &
                 ' RK2 tolerance - gas             : ', ES12.3, /, &
                 '               - point masses    : ', ES12.3, /, &
                 '               - smoothing length: ', ES12.3, //)
!
!--Print out massive point mass details
!
         IF (iptmass.NE.0) THEN
            WRITE(iprint,99014) iptmass, radcrit, ptmcrit, rhocrea
         ELSE
            WRITE(iprint,99015)
         ENDIF
99014    FORMAT (' Point mass creation ALLOWED, type ', I2, /, &
                 '  minimum creation radius      : ', ES12.3, /, &
                 '  creation density (in rhozero): ', ES12.3, /, &
                 '                (in code units): ', ES12.3, /)
99015    FORMAT (' Point mass creation NOT ALLOWED')

         IF (iptmass.NE.0.OR.nptmass.NE.0) THEN
            IF (iptintree) THEN
               WRITE (iprint,99147)
99147                    FORMAT (' Point masses done in TREE')
            ELSEIF (.NOT.iptintree) THEN
               WRITE (iprint,99148)
99148                    FORMAT (' Point masses done in GFORSPT')
            ELSE
               WRITE (iprint,99149)
99149                    FORMAT (' ERROR - Ptmasses in header ')
               CALL quit
            ENDIF

            WRITE(iprint,99016) nptmass
99016       FORMAT (' Number of point masses        : ', I4)

            IF (iaccevol.EQ.'v') THEN
               WRITE (iprint,99150) accfac
99150          FORMAT ('  Accretion radii: VARIABLE ROCHE ',ES12.5)
            ELSEIF (iaccevol.EQ.'s') THEN
               WRITE (iprint,99151) accfac
99151       FORMAT ('  Accretion radii: VARIABLE SEPARATION ',ES12.5)
            ELSE
               WRITE (iprint,99152)
99152          FORMAT ('  Accretion radii: FIXED ')
            ENDIF

            IF (hacc.GT.0.) THEN
               WRITE (iprint,99017) hacc
            ELSE
               CALL error(where2,1)
            END IF
99017       FORMAT ('  Outer accretion radius       : ', ES12.3)

            IF (haccall.GT.0.) THEN
               WRITE (iprint,99018) haccall
            ELSE
               CALL error(where2,2)
            END IF
99018       FORMAT ('  Inner accretion radius       : ', ES12.3)

            WRITE (iprint,99019) iptsoft
            IF (iptsoft.NE.0) WRITE (iprint,99020) ptsoft
99019       FORMAT ('  Gravity ptmass softening     : ', I2)
99020       FORMAT ('  Softening ptmass parameter   : ', ES12.3)

            IF (nptmass.NE.0) THEN
               DO i = 1, nptmass
                  WRITE (iprint,99021) i,iphase(listpm(i)),h(listpm(i))
99021             FORMAT (' Point mass: ',I4, ' type: ',I1, &
                       ' hacc: ',ES12.3)
               END DO
            ENDIF
         ENDIF

      ENDIF
!
!--Smoothing length max min
!
      IF (hmin.NE.0.0) THEN
         WRITE (iprint, 99114) hmin
99114        FORMAT (/,' MINIMUM SMOOTHING LENGTH = ', ES12.3, //)
      ENDIF
      IF (hmaximum.NE.0.0) THEN
         WRITE (iprint, 99115) hmaximum
99115        FORMAT (/,' MAXIMUM SMOOTHING LENGTH = ', ES12.3, //)
      ENDIF
!
!--Write name of file used
!
      WRITE (iprint, 99025) file1
99025 FORMAT (//,' Name of input file : ', A7, //)

      END SUBROUTINE header
