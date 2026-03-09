      SUBROUTINE options
!************************************************************
!                                                           *
!  This subroutine defines all options desired for the run  *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use constants
      use typef
      use units
      use dissi
      use rotat
      use tming
      use integ
      use varet
      use recor
      use rbnd
      use diskbd
      use expan
      use kerne
      use files
      use actio
      use logun
      use debug
      use cgas
      use stepopt
      use init
      use pres
      use xforce
      use numpa
      use ptmass
      use binary
      use ptdump
      use crpart
      use ptbin
      use useles
      use new
      use polyk2
      use sphcom
      use xtorq
      use misali
      use maslos
      use winds
      use split

      implicit none

      INTEGER(I4B) :: i, ibin, icarry, j, namelength, number
      REAL(DP) :: amaxinfl, amininfl, cos1, cos2, cos3, &
               fcutoff, partm0, rd1, rd2, sin1, sin2, sin3, xlog2
!
!--Allow for tracing flow
!
      IF (myrank.eq.0) then
         IF (itrace == 'all') WRITE (iprint, 99001)
      ENDIF
99001 FORMAT (' entry subroutine options')
!
!--Open input file
!
      OPEN (iterm, FILE = inname)
!
!--Determine options for evolution run
!
!--Read name of run
!
      READ (iterm, 99002) namerun
      IF (myrank.eq.0) WRITE (namenextrun, 99002) namerun
      namelength = LEN (namenextrun)
99002 FORMAT (A20)
      DO i = namelength, 1, - 1
         IF (namenextrun(i:i).GE.'0' .AND. namenextrun(i:i).LE.'9') THEN
            READ (namenextrun(i:i), 88001) number
88001       FORMAT(I1)
            icarry = 0
            number = number + 1
            IF (number.GT.9) THEN
               icarry = 1
               number = 0
            ENDIF
            IF (myrank.eq.0) WRITE (namenextrun(i:i), 88001) number
            j = i - 1
            IF ((j.GT.0) .AND. (icarry == 1) .AND. (namenextrun(j:j) &
            .GE.'0' .AND. namenextrun(j:j).LE.'9') ) THEN
               READ (namenextrun(j:j), 88001) number
               number = number + icarry
               IF (myrank.eq.0) WRITE (namenextrun(j:j), 88001) number
            ENDIF
            GOTO 100
         ELSEIF (namenextrun (i:i) .NE.' ') THEN
            GOTO 100
         ENDIF
      ENDDO
!
!--Open output file
!
  100 IF (iprint.NE.6) OPEN (iprint, FILE = namerun)

      CALL labrun
!
!--Read name of file containing physical input
!
      READ (iterm, 99003) file1
99003 FORMAT (A7)
      READ (iterm, 99003) varsta
!
!--Read options
!
      READ (iterm, 99006) encal
99006 FORMAT (A1)
      READ (iterm,*) initialptm
      READ (iterm, 99006) iaccevol
      IF (iaccevol == 'v'.OR.iaccevol == 's') READ (iterm,*) accfac
      READ (iterm,*) iptmass
      READ (iterm,*) igrp
      READ (iterm,*) igphi
      READ (iterm,*) ifsvi, alpha, beta
      READ (iterm,*) ifcor
      READ (iterm,*) ichoc
      READ (iterm,*) iener
      READ (iterm,*) damp
      READ (iterm,*) ibound
      READ (iterm,*) iexf
      READ (iterm,*) iexpan
      READ (iterm,*) nstep
      IF (nstep.LT.1) nstep = 1
      nprout = MAX (nstep/10, 1)

      READ (iterm,*) iptoutnum
      xlog2 = 0.30103 + 0.00001
      ibin = INT(LOG10(DBLE(iptoutnum))/xlog2) + 1
      iptoutnum = 2**ibin

      READ (iterm,*) umass, udist, utime
      CALL unit

      READ (iterm,*) tol, tolptm, tolh
      READ (iterm,*) ipos
      READ (iterm,*) tmax
      READ (iterm,*) tstop
      READ (iterm,*) dtmax
      READ (iterm,*) dtini

      omeg0 = 0.
      IF (ifcor.NE.0) THEN
         READ (iterm,*) omeg0
         omeg0 = omeg0*utime
      ENDIF

      vexpan = 0.
      IF (iexpan.NE.0) THEN
         READ (iterm,*) vexpan
         vexpan = vexpan*utime/udist
      ENDIF

      pext = 0.
      hmaximum = 0.
!!      extu = 0.
!!      extdens = 0.
      IF (ibound == 7) THEN
         READ (iterm,*) hmaximum
         READ (iterm,*) pext
!!         READ (iterm,*) exttemp
!!         extu = 3.0/2.0*exttemp*Rg/gmw/uergg
!!         READ (iterm,*) extdens
!!         extdens = extdens/udens
      ENDIF

      IF (ibound == 8) THEN
         READ (iterm,*) deadbound
         READ (iterm,*) fractan, fracradial, nstop, nfastd
      ENDIF
      IF (ibound.GE.90) THEN
         READ (iterm,*) hmaximum
         READ (iterm,*) deadbound
         READ (iterm,*) fractan, fracradial, nshell, rshell
      ENDIF

      xmass = 0.
      IF (iexf == 5.OR.iexf == 6) THEN
         READ (iterm,*) xmass
      ENDIF

      IF (iptmass.NE.0.OR.initialptm.NE.0) THEN
         READ (iterm,*) hacc
         READ (iterm,*) haccall
      ENDIF

      IF (iptmass.NE.0) THEN
         READ (iterm,*) radcrit
         READ (iterm,*) ptmcrit
      ENDIF
!
!--Read boundaries
!
      rmax = 0.
      rcyl = 0.
      rmind = 0.
      IF (ibound == 1 .OR. ibound == 3 .OR. ibound == 8 &
      .OR.ibound.GE.90) THEN
         READ (iterm,*) igeom, rmind, rmax, xmin, xmax, ymin, ymax,   &
            zmin, zmax
      ELSEIF (ibound == 2) THEN
         READ (iterm,*) igeom, rmind, rcyl, xmin, xmax, ymin, ymax,   &
            zmin, zmax
      ELSE
         READ (iterm,*) igeom, rd1, rd2, xmin, xmax, ymin, ymax, zmin,&
            zmax
         IF (rd1.GE.0.) THEN
            rcyl = rd2
            rmind = rd1
         ELSE
            rmax = rd2
         ENDIF
      ENDIF

      READ (iterm,*) rptmas(1), rptmas(2)
      READ (iterm,*) isphcom
      READ (iterm,*) rangle(1), rangle(2)
!---- The third Euler angle should be zero.
      rangle (3) = 0.0
      SELECT CASE (ibound)
      CASE (92, 93)
         READ (iterm,*) emdot0, partm
      CASE (94, 96)
         READ (iterm,*) nshell1(1), nshell1(2), rshell1, rshell2,   &
            vwind1, vwind2, vrot1, vrot2, emdot0, emdotratio, &
            partmass(1), partmass(2), sinj0(1), sinj0(2), &
            RK21, RK22, RK23
         DO i = 1, 2
            sinj(i) = sinj0(i)
         ENDDO
!----    thermal2 and therma3 are defined here
!        (thermal1 is defined later in "rdump")
         IF (encal == 'i') THEN
            thermal2 = RK22
            thermal3 = RK23
         ELSE IF (encal == 'c') THEN
!----    thermal1 gives the floor temperature for simulations
!        with radiative cooling (in radcool.f)
            thermal2 = RK22
            thermal3 = RK23
         ENDIF
      CASE (95, 97)
         READ (iterm,*) nshell1(1), nshell1(2), rshell1, rshell2,   &
            vwind1, vwind2, vrot1, vrot2, emdot0, emdotratio, &
            partmass(1), partmass(2), sinj0(1), sinj0(2), &
            nshell1(3), partmass(3), sinj0(3), &
            RK21, RK22, RK23
         DO i = 1, 3
            sinj(i) = sinj0(i)
         ENDDO
!----    thermal2 and therma3 are defined here
!        (thermal1 is defined later in "rdump")
         IF (encal == 'i') THEN
            thermal2 = RK22
            thermal3 = RK23
         ELSE IF (encal == 'c') THEN
!----    thermal1 gives the floor temperature for simulations
!        with radiative cooling (in radcool.f)
            thermal2 = RK22
            thermal3 = RK23
         ENDIF
      CASE (99)
         READ (iterm,*) emdot0, sinj0(1), partm
         sinj(1) = sinj0(1)
         READ (iterm, 99007) nameinfl
99007    FORMAT (A20)
         OPEN (UNIT = 19, FILE = nameinfl, FORM = 'formatted')
         READ (19,*) partm0, dphasei, nphasei, &
            (phasei(i), i = 1, nphasei), &
            (eninfl(i), i = 1, nphasei), &
            (xinfl(i), i = 1, nphasei), &
            (yinfl(i), i = 1, nphasei), &
            (zinfl(i), i = 1, nphasei), &
            (sxinfl(i), i = 1, nphasei), &
            (syinfl(i), i = 1, nphasei), &
            (szinfl(i), i = 1, nphasei), &
            (vxinfl(i), i = 1, nphasei), &
            (vyinfl(i), i = 1, nphasei), &
            (vzinfl(i), i = 1, nphasei), &
            (svxinfl(i), i = 1, nphasei), &
            (svyinfl(i), i = 1, nphasei), &
            (svzinfl(i), i = 1, nphasei), &
            (hinfl(i), i = 1, nphasei), &
            (shinfl(i), i = 1, nphasei)
         IF (myrank.eq.0) CLOSE (19)
         amaxinfl = 0.0
         DO i = 1, nphasei
            IF (eninfl(i).GT.amaxinfl) amaxinfl = eninfl(i)
         ENDDO
!----    Set eninfl(i) to be zero if it is less than
!        fcutoff*Max(eninfl). This is for avoiding injection
!        of particles with large h(i), which often causes
!        problems in making neighbor lists.
!!         fcutoff = 0.05
!!         fcutoff = 0.03
!!         fcutoff = 0.01
         fcutoff = 0.0
         amininfl = INT (amaxinfl*fcutoff)
         DO i = 1, nphasei
            IF (eninfl(i).LT.amininfl) eninfl(i) = 0.0
         ENDDO
         IF (myrank.eq.0) WRITE (*,*) (eninfl(i),i=1,nphasei)
      CASE default
      END SELECT

      IF (igeom == 8 .OR. igeom == 9) THEN
         READ (iterm,*) azimuth1(1), azimuth2(1), vangle1(1), vangle2(1)
         READ (iterm,*) azimuth1(2), azimuth2(2), vangle1(2), vangle2(2)
      ENDIF

      SELECT CASE (iexf)
      CASE (3)
         READ (iterm,*) xeps, xbeta
!!      ELSE IF (iexf == 7) THEN
!!         READ (iterm,*) xantgrav(1),xantgrav(2)
      CASE (7, 8)
         READ (iterm,*) xantgrav (1), xantgrav (2)
         READ (iterm,*) akappac (1), akappar (1), akappac (2), &
                        akappar (2)
         READ (iterm,*) vbeta (1), vbeta (2)
      CASE default
      END SELECT

      READ (iterm,*) isplit
      numsplit = 0

!---- Before 31 March 2003, rotation angles about x- and y-axies
!     had been used.
!--   arot: matrix for rotation about x- and y-axes.
!           Rotation about x-axis first, then about y-axis.
!!      sinx = SIN (ranglex)
!!      cosx = COS (ranglex)
!!      siny = SIN (rangley)
!!      cosy = COS (rangley)
!!      arot(1,1) = cosy
!!      arot(2,1) = 0.0
!!      arot(3,1) = -siny
!!      arot(1,2) = sinx*siny
!!      arot(2,2) = cosx
!!      arot(3,2) = sinx*cosy
!!      arot(1,3) = cosx*siny
!!      arot(2,3) = -sinx
!!      arot(3,3) = cosx*cosy
!---- After 31 March 2003, rotation using Euler angles
!     have been used. (A. Okazaki, 31 Mar 2003)
      sin1 = SIN (rangle (1))
      cos1 = COS (rangle (1))
      sin2 = SIN (rangle (2))
      cos2 = COS (rangle (2))
      sin3 = SIN (rangle (3))
      cos3 = COS (rangle (3))
      arot (1, 1) = cos3*cos1 - cos2*sin1*sin3
      arot (2, 1) = cos3*sin1 + cos2*cos1*sin3
      arot (3, 1) = sin2*sin3
      arot (1, 2) = - sin3*cos1 - cos2*sin1*cos3
      arot (2, 2) = - sin3*sin1 + cos2*cos1*cos3
      arot (3, 2) = sin2*cos3
      arot (1, 3) = sin2*sin1
      arot (2, 3) = - sin2*cos1
      arot (3, 3) = cos2

!
!--Check for consistency
!
      CALL chekopt

      if (myrank.eq.0) then
         WRITE (iprint,*) 'tstop=', tstop
         IF (idebug (1:7) == 'options') THEN
            WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, &
               iener, ibound, damp, varsta
99004       FORMAT(1X, 7(I2,1X), E12.5, 1X, A7)
            WRITE (iprint, 99005) file1, ipos, nstep
99005       FORMAT(1X, A7, 1X, I4, 1X, I4)

            SELECT CASE (ibound)
            CASE (1)
               WRITE (iprint,*) rmax, xmin, xmax, ymin, &
                                ymax, zmin, zmax
            CASE (3)
               WRITE (iprint,*) rmind, rmax, xmin, xmax, &
                                ymin, ymax, zmin, zmax
            CASE (2)
               WRITE (iprint,*) rmind, rcyl, xmin, xmax, &
                                ymin, ymax, zmin, zmax
            CASE default
            END SELECT
         ENDIF
      endif

      IF (myrank.eq.0) CLOSE(iterm)

      END SUBROUTINE options
