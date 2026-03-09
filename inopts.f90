      SUBROUTINE inopts
!************************************************************
!                                                           *
!  This subroutine defines all options desired for the run  *
!                                                           *
!************************************************************

      use idims

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
      use stepopt
      use init
      use numpa
      use ptdump
      use xtorq

      implicit none
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine inopts')
!
!--Open input file
!
      OPEN (iterm, FILE='inspho')
!
!--Determine options for insph
!
99002 FORMAT (A20)
99003 FORMAT (A7)
!
!--Read options
!
      READ (iterm, *) igrp
      READ (iterm, *) ifsvi,alpha,beta
      READ (iterm, *) ifcor
      READ (iterm, *) ichoc
      READ (iterm, *) iener
      READ (iterm, *) damp
      READ (iterm, *) iexf
      READ (iterm, *) iexpan
      READ (iterm, *) nstep
!---- Setting for 1st call of "inform". Any nprout /=0 works.
      nprout = MAX(nstep/10,1)

      READ (iterm, *) iptoutnum
      READ (iterm, *) tol, tolptm, tolh
      READ (iterm, *) ipos
      READ (iterm, *) tmax
      READ (iterm, *) tstop
      READ (iterm, *) dtmax
      READ (iterm, *) dtini
      omeg0 = 0.
      IF (ifcor.NE.0) THEN
         READ (iterm, *) omeg0
      ENDIF
      vexpan = 0.
      IF (iexpan.NE.0) THEN
         READ (iterm, *) vexpan
      ENDIF
!
!--Check for consistency
!
      CALL chekopt

      IF (idebug(1:7).EQ.'inopts') THEN
         WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc, iener, &
              ibound, damp, varsta
99004    FORMAT (1X, 7(I2,1X), 2(E12.5,1X), 1X, A7)
         WRITE (iprint, 99005) file1, ipos, nstep
99005    FORMAT (1X, A7, 1X, I4, 1X, I4)
      ENDIF

      CLOSE (iterm)

      END SUBROUTINE inopts
