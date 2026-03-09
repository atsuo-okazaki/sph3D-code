      SUBROUTINE unit
!************************************************************
!                                                           *
!  This routine computes the transformation between the     *
!     physical units (cgs) and the units used in the code.  *
!                                                           *
!************************************************************

      use idims

      use constants
      use units
      use kerne
      use logun
      use debug
      use polyk2
      use vargam
      use physeos

      implicit none

      REAL(DP) :: uerg
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine unit')
!
!--Specify mass unit (g)
!
!c      fnbtot = 1.
!      fnbtot = 1000.
!      fnbtot = 0.1
!      fnbtot = 0.2
!c      umass = DBLE(fnbtot)*solarm
!
!--Specify distance unit (cm)
!
!c      udist = 1.e16
!      udist = 1.e15
!      udist = 0.1*pc
!      udist = pc
!      udist = 1.e14
!      udist = 1.e13
!      udist = 1000.0*pc
!
!--------------------
!--VARIABLE GAMMA EOS
!--------------------
!
!--Critical density for changing gamma from 1 to 1.4 for variable eq. of state
!
!      rhocrit = 1.e-14
!      rhocrit = 2.e-13
      rhocrit = 1.e-13
!      rhocrit = 1.e-12
!      rhocrit = 1.e-16
      gam = 1.4
!      gam = 5./3.
!
!--Critical density for changing gamma from 1.4 to 1.1 for variable e.o.s.
!     i.e. at 2000K assuming trans to 1.4 at 10K
!
!-- gam = 1.4
!
       rhocrit2 = rhocrit * (200.**2.5)
!
!-- gam = 5/3
!
!      rhocrit2 = rhocrit * (200.**1.5)
!
!      rhocrit2 = 1.0e-11
!      rhocrit2 = 1.0e-12
!      gamdh = 1.10
      gamdh = 1.15
!      gamdh = 1.05
!      gamdh = 1.0
!
!--Critical density for changing gamma from 1.1 to 5/3 for variable e.o.s.
!
      rhocrit3 = 1.0e-3
!      rhocrit3 = 1.0e-10
!      rhocrit3 = 1.0e-11
!      rhocrit3 = 1.0e-12
      gamah = 5./3.
!
!--****** For turning off 2nd collapse phase! ******
!
      rhocrit2 = rhocrit3
!
!
!
!-------------------------------------------
!--PHYSICAL EOS (Bodenheimer, Bate, Burkert)
!-------------------------------------------
!
!--Changing gamma from 1 to 5/3 for physical eq. of state
!
      gamphys1 = 5./3.
      rhophys1 = 2.0e-13
      rhoref1 = rhophys1
      rhochange1 = rhophys1*27.
!
!--Changing gamma from 5/3 to 1.4 for physical e.o.s.
!
      gamphys2 = 1.4
      rhoref2 = rhophys1/9.
      rhochange2 = 7.5e-11
!
!--Changing gamma from 1.4 to 2 for physical e.o.s.
!
      gamphys3 = 2.0
      rhoref3 = 2.9e-12
!
!
!
!
!
!--Transformation factor for :
!
!  a) density
!
      udens = DBLE(umass)/DBLE(udist)**3

      rhocrit = rhocrit / udens
      rhocrit2 = rhocrit2 / udens
      rhocrit3 = rhocrit3 / udens

      rhophys1 = rhophys1 / udens
      rhoref1 = rhoref1 / udens
      rhoref2 = rhoref2 / udens
      rhoref3 = rhoref3 / udens
      rhochange1 = rhochange1 / udens
      rhochange2 = rhochange2 / udens
!
!  b) time
!
!c      utime = DSQRT(DBLE(udist)**3/(DBLE(gg)*DBLE(umass)))
!
!  c) ergs
!
      uerg = DBLE(umass)*DBLE(udist)**2/DBLE(utime)**2
!
!  c) ergs per gram
!
      uergg = DBLE(udist)**2/DBLE(utime)**2
!
!  d) ergs per cc
!
      uergcc = DBLE(umass)/(DBLE(udist)*DBLE(utime)**2)

      IF (idebug(1:4).EQ.'unit') THEN
         WRITE (iprint, 99002) umass, udist, udens, utime, uergg, uergcc
99002    FORMAT (1X, 5(1PE12.5,1X), /, 1X, 2(1PE12.5,1X))
      ENDIF

      END SUBROUTINE unit
