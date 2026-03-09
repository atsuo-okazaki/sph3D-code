      SUBROUTINE chekopt
!************************************************************
!                                                           *
!   This subroutine checks for possible incompatibilities   *
!      between the different options selected; it also      *
!      warns the user of possible problems related to the   *
!      present limitations of the code                      *
!                                                           *
!************************************************************

      use idims

      use constants
      use typef
      use units
      use dissi
      use rotat
      use tming
      use integ
      use varet
      use cgas
      use recor
      use rbnd
      use diskbd
      use expan
      use kerne
      use files
      use actio
      use logun
      use debug

      implicit none

      CHARACTER(len=7) :: where='chekopt'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine chekopt')

      IF (iener.EQ.1 .AND. damp.NE.0.) CALL error(where, 2)

      IF (encal.EQ.'i' .AND. (iener.EQ.1 .OR. ichoc.EQ.1)) &
                                          CALL error(where,3)

      IF (encal.EQ.'a' .AND. (iener.EQ.0 .OR. ichoc.EQ.0)) &
                                          CALL error(where,4)

      IF ((encal.EQ.'p' .OR. encal.EQ.'v' .OR. encal.EQ.'x') .AND. &
            (iener.EQ.1 .OR. ichoc.EQ.1)) CALL error(where,5)

      IF (encal.EQ.'c' .AND. (iener.EQ.0 .OR. ichoc.EQ.0)) &
                                          CALL error(where,4)

      END SUBROUTINE chekopt
