      SUBROUTINE place(idisk1, ipos, irec, iflag)
!************************************************************
!                                                           *
!  This routine puts the pointer at beginning of dump ipos  *
!                                                           *
!************************************************************

      use idims, only : I4B
      implicit none

      INTEGER(I4B) :: idisk1, ipos, irec, iflag, i
      CHARACTER(len=7) :: where='place'

      irec = 0
      IF (ipos.NE.9999) THEN
!
!--Go to beginning of dump ipos
!
         DO i = 1, ipos - 1
            READ (idisk1, END=300)
            irec = i
         END DO
         irec = irec + 1
         RETURN
      ELSE
         DO i = 1, 10000
            READ (idisk1, END=200)
            irec = i
         END DO
      ENDIF
!
!--Go back one reccord
!
 200  BACKSPACE idisk1

      IF (iflag.EQ.1) BACKSPACE idisk1
      GOTO 400

 300  CALL error(where, 1)

 400  RETURN

      END SUBROUTINE place
