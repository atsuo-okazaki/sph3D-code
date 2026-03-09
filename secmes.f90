      SUBROUTINE secmes
!************************************************************
!                                                           *
!  This subroutine handles the messages received from the   *
!     operator.                                             *
!                                                           *
!************************************************************

      use idims

      use tming
      use stop
      use logun
      use stepopt
      use init
      use timei
      use part
      use phase
      use gtime
      use secret

      implicit none

      INTEGER(I4B) :: i, iwait, nmes
      REAL(DP) :: attime
      CHARACTER(len=20) :: string(10)
      CHARACTER(len=7) :: where='secret'
      CHARACTER(len=5) :: dummy5
      CHARACTER(len=3) :: dummy3
!
!--Open message file first
!
      iwait = 0
      OPEN (idisk2, FILE='secret', FORM='formatted')
!
!--Read message file (one mes. per line , max=10)
!
      nmes = 0
      DO i = 1, 10
         READ (idisk2, 99000, END=200) string(i)
99000    FORMAT (A20)
         nmes = i
      END DO
!
!--Read messages one ofter another
!
 200  DO 300 i = 1, nmes
!
!--Identify each message
!
!  0) wait before processing messages
!
         IF (string(i)(1:2).EQ.'at') THEN
            READ (string(i), *) dummy3, attime
            IF (gt.LT.attime) THEN
               iwait = 1
               WRITE (*, 99001) i
99001          FORMAT (' message ', I2, ' read. waiting to process.')
               GOTO 400
            ELSE
               GOTO 300
            ENDIF
         ENDIF
!
!  4) change minimum time for keeping GRAPE
!
         IF (string(i)(1:5).EQ.'tkeep') THEN
            READ (string(i), *) dummy5, tkeep
            IF (tkeep.LT.3.0) tkeep = 3.0
            WRITE (*, 99010) i, tkeep
99010    FORMAT (' message ', I2, ' read. tkeep changed to : ',1PE12.5)
            GOTO 300
         ENDIF
!
!--Unexpected messages
!
         WRITE (*, 99010) string(i)
99090    FORMAT (' unexpected message received : ', /, A20)

 300  CONTINUE
!
!--Erase message file
!
 400  CONTINUE
!      IF (iwait.EQ.0) THEN
!         CLOSE (idisk2, STATUS='delete')
!      ELSE
         CLOSE (idisk2)
!      ENDIF

      END SUBROUTINE secmes
