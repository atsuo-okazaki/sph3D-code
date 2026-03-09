      SUBROUTINE mesop
!************************************************************
!                                                           *
!  This subroutine handles the messages received from the   *
!     operator.                                             *
!                                                           *
!************************************************************

      use idims

      use constants
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
      use densi
      use sync
      use polyk2

      implicit none

      INTEGER(I4B) :: i, ifactor, imaxnew, iminnew, ipower, iwait, j, nmes
      REAL(DP) :: attime, contrast, diff, dtmaxsync, fffac, &
               freefalltime, ratio, rhomax, rhomaxsyncold, tempvar, &
               timemax
      CHARACTER(len=20) :: string(10)
      CHARACTER(len=7) :: where='mesop'
      CHARACTER(len=5) :: dummy5
      CHARACTER(len=3) :: dummy3
!
!--Open message file first
!
      iwait = 0
      OPEN (idisk2, FILE='message', FORM='formatted')
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
!  1) stop job
!
         IF (string(i)(1:4).EQ.'stop') THEN
            istop = 1
            WRITE (*, 99002) i
99002       FORMAT (' message ', I2, ' read. please wait for stop.')
            GOTO 300
         ENDIF
!
!  2) change dump frequency
!
         IF (string(i)(1:5).EQ.'nstep') THEN
            READ (string(i), 99005) nstep
99005       FORMAT (6X, I2)
            WRITE (*, 99006) i, nstep
99006       FORMAT (' message ', I2, ' read. nstep set to : ', I2)
            GOTO 300
         ENDIF
!
!  3) change syncronisation time
!
         IF (string(i)(1:5).EQ.'synct') THEN
            READ (string(i), 99007) ipower
99007       FORMAT (6X, I2)
            WRITE (*, 99008) i, ipower
99008       FORMAT (' message ', I2, ' read. synct changed by : ', I2)

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor

            IF (ipower.LT.0) THEN
               dtmax = dtmax/DBLE(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 300
               ENDIF
               dtmax = dtmax*DBLE(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ENDIF

            GOTO 300
         ENDIF
!
!  4) change minimum time for keeping GRAPE
!
         IF (string(i)(1:5).EQ.'tkeep') THEN
            READ (string(i), *) dummy5, tkeep
            WRITE (*, 99010) i, tkeep
99010    FORMAT (' message ', I2, ' read. tkeep changed to : ',1PE12.5)
            GOTO 300
         ENDIF
!
!--Unexpected messages
!
         WRITE (*, 99090) string(i)
99090    FORMAT (' unexpected message received : ', /, A20)

 300  CONTINUE
!
!--Erase message file
!
 400  IF (iwait.EQ.0) THEN
         CLOSE (idisk2, STATUS='delete')
      ELSE
         CLOSE (idisk2)
      ENDIF
!
!--Automatic change of syncronisation time using local free-fall time
!
      OPEN (idisk2, FILE='synctcontrol', FORM='formatted')
      READ (idisk2, 99000, END=500) string(i)

      IF (string(i)(1:4).EQ.'auto') THEN
         READ (string(i)(6:20),*) fffac
         WRITE (*, 99016) fffac
         WRITE (iprint, 99016) fffac
99016    FORMAT (' Synctcontrol : ', 1PE12.5)

         ipower = 0
         rhomax = 0.
         DO j = 1, npart
            IF (rho(j).GT.rhomax) rhomax = rho(j)
         END DO
         freefalltime = SQRT((3 * pi) / (32 * rhomax))
         timemax = fffac*freefalltime
         IF (dtmax.GT.timemax) THEN
            ipower = - (INT(LOG10(dtmax/timemax)/0.30103) + 1)
            WRITE (*, 99009) ipower, freefalltime, rhomax
            WRITE (iprint, 99009) ipower, fffac, freefalltime, rhomax
99009       FORMAT (' Synct autochange : ', I2,' fffac ',1PE12.5, &
              ' free-fall ',1PE12.5,' rhomax ',1PE12.5)

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor

            IF (ipower.LT.0) THEN
               dtmax = dtmax/DBLE(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 500
               ENDIF
               dtmax = dtmax*DBLE(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ENDIF
         ENDIF
      ENDIF
!
!--Automatic change of syncronisation time using density change
!
      IF (string(i)(1:4).EQ.'dens') THEN
         READ (string(i)(6:20),*) contrast, dtmaxsync

         rhomaxsyncold = rhomaxsync
         ipower = 0
         rhomaxsync = 0.
         DO j = 1, npart
            rhomaxsync = MAX(rhomaxsync, rho(j))
         END DO

         freefalltime = SQRT((3 * pi) / (32 * rhozero))

         ratio = LOG10(rhomaxsync/rhomaxsyncold)

         contrast = LOG10(contrast)
         ipower = - ABS(INT(ratio/contrast))
         IF (ABS(ratio/contrast).LT.0.5) ipower = 1

         WRITE (*, 99017) dtmaxsync, 10**contrast, &
              rhomaxsync/rhomaxsyncold, ipower
         WRITE (iprint, 99017) dtmaxsync, 10**contrast, &
              rhomaxsync/rhomaxsyncold, ipower
99017    FORMAT (' Dtmaxsync : ', 1PE12.5,' contrast ',1PE12.5, &
              ' ratio ',1PE12.5,' ipower ',I2)

         IF (ipower.GE.0 .AND. tstep/60.0.GT.720.0) THEN
            WRITE (iprint,99034) tstep/60.0
            WRITE (*,99034) tstep/60.0
99034       FORMAT('  Synct autochange since time ',1PE12.5)
            ipower = -1
         ENDIF
         IF (ipower.EQ.1) THEN
            tempvar = gt/(2.0*dtmax)
            diff = tempvar-INT(tempvar)
            IF (diff.GT.0.25 .AND. diff.LT.0.75) THEN
               WRITE (iprint,99035) diff
               WRITE (*,99035) diff
99035          FORMAT('  Synct autochange attempt, but sync ',1PE12.5)
               ipower = 0
            ENDIF
         ENDIF
         IF (ipower.NE.0) THEN

            ifactor = 2**ABS(ipower)
            imaxnew = imaxstep/ifactor
            iminnew = 2*ifactor

            IF (ipower.LT.0) THEN
               WRITE (*, 99021) ipower
               WRITE (iprint, 99021) ipower
99021          FORMAT ('  Synct autochange : ', I2)
               dtmax = dtmax/DBLE(ifactor)
               istepmin = MIN(istepmin, imaxnew)
               istepmax = MIN(istepmax, imaxnew)
               istepmin = istepmin*ifactor
               istepmax = istepmax*ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     isteps(j) = MIN(isteps(j), imaxnew)
                     isteps(j) = isteps(j)*ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ELSEIF (ipower.GT.0 .AND. &
                    dtmax*DBLE(ifactor).LE.dtmaxsync .AND. &
                    tstep/60.*DBLE(ifactor).LE.720.) THEN
               IF (istepmin/ifactor .LT. 2) THEN
                  CALL error(where, 1)
                  istop = 1
                  GOTO 500
               ENDIF
               WRITE (*, 99021) ipower
               WRITE (iprint, 99021) ipower

               dtmax = dtmax*DBLE(ifactor)
               istepmin = istepmin/ifactor
               istepmax = istepmax/ifactor

               DO j = 1, npart
                  IF (iphase(j).NE.-1) THEN
                     IF (isteps(j)/ifactor .LT. 2) CALL error(where, 2)
                     isteps(j) = isteps(j)/ifactor
                     it1(j) = isteps(j)/2
                  ENDIF
               END DO
            ELSE
               IF (dtmax*DBLE(ifactor).GT.dtmaxsync) THEN
                  WRITE (iprint,99040)
                  WRITE (*,99040)
99040       FORMAT('  Synct autochange attempt, but dtmaxsync reached')
               ELSEIF (tstep/60.*DBLE(ifactor).GT.720.) THEN
                  WRITE (iprint,99041)
                  WRITE (*,99041)
99041       FORMAT('  Synct autochange attempt, but tstepmax reached')
               ENDIF
            ENDIF
         ENDIF
      ENDIF

 500  CLOSE(idisk2)

      END SUBROUTINE mesop
