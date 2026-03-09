      SUBROUTINE labrun
!************************************************************
!                                                           *
!  This routine writes the title page                       *
!                                                           *
!************************************************************

      use idims

      use logun
      use actio
      use infor
      use zzhp

      implicit none

      INTEGER(I4B) :: iday, imon, iyear, ih, im, is
      REAL(DP) :: heure
!
!--Get time and date
!
      CALL getdat(iday, imon, iyear)
      CALL getime(ih, im, is, heure)
!
!--Write title page
!
      WRITE (iprint, 99001)
99001 FORMAT ('1', /, 1X, &
              '**********************************************', /, 1X, &
              '*                                            *', /, 1X, &
              '*      S P H - 3 D - TIMESTEPS - GRAPE       *', /, 1X, &
              '*              WITH POINT MASSES             *', /, 1X, &
              '*                                            *', /, 1X, &
              '**********************************************')
      WRITE (iprint, 99002) version
99002 FORMAT (/, 1X, 'Current version : ', A40)

      IF (igrape.EQ.0) THEN
         WRITE (iprint,*) '                  RUNNING WITH BINARY TREE'
      ELSEIF (igrape.EQ.1) THEN
         WRITE (iprint,*) '                  RUNNING WITH GRAPE'
      ENDIF

      WRITE (iprint, 99003)namerun, job, iday, imon, iyear, ih, im, is
      WRITE (*, 99003)namerun, job, iday, imon, iyear, ih, im, is
99003 FORMAT (//, ' SPH run ', A20, ' is running option : ', A11, /, &
              ' Started on : ', I2, '/', I2, '/', I4, ' at ', I2, &
              ' h. ', I2, ' min. ', I2, ' sec. ')

      END SUBROUTINE labrun
