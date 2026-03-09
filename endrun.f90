      SUBROUTINE endrun
!************************************************************
!                                                           *
!  This subroutine ends the runstream.                      *
!                                                           *
!************************************************************

      use idims

      use logun
      use actio
      use debug
      use typef
      use ptmass
      use ptdump
      use zzhp

      implicit none

      INTEGER(I4B) :: ih, imin, is, ij, im, iy
      REAL(DP) :: heure, tused
!
!--Get time and date
!
      CALL getime(ih, imin, is, heure)
      CALL getdat(ij, im, iy)
      CALL getused(tused)
!
!--Write end page
!
      WRITE (iprint, 99001) namerun, ij, im, iy, ih, imin, is
99001 FORMAT (/////, 1X, 'SPH run ', A20, ' ended on : ', &
                 I2, '/', I2, '/', I4, '  at ', I2, ' h. ', I2, &
                 ' min. ', I2, ' sec.')
      WRITE (iprint, 99002) tused/60.
99002 FORMAT (1X, 'cpu time used for this run :', F10.3, ' min.')

!!      iprint = 6

      IF (job(1:9).EQ.'evolution') CLOSE (iprint)
      IF (job(1:9).EQ.'evolution' .AND. (nptmass.NE.0 .OR. &
                                                 iptmass.NE.0)) THEN
         CLOSE (iptprint)
         CLOSE (iaccpr)
      ENDIF
      IF (job(1:9).EQ.'evolution' .AND. &
                                (ibound.EQ.8 .OR. ibound.GE.90)) THEN
         CLOSE (ikillpr)
         CLOSE (ireasspr)
      ENDIF

      IF (job(1:9).EQ.'evolution') THEN
         WRITE (inotify,*) 'stopped'
         CLOSE (inotify)
      ENDIF

      CALL quit

      END SUBROUTINE endrun
