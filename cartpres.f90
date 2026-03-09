      SUBROUTINE cartpres
!************************************************************
!                                                           *
!     This subroutine gives particles masses in a cartesian *
!              coordinate distribution                      *
!                                                           *
!************************************************************

      use idims

      use constants
      use part
      use rbnd
      use diskbd
      use debug
      use logun
      use maspres

      implicit none

      INTEGER(I4B) :: i, icount, ireg, j
      REAL(DP) :: rxmin(10), rxmax(10), rymin(10), rymax(10)
      REAL(DP) :: rzmin(10), rzmax(10), pres(10)
      REAL(DP) :: fractot, rcyl2, rmax2, rmind2, xmax5, ymax5, zmax5
      CHARACTER(len=1) :: iok
!
!--Allow for tracing flow
!
99004 FORMAT (A1)
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine cartpres ')
!
!--Set Variable Pressure Areas In Cartesian Coordinates
!
      rcyl2 = rcyl * rcyl
      rmind2 = rmind * rmind
      rmax2 = rmax * rmax

 40   WRITE (*, 88002)
88002 FORMAT (' Enter number of different regions (max 10)')
      READ (*,*) ireg
      IF (ireg.GT.10) GOTO 40

      icount = 1
      DO 120 i = 1, ireg
 60      WRITE (*, 88004) icount
88004    FORMAT ('   Enter xmin of region ', I2)
         READ (*,*) rxmin(icount)
         WRITE (*, 88006) icount
88006    FORMAT ('   Enter xmax of region ', I2)
         READ (*,*) rxmax(icount)
         WRITE (*, 88008) icount
88008    FORMAT ('   Enter ymin of region ', I2)
         READ (*,*) rymin(icount)
         WRITE (*, 88010) icount
88010    FORMAT ('   Enter ymax of region ', I2)
         READ (*,*) rymax(icount)
         WRITE (*, 88012) icount
88012    FORMAT ('   Enter zmin of region ', I2)
         READ (*,*) rzmin(icount)
         WRITE (*, 88014) icount
88014    FORMAT ('   Enter zmax of region ', I2)
         READ (*,*) rzmax(icount)

 80      WRITE (*, 88016) icount
88016    FORMAT ('   Enter relative pressure of region ', I2, &
                ' (0.0 to 1.0)')
         READ (*,*) pres(icount)
         IF ((pres(icount).LT.0.).OR.(pres(icount).GT.1.)) GOTO 80

 100     WRITE (*, 88018) icount
88018    FORMAT (' Is region ', I2,' correct (y/n)? ')
         READ (*, 99004) iok
         IF ((iok.NE.'y').AND.(iok.NE.'n')) GOTO 100
         IF (iok.NE.'y') GOTO 60

         icount = icount + 1
 120  CONTINUE

      xmax5 = xmax * 0.5
      ymax5 = ymax * 0.5
      zmax5 = zmax * 0.5
      fractot = 0.

      DO 240 i = 1, npart
         DO 220 j = 1, ireg
            IF ((x(i).GE.rxmin(j)) .AND. (x(i).LE.rxmax(j)) .AND. &
                (y(i).GE.rymin(j)) .AND. (y(i).LE.rymax(j)) .AND. &
                (z(i).GE.rzmin(j)) .AND. (z(i).LE.rzmax(j))) THEN
               disfrac(i) = pres(j)
               fractot = fractot + disfrac(i)
               GOTO 240
            ENDIF
 220     CONTINUE

 240  CONTINUE

      fractot = fractot/DBLE(npart)
      WRITE (*,*) fractot
      DO 300 i = 1, npart
         disfrac(i) = disfrac(i)/fractot
 300  CONTINUE

      END SUBROUTINE cartpres
