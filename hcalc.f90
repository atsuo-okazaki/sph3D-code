      SUBROUTINE hcalc
!************************************************************
!                                                           *
!  This subroutine computes the derivative of the smoothing *
!     length.                                               *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use force
      use densi
      use logun
      use debug
      use rbnd
      use diskbd
      use part
      use carac
      use typef
      use ghost
      use outneigh
      use current
      use neighbor
      use curlist
      use phase

      implicit none

      INTEGER(I4B) :: i, iagain, icount, idumg, idummy, ikount, &
               ioutinf1, ioutinf2, ioutsup1, ioutsup2, icutsup1, icutsup2, &
               ipart, neiinf, neimean, neirange, neisup, ntot, numneigh
      REAL(DP) :: third, acc
!
!--Allow for tracing flow
!
      IF (myrank.eq.0) then
         IF (itrace.EQ.'all') WRITE (iprint, 99001)
      ENDIF
99001 FORMAT (' entry subroutine hcalc')

      third = 1./3.
      neimin = 30
      neimax = 70
!      neimin = 80
!      neimax = 120
      acc = 0.3
      idumg = 0
      idummy = 0

      neimean = (neimax + neimin) / 2
      neirange = (neimax - neimin) / 4
      neisup = neimean + neirange
      neiinf = neimean - neirange

      ikount = 0

      IF (myrank.eq.0) WRITE (*,*) neimax, neimean, neimin, &
          neisup, neiinf, neirange

    5 CONTINUE

      nlst = npart
      DO i = 1, npart
         iscurrent (i) = 1
         llist (i) = i
         hmax (i) = 5.0 * h (i)
      ENDDO

      SELECT CASE (ibound)
      CASE (1)
         CALL ghostp1 (npart, x, y, z, vx, vy, vz, u, h)
      CASE (2)
         CALL ghostp2 (npart, x, y, z, vx, vy, vz, u, h)
      CASE (3, 8, 96, 97)
         CALL ghostp3 (npart, x, y, z, vx, vy, vz, u, h)
      CASE default
!!         IF (ibound.EQ.0.OR.ibound.GE.90) nghost = 0
         nghost = 0
      END SELECT

      ntot = npart + nghost

      IF (igrape.EQ.0) THEN
         IF (myrank.eq.0) WRITE ( * , * ) ' Making tree'
         CALL insulate (1, ntot, npart, x, y, z, pmass, h)
      ENDIF

      icount = 0

   10 CONTINUE

      icount = icount + 1
      ikount = ikount + 1
      ioutinf1 = 0
      ioutsup1 = 0
      ioutinf2 = 0
      ioutsup2 = 0
      iagain = 0
      IF (myrank.eq.0) WRITE ( * , * ) ' Calculating neighbour changes'
!
!--Get neighbours
!
      IF (igrape.EQ.0) THEN
         CALL insulate (5, ntot, npart, x, y, z, pmass, h)
      ELSEIF (igrape.EQ.1) THEN
         CALL insulate (4, ntot, npart, x, y, z, pmass, h)
      ENDIF

      DO ipart = 1, npart

         IF (iphase (ipart) .GE.1) CYCLE

         numneigh = nneigh (ipart)
         IF (DBLE (numneigh) / DBLE (neimean) .LT.0.1) &
            numneigh = neimean / 10
         IF (numneigh.LT.neiinf) THEN
            IF (numneigh.LT.neimin) THEN
               iagain = iagain + 1
               ioutinf1 = ioutinf1 + 1
               h(ipart) = h(ipart)*(DBLE(neimean) &
                           /DBLE(numneigh+1))**third
            ELSE
               ioutinf2 = ioutinf2 + 1
               h(ipart) = h(ipart)*(DBLE(neimean) &
                          /DBLE(numneigh+1))**third
            ENDIF
         ELSEIF (numneigh.GT.neisup) THEN
            IF (numneigh.GT.neimax) THEN
               iagain = iagain + 1
               ioutsup1 = ioutsup1 + 1
               h(ipart) = h(ipart)*(DBLE(neimean) &
                          /DBLE(numneigh+1))**third
            ELSE
               ioutsup2 = ioutsup2 + 1
               h(ipart) = h(ipart)*(DBLE(neimean) &
                          /DBLE(numneigh+1))**third
            ENDIF
         ENDIF
      ENDDO

      IF (myrank.eq.0) then
         IF (ioutinf1.NE.0) WRITE ( * ,  * ) 'h too small ', &
            ioutinf1, ' times'
         IF (ioutsup1.NE.0) WRITE ( * ,  * ) 'h too big ', &
            ioutsup1, ' times'
         IF (ioutinf2.NE.0) WRITE ( * , * ) 'h near lower limit ', &
            ioutinf2, ' times'
         IF (ioutsup2.NE.0) WRITE ( * , * ) 'h near upper limit ', &
            ioutsup2, ' times'
      ENDIF

      IF (ikount.GT.10) GOTO 15
!!      IF (iagain.GT.0 .AND. icount.GT.2) GOTO 5
!!      IF (iagain.GT.0) GOTO 10
      IF (iagain.GT.0) GOTO 5

   15 DO i = 1, npart
         iscurrent (i) = 0
      ENDDO

      END SUBROUTINE hcalc
