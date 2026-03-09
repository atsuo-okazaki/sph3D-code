      SUBROUTINE reduce
!************************************************************
!                                                           *
!  This subroutine compresses binary output files by a      *
!     certain factor and writes the resulting output on a   *
!     single file                                           *
!                                                           *
!************************************************************

      use idims

      use units
      use part
      use densi
      use typef
      use carac
      use cgas
      use kerne
      use gtime
      use bodys
      use ener
      use fracg
      use phase
      use ptmass
      use binary
      use torq
      use timei
      use stepopt

      implicit none

      INTEGER(I4B) :: i, j, k, image, imo, ireduct, nfile
      REAL(DP) :: rhozero, RK2
      CHARACTER(len=20) :: ifile(10), ofile

99001 FORMAT (A20)

      PRINT *, 'number of files to be reduced (max. of 10) ?'
      READ *, nfile

      PRINT *, 'name of files to be reduced ?'

      DO i = 1, nfile
         READ (*, 99001) ifile(i)
      END DO

      PRINT *, 'name of reduced file ?'
      READ (*, 99001) ofile
      OPEN (UNIT = 7, FILE = ofile, FORM = 'unformatted')

      PRINT *, 'compression factor'
      READ *, ireduct

      image = 0
      imo = 0

      DO k = 1, nfile

         OPEN (UNIT = 5, FILE = ifile(k), FORM = 'unformatted')
         PRINT *, 'reading file ', ifile(k)

         DO j = 1, 9999
            image = image + 1
            PRINT *, 'reading image number ', image

            READ (5, END=100) udist, umass, utime, &
                 npart, n1, n2, gt, gamma, rhozero, RK2, &
                 (h(i), i=1, npart), escap, tkin, tgrav, tterm, &
                 (x(i), i=1, npart), (y(i), i=1, npart), &
                 (z(i), i=1, npart), (vx(i), i=1, npart), &
                 (vy(i), i=1, npart), (vz(i), i=1, npart), &
                 (u(i), i=1, npart), (pmass(i), i=1, npart), &
                 (rho(i), i=1, npart), (dgrav(i), i=1, npart), &
                 dtmax, (isteps(i), i=1, npart), &
                 (iphase(i), i=1, npart), &
                 nptmass, (listpm(i), i=1, nptmass), &
                 (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass), &
                 (spinz(i),i=1,nptmass), &
                 (angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass), &
                 (angaddz(i),i=1,nptmass), &
                 anglostx, anglosty, anglostz, &
                 nreassign, naccrete, nkill, specang, ptmassin, &
                 (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass), &
                 (spinadz(i),i=1,nptmass), &
                 (torqt(i), i=1, npart), (torqg(i), i=1, npart), &
                 (torqp(i), i=1, npart),(torqv(i), i=1, npart), &
                 (torqc(i), i=1, npart)


            IF (MOD(image, ireduct).EQ.0 .OR. &
                                         (k.EQ.1 .AND. j.EQ.1)) THEN

               PRINT *, 'writing image just read on output file'
               imo = imo + 1

               WRITE (7) udist, umass, utime, &
                    npart, n1, n2, gt, gamma, rhozero, RK2, &
                    (h(i), i=1, npart), escap, tkin, tgrav, tterm, &
                    (x(i), i=1, npart), (y(i), i=1, npart), &
                    (z(i), i=1, npart), (vx(i), i=1, npart), &
                    (vy(i), i=1, npart), (vz(i), i=1, npart), &
                    (u(i), i=1, npart), (pmass(i), i=1, npart), &
                    (rho(i), i=1, npart), (dgrav(i), i=1, npart), &
                    dtmax, (isteps(i), i=1, npart), &
                (iphase(i), i=1, npart), &
                nptmass, (listpm(i), i=1, nptmass), &
                (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass), &
                (spinz(i),i=1,nptmass), &
                (angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass), &
                (angaddz(i),i=1,nptmass), &
                anglostx, anglosty, anglostz, &
                nreassign, naccrete, nkill, specang, ptmassin, &
                (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass), &
                (spinadz(i),i=1,nptmass), &
                (torqt(i), i=1, npart), (torqg(i), i=1, npart), &
                (torqp(i), i=1, npart),(torqv(i), i=1, npart), &
                (torqc(i), i=1, npart)

            ENDIF

         END DO
 100     CLOSE (5)
      END DO

      PRINT *, 'writing last image on output file'
      imo = imo + 1

      WRITE (7) udist, umass, utime, &
           npart, n1, n2, gt, gamma, rhozero, RK2, &
           (h(i), i=1, npart), escap, tkin, tgrav, tterm, &
           (x(i), i=1, npart), (y(i), i=1, npart), &
           (z(i), i=1, npart), (vx(i), i=1, npart), &
           (vy(i), i=1, npart), (vz(i), i=1, npart), &
           (u(i), i=1, npart), (pmass(i), i=1, npart), &
           (rho(i), i=1, npart), (dgrav(i), i=1, npart), &
           dtmax, (isteps(i), i=1, npart), &
           (iphase(i), i=1, npart), &
           nptmass, (listpm(i), i=1, nptmass), &
           (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass), &
           (spinz(i),i=1,nptmass), &
           (angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass), &
           (angaddz(i),i=1,nptmass), &
           anglostx, anglosty, anglostz, &
           nreassign, naccrete, nkill, specang, ptmassin, &
           (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass), &
           (spinadz(i),i=1,nptmass), &
           (torqt(i), i=1, npart), (torqg(i), i=1, npart), &
           (torqp(i), i=1, npart),(torqv(i), i=1, npart), &
           (torqc(i), i=1, npart)


      CLOSE (7)

      PRINT 99002, ofile, imo
99002 FORMAT ('file ', A10, 'has been created and contains ', I3, &
              ' images')

      END SUBROUTINE reduce
