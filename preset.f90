      SUBROUTINE preset
!************************************************************
!                                                           *
!  This routine makes sure that everything not entered as   *
!     option is defined before integration begins.          *
!                                                           *
!************************************************************

!!      use mpi_mod
      use idims

      use constants
      use part
      use densi
      use force
      use kerne
      use logun
      use files
      use numpa
      use neighbor
      use carac
      use cgas
      use bodys
      use btree
      use stop
      use tming
      use recor
      use typef
      use polyk2
      use ptmass
      use nextmpt
      use phase
      use ptdump
      use actio
      use ptsoftx
      use soft
      use init
      use delay
      use winds
      use cooldata
      use unidis

      implicit none

      INTEGER(I4B) :: i, i1, ichkl, iptcur, nnotacc
      REAL(DP) :: ran1
      REAL(DP) :: hmin, pmassmin, totmas, xmax, xmin, ymax, ymin, &
               zmax, zmin
      CHARACTER(len=7) :: where='preset'
      CHARACTER(len=21) :: ptfile, accfile, killfile, &
               reassfile, antigfile, belongfile
!!      CHARACTER(len=21) :: captfile
      CHARACTER(len=21) :: coolfile
      CHARACTER(len=7) :: inject

!---- Specify the method for particle injection
      inject = 'random'
!!      inject = 'uniform'
!
!--Open file1
!
!!    OPEN (idisk1, FILE=file1, STATUS='unknown', FORM='unformatted', &
!!          RECL=imaxrec)
      OPEN (idisk1, FILE = file1, STATUS = 'unknown', &
            FORM = 'unformatted')

!--Find correct dump to start
!
      CALL place (idisk1, ipos, irec, 1)
!
!--Read dump
!
      CALL rdump (idisk1, ichkl)
      IF (ichkl == 1) CALL error (where, ichkl)

      IF (encal == 'c') THEN
!----    Read the cooling table (ntemp<=MAXDATA) from file
!         1st column: Temperature in Kelvin
!         2nd column: Lambda(T) in erg/cm^3/sec
         OPEN (19, FILE = 'cooltable.dat', STATUS = 'OLD', &
               FORM = 'formatted')
         i = 1
         DO WHILE (i <= MAXDATA)
            READ (19,*, END=20) tempcf(i), alambda(i)
            i = i+1
         ENDDO
   20    CLOSE (19)
         ntemp = i-1

!----    alphacf: slope of the cooling function
         DO i=1,ntemp-1
            i1 = i+1
            alphacf(i) = LOG(alambda(i1)/alambda(i)) &
                         /LOG(tempcf(i1)/tempcf(i))
         ENDDO
         alphacf(ntemp) = alphacf(ntemp-1)

!----    Rich Townsend (2009)'s Y function for a piecewise power-law
!        cooling function [eq.(A6) of Townsend (2009, ApJS, 181, 396)]
         yfunc(ntemp) = 0.0
         DO i=ntemp-1,1,- 1
            i1 = i+1
            IF (alphacf (i) == 1.0) THEN
               yfunc(i) = yfunc(i1)-alambda(ntemp)*tempcf(i) &
                          /(alambda(i)*tempcf(ntemp)) &
                          *LOG(tempcf(i)/tempcf(i1))
            ELSE
               yfunc(i) = yfunc(i1)-alambda(ntemp)*tempcf(i) &
                          /((1-alphacf(i))*alambda(i)*tempcf(ntemp)) &
                          *(1-(tempcf(i)/tempcf(i1))**(alphacf(i)-1))
            ENDIF
         ENDDO
      ENDIF

      IF (inject == 'uniform') THEN
         OPEN (19,FILE='UniDisSphere_3000Points.txt',FORM='formatted')
         DO i=1,idiminj
            READ (19,*) xUD(i), yUD(i), zUD(i)
         ENDDO
         CLOSE (19)
         OPEN (19,FILE='SortedList1_3000Points.txt',FORM='formatted')
         DO i=1,idiminj
            READ (19,*) lUD(i)
         ENDDO
         CLOSE (19)
         nUD = 0
         anginj1 = 2*pi*ran1(1)
         anginj2 = 2*pi*ran1(1)
!         IF (myrank.eq.0) THEN
!            WRITE(*,*) xUD
!            WRITE(*,*)
!            WRITE(*,*) yUD
!            WRITE(*,*)
!            WRITE(*,*) zUD
!            WRITE(*,*)
!            WRITE(*,*) lUD
!            WRITE(*,*)
!            WRITE(*,*) nUD
!         ENDIF
      ENDIF

!--Define rhozerox to calculate the free fall time
!  (13 Oct 2000)
      xmax = -1.0e10
      xmin = 1.0e10
      ymax = -1.0e10
      ymin = 1.0e10
      zmax = -1.0e10
      zmin = 1.0e10
      totmas = 0.0
      DO i = 1, npart
         IF (x(i) > xmax) xmax = x(i)
         IF (x(i) < xmin) xmin = x(i)
         IF (y(i) > ymax) ymax = y(i)
         IF (y(i) < ymin) ymin = y(i)
         IF (z(i) > zmax) zmax = z(i)
         IF (z(i) < zmin) zmin = z(i)
         totmas = totmas + pmass(i)
      ENDDO
      rhozerox = MAX(rhozero,totmas/MAX(xmax-xmin,ymax-ymin, &
                     zmax-zmin)**3)
!
!--Open point mass data output file
!
!!      IF (myrank.eq.0) THEN
         WRITE (ptfile, 99990) namerun
99990    FORMAT ('P', A20)
         WRITE (accfile, 99991) namerun
99991    FORMAT ('A', A20)
!!      ENDIF
      IF (nptmass /= 0.OR.iptmass /= 0) THEN
         OPEN (iptprint, FILE = ptfile, STATUS = 'unknown', FORM = &
         'unformatted')
         OPEN (iaccpr, FILE = accfile, STATUS = 'unknown', FORM = &
         'unformatted')
      ENDIF

!--Open files for killing and reassignment of particles
!
!!      IF (myrank.eq.0) YHEN
         WRITE (killfile, 99992) namerun
99992    FORMAT ('K', A20)
         WRITE (reassfile, 99993) namerun
99993    FORMAT ('R', A20)
!!      ENDIF
      IF (ibound == 8.OR.ibound >= 90) THEN
         OPEN (ikillpr, FILE = killfile, STATUS = 'unknown', &
               FORM = 'unformatted')
         OPEN (ireasspr, FILE = reassfile, STATUS = 'unknown', &
               FORM = 'unformatted')
      ENDIF

!--Open a file for active particles in Roche lobe of secondary
!!      IF (myrank.eq.0) THEN
!!         WRITE (captfile, 99994) namerun
!         WRITE (coolfile, 99994) namerun
!99994    FORMAT ('C', A20)
!!      ENDIF
!!      IF (nptmass == 2) THEN
!!         OPEN (icaptpr, FILE=captfile, STATUS='unknown', &
!!               FORM='unformatted')
!!      ENDIF
!--Open a file for radiatively cooled particles
!         OPEN (icoolpr, FILE=coolfile, STATUS='unknown', &
!               FORM='unformatted')

!--Open notify file
!
      OPEN (inotify, FILE='notify')
!
!--Initialize du
!
      DO i=1,npart
         IF (iphase(i) /= -1) du(i) = 0.
      ENDDO
!
!--Set constant for artificial viscosity
!
      IF (alpha == 0..AND.beta == 0.) THEN
         alpha = 1.0
         beta = 2.0
      ENDIF
!
!--Set accuracy parameter for tree force calculation
!     theoretical limit for 3D tree is 0.57 = 1/SQRT(3)
!
      acc = 0.5
!!      acc = 0.3
!
!--Set stop flag
!
      istop = 0
!
!--Set min and max limit of neighbours the code tries to inforce
!
      neimin = 30
      neimax = 70
      nrange = 12
!
!--Total mass
!
      pmassmin = 1.0E+10
      fmas1 = 0.
      DO i = 1, n1
         IF (iphase (i)  >= 0) THEN
            fmas1 = fmas1 + pmass (i)
            IF (pmassmin > pmass(i)) pmassmin = pmass(i)
         ENDIF
      ENDDO
      fmas2 = 0.
      DO i=n1+1,npart
         IF (iphase(i) >= 0) THEN
            fmas2 = fmas2+pmass(i)
            IF (pmassmin > pmass(i)) pmassmin = pmass(i)
         ENDIF
      ENDDO
!
!--Set mass for particle being partially accreted at which it is
!     completely accreted
!
      pmassleast = pmassmin / 100.
!
!--Point Mass Presets
!
!--Set critical density for point mass creation
!
      rhocrea = rhozero * ptmcrit
!
!--Standardise point mass types and accretion radii
!
      DO i=1,nptmass
         iptcur = listpm(i)
         IF (initialptm < 1.OR.initialptm > 4) CALL error(where,2)
         iphase (iptcur) = initialptm
!!!!         h(iptcur) = hacc
      ENDDO
!
!--For Massive Accretion, only accrete mass and change ptmass properties
!     for certain particles - exclude these:
!
!----    Since 14 May 2004 (for a non-selfgravitating disk) (A. Okazaki)
      nnotacc = 0
      DO i=1,nptmass
         notacc(i) = .FALSE.
      ENDDO
      DO i=nptmass+1,idim
         nnotacc = nnotacc+1
         notacc(i) = .TRUE.
      ENDDO
!      nnotacc = 0
!      IF (myrank.eq.0) WRITE (iprint,*) &
!                   '********* IDELAYACC = 1 *********'
!      OPEN (22, FILE='WANTRUN')
! 50   READ (22,*, END=100) inum
!      nnotacc = nnotacc + 1
!      notacc(inum) = .TRUE.
!      GOTO 50
! 100  CLOSE(22)
!
!--Gravitational softening for ptmass-ptmass interactions
!
      iptsoft = 0
      ptsoft = 1.0E-2
!
!--Set psoft for softening the gravitational potential when the 1/(r+pso
!     potential law is used (i.e. when igrape=1, or isoft=1).
!
      psoft = 0.01
!
!--Value of minimum h in order to save computing time (if hmin is differ
!     from 0 then program does not follow high density regions accuratel
!
      hmin = 0.
!      hmin = 0.01

      IF (isoft == 1 .AND. hmin < psoft) hmin = psoft
!
!--Compute tables for kernel quantities
!
      CALL ktable

      END SUBROUTINE preset
