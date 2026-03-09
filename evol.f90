      SUBROUTINE evol
!************************************************************
!                                                           *
!  This routine drives the time evolution of the system.    *
!                                                           *
!************************************************************

      use idims

      use tming
      use logun
      use part
      use btree
      use useles
      use timei
      use nextmpt
      use phase
      use active
      use ptdump
      use debpt
      use isnpt
      use ptmass
      use accurpt
      use carac
      use torq
      use debugit
      use secret

      implicit none

      INTEGER(I4B) :: i, icount, iseed, j
      REAL(DP) :: ran1
      REAL(DP) :: pmasspt, rnd
      CHARACTER(len=7) :: where='evol'
!
!--Define options of this run
!
!      imax = 2097152
      imax = 1073741824
      imaxstep = imax/2

      CALL options
!
!--Set all quantities needed for the run
!
      CALL preset
!
!--Initialise random number generator
!
      iseed = -4357
      rnd = ran1(iseed)
      WRITE(iprint,*)'Random seed: ',iseed
!
!--Set number of active particles
!
      nactive = 0
      DO i = 1, npart
         IF (iphase(i).NE.-1) nactive = nactive + 1
      END DO

      DO i = 1, nptmass
         j = listpm(i)
         pmasspt = pmass(j)
         ptmsyn(i) = pmasspt
         ptmadd(i) = 0.0
         xmomsyn(i) = pmasspt*vx(j)
         ymomsyn(i) = pmasspt*vy(j)
         zmomsyn(i) = pmasspt*vz(j)
         xmomadd(i) = 0.0
         ymomadd(i) = 0.0
         zmomadd(i) = 0.0
      END DO
!
!--Write all quantities on listings
!
      CALL header(where)
!
!---------------------------
!---- E V O L U T I O N ----
!---------------------------
!
      ncount = 0
      icount = 0
      nbuild = 0
      iaccr = 0
      iptout = 0
      iptdump = 0
      isheld = .FALSE.
      tkeep = 4.0
!cc      tkeep = 0.0

      itest = 0

      DO i = 1, idim
         isnearpt(i) = 0
         torqt(i) = 0.0
         torqg(i) = 0.0
         torqp(i) = 0.0
         torqv(i) = 0.0
         torqc(i) = 0.0
      END DO

      DO i = 1, 1000000
!
!--Evolve one timestep
!
         CALL integs
!
!--Check if saving time has come
!
         CALL save

      END DO

      END SUBROUTINE evol
