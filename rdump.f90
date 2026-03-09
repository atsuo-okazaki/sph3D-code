      SUBROUTINE rdump(idisk1, ichkl)
!************************************************************
!                                                           *
!  This routine reads a dump into memory                    *
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
      use polyk2
      use phase
      use ptmass
      use binary
      use torq
      use timei
      use stepopt
      use new
      use debug
      use ptbin
      use maslos
      use capt
      use winds

      implicit none

      INTEGER(I4B), parameter :: icall=2
      INTEGER(I4B) :: idisk1, ichkl, i, ifactor, imaxnew, iminnew, ipower, j
      REAL(DP) :: umassi, udisti, utimei
      REAL(DP) :: dtmaxdp, gm1
      CHARACTER(len=7) :: where='rdump'

!      DATA icall/2/
!
!--Read
!
      IF (itrace.EQ.'all') WRITE (*, 99001)
99001 FORMAT (' entry subroutine rdump')

      ichkl = 0
      READ (idisk1, END=100) udisti, umassi, utimei, &
           npart, n1, n2, gt, gamma, rhozero, RK21, &
           (h(i), i=1, npart), escap, tkin, tgrav, tterm, &
           (x(i), i=1, npart), (y(i), i=1, npart), &
           (z(i), i=1, npart), (vx(i), i=1, npart), &
           (vy(i), i=1, npart), (vz(i), i=1, npart), &
           (u(i), i=1, npart), (pmass(i), i=1, npart), &
           (rho(i), i=1, npart), (dgrav(i), i=1, npart), &
           dtmaxdp, (isteps(i), i=1, npart) &
           ,(iphase(i), i=1, npart), &
           nptmass, (listpm(i), i=1, nptmass), &
           (spinx(i),i=1,nptmass), (spiny(i),i=1,nptmass), &
           (spinz(i),i=1,nptmass) &
           ,(angaddx(i),i=1,nptmass), (angaddy(i),i=1,nptmass), &
           (angaddz(i),i=1,nptmass), &
           anglostx, anglosty, anglostz, &
           nreassign, naccrete, nkill, specang, ptmassin, &
           (spinadx(i),i=1,nptmass),(spinady(i),i=1,nptmass), &
           (spinadz(i),i=1,nptmass), &
           (torqt(i), i=1, npart), (torqg(i), i=1, npart), &
           (torqp(i), i=1, npart),(torqv(i), i=1, npart), &
           (torqc(i), i=1, npart),(ibelong(i),i=1,npart), &
           (iantigr(i),i=1,npart)

!      DO i = 1, npart
!         iphase(i) = 0
!      END DO
!      nptmass = 0
!      anglostx = 0.
!      anglosty = 0.
!      anglostz = 0.
!      nreassign = 0
!      naccrete = 0
!      nkill = 0
!      specang = 0.
!      ptmassin = 0.

      gtdouble = DBLE(gt)

!-- thermal
      IF (encal.EQ.'i') THEN
         thermal1 = RK21
      ELSE IF (encal.EQ.'c') THEN
!----    thermal gives the floor temperature for simulations
!        with radiative cooling (in radcool.f)
         thermal1 = RK21
      ELSE
         gm1 = gamma - 1.0
         thermal1 = RK21*(rhozero**gm1)
      ENDIF

!
!--Zero torques
!
      DO i = 1, idim
         torqt(i) = 0.
         torqg(i) = 0.
         torqp(i) = 0.
         torqv(i) = 0.
         torqc(i) = 0.
      END DO
!
!--Check units in file the same as in the code!
!
      IF (udisti.LT.0.99999*udist .OR. udisti.GT.1.00001*udist) THEN
         CALL error(where,1)
      ELSEIF (umassi.LT.0.99999*umass .OR.umassi.GT.1.00001*umass) THEN
         CALL error(where,2)
      ENDIF
      IF (npart.GT.idim) THEN
         CALL error(where,3)
      ENDIF
!
!--Check that dtmax times are the same.  If not, modify isteps(i) as in mesop.f
!
!cc      GOTO 50

      IF (gt.NE.0.0 .AND. &
           (dtmaxdp.LT.0.9999*dtmax .OR. dtmaxdp.GT.1.0001*dtmax)) THEN
         ipower = INT(LOG10(dtmax/dtmaxdp)/LOG10(2.0))

         ifactor = 2**ABS(ipower)
         imaxnew = imaxstep/ifactor
         iminnew = 2*ifactor

         IF (ipower.LT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  isteps(j) = MIN(isteps(j), imaxnew)
                  isteps(j) = isteps(j)*ifactor
               ENDIF
            END DO
         ELSEIF (ipower.GT.0) THEN
            DO j = 1, npart
               IF (iphase(j).NE.-1) THEN
                  IF (isteps(j)/ifactor .LT. 2) CALL error(where, 4)
                  isteps(j) = isteps(j)/ifactor
               ENDIF
            END DO
         ENDIF
      ENDIF
!
!--Change reference frame
!
 50   IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2.AND.gt.NE.0.0)) THEN
!      IF (iexpan.NE.0.OR.(ifcor.GT.0.AND.ifcor.LE.2)) THEN
         CALL chanref(icall)
      ELSEIF (ifcor.GT.2) THEN
         ifcor = ifcor - 2
      ENDIF

!      DO i=1,npart
!         IF (iphase(i).EQ.0) pmass(i)=pmass(i)*10.0
!         IF (iphase(i).EQ.0) u(i)=u(i)*8.0
!         IF (iphase(i).EQ.0) h(i)=2.0*h(i)
!      END DO

!---- If nptmass=2 and iaccevol='f', then
!     replace h(i) for point masses with hacc.
!     (This is to enable the code to work when iaccevol is changed
!     from 'v' to 'f' in the middle of the simulation.)
!     (1 September 2003, A. Okazaki)
      IF (nptmass.EQ.2 .AND. iaccevol.EQ.'f') THEN
         h(1) = hacc
         h(2) = hacc
      ENDIF

      IF (nptmass.EQ.2) THEN
         ncapt = 0
         DO i=1,npart
            iinold(i) = 0
         ENDDO
      ENDIF

!c      IF (ibound.EQ.99) THEN
!c         sinj(1) = REAL(nreassign)
!c      ENDIF

      IF (itrace.EQ.'all') WRITE (*, 99002)
99002 FORMAT (' exit subroutine rdump')

      RETURN

 100  ichkl = 1

      END SUBROUTINE rdump

