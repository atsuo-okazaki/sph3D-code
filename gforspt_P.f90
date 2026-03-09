      SUBROUTINE gforspt (x, y, z, vx, vy, vz, h, npart)
!************************************************************
!                                                           *
!  Subroutine by IAB, MRB 1994.  Evaluates forces on point  *
!     mass due to all other particles.                      *
!     Returns list of nearest neighbours to point mass !    *
!     THIS ROUTINE VECTORIZABLE. (?)                        *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use phase
      use ptsoftx
      use logun
      use gravi
      use ener
      use carac
      use ptmass
      use neighbor
      use current
      use isnpt
      use debug
      use perform
      use delay
      use diskbd
      use typef
      use zzhp

      implicit none

      INTEGER(I4B) :: npart, j, k, iptcur, ipt, iptn, neigh
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: vx(idim), vy(idim), vz(idim), h(idim)
      REAL(DP) :: difvx, difvy, difvz, fff, fmass, hipt, hmean2, &
                pmassj, potn, rr, rr05, rr32, rr4, rr4s025, rr54, rrs05, &
                rrx, rry, rrz, difx, dify, difz, vxipt, vyipt, vzipt
      CHARACTER(len=7) :: where='gforspt'

      INTEGER(I4B), parameter :: ip = 64
      INTEGER(I4B) :: i, ii
      INTEGER(I4B) :: l1, l2, ls (ip), le (ip)
      INTEGER(I4B) :: neigh2 (ip)
      INTEGER(I4B) :: nearpt_work (nptmass, npart, ip)
      REAL(DP) :: grav_work (4)

      CALL getused (tgforpt1)
!
!--Allow for tracing flow
!
      IF (myrank.eq.0) then
         IF (itrace.EQ.'all') WRITE (iprint, 99001)
      ENDIF
99001 FORMAT(' entry subroutine gforspt')

      DO k = 1, nptmass
      iptcur = listpm (k)
      IF (iscurrent (iptcur) .EQ.1) THEN
         gravx (iptcur) = 0.
         gravy (iptcur) = 0.
         gravz (iptcur) = 0.
         poten (iptcur) = 0.
         dphit (iptcur) = 0.
      ENDIF
      ENDDO

      IF (nptmass.GE.1) THEN
         DO k = 1, npart
         isnearpt (k) = 0
         ENDDO
      ENDIF

      l1 = npart / ip
      l2 = mod (npart, ip)
      DO i = 1, ip
      IF (i - 1.lt.l2) then
         ls (i) = (l1 + 1) * (i - 1) + 1
         le (i) = (l1 + 1) * i
      ELSE
         ls (i) = l1 * (i - 1 - l2) + 1 + (l1 + 1) * l2
         le (i) = l1 * (i - 1 - l2 + 1) + (l1 + 1) * l2
      ENDIF
      enddo

      DO iptn = 1, nptmass
      ipt = listpm (iptn)
      hipt = h (ipt)
!----    For point mass 1 (= Be star), hipt should be
!        the stellar radius, not the outer accretion radius
!        for the other point mass(es).
!        (14 September 2003, A. Okazaki)
      IF (ipt.EQ.1) THEN
         IF (ibound.EQ.0.OR. (ibound.GE.90.AND.ibound.LE.97) ) THEN
            hipt = rptmas(ipt)
         ENDIF
      ENDIF

      rrx = x (ipt)
      rry = y (ipt)
      rrz = z (ipt)
      vxipt = vx (ipt)
      vyipt = vy (ipt)
      vzipt = vz (ipt)
!!!         neigh = 0
      grav_work (1:4) = 0.d0
!POPTION PARALLEL,INDEP
!POPTION TLOCAL(fff,potn)
!POPTION PSUM(grav_work)
      DO i = 1, ip
      neigh2 (i) = 0
      DO 101 j = ls (i), le (i)
!         DO 101 j = 1, npart
         IF (iphase (j) .GE.1.AND.j.LE.ipt.OR.iphase (j) .EQ. - 1) GOTO &
         101
         pmassj = pmass (j)

         difx = x (j) - rrx
         dify = y (j) - rry
         difz = z (j) - rrz
         difvx = vx (j) - vxipt
         difvy = vy (j) - vyipt
         difvz = vz (j) - vzipt

         rr = difx**2 + dify**2 + difz**2 + tiny
!
!--Check to see if neighbour - NO POINT MASSES AS NEIGHBOURS
!
         hmean2 = ( (hipt + h (j) ) / 2.) **2
         IF (hipt.GT.h (j) .AND.iphase (ipt) .EQ.4) THEN
            hmean2 = ( (hipt + 2.0 * h (j) ) / 2.0) **2
         ENDIF
         IF (rr.LT. (4. * hmean2) .AND.iphase (j) .EQ.0) THEN
!!!               neigh = neigh + 1
            neigh2 (i) = neigh2 (i) + 1
            nearpt_work (iptn, neigh2 (i), i) = j
!!!               IF (neigh.GT.iptneigh) CALL error(where,1)
!!!               nearpt(iptn,neigh) = j
!--Added 19/10/94
            IF (isnearpt (j) .EQ.0) THEN
               isnearpt (j) = iptn
            ELSE
               isnearpt (j) = - 1
            ENDIF
         ENDIF
!
!--Add forces
!
!
!--The force definition:
!
         IF (iphase (j) .EQ.0.OR.iptsoft.EQ.0) THEN
            rr05 = SQRT (rr)
            fff = pmassj / (rr * rr05)
            potn = pmassj / rr05

         ELSEIF (iphase (j) .GE.1) THEN
!
!--Non-pointmass force
!
!
!--Softened potential 1
!
            IF (iptsoft.EQ.1) THEN
               rr = rr + ptsoft * ptsoft
               rrs05 = SQRT (rr)
               rr32 = rr * rrs05
               potn = pmassj / rrs05
               fff = pmassj / rr32
!
!--Softened potential 2
!
            ELSEIF (iptsoft.EQ.2) THEN
               rr4 = rr * rr + ptsoft * ptsoft * ptsoft * ptsoft
               rr4s025 = (rr4) **0.25
               rr54 = rr4 * rr4s025
               potn = pmassj / rr4s025
               fff = rr * pmassj / rr54
            ENDIF
         ENDIF
!
!--End force definition
!
         IF (iscurrent (ipt) .EQ.1.AND. (.NOT. (notacc (j) ) ) ) THEN
!!!               gravx(ipt) = gravx(ipt) + fff*difx
!!!               gravy(ipt) = gravy(ipt) + fff*dify
!!!               gravz(ipt) = gravz(ipt) + fff*difz
!!!               poten(ipt) = poten(ipt) - potn
            grav_work (1) = grav_work (1) + fff * difx
            grav_work (2) = grav_work (2) + fff * dify
            grav_work (3) = grav_work (3) + fff * difz
            grav_work (4) = grav_work (4) - potn
         ENDIF
         IF (iscurrent (j) .EQ.1) THEN
            fmass = pmass (ipt) / pmassj
            gravx (j) = gravx (j) - fmass * fff * difx
            gravy (j) = gravy (j) - fmass * fff * dify
            gravz (j) = gravz (j) - fmass * fff * difz
            poten (j) = poten (j) - fmass * potn
         ENDIF

  101 END DO
      enddo

      gravx (ipt) = gravx (ipt) + grav_work (1)
      gravy (ipt) = gravy (ipt) + grav_work (2)
      gravz (ipt) = gravz (ipt) + grav_work (3)
      poten (ipt) = poten (ipt) + grav_work (4)
!!!         nptlist(iptn) = neigh
      nptlist (iptn) = 0
      DO i = 1, ip
      nptlist (iptn) = nptlist (iptn) + neigh2 (i)
      IF (nptlist (iptn) .GT.iptneigh) CALL error (where, 1)
      enddo

      ii = 0
!POPTION PARALLEL
      DO i = 1, ip
      DO j = 1, neigh2 (i)
      ii = ii + 1
      nearpt (iptn, ii) = nearpt_work (iptn, j, i)
      enddo
      enddo

      ENDDO

      CALL getused (tgforpt2)

      tgforpt = tgforpt + (tgforpt2 - tgforpt1)

      RETURN
      END SUBROUTINE gforspt
