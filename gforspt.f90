      SUBROUTINE gforspt(x, y, z, vx, vy, vz, h, npart)
!************************************************************
!                                                           *
!  Subroutine by IAB, MRB 1994.  Evaluates forces on point  *
!     mass due to all other particles.                      *
!     Returns list of nearest neighbours to point mass !    *
!     THIS ROUTINE VECTORIZABLE. (?)                        *
!                                                           *
!************************************************************

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

      CALL getused(tgforpt1)
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT(' entry subroutine gforspt')

      DO k = 1, nptmass
         iptcur = listpm(k)
         IF (iscurrent(iptcur).EQ.1) THEN
            gravx(iptcur) = 0.
            gravy(iptcur) = 0.
            gravz(iptcur) = 0.
            poten(iptcur) = 0.
            dphit(iptcur) = 0.
         ENDIF
      END DO

      IF (nptmass.GE.1) THEN
         DO k = 1, npart
            isnearpt(k) = 0
         END DO
      ENDIF

      DO iptn = 1, nptmass
         ipt = listpm(iptn)
         hipt = h(ipt)
!----    For point mass 1 (= Be star), hipt should be
!        the stellar radius, not the outer accretion radius
!        for the other point mass(es).
!        (14 September 2003, A. Okazaki)
         IF (ipt.EQ.1) THEN
            IF (ibound.EQ.0 .OR. &
               (ibound.GE.90 .AND. ibound.LE.93)) THEN
!-- Have a margin which is a bit larger than a typical
!   smoothing length of particles near the point mass
!!               hipt = rptmas(ipt)
               hipt = rptmas(ipt) * 0.9
            ENDIF
         ENDIF

         rrx = x(ipt)
         rry = y(ipt)
         rrz = z(ipt)
         vxipt = vx(ipt)
         vyipt = vy(ipt)
         vzipt = vz(ipt)
         neigh = 0
         DO 101 j = 1, npart
            IF(iphase(j).GE.1 .AND. j.LE.ipt &
                                      .OR. iphase(j).EQ.-1) GOTO 101
            pmassj = pmass(j)

            difx = x(j) - rrx
            dify = y(j) - rry
            difz = z(j) - rrz
            difvx = vx(j) - vxipt
            difvy = vy(j) - vyipt
            difvz = vz(j) - vzipt

            rr = difx**2 + dify**2 + difz**2 + tiny
!
!--Check to see if neighbour - NO POINT MASSES AS NEIGHBOURS
!
            hmean2 = ((hipt + h(j))/2.)**2
            IF (hipt.GT.h(j) .AND. iphase(ipt).EQ.4) THEN
               hmean2 = ((hipt + 2.0*h(j))/2.0)**2
            ENDIF
            IF (rr.LT.(4.*hmean2) .AND. iphase(j).EQ.0) THEN
               neigh = neigh + 1
               IF (neigh.GT.iptneigh) CALL error(where,1)
               nearpt(iptn,neigh) = j
!--Added 19/10/94
               IF (isnearpt(j).EQ.0) THEN
                  isnearpt(j) = iptn
               ELSE
                  isnearpt(j) = -1
               ENDIF
            ENDIF
!
!--Add forces
!
!
!--The force definition:
!
            IF(iphase(j).EQ.0 .OR. iptsoft.EQ.0) THEN
               rr05 = SQRT(rr)
               fff = pmassj/(rr*rr05)
               potn = pmassj/rr05

!--This turned out to be helpless. (27/11/00)
!--Added 24/11/00
!--Add a very-steep, short-range potential to prevent particles
!  from accreting onto point mass 1 (=primary star of a binary)
!!               IF (ipt.EQ.1) THEN
!!cc                  sclgth = 0.05*rptmas(ipt)
!!                  sclgth = 0.05*rmind
!!cc                  expr = EXP((rptmas(ipt)-rr05)/sclgth)
!!                  expr = EXP((rmind-rr05)/sclgth)
!!                  potxpl = -pmassj*sclgth * 10.0
!!                  potnadd = potxpl*expr
!!                  fffadd = potxpl/sclgth*expr
!!                  potn = potn + potnadd
!!                  fff = fff+fffadd
!!               ENDIF
            ELSEIF (iphase(j).GE.1) THEN
!
!--Non-pointmass force
!
!
!--Softened potential 1
!
               IF(iptsoft.EQ.1) THEN
                  rr = rr + ptsoft*ptsoft
                  rrs05 = SQRT(rr)
                  rr32 = rr*rrs05
                  potn = pmassj / rrs05
                  fff = pmassj / rr32
!
!--Softened potential 2
!
               ELSEIF (iptsoft.EQ.2) THEN
                  rr4 = rr*rr + ptsoft*ptsoft*ptsoft*ptsoft
                  rr4s025 = (rr4)**0.25
                  rr54 = rr4*rr4s025
                  potn = pmassj / rr4s025
                  fff = rr*pmassj / rr54
               ENDIF
            ENDIF
!
!--End force definition
!
            IF (iscurrent(ipt).EQ.1 .AND. (.NOT.(notacc(j)))) THEN
               gravx(ipt) = gravx(ipt) + fff*difx
               gravy(ipt) = gravy(ipt) + fff*dify
               gravz(ipt) = gravz(ipt) + fff*difz
               poten(ipt) = poten(ipt) - potn
            ENDIF
            IF (iscurrent(j).EQ.1) THEN
               fmass = pmass(ipt)/pmassj
               gravx(j) = gravx(j) - fmass*fff*difx
               gravy(j) = gravy(j) - fmass*fff*dify
               gravz(j) = gravz(j) - fmass*fff*difz
               poten(j) = poten(j) - fmass*potn
            ENDIF

 101     CONTINUE
         nptlist(iptn) = neigh
      END DO


!!      DO j = 1, npart
!!         WRITE (iprint, *) 'poten(',j,')=',poten(j)
!!      ENDDO

      CALL getused(tgforpt2)
      tgforpt = tgforpt + (tgforpt2 - tgforpt1)

      END SUBROUTINE gforspt
