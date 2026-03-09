      SUBROUTINE inform (where)
!************************************************************
!                                                           *
!  This routine computes relevant quantities for print out  *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use tming
      use btree
      use part
      use densi
      use typef
      use carac
      use kerne
      use out
      use bodys
      use logun
      use debug
      use ptdump
      use phase
      use active
      use ptmass
      use maslos
      use winds

      implicit none

      INTEGER(I4B) :: i, imin, j, n1new
      REAL(DP) :: d2, dz, pmassi, rmin2
      CHARACTER (len=7) :: where
!
!--Allow for tracing flow
!
      IF (myrank.eq.0) then
         IF (itrace.EQ.'all') WRITE (iprint, 99001)
      ENDIF
99001 FORMAT (' entry subroutine inform')
!
!--Compute : center of mass, velocity of cm, mean density and dispersion
!     for first object
!
      cmx1 = 0.
      cmy1 = 0.
      cmz1 = 0.
      vcmx1 = 0.
      vcmy1 = 0.
      vcmz1 = 0.
      cmx2 = 0.
      cmy2 = 0.
      cmz2 = 0.
      vcmx2 = 0.
      vcmy2 = 0.
      vcmz2 = 0.

      romean1 = 0.
      romax1 = 0.
      rocen1 = 0.
      romean2 = 0.
      romax2 = 0.
      rocen2 = 0.
      hmi1 = 1.E30
      hma1 = 0.
      hmi2 = 1.E30
      hma2 = 0.

      fmas1 = 0.
      fmas2 = 0.

      n1new = 0
      DO i = 1, n1
         IF (iphase(i) .EQ.0) THEN
            n1new = n1new+1
            pmassi = pmass(i)
            fmas1 = fmas1+pmassi
            cmx1 = cmx1+pmassi*x(i)
            cmy1 = cmy1+pmassi*y(i)
            cmz1 = cmz1+pmassi*z(i)

            vcmx1 = vcmx1+pmassi*vx(i)
            vcmy1 = vcmy1+pmassi*vy(i)
            vcmz1 = vcmz1+pmassi*vz(i)

            hmi1 = MIN(hmi1, h(i) )
            hma1 = MAX(hma1, h(i) )

            romean1 = romean1+rho(i)
            romax1 = MAX(romax1,rho(i))
         ENDIF
      ENDDO
      DO j = 1, nptmass
         i = listpm(j)
         pmassi = pmass(i)
         fmas1 = fmas1+pmassi
         cmx1 = cmx1+pmassi*x(i)
         cmy1 = cmy1+pmassi*y(i)
         cmz1 = cmz1+pmassi*z(i)

         vcmx1 = vcmx1+pmassi*vx(i)
         vcmy1 = vcmy1+pmassi*vy(i)
         vcmz1 = vcmz1+pmassi*vz(i)
      ENDDO

      cmx1 = cmx1/fmas1
      cmy1 = cmy1/fmas1
      cmz1 = cmz1/fmas1

      vcmx1 = vcmx1/fmas1
      vcmy1 = vcmy1/fmas1
      vcmz1 = vcmz1/fmas1

      romean1 = romean1/n1new
!
!--If doing accretion on to binary, then recentre the centre of mass and
!     set the centre of mass velocity to zero
!
      IF (ibound.EQ.8 .OR. ibound.GE.90) THEN
         IF (MOD(ncount,nprout) .EQ.0.OR.MOD(ncount,nstep)  &
         .EQ.0 .OR. iptcreat.EQ.1) THEN
            IF (myrank.eq.0) WRITE (iprint,*) 'ZERO CENTRE OF MASS'
         ENDIF
         DO i = 1, npart
            IF (iphase(i).GE.0) THEN
               x(i) = x(i)-cmx1
               y(i) = y(i)-cmy1
               z(i) = z(i)-cmz1
               vx(i) = vx(i)-vcmx1
               vy(i) = vy(i)-vcmy1
               vz(i) = vz(i)-vcmz1
            ENDIF
         ENDDO
      ENDIF
!
!--Compute : center of mass, velocity of cm, mean density and dispersion
!     for second object (if existing)
!
      IF (n2.NE.0) THEN
         DO i = n1+1, npart
            IF (iphase(i).GE.0) THEN
               pmassi = pmass(i)
               fmas2 = fmas2+pmass(i)
               cmx2 = cmx2+pmassi*x(i)
               cmy2 = cmy2+pmassi*y(i)
               cmz2 = cmz2+pmassi*z(i)

               vcmx2 = vcmx2+pmassi*vx(i)
               vcmy2 = vcmy2+pmassi*vy(i)
               vcmz2 = vcmz2+pmassi*vz(i)

               hmi2 = MIN(hmi2,h(i))
               hma2 = MAX(hma2,h(i))

               romean2 = romean2+rho(i)
               romax2 = MAX(romax2,rho(i))
            ENDIF
         ENDDO

         cmx2 = cmx2/fmas2
         cmy2 = cmy2/fmas2
         cmz2 = cmz2/fmas2

         vcmx2 = vcmx2/fmas2
         vcmy2 = vcmy2/fmas2
         vcmz2 = vcmz2/fmas2

         romean2 = romean2/n2
      ENDIF
!
!--Compute maximum distance
!
      dmax1 = 0.
      zmax1 = 0.
      rmin2 = 1.0E+30
      DO i = 1, n1
         IF (iphase(i).GE.0) THEN
            dz = z(i)-cmz1
            d2 = (x(i)-cmx1)**2+(y(i)-cmy1)**2+dz*dz
            dmax1 = MAX(dmax1,d2)
            zmax1 = MAX(zmax1,ABS(dz))
            IF (d2.LT.rmin2 .AND. iphase(i).EQ.0) THEN
               imin = i
               rmin2 = d2
            ENDIF
         ENDIF
      ENDDO
      dmax1 = SQRT(dmax1)
      rocen1 = rho(imin)

      dmax2 = 0.
      zmax2 = 0.
      rmin2 = 1.0E+30
      DO i = n1+1, npart
         IF (iphase(i).GE.0) THEN
            dz = z(i)-cmz2
            d2 = (x(i)-cmx2)**2+(y(i)-cmy2)**2+dz*dz
            dmax2 = MAX(dmax2,d2)
            zmax2 = MAX(zmax2,ABS(dz))
            IF (d2.LT.rmin2 .AND. iphase(i).EQ.0) THEN
               imin = i
               rmin2 = d2
            ENDIF
         ENDIF
      ENDDO
      dmax2 = SQRT(dmax2)
      IF (npart.GT.n1) rocen2=rho (imin)
!
!--Compute energies
!
      CALL toten
!
!--Compute total angular momentum
!
      CALL angmom

      IF (where(1:6).NE.'newrun') THEN
!
!--Write dump on disk
!
         IF (nstep.LT.1) nstep = 1
!         IF (nbuild.EQ.1 .OR. MOD(ncount,nstep).EQ.0 .OR.
         IF (MOD(ncount,nstep).EQ.0 .OR. iptcreat.EQ.1) THEN
            CALL file
            CALL wdump(idisk1)

            IF (ibound.EQ.94 .OR. ibound.EQ.96) THEN
               DO i=1,nptmass
                  sinj0(i) = sinj(i)
               ENDDO
            ELSEIF (ibound.EQ.95 .OR. ibound.EQ.97) THEN
               !DO i=1,nptmass+1
               DO i=1,3
                  sinj0(i) = sinj(i)
               ENDDO
            ELSEIF (ibound.EQ.99) THEN
               sinj0(1) = sinj(1)
            ENDIF
         ENDIF
!
!--Write global results on listing
!
         IF (MOD(ncount,nprout).EQ.0 .OR. MOD(ncount,nstep) &
         .EQ.0.OR.iptcreat.EQ.1) THEN
            where = 'inform'
            CALL prout(where)
         ENDIF

         iptcreat = 0
!
!--Update input file
!
         CALL wrinsph

      ENDIF

      IF (idebug(1:6).EQ.'inform') THEN
         IF (myrank.eq.0) WRITE (iprint,99002) cmx1, cmy1, cmz1, &
            vcmx1, vcmy1, vcmz1, romean1, zmax1
99002    FORMAT(1X,5(1PE12.5,1X))
      ENDIF

      END SUBROUTINE inform
