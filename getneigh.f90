      SUBROUTINE getneigh(ipart, ntot, hi, x, y, z, nlist, list)

      use idims

      use phase
      use kerne

      implicit none

      INTEGER(I4B) :: list(idim)
      INTEGER(I4B) :: ipart, ntot, nlist, i
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: hi, dx, dy, dz, r2

      write (9,*) 'in getneigh ',ipart, ntot, hi,radkernel
      nlist = 0
      DO i = 1, ntot
         IF (iphase(i).EQ.0) THEN
            dx = x(ipart) - x(i)
            dy = y(ipart) - y(i)
            dz = z(ipart) - z(i)
            r2 = dx**2 + dy**2 + dz**2
            IF (r2.LT.radkernel*radkernel*hi*hi .AND. i.NE.ipart) THEN
               nlist = nlist + 1
               list(nlist) = i
            ENDIF
         ENDIF
      END DO

      END SUBROUTINE getneigh
