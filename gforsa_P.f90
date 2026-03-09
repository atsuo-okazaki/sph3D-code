      SUBROUTINE gforsa(m, nlistga, listga, fsx, fsy, fsz, epot)
!************************************************************
!                                                           *
!  Subroutine by W. Press (11/21/86).  Evaluates force on   *
!     particle M due to a list of other particles.          *
!     THIS ROUTINE VECTORIZABLE.                            *
!                                                           *
!************************************************************

      use idims

!---- treecom.mod <-- treecom.f90 for non-parallel simulation
!---- treecom.mod <-- treecom_P.f90 for parallel simulation
      use treecom
      use logun
      use soft

      implicit none

      INTEGER(I4B) :: listga(idim), m, nlistga, n, j
      REAL(DP) :: fsx, fsy, fsz, epot, rrx, rry, rrz, &
           difx, dify, difz, rr, rr05, fff, potn

      rrx = rx(m)
      rry = ry(m)
      rrz = rz(m)

      DO 101 j = 1, nlistga
         n = listga(j)
         difx = rx(n) - rrx
         dify = ry(n) - rry
         difz = rz(n) - rrz

         IF (isoft.EQ.1) THEN
            rr = difx**2 + dify**2 + difz**2 + psoft**2
         ELSE
            rr = difx**2 + dify**2 + difz**2 + tiny
         ENDIF

         rr05 = SQRT(rr)
!
!--The force definition
!
         fff = em(n)/(rr*rr05)
         potn = em(n)/rr05

         fsx = fsx + fff*difx
         fsy = fsy + fff*dify
         fsz = fsz + fff*difz
         epot = epot + potn
 101  CONTINUE

      END SUBROUTINE gforsa
