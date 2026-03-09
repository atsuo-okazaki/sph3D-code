      SUBROUTINE gforsn(m, nlistgn, listgn, fsx, fsy, fsz, epot)
!************************************************************
!                                                           *
!  Subroutine by W. Press (11/21/86).  Evaluates force on   *
!     particle M due to a list of composite nodes.          *
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

      INTEGER(I4B) :: listgn(idim), m, nlistgn, n, j
      REAL(DP) :: fsx, fsy, fsz, epot, rrx, rry, rrz, &
           difx, dify, difz, rr, rri, rri05, f, fpr, fpprr, &
           tx, ty, tz, rqr, fac, fx, fy, fz, phim, phiq

      rrx = rx(m)
      rry = ry(m)
      rrz = rz(m)

      DO j = 1, nlistgn
         n = listgn(j)
         difx = rx(n) - rrx
         dify = ry(n) - rry
         difz = rz(n) - rrz

         IF (isoft.EQ.1) THEN
            rr = difx**2 + dify**2 + difz**2 + psoft**2
         ELSE
            rr = difx**2 + dify**2 + difz**2 + tiny
         ENDIF
!
!--The force definition
!
         rri = 1./rr
         rri05 = SQRT(rri)
         f = rri*rri05
         fpr = ( - 3.0)*f*rri
         fpprr = ( - 4.0)*fpr*rri

         tx = qxx(n)*difx + qxy(n)*dify + qzx(n)*difz
         ty = qxy(n)*difx + qyy(n)*dify + qyz(n)*difz
         tz = qzx(n)*difx + qyz(n)*dify + qzz(n)*difz
         rqr = difx*tx + dify*ty + difz*tz
         fac = (em(n)*f + 0.5*(rqr*(fpprr-fpr/rr)+ &
                    (qxx(n)+qyy(n)+qzz(n))*fpr))
         fx = fac*difx + fpr*tx
         fy = fac*dify + fpr*ty
         fz = fac*difz + fpr*tz
         fsx = fsx + fx
         fsy = fsy + fy
         fsz = fsz + fz
         phim = em(n)*rri05
         phiq = 0.5*rqr*rri*rri*rri05
         epot = epot + phim + phiq
      END DO

      END SUBROUTINE gforsn
