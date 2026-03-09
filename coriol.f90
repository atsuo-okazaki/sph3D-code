      SUBROUTINE coriol(ipart, ti, x, y, z, vx, vy, vz)
!************************************************************
!                                                           *
!  This routine adds centrifugal and coriolis forces        *
!                                                           *
!************************************************************

      use idims

      use constants
      use rotat
      use force, only: fx, fy, fz, du, dh
      use logun
      use debug
      use typef

      implicit none

      INTEGER(I4B) :: ipart
      REAL(DP) :: x(idim), y(idim), z(idim), vx(idim), &
                        vy(idim), vz(idim)
      REAL(DP) :: ti, a3, a32, aadot, rscale, drdt, dlnrdt, &
               fcentx, fcenty, fcentz, fcorcx, fcorcy, fcorcz, omeg2
!
!--Scaling factors
!
      CALL scaling(ti, rscale, drdt, dlnrdt)
      a3 = rscale**3
      a32 = 2.*omeg0*a3
      aadot = 2.*omeg0*drdt*rscale**2
      omeg2 = omeg0**2
      IF (ifcor.EQ.1) THEN
!
!--Rotation around z axis
!
!--Coriolis
!
         fcorcx = a32*vy(ipart) + aadot*y(ipart)
         fcorcy = -a32*vx(ipart) - aadot*x(ipart)
!
!--Centrifugal
!
         fcentx = a3*omeg2*x(ipart)
         fcenty = a3*omeg2*y(ipart)

         fx(ipart) = (fx(ipart)) + fcorcx + fcentx
         fy(ipart) = (fy(ipart)) + fcorcy + fcenty

         IF (idebug(1:6).EQ.'coriol') THEN
            WRITE (iprint, 99002) fx(ipart)
            WRITE (iprint, 99002) fy(ipart)
         ENDIF

      ELSEIF (ifcor.EQ.2) THEN
!
!--Rotation end over end around x axis
!
!--Coriolis
!
         fcorcy = a32*vz(ipart) + aadot*z(ipart)
         fcorcz = -a32*vy(ipart) - aadot*y(ipart)
!
!--Centrifugal
!
         fcenty = a3*omeg2*y(ipart)
         fcentz = a3*omeg2*z(ipart)

         fy(ipart) = (fy(ipart)) + fcorcy + fcenty
         fz(ipart) = (fz(ipart)) + fcorcz + fcentz

         IF (idebug(1:6).EQ.'coriol') THEN
            WRITE (iprint, 99002) fy(ipart)
            WRITE (iprint, 99002) fz(ipart)
         ENDIF

      ENDIF

99002 FORMAT (1PE12.5)

      END SUBROUTINE coriol

