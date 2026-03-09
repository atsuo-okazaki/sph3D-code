      SUBROUTINE homexp(ipart, ti, vx, vy, vz)
!************************************************************
!                                                           *
!  This subroutine computes the correction to the momentum  *
!     equation due to homologous expansion or contraction.  *
!                                                           *
!************************************************************

      use idims

      use force, only: fx, fy, fz, du, dh
      use logun

      implicit none

      INTEGER(I4B) :: ipart
      REAL(DP) :: ti, rscale, rscale3, drdt, dlnrdt, dlnrdt2
      REAL(DP) :: vx(idim), vy(idim), vz(idim)
!
!--Scaling factors
!
      CALL scaling(ti, rscale, drdt, dlnrdt)

      dlnrdt2 = 2.*dlnrdt
      rscale3 = 1./rscale**3
      fx(ipart) = rscale3*fx(ipart) - vx(ipart)*dlnrdt2
      fy(ipart) = rscale3*fy(ipart) - vy(ipart)*dlnrdt2
      fz(ipart) = rscale3*fz(ipart) - vz(ipart)*dlnrdt2

      END SUBROUTINE homexp
