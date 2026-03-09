      SUBROUTINE scaling(gt, rs, drdt, dlnrdt)
!************************************************************
!                                                           *
!  Subroutine to compute the scaling factor and its various *
!     derivatives.                                          *
!                                                           *
!************************************************************

      use idims

      use expan

      implicit none

      REAL(DP) :: gt, rs, drdt, dlnrdt

      rs = 1.0 + vexpan*gt
      drdt = vexpan
      dlnrdt = drdt/rs

      END SUBROUTINE scaling
