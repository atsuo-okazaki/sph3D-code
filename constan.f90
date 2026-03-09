       SUBROUTINE constan
!************************************************************
!                                                           *
!  This subroutine initializes the mathematical, physical   *
!     and astronomical constants (refs: Handbook of         *
!     Chemistry and Physics, 1975-1976 edition;             *
!     Astrophysical Concepts)                               *
!                                                           *
!************************************************************

      use constants

      implicit none

      REAL(DP) :: solarl
!
!--Mathematical constants
!
      pi = 3.141592654
!
!--Physical constants (in cgs units)
!
      c = 2.997924e10
      gg = 6.672041e-8
      Rg = 8.314e7
!
!--Astronomical constants (in cgs units)
!
!--Solar mass and radius
!
      solarm = 1.991e33
      solarr = 6.959500e10
      solarl = 3.85e33
!
!--Earth mass and radius
!
      earthm = 5.979e27
      earthr = 6.371315e8
!
!--Distance scale
!
      au = 1.496e13
      pc = 3.086e18
!
!--Gas molecular weight, mu
!
!!      gmw = 2.0
!      gmw = 2.46
!  For fully ionized gas with the cosmic abundances
      gmw = 0.6

      END SUBROUTINE constan
