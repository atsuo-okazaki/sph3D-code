      MODULE unidis

      use idims
      implicit none
      save

!---- The following parameters are used to inject disk/wind particles
!     uniformly on (a part of) a sphere
!     (implemented by Chris Russell, 16/08/2013)
      INTEGER(I4B), PARAMETER :: idiminj = 3000
      INTEGER(I4B), DIMENSION(idiminj) :: lUD
      INTEGER(I4B) :: nUD, nUD1, nUD2
      REAL(DP), DIMENSION(idiminj) :: xUD, yUD, zUD
      REAL(DP) :: anginj1, anginj2

      END MODULE unidis
