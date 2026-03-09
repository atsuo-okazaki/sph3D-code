      MODULE table

      use idims
      implicit none
      save

      INTEGER(I4B), parameter :: itable=40000
      REAL(DP) :: wij(0:itable), grwij(0:itable), fmass(0:itable), &
                     fpoten(0:itable), dphidh(0:itable), dvtable

      END MODULE table
