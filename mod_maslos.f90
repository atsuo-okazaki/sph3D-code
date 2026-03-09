      MODULE maslos

      use idims
      implicit none
      save

      INTEGER(I4B), parameter :: MAXPHASEI=1024
      INTEGER(I4B) :: nphasei, iphsi
      REAL(DP) :: phasei(MAXPHASEI), eninfl(MAXPHASEI), &
                     xinfl(MAXPHASEI), yinfl(MAXPHASEI), &
                     zinfl(MAXPHASEI), vxinfl(MAXPHASEI), &
                     vyinfl(MAXPHASEI), vzinfl(MAXPHASEI), &
                     sxinfl(MAXPHASEI), syinfl(MAXPHASEI), &
                     szinfl(MAXPHASEI), svxinfl(MAXPHASEI), &
                     svyinfl(MAXPHASEI), svzinfl(MAXPHASEI), &
                     hinfl(MAXPHASEI), shinfl(MAXPHASEI), &
                     phase0, emdot0, sinj0(3), sinj(3), &
                     dphasei

      END MODULE maslos
