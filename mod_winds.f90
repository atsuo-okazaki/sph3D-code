      MODULE winds

      use idims
      implicit none
      save

!---- vbeta1, vbeta2, akappac,a nd akappar are parameters for
!     the beta-law stellar winds
      INTEGER(I4B) :: nshell1(3), iantigr(idim), ibelong(idim)
      REAL(DP) :: vwind1, vwind2, rshell1, rshell2, &
              vrot1, vrot2, vinf1, vinf2, therm1, therm2, &
              emdotratio, partmass(3), vbeta(2)

      END MODULE winds
