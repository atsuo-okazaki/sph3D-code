      MODULE numpa

      use idims
      implicit none
      save

!---- numpa
      REAL(DP) :: alpha, beta 
!---- Per-particle artificial-viscosity coefficients for ifsvi=6
      REAL(DP) :: alphasph(idim)=0.0_DP, betasph(idim)=0.0_DP

      END MODULE numpa
