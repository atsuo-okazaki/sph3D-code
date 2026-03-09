      MODULE pres

      use idims
      implicit none
      save

!---- pres
      REAL(DP) :: gradpx(idim), gradpy(idim), gradpz(idim)
!---- presb
      REAL(DP) :: pext

      END MODULE pres
