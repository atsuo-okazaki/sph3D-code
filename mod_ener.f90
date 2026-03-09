      MODULE ener

      use idims
      implicit none
      save

      REAL(DP) :: poten(idim), dq(idim), pdv(idim)
      REAL(DP) :: trotz, trotx, tkin, tgrav, tterm
      REAL(DP) :: dphit(idim), dgrav(idim)

      END MODULE ener
