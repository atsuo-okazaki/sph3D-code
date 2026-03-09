      MODULE accnum

      use idims
      implicit none
      save

      INTEGER(I4B) :: nactotal(iptdim), nghtotal(iptdim), &
                 nactotx(iptdim)
      REAL(DP) :: ptmassinner(iptdim), ptmassinx(iptdim)

      END MODULE accnum
