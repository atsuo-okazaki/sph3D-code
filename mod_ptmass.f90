      MODULE ptmass

      use idims
      implicit none
      save

      INTEGER(I4B) :: nptmass, listpm(iptdim), iptmass
      REAL(DP) :: spinx(iptdim), spiny(iptdim), spinz(iptdim), &
                     hacc, haccall, ptmcrit, radcrit, &
                     angaddx(iptdim), angaddy(iptdim), angaddz(iptdim), &
                     spinadx(iptdim), spinady(iptdim), spinadz(iptdim), &
                     rptmas(iptdim)

      END MODULE ptmass
