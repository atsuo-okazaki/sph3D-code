      MODULE ghost

      use idims
      implicit none
      save

      INTEGER(I4B), parameter :: ighost = idim
      INTEGER(I4B) :: nghost, ireal(ighost), hasghost(idim)

      END MODULE ghost
