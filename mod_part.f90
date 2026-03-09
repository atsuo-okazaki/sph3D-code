      MODULE part

      use idims
      implicit none
      save

      INTEGER(I4B) :: npart
      REAL(DP) :: x(idim), y(idim), z(idim), vx(idim), &
                     vy(idim), vz(idim), u(idim), h(idim)

      END MODULE part
