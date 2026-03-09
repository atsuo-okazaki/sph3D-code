      MODULE force

      use idims
      implicit none
      save

!---- force
      REAL(DP) :: fx(idim), fy(idim), fz(idim), du(idim), dh(idim)
!---- f1
      REAL(DP) :: f1vx(idim), f1vy(idim), f1vz(idim), f1u(idim), &
                     f1h(idim)
!---- f2
      REAL(DP) :: f2vx(idim), f2vy(idim), f2vz(idim), f2u(idim), &
                     f2h(idim)

      END MODULE force
