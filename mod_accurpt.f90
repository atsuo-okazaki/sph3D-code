      MODULE accurpt

      use idims
      implicit none
      save

      real(DP) :: ptmsyn(iptdim), ptmadd(iptdim), &
                     xmomsyn(iptdim), ymomsyn(iptdim), &
                     zmomsyn(iptdim), xmomadd(iptdim), &
                     ymomadd(iptdim), zmomadd(iptdim)

      END MODULE accurpt
