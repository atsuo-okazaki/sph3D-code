      MODULE torq

      use idims
      implicit none
      save

      REAL(DP) :: torqt(idim), torqg(idim), torqp(idim), &
                     torqv(idim), torqc(idim)

      END MODULE torq
