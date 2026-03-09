      MODULE timei

      use idims
      implicit none
      save

      INTEGER(I4B) :: imax, imaxstep, it0(idim), it1(idim), isteps(idim), &
              iteighth
      REAL(DP) :: dum2vx(idim), dum2vy(idim), dum2vz(idim)

      END MODULE timei
