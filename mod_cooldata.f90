      MODULE cooldata

      use idims
      implicit none
      save

      INTEGER(I4B), parameter :: MAXDATA=1500
      INTEGER(I4B) :: ntemp
      REAL(DP), dimension(MAXDATA) :: tempcf,alambda, &
                 alphacf,yfunc
      REAL(DP), dimension(idim) :: tempini,tempfin,tcool

      END MODULE cooldata
