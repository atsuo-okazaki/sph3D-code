      FUNCTION ran1(iseed)

!************************************************************
!                                                           *
!  Uniform random deviate (from Numerical Recipes 2)        *
!  This rountine must be called with a negative seed        *
!     and then a POSITIVE number from then on               *
!                                                           *
!************************************************************

      use idims

      use savernd

      implicit none

      INTEGER(I4B), intent(in) :: iseed
      REAL(DP) :: ran1
      INTEGER(I4B), parameter :: IA=16807, IM=2147483647, &
           IQ=127773, IR=2836, NTAB=32
      INTEGER(I4B), parameter :: NDIV=1+(IM-1)/NTAB
      REAL(DP), parameter :: AM=1./IM, EPS=1.2e-7
      REAL(DP), parameter :: RNMX=1.-EPS

      INTEGER(I4B) :: j,k
      INTEGER(I4B), dimension(NTAB) :: iv=0
      INTEGER(I4B) :: iy=0
!      DATA iv /NTAB*0/, iy /0/

      IF (iseed.LE.0) THEN
         idum = iseed
      ENDIF

      IF (idum.LE.0 .OR. iy.EQ.0) THEN
         idum=MAX(-idum, 1)
         DO j = NTAB + 8, 1, -1
            k = idum/IQ
            idum = IA*(idum - k*IQ) - IR*k
            IF (idum.LT.0) idum = idum + IM
            IF (j.LE.NTAB) iv(j) = idum
         END DO
         iy = iv(1)
      ENDIF
      k = idum/IQ
      idum = IA*(idum - k*IQ) - IR*k
      IF (idum.LT.0) idum = idum + IM
      j = 1 + iy/NDIV
      iy = iv(j)
      iv(j) = idum
      ran1 = MIN(AM*iy, RNMX)

      END FUNCTION ran1
