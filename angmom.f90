      SUBROUTINE angmom
!************************************************************
!                                                           *
!  Computes the total angular momentum of the system        *
!                                                           *
!************************************************************

      use idims

      use constants
      use part
      use bodys
      use carac
      use angm
      use kerne
      use logun
      use debug
      use ptmass
      use phase

!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine angmom')
!
!--Initialisation
!
      angx = 0.
      angy = 0.
      angz = 0.
!
!--Compute total angular momentum
!
      DO i = 1, npart
         IF (iphase(i).GE.0) THEN
            angx = angx + pmass(i)*(y(i)*vz(i) - vy(i)*z(i))
            angy = angy + pmass(i)*(vx(i)*z(i) - x(i)*vz(i))
            angz = angz + pmass(i)*(x(i)*vy(i) - vx(i)*y(i))
         ENDIF
      END DO
!
!--Add spin angular momentum of point masses
!
      DO i = 1, nptmass
         angx = angx + spinx(i)
         angy = angy + spiny(i)
         angz = angz + spinz(i)
      END DO

      angto = SQRT(angx**2 + angy**2 + angz**2)

      IF (idebug(1:6).EQ.'angmom') THEN
         WRITE (iprint, 99002) angx, angy, angz, angto
99002    FORMAT (1X, 4(1PE12.5,1X))
      ENDIF

      END SUBROUTINE angmom
