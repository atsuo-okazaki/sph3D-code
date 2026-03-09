      SUBROUTINE energ(ipart, ti, u)
!************************************************************
!                                                           *
!  This routine computes the change in internal energy or   *
!     in entropy.                                           *
!                                                           *
!    varsta = entropy : compute change in specific entropy  *
!    varsta = intener : compute change in specific internal *
!                          energy                           *
!                                                           *
!************************************************************

      use idims

      use constants
      use densi
      use eosq
      use force, only: fx, fy, fz, du, dh
      use cgas
      use ener
      use kerne
      use typef
      use logun
      use varet
      use units
      use phase

      implicit none

      INTEGER(I4B) :: ipart
      REAL(DP) :: u(idim), ti, chi, cnormk05, rscale, drdt, dlnrdt, &
               dqexp, gamma1

      IF (iphase(ipart).NE.0) THEN
            du(ipart) = 0.0
      ELSE
!--      Initialisation
!
         cnormk05 = cnormk*0.5
         gamma1 = gamma - 1.
!
!--      dq from expansion
!
         chi = 4.0 - 3.0*gamma
         CALL scaling(ti, rscale, drdt, dlnrdt)
         dqexp = chi*dlnrdt

!----    versta = 'intener' or 'entropy'
         IF (varsta.NE.'entropy') THEN
!
!--      Compute change in specific internal energy
!
!        a) pdv term first
!
            IF (iexpan.EQ.0) THEN
               du(ipart) = cnormk*pdv(ipart)*pr(ipart)/rho(ipart)**2
!c               WRITE (iprint,*) 'pdv=',pdv(ipart),', pr=',pr(ipart),
!c     &          ', rho=',rho(ipart)
!c               WRITE (iprint,*) 'du(pdV(',ipart,'))=',du(ipart)
            ELSE
              du(ipart) = cnormk*pdv(ipart)*pr(ipart)/rho(ipart)**2 &
                              - u(ipart)*dqexp
            ENDIF
!
!        b) Shock dissipation if appropriate
!
            IF (ichoc.NE.0) THEN
               du(ipart) = du(ipart) + cnormk05*dq(ipart)
!c              WRITE (iprint,*) 'du(shock(',ipart,'))=',
!c     &              cnormk05*dq(ipart)
            ENDIF

         ELSE
!
!--         Compute change in specific entropy
!
!           a) Shock dissipation
!
               IF (ichoc.NE.0) THEN
                  du(ipart) = u(ipart)*gamma1*dq(ipart)*rho(ipart)* &
                       cnormk05/pr(ipart) - pr(ipart)*dqexp
               ENDIF
         ENDIF
      ENDIF

      END SUBROUTINE energ
