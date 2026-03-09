      SUBROUTINE eospg(ipart,u,rho,pr,vsound)
!************************************************************
!                                                           *
!  This routine computes the pressure and sound speed       *
!     according to a perfect gas equation of state          *
!                                                           *
!************************************************************

      use idims

      use cgas
      use varet
      use polyk2
      use units
      use vargam
      use typef
      use logun
      use physeos
      use winds

      implicit none

      INTEGER(I4B) :: ipart
      REAL(DP) :: u(idim), rho(idim), pr(idim), vsound(idim)
      REAL(DP) :: gama1, rhocrt, rhocrt2, rhocrt3, gam1, gamdh1, &
                 gamah1, gam2, gam3, uzero, RK2, RK3, RK4, RK5
      CHARACTER(len=7) :: where='eospg'

      gama1 = gamma - 1.
      rhocrt = rhocrit * udens
      rhocrt2 = rhocrit2 * udens
      rhocrt3 = rhocrit3 * udens
!
!--Variable is internal energy
      IF (ibound >= 94 .AND. ibound <= 97) THEN
         IF (ibelong(ipart) == 1) THEn
            IF (iantigr(ipart) == 0) THEN
                RK2 = RK23
            ELSE
                RK2 = RK21
            END IF
         ELSE
            RK2 = RK22
         END IF
      ELSE
         RK2 = RK21
      END IF
!
!--Variable gamma equation of state
!
      IF (encal.EQ.'v') THEN
          gam1 = gam - 1.
          gamdh1 = gamdh - 1.
          gamah1 = gamah - 1.
          uzero = RK2 * (rhozero ** gama1)
          RK3 = uzero / (rhocrit ** gam1)
          RK4 = RK3 * (rhocrit2 ** gam1) / (rhocrit2 ** gamdh1)
          RK5 = RK4 * (rhocrit3 ** gamdh1) / (rhocrit3 ** gamah1)
          IF (rho(ipart).LT.rhocrit) THEN
             u(ipart) = uzero
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(pr(ipart)/rho(ipart))
          ELSE IF (rho(ipart).LT.rhocrit2) THEN
             u(ipart) = RK3 * (rho(ipart) ** gam1)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gam * pr(ipart) / rho(ipart))
          ELSE IF (rho(ipart).LT. rhocrit3) THEN
             u(ipart) = RK4 * (rho(ipart) ** gamdh1)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamdh * pr(ipart) / rho(ipart))
          ELSE
             u(ipart) = RK5 * (rho(ipart) ** gamah1)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamah * pr(ipart) / rho(ipart))
          ENDIF
!
!--Physical equation of state
!
      ELSEIF (encal.EQ.'x') THEN
          gam1 = gamphys1 - 1.
          gam2 = gamphys2 - 1.
          gam3 = gamphys3 - 1.
          uzero = RK2 * (rhozero ** gama1)
          IF (rho(ipart).LT.rhochange1) THEN
             u(ipart) = uzero*(1.0 + (rho(ipart)/rhoref1)**gam1)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys1*pr(ipart)/rho(ipart))
          ELSE IF (rho(ipart).LT.rhochange2) THEN
             u(ipart) = uzero*(1.0 + (rho(ipart)/rhoref2)**gam2)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys2* pr(ipart) / rho(ipart))
          ELSE
             u(ipart) = uzero*(1.0 + (rho(ipart)/rhoref3)**gam3)
             pr(ipart) = 2./3. * u(ipart) * rho(ipart)
             vsound(ipart) = SQRT(gamphys3* pr(ipart) / rho(ipart))
          ENDIF
!
!--Isothermal equation of state
!
      ELSEIF (encal.EQ.'i') THEN
            pr(ipart) = (2./3.) * u(ipart) * rho(ipart)
            vsound(ipart) = SQRT(pr(ipart)/rho(ipart))
!
!--Adiabatic equation of state
!
      ELSEIF (encal.EQ.'a') THEN
            pr(ipart) =  gama1 * u(ipart) * rho(ipart)
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
!
!--Polytropic equation of state
!
      ELSEIF (encal.EQ.'p') THEN
            u(ipart) = RK2 * (rho(ipart) ** gama1)
            pr(ipart) = (2./3.) * u(ipart) * rho(ipart)
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
!
!--Gas with radiative cooling
!
      ELSEIF (encal.EQ.'c') THEN
            pr(ipart) =  gama1 * u(ipart) * rho(ipart)
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
!
!--Variable of state is entropy
!
      ELSEIF (varsta.EQ.'entropy') THEN
            pr(ipart) = u(ipart)*rho(ipart)**gamma
            vsound(ipart) = SQRT(gamma*pr(ipart)/rho(ipart))
      ELSE
         WRITE (iprint,99500) encal
99500    FORMAT ('encal = ',A1)
         CALL error(where,1)
      ENDIF

      END SUBROUTINE eospg
