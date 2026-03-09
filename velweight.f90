      FUNCTION velweight(nspace, x, y, z, vel)
!************************************************************
!                                                           *
!  This subroutine is called by sphdis.f90                  *
!                                                           *
!************************************************************

      use idims
      use rbnd

      implicit none

      INTEGER(I4B), intent(in) :: nspace
      INTEGER(I4B) :: iposx, iposy, iposz
      REAL(DP) :: velweight
      REAL(DP), intent(in) :: x, y, z, vel(64,64,64)
      REAL(DP) :: deli, delx, dely, delz, &
               velx1, velx2, vely1, vely2

      deli = 1.0/DBLE(nspace/2)

      iposx = INT(x/rmax*(nspace/2)+(nspace/2))
      iposy = INT(y/rmax*(nspace/2)+(nspace/2))
      iposz = INT(z/rmax*(nspace/2)+(nspace/2))

      delx = x - (iposx-(nspace/2))/DBLE(nspace/2)*rmax
      dely = y - (iposy-(nspace/2))/DBLE(nspace/2)*rmax
      delz = z - (iposz-(nspace/2))/DBLE(nspace/2)*rmax
!
!--Find interpolated velocity
!
      velx1 = vel(iposx,iposy,iposz) + delx/deli* &
           (vel(iposx+1,iposy,iposz)-vel(iposx,iposy,iposz))
      velx2 = vel(iposx,iposy+1,iposz) + delx/deli* &
           (vel(iposx+1,iposy+1,iposz)-vel(iposx,iposy+1,iposz))
      vely1 = velx1 + dely/deli*(velx2-velx1)

      velx1 = vel(iposx,iposy,iposz+1) + delx/deli* &
           (vel(iposx+1,iposy,iposz+1)-vel(iposx,iposy,iposz+1))
      velx2 = vel(iposx,iposy+1,iposz+1) + delx/deli* &
           (vel(iposx+1,iposy+1,iposz+1)-vel(iposx,iposy+1,iposz+1))
      vely2 = velx1 + dely/deli*(velx2-velx1)

      velweight = vely1 + delz/deli*(vely2-vely1)

      END FUNCTION velweight
