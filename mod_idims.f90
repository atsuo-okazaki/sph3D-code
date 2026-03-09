      MODULE idims

      implicit none
      save

      INTEGER, parameter :: SP = SELECTED_REAL_KIND(p=6,r=37)
      INTEGER, parameter :: DP = SELECTED_REAL_KIND(p=13,r=200)
      INTEGER, parameter :: I4B = selected_int_kind(r=9)

      INTEGER(I4B), parameter :: idim=1200000,    &
                            iptdim=2,       &
                            iptneigh=200000, &
                            nlmax = 50000

      REAL(DP), parameter :: tiny=1.0E-18

!---- igrape
!---- The parameter `igrape' tells the program whether it is 
!     finding the forces and neighbours of particles using 
!     the BINARY TREE (igrape=0) or the GRAPE board (igrape=1).

      INTEGER(I4B), parameter :: igrape=0
!      INTEGER(I4B), parameter :: igrape=1

      INTEGER(I4B), parameter :: nmaxboards = 1,      &
                            nmaxchips = 8,       &
                            nmaxneighbours = 1024

!---- The parameter `isoft' tells the program whether to soften the 
!     gravitational forces using kernel smoothing of the potential
!     (isoft=0), or using the 1/(r+e) smoothing used by 
!     the GRAPE (isoft=1).
!
!     NOTE:  This only has an effect when the BINARY TREE is used.
!            When igrape=1, the softening is ALWAYS done as 1/(r+e).

      INTEGER(I4B), parameter :: isoft=0
!      INTEGER(I4B), parameter :: isoft=1

!---- Allow or dis-allow timing to be done

!      LOGICAL, parameter :: itiming=.FALSE.
      LOGICAL, parameter :: itiming=.TRUE.

!---- Specify that sink particles are done IN or OUT of the tree

      LOGICAL, parameter :: iptintree=.FALSE.
!      LOGICAL, parameter :: iptintree=.TRUE.

      END MODULE idims
