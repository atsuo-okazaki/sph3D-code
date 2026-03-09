      MODULE logun
!************************************************************
!                                                           *
!  This routine attributes logical unit numbers.            *
!                                                           *
!************************************************************
! iptprint: P.....
! iaccpr: A......
! ikillpr: K....
! ireasspr: R....
!!! icaptpr: C.....
! icoolpr: C.....
! ipgpr: PG....

      use idims
      implicit none
      save

!      INTEGER(I4B), parameter :: iprint = 6, &
      INTEGER(I4B), parameter :: iprint = 9, &
                            iterm = 50, &
                            idisk1 = 11, &
                            idisk2 = 12, &
                            idisk3 = 13, &
                            iptprint = 14, &
                            iaccpr = 15, &
                            ikillpr = 16, &
                            ireasspr = 17, &
                            inotify = 18, &
                            ipgpr = 20, &
                            icoolpr = 21
!!                            icaptpr = 21
      INTEGER(I4B), parameter :: imaxrec = 200 * idim

      END MODULE logun
