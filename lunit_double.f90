      SUBROUTINE lunit
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

      use logun
      use ptdump
      use binfile

      implicit none

      iprint = 9
!      iprint = 6
      iterm = 50
      idisk1 = 11
      idisk2 = 12
      idisk3 = 13
      iptprint = 14
      iaccpr = 15
      ikillpr = 16
      ireasspr = 17
      inotify = 18
      ipgpr = 20
!!      icaptpr = 21
      icoolpr = 21

      imaxrec = 200*idim

      END SUBROUTINE lunit
