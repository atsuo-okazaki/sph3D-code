      subroutine dsortx(x, incx, n, indx)
      implicit none
      integer :: incx, n
      integer :: indx(n)
      real(8) :: x(n)

      integer(2), external :: compar
      type :: sort_t
        real(8) :: val
        integer(8) :: idx
      end type
      type(sort_t) :: wk(n)
      integer :: i
      integer(8) :: nmemb, isize

      if (incx .ne. 1) then
        stop 'dsortx'
      end if

      do i = 1, n
        wk(i)%val = x(i)
        wk(i)%idx = i
      end do

      nmemb = n
      isize = sizeof(wk(1))
      call qsort(wk, nmemb, isize, compar)

      do i = 1, n
        x(i) = wk(i)%val
        indx(i) = wk(i)%idx
      end do
      end

      integer(2) function compar(a, b)
      implicit none
      real(8) :: a, b

      if (a .gt. b) then
        compar = 1
      else if (a .lt. b) then
        compar = -1
      else
        compar = 0
      end if
      end
