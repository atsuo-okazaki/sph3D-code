      MODULE treecom

      use idims
      implicit none
      save

      INTEGER(I4B), parameter :: mmax=2*idim + 2
      INTEGER(I4B), parameter :: nmaxlevel=1000
      INTEGER(I4B) :: natom, nlevel, nactatom
      INTEGER(I4B) :: key(mmax), next(mmax), ihash(mmax), &
                 list(mmax), nay(mmax), level(nmaxlevel), &
                 listmap(idim), nroot, isib(mmax), ipar(mmax), &
                 idau(mmax)
      REAL(DP) :: em(mmax), rx(mmax), ry(mmax), rz(mmax), &
                 qrad(mmax), qxx(mmax), qyy(mmax), qzz(mmax), &
                 qxy(mmax), qyz(mmax), qzx(mmax)
      REAL(DP) :: hnode(mmax)

      END MODULE treecom
