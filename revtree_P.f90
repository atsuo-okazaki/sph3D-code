      SUBROUTINE revtree (nnatom, npart, xx, yy, zz, emm, h)
!************************************************************
!                                                           *
!  Subroutine by W. Press.  Revises the tree structure.     *
!     Assumes that only the positions have been altered     *
!                                                           *
!************************************************************

      use idims

!---- treecom.mod <-- treecom.f90 for non-parallel simulation
!---- treecom.mod <-- treecom_P.f90 for parallel simulation
      use treecom

      implicit none

      INTEGER(I4B) :: nnatom, npart, j, ipart, ilevel, new, l, ll
!!      REAL(DP) :: xx(nnatom), yy(nnatom), zz(nnatom), &
!!                        emm(nnatom), h(nnatom)
      REAL(DP) :: xx(idim), yy(idim), zz(idim), emm(idim), h(idim)
      REAL(DP) :: fl, fll, emred, difx, dify, difz, rr

      natom = nnatom

!$OMP PARALLEL default(none) &
!$OMP shared(natom,npart,rx,ry,rz,em,xx,yy,zz,emm,h,hnode) &
!$OMP shared(nactatom,listmap) &
!$OMP private(j,ipart)

!$OMP DO SCHEDULE(runtime)
      DO j = 1, nactatom
         ipart = listmap(j)
         rx(ipart) = xx(ipart)
         ry(ipart) = yy(ipart)
         rz(ipart) = zz(ipart)
         hnode(ipart) = h(ipart)

         IF (ipart.LE.npart) THEN
            em(ipart) = emm(ipart)
         ELSE
            em(ipart) = 0.
         ENDIF
      END DO
!$OMP END DO
!$OMP END PARALLEL

      DO ilevel = 1, nlevel

!$OMP PARALLEL default(none) &
!$OMP shared(ilevel,level,natom,rx,ry,rz,em,isib,idau) &
!$OMP shared(qrad,qxx,qyy,qzz,qxy,qyz,qzx,hnode) &
!$OMP private(new,l,ll,fl,fll,emred,difx,dify,difz,rr)

!$OMP DO SCHEDULE(runtime)
         DO new = level(ilevel), level(ilevel + 1) - 1

            l = idau(new)
            ll = isib(l)

            em(new) = em(l) + em(ll)

            IF (em(new).NE.0) THEN
               fl = em(l)/em(new)
               fll = em(ll)/em(new)
            ELSE
               fl = 0.5
               fll = 0.5
            ENDIF
            emred = fl*fll*em(new)
            difx = rx(ll) - rx(l)
            dify = ry(ll) - ry(l)
            difz = rz(ll) - rz(l)
            IF (fl.GT.fll) THEN
               rx(new) = rx(l) + fll*difx
               ry(new) = ry(l) + fll*dify
               rz(new) = rz(l) + fll*difz
            ELSE
               rx(new) = rx(ll) - fl*difx
               ry(new) = ry(ll) - fl*dify
               rz(new) = rz(ll) - fl*difz
            ENDIF
!
!--Find radius
!
            rr = SQRT(difx**2 + dify**2 + difz**2) + tiny
            qrad(new) = MAX(fll*rr + qrad(l), fl*rr + qrad(ll))
            hnode(new) = MAX(hnode(l), hnode(ll))
!
!--Find quadrupole moments
!
            qxx(new) = (emred*difx)*difx
            qxy(new) = (emred*difx)*dify
            qzx(new) = (emred*difx)*difz
            qyy(new) = (emred*dify)*dify
            qyz(new) = (emred*dify)*difz
            qzz(new) = (emred*difz)*difz
            IF (l.GT.natom) THEN
               qxx(new) = qxx(new) + qxx(l)
               qyy(new) = qyy(new) + qyy(l)
               qzz(new) = qzz(new) + qzz(l)
               qxy(new) = qxy(new) + qxy(l)
               qyz(new) = qyz(new) + qyz(l)
               qzx(new) = qzx(new) + qzx(l)
            ENDIF
            IF (ll.GT.natom) THEN
               qxx(new) = qxx(new) + qxx(ll)
               qyy(new) = qyy(new) + qyy(ll)
               qzz(new) = qzz(new) + qzz(ll)
               qxy(new) = qxy(new) + qxy(ll)
               qyz(new) = qyz(new) + qyz(ll)
               qzx(new) = qzx(new) + qzx(ll)
            ENDIF
         END DO
!$OMP END DO
!$OMP END PARALLEL

      END DO

      END SUBROUTINE revtree
