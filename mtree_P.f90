      SUBROUTINE mtree(nnatom, npart, xx, yy, zz, emm, h)
!************************************************************
!                                                           *
!  Subroutine by W. Press (11/21/86).  Given the positions  *
!     and masses of NNATOM points, this subroutine          *
!     constructs a nearest neighbor tree (including         *
!     quadrupole moments) for use by TREEFORCE, filling EM, *
!     R, Q, ISIB, IPAR, IDAU, and NROOT, the pointer to the *
!     root of the tree (which is also the largest value     *
!     used in the range 1..MMAX)                            *
!                                                           *
!************************************************************

      use idims

!---- treecom.mod <-- treecom.f90 for non-parallel simulation
!---- treecom.mod <-- treecom_P.f90 for parallel simulation
      use treecom
      use logun
      use debug
      use phase
      use indexx_mod

      implicit none

      INTEGER(I4B) :: nnatom, npart
      INTEGER(I4B), parameter :: ncbrt=150
      INTEGER(I4B) :: next1(mmax), next2(mmax), next3(mmax)
      INTEGER(I4B) :: iused(idim)
      INTEGER(I4B) :: i, j, k, l, m, n, nactive, nactold, newend, ipart, &
               icbrt, icbrt1, nhash, ibin, itemp, ibin1, &
               nj,  np, mm, ix, iy, iz, nc, &
               llx, lly, llz, lllx, llly, lllz, lux, luy, luz, &
               llux, lluy, lluz, jx, jy, jz, icjx, icjy, icjz, nwal, &
               new, newold, ll, nglob
      REAL(DP), parameter :: bigno=1.0d30, hashfac=1.0d0
!!      REAL(DP) :: xx(nnatom), yy(nnatom), zz(nnatom), &
!!               emm(nnatom), h(idim)
      REAL(DP) :: xx(idim), yy(idim), zz(idim), emm(idim), h(idim)
      REAL(DP) :: xmap(ncbrt), ymap(ncbrt), zmap(ncbrt)
      REAL(DP) :: third, ddd, xnp, ynp, znp, d, dmin, fl, fll, &
               emred, difx, dify, difz, rr
      CHARACTER(len=7) :: where='mtree'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine mtree')

      third = 1./3.

!!      WRITE (iprint,*) 'nnatom=',nnatom,', npart=',npart
!!      WRITE (iprint,*) 'xx(1)=',xx(1),', yy(1)=',yy(1),', zz(1)=',zz(1)

      natom = nnatom
      nactatom = 0

      DO j = 1, natom
         IF (iphase(j).EQ.0 .OR. (iphase(j).GE.1 .AND. iptintree)) THEN
            nactatom = nactatom + 1
            listmap(nactatom) = j
         ENDIF
      END DO

      nactive = nactatom
      nactold = 0
      newend = natom
      nlevel = 0

!$OMP PARALLEL default(none) &
!$OMP shared(natom,npart,list,rx,ry,rz,em,xx,yy,zz,emm,h,isib,idau) &
!$OMP shared(qrad,qxx,qyy,qzz,qxy,qyz,qzx,hnode,listmap) &
!$OMP private(j,ipart) &
!$OMP shared(nactatom,nactive,newend,where,iused) &
!$OMP shared(ipar,nay) &
!$OMP shared(ihash,nhash,next,key,nglob) &
!$OMP private(l,new,n,ll,fl,fll,emred,difx,dify,difz,rr) &
!$OMP shared(icbrt,icbrt1,xmap,ymap,zmap,level) &
!$OMP shared(next1,next2,next3,third,nactold,newold,nlevel) &
!$OMP private(m,np,mm,iz,iy,ix,nc,ddd,xnp,ynp,znp) &
!$OMP private(lllx,llux,llly,lluy,lllz,lluz,llx,lux,lly,luy,llz,luz) &
!$OMP private(jx,icjx,jy,icjy,jz,k,d,dmin,nwal) &
!$OMP private(ibin,i,itemp,nj,ibin1)

!$OMP DO SCHEDULE(runtime)
      DO j = 1, nactatom
         ipart = listmap(j)
         list(j) = ipart
         rx(ipart) = xx(ipart)
         ry(ipart) = yy(ipart)
         rz(ipart) = zz(ipart)
         IF (ipart.GT.npart) THEN
            em(ipart) = 0.
         ELSE
            em(ipart) = emm(ipart)
         ENDIF
         isib(ipart) = 0
         idau(ipart) = 0
         qrad(ipart) = 0.
         qxx(ipart) = 0.
         qyy(ipart) = 0.
         qzz(ipart) = 0.
         qxy(ipart) = 0.
         qyz(ipart) = 0.
         qzx(ipart) = 0.
         hnode(ipart) = h(ipart)
      END DO
!$OMP END DO
!
!--Return here to fill in each new level of hierarchy:
!     find nearest neighbors:
!
 50   CONTINUE

!************************************************************
!                                                           *
!  Used to be a subroutine called naybor(nactive)           *
!                                                           *
!************************************************************

!$OMP SINGLE
      icbrt = hashfac*DBLE(nactive)**third
      icbrt1 = icbrt - 1
      nhash = icbrt**3
      IF (nhash.GT.nactive) CALL error(where, 1)
      IF (icbrt.GT.ncbrt) CALL error(where, 2)
!$OMP END SINGLE

!$OMP SECTIONS
      CALL indexx(nactive, list, rx, next1)
!$OMP SECTION
      CALL indexx(nactive, list, ry, next2)
!$OMP SECTION
      CALL indexx(nactive, list, rz, next3)
!$OMP END SECTIONS
!
!--Bin the points by X and fill XMAP
!
!$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) xmap(ibin) = 0.5*(rx(list(next1(i))) &
                                    + rx(list(next1(i + 1))))
         itemp = icbrt*(ibin - 1)
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            key(next1(j)) = itemp
         END DO
      END DO
!$OMP END DO
!
!--Bin the points by Y and fill YMAP
!
!$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) ymap(ibin) = 0.5*(ry(list(next2(i))) &
                                   + ry(list(next2(i + 1))))
         ibin1 = ibin - 1
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            nj = next2(j)
            key(nj) = icbrt*(key(nj) + ibin1)
         END DO
      END DO
!$OMP END DO
!
!--Bin the points by Z and fill ZMAP
!
!$OMP DO SCHEDULE(runtime)
      DO ibin = 1, icbrt
         i = (ibin*nactive)/icbrt
         IF (i.NE.nactive) zmap(ibin) = 0.5*(rz(list(next3(i))) &
                                   + rz(list(next3(i + 1))))
         DO j = ((ibin-1)*nactive)/icbrt+1, i
            nj = next3(j)
            key(nj) = key(nj) + ibin
         END DO
      END DO
!$OMP END DO
!
!--Now fill the head-of-list table and the linked list
!
!$OMP DO SCHEDULE(runtime)
      DO m = 1, nhash
         ihash(m) = 0
      END DO
!$OMP END DO

!$OMP SINGLE
      DO j = 1, nactive
         m = key(j)
         next(j) = ihash(m)
         ihash(m) = j
      END DO
!$OMP END SINGLE
!
!--Now loop over the particles to find nearest neighbors
!
!$OMP DO SCHEDULE(runtime)
      DO np = 1, nactive
!
!--Find its cell
!
         mm = key(np) - 1
         iz = MOD(mm, icbrt)
         mm = mm/icbrt
         iy = MOD(mm, icbrt)
         ix = mm/icbrt
!
!--Set closest distance so far and coordinates
!
         nc = 0
         ddd = bigno
         xnp = rx(list(np))
         ynp = ry(list(np))
         znp = rz(list(np))
!
!--Initialize the limits of the search loop
!
         lllx = ix
         llux = ix
         llly = iy
         lluy = iy
         lllz = iz
         lluz = iz
         llx = ix
         lux = ix
         lly = iy
         luy = iy
         llz = iz
         luz = iz
!
!--Loop for finding closest particle in the volume of new cells
!
 100     DO jx = llx, lux
            icjx = icbrt*jx
            DO jy = lly, luy
               icjy = icbrt*(jy + icjx) + 1
               DO jz = llz, luz
                  m = icjy + jz
                  k = ihash(m)
                  IF (k.NE.0) THEN
 200                 IF (k.NE.np) THEN
                        d = (rx(list(k)) -xnp)**2 +(ry(list(k)) -ynp) &
                            **2 + (rz(list(k)) - znp)**2
                        IF (d.LT.ddd) THEN
                           ddd = d
                           nc = k
                        ENDIF
                     ENDIF
                     k = next(k)
                     IF (k.NE.0) GOTO 200
                  ENDIF
               END DO
            END DO
         END DO
!
!--We must now find the closest wall, if any
!
         dmin = bigno
         nwal = 0
         IF (lllx.GT.0) THEN
            d = xnp - xmap(lllx)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 1
               ENDIF
            ENDIF
         ENDIF
         IF (llux.LT.icbrt1) THEN
            d = xmap(llux + 1) - xnp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 2
               ENDIF
            ENDIF
         ENDIF
         IF (llly.GT.0) THEN
            d = ynp - ymap(llly)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 3
               ENDIF
            ENDIF
         ENDIF
         IF (lluy.LT.icbrt1) THEN
            d = ymap(lluy + 1) - ynp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 4
               ENDIF
            ENDIF
         ENDIF
         IF (lllz.GT.0) THEN
            d = znp - zmap(lllz)
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 5
               ENDIF
            ENDIF
         ENDIF
         IF (lluz.LT.icbrt1) THEN
            d = zmap(lluz + 1) - znp
            IF (d**2.LT.ddd) THEN
               IF (d.LT.dmin) THEN
                  dmin = d
                  nwal = 6
               ENDIF
            ENDIF
         ENDIF
!
!--Reset the search volume to its augmented value
!
         llx = lllx
         lux = llux
         lly = llly
         luy = lluy
         llz = lllz
         luz = lluz
!
!--Augment it and set the next search according to which wall is closest
!
         IF (nwal.NE.0) THEN
            IF (nwal.EQ.1) THEN
               lllx = lllx - 1
               llx = lllx
               lux = lllx
            ELSEIF (nwal.EQ.2) THEN
               llux = llux + 1
               llx = llux
               lux = llux
            ELSEIF (nwal.EQ.3) THEN
               llly = llly - 1
               lly = llly
               luy = llly
            ELSEIF (nwal.EQ.4) THEN
               lluy = lluy + 1
               lly = lluy
               luy = lluy
            ELSEIF (nwal.EQ.5) THEN
               lllz = lllz - 1
               llz = lllz
               luz = lllz
            ELSEIF (nwal.EQ.6) THEN
               lluz = lluz + 1
               llz = lluz
               luz = lluz
            ENDIF
            GOTO 100
         ENDIF
         nay(np) = nc
      END DO
!$OMP END DO

!************************************************************
!                                                           *
!  End of old subroutine called naybor(nactive)             *
!                                                           *
!************************************************************

!$OMP SINGLE
      IF (nactive.EQ.nactold) THEN
         WRITE (iprint,99100) nactive
99100    FORMAT (' MTREE IS IN AN INFINITE LOOP ',I8)
         DO k = 1, nactive
            l = list(k)
            WRITE (iprint,99101) l, rx(l), ry(l), rz(l)
99101       FORMAT (I8, 1F12.7, 1F12.7, 1F12.7)
         END DO
         CALL quit
      ENDIF
!$OMP END SINGLE
!
!--Find new nodes:
!
!$OMP DO SCHEDULE(runtime)
      DO j = 1, nactive
         l = list(j)
         new = newend + j
         IF (new.GT.mmax) CALL error(where, 3)
         iused(j) = 0
         n = nay(j)
         IF (isib(l).EQ.0 .AND. n.GT.j) THEN
            IF (nay(n).EQ.j) THEN
               iused(j) = 1
               ll = list(n)
               isib(ll) = l
               isib(l) = ll
               ipar(ll) = new
               ipar(l) = new
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
               rx(new) = rx(l) + fll*difx
               ry(new) = ry(l) + fll*dify
               rz(new) = rz(l) + fll*difz
!
!--Find radius and maximum h
!
               rr = SQRT(difx**2 + dify**2 + difz**2)
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
               qzz(new) = emred*(difz*difz)
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
               isib(new) = 0
               idau(new) = l
            ENDIF
         ENDIF
      END DO
!$OMP END DO
!
!--Compactify list:
!
!$OMP SINGLE
      newold = newend
      nglob = 1
      DO j = 1, nactive
         IF (iused(j).EQ.1) THEN
            new = newold + j
            newend = newend + 1
            em(newend) = em(new)
            rx(newend) = rx(new)
            ry(newend) = ry(new)
            rz(newend) = rz(new)
            qrad(newend) = qrad(new)
            hnode(newend) = hnode(new)
            qxx(newend) = qxx(new)
            qyy(newend) = qyy(new)
            qzz(newend) = qzz(new)
            qxy(newend) = qxy(new)
            qyz(newend) = qyz(new)
            qzx(newend) = qzx(new)
            isib(newend) = isib(new)
            idau(newend) = idau(new)
            ipar(idau(newend)) = newend
            ipar(isib(idau(newend))) = newend
         ELSEIF (isib(list(j)).EQ.0) THEN
            list(nglob) = list(j)
            nglob = nglob + 1
         ENDIF
      END DO
!$OMP END SINGLE

!$OMP DO SCHEDULE(runtime)
      DO j = newold + 1, newend
         list(nglob+j-(newold + 1)) = j
      END DO
!$OMP END DO

!$OMP SINGLE
      nglob = nglob + (newend - newold)

      nlevel = nlevel + 1
      IF (nlevel.GT.nmaxlevel) CALL error(where,4)
      level(nlevel) = newold + 1
      nactold = nactive
      nactive = nglob - 1
!$OMP END SINGLE

      IF (nactive.GT.1) GOTO 50

!$OMP END PARALLEL

      nlevel = nlevel + 1
      IF (nlevel.GT.nmaxlevel) CALL error(where,4)
      level(nlevel) = newend + 1

      nroot = newend

      END SUBROUTINE mtree
