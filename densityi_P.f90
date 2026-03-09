      SUBROUTINE densityi (jneigh,npart,x,y,z,vx,vy,vz,u,h, &
                  nlst_in,nlst_end,list,itime)
!************************************************************
!                                                           *
!  Subroutine to compute the density and the velocity       *
!     divergence at ipart particle location. This           *
!     subroutine uses the binary tree algorithm to locate   *
!     neighbors.                                            *
!                                                           *
!************************************************************

      use idims
      use mpi_mod

      use constants
      use table
      use tlist
      use carac
      use densi
      use divve
      use eosq
      use kerne
      use typef
      use logun
      use debug
!---- neighbor.mod <-- neighbor.f90 for non-parallel simulation
!---- neighbor.mod <-- neighbor_P.f90 for parallel simulation
      use neighbor
      use current
      use rbnd
      use polyk2
      use phase
      use ptmass
      use nextmpt
      use timei
      use ptdump
      use btree
      use isnpt
      use glrho
      use vbound
      use init
      use debugpt
      use ian
      use visc
      use call
      use xforce
      use winds

      implicit none

      INTEGER(I4B) :: jneigh, npart, nlst_in, nlst_end, itime, ii, &
                 index, index1, indw, indw1, iokay, ipart, ipt, &
                 iptcur, iptden, iptm, iptn, iptnr, j, k, l, ll, n, nl, &
                 nlocalneigh, nmiddle, nsecond, ntanvisc, numneigh1i, &
                 numneigh2i, numneighnew, numneighold
      INTEGER(I4B) :: neighloc(idim), indexs(idim) , numneigh1(idim), &
                 numneigh2(idim), list(idim)
      REAL(DP) :: x(idim), y(idim), z(idim), vx(idim), vy(idim), &
                    vz(idim), h(idim), u(idim)
      REAL(DP) :: curlvx(idim), curlvy(idim), curlvz(idim)
      REAL(DP) :: rr1(idim), temprho(idim)
      REAL(DP) :: akapp(iptdim)
      REAL(DP) :: ck, coordi, curlvxi, curlvyi, curlvzi, ddinv, &
                 densgradrho, dgrwdx, diffr, divvi, dot, dvx, dvy, dvz, &
                 dwdx, dwdxw, dx, dxi, dxj, dxx, dxxw, dy, dyi, dyj, &
                 dz, dzi, dzj, factor, fslope, grwtij, hi, hj, hmean, &
                 hmean2, hmean21, hmean31, hmean41, hi3, pmassi, pmassj, &
                 prmean1, prmean2, procurlvx, procurlvy, procurlvz, &
                 projv, prslope, r2, radcrit2, radi2, radius, rho1dif, &
                 rho2i, rhoi, rhoip2, rhomean1, rhomean2, rhominimum, &
                 rij1, rij2, ript, rjpt, rmean1, rmean2, rr, rr2, &
                 tanv, tanvj, tanvx1j, tanvxi, tanvy1j, tanvyi, &
                 tanvz1j, tanvzi, rho1i, uxip, uxjp, uyip, uyjp, uzip, uzjp, &
                 v2, valuen, vpos, vr, vrj, vsmvfor, vsmvfor2, vxi, vxipt, &
                 vxnl, vyi, vyipt, vynl, vzi, vzipt, vznl, weight, &
                 wtij, xi, xii, xip, xipt, yi, yii, yip, yipt, &
                 zi, zii, zip, zipt
      LOGICAL :: ioutside(idim)

      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT ('entry subroutine densityi')
!
!--Initialise
!
      internum = 0
      rhonext = 0.
      icreate = 0
      radcrit2 = radcrit*radcrit

!$OMP PARALLEL default(none) &
!$OMP shared(nlst_in,nlst_end,list,divv,curlvx,curlvy,curlvz,curlv) &
!$OMP shared(hmax,h,x,y,z,vx,vy,vz,u,pr,vsound,pmass,rho) &
!$OMP shared(ilen,neighb,iscurrent,selfnormkernel) &
!$OMP shared(cnormk,radkernel,dvtable,wij,grwij) &
!$OMP shared(temprho,numneigh1,numneigh2,rho1,rho2) &
!$OMP shared(isnearpt,listpm,iphase,tanvx1,tanvy1,tanvz1) &
!$OMP shared(tanviscfor1,artvix,artviy,artviz,prgradrho,densgrad) &
!$OMP shared(nptmass,iptmass,radcrit2) &
!$OMP shared(vbeta,akappa,iexf,ibelong) &
!$OMP shared(akappac,akappar,rptmas) &
!$OMP private(akapp) &
!$OMP private(iptm,rhoip2,ddinv) &
!$OMP private(n,ipart,j,k,xi,yi,zi,vxi,vyi,vzi,pmassi,hi,hj,rhoi) &
!$OMP private(divvi,curlvxi,curlvyi,curlvzi) &
!$OMP private(numneigh1i,numneigh2i,rho1i,rho2i) &
!$OMP private(iptn,iptcur,xipt,yipt,zipt,vxipt,vyipt,vzipt) &
!$OMP private(xip,yip,zip,nlocalneigh) &
!$OMP private(pmassj,hmean,hmean21,hmean31,hmean41) &
!$OMP private(dx,dy,dz,dvx,dvy,dvz,rij2,rij1,v2) &
!$OMP private(index,dxx,index1,dwdx,wtij,dgrwdx,grwtij) &
!$OMP private(projv,procurlvx,procurlvy,procurlvz) &
!$OMP private(dot,rr2,neighloc,rr1,ii,iptden,xii,yii,zii,dxi,dyi,dzi) &
!$OMP private(dxj,dyj,dzj,vpos,coordi,indexs,ript,uxip,uyip,uzip) &
!$OMP private(vr,tanvxi,tanvyi,tanvzi,tanv,ntanvisc,nmiddle,rmean1) &
!$OMP private(rhomean1,prmean1,ll,nl,rjpt,uxjp,uyjp,uzjp) &
!$OMP private(vxnl,vynl,vznl,vrj,tanvx1j,tanvy1j,tanvz1j,tanvj) &
!$OMP private(factor,rmean2,rhomean2,prmean2,nsecond,diffr) &
!$OMP private(fslope,prslope,hi3,rhominimum,radius,l,radi2) &
!$OMP reduction(+:internum) &
!$OMP reduction(MAX:rhonext)

!$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)
!
!--This counts the number of interactions between 2 particles.  This is
!    necessary, rather than just counting the neighbours of each particle,
!    because, with the GRAPE, some neighbours of `ipart' may be found by
!    a particle `j' and be done when particle `j's neighbours are looped
!    over, but these particles may not appear in `ipart's list of neighbours.
!
         IF (igrape.EQ.1) THEN
            hmax(ipart) = h(ipart)
         ENDIF

         temprho(ipart) = 0.
         divv(ipart) = 0.
         curlvx(ipart) = 0.
         curlvy(ipart) = 0.
         curlvz(ipart) = 0.
         numneigh1(ipart) = 0
         numneigh2(ipart) = 0
         rho1(ipart) = 0.
         rho2(ipart) = 0.
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(runtime)
      DO 50 n = nlst_in, nlst_end
         ipart = list(n)

         IF (iphase(ipart).GE.1) GOTO 50

         xi = x(ipart)
         yi = y(ipart)
         zi = z(ipart)
         vxi = vx(ipart)
         vyi = vy(ipart)
         vzi = vz(ipart)
         pmassi = pmass(ipart)
         hi = h(ipart)

         rhoi = 0.

!----    Computing rhoip2 is needed only if iexf=8, in order to
!        calculate mean opacity in externf.f (A. Okazaki, 03/11/2010)
!        rhoip2: density of secondary wind particles. The following
!        line initializes rhoip2.
         if (iexf.EQ.8) rhoip2 = 0.

         numneigh1i = 0
         numneigh2i = 0
         rho1i = 0.
         rho2i = 0.
         divvi = 0.
         curlvxi = 0.
         curlvyi = 0.
         curlvzi = 0.

         iptn = isnearpt(ipart)
!
!--No corrections if a neighbour of 2 point masses or more
!
         IF (iptn.GT.0) THEN
            iptcur = listpm(iptn)
            IF (iphase(iptcur).EQ.4) THEN
               xipt = x(iptcur)
               yipt = y(iptcur)
               zipt = z(iptcur)
               vxipt = vx(iptcur)
               vyipt = vy(iptcur)
               vzipt = vz(iptcur)
               xip = xi - xipt
               yip = yi - yipt
               zip = zi - zipt
               nlocalneigh = 0
            ENDIF
         ENDIF

         DO 40 k = 1, ilen(ipart)
            j = neighb(k,ipart)
            hj = h(j)

            IF (iptn.GT.0 .AND. iphase(iptcur).EQ.4) THEN
!
!--Find distance and list of all neighbours on the same side
!     of the point mass (not other point masses) and outside hacc
!
               dx = x(j) - xipt
               dy = y(j) - yipt
               dz = z(j) - zipt
               dot = xip*dx + yip*dy + zip*dz
               rr2 = dx*dx + dy*dy + dz*dz
               IF (iphase(j).EQ.0 .AND. dot.GT.0.0 .AND. &
                    rr2.GT.h(j)*h(j)) THEN
                  nlocalneigh = nlocalneigh + 1
                  neighloc(nlocalneigh) = j
                  rr1(nlocalneigh) = SQRT(rr2)
               ENDIF
            ENDIF

!!!CCCCCCC            IF (j.EQ.ipart .OR.
!!!CCCCCCC     &           ((hi.LT.hj .OR. (hi.EQ.hj .AND. j.LT.ipart)) .AND.
!!!CCCCCCC     &           iscurrent(j).EQ.1)) GOTO 40

            IF (iphase(ipart).EQ.-1 .OR. iphase(j).EQ.-1) GOTO 40
            IF (iphase(j).GE.1) THEN
               WRITE(iprint,*)'ERROR - ptmass in neighbour list ', j
               CALL quit
            ENDIF
!
!--No density contribution between particles across a point mass of type 4
!
            DO ii = 1, nptmass
               iptden = listpm(ii)
               IF (iphase(iptden).EQ.4) THEN
                  xii = x(iptden)
                  yii = y(iptden)
                  zii = z(iptden)
                  dxi = xii - xi
                  dyi = yii - yi
                  dzi = zii - zi
                  dxj = xii - x(j)
                  dyj = yii - y(j)
                  dzj = zii - z(j)
                  vpos = dxi*dxj + dyi*dyj + dzi*dzj
                  IF (vpos.LT.0.0) GOTO 40
               ENDIF
            END DO
!
!--Use mean h
!
            pmassj = pmass(j)
            hmean = 0.5*(hi + hj)
            hmean21 = 1./(hmean*hmean)
            hmean31 = hmean21/hmean
            hmean41 = hmean21*hmean21

            dx = xi - x(j)
            dy = yi - y(j)
            dz = zi - z(j)

            dvx = vxi - vx(j)
            dvy = vyi - vy(j)
            dvz = vzi - vz(j)

            rij2 = dx*dx + dy*dy + dz*dz + tiny
            v2 = rij2*hmean21

            IF (v2.LT.radkernel**2) THEN
               IF (igrape.EQ.1) THEN
                  internum = internum + 1
                  hmax(ipart) = MAX(hmax(ipart), hj)
               ENDIF

               rij1 = SQRT(rij2)
!
!--Get kernel quantities from interpolation in table
!
               index = v2/dvtable
               dxx = v2 - index*dvtable
               index1 = index + 1
               IF (index1.GT.itable) index1 = itable
               dwdx = (wij(index1) - wij(index))/dvtable
               wtij = (wij(index) + dwdx*dxx)*hmean31
               dgrwdx = (grwij(index1) - grwij(index))/dvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hmean41/rij1
!
!--Compute density
!
               rhoi = rhoi + pmassj*wtij
               IF (iexf.EQ.8 .AND. ibelong(j).EQ.2) THEN
                  rhoip2 = rhoip2 + pmassj*wtij
               ENDIF
               IF (iptn.GT.0 .AND. iphase(iptcur).EQ.4) THEN
                  coordi = xip*dx + yip*dy + zip*dz
                  IF (coordi.LT.0) THEN
                     numneigh2i = numneigh2i + 1
                     rho2i = rho2i + pmassj*wtij
                  ELSE
                     numneigh1i = numneigh1i + 1
                     rho1i = rho1i + pmassj*wtij
                  ENDIF
               ENDIF
!
!--Velocity divergence times density
!
               projv = grwtij*(dvx*dx + dvy*dy + dvz*dz)
               divvi = divvi - pmassj*projv
!
!--Velocity curl in 3D times density
!
               procurlvz = grwtij*(dvy*dx - dvx*dy)
               procurlvy = grwtij*(dvx*dz - dvz*dx)
               procurlvx = grwtij*(dvz*dy - dvy*dz)

               curlvxi = curlvxi - pmassj*procurlvx
               curlvyi = curlvyi - pmassj*procurlvy
               curlvzi = curlvzi - pmassj*procurlvz

!!!CCCCCCC               rhoij = (rho(ipart) + rho(j))/2.0
!!!CCCCCCC               IF (rhoij.NE.0.) THEN
!!!CCCCCCC                  wri = pmassj*wtij/rhoij
!!!CCCCCCC               ELSE
!!!CCCCCCC                  wri = 0.0
!!!CCCCCCC               ENDIF
!!!CCCCCCC
!!!CCCCCCC               IF (iscurrent(j).GT.0) THEN
!!!CCCCCCC                  IF (igrape.EQ.1) THEN
!!!CCCCCCC                     internum = internum + 1
!!!CCCCCCC                     hmax(j) = MAX(hmax(j), hi)
!!!CCCCCCC                  ENDIF
!!!CCCCCCC
!!!CCCCCCC                  temprho(j) = temprho(j) + pmassi*wtij
!!!CCCCCCC                  IF (isnearpt(j).GT.0) THEN
!!!CCCCCCC                     iptj = listpm(isnearpt(j))
!!!CCCCCCC                     xjp = x(j) - x(iptj)
!!!CCCCCCC                     yjp = y(j) - y(iptj)
!!!CCCCCCC                     zjp = z(j) - z(iptj)
!!!CCCCCCC                     coordj = xjp*dx + yjp*dy + zjp*dz
!!!CCCCCCC                     IF (coordj.GT.0) THEN
!!!CCCCCCC                        numneigh2(j) = numneigh2(j) + 1
!!!CCCCCCC                        rho2(j) = rho2(j) + pmassi*wtij
!!!CCCCCCC                     ELSE
!!!CCCCCCC                        numneigh1(j) = numneigh1(j) + 1
!!!CCCCCCC                        rho1(j) = rho1(j) + pmassi*wtij
!!!CCCCCCC                     ENDIF
!!!CCCCCCC                  ENDIF
!!!CCCCCCC                  divv(j) = divv(j) - pmassi*projv
!!!CCCCCCC                  curlvx(j) = curlvx(j) - pmassi*procurlvx
!!!CCCCCCC                  curlvy(j) = curlvy(j) - pmassi*procurlvy
!!!CCCCCCC                  curlvz(j) = curlvz(j) - pmassi*procurlvz
!!!CCCCCCC                  wrj = wri*pmassi/pmassj
!!!CCCCCCC               ENDIF
            ENDIF
 40      CONTINUE

!----    Calculate the mean opacity if iext=8.
!        Use akappa(1) or akappa(2) if iexf=7.
         IF (iexf.EQ.7 .OR. iexf.EQ.8) THEN
            IF (iphase(ipart).EQ.0) THEN
               DO iptm=1,nptmass
                  IF (akappar(iptm).NE.0.0) THEN
                     ddinv = rptmas(iptm) &
                             /SQRT((xi-x(iptm))*(xi-x(iptm)) &
                             +(yi-y(iptm))*(yi-y(iptm)) &
                             +(zi-z(iptm))*(zi-z(iptm)))
                     akapp(iptm) = akappac(iptm) &
                          +akappar(iptm) &
                          *(1.0d0-ddinv) &
                          **(2.0d0*vbeta(iptm)-1.0d0)
                  ELSE
                     akapp(iptm) = akappac(iptm)
                  ENDIF
               ENDDO
               IF (iexf.EQ.7) THEN
                  akappa(ipart) = akapp(ibelong(ipart))
               ELSE
                  IF (rhoi.NE.0.0) THEN
                     akappa(ipart) = akapp(1) &
                            +(akapp(2)-akapp(1))*rhoip2/rhoi
                  ELSE
                     akappa(ipart) = 0.0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF

         IF (iptn.GT.0 .AND. iphase(iptcur).EQ.4) THEN
!
!--Compute densities, and pressure gradients for the neighbours of a
!     point mass using a linear fit from each particle's neighbours
!     which are on the same side of the point mass and which are outside
!     hacc.
!
!--Sort in increasing radius
!
            IF (nlocalneigh.GT.1) THEN
               CALL indexx2(nlocalneigh, rr1, indexs)
!!!CCCCCCC            ELSE
!!!CCCCCCC               WRITE(iprint,*)'nlocal=',nlocalneigh,nneigh(ipart),ipart
            ENDIF
!
!--Unit vector in direction of mean tangential velocity of neighbours of ipart
!    and ipart itself.
!
            ript = SQRT(xip*xip + yip*yip + zip*zip)
            uxip = xip/ript
            uyip = yip/ript
            uzip = zip/ript
            vr = vxi*uxip + vyi*uyip + vzi*uzip
            tanvxi = vxi - vr*uxip
            tanvyi = vyi - vr*uyip
            tanvzi = vzi - vr*uzip
            tanv = SQRT(tanvxi*tanvxi + tanvyi*tanvyi + tanvzi*tanvzi)
            IF (tanv.NE.0.0) THEN
               tanvxi = tanvxi/tanv
               tanvyi = tanvyi/tanv
               tanvzi = tanvzi/tanv
            ENDIF
            tanvx1(ipart) = tanvxi
            tanvy1(ipart) = tanvyi
            tanvz1(ipart) = tanvzi

            tanviscfor1(ipart) = 0.
            ntanvisc = 0

            nmiddle = (1+nlocalneigh)/2
            rmean1 = 0.
            rhomean1 = 0.
            prmean1 = 0.
            DO ll = 1, nmiddle
               nl = neighloc(indexs(ll))
               rmean1 = rmean1 + rr1(indexs(ll))
               rhomean1 = rhomean1 + rho(nl)
               prmean1 = prmean1 + pr(nl)
               rjpt = rr1(indexs(ll))
               IF (rjpt.GT.ript) THEN
                  uxjp = (x(nl)-xipt)/rjpt
                  uyjp = (y(nl)-yipt)/rjpt
                  uzjp = (z(nl)-zipt)/rjpt
                  vxnl = vx(nl)
                  vynl = vy(nl)
                  vznl = vz(nl)
                  vrj = vxnl*uxjp + vynl*uyjp + vznl*uzjp
                  tanvx1j = vxnl - vrj*uxjp
                  tanvy1j = vynl - vrj*uyjp
                  tanvz1j = vznl - vrj*uzjp
                  tanvj = SQRT(tanvx1j*tanvx1j + tanvy1j*tanvy1j + &
                       tanvz1j*tanvz1j)
                  IF (tanvj.NE.0.0) THEN
                     tanvx1j = tanvx1j/tanvj
                     tanvy1j = tanvy1j/tanvj
                     tanvz1j = tanvz1j/tanvj
                  ENDIF

                  ntanvisc = ntanvisc + 1
!                  IF (ript.EQ.0.0) THEN
                     factor = 1.0
!                  ELSE
!                     factor = (rjpt/ript)**1.5
!                     factor = (rjpt/ript)**3.0
!                  ENDIF
                  tanviscfor1(ipart) = tanviscfor1(ipart) + &
                       (artvix(nl)*tanvx1j + artviy(nl)*tanvy1j + &
                       artviz(nl)*tanvz1j) * factor
               ENDIF
            END DO
            IF (nmiddle.NE.0) THEN
               rmean1 = rmean1/nmiddle
               rhomean1 = rhomean1/nmiddle
               prmean1 = prmean1/nmiddle
            ENDIF

            rmean2 = 0.
            rhomean2 = 0.
            prmean2 = 0.
            DO ll = nmiddle + 1, nlocalneigh
               nl = neighloc(indexs(ll))
               rmean2 = rmean2 + rr1(indexs(ll))
               rhomean2=rhomean2 + rho(nl)
               prmean2 = prmean2 + pr(nl)
               rjpt = rr1(indexs(ll))
               IF (rjpt.GT.ript) THEN
                  uxjp = (x(nl)-xipt)/rjpt
                  uyjp = (y(nl)-yipt)/rjpt
                  uzjp = (z(nl)-zipt)/rjpt
                  vxnl = vx(nl)
                  vynl = vy(nl)
                  vznl = vz(nl)
                  vrj = vxnl*uxjp + vynl*uyjp + vznl*uzjp
                  tanvx1j = vxnl - vrj*uxjp
                  tanvy1j = vynl - vrj*uyjp
                  tanvz1j = vznl - vrj*uzjp
                  tanvj = SQRT(tanvx1j*tanvx1j + tanvy1j*tanvy1j + &
                       tanvz1j*tanvz1j)
                  IF (tanvj.NE.0.0) THEN
                     tanvx1j = tanvx1j/tanvj
                     tanvy1j = tanvy1j/tanvj
                     tanvz1j = tanvz1j/tanvj
                  ENDIF

                  ntanvisc = ntanvisc + 1
!                  IF (ript.EQ.0.0) THEN
                     factor = 1.0
!                  ELSE
!                     factor = (rjpt/ript)**1.5
!                     factor = (rjpt/ript)**3.0
!                  ENDIF
                  tanviscfor1(ipart) = tanviscfor1(ipart) + &
                       (artvix(nl)*tanvx1j + artviy(nl)*tanvy1j + &
                       artviz(nl)*tanvz1j) * factor
               ENDIF
            END DO
            nsecond = nlocalneigh-nmiddle
            IF (nsecond.NE.0) THEN
               rmean2 = rmean2/nsecond
               rhomean2 = rhomean2/nsecond
               prmean2 = prmean2/nsecond
            ENDIF
            IF (ntanvisc.NE.0) THEN
               tanviscfor1(ipart) = tanviscfor1(ipart)/ntanvisc
!               IF (tanviscfor1(ipart).LT.0.0) THEN
!                  tanviscfor1(ipart)=0.0
!               ENDIF
            ENDIF

            diffr = rmean1 - rmean2
            IF(diffr.NE.0. .AND. nlocalneigh.GT.1) THEN
               fslope = (rhomean1 - rhomean2)/diffr
               prslope = (prmean1 - prmean2)/diffr
            ELSE
               fslope = 0.
               prslope = 0.
            ENDIF

            IF (rho(ipart).NE.0.) THEN
               prgradrho(ipart) = prslope/rho(ipart)
               densgrad(ipart) = fslope/rho(ipart)
            ELSE
               prgradrho(ipart) = 0.
               densgrad(ipart) = 0.
            ENDIF
!            IF (densgrad(ipart).GT.0.) densgrad(ipart) = 0.

         ENDIF

         hi3 = hi*hi*hi
         temprho(ipart) = temprho(ipart) + rhoi + &
              selfnormkernel*pmassi/hi3
         numneigh1(ipart) = numneigh1(ipart) + numneigh1i
         numneigh2(ipart) = numneigh2(ipart) + numneigh2i
         rho1(ipart) = rho1(ipart) + rho1i + &
              selfnormkernel*pmassi/hi3/2.0
         rho2(ipart) = rho2(ipart) + rho2i + &
              selfnormkernel*pmassi/hi3/2.0
         divv(ipart) = divv(ipart) + divvi
         curlvx(ipart) = curlvx(ipart) + curlvxi
         curlvy(ipart) = curlvy(ipart) + curlvyi
         curlvz(ipart) = curlvz(ipart) + curlvzi

 50   CONTINUE
!$OMP END DO

!$OMP DO SCHEDULE(runtime)
      DO 100 n = nlst_in, nlst_end
         ipart = list(n)
         rho(ipart) = cnormk*temprho(ipart)

!         IF (myrank == 0) THEN
!            IF (MOD(ipart,1000) == 1) THEN
!               WRITE (iprint,*) 'rho(',ipart,') =',rho(ipart)
!            END IF
!         END IF

         rhominimum = cnormk*selfnormkernel*pmass(ipart)/(h(ipart))**3
         IF (iphase(ipart).EQ.0 .AND. rho(ipart).LT.rhominimum) THEN
!!!CCCCCCC            IF (rho(ipart).LT.0.) THEN
!!!CCCCCCC               WRITE (iprint,*) 'Density < 0'
!!!CCCCCCC            ELSE
!!!CCCCCCC               WRITE (iprint,*) 'Density < rhomin'
!!!CCCCCCC            ENDIF
            xi = x(ipart)
            yi = y(ipart)
            zi = z(ipart)
            radius = SQRT(xi*xi + yi*yi + zi*zi)
!!!CCCCCCC            WRITE (iprint,*) '   ',radius, rho(ipart), xi, yi, zi
            rho(ipart) = rhominimum
         ENDIF

!         WRITE (iprint,88888) x(ipart),y(ipart),z(ipart),h(ipart),
!     &        pmass(ipart), rho(ipart)
!         WRITE (iprint,*) ilen(ipart)
!88888    FORMAT(6(1PE12.5,1X))
         rho1(ipart) = cnormk*rho1(ipart)
         rho2(ipart) = cnormk*rho2(ipart)
         divv(ipart) = cnormk*divv(ipart)
         curlv(ipart) = cnormk* &
              SQRT(curlvx(ipart)**2+curlvy(ipart)**2+curlvz(ipart)**2)

         IF (iphase(ipart).NE.0) GOTO 100
!
!--Pressure and sound velocity from ideal gas law...
!
         CALL eospg(ipart, u, rho, pr, vsound)
!
!--Find particle with highest density outside radcrit of point mass
!
         IF (iptmass.EQ.0) GOTO 100
         DO l = 1, nptmass
            iptcur = listpm(l)
            radi2 = (x(ipart) - x(iptcur))**2 + &
                 (y(ipart) - y(iptcur))**2 + (z(ipart) - z(iptcur))**2
            IF (radi2.LT.radcrit2) GOTO 100
         END DO
         rhonext = MAX(rhonext, rho(ipart))
 100  CONTINUE
!$OMP END DO
!$OMP END PARALLEL
!
!--Find particle with highest density outside radcrit of point mass
!
      IF (iptmass.NE.0) THEN
         DO n = nlst_in, nlst_end
            ipart = list(n)
            IF (rho(ipart).EQ.rhonext) THEN
               irhonex = ipart
            ENDIF
         END DO
      ENDIF
!
!--Possible to create a point mass
!
      IF (rhonext.GT.rhocrea .AND. nptmass.LT.iptdim .AND. &
                           iptmass.NE.0 .AND. icall.EQ.3) THEN
!
!--Make sure that all neighbours of point mass candidate are being
!     done on this time step. Otherwise, not possible to accrete
!     them to form a point mass and it may create a point mass without
!     accreting many particles!
!
         IF ((2.0*h(irhonex)).LT.hacc) THEN
            WRITE(iprint,*)'Ptmass creation passed h ', h(irhonex)

            CALL getneigh(irhonex,npart,h(irhonex),x,y,z,nlist,nearl)

            iokay = 1
            DO n = 1, nlist
               j = nearl(n)
               IF (it0(j).NE.itime) iokay = 0
            END DO
!
!--Set creation flag to true. Other tests done in accrete.f
!
            IF (iokay.EQ.1) THEN
               icreate = 1
               WRITE(iprint,*) ' and all particles on step'
            ELSE
               WRITE(iprint,*) ' but not all particles on step'
            ENDIF
         ELSE
            WRITE(iprint,*) 'Ptmass creation failed on h ',h(irhonex)
         ENDIF
      ENDIF
!
!--Point Mass Boundaries
!
      IF ((initialptm.EQ.4 .OR. iptmass.EQ.4) .AND. nptmass.GE.1) THEN
!$OMP PARALLEL DO SCHEDULE(runtime) default(none) &
!$OMP shared(nlst_in,nlst_end,list,numneighadd,isnearpt,iphase) &
!$OMP shared(listpm,x,y,z,h,wij,dvtable,pi,rho2,rho1) &
!$OMP shared(numneigh1,numneigh2,rho,u,pr,vsound,densgrad) &
!$OMP private(n,ipart,numneighold,iptnr,ipt,dx,dy,dz,r2,hmean2,rr) &
!$OMP private(vsmvfor,weight,vsmvfor2,indw,dxxw,indw1,dwdxw) &
!$OMP private(densgradrho,ck,valuen,rho1dif,numneighnew)
         DO n = nlst_in, nlst_end
            ipart = list(n)
            IF (iphase(ipart).EQ.0) THEN
               numneighold = numneighadd(ipart)
               numneighadd(ipart) = 0
!               rhoold(ipart) = rho(ipart)

               iptnr = isnearpt(ipart)
!--No corrections if a neighbour of 2 point masses or more
               IF (iptnr.GT.0) THEN
                  IF (iphase(listpm(iptnr)).EQ.4) THEN
                     ipt = listpm(iptnr)
!
!--Find what the contribution to a particle's density should be assuming
!     globally constant density gradient divided by density of particle
!     (which is making the assumption the radial pressure acceleration on
!     all particles is the same if isothermal).
!
                     dx = x(ipart)-x(ipt)
                     dy = y(ipart)-y(ipt)
                     dz = z(ipart)-z(ipt)
                     r2 = dx*dx + dy*dy + dz*dz
                     hmean2 = ((h(ipt)+2.0*h(ipart))/2.)**2
!
!--Only correct density if particle overlaps accretion radius by more
!     than h(ipart).  Because of kernel, density contribution for particles
!     that overlap less than this not really needed.  However, they DO
!     still need a pressure correction (because of shape of grad kernel).
!
                     IF (r2.LT.(4.*hmean2)) THEN
                        rr = SQRT(r2)
!!!                  htemp = MIN(h(ipart), h(ipt))
!!!                  vsmvfor = 2.0*(rr - h(ipt) - htemp)/htemp
                        vsmvfor = 2.0*(rr - h(ipt) - h(ipart))/h(ipart)
                        IF (vsmvfor.LE.0.0) THEN
                           weight = 1.0
                        ELSEIF (vsmvfor.GE.2.0) THEN
                           weight = 0.0
                        ELSE
                           vsmvfor2 = vsmvfor*vsmvfor
                           indw = vsmvfor2/dvtable
                           dxxw = vsmvfor2 - indw*dvtable
                           indw1 = indw + 1
                           IF (indw1.GT.itable) indw1 = itable
                           dwdxw = (wij(indw1) - wij(indw))/dvtable
                           weight = (wij(indw) + dwdxw*dxxw)
                        ENDIF

                        densgradrho = densgrad(ipart)

!                       IF (densgradrho .LT. -0.5/h(ipart))
!     &                 densgradrho = -0.5/h(ipart)
                        IF (densgradrho .LT. -1.0/h(ipart)) THEN
                           densgradrho = -1.0/h(ipart)
                        ELSEIF (densgradrho .GT. 1.0/h(ipart)) THEN
                           densgradrho = 1.0/h(ipart)
                        ENDIF
                        ck = 8.0/pi*31.0/140.0*h(ipart)*densgradrho
                        valuen = (1.0 - ck)/(1.0 + ck)
                        rho1dif = (rho2(ipart)*valuen - rho1(ipart))* &
                             weight
                        IF (rho1dif.LT.0.) rho1dif = 0.
!
!--Correct number of neighbours for boundary of point mass
!
!                       IF (valuen.LT.1.0) valuen = 1.0
                        numneighnew = INT((valuen* &
                        numneigh2(ipart) - numneigh1(ipart))*weight)
                        IF (numneighnew.LT.0) numneighnew=0

                        IF (numneighnew .LT. numneighold - 1) THEN
                           numneighadd(ipart) = numneighold - 1
                        ELSEIF (numneighnew .GT. numneighold + 1) THEN
                           numneighadd(ipart) = numneighold + 1
                        ELSE
                           numneighadd(ipart) = numneighnew
                        ENDIF

!                  numneighadd(ipart)=0

!
!--Correct density of particle for boundary of point mass
!
                        rho(ipart) = rho(ipart) + rho1dif
                        CALL eospg(ipart,u,rho,pr,vsound)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         END DO
!$OMP END PARALLEL DO
      ENDIF

      IF (itrace.EQ.'all') WRITE (iprint,300)
  300 FORMAT ('exit subroutine densityi')

      END SUBROUTINE densityi
