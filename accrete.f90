      SUBROUTINE accrete(dt, realtime, isave)
!************************************************************
!                                                           *
!  This routine does accretion for pt masses                *
!                                                           *
!************************************************************


      use mpi_mod
      use idims

      use constants
      use units
      use logun
      use part
      use densi
      use typef
      use carac
      use cgas
      use kerne
      use ener
      use divve
      use eosq
      use bodys
      use force, only: f1vx, f1vy, f1vz, f1u, f1h, &
                 f2vx, f2vy, f2vz, f2u, f2h
      use btree
      use fracg
      use polyk2
      use phase
      use ptmass
      use nextmpt
      use neighbor
      use tlist
      use ptdump
      use curlist
      use timei
      use init
      use ghost
      use table
      use active
      use accnum
      use varet
      use binary
      use accurpt
      use delay
      use new
      use sphcom
      use crpart
      use misali
      use capt

      implicit none

      INTEGER(I4B) :: isave
      INTEGER(I4B) :: iremove(idim), numberacc(iptdim)
      INTEGER(I4B) :: i, j, iii, jjj, i1list, i2list, icount, ikk, &
                imerge, index, index1, ipt, iptm, iptm1, iptm2, &
                itemp, j1, j2, jj, k, k1, k2, kk, nboundto, &
                nvalid
      REAL(DP) :: dt, realtime
      REAL(DP) :: ptminner(iptdim)
      REAL(DP) :: xi, yi, zi, hi, pmassi, vxi, vyi, vzi, &
                f1vxi, f1vyi, f1vzi, f1hi, f1ui, f2vxi, f2vyi, &
                f2vzi, f2hi, f2ui, dum2vxi, dum2vyi, dum2vzi, &
                ui, rhoi, dphiti, dgravi, poteni, pri, divvi, &
                spinxi, spinyi, spinzi, divai, gama1
      INTEGER(I4B) :: it0i, it1i, istepsi
      REAL(DP) :: alphapt, alphatot,  betatotal, cnorm2, &
                hacc2, hacc11, hacc21, hacc3
      REAL(DP) :: d12, dax, day, daz, dfptdx, dgrwdx, &
                divvr2, dtaccj, dx, dy, dz, dxx, &
                dvx, dvy, dvz, dvxtemp, dvytemp, &
                dvztemp, d2vxtemp, d2vytemp, d2vztemp, &
                dwdx, f1vxtemp, f1vytemp, f1vztemp, &
                f2vxtemp, f2vytemp, f2vztemp, grwtij
      REAL(DP) :: haccall2, hacccur, hacccur2, haccmin2, &
                hj3, hmean, hmean21, hmean31, hmean41, hratio, &
                haccmin
      REAL(DP) :: phi, pmass1, pmass2, pmassj, pmassj1, &
                pmassj2, pmassjn, pmassnew, pmiold, &
                proja, q13, q23, r1, r1temp, r2, r2xy, r2xz, &
                r2yz, rij, rij1, rij2, roche, roche2, &
                rvx, rvy, rvz, rx, ry, rz, rxtemp, rytemp, rztemp
      REAL(DP) :: specangmom2, specangmomhacc2, spinm, &
                spinxtemp, spinytemp, spinztemp, tcomp, total, &
                totalmass, trot, troty, utemp, v, v2, v2temp, &
                vkep2, vpotdif, vpottemp, vrad, vrad2, vtan2, &
                vtot2, x1, y1, z1, xiold, xj1, xmjeans, xtemp, &
                yiold, yj1, ytemp, ziold, zj1, ztemp, &
                vx1, vy1, vz1, wptj
      CHARACTER(len=7) :: where='accrete'

!!    LOGICAL ifirst
!!    DATA ifirst/.true./
!!    SAVE rscalebnd,rscalebnd2,vk0,hr0,height0,height5,gg2

      cnorm2 = 0.75 - haccall + 0.5*haccall**3 - 0.1875*haccall**4
      cnorm2 = 2. / cnorm2
      hacc2  = hacc*hacc
      hacc11 = 2./hacc
      hacc21 = 4./hacc2
      hacc3  = hacc*hacc2
!
!--Only allow accretion of GAS particles evaluated at the CURRENT timestep
!
      DO i = 1, npart + nghost
         iremove(i) = -1
      END DO
      DO i = 1, nlst0
         j = llist(i)
         IF (iphase(j).EQ.0) iremove(j) = 0
      END DO

      DO iii = 1, nptmass
         numberacc(iii) = 0
         ptminner(iii) = 0.
      END DO
!
!--CREATION OF NEW POINT MASS, by accretion of particles to form it
!
      IF (icreate.EQ.1 .AND. iremove(irhonex).EQ.0) THEN
!
!--Test before creating a new point mass:
!
!     (a) whether the particle(irhonex) contains a jeans mass within 2h
!           =>calculate alpha
!
!     (b) calculate beta ratio of rotational to gravitational pot. energy
!
!     (c) calculate the divergence of the acceleration
!           -ve => self gravitating/collapsing
!           +ve => in process of tidal/disruption or core bounce
!
!
!--Fraction of hacc for accreting a particle regardless of tests
!
         IF(iptmass.EQ.1) THEN
            haccmin2 = hacc2
         ELSE
            haccmin2 = haccall*haccall
         ENDIF

         xi = x(irhonex)
         yi = y(irhonex)
         zi = z(irhonex)
         hi = h(irhonex)
         pmassi = pmass(irhonex)
         vxi = vx(irhonex)
         vyi = vy(irhonex)
         vzi = vz(irhonex)
         f1vxi = f1vx(irhonex)
         f1vyi = f1vy(irhonex)
         f1vzi = f1vz(irhonex)
         f1hi = f1h(irhonex)
         f1ui = f1u(irhonex)
         f2vxi = f2vx(irhonex)
         f2vyi = f2vy(irhonex)
         f2vzi = f2vz(irhonex)
         f2hi = f2h(irhonex)
         f2ui = f2u(irhonex)
         dum2vxi = dum2vx(irhonex)
         dum2vyi = dum2vy(irhonex)
         dum2vzi = dum2vz(irhonex)
         ui  =  u(irhonex)
         rhoi = rho(irhonex)
         dphiti = dphit(irhonex)
         dgravi = dgrav(irhonex)
         poteni = poten(irhonex)
         pri = pr(irhonex)
         divvi = divv(irhonex)
         it0i = it0(irhonex)
         it1i = it1(irhonex)
         istepsi = isteps(irhonex)
         spinxi = 0.
         spinyi = 0.
         spinzi = 0.

         tkin = 0.
         trotx = 0.
         troty = 0.
         trotz = 0.
         tgrav = 0.
         divai = 0.
         gama1 = gamma - 1.0
         IF ( varsta.NE.'entropy' ) THEN
            tterm = pmassi*ui
         ELSEIF (gama1.EQ.0.) THEN
            tterm = 1.5*ui
         ELSE
            tterm = pmassi*ui*rhoi**gama1/gama1
         ENDIF

         CALL getneigh(irhonex, npart, hi, x, y, z, nlist, nearl)
         if(myrank.eq.0) &
         write (iprint,*) 'getneigh ',irhonex,npart,hi,nlist

         IF (nlist.GT.iptneigh) CALL error(where,2)
         nptlist(nptmass+1) = nlist

         nvalid = 0
         DO k = 1, nlist
            j = nearl(k)
            nearpt(nptmass+1,k) = j
            dx = x(j) - xi
            dy = y(j) - yi
            dz = z(j) - zi
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            IF (iremove(j).EQ.0 .AND. rij2.LT.hacc2) THEN
               nvalid = nvalid + 1
               hmean = 0.5*(hi + h(j))
               hmean21 = 1./(hmean*hmean)
               hmean31 = hmean21/hmean
               hmean41 = hmean21*hmean21
               dvx = vx(j) - vxi
               dvy = vy(j) - vyi
               dvz = vz(j) - vzi
!
!--Relative kinetic energy, tkin
!
               vtot2 = dvx*dvx + dvy*dvy + dvz*dvz
               tkin = tkin + pmass(j)*vtot2
!
!--Relative rotational energy around x, trotx
!
               r2yz = dz*dz + dy*dy
               rvx = dy*dvz - dz*dvy
               IF(r2yz.NE.0.) trotx = trotx + pmass(j)*rvx*rvx/r2yz
!
!--Relative rotational energy around y, troty
!
               r2xz = dx*dx + dz*dz
               rvy = dz*dvx - dx*dvz
               IF(r2xz.NE.0.) troty = troty + pmass(j)*rvy*rvy/r2xz
!
!--Relative rotational energy around z, trotz
!
               r2xy = dx*dx + dy*dy
               rvz = dx*dvy - dy*dvx
               IF(r2xy.NE.0.) trotz = trotz + pmass(j)*rvz*rvz/r2xy

               v2 = rij2*hmean21
               rij = SQRT(rij2)
               rij1 = 1.0/rij
               v = rij/hmean
!
!--Get kernel quantities from interpolation in table
!
               IF (v.LT.radkernel) THEN
                  index = v2/dvtable
                  dxx = v2 - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dgrwdx = (grwij(index1) - grwij(index))/dvtable
                  grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phi = (fpoten(index) + dfptdx*dxx)/hmean
                  IF (v2.GT.1.) phi = phi + rij1/15.0
               ELSE
                  grwtij = 0.
                  phi = -rij1
               ENDIF
!
!--Acceleration divergence times density...
!
               dax = f1vx(j) - f1vxi
               day = f1vy(j) - f1vyi
               daz = f1vz(j) - f1vzi
               proja = grwtij*(dax*dx + day*dy + daz*dz)/rij
               divai = divai - pmass(j)*proja
!
!--Gravitational energy of particles
!
               tgrav = tgrav + phi*pmass(j)*pmassi
!
!--Thermal energy, tterm
!
               IF ( varsta.NE.'entropy' ) THEN
                  tterm = tterm + pmass(j)*u(j)
               ELSEIF (gama1.EQ.0.) THEN
                  tterm = tterm + 1.5*u(j)
               ELSE
                  tterm = tterm + pmass(j)*u(j)*rho(j)**gama1/gama1
               ENDIF
            ENDIF
         END DO
!
!--Normalise acceleration divergence, divai
!
         if(myrank.eq.0) &
         WRITE(iprint,*) 'Ptmass nvalid, nlist ',nvalid, nlist, &
              irhonex,h(irhonex),rho(irhonex),iphase(irhonex), &
              nneigh(irhonex)
         divai = cnormk*divai
         IF (divai.GE.0) THEN
            if(myrank.eq.0) &
            WRITE(iprint,*)'Divai +ve => no pt mass creation yet ', &
                 divai
            icreate = 0
            GOTO 100
         ENDIF
!
!--Other normalisations
!
         tkin = 0.5*tkin
         trotx = 0.5*trotx
         troty = 0.5*troty
         trotz = 0.5*trotz
         trot = SQRT(trotx*trotx + troty*troty + trotz*trotz)
!
!--Now calculate tgrav, the potential energy => alpha
!
         DO k1 = 1, nlist
            j1 = nearl(k1)
            xj1 = x(j1)
            yj1 = y(j1)
            zj1 = z(j1)
            pmassj1 = pmass(j1)

            rx = xj1 - xi
            ry = yj1 - yi
            rz = zj1 - zi
            r2 = rx*rx + ry*ry + rz*rz
            IF (iremove(j1).NE.0 .OR. r2.GE.hacc2) GOTO 20

            DO k2 = k1+1, nlist
               j2 = nearl(k2)

               rx = x(j2) - xi
               ry = y(j2) - yi
               rz = z(j2) - zi
               r2 = rx*rx + ry*ry + rz*rz
               IF (iremove(j2).NE.0 .OR. r2.GE.hacc2) GOTO 10

               dx = x(j2) - xj1
               dy = y(j2) - yj1
               dz = z(j2) - zj1
               rij2 = dx*dx + dy*dy + dz*dz + tiny
               rij = SQRT(rij2)
               rij1 = 1./rij
               pmassj2 = pmass(j2)
!
!--Define mean h ...
!
               hmean = 0.5*(h(j1) + h(j2))
               hmean21 = 1./(hmean*hmean)
               v2 = rij2*hmean21
               IF (v2.LT.radkernel*radkernel) THEN
                  index = v2/dvtable
                  dxx = v2 - index*dvtable
                  index1 = index + 1
                  IF (index1.GT.itable) index1 = itable
                  dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                  phi = (fpoten(index) + dfptdx*dxx)/hmean
                  IF (v2.GT.1.) phi = phi + rij1/15.0
               ELSE
                  phi = -rij1
               ENDIF
               tgrav = tgrav + phi*pmassj1*pmassj2
 10         END DO
 20      END DO
!
!--Now test to see if the particles to be turned into a point mass
!     satisfy all of the criteria
!
         xmjeans = ABS(tgrav/tterm)
         alphapt = 1./xmjeans
         betatotal = ABS(trot/tgrav)
         alphatot = alphapt + betatotal
         total = tterm + tgrav + tkin
         IF(alphapt.GT.0.5) THEN
          if(myrank.eq.0)then
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on alpha = ',alphapt
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans, &
                 alphapt,betatotal,rho(irhonex),h(irhonex),irhonex
            FLUSH(iprint)
          endif
            icreate = 0
            GOTO 100
         ELSE IF (alphatot.GT.1.0) THEN
           if(myrank.eq.0)then
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on alpha + beta = ', &
                           alphapt,betatotal
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans, &
                 alphapt,betatotal,rho(irhonex),h(irhonex),irhonex
           endif
            icreate = 0
            if(myrank.eq.0)FLUSH(iprint)
            GOTO 100
         ELSE IF(total.GE.0) THEN
          if(myrank.eq.0)then
            WRITE(iprint,99001)
            WRITE(iprint,*)'Ptmass failed on total energy (pos) = ', &
                 total
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans, &
                 alphapt,betatotal,rho(irhonex),h(irhonex),irhonex
           endif
            icreate = 0
            if(myrank.eq.0)FLUSH(iprint)
            GOTO 100
         ELSE
          if(myrank.eq.0)then
            WRITE(iprint,99002)
            WRITE(iprint,99003)total,tterm,tgrav,tkin,xmjeans, &
                 alphapt,betatotal,rho(irhonex),h(irhonex),irhonex
            FLUSH(iprint)
          endif
99001       FORMAT(' PROTOSTAR FORMATION UNSUCCESSFUL !!!')
99002       FORMAT(' PROTOSTAR FORMATION SUCCESSFUL !!!')
99003       FORMAT(' total energy                   :',1PE14.5,/, &
           ' thermal energy                 :',1PE14.5,/, &
           ' gravitational potential energy :',1PE14.5,/, &
           ' kinetic energy                 :',1PE14.5,/, &
           ' Jeans no.                      :',1PE14.5,/, &
           ' alpha                          :',1PE14.5,/, &
           ' beta total                     :',1PE14.5,/, &
           ' rho(irhonex)                   :',1PE14.5,/, &
           ' h(irhonex)                     :',1PE14.5,/, &
           ' irhonex                        :',I6)
         ENDIF
!
!--Create point mass from particles
!
!--Make a dump next time save is called
!
         iptcreat = 1
         nptmass = nptmass + 1
         IF (nptmass.GE.iptdim) CALL error(where,3)

         numberacc(nptmass) = 0
         ptminner(nptmass) = 0.0
         DO jj = 1,nptlist(nptmass)
            j = nearpt(nptmass,jj)
            rx = x(j)-xi
            ry = y(j)-yi
            rz = z(j)-zi
            dvx = vx(j)-vxi
            dvy = vy(j)-vyi
            dvz = vz(j)-vzi

            r2 = rx*rx+ry*ry+rz*rz
            IF (r2.LT.hacc2.AND.iremove(j).EQ.0) THEN
               iremove(j) = 1
               iphase(j) = -1
               iaccr = 1
               numberacc(nptmass) = numberacc(nptmass) + 1
               ptminner(nptmass) = ptminner(nptmass) + pmass(j)
!
!--Accrete particle if it lies within hacc
!      accrete mass
!      angular momentum
!      linear momentum
!
               totalmass = pmassi + pmass(j)
               spinm = pmass(j)*pmassi/totalmass
               spinxi = spinxi + spinm*(ry*dvz - dvy*rz)
               spinyi = spinyi + spinm*(dvx*rz - rx*dvz)
               spinzi = spinzi + spinm*(rx*dvy - dvx*ry)
               vxi = (pmassi*vxi + pmass(j)*vx(j))/totalmass
               vyi = (pmassi*vyi + pmass(j)*vy(j))/totalmass
               vzi = (pmassi*vzi + pmass(j)*vz(j))/totalmass
               xi = (pmassi*xi + pmass(j)*x(j))/totalmass
               yi = (pmassi*yi + pmass(j)*y(j))/totalmass
               zi = (pmassi*zi + pmass(j)*z(j))/totalmass
               f1vxi = (pmassi*f1vxi + pmass(j)*f1vx(j))/totalmass
               f1vyi = (pmassi*f1vyi + pmass(j)*f1vy(j))/totalmass
               f1vzi = (pmassi*f1vzi + pmass(j)*f1vz(j))/totalmass
               f2vxi = (pmassi*f2vxi + pmass(j)*f2vx(j))/totalmass
               f2vyi = (pmassi*f2vyi + pmass(j)*f2vy(j))/totalmass
               f2vzi = (pmassi*f2vzi + pmass(j)*f2vz(j))/totalmass
               dum2vxi=(pmassi*dum2vxi + pmass(j)*dum2vx(j))/totalmass
               dum2vyi=(pmassi*dum2vyi + pmass(j)*dum2vy(j))/totalmass
               dum2vzi=(pmassi*dum2vzi + pmass(j)*dum2vz(j))/totalmass
               ui = (pmassi*ui + pmass(j)*u(j))/totalmass
               if(myrank.eq.0)WRITE(iprint,*)'add = ', jj, r2, pmassi
               pmassi = totalmass
            ENDIF
         END DO
        if(myrank.eq.0)then
         WRITE(iprint,*)'PROTOSTAR CREATION, mass = ', pmassi, &
                          ' time = ', realtime
         WRITE(iprint,*)'   accretion radius = ', hacc
        endif
!
!--Set spin arrays for ptmass and change the number of point masses
!
         spinx(nptmass) = spinxi
         spiny(nptmass) = spinyi
         spinz(nptmass) = spinzi
         spinadx(nptmass) = spinxi
         spinady(nptmass) = spinyi
         spinadz(nptmass) = spinzi
!
!--Finally set new point mass's properties
!
         x(irhonex) = xi
         y(irhonex) = yi
         z(irhonex) = zi
         vx(irhonex) = vxi
         vy(irhonex) = vyi
         vz(irhonex) = vzi
         xmomsyn(nptmass) = pmassi*vxi
         ymomsyn(nptmass) = pmassi*vyi
         zmomsyn(nptmass) = pmassi*vzi
         xmomadd(nptmass) = 0.0
         ymomadd(nptmass) = 0.0
         zmomadd(nptmass) = 0.0

         f1vx(irhonex) = f1vxi
         f1vy(irhonex) = f1vyi
         f1vz(irhonex) = f1vzi
         f1h(irhonex) = f1hi
         f1u(irhonex) = f1ui
         f2vx(irhonex) = f2vxi
         f2vy(irhonex) = f2vyi
         f2vz(irhonex) = f2vzi
         f2h(irhonex) = f2hi
         f2u(irhonex) = f2ui
         dum2vx(irhonex) = dum2vxi
         dum2vy(irhonex) = dum2vyi
         dum2vz(irhonex) = dum2vzi
         pmass(irhonex) = pmassi
         ptmsyn(nptmass) = pmassi
         ptmadd(nptmass) = 0.0
         angaddx(nptmass) = 0.0
         angaddy(nptmass) = 0.0
         angaddz(nptmass) = 0.0

         h(irhonex) = hacc

         u(irhonex) = ui
         rho(irhonex) = rhoi
         dphit(irhonex) = 0.
         dgrav(irhonex) = 0.
         poten(irhonex) = poteni
         pr(irhonex) = pri
         divv(irhonex) = divvi
         IF (istepmin.LT.istepsi) THEN
            it0(irhonex) = it0i
            it1(irhonex) = it0(irhonex) + istepmin/2
            isteps(irhonex) = istepmin
         ELSE
            it0(irhonex) = it0i
            it1(irhonex) = it1i
            isteps(irhonex) = istepsi
         ENDIF

         CALL getneigh(irhonex,npart,h(irhonex),x,y,z,nlist,nearl)

         IF (nlist.GT.iptneigh) CALL error(where,2)
         nptlist(nptmass) = nlist

         DO k = 1, nlist
            nearpt(nptmass,k) = nearl(k)
         END DO

         iphase(irhonex) = iptmass
         IF (initialptm.EQ.0) initialptm = iptmass
         listpm(nptmass) = irhonex
         hasghost(irhonex) = 0
      ENDIF
!
!--ACCRETION OF PARTICLES NEAR AN EXISTING POINT MASS
!
!--Method for accreting a particle
!     iphase = 1  Point mass that accretes everything regardless of tests.
!              2  Point mass that accretes particles if they pass tests.
!              3  Point mass that accretes part of a particle until all gone,
!                    they must also pass tests.
!              4  Point mass with accretion radius boundary corrections
!                   (a) smoothing length corrections
!                   (b) local density gradient correction of density
!                   (c) local pressure gradient correction to pressure force
!                   (d) local shear viscosity correction
!
!     hacc    = outer radius at which particle's accretion begins
!     haccall = radius at which all particles are accreted without test
!
!
 100  DO iii = 1, nptmass
         i = listpm(iii)

!--Don't allow accretion onto point mass 1 (=Be star) unless
!  a particle is not inside haccall.
         hacccur = h(i)
         IF (i.EQ.1) THEN
            IF (ibound.EQ.0 .OR. ibound.GE.90) THEN
               hacccur = rptmas(1)
            ELSE
               hacccur = hacc
            ENDIF
         ENDIF

         IF (iphase(i).EQ.3) THEN
            haccmin = haccall
         ELSE
            haccmin = hacccur
         ENDIF
         hacccur2 = hacccur*hacccur
         haccmin2 = haccmin*haccmin
         haccall2 = haccall*haccall

!----    For point mass 1 (=Be star), haccall2 should be
!        rptmas(i)**2, not haccall**2
!        (14 September 2003, A. Okazaki)
         IF (i.EQ.1) THEN
            IF (ibound.EQ.0 .OR. ibound.GE.90) THEN
               haccall2 = hacccur2
            ENDIF
         ENDIF

         xi = x(i)
         yi = y(i)
         zi = z(i)
         xiold = xi
         yiold = yi
         ziold = zi
         hi = h(i)
         pmassi = pmass(i)
         vxi = vx(i)
         vyi = vy(i)
         vzi = vz(i)
         f1vxi = f1vx(i)
         f1vyi = f1vy(i)
         f1vzi = f1vz(i)
         f2vxi = f2vx(i)
         f2vyi = f2vy(i)
         f2vzi = f2vz(i)
         dum2vxi = dum2vx(i)
         dum2vyi = dum2vy(i)
         dum2vzi = dum2vz(i)
         ui = u(i)

         xtemp = 0.0
         ytemp = 0.0
         ztemp = 0.0
         f1vxtemp = 0.0
         f1vytemp = 0.0
         f1vztemp = 0.0
         f2vxtemp = 0.0
         f2vytemp = 0.0
         f2vztemp = 0.0
         d2vxtemp = 0.0
         d2vytemp = 0.0
         d2vztemp = 0.0
         utemp = 0.0
         spinxtemp = 0.0
         spinytemp = 0.0
         spinztemp = 0.0

!----    Calculation of the Roche lobe radius to test whether
!        particles enter the Roche radii of the companion.
!        (A. Okazaki, 02/02/2007)
         IF (nptmass.EQ.2 .AND. i.NE.isphcom) THEN
            q13 = (pmass(listpm(2))/pmass(listpm(1)))**(1.0/3.0)
            q23 = q13*q13
            d12 = SQRT((x(listpm(2))-x(listpm(1))) &
                      *(x(listpm(2))-x(listpm(1))) &
                     + (y(listpm(2))-y(listpm(1))) &
                      *(y(listpm(2))-y(listpm(1))) &
                     + (z(listpm(2))-z(listpm(1))) &
                      *(z(listpm(2))-z(listpm(1))))
!----       Eggleton (1983)
            roche = d12*0.49*q23/(0.6*q23+LOG(1.0+q13))
            roche2 = roche*roche
         ENDIF

         DO jj = 1, nptlist(iii)
            j = nearpt(iii,jj)
            rx = x(j) - xi
            ry = y(j) - yi
            rz = z(j) - zi
            r2 = rx*rx + ry*ry + rz*rz

!--         Is point mass neighbour inside Roche radius?
!           (A. Okazaki, 02/02/2007)
            IF (nptmass.EQ.2 .AND. i.NE.isphcom) THEN
               iin(j) = 0
               IF (iphase(j).EQ.0) THEN
                  IF (r2.LT.roche2) THEN
                     IF (iinold(j).EQ.0 &
                       .AND. iremove(j).EQ.0) THEN
                        ncapt = ncapt + 1
                       if(myrank.eq.0)then
!!                        WRITE(icaptpr) j,realtime,x(j),y(j),z(j),
!!     &                      vx(j),vy(j),vz(j),h(j),u(j),
!!     &                      poten(j),dgrav(j),i,roche,
!!     &                      (x(ipt),y(ipt),z(ipt),
!!     &                      vx(ipt),vy(ipt),vz(ipt),ipt=1,nptmass)
!----                   Since 10/29/2011 (A. Okazaki)
!                        WRITE(icaptpr) j,realtime,x(j),y(j),z(j), &
!                            vx(j),vy(j),vz(j),h(j),u(j), &
!                            poten(j),dgrav(j),pmass(j),i,roche, &
!                            (x(ipt),y(ipt),z(ipt), &
!                            vx(ipt),vy(ipt),vz(ipt),ipt=1,nptmass)
!                        FLUSH (icaptpr)
                       endif
                     ENDIF
                     iin(j) = -1
                  ENDIF
               ELSE
                  iin(j) = -1
               ENDIF
               iinold(j) = iin(j)
            ENDIF
!
!--Is point mass neighbour inside accretion radius and accretable?
!
            IF (r2.LT.hacccur2 .AND. iremove(j).EQ.0) THEN
!
!--Check to ensure that the particle is actually bound to the
!     point mass and not just passing through its neighbourhood
!      (a) particle must be bound
!      (b) particle must be more bound to current point mass than any other
!      (c) specific angular momentum of particle must be less than that
!             required for it to form a circular orbit at hacc (8/9/94)
!
               dvx = vx(j) - vxi
               dvy = vy(j) - vyi
               dvz = vz(j) - vzi

               r1 = SQRT(r2)
               divvr2 = dvx*dvx + dvy*dvy + dvz*dvz
               vrad = (dvx*rx + dvy*ry + dvz*rz)/r1
               vrad2 = vrad*vrad
               vkep2 = pmassi/r1
               vpotdif = -vkep2 + divvr2/2.
               vtan2 = divvr2 - vrad2
               specangmom2 = vtan2*r2
               specangmomhacc2 = pmassi*hacccur
!
!--Check which point mass the particle is MOST bound to (6/8/94)
!
               nboundto = iii
               DO kk = 1, nptmass
                  IF (kk.NE.iii) THEN
                     ikk = listpm(kk)
                     rxtemp = x(j) - x(ikk)
                     rytemp = y(j) - y(ikk)
                     rztemp = z(j) - z(ikk)
                     r1temp = SQRT(rxtemp*rxtemp + rytemp*rytemp + &
                          rztemp*rztemp)
                     dvxtemp = vx(j) - vx(ikk)
                     dvytemp = vy(j) - vy(ikk)
                     dvztemp = vz(j) - vz(ikk)
                     v2temp = dvxtemp*dvxtemp + dvytemp*dvytemp + &
                          dvztemp*dvztemp
                     vpottemp = -pmass(ikk)/r1temp + v2temp/2.
                     IF (vpottemp.LT.vpotdif) nboundto = kk
                  ENDIF
               END DO
               IF (nboundto.NE.iii) vpotdif = 1.0E+10

               IF((vpotdif.LE.0. .AND. specangmom2.LT.specangmomhacc2) &
                    .OR. r2.LT.haccall2 .OR. iphase(i).EQ.1) THEN

                  IF (r2.LT.haccmin2 .OR. (iphase(i).EQ.3 .AND. &
                       pmassjn.LT.pmassleast)) THEN
!
!--Accrete all of particle's mass, and set particle accreted (iphase=-1).
!
                     iphase(j) = -1
                     iremove(j) = 1
                     iaccr = 1
                     numberacc(iii) = numberacc(iii) + 1
                     pmassj = pmass(j)
                     pmassjn = pmass(j)
                     ptminner(iii) = ptminner(iii) + pmassj
                     IF (iphase(i).EQ.3) THEN
                        if(myrank.eq.0)WRITE(iprint,*)'accreted ',j,r1, &
                          pmassj/pmassleast, pmassleast
                     ENDIF
                    if(myrank.eq.0)then
                     WRITE(iaccpr) j,realtime,x(j),y(j),z(j), &
                            vx(j),vy(j),vz(j),h(j),u(j), &
                            poten(j),dgrav(j),pmass(j),i,hacccur, &
                            (x(ipt),y(ipt),z(ipt), &
                            vx(ipt),vy(ipt),vz(ipt),ipt=1,nptmass)
                     FLUSH (iaccpr)
                    endif
                  ELSEIF (iphase(i).EQ.3) THEN
!
!--Point mass accretes part of a particle's mass at time
!
!--Use timestep for relevant accretion timestep
!
                     dtaccj = dt*isteps(j)/imaxstep

                     hj3 = h(j)*h(j)*h(j)
                     hratio = hacc3/hj3/8.
                     IF (hratio.GT.1.) hratio = 1.0
                     v2 = r2*hacc21
!
!--Use smoothing kernel of point mass renormalised as accretion kernel
!
                     index = v2/dvtable
                     dxx = v2 - index*dvtable
                     index1 = index + 1
                     IF (index1.GT.itable) index1 = itable
                     dwdx = (wij(index1) - wij(index))/dvtable
                     wptj = 0.333333*(wij(index) + dwdx*dxx)* &
                                                      hacc11*cnorm2
                     pmassj = -wptj*pmass(j)*vrad*hratio*dtaccj/hacc
                     pmassjn = pmass(j) - pmassj
                  ELSE
                     CALL error(where,1)
                  ENDIF
!
!--accrete particle if it lies within h and above condtions okay
!      accrete mass
!      angular momentum
!      linear momentum
!
                  IF (notacc(j)) GOTO 888

                  pmiold = ptmsyn(iii) + ptmadd(iii)
                  ptmadd(iii) = ptmadd(iii) + pmassj
                  totalmass = ptmsyn(iii) + ptmadd(iii)
                  spinm = pmassj*pmiold/totalmass
                  spinxtemp = spinxtemp + spinm*(ry*dvz - dvy*rz)
                  spinytemp = spinytemp + spinm*(dvx*rz - rx*dvz)
                  spinztemp = spinztemp + spinm*(rx*dvy - dvx*ry)

                  xmomadd(iii) = xmomadd(iii) + pmassj*vx(j)
                  ymomadd(iii) = ymomadd(iii) + pmassj*vy(j)
                  zmomadd(iii) = zmomadd(iii) + pmassj*vz(j)
                  angaddx(iii) = angaddx(iii) + pmassj* &
                       (yi*vz(j) - vy(j)*zi)
                  angaddy(iii) = angaddy(iii) + pmassj* &
                       (vx(j)*zi - xi*vz(j))
                  angaddz(iii) = angaddz(iii) + pmassj* &
                       (xi*vy(j) - vx(j)*yi)

                  xtemp = xtemp + pmassj*x(j)
                  ytemp = ytemp + pmassj*y(j)
                  ztemp = ztemp + pmassj*z(j)

                  xi = (xiold*pmassi + xtemp)/totalmass
                  yi = (yiold*pmassi + ytemp)/totalmass
                  zi = (ziold*pmassi + ztemp)/totalmass
                  vxi = (xmomsyn(iii) + xmomadd(iii))/totalmass
                  vyi = (ymomsyn(iii) + ymomadd(iii))/totalmass
                  vzi = (zmomsyn(iii) + zmomadd(iii))/totalmass
                  f1vxtemp = f1vxtemp + pmassj*f1vx(j)
                  f1vytemp = f1vytemp + pmassj*f1vy(j)
                  f1vztemp = f1vztemp + pmassj*f1vz(j)
                  f2vxtemp = f2vxtemp + pmassj*f2vx(j)
                  f2vytemp = f2vytemp + pmassj*f2vy(j)
                  f2vztemp = f2vztemp + pmassj*f2vz(j)
                  d2vxtemp = d2vxtemp + pmassj*dum2vx(j)
                  d2vytemp = d2vytemp + pmassj*dum2vy(j)
                  d2vztemp = d2vztemp + pmassj*dum2vz(j)
                  utemp = utemp + pmassj*u(j)
                  pmass(j) = pmassjn

 888              CONTINUE
               ENDIF
            ENDIF
         END DO
         pmass(i) = ptmsyn(iii) + ptmadd(iii)
         spinx(iii) = spinx(iii) + spinxtemp
         spiny(iii) = spiny(iii) + spinytemp
         spinz(iii) = spinz(iii) + spinztemp
         spinadx(iii) = spinadx(iii) + spinxtemp
         spinady(iii) = spinady(iii) + spinytemp
         spinadz(iii) = spinadz(iii) + spinztemp

         pmassnew = pmass(i)
         vx(i) = (xmomsyn(iii) + xmomadd(iii))/pmassnew
         vy(i) = (ymomsyn(iii) + ymomadd(iii))/pmassnew
         vz(i) = (zmomsyn(iii) + zmomadd(iii))/pmassnew

         x(i) = (xiold*pmassi + xtemp)/pmassnew
         y(i) = (yiold*pmassi + ytemp)/pmassnew
         z(i) = (ziold*pmassi + ztemp)/pmassnew
         f1vx(i) = (f1vxi*pmassi + f1vxtemp)/pmassnew
         f1vy(i) = (f1vyi*pmassi + f1vytemp)/pmassnew
         f1vz(i) = (f1vzi*pmassi + f1vztemp)/pmassnew
         f2vx(i) = (f2vxi*pmassi + f2vxtemp)/pmassnew
         f2vy(i) = (f2vyi*pmassi + f2vytemp)/pmassnew
         f2vz(i) = (f2vzi*pmassi + f2vztemp)/pmassnew
         dum2vx(i) = (dum2vxi*pmassi + d2vxtemp)/pmassnew
         dum2vy(i) = (dum2vyi*pmassi + d2vytemp)/pmassnew
         dum2vz(i) = (dum2vzi*pmassi + d2vztemp)/pmassnew
         u(i) = (ui*pmassi + utemp)/pmassnew
         h(i) = hi
      END DO
!
!--MERGER OF TWO POINT MASSES
!
      imerge = 0

      GOTO 1000

      DO iii = 1, nptmass
         iptm1 = listpm(iii)
         x1 = x(iptm1)
         y1 = y(iptm1)
         z1 = z(iptm1)
         DO jjj = iii + 1, nptmass
            iptm2 = listpm(jjj)
            rx = x1 - x(iptm2)
            ry = y1 - y(iptm2)
            rz = z1 - z(iptm2)
            r2 = rx*rx + ry*ry + rz*rz
            IF (r2.LT.(MAX(h(iptm1),h(iptm2))**2.0)) THEN
               imerge = 1
              if(myrank.eq.0)then
               WRITE(iprint,*) 'POINT MASSES MERGED, TIME=',realtime
               WRITE(iprint,*) '   Radius**2 ',r2
               WRITE(iprint,*) '   Point masses ',iii,jjj,iptm1,iptm2
              endif
               i1list = iii
               i2list = jjj
               IF (pmass(iptm2).GT.pmass(iptm1)) THEN
                  itemp = iptm1
                  iptm1 = iptm2
                  iptm2 = itemp
                  itemp = i1list
                  i1list = i2list
                  i2list = itemp
                  rx = -rx
                  ry = -ry
                  rz = -rz
                  x1 = x(iptm1)
                  y1 = y(iptm1)
                  z1 = z(iptm1)
               ENDIF
               iphase(iptm2) = -1
               numberacc(i1list) = numberacc(i1list) + 1
               pmass1 = pmass(iptm1)
               pmass2 = pmass(iptm2)
               if(myrank.eq.0)then
               WRITE(iprint,*) '   Point masses ',iii,jjj,iptm1,iptm2
               WRITE(iprint,*) '   Mases ',pmass1, pmass2
               endif
               dvx = x(iptm1) - x(iptm2)
               dvy = y(iptm1) - y(iptm2)
               dvz = z(iptm1) - z(iptm2)
               vx1 = vx(iptm1)
               vy1 = vy(iptm1)
               vz1 = vz(iptm1)

               totalmass = pmass1 + pmass2
               spinx(i1list) = spinx(i1list) + spinx(i2list)
               spiny(i1list) = spiny(i1list) + spiny(i2list)
               spinz(i1list) = spinz(i1list) + spinz(i2list)
               spinadx(i1list) = spinadx(i1list) + spinadx(i2list)
               spinady(i1list) = spinady(i1list) + spinady(i2list)
               spinadz(i1list) = spinadz(i1list) + spinadz(i2list)
               spinm = pmass2*pmass1/totalmass
               spinx(i1list) = spinx(i1list) + spinm*(ry*dvz - dvy*rz)
               spiny(i1list) = spiny(i1list) + spinm*(dvx*rz - rx*dvz)
               spinz(i1list) = spinz(i1list) + spinm*(rx*dvy - dvx*ry)
               spinadx(i1list) = spinadx(i1list) + &
                    spinm*(ry*dvz - dvy*rz)
               spinady(i1list) = spinady(i1list) + &
                    spinm*(dvx*rz - rx*dvz)
               spinadz(i1list) = spinadz(i1list) + &
                    spinm*(rx*dvy - dvx*ry)
               vx(iptm1) = (pmass1*vx(iptm1) + &
                    pmass2*vx(iptm2))/totalmass
               vy(iptm1) = (pmass1*vy(iptm1) + &
                    pmass2*vy(iptm2))/totalmass
               vz(iptm1) = (pmass1*vz(iptm1) + &
                    pmass2*vz(iptm2))/totalmass
               x(iptm1) = (pmass1*x(iptm1)+pmass2*x(iptm2))/totalmass
               y(iptm1) = (pmass1*y(iptm1)+pmass2*y(iptm2))/totalmass
               z(iptm1) = (pmass1*z(iptm1)+pmass2*z(iptm2))/totalmass
               f1vx(iptm1) = (pmass1*f1vx(iptm1) + &
                    pmass2*f1vx(iptm2))/totalmass
               f1vy(iptm1) = (pmass1*f1vy(iptm1) + &
                    pmass2*f1vy(iptm2))/totalmass
               f1vz(iptm1) = (pmass1*f1vz(iptm1) + &
                    pmass2*f1vz(iptm2))/totalmass
               f2vx(iptm1) = (pmass1*f2vx(iptm1) + &
                    pmass2*f2vx(iptm2))/totalmass
               f2vy(iptm1) = (pmass1*f2vy(iptm1) + &
                    pmass2*f2vy(iptm2))/totalmass
               f2vz(iptm1) = (pmass1*f2vz(iptm1) + &
                    pmass2*f2vz(iptm2))/totalmass
               dum2vx(iptm1) = (pmass1*dum2vx(iptm1) + &
                    pmass2*dum2vx(iptm2))/totalmass
               dum2vy(iptm1) = (pmass1*dum2vy(iptm1) + &
                    pmass2*dum2vy(iptm2))/totalmass
               dum2vz(iptm1) = (pmass1*dum2vz(iptm1) + &
                    pmass2*dum2vz(iptm2))/totalmass
               u(iptm1) = (pmass1*u(iptm1)+pmass2*u(iptm2))/totalmass
               numberacc(i1list) = numberacc(i1list)+numberacc(i2list)
               ptminner(i1list) = ptminner(i1list) + ptminner(i2list)

               ptmsyn(i1list) = ptmsyn(i1list) + ptmsyn(i2list)
               ptmadd(i1list) = ptmadd(i1list) + ptmadd(i2list)
               pmass(iptm1) = ptmsyn(i1list) + ptmadd(i1list)
               h(iptm1) = MAX(h(iptm1),h(iptm2))

               xmomsyn(i1list) = pmass(iptm1)*vx(iptm1)
               ymomsyn(i1list) = pmass(iptm1)*vy(iptm1)
               zmomsyn(i1list) = pmass(iptm1)*vz(iptm1)
               xmomadd(i1list) = 0.0
               ymomadd(i1list) = 0.0
               zmomadd(i1list) = 0.0

               GOTO 1000
            ENDIF
         END DO
      END DO
!
!--Compactify list of point masses if merger
!
 1000 IF (imerge.EQ.1) THEN
         icount = 0
         iaccr = 1
         DO iii = 1, nptmass
            iptm = listpm(iii)
            IF (iphase(iptm).GE.1) THEN
               icount = icount + 1
               listpm(icount) = listpm(iii)
               nactotal(icount) = nactotal(iii)
               nghtotal(icount) = nghtotal(iii)
               ptmassinner(icount) = ptmassinner(iii)
               spinx(icount) = spinx(iii)
               spiny(icount) = spiny(iii)
               spinz(icount) = spinz(iii)
               spinadx(icount) = spinadx(iii)
               spinady(icount) = spinady(iii)
               spinadz(icount) = spinadz(iii)

               ptmsyn(icount) = ptmsyn(iii)
               ptmadd(icount) = ptmadd(iii)
               xmomsyn(icount) = xmomsyn(iii)
               ymomsyn(icount) = ymomsyn(iii)
               zmomsyn(icount) = zmomsyn(iii)
               xmomadd(icount) = xmomadd(iii)
               ymomadd(icount) = ymomadd(iii)
               zmomadd(icount) = zmomadd(iii)

               nptlist(icount) = nptlist(iii)
               numberacc(icount) = numberacc(iii)
               ptminner(icount) = ptminner(iii)
               DO jjj = 1, nptlist(icount)
                  nearpt(icount,jjj) = nearpt(iii,jjj)
               END DO
            ENDIF
         END DO
         nptmass = icount
      ENDIF
!
!--Consider accreted particles' effect on ghosts
!
      DO i = 1, nghost
         j = ireal(i + npart)
         IF (iremove(j).EQ.1 .OR. &
              (icreate.EQ.1 .AND. j.EQ.irhonex)) THEN
            iphase(npart + i) = -1
         ENDIF
      END DO
!
!--Reset each type of particle count
!
      DO i = 1, nptmass
         nactotal(i) = nactotal(i) + numberacc(i)
         naccrete = naccrete + numberacc(i)
         ptmassinner(i) = ptmassinner(i) + ptminner(i)
         nactive = nactive - numberacc(i)
         nactotx(i) = nactotx(i) + numberacc(i)
         ptmassinx(i) = ptmassinx(i) + ptminner(i)
      END DO
!
!--Dump point mass details to ptprint file
!
      IF (isave.EQ.1 .OR. icreate.EQ.1) THEN
         tcomp = SQRT((3 * pi) / (32* rhozerox))
!!         tcomp = SQRT((3 * pi) / (32 * rhozero))
         if(myrank.eq.0) &
         WRITE (iptprint) realtime/tcomp, realtime, nptmass

         DO i = 1, nptmass
            j = listpm(i)
            if(myrank.eq.0)then
            IF (ibound.EQ.8 .OR. ibound.GE.90) THEN
               WRITE (iptprint)j,x(j),y(j),z(j),vx(j),vy(j), &
                    vz(j),pmass(j),rho(j),nactotal(i),ptmassinner(i), &
                    spinx(i),spiny(i),spinz(i),angaddx(i),angaddy(i), &
                    angaddz(i),spinadx(i),spinady(i),spinadz(i), &
                    naccrete,anglostx,anglosty,anglostz,nkill
            ELSE
               WRITE (iptprint)j,x(j),y(j),z(j),vx(j),vy(j), &
                    vz(j),pmass(j),rho(j),nactotal(i),ptmassinner(i), &
                    spinx(i),spiny(i),spinz(i),angaddx(i),angaddy(i), &
                    angaddz(i),spinadx(i),spinady(i),spinadz(i), &
                    naccrete
            ENDIF
            endif
            nactotal(i) = 0
            ptmassinner(i) = 0.
         END DO
      ENDIF

      IF (icreate.EQ.1) THEN
         iaccr = 1
         icreate = 0
      ENDIF

      END SUBROUTINE accrete
