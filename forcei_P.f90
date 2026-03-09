      SUBROUTINE forcei(jneigh,nlst_in,nlst_end,list,dt,itime, &
            npart,x,y,z,vx,vy,vz,u,h,delvx,delvy,delvz,delhp,delup, &
            trho,pr,vsound)
!************************************************************
!                                                           *
!  This subroutine computes the forces on particle ipart    *
!                                                           *
!************************************************************

      use idims

      use constants
      use table
      use force, only: fx, fy, fz, du, dh
      use tlist
      use btree
      use carac
      use numpa
      use gravi
      use ener
      use kerne
      use tming
      use typef
      use timei
      use logun
      use gtime
      use debug
      use divve
      use cgas
!---- neighbor.mod <-- neighbor.f90 for non-parallel simulation
!---- neighbor.mod <-- neighbor_P.f90 for parallel simulation
      use neighbor
      use dissi
      use current
      use useles
      use densi
      use rbnd
      use xforce
      use phase
      use ptbin
      use ptmass
      use debugpt
      use isnpt
      use glrho
      use vbound
      use artvb
      use ian
      use visc
      use pres
      use soft
      use call
      use outneigh
      use units
      use new
      use diskbd
      use sphcom
      use misali

      implicit none

      INTEGER(I4B) :: list(idim)
      INTEGER(I4B) :: jneigh, nlst_in, nlst_end, itime, npart, n, ipart, &
              numneigh, iptn, iptcur, nlocalneigh, j, k, &
              index, index1, ii, iptcurv, isphcom2, ipt, ilenipart, &
              indw, indw1
      REAL(DP) :: x(idim), y(idim), z(idim)
      REAL(DP) :: vx(idim), vy(idim), vz(idim), u(idim), h(idim)
      REAL(DP) :: trho(idim), pr(idim), vsound(idim)
      REAL(DP) :: delvx(idim), delvy(idim), delvz(idim)
      REAL(DP) :: delhp(idim), delup(idim)
      REAL(DP) :: dt, rmax2, stepsi, xi, yi, zi, vxi, vyi, vzi, &
              pmassi, dhi, hi, gravxi, gravyi, gravzi, poteni, &
              dphiti, gradxi, gradyi, gradzi, artxi, artyi, artzi, &
              pdvi, dqi, rhoi, pro2i, vsoundi, xipt, yipt, zipt, &
              vxipt, vyipt, vzipt, xip, yip, zip, hj, dx, dy, dz, &
              rij, rij1, rij2, pmassj, runix, runiy, runiz, &
              hmean, hmean21, hmean41, dhmean, v, v2, dxx, rij2grav, &
              fm, phi, dphi, dfmassdx, dfptdx, dpotdh, xmasj, &
              xii, yii, zii, vpos, rhoj, robar, dgrwdx, grwtij, grpm, &
              poro2j, dvx, dvy, dvz, projv, vsbar, f, adivi, acurlvi, &
              fi, adivj, acurlvj, fj, t12j, alphafac, ri1, ri2, &
              rj1, rj2, hcoeff, rangmx, hdiski, hdiskj, hfac, &
              vlowcorrection, qi, qj, pxtemp, pytemp, pztemp, &
              artxtemp, artytemp, artztemp, spinxtemp, spinytemp, spinztemp, &
              rr, radpgrad, prgradrhoj, acceltrue, vsmvfor, weight, &
              vsmvfor2, dxxw, dwdxw, acceladd, vnormj, vnormpt, &
              tanvx1j, tanvy1j, tanvz1j, tvf1, acc1diff, addx, addy, addz, &
              addspinx, addspiny, addspinz, pmassrat, stepsipt, addspin, &
              realtime, thermal
      REAL(DP), parameter :: epsil=1.E-2
      REAL(DP), parameter :: epsil2=1.E-4
      CHARACTER(len=7) :: where='forcei'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 250)
  250 FORMAT(' entry subroutine forcei')
!
!--Initialize
!
      rmax2 = rmax*rmax
      IF (ibound == 95 .OR. ibound == 97) THEN
         thermal = thermal3
      ELSE
         thermal = thermal1
      END IF

!$OMP PARALLEL default(none) &
!$OMP shared(nlst_in,nlst_end,npart,hmin,list,nneigh,neimin,neimax) &
!$OMP shared(icall,dt,imaxstep,isteps) &
!$OMP shared(h,gradpx,gradpy,gradpz,artvix,artviy,artviz,pdv,dq) &
!$OMP shared(x,y,z,vx,vy,vz,pmass,dh,u,trho,pr,vsound) &
!$OMP shared(ilen,neighb,dvtable,igphi,psoft,radkernel) &
!$OMP shared(fmass,fpoten,dphidh,part2kernel,part1kernel) &
!$OMP shared(part2potenkernel,part1potenkernel,grwij,ifsvi) &
!$OMP shared(isphcom,udist,utime,uergg,gamma,thermal) &
!$OMP shared(divv,curlv,beta,alpha,gravx,gravy,gravz,poten,dphit) &
!$OMP shared(cnormk,where,pext) &
!$OMP shared(iphase,isnearpt,listpm,nptmass) &
!$OMP shared(rangle,alphasph,betasph,encal,iaccevol) &
!$OMP private(n,ipart,stepsi,numneigh) &
!$OMP private(xi,yi,zi,vxi,vyi,vzi,pmassi,dhi,hi,gravxi,gravyi,gravzi) &
!$OMP private(poteni,dphiti,gradxi,gradyi,gradzi,artxi,artyi,artzi) &
!$OMP private(pdvi,dqi,rhoi,pro2i,vsoundi,k,j,hj,dx,dy,dz) &
!$OMP private(rij2,rij,rij1,pmassj,runix,runiy,runiz,hmean,hmean21) &
!$OMP private(hmean41,dhmean,v2,v,index,dxx,index1,rij2grav,fm) &
!$OMP private(phi,dphi,dfmassdx,dfptdx,dpotdh,xmasj,rhoj,robar) &
!$OMP private(dgrwdx,grwtij,grpm,poro2j,dvx,dvy,dvz,projv,vsbar) &
!$OMP private(f,adivi,acurlvi,fi,adivj,acurlvj,fj,t12j) &
!$OMP private(alphafac,hfac,hcoeff,hdiski,hdiskj,ri1,ri2,rj1,rj2) &
!$OMP private(isphcom2,vlowcorrection,qi,qj,rangmx) &
!$OMP private(iptn,iptcur,xipt,yipt,zipt,vxipt,vyipt,vzipt) &
!$OMP private(xip,yip,zip,nlocalneigh,ii,iptcurv,xii,yii,zii,vpos) &
!$OMP reduction(+:ioutmin,ioutsup,ioutinf) &
!$OMP reduction(MIN:inmin,inminsy) &
!$OMP reduction(MAX:inmax,inmaxsy)

!$OMP DO SCHEDULE(runtime)
      DO n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).EQ.-1) THEN
            WRITE(iprint,*) 'Error: Force for non-existant particle'
            CALL quit
         ENDIF
         gradpx(ipart) = 0.
         gradpy(ipart) = 0.
         gradpz(ipart) = 0.
         artvix(ipart) = 0.
         artviy(ipart) = 0.
         artviz(ipart) = 0.
         pdv(ipart) = 0.
         dq(ipart) = 0.
!
!--Derivative of smoothing length
!
         IF (icall.EQ.3) THEN
            numneigh = nneigh(ipart)
            inmin = MIN(inmin,numneigh)
            inmax = MAX(inmax,numneigh)
            inminsy = MIN(inminsy,numneigh)
            inmaxsy = MAX(inmaxsy,numneigh)
            IF (h(ipart).LT.hmin .AND. numneigh.GT.neimin) &
                 ioutmin = ioutmin + 1
            IF (numneigh.GT.neimax) ioutsup = ioutsup + 1
            IF (numneigh.LT.neimin) ioutinf = ioutinf + 1
         ENDIF
         stepsi = dt*isteps(ipart)/imaxstep
         CALL hdot(npart, ipart, stepsi, h)
      END DO
!$OMP END DO

!$OMP DO SCHEDULE(runtime)
      DO 80 n = nlst_in, nlst_end
         ipart = list(n)
         IF (iphase(ipart).GE.1) GOTO 80
!
!--Compute forces on particle ipart
!
         xi = x(ipart)
         yi = y(ipart)
         zi = z(ipart)
         vxi = vx(ipart)
         vyi = vy(ipart)
         vzi = vz(ipart)
         pmassi = pmass(ipart)
         dhi = dh(ipart)
         hi = h(ipart)

         gravxi = 0.
         gravyi = 0.
         gravzi = 0.
         poteni = 0.
         dphiti = 0.

         gradxi = 0.
         gradyi = 0.
         gradzi = 0.

         artxi = 0.
         artyi = 0.
         artzi = 0.

         pdvi = 0.
         dqi = 0.

         rhoi = trho(ipart)
         IF (iphase(ipart).GE.1) THEN
            pro2i = 0.0
         ELSE
            pro2i = (pr(ipart) - pext)/(rhoi*rhoi)
         ENDIF
         vsoundi = vsound(ipart)

         stepsi = dt*isteps(ipart)/imaxstep

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
!
!--Loop over neighbors
!
         DO 70 k = 1, ilen(ipart)
            j = neighb(k,ipart)

            IF (iphase(j).GE.1) GOTO 70

            hj = h(j)

            IF (iphase(j).EQ.-1) THEN
               WRITE(iprint,*)'ERROR - Accreted particle as neighbour!'
               WRITE(iprint,*) j,x(j),y(j),z(j),vx(j),vy(j),vz(j)
               WRITE(iprint,*) ipart,icall,xi,yi,zi
               CALL quit
            ENDIF
!
!--Gravity and potential energy
!
            dx = xi - x(j)
            dy = yi - y(j)
            dz = zi - z(j)
            rij2 = dx*dx + dy*dy + dz*dz + tiny
            rij = SQRT(rij2)
            rij1 = 1./rij
            pmassj = pmass(j)
!
!--Unit vectors
!
            runix = dx*rij1
            runiy = dy*rij1
            runiz = dz*rij1
!
!--Define mean h
!
            IF (iphase(ipart).GE.1) THEN
               hmean = hj/2.0
            ELSEIF (iphase(j).GE.1) THEN
               hmean = hi/2.0
            ELSE
               hmean = 0.5*(hi + hj)
            ENDIF
            hmean21 = 1./(hmean*hmean)
            hmean41 = hmean21*hmean21
!
!--dhmean uses old dh(j) if particle j is not currently being evaluated
!
            IF (iphase(ipart).GE.1) THEN
               dhmean = dh(j)*stepsi
            ELSE IF (iphase(j).GE.1) THEN
               dhmean = dhi*stepsi
            ELSE
               dhmean = 0.5*(dhi + dh(j))*stepsi
            ENDIF
            v2 = rij2*hmean21
            v = rij/hmean

            index = v2/dvtable
            dxx = v2 - index*dvtable
            index1 = index + 1
            IF (index1.GT.itable) index1 = itable

            IF (igrape.EQ.0 .AND. igphi.NE.0) THEN
               IF (isoft.EQ.1) THEN
                  rij2grav = dx*dx + dy*dy + dz*dz + psoft**2
                  fm = 1.0
                  phi = - 1./SQRT(rij2grav)
                  dphi = 0.0
               ELSEIF (isoft.EQ.0) THEN
                  rij2grav = rij2
                  IF (v.GE.radkernel) THEN
                     fm = 1.0
                     phi = -rij1
                     dphi = 0.0
                  ELSE
                     dfmassdx = (fmass(index1) - fmass(index))/dvtable
                     fm = (fmass(index) + dfmassdx*dxx)
                     dfptdx = (fpoten(index1) - fpoten(index))/dvtable
                     phi = (fpoten(index) + dfptdx*dxx)/hmean
                     dpotdh = (dphidh(index1) - dphidh(index))/dvtable
                     dphi = (dphidh(index) + dpotdh*dxx)*hmean21*dhmean
                     IF (v.GT.part2kernel) THEN
                        phi = phi + rij1*part2potenkernel
                     ELSEIF (v.GT.part1kernel) THEN
                        phi = phi + rij1*part1potenkernel
                     ENDIF
                  ENDIF
               ELSE
                  CALL error(where,1)
               ENDIF
!
!--Gravitational force calculation
!
               IF (j.LE.npart) THEN
                  xmasj = fm*pmassj/rij2grav
                  gravxi = gravxi - xmasj*runix
                  gravyi = gravyi - xmasj*runiy
                  gravzi = gravzi - xmasj*runiz
                  poteni = poteni + phi*pmassj
                  dphiti = dphiti + pmassj*dphi
               ENDIF
            ENDIF
!
!--Pressure and artificial viscosity
!     There is no pressure between point masses and particles
!     There is no viscosity between point masses and particles
!     There is no viscosity or pressure between two point masses
!
!--No artificial viscosity or pressure between particles across a point mass
!
            DO ii = 1, nptmass
               iptcurv = listpm(ii)
               xii = x(iptcurv)
               yii = y(iptcurv)
               zii = z(iptcurv)
               vpos = (xii-xi)*(xii-x(j)) + (yii-yi)*(yii-y(j)) + &
                    (zii-zi)*(zii-z(j))
               IF (vpos.LT.0.0) GOTO 70
            END DO

            IF(iphase(ipart).GE.1 .OR. iphase(j).GE.1) GOTO 70

            IF (v.LT.radkernel) THEN
               rhoj = trho(j)
               robar = 0.5*(rhoi + rhoj)
!
!--Get kernel quantities from interpolation in table
!
               dgrwdx = (grwij(index1) - grwij(index))/dvtable
               grwtij = (grwij(index) + dgrwdx*dxx)*hmean41
               grpm = pmassj*grwtij
!
!--Pressure gradient and pdv
!
               poro2j = grpm*(pro2i + (pr(j) - pext)/(rhoj**2))
               gradxi = gradxi + poro2j*runix
               gradyi = gradyi + poro2j*runiy
               gradzi = gradzi + poro2j*runiz

               dvx = vxi - vx(j)
               dvy = vyi - vy(j)
               dvz = vzi - vz(j)
               projv = dvx*runix + dvy*runiy + dvz*runiz

               pdvi = pdvi + grpm*projv
!
!--Artificial viscosity and energy dissipation
!
!             IF (ifsvi.NE.0 .AND. projv.LT.0.) THEN
               IF (ifsvi.NE.0 .AND. projv.LT.0. .AND. j.LE.npart) THEN
!
!--Calculate artificial viscosity:
!     If ifsvi=1 then normal viscosity
!     If ifsvi=2 then divv/curl weighted viscosity
!     If ifsvi=3 then viscosity reduced linearly to zero below vsound/2
!     If ifsvi=4 then Balsara viscosity (divv/curl weighted,but times pressure)
!     If ifsvi=5 then viscosity in Hernquist and Katz, ApJS 70, 424
!     If ifsvi=6 then constant alpha (Shakura-Sunyaev) viscosity
!
                  IF (ifsvi.NE.5) THEN
                     vsbar = 0.5*(vsoundi + vsound(j))
                     f = projv*v/(v2 + epsil)
                     IF (ifsvi.EQ.2 .OR. ifsvi.EQ.4) THEN
                        adivi = ABS(divv(ipart)/rhoi)
                        acurlvi = ABS(curlv(ipart)/rhoi)
                        fi = adivi/(adivi+acurlvi+epsil2*vsoundi/hi)
                        adivj = ABS(divv(j)/rhoj)
                        acurlvj = ABS(curlv(j)/rhoj)
                        fj = adivj/(adivj+acurlvj+epsil2*vsound(j)/hj)

                        IF (ifsvi.EQ.2) THEN
                           f = f*(fi+fj)/2.0

                           t12j = grpm*f*(beta*f - alpha*vsbar)/robar
                        ELSEIF (ifsvi.EQ.4) THEN
                           f = f*(fi+fj)/(vsoundi+vsound(j))

                           t12j = poro2j*f*(beta*f - alpha)
                        ENDIF
                     ELSEIF (ifsvi.EQ.6) THEN
!----                   The following method may be used only for either
!                       isothermal or polytropic disks in coplanar
!                       binaries. Never use it for misaligned systems.
!                       (A. Okazaki, 01/23/2007)
                        alphafac = 0.1
                        isphcom2 = 3-isphcom
                        ri1 = SQRT((xi-x(isphcom))*(xi-x(isphcom)) &
                                  +(yi-y(isphcom))*(yi-y(isphcom)) &
                                  +(zi-z(isphcom))*(zi-z(isphcom)))
                        ri2 = SQRT((xi-x(isphcom2))*(xi-x(isphcom2)) &
                                  +(yi-y(isphcom2))*(yi-y(isphcom2)) &
                                  +(zi-z(isphcom2))*(zi-z(isphcom2)))
                        rj1 = SQRT((x(j)-x(isphcom))*(x(j)-x(isphcom)) &
                                  +(y(j)-y(isphcom))*(y(j)-y(isphcom)) &
                                  +(z(j)-z(isphcom))*(z(j)-z(isphcom)))
                        rj2 = SQRT((x(j)-x(isphcom2))*(x(j)-x(isphcom2)) &
                                +(y(j)-y(isphcom2))*(y(j)-y(isphcom2)) &
                                +(z(j)-z(isphcom2))*(z(j)-z(isphcom2)))
!---- alphasph, betasph: parameters of numerical voscosity if ifsvi=6
                        IF (encal.EQ.'i') THEN
!                          By definition, SQRT(uergg)*utime/udist = 1.
!!                           hcoeff = SQRT(gamma*thermal/(3.0/2.0/uergg))
!!     &                           /(SQRT(pmass(isphcom))*udist/utime)
                           hcoeff = SQRT(gamma*thermal/(3.0/2.0) &
                                        /pmass(isphcom))
                           rangmx = MAX(ABS(rangle(1)),ABS(rangle(2)), &
                                        ABS(rangle(3)))
                           IF (iaccevol.EQ.'v' .OR. rangmx.NE.0.0) THEN
                              hdiski = hcoeff*ri1**1.5
                              hdiskj = hcoeff*rj1**1.5
                           ELSEIF (iaccevol.EQ.'f') THEN
                              hdiski = hcoeff*ri1**1.5 &
                                /SQRT(1.0+pmass(isphcom2)/pmass(isphcom) &
                                      *(ri1/ri2)**3)
                              hdiskj = hcoeff*rj1**1.5 &
                                /SQRT(1.0+pmass(isphcom2)/pmass(isphcom) &
                                      *(rj1/rj2)**3)
                           ELSE
                              WRITE (iprint,*) 'ERROR - iaccevol:', &
                                               iaccevol
                              CALL quit
                           ENDIF
                           alphasph = alpha/(alphafac*hmean &
                                          /(0.5*(hdiski+hdiskj)))
                        ELSE IF (encal.EQ.'p' .OR. encal.EQ.'a' &
                                .OR. encal.EQ.'c') THEN
!----                      The following implementation calculates
!                          the scale-height of the disk from rho,
!                          u, and z of individual particles
!                          (not the z-distribution of particles).
!                          It uses the scale-height of the disk,
!                          not the half thickness of the disk, in
!                          order to give a smooth transition from the
!                          isothermal case.
!                          Be careful to adopt this option.
!                          (A. Okazaki, 02/07/2007)
                           hfac = SQRT(1.0-EXP(-0.5*(gamma-1.0)))
                           hcoeff = 4.0*gamma/(3.0*(gamma-1.0)) &
                                    /pmass(isphcom)
                           hdiski = hfac*SQRT(hcoeff*u(ipart)*ri1**3 &
                                 /(1.0+pmass(isphcom2)/pmass(isphcom) &
                                   *(ri1/ri2)**3) &
                                +zi*zi)
                           hdiskj = hfac*SQRT(hcoeff*u(j)*rj1**3 &
                                 /(1.0+pmass(isphcom2)/pmass(isphcom) &
                                   *(rj1/rj2)**3) &
                                +z(j)*z(j))
                           alphasph = alpha/(alphafac*hmean &
                                          /(0.5*(hdiski+hdiskj)))
                        ENDIF

!---- betasph=2*alphasph is adopted if beta given initially is not zero
                        IF (beta.NE.0.0) THEN
                           betasph = 2.0*alphasph
                        ELSE
                           betasph = beta
                        ENDIF
                        t12j = grpm*f*(betasph*f - alphasph*vsbar) &
                               /robar
                     ELSE
                        t12j = grpm*f*(beta*f - alpha*vsbar)/robar
                     ENDIF
                     IF (ifsvi.EQ.3) THEN
                        vlowcorrection = ABS(projv/vsbar)
                        IF (vlowcorrection.LT.0.5) &
                             t12j = 2.0*vlowcorrection*t12j
                     ENDIF
                  ELSE
!
!--Hernquist and Katz
!
                     IF (divv(ipart).LT.0) THEN
                        adivi = ABS(divv(ipart)/rhoi)
                        qi = hi*rhoi*adivi*(alpha*vsoundi + &
                             beta*hi*adivi)
                     ELSE
                        qi = 0.0
                     ENDIF
                     IF (divv(j).LT.0) THEN
                        adivj = ABS(divv(j)/rhoj)
                        hj = h(j)
                        qj = hj*rhoj*adivj*(alpha*vsound(j) + &
                             beta*hj*adivj)
                     ELSE
                        qj = 0.0
                     ENDIF
                     t12j = grpm*(qi/(rhoi**2) + qj/(rhoj**2))
                  ENDIF

                  artxi = artxi + t12j*runix
                  artyi = artyi + t12j*runiy
                  artzi = artzi + t12j*runiz
                  dqi = dqi + t12j*projv
               ENDIF
            ENDIF
 70      CONTINUE
!
!--Store quantities
!
         IF (igrape.EQ.0) THEN
            gravx(ipart) = gravx(ipart) + gravxi
            gravy(ipart) = gravy(ipart) + gravyi
            gravz(ipart) = gravz(ipart) + gravzi
            poten(ipart) = poten(ipart) + poteni
            dphit(ipart) = dphit(ipart) + dphiti
!!            IF (ipart.LE.10) THEN
!!               WRITE (iprint,'(1h , ''forcei: poten('',i5,'')='',1pd12.4)') &
!!                  ipart,poten(ipart)
!!            ENDIF
         ENDIF
         gradpx(ipart) = gradpx(ipart) + gradxi
         gradpy(ipart) = gradpy(ipart) + gradyi
         gradpz(ipart) = gradpz(ipart) + gradzi
         artvix(ipart) = artvix(ipart) + artxi
         artviy(ipart) = artviy(ipart) + artyi
         artviz(ipart) = artviz(ipart) + artzi
         pdv(ipart) = pdv(ipart) + pdvi
         dq(ipart) = dq(ipart) + dqi
 80   CONTINUE
!$OMP END DO
!$OMP END PARALLEL

      IF (igrp.NE.0) THEN
         DO iptn = 1, nptmass
            ipt = listpm(iptn)
            IF (iphase(ipt).EQ.4 .AND. iscurrent(ipt).EQ.1) THEN
!
!--Loop over neighbours
!
               pxtemp = 0.0
               pytemp = 0.0
               pztemp = 0.0
               artxtemp = 0.0
               artytemp = 0.0
               artztemp = 0.0
               spinxtemp = 0.0
               spinytemp = 0.0
               spinztemp = 0.0

               ilenipart = nptlist(iptn)
               DO 400 k = 1, ilenipart
                  j = nearpt(iptn,k)
                  IF (isnearpt(j).NE.iptn) GOTO 400

                  dx = x(j) - x(ipt)
                  dy = y(j) - y(ipt)
                  dz = z(j) - z(ipt)
                  rr = SQRT(dx*dx + dy*dy + dz*dz)
!
!--Radial pressure gradient from point mass
!
                  radpgrad = gradpx(j)*dx+ gradpy(j)*dy+ gradpz(j)*dz
                  IF (rr.EQ.0.0) THEN
                     WRITE(iprint,*)'ERROR - Part. on ptmass'
                     CALL quit
                  ELSE
                     radpgrad = -cnormk*radpgrad/rr
                  ENDIF
!
!--Find radial pressure gradient approximation
!
                  rhoi = trho(j)
                  prgradrhoj = prgradrho(j)
                  IF (prgradrhoj .LT. -pr(j)/rhoi/h(j)) THEN
                     prgradrhoj = -pr(j)/rhoi/h(j)/2.0
                  ELSEIF (prgradrhoj .GT. pr(j)/rhoi/h(j)) THEN
                     prgradrhoj = pr(j)/rhoi/h(j)/2.0
                  ENDIF
                  acceltrue = -prgradrhoj
!
!--Smooth the correction to viscous forces outside 2*hacc
!
!!c                  htemp = MIN(h(j), h(ipt))
!!c                  vsmvfor = 4.0*(rr - h(ipt) - 1.5*htemp)/htemp
                  vsmvfor = 4.0*(rr - h(ipt) - 1.5*h(j))/h(j)
!!c                  vsmvfor = 2.0*(rr - h(ipt) - h(j))/h(j)
!!c                  vsmvfor = (rr - h(ipt))/h(j)
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

                  IF(radpgrad.LT.acceltrue) THEN
                     acceladd = (acceltrue - radpgrad)*weight
                     vnormj = 1.0/(rr*cnormk)
                     IF (iscurrent(j).EQ.1) THEN
                        gradpx(j) = gradpx(j) - acceladd*dx*vnormj
                        gradpy(j) = gradpy(j) - acceladd*dy*vnormj
                        gradpz(j) = gradpz(j) - acceladd*dz*vnormj
                     ENDIF
                     vnormpt = vnormj*pmass(j)/pmass(ipt)
                     pxtemp = pxtemp + acceladd*dx*vnormpt
                     pytemp = pytemp + acceladd*dy*vnormpt
                     pztemp = pztemp + acceladd*dz*vnormpt
                  ENDIF
!
!--Viscosity boundaries.  First, find viscous forces on particle in
!     direction of smoothed velocity of particle (from densityi) (the
!     total viscous force, and the bits from part 1 and part 2.
!     Then compare this to the expected force by extrapolating that
!     from part 2.  Finally add, or subtract, the necessary extra force
!     in that direction, and put a backward force on the point mass, to
!     conserve linear momentum, and angular momentum.
!
                  tanvx1j = tanvx1(j)
                  tanvy1j = tanvy1(j)
                  tanvz1j = tanvz1(j)

                  IF (iscurrent(j).EQ.1) THEN
                     accgot1(j) = artvix(j)*tanvx1j + &
                          artviy(j)*tanvy1j + artviz(j)*tanvz1j
                  ENDIF

                  tvf1 = tanviscfor1(j)
                  IF (tvf1.LT.0.) tvf1 = 0.0
                  acc1diff = (tvf1 - accgot1(j))*weight

                  addx = acc1diff*tanvx1j
                  addy = acc1diff*tanvy1j
                  addz = acc1diff*tanvz1j

                  addspinx = (dy*addz - addy*dz)
                  addspiny = (addx*dz - dx*addz)
                  addspinz = (dx*addy - addx*dy)
!
!--Add corrections to viscous forces
!
                  IF (iscurrent(j).EQ.1) THEN
                     artvix(j) = artvix(j) + addx
                     artviy(j) = artviy(j) + addy
                     artviz(j) = artviz(j) + addz
                  ENDIF
                  pmassrat = pmass(j)/pmass(ipt)
                  artxtemp = artxtemp - addx*pmassrat
                  artytemp = artytemp - addy*pmassrat
                  artztemp = artztemp - addz*pmassrat
                  IF (icall.EQ.3) THEN
                     stepsipt = dt*isteps(ipt)/imaxstep
                     addspin = pmass(j)*stepsipt*cnormk
                     spinxtemp = spinxtemp + addspin*addspinx
                     spinytemp = spinytemp + addspin*addspiny
                     spinztemp = spinztemp + addspin*addspinz
                  ENDIF
 400           CONTINUE

               gradpx(ipt) = gradpx(ipt) + pxtemp
               gradpy(ipt) = gradpy(ipt) + pytemp
               gradpz(ipt) = gradpz(ipt) + pztemp
               artvix(ipt) = artvix(ipt) + artxtemp
               artviy(ipt) = artviy(ipt) + artytemp
               artviz(ipt) = artviz(ipt) + artztemp
               spinx(iptn) = spinx(iptn) + spinxtemp
               spiny(iptn) = spiny(iptn) + spinytemp
               spinz(iptn) = spinz(iptn) + spinztemp
            ENDIF
         END DO
      ENDIF

!!      OPEN (30,FILE='forcetree', FORM = 'unformatted')
!!      DO i = 1, npart
!!      WRITE (30) gravx(i), gravy(i), gravz(i)
!!      END DO
!!      CLOSE(30)
!!      CALL quit

      realtime = dt*itime/imaxstep + gt

!$OMP PARALLEL default(none) &
!$OMP shared(nlst_in,nlst_end,nptmass,list,fx,fy,fz,cnormk) &
!$OMP shared(gradpx,gradpy,gradpz,igrp,igphi,ifsvi,iexf) &
!$OMP shared(ifcor,iexpan,iener,damp) &
!$OMP shared(gravx,gravy,gravz,artvix,artviy,artviz) &
!$OMP shared(x,y,z,vx,vy,vz,u,du,dh,realtime) &
!$OMP shared(delvx,delvy,delvz,delhp,delup) &
!$OMP shared(encal,dt) &
!$OMP private(n,ipart)

!$OMP DO SCHEDULE(runtime)
      DO n = nlst_in,nlst_end
         ipart = list(n)
!
!--Normalise and add forces
!
!--Pressure gradients or initialize forces
!
         IF (igrp.NE.0) THEN
            fx(ipart) = -gradpx(ipart)*cnormk
            fy(ipart) = -gradpy(ipart)*cnormk
            fz(ipart) = -gradpz(ipart)*cnormk
         ELSE
            fx(ipart) = 0.
            fy(ipart) = 0.
            fz(ipart) = 0.
         ENDIF
!
!--Gravity
!
         IF (igphi.NE.0 .OR. nptmass.NE.0 ) THEN
            fx(ipart) = fx(ipart) + gravx(ipart)
            fy(ipart) = fy(ipart) + gravy(ipart)
            fz(ipart) = fz(ipart) + gravz(ipart)
         ENDIF
!
!--Artificial viscosity
!
         IF (ifsvi.NE.0) THEN
            fx(ipart) = fx(ipart) - artvix(ipart)*cnormk
            fy(ipart) = fy(ipart) - artviy(ipart)*cnormk
            fz(ipart) = fz(ipart) - artviz(ipart)*cnormk
         ENDIF
!
!--External forces
!
         IF (iexf.GE.1) CALL externf(ipart,x,y,z,iexf)
!
!--Coriolis and centrifugal forces
!
         IF (ifcor.NE.0) CALL coriol(ipart,realtime,x,y,z,vx,vy,vz)
!
!--Homologous expansion or contraction
!
         IF (iexpan.GT.0) CALL homexp(ipart,realtime,vx,vy,vz)
!
!--Energy conservation
!
         IF (iener.NE.0) THEN
            CALL energ(ipart,realtime,u)

!----    Don't call 'radcool' here. Calling it here in the Runge-
!        Kutta-Fehlberg part is most likely incompatible with
!        Rich Townsend (2009)'s Exact Integration Scheme inplemented
!        in 'radcool'. Instead, call it in step_P.f90 after the
!        R-K-F part is finished.
!----       Correction of u(i) due to optically thin,
!           radiative cooling. This uses Townsend (2009)'s
!           Exact Integration scheme.
!           (Note that this scheme can't be used in energ.f,
!           because the EI scheme does not fit in an explicit method.
!           (A. Okazaki, 05/05/2010)
!
!!         IF (encal.EQ.'c') THEN
!!            CALL radcool (ipart, dt, u(ipart) )
!!         ENDIF
         END IF
!
!--Damp velocities if appropiate
!
         IF (damp.NE.0.) THEN
            fx(ipart) = fx(ipart) - damp*vx(ipart)
            fy(ipart) = fy(ipart) - damp*vy(ipart)
            fz(ipart) = fz(ipart) - damp*vz(ipart)
         ENDIF

!      IF (ipart.EQ.1) THEN
!         write(iprint,99200) gradpx(ipart),gradpy(ipart),gradpz(ipart)
!99200    FORMAT('pres ',1F15.10,1F15.10,1F15.10)
!         write(iprint,99201) gravx(ipart),gravy(ipart),gravz(ipart)
!99201    FORMAT('grav ',1F15.10,1F15.10,1F15.10)
!         write(iprint,99202) artvix(ipart),artviy(ipart),artviz(ipart)
!99202    FORMAT('art ',1F15.10,1F15.10,1F15.10)
!         write(iprint,99203) pdv(ipart)
!99203    FORMAT('pdv ',1F15.10)
!         write(iprint,99204) dq(ipart)
!99204    FORMAT('dq ',1F15.10)
!         write(iprint,99207) fx(ipart), fy(ipart), fz(ipart)
!99207    FORMAT('fx,fy,fz ',1F15.10,1F15.10,1F15.10)
!      END IF

         delvx(ipart) = fx(ipart)
         delvy(ipart) = fy(ipart)
         delvz(ipart) = fz(ipart)
         delhp(ipart) = dh(ipart)
         delup(ipart) = du(ipart)

      END DO
!$OMP END DO
!$OMP END PARALLEL
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 255)
 255  FORMAT(' exit subroutine forcei')

      END SUBROUTINE forcei
