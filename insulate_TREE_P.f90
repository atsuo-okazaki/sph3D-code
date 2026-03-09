      SUBROUTINE insulate(itype, ntot, npart, x, y, z, pmass, h)
!************************************************************
!                                                           *
!  The purpose of this subroutine is to insulate the main   *
!     SPH code from being able to tell whether it uses a    *
!     binary tree or the GRAPE to get gravity forces and    *
!     neighbours.                                           *
!                                                           *
!     Written by M R Bate 2/2/96                            *
!                                                           *
!************************************************************

      use idims

      use btree
      use gravi
      use ener
      use typef
      use curlist
      use phase
      use perform
      use useles
!---- neighbor.mod <-- neighbor.f90 for non-parallel simulation
!---- neighbor.mod <-- neighbor_P.f90 for parallel simulation
      use neighbor
!---- treecom.mod <-- treecom.f90 for non-parallel simulation
!---- treecom.mod <-- treecom_P.f90 for parallel simulation
      use treecom
      use logun
      use zzhp

      implicit none

      INTEGER(I4B) :: itype, ntot, npart, i, ipart, ichg
      REAL(DP) :: x(idim), y(idim), z(idim), pmass(idim), h(idim)
      REAL(DP) :: epot, fsx, fsy, fsz
      CHARACTER(len=7) :: where='insul'

!!      WRITE (iprint,*) 'itype=',itype,', ntot=',ntot,
!!     &       ', npart=',npart

      IF (itiming) CALL getused(tins1)

!------------------------------------------------------------
!--Call mtree to make the binary tree
!------------------------------------------------------------

      IF (itype.EQ.1) THEN

         IF (itiming) CALL getused(tmtree1)

!!         WRITE (iprint,*) 'x(1)=',x(1),', y(1)=',y(1),', z(1)=',z(1)
         CALL mtree(ntot, npart, x, y, z, pmass, h)

         IF (itiming) THEN
            CALL getused(tmtree2)
            tmtree = tmtree + (tmtree2 - tmtree1)
         ENDIF

!------------------------------------------------------------
!--Call revtree to revise the binary tree
!------------------------------------------------------------

      ELSEIF (itype.EQ.2) THEN

         IF (itiming) CALL getused(trevt1)

         CALL revtree(ntot, npart, x, y, z, pmass, h)

         IF (itiming) THEN
            CALL getused(trevt2)
            trevt = trevt + (trevt2 - trevt1)
         ENDIF

!------------------------------------------------------------
!--Get gravity forces and neighbours using the binary tree
!------------------------------------------------------------

      ELSEIF (itype.EQ.3) THEN
         IF (itiming) CALL getused(ttreef1)
!
!--Compute the neighbour indexes & gravitational forces of the distant
!     particles for all the particles in the list
!
 110     CONTINUE
!$OMP PARALLEL default(none), shared(nlst,npart,acc,igphi) &
!$OMP shared(llist,h,gravx,gravy,gravz,poten,dphit,iphase) &
!$OMP private(i,ipart,fsx,fsy,fsz,epot)
!$OMP DO SCHEDULE(runtime)
         DO 100 i = 1, nlst
            ipart = llist(i)
            IF (iphase(ipart).GE.1 .AND. (.NOT.(iptintree))) GOTO 100
!
!--Walk through the tree to get the neighbours, the acceleration
!     and the potential due to outside 2h particles
!
            CALL treef(ipart,npart,h,acc,igphi,fsx,fsy,fsz,epot)

            gravx(ipart) = fsx
            gravy(ipart) = fsy
            gravz(ipart) = fsz
            IF (igphi.NE.0) THEN
               poten(ipart) = epot
            ELSE
               poten(ipart) = 0.
            ENDIF
            dphit(ipart) = 0.

 100     CONTINUE
!$OMP END DO
!$OMP END PARALLEL

!---------------------------------------------------------------
!        Reduce h(i) if the number of neighbors of any particle
!        exceeds nlmax, and execute the above DO loop again
!        (A. Okazaki, 16 Feb 2006)
!---------------------------------------------------------------
         ichg = 0
         DO i=1,nlst
            ipart = llist(i)
            IF (ilen(ipart).GT.nlmax) THEN
               IF (h(ipart).LT.0.0) THEN
                  WRITE(iprint,*) 'Error in h: h(',ipart,')=', &
                                  h(ipart)
                  h(ipart) = hmaximum
                  WRITE(iprint,*) 'h is set to',hmaximum
               ELSE
                  WRITE (iprint,1000) ilen(ipart),ipart,nlmax
 1000             FORMAT (1H ,'NO. OF NEIGH (',I8,') for ipart =', &
                          I8,' EXCEEDED NLMAX (',I8,')!')
                  WRITE (iprint,1010) x(ipart),y(ipart), &
                         z(ipart),h(ipart)
 1010             FORMAT (1H ,'(x=',1PE14.5,', y=',1PE14.5, &
                          ', z=',1PE14.5,', h=',1PE17.8,')')
                  h(ipart) = h(ipart)*0.5
                  WRITE (iprint,1020) h(ipart)
 1020             FORMAT (1H ,'h is reduced to:',1PE17.8)
               ENDIF
               ichg = ichg+1
            ENDIF
         ENDDO
         IF (ichg.NE.0) GOTO 110

         IF (itiming) THEN
            CALL getused(ttreef2)
            ttreef = ttreef + (ttreef2 - ttreef1)
         ENDIF

!------------------------------------------------------------
!--Get gravity forces and neighbours using the binary tree
!------------------------------------------------------------

      ELSEIF (itype.EQ.5) THEN
         IF (itiming) CALL getused(ttreef1)
!
!--Compute the neighbour indexes & gravitational forces of the distant
!     particles for all the particles in the list
!
 120     CONTINUE
!$OMP PARALLEL default(none), shared(nlst,npart,acc,igphi) &
!$OMP shared(llist,h,iphase) &
!$OMP private(i,ipart,fsx,fsy,fsz,epot)
!$OMP DO SCHEDULE(runtime)
         DO 200 i = 1, nlst
            ipart = llist(i)
            IF (iphase(ipart).GE.1 .AND. (.NOT.(iptintree))) GOTO 200
!
!--Walk through the tree to get the neighbours, the acceleration
!     and the potential due to outside 2h particles
!
            CALL treef(ipart,npart,h,acc,igphi,fsx,fsy,fsz,epot)

 200     CONTINUE
!$OMP END DO
!$OMP END PARALLEL

!---------------------------------------------------------------
!        Reduce h(i) if the number of neighbors of any particle
!        exceeds nlmax, and execute the above DO loop again
!        (A. Okazaki, 16 Feb, 2006)
!---------------------------------------------------------------
         ichg = 0
         DO i=1,nlst
            ipart = llist(i)
            IF (ilen(ipart).GT.nlmax) THEN
               IF (h(ipart).LT.0.0) THEN
                  WRITE(iprint,*) 'Error in h: h(',ipart,')=', &
                                  h(ipart)
                  h(ipart) = hmaximum
                  WRITE(iprint,*) 'h is set to',hmaximum
               ELSE
                  WRITE (iprint,1000) ilen(ipart),ipart,nlmax
                  WRITE (iprint,1010) x(ipart),y(ipart), &
                         z(ipart),h(ipart)
                  h(ipart) = h(ipart)*0.5
                  WRITE (iprint,1020) h(ipart)
               ENDIF
               ichg = ichg+1
            ENDIF
         ENDDO
         IF (ichg.NE.0) GOTO 120

      ELSE
         CALL error(where,2)
      ENDIF

      IF (itiming) THEN
         CALL getused(tins2)
         tins = tins + (tins2 - tins1)
      ENDIF

      END SUBROUTINE insulate
