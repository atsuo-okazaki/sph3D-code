      SUBROUTINE integs
!************************************************************
!                                                           *
!  This subroutine integrates the systeme of differential   *
!     equations over one timestep and measures time neded.  *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use part
      use densi
      use kerne
      use ghost
      use integ
      use ener
      use tming
      use typef
      use gtime
      use logun
      use stepopt
      use init
      use timei
      use outneigh
      use useles
      use phase
      use active
      use g3monit
      use neighbor
      use perform
      use secret
      use ptmass
      use zzhp

      implicit none

      INTEGER(I4B) :: nsteplist(30), j, ipos
      REAL(DP), save :: tafter,tbefor
      REAL(DP) :: dt, dttime
      CHARACTER(len=7) :: where='integs'

      ifail = 0
      ioutinf = 0
      ioutsup = 0
      ioutmin = 0
      ioutmax = 0
      inmin = 10000
      inmax = -1
      inminsy = 10000
      inmaxsy = 0
      IF (igrape.EQ.1) THEN
         searchfacmax = 0.
         searchfacmin = 1.0E+10
         neighover = 0
         neighretry = 0
      ENDIF
!
      IF (ncount.EQ.0) THEN
!--   Get the intial time
!
         CALL getused(tbefor) !&
      ENDIF
!
!--Advance system one hydro time-step
!
      dt = dtmax
      CALL step(dt)
      gtdouble = gtdouble + DBLE(dtmax)
!c      gt = REAL(gtdouble)
      gt = gtdouble

      IF (MOD(ncount, nprout).EQ.0 .OR. &
          MOD(ncount, nstep).EQ.0) THEN
      if(myrank.eq.0)then
      IF (ioutinf.NE.0) THEN
         WRITE (iprint, *) 'h too small', ioutinf, ' times'
         WRITE (iprint, *) 'minimum no. of neigh ', inmin
      ENDIF
      WRITE (iprint, *) 'minimum no. neigh for sync step', inminsy, &
           ' and ilen(i) was', ilenmin, ' with h of ',hilenmin
      WRITE (iprint, *) 'maximum no. neigh for sync step', inmaxsy
      IF (ioutsup.NE.0) THEN
         WRITE (iprint, *) 'h too big  ', ioutsup, ' times'
         WRITE (iprint, *) 'maximum no. of neigh (not hmin) ', inmax
      ENDIF
      IF (ioutmin.NE.0) &
           WRITE (iprint, *) 'h less than ', hmin, ioutmin,' times'
      IF (hmaximum.NE.0.0 .AND. ioutmax.NE.0) &
           WRITE (iprint, *) 'h greater than ', hmaximum, ioutmax, &
              ' times'
!
!--If GRAPE, report number of neighbour lists overflows and retries
!
      IF (igrape.EQ.1) THEN
         WRITE (iprint, 99008) searchfacmax
         WRITE (iprint, 99009) searchfacmin
         WRITE (iprint, 99010) neighover
         WRITE (iprint, 99011) neighretry
         WRITE (iprint, 99012) REAL(ilenlist)/REAL(nactive-nptmass)
         WRITE (iprint, 99014) REAL(nneightot)/REAL(nactive-nptmass), &
              REAL(internum)/REAL(nactive-nptmass)
         WRITE (iprint, 99015) tkeep

99008    FORMAT('Maximum value of the search factor       = ',F8.4)
99009    FORMAT('Minimum value of the search factor       = ',F8.4)
99010    FORMAT('Number of GRAPE-neighbour-list overflows = ', I6)
99011    FORMAT('Number of GRAPE-neighbour-list retries   = ', I6)
99012    FORMAT('Mean number of GRAPE-list neighbours     = ',F9.4)
99014    FORMAT('Mean number of nneigh(i) = ',F8.4, &
              ', Mean number of interactions = ',F8.4)
99015    FORMAT('Keep GRAPE for ',F6.1,' seconds')
      ELSE
         WRITE (iprint, 99016) REAL(ilenlist)/REAL(nactive-nptmass)
         WRITE (iprint, 99018) REAL(nneightot)/REAL(nactive-nptmass)

99016    FORMAT('Mean number of TREE-list neighbours = ',F9.4)
99018    FORMAT('Mean number of nneigh(i) = ',F8.4)
      ENDIF
      endif
      ENDIF
!
!--Report number of particles with each timestep
!
      DO j = 1, 30
         nsteplist(j) = 0
      END DO

      DO j = 1, npart
         IF (iphase(j).GE.0) THEN
            dttime = dtmax*isteps(j)/imaxstep
            ipos = INT((LOG10(dtmax/dttime) * 3.3219281) + 1.01)
            IF ((ipos.GE.1) .AND. (ipos.LE.29)) THEN
               nsteplist(ipos) = nsteplist(ipos) + 1
            ELSE
               nsteplist(30) = nsteplist(30) + 1
            ENDIF
         ENDIF
      END DO

      IF (MOD(ncount, nprout).EQ.0 .OR. &
          MOD(ncount, nstep).EQ.0) THEN
      if(myrank.eq.0)then
      WRITE (iprint,99100) npart, nactive
99100 FORMAT ('Timestep distribution ',I8,' number active ',I8)

      DO j = 0, 29
         WRITE (iprint,99101) j, dtmax/2**j, nsteplist(j+1)
      END DO
99101 FORMAT (I2,1PE12.4,I8)
      endif
!
!--Get time and compute time needed for this timestep
!
      CALL getused(tafter) !&
      tstep = (tafter - tbefor)
      tbefor = tafter
!
!--Check for messages
!
      CALL mesop
!
!--Check for errors during integration
!
      IF (ifail.NE.0) CALL error(where, ifail)
!
!--Write time, and number of ghosts
!
      if(myrank.eq.0)WRITE (iprint, 98001) tstep/60., nghost
98001 FORMAT (' time used  : ', F8.3, ' min. , ghosts : ', I8)

      IF (itiming) THEN
         IF (igrape.EQ.1) THEN
            if(myrank.eq.0)then
            WRITE (iprint, 98002)
            WRITE (iprint, 98008) tdens/60.,   tdens/tstep*100.
            WRITE (iprint, 98009) tforce/60.,  tforce/tstep*100.
            WRITE (iprint, 98016) tgforpt/60., tgforpt/tstep*100.
            WRITE (iprint, 98003) tins/60.,    tins/tstep*100.
            WRITE (iprint, 98004) tginit/60.,  tginit/tins*100.
            WRITE (iprint, 98005) tgsort/60.,  tgsort/tins*100.
            WRITE (iprint, 98006) tgload/60.,  tgload/tins*100.
            WRITE (iprint, 98007) tgcall1/60., tgcall1/tins*100.
            WRITE (iprint, 98010) tggrav/60.,  tggrav/tins*100.
            WRITE (iprint, 98012) tgcall2/60., tgcall2/tins*100.
            WRITE (iprint, 98011) tgnei/60.,   tgnei/tins*100.
            WRITE (iprint, 98015) tgnmisc/60., tgnmisc/tins*100.
            WRITE (iprint, 98014) tgnstor/60., tgnstor/tins*100.
            WRITE (iprint, 98013) tgcall3/60., tgcall3/tins*100.
            WRITE (iprint,*) 'ninit = ',ninit1,ninit2,ninit3,ninit4, &
                           ninit5
!            WRITE (iprint,*) 'ttest = ',ttest/60., ttest/tins*100.
            endif
            tdens = 0.
            tforce = 0.
            tgforpt = 0.
            tins = 0.
            tginit = 0.
            tgsort = 0.
            tgload = 0.
            tgcall1 = 0.
            tggrav = 0.
            tgcall2 = 0.
            tgnei = 0.
            tgnmisc = 0.
            tgnstor = 0.
            tgcall3 = 0.

            ninit1 = 0
            ninit2 = 0
            ninit3 = 0
            ninit4 = 0
            ninit5 = 0
            ttest = 0.

98002       FORMAT (' GRAPE Timing:')
98008       FORMAT (' densityi     time: ', F8.3, ' min. ',F5.2,' %')
98009       FORMAT (' forcei       time: ', F8.3, ' min. ',F5.2,' %')
98016       FORMAT (' gforspt      time: ', F8.3, ' min. ',F5.2,' %')
98003       FORMAT (' insulate     time: ', F8.3, ' min. ',F5.2,' %')
98004       FORMAT ('  GRP init    time: ', F8.3, ' min. ',F5.2,' %')
98005       FORMAT ('  GRP sort    time: ', F8.3, ' min. ',F5.2,' %')
98006       FORMAT ('  GRP load    time: ', F8.3, ' min. ',F5.2,' %')
98007       FORMAT ('  GRP call1   time: ', F8.3, ' min. ',F5.2,' %')
98010       FORMAT ('  GRP gravity time: ', F8.3, ' min. ',F5.2,' %')
98011       FORMAT ('  GRP nei get time: ', F8.3, ' min. ',F5.2,' %')
98012       FORMAT ('  GRP call2   time: ', F8.3, ' min. ',F5.2,' %')
98013       FORMAT ('  GRP call3   time: ', F8.3, ' min. ',F5.2,' %')
98014       FORMAT ('  GRP nstore  time: ', F8.3, ' min. ',F5.2,' %')
98015       FORMAT ('  GRP nmisc   time: ', F8.3, ' min. ',F5.2,' %')
98017       FORMAT ('  TRE mtree   time: ', F8.3, ' min. ',F5.2,' %')
98018       FORMAT ('  TRE treef   time: ', F8.3, ' min. ',F5.2,' %')
98019       FORMAT ('  TRE revtree time: ', F8.3, ' min. ',F5.2,' %')
         ELSE
            if(myrank.eq.0)then
            WRITE (iprint, 98020)
            WRITE (iprint, 98008) tdens/60.,   tdens/tstep*100.
            WRITE (iprint, 98009) tforce/60.,  tforce/tstep*100.
            WRITE (iprint, 98016) tgforpt/60., tgforpt/tstep*100.
            WRITE (iprint, 98003) tins/60.,    tins/tstep*100.
            WRITE (iprint, 98017) tmtree/60.,  tmtree/tins*100.
            WRITE (iprint, 98018) ttreef/60.,  ttreef/tins*100.
            WRITE (iprint, 98019) trevt/60.,   trevt/tins*100.
            endif
            tdens = 0.
            tforce = 0.
            tgforpt = 0.
            tins = 0.
            tmtree = 0.
            ttreef = 0.
            trevt = 0.

98020       FORMAT (' TREE Timing:')
         ENDIF
      ENDIF
      ENDIF

      END SUBROUTINE integs
