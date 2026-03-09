      SUBROUTINE save
!************************************************************
!                                                           *
!  This routine determines whether or not the time has come *
!     to write a dump on disk and print out a detailed      *
!     state of the run                                      *
!                                                           *
!************************************************************

      use mpi_mod
      use tming
      use gtime
      use stop
      use btree
      use integ
      use typef
      use logun
      use debug
      use binary
      use zzhp

      implicit none

      INTEGER(I4B) :: nleft
      REAL(DP) :: tleft, tneed, tused
      CHARACTER(len=7) :: where='save'
!
!--Allow for tracing flow
!
      IF (itrace.EQ.'all') WRITE (iprint, 99001)
99001 FORMAT (' entry subroutine save')
!
!--Increment time step counter
!
      ncount = ncount + 1
      nbuild = nbuild + 1
!
!--Get remaining time for the job
!
      CALL getused(tused)
      tleft = 60.*tmax - tused
!
!--Evaluate time needed for next timestep
!
      nleft = nstep - ncount
      tneed = 1.0*tstep + 1.0
!
!--If enough time and no dump required go on with integration
!
      IF (tleft.LT.tneed) istop = 1
      IF (gt.GT.tstop) istop = 1
      IF (ibound.EQ.8 .AND. naccrete.GT.nstop) istop = 1
      IF (nleft.LE.0) THEN
         ncount = 0
         IF (ibound.EQ.8 .AND. naccrete.GT.nstop-nfastd) nstep = 1
      ENDIF
!
!--Transform into original frame of reference
!
      CALL chanref(1)
!
!--Compute and print out the state of the system
!
      CALL inform(where)
!
!--Transform into expanding frame of reference
!
      CALL chanref(2)

      IF (idebug(1:4).EQ.'save') THEN
      if(myrank.eq.0)then
         WRITE (iprint, 99002) nstep, ncount, nleft
99002    FORMAT (1X, 3(I4,1X))
         WRITE (iprint, 99003) tmax, tleft, tstep, tneed
99003    FORMAT (1X, 5(1PE12.5,1X))
      endif
      ENDIF
!
!--Check for run termination
!
!c      IF (istop.EQ.1) CALL endrun
      IF (istop.EQ.1) THEN
         if(myrank.eq.0)then
            WRITE (iprint,*) 'tleft=',tleft,', tneed=',tneed
            WRITE (iprint,*) 'gt=',gt,', tstop=',tstop
            WRITE (iprint,*) 'ibound=',ibound,', naccrete=', &
                             naccrete,', nstop=',nstop
            WRITE (iprint,*) 'tstep=',tstep,', tmax=',tmax, &
                  ', nstep=',nstep,', ncount=',ncount, &
                  ', nprout=',nprout,', tstop=',tstop
         endif
         CALL endrun
      ENDIF

      END SUBROUTINE save
