      SUBROUTINE prout (where)
!************************************************************
!                                                           *
!  This routine prints out all interesting quantities at    *
!     the present time                                      *
!                                                           *
!************************************************************

      use mpi_mod
      use idims

      use constants
      use tming
      use recor
      use gtime
      use bodys
      use ener
      use angm
      use fracg
      use cgas
      use typef
      use latti
      use expan
      use trans
      use kerne
      use out
      use new
      use files
      use units
      use logun
      use task
      use varet
      use debug
      use polyk2
      use btree
      use active
      use ptdump
      use binary
      use phase
      use ptmass
      use ptbin
      use rbnd
      use part
      use physeos
      use accnum
      use carac
      use ghost
      use crpart
      use sphcom
      use capt
      use maslos
      use winds
      use split

      implicit none

      INTEGER(I4B) :: i, ipt
      REAL(DP) :: ajeans, alph, betl, betr, tcomp, total
      CHARACTER(len=24) :: sentenc
      CHARACTER(len=7) :: where
!
!--Allow for tracing flow
!
      IF (myrank.eq.0) then

         IF (itrace == 'all') WRITE (iprint, 99001)
99001    FORMAT (' entry subroutine prout')
!
!--Check for output type
!
         IF (where (1:6)  == 'inform') THEN
!
!--Write integration output :
!
!--Number of cycle
!
!         IF (nbuild == 1 .OR. MOD(ncount, nstep) == 0 .OR.
            IF (MOD (ncount, nstep)  == 0.OR.iptcreat == 1) THEN
               WRITE (iprint, 99002, ERR = 100) irec, &
                      MOD(ncount,nstep) + 1
            ELSE
               WRITE (iprint, 99202, ERR = 100) irec, &
                      MOD (ncount,nstep) + 1
            ENDIF
99002       FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I4,  &
                 '  W R I T T E N   O N   D I S K <------', /)
99202       FORMAT (//, ' ------> C Y C L E   N O  : ', I4,'/', I4, &
                 ' <------', /)
!
!--Write present time
!
            WRITE (iprint, 99003, ERR = 100) gt
99003       FORMAT (1X, 'TIME  : ', 1PE17.10)
            tcomp = SQRT ( (3 * pi) / (32 * rhozerox) )
!!            tcomp = SQRT((3 * pi) / (32 * rhozero))
            IF (nptmass == 2) THEN
               WRITE (iprint, 99204, ERR = 100) gt / (2.0 * pi)
99204          FORMAT (1X, '( in units of orbital period : ', &
                    1PE17.10,' )')
            ELSE
               WRITE (iprint, 99203, ERR = 100) gt / tcomp
99203          FORMAT (1X, '( in free-fall time unit : ', 1PE17.10,' )')
            ENDIF
!
!--Energies + total angular momentum
!
            total = tkin + tgrav + tterm
            WRITE (iprint, 99004, ERR = 100) total, tkin, trotz, trotx, &
                   tgrav, tterm, angto
99004       FORMAT (' General properties of system : ', /,                 &
                 ' total energy                       : ', 1PE14.5, /,  &
                 ' kinetic energy                     : ', 1PE14.5, /,  &
                 ' rotational energy around z         : ', 1PE14.5, /,  &
                 ' rotational energy around x         : ', 1PE14.5, /,  &
                 ' potential energy                   : ', 1PE14.5, /,  &
                 ' internal energy                    : ', 1PE14.5, /,  &
                 ' total angular momentum             : ', 1PE14.5)
            IF (tgrav == 0.0) THEN
               WRITE (iprint, 99020) 1, pmass (1), 2, pmass (2)
99020          FORMAT ('pmass(',I3,')=',1PE12.5, &
                   ' pmass(',I3,')=',1pE12.5)
            ENDIF
            alph = tterm/ABS(tgrav)
            betl = trotz/ABS(tgrav)
            betr = trotx/ABS(tgrav)
            ajeans = 1/alph
            WRITE (iprint, 88001, ERR = 100) alph, betl, betr, ajeans
88001       FORMAT (' evolutionary energy parmeters : ', /, &
                 '                        alpha  : ', 1PE14.5, /, &
                 '       (z)       beta parallel : ', 1PE14.5, /, &
                 '       (x)       beta perpend  : ', 1PE14.5, /, &
                 '                 Jeans number  : ', 1PE14.5)
!
!--Escapors
!
            WRITE (iprint, 99005, ERR = 100) escap
99005       FORMAT (' escapors (total mass)         : ', 1PE14.5)
!
!--Accretion
!
            WRITE (iprint, 88002, ERR = 100) nactive
88002       FORMAT (' number of active particles    : ', I8)
            DO i = 1, nptmass
               IF (i == 1) THEN
                  WRITE (iprint, 88009, ERR = 100)
88009             FORMAT (' number of accreted part.', &
                       ' onto point masses :')
               ENDIF
               WRITE (iprint, 88010, ERR = 100) i, nactotx (i)
88010          FORMAT (24X, I6,' : ', I9)
            ENDDO
            IF (nptmass == 2) THEN
               IF (isplit /= 0) THEN
                  WRITE (iprint, 88013, ERR=100) numsplit
88013             FORMAT (' number of split particles',15x, &
                       ' : ', I8)
               ENDIF
               WRITE (iprint, 88014, ERR = 100) 3-isphcom, ncapt
88014          FORMAT (' number of part. captured', &
                       ' by point mass',I2,' : ', I8)
            ENDIF
            IF (ibound == 8 .OR. ibound >= 90) THEN
               WRITE (iprint, 88003, ERR = 100) nreassign
               IF (ibound == 94 .OR. ibound == 96) THEN
                  DO i = 1, nptmass
                     WRITE (iprint, 88010, ERR = 100) i, NINT(sinj(i))
                  ENDDO
               ELSEIF (ibound == 95 .OR. ibound == 97) THEN
!!                  DO i = 1, nptmass + 1
                  DO i = 1, 3
                     WRITE (iprint, 88010, ERR = 100) i, NINT(sinj(i))
                  ENDDO
               ENDIF
88003          FORMAT (' number of reassigned part.    : ', I9)
               WRITE (iprint, 88004, ERR = 100) naccrete
88004          FORMAT (' number of accreted part.      : ', I9)
               WRITE (iprint, 88005, ERR = 100) nkill
88005          FORMAT (' number of killed particles    : ', I9)
               WRITE (iprint, 88011, ERR = 100) nghost
88011          FORMAT (' number of ghost particles     : ', I9)
               WRITE (iprint, 88012, ERR = 100) inshell
88012          FORMAT (' number of part. within rshell : ', I9)
            ENDIF
            IF (iaccevol == 'v'.OR.iaccevol == 's') THEN
               WRITE (iprint, 88006, ERR = 100)
88006          FORMAT (' point mass data:')
               DO i = 1, nptmass
                  ipt = listpm (i)
                  WRITE (iprint, 88007, ERR = 100) i, iphase (ipt), &
                         h (ipt)
88007             FORMAT ('  point mass: ',I2,' type: ',I1, &
                       ' hacc: ',1PE12.3)
               ENDDO
               WRITE (iprint, 88008) rmax
88008          FORMAT (' new boundary radius: ',1PE12.5)
            ENDIF
!
!--Object no 1
!
            WRITE (iprint, 99006, ERR = 100) n1, cmx1, cmy1, cmz1, &
                vcmx1, vcmy1, vcmz1, hmi1, hma1, dmax1, zmax1, romean1, &
                romax1, rocen1
99006       FORMAT (/, ' Object number 1 (', I8, ' particles ) : ', /, &
                 ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,     &
                 '  z :', 1PE14.5, /, ' velocity cm    vx :', 1PE14.5,  &
                 ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,                 &
                 ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5, /,  &
                 ' max. dist. cm   r :', 1PE14.5, '  z :', 1PE14.5, /,  &
                 ' density      mean :', 1PE14.5, ' max:', 1PE14.5,     &
                 ' cen:', 1PE14.5, /)
!
!--Object no 2
!
            IF (n2 /= 0) THEN
               WRITE (iprint, 99007, ERR = 100) n2, cmx2, cmy2, cmz2, &
                  vcmx2, vcmy2, vcmz2, hmi2, hma2, dmax2, zmax2, romean2, &
                  romax2, rocen2
99007          FORMAT (/, ' Object number 2 (', I8, ' particles ) : ', /, &
                    ' center of mass  x :', 1PE14.5, '  y :', 1PE14.5,  &
                    '  z :', 1PE14.5, /, ' velocity cm    vx :',        &
                    1PE14.5, ' vy :', 1PE14.5, ' vz :', 1PE14.5, /,     &
                    ' smoothing l.  min :', 1PE14.5, ' max:', 1PE14.5,  &
                    /, ' max. dist. cm   r :', 1PE14.5, '  z :',        &
                    1PE14.5, /, ' density      mean :', 1PE14.5,        &
                    ' max:', 1PE14.5,' cen:', 1PE14.5, /)
            ENDIF
!
!--Write transfer output
!
         ELSEIF (where (1:5)  == 'trans') THEN
            WRITE (iprint, 99008) ibegin, iend, file1, file2
99008       FORMAT (//, ' dumps no ', I4, ' to ', I4, &
                 ' copied from file ', A7, ' to file ', A7, /)
            IF (ichang == 0) WRITE (iprint, 99009)
99009       FORMAT (' no change done. ')
            IF (ichang == 1) WRITE (iprint, 99010) frac, energc
99010       FORMAT (' the inner ', F6.3, &
                 ' of the total number of particles', &
                 ' have had their internal energy multiplied by: ', &
                 1PE12.5)
            IF (ichang == 2) THEN
               WRITE (iprint, 99010) frac, energc
               WRITE (iprint, 99011) vexpan * udist / utime, rnorm *    &
                  udist
99011          FORMAT (' a homologous expansion of ', 1PE12.5, ' cm/s', &
                    ' at ', 1PE12.5, ' cm has been added')
            ENDIF
!
!--Write output for initialisation
!
         ELSEIF (where (1:6)  == 'newrun') THEN
            IF (varsta == 'entropy') sentenc = 'specific entropy'
            IF (varsta == 'intener') sentenc = 'specific internal energy'

            IF (what == 'scratch'.OR.what == 's') THEN
               WRITE (iprint, 99012)
99012          FORMAT (' New initial conditions defined from scratch')
               WRITE (iprint, 99013) idist
99013          FORMAT (' Particles distributed according to ', &
                 'distribution no :', I2)
               IF (ibound >= 94 .AND. ibound <= 97) THEN
                  WRITE (iprint, "(' Values of ', A24, ' : ', 1P2E12.5)") &
                     sentenc, thermal1, thermal2
               ELSE
                  WRITE (iprint, "(' Values of ', A24, ' : ', 1PE12.5)") &
                     sentenc, thermal1
               END IF

               IF (encal == 'i') WRITE (iprint, 99100)
99100          FORMAT (' Isothermal equation of state')
               IF (encal == 'a') WRITE (iprint, 99103)
99103          FORMAT (' Adiabatic equation of state')
               IF (encal == 'p') WRITE (iprint, 99101) gamma
99101          FORMAT (' Polytropic equation of state gamma is ',1PE12.5)
               IF (encal == 'v') WRITE (iprint, 99102)
99102          FORMAT (' Gamma variable equation of state')
               IF (encal == 'x') WRITE (iprint, 99104)
99104          FORMAT (' Physical equation of state')
               IF (encal == 'c') WRITE (iprint, 99105)
99105          FORMAT (' Gas with shock heating and radiative cooling')

               IF (idist >= 1) WRITE (iprint, 99015) xlmax, xlmin, &
                     ylmax, ylmin, zlmax, zlmin
99015          FORMAT (/,' Particles are in volume given by :', /, &
                    ' xmax : ', 1PE12.5, ' xmin : ', 1PE12.5, /, &
                    ' ymax : ', 1PE12.5, ' ymin : ', 1PE12.5, /, &
                    ' zmax : ', 1PE12.5, ' zmin : ', 1PE12.5)
            ENDIF

            IF (what (1:5)  == 'exist'.OR.what == 'e') THEN
               WRITE (iprint, 99016) file1, file2
99016          FORMAT (//, ' new initial conditions made from file ', A7, &
                    ' and ', A7)
               WRITE (iprint, 99017) n1, xx1, yy1, zz1, vvx1, vvy1,     &
                  vvz1
99017          FORMAT (/, ' first   object is made of ', I4, &
                    ' particles : ', /, '  x : ', 1PE12.5, '   y : ',   &
                    1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,   &
                    '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
               WRITE (iprint, 99018) n2, xx2, yy2, zz2, vvx2, vvy2,     &
                  vvz2
99018          FORMAT (/, ' second  object is made of ', I4, &
                    ' particles : ', /, '  x : ', 1PE12.5, '   y : ',   &
                    1PE12.5, '   z :', 1PE12.5, /, ' vx : ', 1PE12.5,   &
                    '  vy : ', 1PE12.5, '  vz :', 1PE12.5)
               WRITE (iprint, 99019) angto
99019          FORMAT (/, ' total angular momentum of the system : ', &
                    1PE12.5)
            ENDIF
         ENDIF

      ENDIF

      GOTO 200
!
!--Handle errors during writing
!
  100 where = 'prout'
      CALL error (where, 1)

  200 CONTINUE

      END SUBROUTINE prout
