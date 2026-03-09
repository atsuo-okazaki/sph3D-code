      SUBROUTINE wrinsph 
!************************************************************           
!                                                           *           
!  This subroutine writes a new file inname from parameters *           
!     taken from setpart                                    *           
!                                                           *           
!************************************************************           
                                                                        
!!      use mpi_mod 
      use idims

      use constants
      use typef
      use units
      use dissi
      use rotat
      use tming
      use integ
      use varet
      use recor
      use rbnd
      use diskbd
      use expan
      use kerne
      use files
      use actio
      use logun
      use debug
      use cgas
      use stepopt
      use init
      use pres
      use xforce
      use numpa
      use ptmass
      use binary
      use ptdump
      use crpart
      use ptbin
      use useles
      use new
      use xtorq
      use sphcom
      use misali
      use maslos
      use winds
      use split
      use polyk2

      implicit none

!      REAL(DP) :: omeg0r, rminus1, rzero
      REAL(DP) :: omeg0r, rminus1
!                                                                       
!--Allow for tracing flow                                               
!                                                                       
!!      IF (myrank.eq.0) then 
         IF (itrace == 'all') WRITE (iprint, 99001) 
!!      ENDIF 
99001 FORMAT (' entry subroutine wrinsph') 
!                                                                       
!--Open input file                                                      
!                                                                       
      OPEN (iterm, FILE = inname) 
!                                                                       
!--Determine options for evolution run                                  
!                                                                       
!--Write name of run                                                    
!                                                                       
!!      IF (myrank.eq.0) then 
                                                                        
         WRITE (iterm, 99002) namenextrun 
99002 FORMAT (A20) 
!                                                                       
!--Write name of file containing physical input                         
!                                                                       
         WRITE (iterm, 99003) file1 
99003 FORMAT (A7) 
         WRITE (iterm, 99003) varsta 
!                                                                       
!--Write options                                                        
!                                                                       
         WRITE (iterm, 99006) encal 
99006 FORMAT (A1) 
                                                                        
         WRITE (iterm, 89000) initialptm 
89000 FORMAT(2X,I2,'  Point Masses Initially') 
         WRITE (iterm, 89005) iaccevol 
89005 FORMAT(A1,'  Variable:Roche,Sep/Fixed Accretion Radii (v/s/f)') 
         IF (iaccevol == 'v'.OR.iaccevol == 's') WRITE (iterm, 89010)   &
         accfac                                                         
89010 FORMAT(2X,1PE12.5,'  Fraction of Roche lobe size') 
         WRITE (iterm, 89015) iptmass 
89015 FORMAT(2X,I2,'  Point Mass Creation') 
         WRITE (iterm, 88000) igrp 
88000 FORMAT(2X,I2,'  Pressure') 
         WRITE (iterm, 88001) igphi 
88001 FORMAT(2X,I2,'  Gravity') 
         WRITE (iterm, 88002) ifsvi, alpha, beta 
88002 FORMAT(2X,I2,1X,1F6.4,1X,1F6.4,'  Artificial Viscosity') 
         WRITE (iterm, 88003) ifcor 
88003 FORMAT(2X,I2,'  Coriolis Forces') 
         WRITE (iterm, 88004) ichoc 
88004 FORMAT(2X,I2,'  Heating from Shocks') 
         WRITE (iterm, 88005) iener 
88005 FORMAT(2X,I2,'  PdV') 
         WRITE (iterm, 88006) damp 
88006 FORMAT(2X,1PE12.5,'  Damping') 
         WRITE (iterm, 88007) ibound 
88007 FORMAT(2X,I2,'  Boundry') 
         WRITE (iterm, 88008) iexf 
88008 FORMAT(2X,I2,'  External Forces') 
         WRITE (iterm, 88009) iexpan 
88009 FORMAT(2X,I2,'  Expansion') 
         WRITE (iterm, 88010) nstep 
88010 FORMAT(I4,'  Binary Dump Every N Max Time Steps') 
         WRITE (iterm, 89020) iptoutnum 
89020 FORMAT(I4,'  N Ptmass Dumps Per Max Timestep') 
         WRITE (iterm, 88017) umass, udist, utime 
88017 FORMAT(2X,3(1PD22.15,1X),'Units (umass, udist, utime)') 
         WRITE (iterm, 88011) tol, tolptm, tolh 
88011 FORMAT (2X,3(1PE12.5,1X),'Tolerance (Gas, Ptmass, Smooth Len)') 
         WRITE (iterm, 88012) ipos 
88012 FORMAT (I4,'  File Position') 
         WRITE (iterm, 88013) tmax 
88013 FORMAT (2X,1PE12.5,'  Max Time (min) (Status Written)') 
         WRITE (iterm, 88014) tstop 
88014 FORMAT (2X,1PE12.5,'  Max Dynamic Time') 
         WRITE (iterm, 88015) dtmax 
88015 FORMAT (2X,1PE14.7,'  Max Timestep') 
         WRITE (iterm, 88016) dtini 
88016 FORMAT (2X,1PE12.5,'  Initial Timestep') 
                                                                        
!!      ENDIF 
                                                                        
      IF (ifcor.NE.0) THEN 
         IF (job (1:9)  == 'evolution') THEN 
            omeg0r = omeg0 / utime 
         ELSE 
            omeg0r = omeg0 
         ENDIF 
!!         IF (myrank.eq.0) THEN
            WRITE (iterm, 88020) omeg0r 
!!         END IF
88020 FORMAT   (2X,1PE12.5,'  Omega') 
      ENDIF 
                                                                        
!!      IF (myrank.eq.0) then 
                                                                        
         IF (iexpan.NE.0) THEN 
            WRITE (iterm, 88025) vexpan 
88025 FORMAT   (2X,1PE12.5,'  Expansion Velocity') 
         ENDIF 
                                                                        
         IF (ibound == 7) THEN 
            WRITE (iterm, 88035) hmaximum 
88035 FORMAT   (2X,1PE12.5,'  Maximum h (0 to disable)') 
            WRITE (iterm, 88030) pext 
88030 FORMAT   (2X,1PE12.5,'  Const. External Pressure') 
         ENDIF 
                                                                        
         IF (ibound == 8) THEN 
            WRITE (iterm, 88036) deadbound 
88036 FORMAT   (2X,1PE12.5,'  Dead Particle Boundary') 
            WRITE (iterm, 88037) fractan, fracradial, nstop, nfastd 
88037 FORMAT   (2X,2(1PE12.5,1X),I6,1X,I6,                              &
              '  New Particle Velocities And Nstop')                    
         ENDIF 
         IF (ibound.GE.90) THEN 
            WRITE (iterm, 88039) hmaximum 
88039 FORMAT   (2X,1PE12.5,'  Maximum smoothing length') 
            WRITE (iterm, 88036) deadbound 
            WRITE (iterm, 88038) fractan, fracradial, nshell, rshell 
88038 FORMAT   (2X,2(1PE12.5,1X),I6,1X,1PE12.5,                         &
              '  New Part. Vels., Nshell, and Rshell')                  
         ENDIF 
                                                                        
         IF (iexf == 5.OR.iexf == 6) THEN 
            WRITE (iterm, 88040) xmass 
88040 FORMAT   (2X,1PE12.5,'  External Forces Mass') 
         ENDIF 
                                                                        
         IF (iptmass.NE.0.OR.nptmass.NE.0) THEN 
            WRITE (iterm, 88042) hacc 
88042 FORMAT   (2X,1PE12.5,'  Outer Accretion Radius') 
            WRITE (iterm, 88044) haccall 
88044 FORMAT   (2X,1PE12.5,'  Inner Accretion Radius') 
         ENDIF 
                                                                        
         IF (iptmass.GE.1) THEN 
            WRITE (iterm, 88046) radcrit 
88046 FORMAT   (2X,1PE12.5,'  Critical Radius') 
            WRITE (iterm, 88048) ptmcrit 
88048 FORMAT   (2X,1PE12.5,'  Critical Density') 
                                                                        
         ENDIF 
                                                                        
!!      ENDIF 
                                                                        
!      rzero = 0.0 
      rminus1 = - 1.0 
                                                                        
!!      IF (myrank.eq.0) then 
                                                                        
         IF (ibound == 1.OR.ibound == 3.OR.ibound == 8.OR.ibound.GE.90) &
         THEN                                                           
            WRITE (iterm, 88050) igeom, rmind, rmax, xmin, xmax, ymin,  &
            ymax, zmin, zmax                                            
         ELSEIF (ibound == 2) THEN 
            WRITE (iterm, 88050) igeom, rmind, rcyl, xmin, xmax, ymin,  &
            ymax, zmin, zmax                                            
         ELSE 
            IF (rmax == rcyl) THEN 
               WRITE (iterm, 88050) igeom, rminus1, rmax, xmin, xmax,   &
               ymin, ymax, zmin, zmax                                   
            ELSE 
               WRITE (iterm, 88050) igeom, rmind, rcyl, xmin, xmax,     &
               ymin, ymax, zmin, zmax                                   
            ENDIF 
         ENDIF 
88050 FORMAT(1X,I3,1P8E13.5) 
                                                                        
         WRITE (iterm, 88060) rptmas(1), rptmas(2) 
88060 FORMAT(2X,2(1PE12.5,1X),                                          &
              '  Radii of point masses 1 and 2')                        
         WRITE (iterm, 88062) isphcom 
88062 FORMAT(2X,I2,'  Centre of mass of SPH particles') 
         WRITE (iterm, 88064) rangle(1), rangle(2) 
88064 FORMAT(2X,2(1PE12.5,1X),                                          &
           ' rangle(1), rangle(2): Two Euler angles for misalignment')  
         IF (ibound == 93) THEN 
            WRITE (iterm, 88066) emdot0, partm 
88066 FORMAT   (2X,1P2E15.8,                                            &
                '  mass injection rate in units of nshell/utime')       
         ELSEIF (ibound == 94 .OR. ibound == 96) THEN 
            WRITE (iterm, 88065) nshell1(1), nshell1(2), rshell1,     &
            rshell2, vwind1, vwind2, vrot1, vrot2, emdot0, emdotratio,  &
            partmass(1), partmass(2), sinj0(1), sinj0(2), &
            RK21, RK22, RK23
88065 FORMAT   (2X,2I6,1P15E15.8,                                       &
                '  parameters for stellar winds')                       
         ELSEIF (ibound == 95 .OR. ibound == 97) THEN 
            WRITE (iterm, 88072) nshell1(1), nshell1(2), rshell1,     &
            rshell2, vwind1, vwind2, vrot1, vrot2, emdot0, emdotratio,  &
            partmass(1), partmass(2), sinj0(1), sinj0(2), nshell1(3), &
            partmass(3), sinj0(3), &
            RK21, RK22, RK23
88072 FORMAT   (2X,2I6,1P12E15.8,I6,1P5E15.8,                           &
                '  parameters for stellar winds')                       
         ELSEIF (ibound == 99) THEN 
            WRITE (iterm, 88067) emdot0, sinj0(1), partm 
88067 FORMAT   (2X,3(1PE15.8,1X),                                       &
            '  injection factor, # injectd particles (model), partm')   
!--      Write name of file containing infall particles                 
            WRITE (iterm, 99002) nameinfl 
         ENDIF 
                                                                        
         IF (igeom == 8 .OR. igeom == 9) THEN 
            WRITE (iterm, "(2X,4(1PE12.5,1X), &
      &           ' range of angles for mass injection for star 1')") &
      &           azimuth1(1), azimuth2(1), vangle1(1), vangle2(1) 
            WRITE (iterm, "(2X,4(1PE12.5,1X), &
      &           ' range of angles for mass injection for star 2')") &
      &           azimuth1(2), azimuth2(2), vangle1(2), vangle2(2) 
         ENDIF 
                                                                        
         IF (iexf == 3) THEN 
            WRITE (iterm, 88069) xeps, xbeta 
88069 FORMAT   (2X,2(1PE12.5,1X),                                       &
                '  parameters for ang. mom. injection rate')            
!!      ELSEIF (iexf == 7) THEN                                         
!!         WRITE (iterm, 88070) xantgrav(1),xantgrav(2)                 
         ELSEIF (iexf == 7.OR.iexf == 8) THEN 
            WRITE (iterm, 88070) xantgrav(1), xantgrav(2) 
88070 FORMAT   (2X,1P2E12.5,1X,                                         &
                '  antigravity parameters for two winds')               
            WRITE (iterm, 88080) akappac(1), akappar(1), akappac(2), &
            akappar(2)                                                 
88080 FORMAT   (2X,1P4E12.5,1X,                                         &
                '  opacities for two winds')                            
            WRITE (iterm, 88082) vbeta(1), vbeta(2) 
88082 FORMAT   (2X,1P2E12.5,1X,                                         &
                '  beta parameters for two winds')                      
         ENDIF 

         WRITE (iterm, 88084) isplit
88084    FORMAT(2X,I2,'  Particle splitting (0:off, 1:on)')
                                                                        
!!      ENDIF 
                                                                        
!                                                                       
!--Check for consistency                                                
!                                                                       
      CALL chekopt 
                                                                        
!!      IF (myrank.eq.0) then 
         IF (idebug (1:7) == 'wrinsph') THEN 
            WRITE (iprint, 99004) igrp, igphi, ifsvi, ifcor, ichoc,     &
            iener, ibound, damp, varsta                                 
99004 FORMAT    (1X, 7(I2,1X), 2(E12.5,1X), 1X, A7) 
            WRITE (iprint, 99005) file1, ipos, nstep 
99005 FORMAT    (1X, A7, 1X, I4, 1X, I4) 
         ENDIF 
!!      ENDIF 
                                                                        
!!      IF (myrank.eq.0) THEN
         CLOSE (iterm) 
!!      ENDIF

                                                                        
      END SUBROUTINE wrinsph                        
