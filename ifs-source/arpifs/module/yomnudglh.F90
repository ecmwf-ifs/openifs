! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMNUDGLH
! define namelist switches for latent heat nudging and related values 
! Author: Florian Meier *ZAMG* based on other such modules 
! Original nov 2015

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE
! define namelist switches for latent heat nudging and related values 

! Logical to activate LH nudging
LOGICAL :: LNUDGLH
! LOGICAL to activate moisture compensation for LHN increment (=keep RH constant)
LOGICAL :: LNUDGLHCOMPT
! LOGICAL to activate replacement of 0-LH profile by horizontal mean profile
LOGICAL :: LNUDGLHCLIM
! LOGICAL to activate usage of environment profile in case of rain obs but not rain model=0
LOGICAL :: LNUDGLHENVI
! in case of model but no obs rain nudging against critical humidity profile
! not saturated
LOGICAl :: LNUDGLHSUBSAT
! logical if true use Stephan 2008 LHN factor: ln(RRobs/RRmodel) instead of (RRobs-RRmodel)/RRmodel
LOGICAL :: LNUDGLHSTEPH
!
LOGICAL :: LNUDGLHREAD
! Integer
INTEGER(KIND=JPIM) :: NSTARTNUDGLH  ! fist time step of the lh nudging
INTEGER(KIND=JPIM) :: NSTOPNUDGLH   ! last time step of the lh nudging
INTEGER(KIND=JPIM) :: NINTNUDGLH    ! interval of lh nudging steps
INTEGER(KIND=JPIM) :: NTAUNUDGLH    ! number of timesteps after 
                                    !  LHN-observation step, where LHN 
                                    !  is still applied with old factor
INTEGER(KIND=JPIM) :: NTOPNUDGLH    ! top level where LHN is not applied above
INTEGER(KIND=JPIM) :: NBOTTOMNUDGLH ! bottom level where LHN is not applied below
! REAL
REAL(KIND=JPRB) :: TSTARTNUDGLH ! fist time (sec.) of the lh nudging
REAL(KIND=JPRB) :: TSTOPNUDGLH  ! last time (sec.) of the lh nudging
REAL(KIND=JPRB) :: TINTNUDGLH   ! time interval of nudging steps (sec.)
REAL(KIND=JPRB) :: TTAUNUDGLH   ! seconds after LHN-observation step, where LHN
                                ! is still applied with old factor
REAL(KIND=JPRB) :: RDAMPNUDGLH  ! damping of LHN factor for following LHN steps in case of NTAUNUDGLH>1
REAL(KIND=JPRB) :: RPRTOPSNUDGLH ! pressure Pa of clim profile top (stratiform)
REAL(KIND=JPRB) :: RPRTOPANUDGLH !
REAL(KIND=JPRB) :: RPRTOPBNUDGLH !
REAL(KIND=JPRB) :: RPRTOPCNUDGLH !
REAL(KIND=JPRB) :: RPRTOPDNUDGLH ! pressure Pa of clim profile top (convective)
REAL(KIND=JPRB) :: RPRTOPTNUDGLH !
REAL(KIND=JPRB) :: RPRTOPWNUDGLH !
REAL(KIND=JPRB) :: RPRTOPSTABNUDGLH ! pressure Pa of clim profile top (stable)
REAL(KIND=JPRB) :: RPRMAXSNUDGLH ! pressure Pa of clim profile LH max (stratiform)
REAL(KIND=JPRB) :: RPRMAXANUDGLH !
REAL(KIND=JPRB) :: RPRMAXBNUDGLH !
REAL(KIND=JPRB) :: RPRMAXCNUDGLH !
REAL(KIND=JPRB) :: RPRMAXDNUDGLH ! pressure Pa of clim profile LH max (convective)
REAL(KIND=JPRB) :: RPRMAXTNUDGLH !
REAL(KIND=JPRB) :: RPRMAXWNUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMSNUDGLH ! pressure Pa of clim profile bottom (stratiform)
REAL(KIND=JPRB) :: RPRBOTTOMANUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMBNUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMCNUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMDNUDGLH ! pressure Pa of clim profile bottom (convective)
REAL(KIND=JPRB) :: RPRBOTTOMTNUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMWNUDGLH !
REAL(KIND=JPRB) :: RPRBOTTOMSTABNUDGLH ! pressure Pa of clim profile bottom (stable)
REAL(KIND=JPRB) :: RHEATSNUDGLH ! max value of LH K in clim profile (stratiform)
REAL(KIND=JPRB) :: RHEATCNUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATANUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATBNUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATDNUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATTNUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATWNUDGLH ! max value of LH K in clim profile (convective)
REAL(KIND=JPRB) :: RHEATSTABNUDGLH ! max value of LH K in clim profile (stable)
REAL(KIND=JPRB) :: RAMPLIFY ! factor multiflied with LHN increment default 1.0
REAL(KIND=JPRB) :: RMAXNUDGLH ! upper limit for LHN-tendency for security reason default 3.0
REAL(KIND=JPRB) :: RMINNUDGLH ! lower limit for LHN-tendency for security reason default -0.001
REAL(KIND=JPRB) :: RORONUDGLH ! reduce LHN over mountains
REAL(KIND=JPRB) :: RALPHANUDGLH ! constant from Macpherson 1997
REAL(KIND=JPRB) :: REPSILONNUDGLH ! constant from Macpherson 1997
REAL(KIND=JPRB) :: RFACNUDGLH ! nudging factor for subsaturated profile LNUDGLHSUBSAT=T
INTEGER(KIND=JPIM) :: NTIMESPLITNUDGLH ! segregate nudging interval for LHN from observation interval
                                    ! by using the same observation divided by a factor several times
                                    ! to reduce time lag between observation and effect on the model
                                
END MODULE YOMNUDGLH
