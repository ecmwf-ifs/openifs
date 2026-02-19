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

MODULE TESTVAR_MIX

!     Purpose.
!     --------
!       Setup for tests of variational assimilation components.

!     Author.
!     -------
!       Y. Tremolet
!       Original    20-11-03

!     Modifications.
!     --------------
!        Y.Tremolet    18-Feb-2005 Identical twin experiments
!        Y.Tremolet    30-Mar-2005 Improve incremental 4D-Var test
!        A.Geer        27-Jun-2013 Per-observation adjoint test for radiances
!        K. Yessad (July 2014): Move some variables.
!        A.Geer        20-Aug-2015 Switches in support of offline HOP_DRIVER
!        A. Geer       10-Sep-2015 Use checksums for bit-reproducibility tracing
!        Y. Michel, MF,  Feb 2020 Jb test (former CV test)
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB
USE YOMHOOK , ONLY: LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY: NULOUT, NULERR, NULNAM
USE IOSTREAM_MIX, ONLY: SETUP_IOSTREAM,SETUP_IOREQUEST,IO_INQUIRE, &
 & CLOSE_IOSTREAM, TYPE_IOSTREAM, TYPE_IOREQUEST, CLOSE_IOREQUEST

IMPLICIT NONE
SAVE

PRIVATE
PUBLIC SETUP_TESTVAR, &
 & NTESTVAR,LINTEST,LTLTEST,LADTEST,LGRTEST,LJBTEST, &
 & NSTOP_TLTEST, NPHYS_TLTEST, LTESTADJM, LTESTADJH, LTESTANY, &
 & LTESTINC, LSTEP_TESTINC, INTERP_TEST, NGRB_TESTINC, LTESTSPPA, &
 & LTWINTRUTH, LTWINEXP, TWINDEPART, LRTTOV_ADTEST, &
 & LGOM_PLUS_DUMP, LGOM_PLUS_NETCDF, LHOP_RESULTS, LGOM_CHECKSUM

LOGICAL :: LTLTEST   ! Test of tangent linear model
LOGICAL :: LINTEST   ! Test linearity of tangent linear model
LOGICAL :: LADTEST   ! Test of adjoint
LOGICAL :: LGRTEST   ! Gradient test
LOGICAL :: LJBTEST   ! Test of adjoint for Jb and control variable

LOGICAL :: LRTTOV_ADTEST ! Per-observation adjoint test for RTTOV_EC inside RADTR_ML_TL
                         ! This can be run independently or during 4D-Var adjoint test (LADTEST)

INTEGER(KIND=JPIM) :: NTESTVAR       ! Controls when tests are performed:
                            ! 0 - No tests
                            ! 1 - Tests run before minimisation
                            ! 2 - Tests run after  minimisation
                            ! 3 - Tests run before and after minimisation

INTEGER(KIND=JPIM) :: NSTOP_TLTEST   ! Number of time-steps to run TL test for
INTEGER(KIND=JPIM) :: NPHYS_TLTEST   ! Type of physics used in TL test:
                            ! 0 - No physics in NL or TL models
                            ! 1 - Only physics available in TL model is used
                            !     NL model
                            ! 2 - All physics available in each model is used
LOGICAL :: LTESTADJM        ! Test of M adjoint is being run (internal switch)
LOGICAL :: LTESTADJH        ! Test of H adjoint is being run (internal switch)
LOGICAL :: LTESTANY         ! Any of the tests is being run
LOGICAL :: LTESTINC                  ! Save data for test of incremental 4D-Var
REAL(KIND=JPRB) :: TSTEP_TESTINC     ! Interval (seconds) for incr. 4D-Var test
                                     ! save for initial time only if set to 0
REAL(KIND=JPRB) :: TSTEP_TESTBGN     ! Interval of additional incr. 4D-Var test
REAL(KIND=JPRB) :: TSTOP_TESTBGN     ! Stop time   additional incr. 4D-Var test
INTEGER(KIND=JPIM) :: INTERP_TEST    ! Interpolation algo for incr. 4D-Var test
LOGICAL, ALLOCATABLE :: LSTEP_TESTINC(:)   ! Active steps for incr. 4D-Var test
INTEGER(KIND=JPIM) :: NGRB_TESTINC
LOGICAL :: LTESTSPPA                 ! Surf pressure in Pa (LNSP otherwise)

LOGICAL :: LTWINTRUTH          ! Running "truth" for identical twin experiment
LOGICAL :: LTWINEXP            ! Running identical twin experiment
REAL(KIND=JPRB) :: TWINDEPART  ! Obs departure from truth in twin exp

! For use with the HOP_DRIVER offline test harness for the observation operator
LOGICAL :: LGOM_PLUS_DUMP  ! Dump GOM plus to file (model fields at obs locations)
LOGICAL :: LGOM_PLUS_NETCDF! Deafult format of this file is F77 binary; make it NetCDF (incompatible with offline test harness)
LOGICAL :: LHOP_RESULTS    ! Write text files containing the HOP output     

! Compute & print a checksum for all GOM variables, useful for trakcing bit-reproducility issues:
LOGICAL :: LGOM_CHECKSUM       

! ------------------------------------------------------------------
NAMELIST/NAMTESTVAR/NTESTVAR,LINTEST,LTLTEST,LADTEST,LGRTEST,LJBTEST, &
                  & NSTOP_TLTEST,NPHYS_TLTEST, LTESTINC, &
                  & INTERP_TEST,TSTEP_TESTINC,TSTEP_TESTBGN,TSTOP_TESTBGN, &
                  & LTWINTRUTH, LTWINEXP, TWINDEPART, LRTTOV_ADTEST, &
                  & LGOM_PLUS_DUMP, LGOM_PLUS_NETCDF, LHOP_RESULTS, LGOM_CHECKSUM
! ------------------------------------------------------------------

#include "posnam.intfb.h"

! ------------------------------------------------------------------
CONTAINS
! ------------------------------------------------------------------

SUBROUTINE SETUP_TESTVAR(YDRIP)

USE YOMRIP , ONLY : TRIP
TYPE(TRIP),INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM) :: ISTEP_TESTINC, ISTEP_TESTBGN, ISTOP_TESTBGN, II, JJ
CHARACTER(LEN=12) :: CLFILE
TYPE(TYPE_IOSTREAM) :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('TESTVAR_MIX:SETUP_TESTVAR',0,ZHOOK_HANDLE)
ASSOCIATE(NSTART=>YDRIP%NSTART, NSTOP=>YDRIP%NSTOP, TSTEP=>YDRIP%TSTEP)
! Internal consistency tests (inner loop)
NTESTVAR=0
LTLTEST =.FALSE.
LINTEST =.FALSE.
LADTEST =.FALSE.
LGRTEST =.FALSE.
LJBTEST =.FALSE.
LRTTOV_ADTEST =.FALSE.
NSTOP_TLTEST=-1
NPHYS_TLTEST=0

! Incremental 4D-Var tests (inner vs. outer loops)
LTESTINC=.FALSE.
TSTEP_TESTINC=0.0_JPRB
TSTEP_TESTBGN=0.0_JPRB
TSTOP_TESTBGN=0.0_JPRB
INTERP_TEST=1          ! 3 for conserving interpolation
LTESTSPPA=.TRUE.

! Identical twin experiments
LTWINTRUTH=.FALSE.
LTWINEXP=.FALSE.
TWINDEPART=0.0_JPRB

! Offline hop_driver support
LGOM_PLUS_DUMP=.FALSE.
LGOM_PLUS_NETCDF=.FALSE.
LHOP_RESULTS=.FALSE.

LGOM_CHECKSUM=.FALSE.

CALL POSNAM(NULNAM,'NAMTESTVAR')
READ(NULNAM,NAMTESTVAR)

WRITE(NULOUT,*)'TESTVAR_MIX: NTESTVAR=',NTESTVAR
WRITE(NULOUT,*)'TESTVAR_MIX: LINTEST, LTLTEST, LADTEST, LGRTEST, LJBTEST= ',&
                           & LINTEST, LTLTEST, LADTEST, LGRTEST, LJBTEST  
WRITE(NULOUT,*)'TESTVAR_MIX: LRTTOV_ADTEST= ', LRTTOV_ADTEST 
WRITE(NULOUT,*)'TESTVAR_MIX: LGOM_PLUS_DUMP, LGOM_PLUS_NETCDF, LHOP_RESULTS= ', LGOM_PLUS_DUMP, LGOM_PLUS_NETCDF, LHOP_RESULTS
WRITE(NULOUT,*)'TESTVAR_MIX: LGOM_CHECKSUM= ', LGOM_CHECKSUM
IF (LTLTEST) THEN
  WRITE(NULOUT,*)'TESTVAR_MIX: NSTOP_TLTEST=',NSTOP_TLTEST,&
                           & ' NPHYS_TLTEST=',NPHYS_TLTEST
ENDIF

WRITE(NULOUT,*)'TESTVAR_MIX: LTESTINC=',LTESTINC
IF (LTESTINC) THEN
  ALLOCATE(LSTEP_TESTINC(NSTART:NSTOP))
  LSTEP_TESTINC(:)=.FALSE.
  IF (TSTEP_TESTINC==0.0_JPRB) LSTEP_TESTINC(NSTART)=.TRUE.
  ISTEP_TESTINC=NINT(TSTEP_TESTINC/TSTEP)
  IF (ISTEP_TESTINC>0) THEN
    DO JJ=NSTART,NSTOP
      II=JJ-NSTART
      LSTEP_TESTINC(JJ)=(MOD(II,ISTEP_TESTINC)==0)
    ENDDO
  ENDIF
  ISTEP_TESTBGN=NINT(TSTEP_TESTBGN/TSTEP)
  ISTOP_TESTBGN=NINT(TSTOP_TESTBGN/TSTEP)
  IF (ISTEP_TESTBGN>0 .AND. ISTOP_TESTBGN>0) THEN
    DO JJ=NSTART,NSTART+ISTOP_TESTBGN
      II=JJ-NSTART
      LSTEP_TESTINC(JJ)=LSTEP_TESTINC(JJ).OR.(MOD(II,ISTEP_TESTBGN)==0)
    ENDDO
  ENDIF
  WRITE(NULOUT,*)'TESTVAR_MIX: TSTEP_TESTINC=',TSTEP_TESTINC
  WRITE(NULOUT,*)'TESTVAR_MIX: TSTEP_TESTBGN, TSTOP_TESTBGN=',&
                         & TSTEP_TESTBGN, TSTOP_TESTBGN
  WRITE(NULOUT,*)'TESTVAR_MIX: LSTEP_TESTINC=',LSTEP_TESTINC(:)
  CLFILE='test4dinc_gg'
  CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',CLFILE,'r',KIOMASTER=1)
  CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,&
   & LDINTONLY=.TRUE.,KCHUNKSIZE=1)
  CALL IO_INQUIRE(YL_IOSTREAM,YL_IOREQUEST,KGRIB_HANDLE=NGRB_TESTINC)
  CALL CLOSE_IOREQUEST(YL_IOREQUEST)
  CALL CLOSE_IOSTREAM(YL_IOSTREAM)
ENDIF

WRITE(NULOUT,*)'TESTVAR_MIX: LTWINTRUTH= ',LTWINTRUTH,', LTWINEXP= ',LTWINEXP
WRITE(NULOUT,*)'TESTVAR_MIX: TWINDEPART=',TWINDEPART
IF (TWINDEPART/=0.0_JPRB) THEN
  WRITE(NULOUT,*)'TESTVAR_MIX: WARNING: non zero departure specified'
  WRITE(NULERR,*)'TESTVAR_MIX: WARNING: non zero departure specified'
ENDIF

LTESTADJM=.FALSE.
LTESTADJH=.FALSE.
LTESTANY=.FALSE.

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TESTVAR_MIX:SETUP_TESTVAR',1,ZHOOK_HANDLE)
END SUBROUTINE SETUP_TESTVAR

! ------------------------------------------------------------------

END MODULE TESTVAR_MIX
