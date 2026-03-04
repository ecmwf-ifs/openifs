! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WAVEMDL (CBEGDAT, PSTEP, KSTOP, KSTPW,                 &
     &              NFIELDS, NGPTOTG, NC, NR,                     &
     &              IGRIB_HANDLE, RMISS, ZRCHAR, FIELDS,          &
     &              NATMFLX,                                      &
     &              LWCUR, LWSTOKES,                              &
     &              NWVFIELDS, WVFLDG,                            &
     &              NLONW, NLATW, LDSTOP, LDWRRE,                 &
     &              LDRESTARTED, ZDELATM,                         &
     &              LDWCOUNORMS, LDNORMWAMOUT_GLOBAL,             &
     &              MASK_IN, MASK_OUT,                            &
     &              FRSTIME, NADV, PRPLRADI, PRPLRG,              &
     &              RNU_ATM, RNUM_ATM,                            &
     &              IDATE_TIME_WINDOW_END, NSTEP,                 &
     &              LDIFS_IO_SERV_ENABLED )

!****  *WAVEMDL* - SUPERVISES EXECUTION OF THE WAVE MODEL

!      LIANA ZAMBRESKY    GKSS/ECMWF    OCTOBER 1988

!      MODIFICATION.
!      -------------
!         H. GUNTHER   ECMWF  MARCH    1990
!         P. LIONELLO  ECMWF  APRIL    1990  DATA ASSIMILATION
!                                            MODULE WAMASSI ADDED.
!         J. BIDLOT    ECMWF  FEBRUARY 1996  MESSAGE PASSING.
!         J. DOYLE     ECMWF  OCTOBER  1996  ATMOSPHERIC COUPLING  .
!         J. BIDLOT    ECMWF  FEBRUARY 1997  MESSAGE PASSING.
!         B. HANSEN    ECMWF  MARCH    1997  SIGNAL HANDLING.
!         J. BIDLOT    ECMWF  April    1997  ADD ZDELATM IN PARAM. LIST
!         S. ABDALLA   ECMWF  OCTOBER  2001  INCLUSION OF AIR DENSITY & Zi/L
!                                            GENERALIZE THE INTERFACE WITH
!                                            THE ATMOSPHERIC MODEL
!         G.Mozdzynski ECMWF  January  2005  OPTIMISE COUPLING COMMS
!         J BIDLOT     ECMWF  August 2006  PASS RMISS TO SPECIFY MISSING DATA.
!         M. Drusch    ECMWF  Sep 2007  Re-initialize through FRSTIME and NADV
!         J BIDLOT     ECMWF  August 2008 ADD LWCUR, MOVE PREWIND.
!         J BIDLOT     ECMWF  June 2009 ADD LWSTOKES.
!         J BIDLOT     ECMWF  August 2010 ADD IGRIB_HANDLE.
!         P Bechtold   ECMWF  March 2012 ADD SMALL PLANET RADIUS&GRAVITY FACTOR
!

!     PURPOSE.
!     --------

!          THIS SUBROUTINE SUPERVISES THE EXECUTION OF
!          MAIN MODULES FOR WAM MODEL INITIALIZATION,
!          WIND FIELD PREPROCESSING, WAM MODEL EXECUTION,
!          AND WAVE DATA ASSIMILATION.

!*    INTERFACE.
!     ----------

!          SEE MAIN MODULES SUB INITMDL, PREWIND, WAMODEL, WAMASSI.

!     METHOD.
!     -------

!          THE FIRST TIME WAVEMDL IS CALLED, THE WAM MODEL IS
!          INITIALIZED. THIS INITIALIZATION INCLUDES GETTING
!          FROM ECFILE THE INITIAL SEA STATE FILES, FILLING
!          COMMON BLOCKS DEFINING THE GRID AND SETTING GENERAL
!          PARAMETERS. IN THE FIRST AND ALL SUBSEQUENT CALLS TO
!          WAVEMDL  PREWIND REFORMATS THE WINDS INTO THE WAM
!          MODEL BLOCKED STRUCTURE AND THE WAM MODEL IS EXECUTED.
!          EACH CALL TO WAMODEL INTEGRATES THE WAVE SPECTRA FORWARD
!          IN TIME BY ONE INPUT WIND TIME STEP OR PROPAGATION TIME
!          STEP, WHAT EVER IS GREATER.

!     EXTERNALS.
!     ----------

!          INITMDL  -  INITIALIZES THE WAM MODEL.
!                      GETS RECOVERY FILES OUT OF ECFILE,
!                      SETS COMMON BLOCKS NECESSARY TO DEFINE
!                      THE GRID AND BLOCKING STRUCTURE.
!                      DEFINES GENERAL PARAMETERS.

!          PREWIND  -  REFORMATS WINDS ON THE GAUSSIAN GRID
!                      INTO THE WAM MODEL BLOCKED STRUCTURE.

!          WAMODEL  -  INTEGRATES THE WAVE SPECTRA FORWARD IN TIME BY
!                      ONE WIND INPUT TIME STEP OR ONE PROPAGATION
!                      TIME STEP, WHATEVER IS GREATER.

!          WAMASSI  -  SUPERVISES DATA ASSIMILATION:
!                      PREPROCESSES DATA; PRODUCES ANALYSED
!                      INTEGRATED QUANTITIES BY OPTIMAL INTERPOLATION;
!                      ANALYSES WAVE SPECTRA ; SAVE ANALYSIS FOR
!                      OUTPUT AND NEXT TIME STEP MODEL COMPUTATION.

!     REFERENCES.
!     -----------

!          NONE

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : CASS     ,NASS
      USE YOWCOUP  , ONLY : LWCOU, LWCOU2W, LWCOUHMF, LWFLUX,           &
     &         LWCOUNORMS, LLNORMWAMOUT_GLOBAL, LLNORMWAM2IFS,          &
     &         KCOUSTEP, LMASK_OUT_NOT_SET, LMASK_TASK_STR,             &
     &         I_MASK_OUT, J_MASK_OUT, N_MASK_OUT, LWNEMOCOU,           &
     &         IFSTSTEP, IFSNSTEP,                                      &
     &         LIFS_IO_SERV_ENABLED
      USE YOWGRIBHD, ONLY : NDATE_TIME_WINDOW_END
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK 
      USE YOWCURR  , ONLY : IDELCUR  ,LLCHKCFL
      USE YOWFRED  , ONLY : FR
      USE YOWGRID  , ONLY : IJS      ,IJL
      USE YOWICE   , ONLY : FLMIN
      USE YOWMAP   , ONLY : BLK2GLO  ,BLK2LOC  ,ZDELLO   ,IQGAUSS,   AMONOP
      USE YOWMEAN  , ONLY : INTFLDS
      USE YOWNEMOFLDS , ONLY : WAM2NEMO, NEMO2WAM
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NANG     ,NFRE    ,     &
     &                      LL1D
      USE YOWPCONS , ONLY : ZMISS    ,G        ,GM1
      USE YOWPHYS  , ONLY : RNU      ,RNUM     ,PRCHAR
      USE YOWSTAT  , ONLY : MARSTYPE ,CDATEA   ,CDATEE   ,CDATEF   ,    &
     &            CDTPRO   ,IDELPRO  ,IDELWI   ,IDELWO   ,IASSI    ,    &
     &            LSMSSIG_WAM,CMETER ,CEVENT   ,                        &
     &            IDELWI_LST,IDELWO_LST,CDTW_LST,NDELW_LST
      USE YOWSHAL  , ONLY : WVENVI   ,WVPRPT
      USE YOWTEST  , ONLY : IU06
      USE YOWWNDG  , ONLY : ICODE_CPL
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FF_NOW   ,FL1 
      USE YOWWIND  , ONLY : CDAWIFL  ,IUNITW   ,CDATEWO  ,CDATEFL ,     &
     &                      FF_NEXT  ,                                  &
     &                      NXFFS    ,NXFFE    ,NYFFS    ,NYFFE,        &
     &                      NXFFS_LOC,NXFFE_LOC,NYFFS_LOC,NYFFE_LOC
      USE YOWUNPOOL,ONLY  : LLUNSTR

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOWASSI  , ONLY : WAMASSI
      USE YOWGRIB  , ONLY : IGRIB_GET_VALUE
      USE MPL_MODULE, ONLY : MPL_BARRIER, MPL_GATHERV
! ---------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "chkoops.intfb.h"
#include "difdate.intfb.h"
#include "incdate.intfb.h"
#include "initmdl.intfb.h"
#include "mpfldtoifs.intfb.h"
#include "outbeta.intfb.h"
#include "prewind.intfb.h"
#include "setmarstype.intfb.h"
#include "wamodel.intfb.h"

!     INITIAL DATE OF THE RUN
      CHARACTER(LEN=14), INTENT(IN) :: CBEGDAT
!     ATMOSPHERIC TIME STEP
      REAL(KIND=JWRB), INTENT(IN) :: PSTEP
!     NUMBER OF ATMOSPHERIC TIME STEPS UNTIL THE END OF THE RUN
      INTEGER(KIND=JWIM), INTENT(IN) :: KSTOP
!     NUMBER OF ATMOSPHERIC TIME STEP(S) FOR EACH WAVE MODEL CALL
      INTEGER(KIND=JWIM), INTENT(IN) :: KSTPW
!     NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
!     NUMBER OF ATMOSPHERIC GRID POINTS IN THOSE FIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
!     NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
!     NUMBER OF ATM. ROWS OF LATITUDES
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
!     IFS GRIB HANDLE
      INTEGER(KIND=JWIM), INTENT(IN) :: IGRIB_HANDLE
!     GRIB MISSING DATA INDICATOR
      REAL(KIND=JWRB), INTENT(IN) :: RMISS
!     DEFAULT VALUE FOR CHARNOCK
      REAL(KIND=JWRB), INTENT(IN) :: ZRCHAR
!     FIELDS CONTAINING ATMOSPHERIC DATA
      REAL(KIND=JWRB), INTENT(INOUT) :: FIELDS(NGPTOTG,NFIELDS)
!     INDICATES WHICH PHYSICAL QUANTITY IS USED FOR WIND INPUT
!      0: NOT DECIDED WILL BE DETERMINED FROM INPUT
!      1: FRICTION VELOCITY, 2: SURFACE STRESS, 3: 10m WIND
      INTEGER(KIND=JWIM), INTENT(IN) :: NATMFLX
!     INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS FROM IFS
      LOGICAL, INTENT(INOUT) :: LWCUR
!     INDICATES WHETHER THE PRODUCTION OF THE STOKES DRIFT IS REQUIRED FOR THE IFS
      LOGICAL, INTENT(INOUT) :: LWSTOKES
!     NUMBER OF FIELDS RETURNED TO ATMOSPHERIC MODEL
      INTEGER(KIND=JWIM), INTENT(IN) :: NWVFIELDS
!     FIELDS RETURNED TO ATMOSPHERIC MODEL
      REAL(KIND=JWRB), INTENT(INOUT) :: WVFLDG(:,:,:)
!     FIRST DIMENSION OF WVFLDG
      INTEGER(KIND=JWIM), INTENT(IN) :: NLONW
!     SECOND DIMENSION OF WVFLDG
      INTEGER(KIND=JWIM), INTENT(IN) :: NLATW
!     SET .TRUE. IF STOP SIGNAL RECEIVED.
      LOGICAL, INTENT(INOUT) :: LDSTOP
!     SET .TRUE. IF RESTART FILE SIGNAL RECEIVED.
      LOGICAL, INTENT(INOUT) :: LDWRRE
!     TELLS ATMOSPHERIC MODEL THAT IT WAS A RESTART
      LOGICAL, INTENT(OUT) :: LDRESTARTED
!     WAVE MODEL GRID SPACING FOR EACH LATITUDE
      REAL(KIND=JWRB), INTENT(OUT) :: ZDELATM(NLATW)
!     TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
      LOGICAL, INTENT(INOUT) :: LDWCOUNORMS
!     TELL WAM TO PROCUCE REPRODUCIBLE GLOBAL NORMS
      LOGICAL, INTENT(IN) :: LDNORMWAMOUT_GLOBAL
!     MASK TO INDICATE WHICH PART OF ARRAY FIELDS IS RELEVANT
      INTEGER(KIND=JWIM), INTENT(INOUT) :: MASK_IN(NGPTOTG)
!     MASK TO INDICATE WHICH PART OF ARRAY WVFLDG IS RELEVANT
      INTEGER(KIND=JWIM), INTENT(INOUT) :: MASK_OUT(NLONW,NLATW)
!     CONTROLS FIRST CALL
      LOGICAL, INTENT(INOUT) :: FRSTIME
!     NUMBER OF ADVECTION STEPS
      INTEGER(KIND=JWIM), INTENT(INOUT) :: NADV
!     MODIFICATION FACTOR FOR EARTH RADIUS FOR SMALL PLANET RUNS
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI
!     MODIFICATION FACTOR FOR GRAVITY FOR SMALL PLANET RUNS
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRG
!     KINEMATIC AIR DENSITY
      REAL(KIND=JWRB), INTENT(IN) :: RNU_ATM
!     REDUCED KINEMATIC AIR DENSITY FOR MOMENTUM TRANSFER
      REAL(KIND=JWRB), INTENT(IN) :: RNUM_ATM 
!     USED TO SPECIFY THE END OF THE 4DVAR ANALYSIS WINDOW WHEN COUPLED 
      INTEGER(KIND=JWIM), INTENT(IN) :: IDATE_TIME_WINDOW_END
!     ATMOSPHERIC NSTEP (CT3)
      INTEGER(KIND=JWIM), INTENT(IN) :: NSTEP
!     IFS IO SERVER ENABLED
      LOGICAL, INTENT(IN) :: LDIFS_IO_SERV_ENABLED

      INTEGER(KIND=JWIM) :: IJ, I, J, K, ICPLEN,ICPLEN_ECF
      INTEGER(KIND=JWIM) :: KDELWI, IDURAT
      INTEGER(KIND=JWIM) :: NDUR, KSTOP_BY, ISTOP
      INTEGER(KIND=JWIM) :: N_MASK_OUT_GLOBAL
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJSB, IJLB
      INTEGER(KIND=JWIM) :: IREAD, ISTAT, IRECV
      INTEGER(KIND=JWIM) :: IFLDOFFSET
      INTEGER(KIND=JWIM) :: IYYYYMMDD, IHHMM, ISTEP, IRET, KRET
      INTEGER(KIND=JWIM) :: IC, IL, IST, IED, ICOUNT, JF, IP, IFLD
      INTEGER(KIND=JWIM) :: NCOMBUF, NCOMLOC, NTOT, NMASK
      INTEGER(KIND=JWIM) :: IFCST, IFCSTEP_HOUR
      INTEGER(KIND=JWIM) :: NXS, NXE, NYS, NYE
      INTEGER(KIND=JWIM), ALLOCATABLE :: ZCOMCNT(:)

      REAL(KIND=JWRB) :: VAL
      REAL(KIND=JWRB) :: STEP
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB) :: DURATION, DURATION_MAX, PSTEP8

      REAL(KIND=JWRB), DIMENSION(NWVFIELDS) :: FAVG,FMIN,FMAX
      REAL(KIND=JWRB), DIMENSION(NWVFIELDS) :: DEFVAL
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NCHNK) :: BETAM, BETAB
      REAL(KIND=JWRB), ALLOCATABLE :: ZCOMBUFS(:), ZCOMBUFR(:)
      REAL(KIND=JWRB), ALLOCATABLE :: WVBLOCK(:,:)

      CHARACTER(LEN=7), PARAMETER :: CL_CPENV="SMSNAME"
      CHARACTER(LEN=8), PARAMETER :: CL_CPENV_ECF="ECF_NAME"

      CHARACTER(LEN=2) :: CCLASS, CTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=4) :: CEXPVER
      CHARACTER(LEN=12) :: C12
      CHARACTER(LEN=12) :: FLABEL(NWVFIELDS)
      CHARACTER(LEN=14) :: CDUM
      CHARACTER(LEN=14) :: CDTASS
      CHARACTER(LEN=40) :: CLSETEV
      CHARACTER(LEN=256) :: CLSMSNAME,CLECFNAME

      LOGICAL, SAVE :: LFRST
      LOGICAL, SAVE :: LFRSTCHK
      LOGICAL, SAVE :: LLGRAPI
      LOGICAL :: LLGLOBAL_WVFLDG
      LOGICAL :: LLINIT
      LOGICAL :: LLINIT_FIELDG

      DATA LFRST /.TRUE./
      DATA LFRSTCHK /.TRUE./
      DATA LLGRAPI /.TRUE./

! ---------------------------------------------------------------------

!*    1.  THE FIRST CALL TO WAVEMDL PERFORMS INITIALIZATION.
!         --------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAVEMDL',0,ZHOOK_HANDLE)

      NDATE_TIME_WINDOW_END=IDATE_TIME_WINDOW_END

      IFSTSTEP = PSTEP
      IFSNSTEP = NSTEP
      LIFS_IO_SERV_ENABLED = LDIFS_IO_SERV_ENABLED
      ZMISS=RMISS
      G=G*PRPLRG ! modified for small planet
      GM1 = 1.0_JWRB/G

      RNU=RNU_ATM
      RNUM=RNUM_ATM
      PRCHAR=ZRCHAR

      KDELWI=NINT(PSTEP)*KSTPW

      ICODE_CPL=NATMFLX

      IF (LWCOU) THEN
        IF ( IQGAUSS /= 1 ) THEN
          IF ( AMONOP < 90._JWRB ) THEN
              WRITE (IU06,*) ' *********************************'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' * PROBLEM IN WAVEMDL..........  *'
              WRITE (IU06,*) ' *   *'
              WRITE (IU06,*) ' * AMONOP SHOULD NOT BE < 90 IF  *'
              WRITE (IU06,*) ' * COUPLED AND ON LAT LON GRID   *'
              WRITE (IU06,*) ' * ============================= *'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' * AMONOP=', AMONOP
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' *********************************'
              CALL FLUSH(IU06)
              CALL ABORT1
          ENDIF
        ENDIF

        NGRIB_HANDLE_IFS=IGRIB_HANDLE


        IF (NGRIB_HANDLE_IFS < 0 ) THEN
          WRITE(IU06,*)' SUB: WAVEMDL:  NGRIB_HANDLE_IFS < 0 !'
          WRITE(IU06,*)' CALL ABORT1 '
          WRITE(IU06,*)'  '
          CALL ABORT1
        ENDIF

        IF (LLGRAPI) THEN
           WRITE(IU06,*)'  '
           WRITE (IU06,*) ' WAVEMDL: GRIB HANDLE FROM IFS'
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'dataDate',IYYYYMMDD)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'time',IHHMM)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'step',STEP)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'endStep',ISTEP)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'expver',C12,KRET=IRET)
           IF (IRET /= 0) THEN
             CEXPVER='****'
           ELSE
             CEXPVER=C12(1:4)
           ENDIF
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'class',C12)
           CCLASS=C12(1:2)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'stream',C12)
           CSTREAM=C12(1:4)
           CALL IGRIB_GET_VALUE(NGRIB_HANDLE_IFS,'type',C12)
           CTYPE=C12(1:2)
           WRITE(IU06,*)' EXPVER=', CEXPVER,   &
     &                  ' CLASS=', CCLASS,     &
     &                  ' STREAM=', CSTREAM,   &
     &                  ' TYPE=', CTYPE
           LLGRAPI=.FALSE.
        ENDIF

        LLNORMWAMOUT_GLOBAL = LLNORMWAMOUT_GLOBAL .OR. LDNORMWAMOUT_GLOBAL

      ENDIF

      KCOUSTEP=KDELWI

      IREAD=1
      IF (NPROC == 1) IREAD=1


      IF (FRSTIME) THEN

!     !!!! FIRST TIME AROUND !!!!

        LMASK_OUT_NOT_SET=.TRUE.
        LMASK_TASK_STR=.TRUE.

        IF (LWCOU) THEN
          IU06=20
          IL = LEN_TRIM(CBEGDAT)
          IF (IL /= 14) THEN
            WRITE (IU06,*) ' NON-Y2K COMPLIANT DATE IN CALLING WAVEMDL'
            WRITE (IU06,*) ' CBEGDAT =  ',CBEGDAT
            WRITE (IU06,*) ' PROGRAM WILL ABORT '
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF
          CDATEA = CBEGDAT
          CDATEE = CDATEA

          DURATION_MAX=HUGE(KSTOP)
          PSTEP8=PSTEP
          DURATION=REAL(KSTOP,KIND(DURATION))*PSTEP8
          NDUR=INT(DURATION/DURATION_MAX)+1
          KSTOP_BY=KSTOP/NDUR
          DO ISTOP=1,KSTOP,KSTOP_BY
            IDURAT=MIN(KSTOP_BY,KSTOP-ISTOP+1)*NINT(PSTEP)
           CALL INCDATE (CDATEE,IDURAT)
          ENDDO

          WRITE (IU06,1011)
          WRITE (IU06,1012)
          WRITE (IU06,1013)
          IF (LWCOU2W) THEN
            WRITE (IU06,1014)
          ELSE
            WRITE (IU06,1015)
          ENDIF
          WRITE (IU06,1011)
          WRITE (IU06,1016) CDATEA
          WRITE (IU06,1018) CDATEE
          WRITE (IU06,1019) NINT(PSTEP)
          WRITE (IU06,1020) KSTOP
          WRITE (IU06,1021) KDELWI
          WRITE (IU06,1011)
          WRITE (IU06,1010)
          CALL FLUSH(IU06)

 1010 FORMAT ('  **************************************************')
 1011 FORMAT ('  *                                                *')
 1012 FORMAT ('  *    WAVEMDL: INITMDL STARTS FOR (RE-)START.     *')
 1013 FORMAT ('  *    ========================================    *')
 1014 FORMAT ('  *    TWO-WAY INTERACTION WIND AND WAVES          *')
 1015 FORMAT ('  *    ONE-WAY INTERACTION WIND AND WAVES          *')
 1016 FORMAT ('  * START DATE OF RUN      ', A14,  ' (CDATEA) *')
 1018 FORMAT ('  * END DATE OF RUN        ', A14,  ' (CDATEE) *')
 1019 FORMAT ('  * IFS MODEL TIME STEP    ', I10,  ' (PSTEP)     *')
 1020 FORMAT ('  * TOTAL NUMBER OF PSTEP  ', I10,  ' (KSTOP)     *')
 1021 FORMAT ('  * INTERACTION TIME STEP  ', I10,  ' (KDELWI)    *')

        ENDIF

!       INQUIRE IF IUNITW IS ALREADY OPEN THEN CLOSE IT
        IF (IUNITW /= 0) CLOSE(IUNITW)

        CALL INITMDL (NADV,                                    &
     &                IREAD,                                   &
     &                BLK2GLO, BLK2LOC,                        &
     &                WVENVI, WVPRPT, FF_NOW,                  &
     &                FL1,                                     &
     &                NFIELDS, NGPTOTG, NC, NR,                &
     &                FIELDS, LWCUR, MASK_IN, PRPLRADI,        &
     &                NEMO2WAM) 


        LLCHKCFL=.FALSE.

        FRSTIME = .FALSE.                                              

        IF (LWCOU) THEN
          IF (NLONW /= NGX .OR. NLATW /= NGY) THEN
            WRITE (IU06,*) ' *********************************'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * PROBLEM IN WAVEMDL..........  *'
            WRITE (IU06,*) ' * PROBLEM WITH NLONW AND NLATW  *'
            WRITE (IU06,*) ' * NOT EQUAL TO NGX   AND NGY  : *'
            WRITE (IU06,*) ' * ============================= *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' * NLONW=',NLONW
            WRITE (IU06,*) ' * NLATW=',NLATW
            WRITE (IU06,*) ' * NGX=',NGX
            WRITE (IU06,*) ' * NGY=',NGY
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *                               *'
            WRITE (IU06,*) ' *********************************'
            CALL ABORT1
          ENDIF
          IF (LRESTARTED)  THEN
            LDRESTARTED = .TRUE.
!            RETURN
          ELSE
            LDRESTARTED = .FALSE.
          ENDIF
        ENDIF

      ELSE  ! .NOT. FRSTIME

!     !!!! ANY OTHER TIMES !!!!

!       CHECK IF THE INTERACTION/WIND-INPUT TIME STEP HAS NOT CHANGED
!       ALSO FOR SURFACE CURRENT INPUT,
!       OTHERWISE RESET NADV
        IF (LWCOU) THEN
          IF (IDELWI /= KDELWI) THEN
            IDELWO = KDELWI
            IDELWI = KDELWI
            IDELCUR= KDELWI
            NADV = IDELWI/IDELPRO
            CDAWIFL=CDTPRO
            CALL INCDATE(CDAWIFL,IDELWI)
            WRITE (IU06,*) ' ***********************************'
            WRITE (IU06,*) ' *                                 *'
            WRITE (IU06,*) ' * IN WAVEMDL :                    *'
            WRITE (IU06,*) ' * INTERACTION TIME STEP WAS RESET *'
            WRITE (IU06,*) ' * KDELWI = ', KDELWI
            WRITE (IU06,*) ' *                                 *'
            WRITE (IU06,*) ' ***********************************'
            CALL FLUSH(IU06)
            IF (IDELPRO*NADV /= IDELWI) THEN
              WRITE (IU06,*) ' ***************************************'
              WRITE (IU06,*) ' *                                     *'
              WRITE (IU06,*) ' * PROBLEM IN WAVEMDL :                *'
              WRITE (IU06,*) ' * THE NEW INTERACTION TIME STEP IS NOT*'
              WRITE (IU06,*) ' * A MULTIPLE OF IDELPRO  !!!          *'
              WRITE (IU06,*) ' * KDELWI = ', KDELWI
              WRITE (IU06,*) ' * IDELPRO = ', IDELPRO
              WRITE (IU06,*) ' *                                     *'
              WRITE (IU06,*) ' ***************************************'
              CALL ABORT1
            ENDIF
          ENDIF
        ELSE
          IF (NDELW_LST > 0) THEN
            DO IC=1,NDELW_LST
              IF (CDTPRO < CDTW_LST(IC)) THEN
                IF (IDELWI /= IDELWI_LST(IC) .OR. IDELWO /= IDELWO_LST(IC)) THEN

                  CALL INCDATE(CDATEWO,-IDELWO/2)
                  IDELWI=IDELWI_LST(IC)
                  IDELWO=IDELWO_LST(IC)
                  CALL INCDATE(CDATEWO,IDELWO/2)
                  CDATEFL=CDATEWO
                  NADV = IDELWI/IDELPRO
                  CDAWIFL=CDTPRO
                  CALL INCDATE(CDAWIFL,IDELWI)
                  WRITE (IU06,*) ' **********************************'
                  WRITE (IU06,*) ' *                                *'
                  WRITE (IU06,*) ' * IN WAVEMDL :                   *'
                  WRITE (IU06,*) ' * WIND INPUT TIME STEP WAS RESET *'
                  WRITE (IU06,*) ' * IDELWI = ', IDELWI
                  WRITE (IU06,*) ' * IDELWO = ', IDELWO
                  WRITE (IU06,*) ' *                                *'
                  WRITE (IU06,*) ' **********************************'
                  CALL FLUSH(IU06)
                ENDIF

                EXIT

              ENDIF
            ENDDO
          ENDIF

        ENDIF


!*      REFORMAT FORCING FIELDS FROM INPUT GRID TO BLOCKED.
!       ---------------------------------------------------
        IF (LWCOU) THEN
          NXS = NXFFS_LOC
          NXE = NXFFE_LOC
          NYS = NYFFS_LOC
          NYE = NYFFE_LOC

          IF (LFRSTCHK) THEN
!           CHECK THAT THE STRUCTURE OF FIELDS IS IN AGREEMENT WITH IFROMIJ AND JFROMIJ
            DO ICHNK = 1, NCHNK
              DO IJ = 1, KIJL4CHNK(ICHNK) 
                I = BLK2LOC(IJ,ICHNK)%IFROMIJ 
                J = BLK2LOC(IJ,ICHNK)%JFROMIJ
                IF (I < NXS .OR. I > NXE .OR. J < NYS .OR. J > NYE) THEN
                  WRITE(IU06,*) '*************ERROR********************'
                  WRITE(IU06,*) '* WAVEMDL: IJS to IJL SPAN TOO MUCH !'
                  WRITE(IU06,*) '* ICHNK, IJ = ', ICHNK, IJ
                  WRITE(IU06,*) '* I, J = ', I, J
                  WRITE(IU06,*) '* NXS, NXE  = ', NXS, NXE
                  WRITE(IU06,*) '* NYS, NYE  = ', NYS, NYE
                  WRITE(IU06,*) '*************ERROR********************'
                  CALL ABORT1
                ENDIF
              ENDDO
            ENDDO

            LFRSTCHK = .FALSE.
          ENDIF

        ELSE
          NXS = NXFFS
          NXE = NXFFE
          NYS = NYFFS
          NYE = NYFFE
        ENDIF

        LLINIT = .FALSE.
        LLINIT_FIELDG = .NOT. LWCOU
!       !!!! PREWIND IS CALLED THE FIRST TIME IN INITMDL !!!!
        CALL PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT,    &
                      NXS, NXE, NYS, NYE, LLINIT_FIELDG,   &
     &                LLINIT, IREAD,                       &
     &                NFIELDS, NGPTOTG, NC, NR,            &
     &                FIELDS, LWCUR, MASK_IN,              &
     &                NEMO2WAM) 

      ENDIF


! --------------------------------------------------------------------

!*    2.1  INTEGRATE THE WAVE SPECTRA FORWARD IN TIME.
!          -------------------------------------------

      IF (LWCOU) THEN
        CALL CHKOOPS(LDUPDATEOOPS=.TRUE.)
      ENDIF

      CALL SETMARSTYPE

      CALL WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,            &
     &              WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS, &
     &              WAM2NEMO, NEMO2WAM, FL1)


!*    2.2  DATA ASSIMILATION
!*         IF REQUESTED AND MODEL IS IN ANALYSIS PERIOD.
!          THE DATA ASSIMILATION SOFTWARE IS NOT AVAILABLE
!          FOR GENERAL DISSIMINATION !
!          ---------------------------------------------

      IF (IASSI == 1) THEN

        MARSTYPE = 'an'

!       UPDATE ANALYSIS TIME
        DO J=1,NASS
          IF (CDTPRO == CASS(J)) THEN
            CDTASS=CDTPRO
            EXIT
          ENDIF
        ENDDO

        IF (NASS > 0 ) THEN
          IF ( CDTPRO == CDTASS ) THEN
            CALL WAMASSI (LDSTOP, LDWRRE, BLK2GLO,          &
 &                        WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                        WAM2NEMO, NEMO2WAM, FL1)
          ENDIF
        ELSEIF ( (.NOT.LWCOU .AND. CDTPRO <= CDATEF ) .OR. (LWCOU .AND. CDTPRO == CDATEF) ) THEN
          CALL WAMASSI (LDSTOP, LDWRRE, BLK2GLO,          &
 &                      WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                      WAM2NEMO, NEMO2WAM, FL1)
        ENDIF
      ENDIF


!     ECFLOW METER:
      IF (LSMSSIG_WAM) THEN

        CALL MPL_BARRIER(CDSTRING='WAVEMDL:')

        CALL DIFDATE (CDATEF, CDTPRO, IFCST)
        IFCSTEP_HOUR=IFCST/3600
        IF (IRANK == 1) THEN
          WRITE(CLSETEV,' (A25,'' step '',I8,''&'') ') CMETER,IFCSTEP_HOUR
          CLSMSNAME="                                             "
          CLECFNAME="                                             "
          CALL GET_ENVIRONMENT_VARIABLE(NAME=CL_CPENV,     VALUE=CLSMSNAME, LENGTH=ICPLEN)
          CALL GET_ENVIRONMENT_VARIABLE(NAME=CL_CPENV_ECF, VALUE=CLECFNAME, LENGTH=ICPLEN_ECF)
          IF( ICPLEN     == 0 ) CLSMSNAME = 'NOSMS'
          IF( ICPLEN_ECF == 0 ) CLECFNAME = 'NOECF'
          IF ((ICPLEN > 0.AND.CLSMSNAME(1:5) /= 'NOSMS') .OR.           &
     &        (ICPLEN_ECF > 0.AND.CLECFNAME(1:5) /= 'NOECF') ) THEN
            CALL SYSTEM(CLSETEV)
            WRITE(IU06,'(2X,A25,I8,'' posted '')') CMETER,IFCSTEP_HOUR
          ELSE
            WRITE(IU06,'(A25,I8)') CMETER,IFCSTEP_HOUR
            WRITE(IU06,*) 'not posted  because neither SMSNAME'
            WRITE(IU06,*) ICPLEN, CLSMSNAME
            WRITE(IU06,*) 'nor ECF_NAME  is defined. '
            WRITE(IU06,*) ICPLEN_ECF, CLECFNAME
          ENDIF
        ENDIF
      ENDIF

!----------------------------------------------------------------------

!     3. PREPARE FIELDS THAT ARE RETURNED TO IFS (if coupled).
!        -----------------------------------------------------

      IF (LWCOU) THEN

        ALLOCATE(WVBLOCK(IJS:IJL,NWVFIELDS))

!       GRID LAYOUT
        DO K=1,NLATW
          ZDELATM(K) = ZDELLO(K)
        ENDDO

!       FIELDS TO BE PASSED TO THE ATMOSPHERIC MODEL ARE:

!       1. CHARNOCK FIELD(S)
        FLABEL(1)=' Charnock'
        DEFVAL(1)=PRCHAR ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                         ! THE WAVE MODEL ICE FREE SEA POINTS.
        IFLDOFFSET=1

        IF (LWCOUHMF) THEN
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' Eqv Chnk'
          DEFVAL(IFLDOFFSET)=0.25_JWRB*PRCHAR
        ENDIF

        IF (LWSTOKES) THEN
!         2. U-STOKESDRIFT
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' U-Stokes'
          DEFVAL(IFLDOFFSET)=0.0_JWRB  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.

!         3. V-STOKESDRIFT
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' V-Stokes'
          DEFVAL(IFLDOFFSET)=0.0_JWRB  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

        IF (LWFLUX) THEN
!         4. ENERGY FLUX TO OCEAN (dimensional)
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' Phi_ocd'
          DEFVAL(IFLDOFFSET)=0.0_JWRB  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.

!         5. X-COMPONENT OF THE MOMENTUM FLUX TO OCEAN (dimensional)
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' Tau_ocx'
          DEFVAL(IFLDOFFSET)=0.0_JWRB  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.

!         6. X-COMPONENT OF THE MOMENTUM FLUX TO OCEAN (dimensional)
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' Tau_ocy'
          DEFVAL(IFLDOFFSET)=0.0_JWRB  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.

!         7. WINDSEA VARIANCE
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' WSEmean'
          DEFVAL(IFLDOFFSET)=FLMIN  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                    ! THE WAVE MODEL ICE FREE SEA POINTS.

!         8. WINDSEA MEAN FREQUENCY
          IFLDOFFSET=IFLDOFFSET+1
          FLABEL(IFLDOFFSET)=' WSFmean'
          DEFVAL(IFLDOFFSET)=FR(NFRE)  ! DEFAULT VALUE FOR GRID POINTS NOT COVERED BY
                                       ! THE WAVE MODEL ICE FREE SEA POINTS.
        ENDIF

        CALL GSTATS(1443,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK, KIJS, KIJL, IJSB, IJLB, IFLDOFFSET, IFLD)
        DO ICHNK = 1, NCHNK
          CALL OUTBETA (1, NPROMA_WAM,                                                                          &
     &                 FF_NOW(:,ICHNK)%WSWAVE, FF_NOW(:,ICHNK)%UFRIC, FF_NOW(:,ICHNK)%Z0M, FF_NOW(:,ICHNK)%Z0B, &
     &                 FF_NOW(:,ICHNK)%CHRNCK, BETAM(:,ICHNK), BETAB(:,ICHNK))

          KIJS = 1
          IJSB = IJFROMCHNK(KIJS,ICHNK)
          KIJL = KIJL4CHNK(ICHNK)
          IJLB = IJFROMCHNK(KIJL,ICHNK) 

          IFLDOFFSET=1
          WVBLOCK(IJSB:IJLB,IFLDOFFSET)=BETAM(KIJS:KIJL,ICHNK)

          IF (LWCOUHMF) THEN
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=BETAB(KIJS:KIJL,ICHNK)
          ENDIF

!         SURFACE STOKES DRIFT NEEDED FOR THE IFS
!         IT MIGHT ALSO BE USED FOR NEMO !!!!!!!!!!

          IF (LWSTOKES) THEN
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%USTOKES
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%VSTOKES
          ENDIF

          IF (LWFLUX) THEN
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%PHIOCD
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%TAUOCXD
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%TAUOCYD
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%WSEMEAN
            IFLDOFFSET=IFLDOFFSET+1
            WVBLOCK(IJSB:IJLB,IFLDOFFSET)=INTFLDS(KIJS:KIJL,ICHNK)%WSFMEAN
          ENDIF

          DO IFLD=IFLDOFFSET+1,NWVFIELDS
            WVBLOCK(IJSB:IJLB,IFLD)=0.0_JWRB
          ENDDO

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1443,1)


!       GRIDDED FIELDS ARE NEEDED FOR IFS, GATHERING OF THE NECESSARY
!       INFORMATION
        CALL MPFLDTOIFS(IJS, IJL, BLK2GLO, NWVFIELDS, WVBLOCK,              &
     &                  WVFLDG, DEFVAL, MASK_OUT, LLGLOBAL_WVFLDG)


        DEALLOCATE(WVBLOCK)


!       COMPUTATION OF THE NORMS OF OUTPUT FIELDS
        IF (LLGLOBAL_WVFLDG) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' GLOBAL NORM OF FIELDS RETURNED TO IFS :'

          N_MASK_OUT_GLOBAL=NLATW*NLONW
          CALL GSTATS(1443,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF, J, I)
          DO JF = 1, NWVFIELDS
            FAVG(JF) = 0._JWRB
            FMIN(JF) = WVFLDG(1,1,JF)
            FMAX(JF) = WVFLDG(1,1,JF)
            DO J=1, NLATW
              DO I=1, NLONW
                FAVG(JF) = FAVG(JF) + WVFLDG(I,J,JF)
                FMIN(JF) = MIN(FMIN(JF), WVFLDG(I,J,JF))
                FMAX(JF) = MAX(FMAX(JF), WVFLDG(I,J,JF))
              ENDDO
            ENDDO
            FAVG(JF)=FAVG(JF)/N_MASK_OUT_GLOBAL
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)

          DO IFLD = 1, NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD), FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),N_MASK_OUT_GLOBAL
            WRITE(IU06,111) FAVG(IFLD),FMIN(IFLD),FMAX(IFLD)
          ENDDO
111       FORMAT(14x,'HEX: ',3(Z16.16,2x))
          WRITE(IU06,*) ' '
          CALL FLUSH(IU06)

        ELSEIF (LLNORMWAM2IFS) THEN
!         Local NORM OF FIELDS RETURNED TO IFS :
!         INITIALISE I,J POINTER TO MASK_OUT
          IF (LFRST) THEN
            LFRST=.FALSE.
            JF=1
            N_MASK_OUT=0
            DO J=1,NLATW
              DO I=1,NLONW
                IF (MASK_OUT(I,J) == 1) THEN
                  N_MASK_OUT=N_MASK_OUT+1
                ENDIF
              ENDDO
            ENDDO
            ALLOCATE(I_MASK_OUT(MAX(1,N_MASK_OUT)))
            I_MASK_OUT(1)=1
            ALLOCATE(J_MASK_OUT(MAX(1,N_MASK_OUT)))
            J_MASK_OUT(1)=1
            JF=1
            IC=0
            DO J=1,NLATW
              DO I=1,NLONW
                IF (MASK_OUT(I,J) == 1) THEN
                  IC=IC+1
                  I_MASK_OUT(IC)=I
                  J_MASK_OUT(IC)=J
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!         COMPUTE local NORMS
          CALL GSTATS(1443,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF, VAL, IC)
          DO JF = 1, NWVFIELDS
            VAL = WVFLDG(I_MASK_OUT(1), J_MASK_OUT(1),JF)
            FAVG(JF) = VAL
            FMIN(JF) = VAL
            FMAX(JF) = VAL
            DO IC = 2, N_MASK_OUT
              VAL = WVFLDG(I_MASK_OUT(IC), J_MASK_OUT(IC),JF)
              FAVG(JF) = FAVG(JF) + VAL
              FMIN(JF) = MIN(FMIN(JF), VAL)
              FMAX(JF) = MAX(FMAX(JF), VAL)
            ENDDO
            FAVG(JF) = FAVG(JF)/MAX(N_MASK_OUT,1)
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1443,1)

!         FOR PRIMARY PE (IRECV), COLLECT ALL THE NORMS TO PRODUCE A
!         PSEUDO GLOBAL NORM.

          IRECV=1
          NCOMLOC=1+3*NWVFIELDS
          NCOMBUF=NCOMLOC*NPROC
          ALLOCATE(ZCOMBUFS(NCOMBUF))

          IST=1+(IRANK-1)*NCOMLOC
          IED=IST+NCOMLOC-1
          ICOUNT=IST
          ZCOMBUFS(ICOUNT)=N_MASK_OUT
          DO IFLD=1,NWVFIELDS
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FAVG(IFLD)
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FMIN(IFLD)
            ICOUNT=ICOUNT+1
            ZCOMBUFS(ICOUNT)=FMAX(IFLD)
          ENDDO

          IF (IRANK == IRECV) THEN
            ALLOCATE(ZCOMBUFR(NCOMBUF))
          ENDIF
          ALLOCATE(ZCOMCNT(NPROC))
          ZCOMCNT=NCOMLOC
          CALL MPL_GATHERV(PSENDBUF=ZCOMBUFS(IST:IED),KROOT=IRECV,      &
     &                    PRECVBUF=ZCOMBUFR(:),KRECVCOUNTS=ZCOMCNT,     &
     &                    CDSTRING='WAVEMDL:')
          DEALLOCATE(ZCOMCNT)

!         COMPUTE PSEUDO GLOBAL NORM
          IF (IRANK == IRECV) THEN
            WRITE(IU06,*) ' '
            WRITE(IU06,*) ' Local NORM OF global FIELDS RETURNED TO IFS :'
            IST=1+(IRANK-1)*NCOMLOC
            ICOUNT=IST
            NMASK=ZCOMBUFR(ICOUNT)
            NTOT=NMASK
            DO JF=1,NWVFIELDS
              ICOUNT=ICOUNT+1
              FAVG(JF)=NMASK*ZCOMBUFR(ICOUNT)
              ICOUNT=ICOUNT+2
            ENDDO
            DO IP=1,NPROC
              IF (IP /= IRECV) THEN
                IST=1+(IP-1)*NCOMLOC
                ICOUNT=IST
                NMASK=ZCOMBUFR(ICOUNT)
                NTOT=NTOT+NMASK
                DO JF=1,NWVFIELDS
                  ICOUNT=ICOUNT+1
                  FAVG(JF)=FAVG(JF)+NMASK*ZCOMBUFR(ICOUNT)
                  ICOUNT=ICOUNT+1
                  FMIN(JF)=MIN(FMIN(JF),ZCOMBUFR(ICOUNT))
                  ICOUNT=ICOUNT+1
                  FMAX(JF)=MAX(FMAX(JF),ZCOMBUFR(ICOUNT))
                ENDDO
              ENDIF
            ENDDO
            DEALLOCATE(ZCOMBUFR)
            DO JF=1,NWVFIELDS
              FAVG(JF)=FAVG(JF)/MAX(NTOT,1)
            ENDDO
          ELSE
            NTOT=N_MASK_OUT
            WRITE(IU06,*) ' '
            WRITE(IU06,*) ' Local NORM OF FIELDS RETURNED TO IFS :'
          ENDIF

          DEALLOCATE(ZCOMBUFS)

          DO IFLD=1,NWVFIELDS
            WRITE(IU06,*) FLABEL(IFLD), FAVG(IFLD),FMIN(IFLD),FMAX(IFLD),NTOT, &
     &                    IRANK, NPROC, LL1D
            WRITE(IU06,111) FAVG(IFLD),FMIN(IFLD),FMAX(IFLD)
          ENDDO

        ENDIF  ! end of norm production


!       TELL ATMOS MODEL THE WAM REQUIREMENT FOR GLOBAL NORMS
        LDWCOUNORMS=LWCOUNORMS
      ENDIF ! end lwcou

!     4. END OF RUN ?
!        -----------
      IF (CDATEE == CDTPRO) THEN
        CALL MPL_BARRIER(CDSTRING='WAVEMDL: END')
        CALL FLUSH(IU06)
      ENDIF

      IF (LHOOK) CALL DR_HOOK('WAVEMDL',1,ZHOOK_HANDLE)

END SUBROUTINE WAVEMDL
