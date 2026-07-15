! (C) Copyright 1989- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,             &
 &                  WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,  &
 &                  WAM2NEMO, NEMO2WAM, FL1)

! ----------------------------------------------------------------------

!**** *WAMODEL* - 3-G WAM MODEL - WRAPPER FOR TIME INTEGRATION OF WAVE FIELDS
!                                 AND OUTPUTS

!*    PURPOSE.
!     --------

!       COMPUTATION OF THE 2-D FREQUENCY-DIRECTION WAVE SPECTRUM AT ALL
!       GRID POINTS FOR A GIVEN INITIAL SPECTRUM AND FORCING SURFACE
!       STRESS FIELD.

!**   INTERFACE.
!     ----------

!     *CALL* *WAMODEL (NADV, LDSTOP, LDWRRE, BLK2GLO,
!    &                 WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,
!    &                 WAM2NEMO, NEMO2WAM, FL1)
!        *NADV*      NUMBER OF ADVECTION ITERATIONS
!                    PER CALL OF WAMODEL, OUTPUT PARAMETER.
!        *LDSTOP*    SET .TRUE. IF STOP SIGNAL RECEIVED.
!        *LDWRRE*    SET .TRUE. IF RESTART SIGNAL RECEIVED.
!        *BLK2GLO*   BLOCK TO GRID TRANSFORMATION
!        *WVENVI*    WAVE ENVIRONMENT FIELDS
!        *WVPRPT*    WAVE PROPERTIES FIELDS
!        *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!        *FF_NEXT*   DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!        *INTFLDS*   INTEGRATED/DERIVED PARAMETERS
!        *WAM2NEMO*  WAVE FIELDS PASSED TO NEMO
!        *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, WAVE2OCEAN, OCEAN2WAVE

      USE YOWCPBO  , ONLY : IBOUNC   ,GBOUNC  , IPOGBO  , CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,                                  &
     &                      LWNEMOCOU,                                  &
     &                      NEMOWSTEP, NEMOFRCO     ,                   &
     &                      NEMOCSTEP, NEMONSTEP
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,FFLAG20  ,GFLAG20  ,    &
     &                      NGOUT    ,                                  &
     &                      NIPRMOUT ,                                  &
     &                      LFDB     ,NOUTT    ,NOUTS    ,              &
     &                      CASS     ,NASS     ,LOUTINT  ,              &
     &                      LRSTPARALW, LRSTINFDAT,                     &
     &                      LRSTST0  ,LWAMANOUT
      USE YOWCURR  , ONLY : CDTCUR
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE
      USE YOWMESPAS, ONLY : LFDBIOOUT,LGRIBOUT ,LNOCDIN  ,LWAVEWIND
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &                      CDATER   ,CDATES   ,CDTINTT  ,IDELPRO  ,IDELT    ,    &
     &                      IDELWI   ,IREST    ,IDELRES  ,IDELINT  ,              &
     &                      CDTBC    ,IDELBC   ,                                  &
     &                      IASSI    ,MARSTYPE ,                                  &
     &                      LLSOURCE ,LANAONLY ,LFRSTFLD ,IREFDATE
      USE YOWSPEC, ONLY   : NBLKS    ,NBLKE
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : ICPLEN   ,CPATH    ,CWI      ,LRESTARTED
      USE YOWUNIT  , ONLY : IU02     ,IU19     ,IU20
      USE YOWWAMI  , ONLY : CBPLTDT  ,CEPLTDT  ,IANALPD  ,IFOREPD  ,    &
     &                      IDELWIN  ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CDATEWO

      USE MPL_MODULE, ONLY : MPL_BARRIER
      USE WAM_MULTIO_MOD, ONLY : WAM_MULTIO_FLUSH
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "outwint.intfb.h"
#include "outwpsp.intfb.h"
#include "abort1.intfb.h"
#include "bouinpt.intfb.h"
#include "difdate.intfb.h"
#include "gsfile_new.intfb.h"
#include "headbc.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "incdate.intfb.h"
#include "outbc.intfb.h"
#include "outbs.intfb.h"
#include "outspec.intfb.h"
#include "outstep0.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"
#include "unsetice.intfb.h"
#include "updnemofields.intfb.h"
#include "updnemostress.intfb.h"
#include "wamintgr.intfb.h"
#include "writsta.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN)                                           :: NADV
      LOGICAL, INTENT(INOUT)                                                   :: LDSTOP, LDWRRE
      TYPE(WVGRIDGLO), INTENT(IN)                                              :: BLK2GLO
      TYPE(ENVIRONMENT), INTENT(INOUT)                                         :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT)                                           :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)                                      :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(IN)                                         :: FF_NEXT
      TYPE(INTGT_PARAM_FIELDS), INTENT(INOUT)                                  :: INTFLDS
      TYPE(WAVE2OCEAN), INTENT(INOUT)                                          :: WAM2NEMO
      TYPE(OCEAN2WAVE), INTENT(IN)                                             :: NEMO2WAM
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1


      INTEGER(KIND=JWIM) :: IJ, K, M, J, IRA, KADV, ICH
      INTEGER(KIND=JWIM) :: IFIL, IC, ICL, ICR, II, ILOOP
      INTEGER(KIND=JWIM) :: ICHNK
      INTEGER(KIND=JWIM) :: JSTPNEMO, IDATE, ITIME
      INTEGER(KIND=JWIM) :: IU04
      INTEGER(KIND=JWIM), DIMENSION(NPROMA_WAM, NCHNK) :: MIJ

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, MAX(NIPRMOUT,1), NCHNK) :: BOUT
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK) :: XLLWS

      CHARACTER(LEN= 2) :: MARSTYPEBAK
      CHARACTER(LEN=14) :: CDATEWH, CZERO
      CHARACTER(LEN=14) :: CDATE, CDTPRA, CDTIMP, CDTIMPNEXT, CDTRCF

      LOGICAL :: LLFLUSH
      LOGICAL :: LSV, LRST, LOUT
      LOGICAL :: LLNONASSI

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WAMODEL',0,ZHOOK_HANDLE)

!     0.0 INITIALISATION
!         --------------

      CZERO = ' '
      LLFLUSH = .FALSE.
      KTAG = 200
      LRSTST0 = .FALSE.

!     TIME FOR THE NEXT SOURCE TERM INTEGRATION
      CDTIMPNEXT = CDTPRO
      CALL INCDATE(CDTIMPNEXT, IDELT)
!     TIME FOR WIND INPUT UPDATE (SEE NEWWIND)
      CDTIMP = CDTPRO

!     0.1 MINIMUM ENERGY
!         --------------
      IF (CDTPRO == CDATEA .AND. LLSOURCE ) THEN
!       INSURE THERE IS SOME WAVE ENERGY FOR GRID POINTS THAT HAVE BEEN
!       FREED FROM SEA ICE (ONLY DONE INITIALLY AND IF THE MODEL IS NOT
!       RESTARTED).
!       IT ALSO RESETS THE MIMIMUM ENERGY LEVEL THAT MIGHT HAVE BEEN LOST
!       WHEN GETTING THE DATA FROM GRIB.
        CALL GSTATS(1236,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          CALL UNSETICE(1, NPROMA_WAM, WVENVI%DEPTH(:,ICHNK), WVENVI%EMAXDPT(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK), &
 &                      FF_NOW%WSWAVE(:,ICHNK), FF_NOW%CICOVER(:,ICHNK), FL1(:,:,:,ICHNK) )
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1236,1)
      ENDIF


!     0.2 OUTPUT INITIAL CONDITION AND/OR FORECAST STEP 0
!         -----------------------------------------------
      IF (CDTPRO == CDATEA .OR. CDTPRO == CDATEF) THEN
         CALL OUTSTEP0 (WVENVI, WVPRPT, FF_NOW, INTFLDS,  &
 &                      WAM2NEMO, NEMO2WAM, FL1)
      ENDIF




!*    1. ADVECTION/PHYSICS TIME LOOP.
!        ----------------------------

      ADVECTION : DO KADV = 1,NADV

!*      1.1 FIX END DATE OF THIS PROPAGATION STEP AND OUTPUT TIMES.
!           -------------------------------------------------------

        CDTPRA = CDTPRO
        CALL INCDATE(CDTPRO, IDELPRO)

!       UPDATE OUTPUT TIMES.
        IF (NOUTT > 0) THEN
          CDTINTT = CZERO
          DO J=1,NOUTT
            IF (CDTPRO == COUTT(J)) THEN
              IF (FFLAG20 .OR. GFLAG20) CDTINTT = COUTT(J)
            ENDIF
          ENDDO
        ELSE
          IF ((FFLAG20.OR.GFLAG20) .AND. CDTINTT.LT.CDTPRO) CALL INCDATE (CDTINTT, IDELINT)
        ENDIF

!       UPDATE SPECTRA OUTPUT DATE
        IF (NOUTS > 0) THEN
!         reset CDATES to insure that spectra output is only controlled
!         by list COUTS
          CDATES='000000000000'
          DO J=1,NOUTS
            IF (CDTPRO == COUTS(J)) THEN
              CDTRES=CDTPRO
              CDATER=CDTRES
              CDATES=CDTRES
              EXIT
            ENDIF
          ENDDO
        ELSE
          IF (CDTRES < CDTPRO) CALL INCDATE(CDTRES, IDELRES)
        ENDIF

!NEST (not used at ECMWF)
        IF ((IBOUNC == 1 .OR. IBOUNF == 1) .AND. CDTBC < CDTPRO) CALL INCDATE(CDTBC, IDELBC)
!NEST



!            -------------------------------------------------------
!*      1.2  THIS IS THE CORE OF THE WAVE MODEL:
!*           COMPUTATION OF PROPAGATION
!*           INTEGRATION OF SOURCE TERMS OVER SUB TIME STEPS BETWEEN
!*           PROPAGATION TIME STEPS.
!            -------------------------------------------------------
!       SET TIME COUNTER.
        CDATE   = CDTPRA
        CDATEWH = CDATEWO
        ILOOP = 1
        DO WHILE ( ILOOP == 1 .OR. CDTIMPNEXT <= CDTPRO)
          CALL WAMINTGR (CDTPRA, CDATE, CDATEWH, CDTIMP, CDTIMPNEXT, &
 &                       BLK2GLO,                                    &
 &                       WVENVI, WVPRPT, FF_NOW, FF_NEXT, INTFLDS,   &
 &                       WAM2NEMO, MIJ, FL1, XLLWS)
          ILOOP = ILOOP +1
        ENDDO


!       1.3 CHECK WHETHER OUTPUT(s) NEEDED
!           ------------------------------
        LRST = (LDWRRE .AND. KADV == NADV )
        IF (LRST) THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  ******************************************'
          IF (LDSTOP) THEN
          WRITE(IU06,*) '  AN INTERRUPT SIGNAL HAS BEEN RECEIVED '
          ENDIF
          WRITE(IU06,*) '  THE NECESSARY BINARY RESTART FILES WILL BE'
          WRITE(IU06,*) '  GENERATED.'
          WRITE(IU06,*) '  ******************************************'
          WRITE(IU06,*) ' '
          CALL FLUSH (IU06)
        ENDIF


!NEST (not used at ECMWF)
!*      1.4.1 INPUT OF BOUNDARY VALUES.
!           -------------------------
        IF (IBOUNF == 1) CALL BOUINPT (IU02, FL1, NBLKS, NBLKE)
!*      1.4.2 OUTPUT OF BOUNDARY POINTS.
!           --------------------------
        IF (IBOUNC == 1) CALL OUTBC (FL1, BLK2GLO, IU19)
!NEST


!*      1.5 POINT OUTPUT (not usually used at ECMWF)
!           ----------------------------------------
        IF ( NGOUT > 0 .AND. (CDTINTT == CDTPRO .OR. LRST) ) THEN
!           OUTPUT POINT SPECTRA (not usually used at ECMWF)
            CALL OUTWPSP (FL1, FF_NOW)
        ENDIF


!       1.6 COMPUTE OUTPUT PARAMETERS FIELDS AND PRINT OUT NORMS
!           ----------------------------------------------------
        IF ( (CDTINTT == CDTPRO .OR. LRST) .AND. NIPRMOUT > 0 ) THEN
          CALL OUTBS (MIJ, FL1, XLLWS,                             &
     &                WVPRPT, WVENVI, FF_NOW, INTFLDS, NEMO2WAM,   &
     &                BOUT)
        ENDIF


!*      1.7 ONE PROPAGATION TIMESTEP DONE
!           -----------------------------
        WRITE(IU06,*) ' !!!!!!!!!!!!!! WAVE FIELDS INTEGRATED FOR DATE : ', CDTPRO

!       IS THIS AN ANALYSIS TIME FROM THE INPUT LIST ?
        LLNONASSI=.TRUE.
        IF (IASSI == 1) THEN
          DO J=1,NASS
            IF (CDTPRO == CASS(J)) THEN
              LLNONASSI=.FALSE.
              EXIT
            ENDIF
          ENDDO
        ENDIF


!*      1.8 SAVE FIELD OF 2D SPECTRA AND/OR BINARY RESTART FILES.
!           -----------------------------------------------------
!         THIS WILL HAPPEN WHEN IT IS NOT IN DATA ASSIMILATION MODE AND
!         IT IS EITHER A DETERMINED OUTPUT TIME
!         OR THE INTERUPT SIGNAL HAS BEEN TRIGGERED and it will wait
!         until the end of the advection loop.
!         OTHERWISE THE OUTPUT WILL OCCUR IN WAMASSI.

        LSV=(CDTRES == CDTPRO .OR. CDATEE == CDTPRO .OR. CDTPRO == CDATER)

        IF (LSV .OR. LRST) THEN

          LOUT = ((IREST == 1) .AND. (CDTPRO == CDATER .OR. CDTPRO <= CDATES)) .AND. LSV .AND. LWAMANOUT

          IF ( LOUT .OR. LRST ) THEN

!           1.8.1 SAVE SPECTRUM IN GRIB
!                 ---------------------
            IF (LOUT .AND. LGRIBOUT) THEN
!             we have insured that the spectra will be written to FDB
!             even when the restart option is triggered and it is an
!             output step for the spectra.

!             IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!             BECOMES TYPE AN (i.e. pseudo analysis)
              MARSTYPEBAK=MARSTYPE
              IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.             &
     &            (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
                MARSTYPE='an'
              ENDIF

              CALL OUTSPEC(FL1, FF_NOW)
              LLFLUSH = .TRUE.

              MARSTYPE=MARSTYPEBAK

              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  GRIB WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
            ENDIF

!           1.8.2 SAVE RESTART FILES IN PURE BINARY FORM (in needed)
!                 --------------------------------------
            IF ( .NOT.LGRIBOUT .OR. LDWRRE ) THEN

              CALL SAVSTRESS(WVENVI, FF_NOW, NBLKS, NBLKE, CDTPRO, CDATEF)
              WRITE(IU06,*) ' '
              WRITE(IU06,*) '  BINARY STRESS FILE DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '

              CALL SAVSPEC(FL1, NBLKS, NBLKE, CDTPRO, CDATEF, CDATER)
              WRITE(IU06,*) '  BINARY WAVE SPECTRA DISPOSED AT........ CDTPRO  = ', CDTPRO
              WRITE(IU06,*) ' '
              CALL FLUSH(IU06)
            ENDIF


!*          1.8.3 UPDATE, WRITE AND SAVE WAMINFO FILE.
!                 -----------------------------------
            IF (LRST .AND. IRANK == 1) THEN
              ICH = 7
              CALL DIFDATE (CDATEF, CDATEE, IFOREPD)
              IF (CDTPRO <= CDATEF) THEN
                CALL DIFDATE (CDTPRO, CDATEF, IANALPD)
                CBPLTDT = CDTPRO
                NFCST = 1
              ELSE
                NFCST = 0
                IANALPD = 0
                CBPLTDT = CDATEF
                CALL DIFDATE (CDTPRO, CDATEE, IFOREPD)
              ENDIF
              ISTAT(:) = 0
              IF (CDATE == CDATEE) ISTAT(1) = 1
              IDELWIN = IDELWI

              CEPLTDT = CDATEF

              IU04 = IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8) , 'w', 'f', 0, 'READWRITE')

              CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,     &
     &                      IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,  &
     &                      IASSI, NFCST, ISTAT, CDTCUR,                &
     &                      LRSTPARALW, NPROC)

              CLOSE (IU04)
              WRITE(IU06,*) ' WAMINFO FILE WRITTEN FOR RESTART... CDTPRO  = ', CDTPRO
              WRITE(IU06,*) '                                     CDATEF  = ', CDATEF
              WRITE(IU06,*) ' TO ', CWI(1:ICPLEN+8)
              CALL FLUSH(IU06)

              IF (LRSTINFDAT) THEN
!               WRITE AN ADDITIONAL wamfile WITH DATE/TIME INFO added to filename
                CDTRCF=CDTPRO
                CALL INCDATE(CDTRCF, IDELPRO)
                IU04 =  IWAM_GET_UNIT (IU06,CWI(1:ICPLEN+8)//'.'//      &
     &                              CDTRCF(1:8)//'_'//CDTRCF(9:14),     &
     &                              'w', 'f', 0, 'READWRITE')
                CALL WRITSTA (IU04, CDTPRO, CDATEE, IANALPD, IFOREPD,   &
     &                        IDELWIN, CDATER, CDATES, CBPLTDT, CEPLTDT,&
     &                        IASSI, NFCST, ISTAT, CDTCUR,              &
     &                        LRSTPARALW, NPROC)
                CLOSE (IU04)
              ENDIF

            ENDIF  ! waminfo

          ENDIF  ! END SAVE RESTART FIELDS

        ENDIF  ! END SAVE 1.8


!       1.9 WRITE INTEGRATED PARAMETER DATA TO FDB OR TO FILE
!           -------------------------------------------------
        IF (LWAMANOUT .AND. CDTINTT == CDTPRO .AND. NIPRMOUT > 0 ) THEN
!         IF THE OUTPUT TIME IS NOT AN ANALYSIS TIME THEN TYPE FG or 4V
!         BECOMES TYPE AN (i.e. speudo analysis)
          MARSTYPEBAK=MARSTYPE
          IF ((MARSTYPE == 'fg' .AND. KADV < NADV) .OR.                 &
     &        (MARSTYPE == '4v' .AND. LLNONASSI) ) THEN
            MARSTYPE='an'
          ENDIF

          CALL OUTWINT(BOUT)
          LLFLUSH = .TRUE.

          MARSTYPE=MARSTYPEBAK

          CALL GSTATS(753,0)
          CALL MPL_BARRIER(CDSTRING='WAMODEL:')
          CALL GSTATS(753,1)
        ENDIF


!*      1.10 FLUSH FDB IF IT HAS BEEN USED AND IT IS NOT AN ANALYSIS (it will be done in *wamassi*)
!            -------------------------------------------------------

        IF ( LFDB .AND. LLFLUSH .AND. (IASSI /= 1 .OR. CDTPRO > CDATEF) ) THEN
          CALL GSTATS(1976,0)
          CALL WAM_MULTIO_FLUSH()
          CALL GSTATS(1976,1)
          WRITE(IU06,*) ' '
          WRITE(IU06,*) '  FDB FLUSHED AT ',  CDTPRO, ' FROM WAMODEL. '
          CALL FLUSH (IU06)
          LLFLUSH=.FALSE.
        ENDIF


!       1.11 FOR PURE BINARY DATA (obsolete option at ECMWF !!!!!!!):
!*           OUTPUT FILES AND RECOVERY FILES ARE DISPOSED WHEN
!            TIME REACHES THE DISPOSE DATE OR WHEN THE MODEL
!            HAS BEEN SIGNALLED TO DO SO.
!            -------------------------------------------------
        IF (FFLAG20) THEN
          IF (CDATEE == CDTPRO .AND. LOUTINT .AND. IRANK == 1  ) THEN
            CALL GSFILE (IU06, IU20, 0, CDTPRO, CDATEF, 'MAP', 'S')
          ENDIF
        ENDIF

!NEST (not used at ECMWF)
!       SAVE BOUNDARY VALUE FILE.
        IF (CDTBC == CDTPRO) THEN
          IF (IBOUNC == 1 .AND. IRANK == 1 ) THEN
            DO II=1,GBOUNC
            CALL GSFILE(IU06, IU19(II), 0, CDTBC, CDTBC,                &
     &        CBCPREF(II), 'S')
            IF (CDTBC < CDATEE)                                         &
     &        CALL HEADBC (IPOGBO(II)-IPOGBO(II-1), IDELPRO,            &
     &                     TH(1), FR(1), IU19(II), IU06)
            ENDDO
          ENDIF
        ENDIF
!NEST


!*      1.12 WAM-NEMO COUPLING (!!!!! WHEN NO atmospheric model !!!!!!) (currently not used at ECMWF)
!       (when coupled see cnt4 in ifs)
!       ----------------------------------------------------
        IF (LWNEMOCOU .AND. (.NOT.LWCOU)) THEN
          NEMOWSTEP=NEMOWSTEP+1

          IF (MOD(NEMOWSTEP,NEMOFRCO) == 0) THEN
            CALL UPDNEMOFIELDS
            CALL UPDNEMOSTRESS

#ifdef WITH_NEMO
            DO JSTPNEMO = NEMOCSTEP, NEMOCSTEP+NEMONSTEP-1
               ! Advance the NEMO model 1 time step
               CALL NEMOGCMCOUP_STEP( JSTPNEMO, IDATE, ITIME )
               WRITE(IU06,*)'NEMO TIME IS : ',JSTPNEMO, IDATE, ITIME
            ENDDO
#endif
            NEMOCSTEP = NEMOCSTEP + NEMONSTEP
          ENDIF
        ENDIF


!*    BRANCHING BACK TO 1.0 FOR NEXT PROPAGATION STEP.
      ENDDO ADVECTION

IF (LHOOK) CALL DR_HOOK('WAMODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WAMODEL
