! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INITMDL (NADV,                                 &
 &                  IREAD,                                &
 &                  BLK2GLO,  BLK2LOC,                    &
 &                  WVENVI, WVPRPT, FF_NOW,               &
 &                  FL1,                                  &
 &                  NFIELDS, NGPTOTG, NC, NR,             &
 &                  FIELDS, LWCUR, MASK_IN, PRPLRADI,     &
 &                  NEMO2WAM)

! ----------------------------------------------------------------------

!**** *INITMDL* - INITIALIZES THE WAM MODEL.

!     L. ZAMBRESKY   GKSS/ECMWF    JULY 1988

!     MODIFIED BY:   H. GUNTHER    NOVEMBER 1989
!                    J. BIDLOT     FEBRUARY 1996-1997
!                    J. DOYLE      SEPTEMBER 1996
!                    J. BIDLOT     APRIL 1997
!                    J. BIDLOT     SEPTEMBER : MODIFY PARALLEL INPUT
!                    S. ABDALLA    OCTOBER 2001: INCLUSION OF AIR
!                                                DENSITY AND Zi/L
!                    J. BIDLOT     AUGUST 2008. ADD CALL TO PREWIND
!                                  TO GET CURRENTS EARLIER
!                                  (i.e BEFORE CALL TO GETSPEC).

!SHALLOW
!          DIFFERENCES FOR SHALLOW WATER RUNS TO DEEP WATER RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CSHALLOW'.
!SHALLOW
!NEST
!          DIFFERENCES FOR NESTED GRID RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CNEST'.
!NEST
!REFRA
!          DIFFERENCES FOR REFRACTION RUNS TO NORMAL RUNS
!          ARE ENCLOSED IN COMMENT LINES : 'CREFRA'.
!REFRA

!*    PURPOSE.
!     --------

!       INITIALIZE THE WAM MODEL.

!**   INTERFACE.
!     ----------

!          ---- FORMAL PARAMETERS ----

!    *CALL* *INITMDL (NADV,
!    &                IREAD,
!    &                BLK2GLO, BLK2LOC, 
!    &                WVENVI, WVPRPT, FF_NOW,
!    &                FL1,                  
!    &                NFIELDS, NGPTOTG, NC, NR,
!    &                FIELDS, LWCUR, MASK_IN,
!    &                NEMO2WAM)*

!      *NADV*      NUMBER OF ADVECTION ITERATIONS
!      *IREAD*     PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                  (IF NEEDED).
!      *BLK2GLO*   BLOCK TO GRID TRANSFORMATION
!      *BLK2LOC*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *WVENVI*    WAVE ENVIRONMENT FIELDS
!      *WVPRPT*    WAVE PROPERTIES FIELDS
!      *FF_NOW*    FORCING FIELDS AT CURRENT TIME.
!      *FL1*       SPECTRUM
!      *NFIELDS*   NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!      *NGPTOTG*   NUMBER OF ATMOSPHERIC GRID POINTS
!      *NC*        NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!      *NR*        NUMBER OF ATM. ROWS OF LATITUDES
!      *FIELDS*    ATM. FIELDS (U10, V10, AIR DENSITY, Zi/L, U and V CURRENTS)
!      *LWCUR*     INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!      *MASK_IN*   MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.
!      *PRPLRADI*  EARTH RADIUS REDUCTION FACTOR FOR SMALL PLANET
!      *NEMO2WAM*  FIELDS FRON OCEAN MODEL to WAM

!          ---- INPUT/OUTPUT UNITS ---

!          THE NAMES ARE DEFINED IN SECTION 1. OF THIS PROGRAM,
!          IF IT IS NOT MENTIONED OTHERWISE.

!           *IU02*   - INPUT  UNIT OF BOUDARY VALUES FROM A PREVIOUS
!                      COARSE GRID IF THIS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'FBI'
!                      (OUTPUT OF BOUINT).
!           *IU06*   - PRINTER OUTPUT.
!           *IU07*   - INPUT  UNIT OF PRECOMPUTED GRID PARAMETERS.
!                      (OUTPUT OF PREPROC).
!           *IU08*   - INPUT  UNIT OF MODULE YOWUBUF.
!                      (OUTPUT OF PREPROC).
!NEST
!           *IU09*   - INPUT  UNIT MODULE YOWCPBO (OUTPUT OF PREPROC).
!           *IU10*   - INPUT  UNIT MODULE YOWFPBO (OUTPUT OF PREPROC).
!NEST
!           *IU11*   - INPUT UNIT OF SPECTRA AT ALL GRID POINTS.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU11 AND IU12 ARE INTERCHANGED.
!           *IU12*   - OUTPUT UNIT BLOCKS OF SPECTRA (SEE IU11).
!           *IU13*   - INPUT UNIT OF SPECTRA ON LAST LATITUDE
!                      OF A BLOCK. SPECTRA ARE SAVED FROM THE
!                      SECOND LATITUDE OF NEXT BLOCK.
!                      EACH PROPAGATION STEP THE FILES CONNECTED
!                      TO IU13 AND IU14 ARE INTERCHANGED.
!           *IU14*   - OUTPUT UNIT SECOND LATITUDES (SEE IU13).
!           *IU15*   - OUTPUT UNIT LAST WINDFIELDS.
!NEST
!           *IU19*   - OUTPUT UNIT OF BOUNDARY VALUES IF
!                      THIS IS A FINE GRID RUN.
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'CBO'
!NEST
!           *IU20*   - OUTPUT UNIT OF INTEGRATED PARAMETERS
!                      OF THE TOTAL SPECTRUM (FIRST GUESS).
!                      THIS FILE IS DYNAMICALLY ASSIGNED FILEID = 'MAP'

!          FOR A START THE RESTART FILES WILL DYNAMICALLY BE ASSIGNED
!          AND COPIED TO THE *OUTPUT* UNITS (IU12, IU14, IU15).
!          IF IT IS REQUESTED TO SAVE RESTART FILES THESE WILL BE COPIED
!          IN REGULAR INTERVALS FROM THE UNIT ALIAS FILES OF IU11 OR
!          IU12, IU13 OR IU14, AND IU15 TO THE PREMANENT RESTART FILES.
!          FOR DETAILS OF THE FILE NAMES SEE SUB GSFILE.

!          THE PROGRAM CLOSES AND DELETES ALL WORK FILES.

!          ALL PARAMETERS HAVE TO BE THE VALUES GIVEN AT THE END
!          OF THE PREPROC OUTPUT IN COLUMN 'REQUIRED'.

!     METHOD.
!     -------

!          THIS ROUTINE INITIALISES THE WAVEMODEL:
!            -  DEFINES THE UNITS FOR INPUT/OUTPUT,
!            -  READS THE USER INPUT FILE,
!            -  INITIALIZES SOME MODEL CONSTANTS,
!            -  GETS THE RECOVERY FILES,
!            -  READS THE YOWMON BLOCKS PRECOMPUTED BY PROG PREPROC,
!            -  DOES SOME GENERAL BOOKEEPING REGARDING
!               DATES, INTEGRATION TIME STEPS AND OUTPUT TIME STEPS.
!            -  READ MODULE YOWUBUF AND SPECTRA IF ONE BLOCK VERSION.

!     REFERENCE
!     ---------

!          A MORE DETAILED  DISCUSSION MAY BE FOUND IN SUB WAMODEL.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO, WVGRIDLOC,                    &
     &                         ENVIRONMENT, FREQUENCY, FORCING_FIELDS,  &
     &                         INTGT_PARAM_FIELDS, OCEAN2WAVE

      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC   ,                        &
     &                      GBOUNC  , IPOGBO   ,CBCPREF
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP ,LWFLUX   ,LWNEMOCOU
      USE YOWCOUT  , ONLY : COUTT    ,COUTLST  ,FFLAG20  ,GFLAG20  ,    &
     &                      NGOUT    ,NOUTT    ,LOUTINT  ,LSECONDORDER
      USE YOWCURR  , ONLY : CDTCUR   ,IDELCUR  ,CDATECURA
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I, NGRIB_HANDLE_WAM_I2 ,NGRIB_HANDLE_WAM_S
      USE YOWFRED  , ONLY : FR       ,TH       ,DELTH   ,FR5      ,     &
     &            ZPIFR    ,                                            &
     &            FRM5     ,COFRM4   ,COEF4    ,FRATIO  ,FLOGSPRDM1,    &
     &            FLMAX    ,RHOWG_DFIM,                                 &
     &            DFIM     ,DFIMOFR  ,DFIMFR  ,DFIMFR2   ,              &
     &            DFIM_SIM ,DFIMOFR_SIM ,DFIMFR_SIM ,DFIMFR2_SIM ,      &
     &            DFIM_END_L, DFIM_END_U,                               &
     &            WVPRPT_LAND
      USE YOWGRIBHD, ONLY : LGRHDIFS
      USE YOWGRID  , ONLY : DELPHI, DELLAM, COSPH,                      &
     &                      NPROMA_WAM, NCHNK, IJFROMCHNK  
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,                        &
     &            AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,    &
     &            NIBLO    ,NGX      ,NGY      ,                        &
     &            KMNOP    ,KMSOP    ,IPER     ,IRGG     ,IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,KTAG
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,NFRE_ODD ,    & 
     &                      LLUNSTR
      USE YOWPCONS , ONLY : G        ,CIRC     ,PI       ,ZPI      ,    &
     &                      RAD      ,ROWATER  ,ZPI4GM2  ,FM2FP
      USE YOWPHYS  , ONLY : ALPHAPMAX, ALPHAPMINFAC, FLMINFAC
      USE YOWREFD  , ONLY : LLUPDTTD
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTHA   ,DEPTHD   ,TOOSHALLOW
      USE YOWSPEC  , ONLY : NBLKS    ,NBLKE    ,KLENTOP  ,KLENBOT
      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDTPRO   ,CDTRES   ,    &
     &            CDTINTT  ,CDTBC    ,                                  &
     &            IFRELFMAX, DELPRO_LF, IDELPRO  ,IDELT ,               &
     &            IDELWI   ,IDELWO   ,IDELRES  ,IDELINT  ,              &
     &            IREFRA   ,                                            &
     &            IPHYS    ,                                            &
     &            CDATEA   ,MARSTYPE ,LANAONLY ,ISNONLIN ,IPROPAGS ,    &
     &            IDELWI_LST,IDELWO_LST,CDTW_LST,NDELW_LST
      USE YOWTABL  , ONLY : FAC0     ,FAC1     ,FAC2     ,FAC3     ,    &
     &                      FAK      ,FRHF     ,DFIMHF
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWWNDG  , ONLY : ICODE
      USE YOWUBUF  , ONLY : LUPDTWGHT
      USE YOWUNIT  , ONLY : IU02     ,IU11     ,IU12     ,IU13     ,    &
     &                      IU14     ,IU15     ,IU19     ,IU20
      USE YOWWAMI  , ONLY : CBPLTDT
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEWO  ,CDATEFL  ,    &
     &                      NXFFS    , NXFFE   ,NYFFS    , NYFFE   ,    &
     &                      LLNEWCURR,LLWSWAVE ,LLWDWAVE ,FF_NEXT

#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL, ONLY : OUT_METHOD
      USE OUTPUT_STRUCT, ONLY : INITIAL_OUTPUT_INITS
      USE UNSTRUCT_BOUND , ONLY : IOBP
#endif
      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! -------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "cigetdeac.intfb.h"
#include "cireduce.intfb.h"
#include "getspec.intfb.h"
#include "getstress.intfb.h"
#include "headbc.intfb.h"
#include "incdate.intfb.h"
#include "inisnonlin.intfb.h"
#include "init_sdiss_ardh.intfb.h"
#include "init_x0tauhf.intfb.h"
#include "initdpthflds.intfb.h"
#include "initnemocpl.intfb.h"
#include "iniwcst.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "mcout.intfb.h"
#include "secondhh_gen.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "prewind.intfb.h"
#include "readbou.intfb.h"
#include "setmarstype.intfb.h"
#include "tabu_swellft.intfb.h"
#include "userin.intfb.h"

      INTEGER(KIND=JWIM), INTENT(OUT) :: NADV
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      TYPE(WVGRIDGLO), INTENT(IN) :: BLK2GLO
      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT) :: WVENVI
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NANG, NFRE, NCHNK), INTENT(INOUT) :: FL1
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      REAL(KIND=JWRB),DIMENSION(NGPTOTG,NFIELDS), INTENT(IN) :: FIELDS
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT) :: MASK_IN
      REAL(KIND=JWRB), INTENT(IN) :: PRPLRADI
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM


      INTEGER(KIND=JWIM) :: IFORCA
      INTEGER(KIND=JWIM) :: IJ, I, II, K, M, IP, LFILE, IX, IY, KX, ID
      INTEGER(KIND=JWIM) :: IC, ICR
      INTEGER(KIND=JWIM) :: JD
      INTEGER(KIND=JWIM) :: IDELWH
      INTEGER(KIND=JWIM) :: IU05, IU09, IU10
      INTEGER(KIND=JWIM) :: ICHNK, IPRM, KIJS, KIJL
      INTEGER(KIND=JWIM) :: MTHREADS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS


      REAL(KIND=JWRB) :: FCRANGE, XD
      REAL(KIND=JWRB) :: GVE, DPH, CFLP, CFLL, DLH, DLH_KX
      REAL(KIND=JWRB) :: FAC, SCDF_L, SCDF_U
      REAL(KIND=JWRB) :: D, OM, XK
      REAL(KIND=JWRB) :: XLOGFRATIO
      REAL(KIND=JWRB) :: GAM
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      REAL(KIND=JWRB) :: XLA, XLO

      CHARACTER(LEN=14) :: ZERO, CDUM
      CHARACTER(LEN=24) :: FILNM
      CHARACTER(LEN=80) :: FILENAME

      LOGICAL :: LLEXIST
      LOGICAL :: LLINIT
      LOGICAL :: LLINIT_FIELDG

!----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INITMDL',0,ZHOOK_HANDLE)


      CALL INIWCST(PRPLRADI)

!*    1. DEFINITION OF MODEL PARAMETERS.
!        -------------------------------

      ZERO = ' '

!*    1.1  DEFINE UNIT NAMES.
!          ------------------

!NEST
      IU02 = 0
!NEST


! GRID AND UBUF FILES UNITS
! -------------------------

!    DEFINITION HAS BEEN MOVED TO MPDECOMP

      IU11 = 11


! RESTART FILES UNITS
! -------------------

! NOTE IF MESSAGE PASSING
! NOTE IU12 WILL ONLY BE CONNECTED TO UNIT 12 FOR PURE BINARY INPUT
! NOTE OTHERWISE IT WILL CONNECTED TO THE FILENAME BLS WITH
! NOTE AN EXTENSION FUNCTION OF THE PE.

      IF (LWCOU) THEN
        IU12 = 36
        IU15 = 37
      ELSE
        IU12 = 12
        IU15 = 15
      ENDIF

      IU13 = 13
      IU14 = 14

      IU20 = 20


      LOUTINT=.FALSE.

      LUPDTWGHT=.TRUE.

      LLNEWCURR=.TRUE.

      ICODE=0

!*    1.2 INPUT OF USER PARAMETER.
!         ------------------------

      CALL USERIN (IFORCA, LWCUR)

!     DEFINE COEFFICIENT FOR MEAN PERIODS CALCULATION
      IF (ALLOCATED(DFIMOFR)) DEALLOCATE(DFIMOFR)
      ALLOCATE(DFIMOFR(NFRE))
      IF (ALLOCATED(DFIMFR)) DEALLOCATE(DFIMFR)
      ALLOCATE(DFIMFR(NFRE))
      IF (ALLOCATED(DFIMFR2)) DEALLOCATE(DFIMFR2)
      ALLOCATE(DFIMFR2(NFRE))

      DO M=1,NFRE
        DFIMOFR(M) = DFIM(M)/FR(M)
        DFIMFR(M)  = DFIM(M)*FR(M)
        DFIMFR2(M) = DFIM(M)*FR(M)**2
      ENDDO

!     DEFINE A FEW CONSTANTS FOR USE IN IMPLSCH

      IF (.NOT.ALLOCATED(ZPIFR)) ALLOCATE(ZPIFR(NFRE))
      IF (.NOT.ALLOCATED(FR5)) ALLOCATE(FR5(NFRE))
      IF (.NOT.ALLOCATED(FRM5)) ALLOCATE(FRM5(NFRE))
      IF (.NOT.ALLOCATED(COFRM4)) ALLOCATE(COFRM4(NFRE))
      IF (.NOT.ALLOCATED(FLMAX)) ALLOCATE(FLMAX(NFRE))
      IF (.NOT.ALLOCATED(DFIM_END_L)) ALLOCATE(DFIM_END_L(NFRE))
      IF (.NOT.ALLOCATED(DFIM_END_U)) ALLOCATE(DFIM_END_U(NFRE))

      SCDF_L = 0.5_JWRB*DELTH*(FRATIO-1.0_JWRB)
      SCDF_U = 0.5_JWRB*DELTH*(1.0_JWRB-1.0_JWRB/FRATIO)

      DO M=1,NFRE
        ZPIFR(M) = ZPI*FR(M)
        FR5(M) = FR(M)**5
        FRM5(M) = 1.0_JWRB/FR5(M)
        COFRM4(M) = COEF4*G/FR(M)**4
        FLMAX(M) = (ALPHAPMAX/PI)/(ZPI4GM2*FR5(M))
        DFIM_END_L(M) = SCDF_L*FR(M)
        DFIM_END_U(M) = SCDF_U*FR(M)
      ENDDO

      FLMINFAC = ALPHAPMINFAC*FM2FP*G/(PI*ZPI**3*FR(NFRE)**5)

      FLOGSPRDM1=1.0_JWRB/LOG10(FRATIO)

      XLOGFRATIO=LOG(FRATIO)

      IF (.NOT.ALLOCATED(RHOWG_DFIM)) ALLOCATE(RHOWG_DFIM(NFRE))
      RHOWG_DFIM(1) = 0.5_JWRB*ROWATER*G*DELTH*XLOGFRATIO*FR(1)
      DO M=2,NFRE-1
        RHOWG_DFIM(M) = ROWATER*G*DELTH*XLOGFRATIO*FR(M)
      ENDDO
      RHOWG_DFIM(NFRE) = 0.5_JWRB*ROWATER*G*DELTH*XLOGFRATIO*FR(NFRE)

      IF (.NOT.ALLOCATED(DFIM_SIM)) ALLOCATE(DFIM_SIM(NFRE))
      NFRE_ODD=NFRE-1+MOD(NFRE,2)
      DFIM_SIM(NFRE)=0.0_JWRB
      DFIM_SIM(1)=DELTH*XLOGFRATIO*FR(1)/3.0_JWRB
      DO M=2,NFRE_ODD-1,2
        DFIM_SIM(M)=4.0_JWRB*DELTH*XLOGFRATIO*FR(M)/3.0_JWRB
        DFIM_SIM(M+1)=2.0_JWRB*DELTH*XLOGFRATIO*FR(M+1)/3.0_JWRB
      ENDDO
      DFIM_SIM(NFRE_ODD)=DELTH*XLOGFRATIO*FR(NFRE_ODD)/3.0_JWRB

      IF (.NOT.ALLOCATED(DFIMOFR_SIM)) ALLOCATE(DFIMOFR_SIM(NFRE))
      IF (.NOT.ALLOCATED(DFIMFR_SIM)) ALLOCATE(DFIMFR_SIM(NFRE))
      IF (.NOT.ALLOCATED(DFIMFR2_SIM)) ALLOCATE(DFIMFR2_SIM(NFRE))
      DO M=1,NFRE
        DFIMOFR_SIM(M) = DFIM_SIM(M)/FR(M)
        DFIMFR_SIM(M)  = DFIM_SIM(M)*FR(M)
        DFIMFR2_SIM(M) = DFIM_SIM(M)*FR(M)**2
      ENDDO


      CALL TABU_SWELLFT

      CALL INIT_X0TAUHF

      KTAG=100

!     1.5 DETERMINE LAST OUTPUT DATE
!         --------------------------

      IF (NOUTT > 0) THEN
        COUTLST=COUTT(NOUTT)
      ELSE
        COUTLST=CDATEA
        IF (GFLAG20) THEN
          DO WHILE (COUTLST <= CDATEE .AND. IDELINT > 0)
            CALL INCDATE (COUTLST, IDELINT)
          ENDDO
          CALL INCDATE (COUTLST, -IDELINT)
        ENDIF

        IF (COUTLST > CDATEE) COUTLST=CDATEE
      ENDIF

! ----------------------------------------------------------------------

!*    2. INPUT FROM PREPROCESSING PROGRAMS.
!        ----------------------------------

!     THE ACTUAL READING HAS BEEN MOVED TO MPDECOMP.


!*    2.0 COMPUTE OUTPUT POINT INDICES (MODULE YOWCOUT).
!         ----------------------------------------------

      IF ( NGOUT > 0 .AND. .NOT. LLUNSTR ) CALL MCOUT


!*    2.1 READ MODULE YOWCPBO AND YOWFPBO.
!     ------------------------------------

      GBOUNC = 1

      IF (.NOT. LLUNSTR) THEN

      IF (IBOUNC == 1 .OR. IBOUNF == 1) THEN
        IF (IBOUNC == 1) THEN
!         READ INFORMATION ABOUT WHERE THE FINE GRID(S) ARE
          FILENAME='wam_nested_grids_info'
          LFILE=0
          LLEXIST=.FALSE.
          IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
          INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            WRITE(IU06,*) '************************************'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITMDL     *'
            WRITE(IU06,*) '*  =============================   *'
            WRITE(IU06,*) '*  WITH OPTION IBOUNC = 1          *'
            WRITE(IU06,*) '*  YOU MUST PROVIDE THE FILE:      *'
            WRITE(IU06,*) '* ',FILENAME(1:LFILE)
            WRITE(*,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),   &
     &        ' IS MISSING !!!!'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '************************************'
            CALL ABORT1
          ENDIF
          IU09 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0, 'READWRITE')
        ELSE
          IU09 = 9
        ENDIF

        IF (IBOUNF == 1) THEN
!         READ INFORMATION ABOUT WHERE THE BOUNDARY VALUES ARE
          FILENAME='wam_boundary_grid_info'
          LFILE=0
          LLEXIST=.FALSE.
          IF (FILENAME.NE. ' ') LFILE=LEN_TRIM(FILENAME)
          INQUIRE(FILE=FILENAME(1:LFILE),EXIST=LLEXIST)
          IF (.NOT. LLEXIST) THEN
            WRITE(IU06,*) '************************************'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITMDL     *'
            WRITE(IU06,*) '*  =============================   *'
            WRITE(IU06,*) '*  WITH OPTION IBOUNF = 1          *'
            WRITE(IU06,*) '*  YOU MUST PROVIDE THE FILE:      *'
            WRITE(IU06,*) '* ',FILENAME(1:LFILE)
            WRITE(NULERR,*) '*  WAVE MODEL INPUT FILE ',FILENAME(1:LFILE),' IS MISSING !!!!'
            WRITE(IU06,*) '*                                  *'
            WRITE(IU06,*) '************************************'
            CALL ABORT1
          ENDIF
          IU10 = IWAM_GET_UNIT(IU06, FILENAME(1:LFILE) , 'r', 'u', 0, 'READWRITE')
        ELSE
          IU10 = 10
        ENDIF

        CALL READBOU (IU09, IU10, IU06)

        CLOSE (UNIT=IU09, STATUS='KEEP')
        CLOSE (UNIT=IU10, STATUS='KEEP')
      ENDIF

      ENDIF ! .NOT. LLUNSTR

      IF (.NOT.ALLOCATED(IU19)) ALLOCATE(IU19(GBOUNC))


!*    2.2.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!          ------------------------------------
      LANAONLY=.FALSE.
      IF ((CDATEA == CDATEE) .AND. (CDATEA == CDATEF)) LANAONLY=.TRUE.

      CALL SETMARSTYPE

      IF (.NOT. LGRHDIFS) THEN
!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I2, NGRIBV=2)
!       FOR SPECTRA
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)
      ENDIF

      IF (MARSTYPE == 'cf' .OR. MARSTYPE == 'pf' .OR. MARSTYPE == 'fc' ) THEN
        IF (ALLOCATED(FAC0))   DEALLOCATE(FAC0)
        IF (ALLOCATED(FAC1))   DEALLOCATE(FAC1)
        IF (ALLOCATED(FAC2))   DEALLOCATE(FAC2)
        IF (ALLOCATED(FAC3))   DEALLOCATE(FAC3)
        IF (ALLOCATED(FAK))    DEALLOCATE(FAK)
        IF (ALLOCATED(FRHF))   DEALLOCATE(FRHF)
        IF (ALLOCATED(DFIMHF)) DEALLOCATE(DFIMHF)
      ENDIF

! ----------------------------------------------------------------------

!*    3. PRINT INITIAL CONDITIONS AS READ FROM PERPROCESSING.
!        ----------------------------------------------------

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' WAVE MODEL GRID ORGANISATION:'
      WRITE(IU06,3002) ' SOUTHERNMOST LATITUDE IN GRID IS .......: ', AMOSOP, ' DEGREE'
      WRITE(IU06,3002) ' NORTHERNMOST LATITUDE IN GRID IS .......: ', AMONOP, ' DEGREE'
      WRITE(IU06,3002) ' WESTERNMOST LONGITUDE IN GRID IS .......: ', AMOWEP, ' DEGREE'
      WRITE(IU06,3002) ' EASTERNMOST LONGITUDE IN GRID IS .......: ', AMOEAP, ' DEGREE'
      IF ( IQGAUSS == 1 ) THEN
        WRITE(IU06,*) ' GAUSSIAN GRID ..........................: '
        WRITE(IU06,3002) ' APPROXIMATE LATITUDE INCREMENT IS ......: ', XDELLA, ' DEGREE'
      ELSE
        IF ( IRGG == 1 ) THEN
          WRITE(IU06,*) ' IRREGULAR LAT/LON  GRID ................: '
          WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ', XDELLA, ' DEGREE'
        ELSE
          WRITE(IU06,*) ' LAT/LON  GRID ..........................: '
          WRITE(IU06,3002) ' LATITUDE INCREMENT IS ..................: ', XDELLA, ' DEGREE'
          WRITE(IU06,3002) ' LONGITUDE INCREMENT IS .................: ', XDELLO, ' DEGREE'
        ENDIF
      ENDIF
      WRITE(IU06,*) '  '
      WRITE(IU06,3003) ' TOTAL LENGTH OF EACH BLOCK .............: ', NIBLO
      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' SPECTRAL RESOLUTION:'
      WRITE(IU06,3003) ' TOTAL NUMBER OF DIRECTIONS .............: ', NANG 
      WRITE(IU06,3003) ' TOTAL NUMBER OF FREQUENCIES I/O & PROPAG: ', NFRE_RED
      WRITE(IU06,3003) ' TOTAL NUMBER OF FREQUENCIES FOR  PHYSICS: ', NFRE

      MTHREADS=1
!$    MTHREADS=OMP_GET_MAX_THREADS()

      WRITE(IU06,*) '  '
      WRITE(IU06,*) ' PARALLEL ORGANISATION : '
      WRITE(IU06,*) ' NPROC      : ', NPROC
      WRITE(IU06,*) ' MTHREADS   : ', MTHREADS
      WRITE(IU06,*) ' NPROMA_WAM : ', NPROMA_WAM
      WRITE(IU06,*) ' NCHNK      : ', NCHNK
      WRITE(IU06,*) '  '
      CALL FLUSH (IU06)

 3002 FORMAT(3x,a,f9.3,a)
 3003 FORMAT(3x,a,i8,  a)

!NEST
      IF (IBOUNC == 1) THEN
        WRITE(IU06,*)
        WRITE(IU06,*) ' COARSE GRID: BOUNDARY OUTPUT POINTS :'
        WRITE(IU06,*) ' TOTAL NUMBER OF BOUNDARY POINTS IS: ',NBOUNC
      ENDIF
!NEST

!*    3.1 COMPUTE SHALLOW WATER TABLE INDICES.
!         ------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IPRM, IJ, XD, ID)
      DO ICHNK = 1, NCHNK
        KIJS = 1
        KIJL = NPROMA_WAM

!       BUILD DEPTH POINTER AND RESET DEPTH TO ALL POSITIVE
        DO IPRM = KIJS, KIJL
          IF (WVENVI%DEPTH(IPRM,ICHNK) <= TOOSHALLOW) THEN
            WVENVI%IODP(IPRM,ICHNK) = 0
          ELSE
            WVENVI%IODP(IPRM,ICHNK) = 1
          ENDIF
        ENDDO

        IF (.NOT. LLUNSTR) THEN
          WVENVI%IOBND(KIJS:KIJL,ICHNK) = 1
        ELSE
#ifdef WAM_HAVE_UNWAM
          DO IPRM = KIJS, KIJL
            IJ = IJFROMCHNK(IPRM, ICHNK)
            IF (IJ /= 0) THEN
              IF (IOBP(IJ) /= 0) THEN
                WVENVI%IOBND(IPRM,ICHNK) = 0
              ELSE
                WVENVI%IOBND(IPRM,ICHNK) = 1
              ENDIF
            ENDIF
          ENDDO
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
        ENDIF

!       SET DEPTH MINIMUM
        DO IPRM = KIJS, KIJL
          WVENVI%DEPTH(IPRM,ICHNK) = MAX(WVENVI%DEPTH(IPRM,ICHNK), DEPTHA)
        ENDDO

!       COMPUTE INDEP
         DO IPRM = KIJS, KIJL
          XD = LOG(WVENVI%DEPTH(IPRM,ICHNK)/DEPTHA)/LOG(DEPTHD)+1.0_JWRB
          ID = NINT(XD)
          ID = MAX(ID,1)
          WVENVI%INDEP(IPRM,ICHNK) = MIN(ID, NDEPTH)
         ENDDO
      ENDDO
!$OMP END PARALLEL DO


!     INITIALISE GRID POINT FIELDS DEPENDENT ON WATER DEPTH AND FREQUENCY
      IF (.NOT.ALLOCATED(WVPRPT_LAND%WAVNUM)) CALL WVPRPT_LAND%ALLOC(NFRE)
      CALL INITDPTHFLDS(WVENVI, WVPRPT, WVPRPT_LAND)


!     INITIALISATION FOR SDISS_ARD
      IF (IPHYS == 1) CALL INIT_SDISS_ARDH

      WRITE(IU06,*) '  '
      WRITE(IU06,*) '  SUB. INITMDL: end arrays initialisation'
      WRITE(IU06,*) '  '
      CALL FLUSH (IU06)

! ----------------------------------------------------------------------

!*    4. CONNECT RESTART FIELDS TO OUTPUT UNITS (IF PBIO SOFTWARE NOT
!        USED). AND READ IN LAST WINDFIELDS FROM RESTARTFILE.
!        ---------------------------------------------------------------

      CALL GETSTRESS(BLK2LOC, WVENVI, FF_NOW, NEMO2WAM,   &
     &               NBLKS, NBLKE, IREAD) 

      WRITE(IU06,*) '  SUB. INITMDL: integrated paramaters READ IN'
      WRITE(IU06,*) '  '
      CALL FLUSH (IU06)

      CDATEWL = CDTPRO

      IF (LWCOU) LLWSWAVE = .FALSE.
      IF (LWCOU) LLWDWAVE = .FALSE.

      IF (CDTPRO /= ZERO .OR. LRESTARTED) THEN

!*    4.1 MODEL STARTS FROM FILES OUT OF A PREVIOUS MODEL RUN.
!         ----------------------------------------------------

        IF (CDTPRO < CDATEA .OR. CDTPRO > CDATEE .OR.                 &
     &      (IFORCA == 1 .AND. CDTPRO > CDATEF) .OR.                  &
     &      (IFORCA /= 1 .AND. CDTPRO <= CDATEF)       ) THEN
          WRITE(IU06,*) ' *******************************************'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. INITMDL          *'
          WRITE(IU06,*) ' *    ===========================          *'
          WRITE(IU06,*) ' * START DATE FROM RESTART FIELD IS NOT    *'
          WRITE(IU06,*) ' * MODEL PERIOD.                           *'
          IF (IFORCA == 1) THEN
            WRITE(IU06,*) ' *  IN ANALYSIS PERIOD AS REQUESTED.       *'
          ELSE
            WRITE(IU06,*) ' *  IN FORECAST PERIOD AS REQUESTED.       *'
          ENDIF
          WRITE(IU06,*) ' * START DATE OF RUN       IS CDATEA = ', CDATEA
          WRITE(IU06,*) ' * START DATE OF FORECAST  IS CDATEF = ', CDATEF
          WRITE(IU06,*) ' * END   DATE OF RUN       IS CDATEE = ', CDATEE
          WRITE(IU06,*) ' * START DATE FROM RESTART IS CDTPRO = ', CDTPRO
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' * PROGRAM ABORTS     PROGRAM ABORTS       *'
          WRITE(IU06,*) ' *                                         *'
          WRITE(IU06,*) ' *******************************************'
          CALL ABORT1
        ENDIF
      ELSE

!*    4.2 MODEL STARTS FROM FIELDS CREATED BY PRESET.
!         -------------------------------------------

        CDTPRO = CDATEA
      ENDIF

      IF (LWCOU) THEN
        IDELWO = KCOUSTEP
        IDELWI = KCOUSTEP
        IDELCUR= KCOUSTEP
      ELSE
        IF (NDELW_LST > 0) THEN
          DO IC=1,NDELW_LST
            IF (CDTPRO < CDTW_LST(IC)) THEN
              IDELWI=IDELWI_LST(IC)
              IDELWO=IDELWO_LST(IC)
              EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      CDATEWO = CDTPRO
      CDTBC= CDATEA
      IF (IDELT < IDELWO) CALL INCDATE(CDATEWO, IDELWO/2)
      CDAWIFL = CDTPRO
      IDELWH = MAX(IDELWI, IDELPRO)
      CALL INCDATE(CDAWIFL, IDELWH)
      CDATEFL = CDATEWO

! ----------------------------------------------------------------------

!*    5. INITIALIZE MODEL TIME VARIABLES
!        -------------------------------

!*    5.1 OUTPUT TIME VARIABLES.
!         ----------------------

      IF (FFLAG20 .OR. GFLAG20) THEN
        CDTINTT= CBPLTDT
        IF (LANAONLY) THEN
          CDTINTT=CDATEA
        ELSE
          DO WHILE (CDTINTT <= CDTPRO)
            CALL INCDATE (CDTINTT, IDELINT)
          ENDDO
        ENDIF
      ELSE
        CDTINTT = ZERO
      ENDIF

!*    5.2 FILE DISPOSE TIME AND RESTART TIME.
!         -----------------------------------

      CDTRES = CBPLTDT
      CDTBC  = CDATEA
      CALL INCDATE(CDTRES, IDELRES)

! ----------------------------------------------------------------------

!*    6. CONSISTENCY CHECK ACCORDING TO CFL CRITERION.
!        ---------------------------------------------
      IF (.NOT. LLUNSTR) THEN

        IF (IFRELFMAX > NFRE_RED ) THEN
          WRITE(IU06,*)'+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! +'
          WRITE(IU06,*)'+   IFRELFMAX SHOULD BE LESS THAN NFRE ' 
          WRITE(IU06,*)'+   IFRELFMAX = ', IFRELFMAX
          WRITE(IU06,*)'+   NFRE_RED = ', NFRE_RED
          WRITE(IU06,*)'+   ABORT SERVICE ROUTINE CALLED BY INITMDL +'
          WRITE(IU06,*)'+   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! +'
          CALL ABORT1
        ENDIF

        GVE = G/(ZPI*FR(1)*2.0_JWRB)
        DPH = DELPHI
        IF (DPH == 0.0_JWRB) THEN
          CFLP= 0.0_JWRB
        ELSE
          CFLP= DELPRO_LF*GVE/DPH
        ENDIF
        DLH = CIRC 
!       FIND THE SMALLEST DLH
        DO KX=KMSOP,KMNOP
          DLH_KX = DELLAM(KX)*COSPH(KX)
          DLH = MIN(DLH_KX,DLH)
        ENDDO

        IF (DLH == 0.0_JWRB) THEN
          CFLL= 0.0_JWRB
        ELSE
          CFLL= DELPRO_LF*GVE/DLH
        ENDIF
        IF (CFLP > 1.0_JWRB .OR. CFLL > 1.0_JWRB) THEN
          WRITE(IU06,*) ' **********************************************'
          WRITE(IU06,*) ' *                                            *'
          WRITE(IU06,*) ' *       FATAL ERROR IN SUB. INITMDL          *'
          WRITE(IU06,*) ' *       ===========================          *'
          WRITE(IU06,*) ' * CFL-CRITERION NOT FULFILLED.               *'
          WRITE(IU06,*) ' * CFLP: ',CFLP,'  GROUP VELOCITY: ',GVE
          WRITE(IU06,*) ' * CFLL: ',CFLL,'  GROUP VELOCITY: ',GVE
          WRITE(IU06,*) ' * PROPAGATION TIME DELPRO_LF : ',DELPRO_LF
          WRITE(IU06,*) ' * PROPAGATION TIME IDELPRO    : ',IDELPRO
          WRITE(IU06,*) ' * LOW FREQUENCY FREQUENCY INDEX : ',IFRELFMAX
          WRITE(IU06,*) ' * GRID Y-DISTANCE AT MOST NORTH LATITUDE: ',DPH
          WRITE(IU06,*) ' * GRID X-DISTANCE AT MOST NORTH LATITUDE: ',DLH
          WRITE(IU06,*) ' *     PROGRAM ABORTS   PROGRAM ABORTS        *'
          WRITE(IU06,*) ' *                                            *'
          WRITE(IU06,*) ' **********************************************'
          CALL ABORT1
        ENDIF

      ENDIF ! .NOT. LLUNSTR

! ----------------------------------------------------------------------

!*    7. NUMBER OF PROPAGATION TIME STEPS PER CALL.
!        ------------------------------------------

      NADV = IDELWI/IDELPRO
      NADV = MAX(NADV,1)
      IF (LANAONLY) THEN
        NADV = 0
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' **** NOTE ****'
        WRITE(IU06,*) ' NO FORWARD TIME INTEGRATION WILL BE PERFORMED'
        WRITE(IU06,*) ' ONLY THE DATA ASSIMILATION.'
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' '
      ENDIF

! ----------------------------------------------------------------------


!     INITIALISE CDTCUR
      CDTCUR=CDATECURA
      IF (.NOT.LWCOU .AND. .NOT.LRESTARTED) THEN
        IF (IREFRA == 2 .OR. IREFRA == 3) CALL INCDATE(CDTCUR, -IDELCUR)
      ENDIF

!     COMPUTE BOTTOM REFRACTION TERMS (SEE *PROPAGS_WAM*)
      IF (IREFRA /= 0) THEN
        LLUPDTTD = .TRUE.
      ELSE
        LLUPDTTD = .FALSE.
      ENDIF

!     INITIALIZE THE NEMO COUPLING
      IF (LWNEMOCOU) CALL INITNEMOCPL(BLK2LOC)


      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (OUT_METHOD == 1) THEN
          CALL INITIAL_OUTPUT_INITS
        END IF
#else
          CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      ENDIF

!     NOTE THAT CURRENT REFRACTION TERMS ARE NOW COMPUTED WHEN
!     *GETCURR* IS CALLED (SEE *PREWIND*).
!     A CALL TO PREWIND IS NEEDED TO GET THE SURFACE CURRENTS
!     BEFORE A CALL TO GETSPEC DUE TO THE TRANSFORMATION FROM
!     ABSOLUTE TO RELATIVE FRAME OF REFERENCE.

      LLINIT = .NOT.LRESTARTED
      LLINIT_FIELDG = .TRUE.

      CALL PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT,            &
     &              NXFFS, NXFFE, NYFFS, NYFFE, LLINIT_FIELDG,   &
     &              LLINIT, IREAD,                               &
     &              NFIELDS, NGPTOTG, NC, NR,                    &
     &              FIELDS, LWCUR, MASK_IN,                      &
     &              NEMO2WAM)

      WRITE(IU06,*) ' SUB. INITMDL: PREWIND DONE'
      CALL FLUSH (IU06)


!    GET SEA ICE DIMENSIONLESS ENERGY ATTENUATION COEFFICIENT
!!!! might need to restrict call when needed !!!
      CALL CIGETDEAC

!     DETERMINE THE SEA ICE REDUCTION FACTOR
      CALL CIREDUCE (WVPRPT, FF_NOW)

      WRITE(IU06,*) ' SUB. INITMDL: CIREDUCE DONE'                   

! ----------------------------------------------------------------------

!*    9.1 READ SPECTRA
!         ------------

      CALL GETSPEC(FL1, BLK2GLO, BLK2LOC, WVENVI, NBLKS, NBLKE, IREAD)

      WRITE(IU06,*) '    SUB. INITMDL: SPECTRA READ IN'
      WRITE(IU06,*) ' '
      CALL FLUSH (IU06)


!     9.2 COMPUTE FREQUENCY DEPENDENT INDICES AND COEFFICIENTS FOR SNONLIN
!         AND THE FREQUENCY FRONT TAIl REDUCTION COEFFICIENTS.
!         ------------------------------------------------------------

      CALL INISNONLIN


!*    9.3 HIGHER HARMONICS TABLES
!         -----------------------

      IF (LSECONDORDER) CALL SECONDHH_GEN

! ----------------------------------------------------------------------
!NEST
!     10. WRITE BOUNDARY VALUE FILE HEADER.
!         ------------------------------
      IF (IBOUNC == 1 .AND. .NOT. LLUNSTR) THEN
        IF (IRANK == 1) THEN
          DO II=1,GBOUNC
            IU19(II)=IWAM_GET_UNIT(IU06, CBCPREF(II), 'w', 'u', 0, 'READWRITE')
!           make the unit available for a silly fort.unit output
!           we will need to recode this a bit better !!!

            CLOSE(IU19(II))
            CALL HEADBC (IPOGBO(II)-IPOGBO(II-1),IDELPRO,TH(1),FR(1),   &
     &                   IU19(II), IU06)
          ENDDO
        ENDIF
      ENDIF
!NEST

IF (LHOOK) CALL DR_HOOK('INITMDL',1,ZHOOK_HANDLE)

END SUBROUTINE INITMDL
