! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT,       &
 &                  NXS, NXE, NYS, NYE, LLINIT_FIELDG,      &
 &                  LLINIT, IREAD,                          &
 &                  NFIELDS, NGPTOTG, NC, NR,               &
 &                  FIELDS, LWCUR, MASK_IN,                 &
 &                  NEMO2WAM)

! ----------------------------------------------------------------------

!**** *PREWIND* - PREPARES WIND/FORCING DATA FOR WAVE MODEL.

!*    PURPOSE.
!     --------

!     EVALUATE THE FORCING AT WAVE MODEL GRID POINTS.

!**   INTERFACE.
!     ----------

!     *CALL* *PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT, 
!    &                 LLINIT, LLINIT_FIELDG,
!    &                 IREAD,
!    &                 NFIELDS, NGPTOTG, NC, NR,
!    &                 FIELDS, LWCUR, MASK_IN,
!    &                 NEMO2WAM)*

!      *BLK2LOC*             POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *WVENVI*              WAVE ENVIRONMENT.
!      *FF_NOW*    REAL      DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!      *FF_NEXT*   REAL      DATA STRUCTURE WITH THE NEXT FORCING FIELDS
!      *NXS:NXE*   INTEGER   FIRST DIMENSION OF FIELDG
!      *NYS:NYE*   INTEGER   SECOND DIMENSION OF FIELDG
!      *LLINIT_FIELDG  LOGICAL IF TRUE THEN FIELDG DATA STRUCTURE WILL BE INITIALISED
!      *LLINIT*    LOGICAL   TRUE IF WAMADSWSTAR NEEDS TO BE CALLED.   
!      *IREAD*     INTEGER   PROCESSOR WHICH WILL ACCESS THE FILE ON DISK
!                            (IF NEEDED).
!      *NFIELDS*   INTEGER   NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!      *NGPTOTG*   INTEGER   NUMBER OF ATMOSPHERIC GRID POINTS
!      *NC*        INTEGER   NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!      *NR*        INTEGER   NUMBER OF ATM. ROWS OF LATITUDES
!      *FIELDS*    REAL      ATM. FIELDS (U10, V10, AIR DENSITY, w*, U and V CURRENTS)
!      *LWCUR*     LOGICAL   INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!      *MASK_IN*   INTEGER   MASK TO INDICATE WHICH PART OF FIELDS IS RELEVANT.
!      *NEMO2WAM*            FIELDS FRON OCEAN MODEL to WAM

!     METHOD.
!     -------

!       INPUT WIND FIELDS WHICH CAN BE YOWPONENTS OF
!                USTAR, U10, USTRESS
!       ARE TRANSFORMED TO FRICTION VELOCITIES.
!       THE INPUT FIELDS HAVE TO BE ON A LAT /LONG GRID.
!       SEE SUB READWIND FOR FORMATS AND HEADER INFORMATION,
!       WHICH HAVE TO BE GIVEN TO THE PROGRAM.

!       A DOUBLE LINEAR INTERPOLATION IN SPACE IS PERFORMED
!       ONTO THE MODEL BLOCKS.
!       IF THE WIND OUTPUT TIMSTEP IS LESS THAN THE INPUT TIMESTEP
!       A LINEAR INTERPOLATION IN TIME IS PERFORMED.

!       THERE ARE TWO POSSIBILITIES WITH RESPECT TO THE WIND
!       OUTPUT FILES:

!           1. PROPAGATION TIMESTEP >= WIND INPUT STEP
!              ONE OUTPUT FILE CONTAINS IDELPRO/IDELWO WINDFIELDS
!              I.E. INFORMATION FOR ONE PROPAGATION TIMESTEP.
!              TIME FILE(I+1)= TIME FILE(I)+ IDELPRO

!           2. PROPAGATION TIMESTEP < INPUT WIND TIMESTEP
!              ONE OUTPUT FILE CONTAINS IDELWI/IDELWO WINDFIELDS
!              I.E. INFORMATION FOR ONE WIND INPUT TIMESTEP.
!              TIME FILE(I+1)= TIME FILE(I) + IDELWI

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, ENVIRONMENT, FORCING_FIELDS, OCEAN2WAVE

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOUSAMEGRID, LWNEMOCOU, LWNEMOCOURECV, IFSTOWAM
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWMAP   , ONLY : NGX      ,NGY
      USE YOWSTAT  , ONLY : CDATEA   ,CDATEE   ,IDELPRO  ,IDELWI   ,    &
     &                      IDELWO   ,LANAONLY, IDELT
      USE YOWTEST  , ONLY : IU06
      USE YOWTEXT  , ONLY : LRESTARTED
      USE YOWWIND  , ONLY : CDATEWL   ,CDAWIFL

      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "getcurr.intfb.h"
#include "getfrstwnd.intfb.h"
#include "incdate.intfb.h"
#include "init_fieldg.intfb.h"
#include "notim.intfb.h"
#include "recvnemofields.intfb.h"
#include "wamadswstar.intfb.h"

      TYPE(WVGRIDLOC), INTENT(IN)         :: BLK2LOC
      TYPE(ENVIRONMENT), INTENT(INOUT)    :: WVENVI
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FF_NEXT
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      LOGICAL, INTENT(IN) :: LLINIT_FIELDG
      LOGICAL, INTENT(IN) :: LLINIT
      INTEGER(KIND=JWIM), INTENT(IN) :: IREAD
      INTEGER(KIND=JWIM), INTENT(IN) :: NFIELDS
      INTEGER(KIND=JWIM), INTENT(IN) :: NGPTOTG
      INTEGER(KIND=JWIM), INTENT(IN) :: NC
      INTEGER(KIND=JWIM), INTENT(IN) :: NR
      REAL(KIND=JWRB),DIMENSION(NGPTOTG, NFIELDS), INTENT(IN) :: FIELDS
      LOGICAL, INTENT(IN) :: LWCUR
      INTEGER(KIND=JWIM),DIMENSION(NGPTOTG), INTENT(INOUT)  :: MASK_IN
      TYPE(OCEAN2WAVE), INTENT(INOUT) :: NEMO2WAM


      INTEGER(KIND=JWIM) :: IDELWH
      INTEGER(KIND=JWIM) :: ISTORE, IJ, ICHNK

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      TYPE(FORCING_FIELDS) :: FIELDG

      CHARACTER(LEN=14) :: CDTWIE, CDTWIS, ZERO

      LOGICAL, SAVE :: LLFRSTNEMO
      LOGICAL :: LLINIALL, LLOCAL, LLMORE
      LOGICAL :: LLNREST, LLNINIT

      DATA LLFRSTNEMO / .TRUE. /

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('PREWIND',0,ZHOOK_HANDLE)


!*    1. BEGIN AND END DATES OF WIND FIELDS TO BE PROCESSED.
!        ---------------------------------------------------

      ZERO = ' '

      IF (CDATEWL == ZERO .OR. LANAONLY) THEN
!       IF START FROM PRESET FIELDS DO FIRST FIELD IN ADDITION.
        CDTWIS = CDATEA
      ELSE
        CDTWIS = CDAWIFL
        IDELWH = -MAX(IDELPRO, IDELWI) + IDELWO
        CALL INCDATE (CDTWIS, IDELWH)
      ENDIF

      IF (CDAWIFL < CDATEE) THEN
        CDTWIE = CDAWIFL
      ELSE
        CDTWIE = CDATEE
      ENDIF
! ----------------------------------------------------------------------

!*    2. PROCESS FORCING FIELDS.
!        -----------------------

!*    2.0 GLOBAL FIELD FOR THE INPUTS
!         ---------------------------

      CALL FIELDG%ALLOC(NXS, NYS, NXE, NYE)

      IF (LLINIT_FIELDG) THEN
        LLINIALL=.TRUE.
        LLOCAL=.TRUE.
        CALL INIT_FIELDG(BLK2LOC, LLINIALL, LLOCAL,         &
     &                   NXS, NXE, NYS, NYE, FIELDG)
      ENDIF


!     2.1 IN COUPLED RUNS, TRANSFORM INPUT FORCING FIELDS TO WAM GRID.
!         -----------------------------------------------------------

      IF (LWCOU) THEN
        CALL IFSTOWAM (BLK2LOC,                    &
 &                     NFIELDS, NGPTOTG, NC, NR,   &
 &                     FIELDS, LWCUR, MASK_IN,     &
 &                     NXS, NXE, NYS, NYE, FIELDG)
      ELSE
        LWCOUSAMEGRID = .FALSE.
      ENDIF

!     2.1.1 GET DATA FROM NEMO (OR BINARY RESTART).
!           -------------------------------------
 
      IF (LWNEMOCOU .AND. LWNEMOCOURECV) THEN
        LLNREST = LLFRSTNEMO .AND. LRESTARTED
        LLNINIT = LLFRSTNEMO .AND. .NOT.LRESTARTED
        CALL RECVNEMOFIELDS(BLK2LOC, WVENVI, NEMO2WAM,  &
 &                          NXS, NXE, NYS, NYE, FIELDG, &
 &                          FF_NOW, LLNREST, LLNINIT)
        LLFRSTNEMO=.FALSE.
      ENDIF

!!!   WHEN COUPLED, IT IS ALSO NEEDED TO INITIALISE THE AIR DENSITY AND
!!!   THE CONVECTIVE VELOCITY SCALE ARRAYS SINCE THESE ARE NOT (YET) PROVIDED
!!!   AS PART OF THE GRIB RESTART FILES.
      IF (LLINIT .AND. LWCOU) THEN
        CALL WAMADSWSTAR(NXS, NXE, NYS, NYE, FIELDG, BLK2LOC, FF_NOW) 
      ENDIF

!     2.2 GET SURFACE CURRENTS TO WAM BLOCK STRUCTURE (if needed) 
!         -------------------------------------------

      CALL GETCURR(LWCUR, IREAD, BLK2LOC,          &
     &             NXS, NXE, NYS, NYE, FIELDG,     &
     &             NEMO2WAM, WVENVI)


!*    PROCESS THE OTHER FORCING FIELDS.
!     ---------------------------------
      IF (IDELWO >= IDELWI) THEN

!*      2.2 NO TIME INTERPOLATION.
!       ----------------------


        IF (CDATEWL == ZERO) THEN
!         Initialisation (either first time or following a restart)
          CALL GETFRSTWND (CDTWIS, CDTWIE,                 &
     &                     NXS, NXE, NYS, NYE, FIELDG,     &
     &                     BLK2LOC, WVENVI, FF_NOW,        &
     &                     IREAD, LWCUR, NEMO2WAM,         &
     &                     LLMORE)
        ELSE
          LLMORE = .TRUE.
        ENDIF


        IF (LLMORE) THEN
!         Update forcing
          CALL NOTIM (CDTWIS, CDTWIE,             &
     &                NXS, NXE, NYS, NYE, FIELDG, &
     &                BLK2LOC, WVENVI, FF_NEXT,   &
     &                IREAD, LWCUR, NEMO2WAM)
        ENDIF
      ELSE

!*      2.3 TIME INTERPOLATION.
!       -------------------

              WRITE (IU06,*) ' *********************************'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' * PROBLEM IN PREWIND .........  *'
              WRITE (IU06,*) ' * NOT AN OPTION ANY MORE        *'
              WRITE (IU06,*) ' *                               *'
              WRITE (IU06,*) ' *********************************'
              CALL FLUSH(IU06)
              CALL ABORT1

      ENDIF


!*    2.5 DEALLOCATE GRID ARRAYS FOR INPUT FORCING FIELDS
!         -----------------------------------------------
      CALL FIELDG%DEALLOC()

IF (LHOOK) CALL DR_HOOK('PREWIND',1,ZHOOK_HANDLE)

END SUBROUTINE PREWIND
