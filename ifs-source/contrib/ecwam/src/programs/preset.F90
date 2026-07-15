! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

PROGRAM preset

! ----------------------------------------------------------------------

!**** *PRESET* - GENERATES ALL BINARY FILES REQUIRED FOR A WAMODEL START
!             OR GENERATES SPECTRAL GRIB FILE REQUIRED FOR A WAMODEL
!             START. IN THAT CASE NO DRAG COEFFICIENT FIELD IS REQUIRED
!             TO RUN WAMODEL (SEE LNOCDIN IN WAMODEL).
!             OR CONVERTS AND INTERPOLATES GRIB SPECTRA PARAMETER 250 
!             INTO PARAMETER 251. THE CORRESPONDING DRAG COEFFICIENT 
!             FIELD SHOULD BE PROVIDED WITH AN APPROPRIATE MARS REQUEST.

!     SUSANNE HASSELMANN  MPI     JULY 1986.
!     ANNEGRET SPEIDEL    MPI      MAY 1988 PARAMETER STATEMENTS.
!     ANNEGRET SPEIDEL    MPI NOVEMBER 1988 CRAY-2 VERSION.
!     CYCLE_3 MODICIFATIONS:
!     ----------------------
!     RENATE PORTZ       MPI      JUNE 1990 COMPUTATION OF INITIAL
!                                           JONSWAP SPECTRA FROM
!                                           INITIAL WIND FIELD.
!     CYCLE_4 MODIFICATIONS:
!     ----------------------
!     H. GUNTHER  GKSS/ECMWF  DECEMBER 1990
!     J. BIDLOT    ECMWF   FEBRUARY 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   MARCH 1997  MODIFY ROUTINES FOR OUTPUT OF
!                          RESTART FILE 
!     B. HANSEN    ECMWF   APRIL    1997  RESTART FACILITY.
!     B. HANSEN    ECMWF   JANUARY  1998  NAMELIST INPUT.
!     B. HANSEN    ECMWF   FEBRUARY 1998  WRITE 2DSP TO FDB.
!     J. BIDLOT    ECMWF   MARCH 1998 PRODUCE GRIB SPECTRA (215) AND
!                                     INTERPOLATE TO OUTPUT GRID.
!                                     NO LAW FILE IS PRODUCED THEN.
!     J. BIDLOT    ECMWF   OCTOBER 1998   MODULES.

!*    PURPOSE.
!     --------
!       TO INITIALISE ALL FILES REQUESTED BY THE WAMODEL.
!**   INTERFACE.
!     ----------
!       *IU05*   INTEGER    USER INPUT UNIT.
!       *IU06*   INTEGER    PRINTER OUTPUT.
!       *IU07*   INTEGER    INPUT  UNIT PREPROC GRID OUTPUT.
!       *IU12*   INTEGER    OUTPUT UNIT BLOCKS OF SPECTRA.
!       *IU14*   INTEGER    OUTPUT UNIT SECOND LAT OF BLOCKS.
!       *IU15*   INTEGER    OUTPUT UNIT LAST WINDFIELDS.
!                 !OR!      OUTPUT UNIT LAST DRAG COEFFICIENT.
!     METHOD.
!     -------
!       IT IS USED TO PRODUCED A COLD START:

!       A JONSWAP SPECTRAL SHAPE IS ASSUMED.
!       JONSWAP PARAMETERS ARE DEFINED EITHER BY USER INPUT (IOPTI=0) OR
!       BY FETCH LAWS (IOPTI=1, or 2). IOPTI=3 CAN BE USED TO SPECIFY
!       INITIAL LOCALISED SWELL SYSTEMS (SEE MSWELL) 
!       THE 2-D SPECTRA ARE COMPUTED FOR
!       EACH POINT IN A BLOCK, THE WAMODEL BLOCKS ARE INITIALISED BY
!       THESE SPECTRA AND ALL BLOCKS AND OVERLAPPING LATITUDES ARE SAVED
!       IF FETCH LAWS ARE USED TO DEFINE PARAMETERS THE FIRST WIND
!       FIELD IS GENERATED OTHERWISE THE WIND FIELD IS INITIALISED
!       WITH ZEROS. THE MAIN MODEL WILL RECONSTRUCT THE WIND ANY HOW.
!       THE FILE HANDLING OF THE RESTART FILES IS COMPUTER DEPENDENT.
!       SUB GSFILE HAS TO BE MODIFIED, TO COPY THE UNIT ALIAS FILES
!       (UNITS IU12,IU14, AND IU15) TO PERMANENT FILES.


! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWFRED  , ONLY : FR       ,TH       ,IFRE1    , FR1
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWGRIBHD, ONLY : PPMISS   ,PPEPS    ,PPREC    ,NTENCODE ,    &
     &            NGRBRESS ,HOPERS   ,PPRESOL  ,LGRHDIFS ,LNEWLVTP ,    &
     &            NGRIB_VERSION
      USE YOWGRID  , ONLY : DELPHI   ,IJS      , IJL     , NTOTIJ  ,    &
     &            NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK,             & 
     &            IJSLOC   ,IJLLOC   ,IJGLOBAL_OFFSET
      USE YOWMAP   , ONLY : CLDOMAIN ,BLK2GLO   ,IRGG    ,AMOWEP   ,    &
     &            AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,XDELLO   ,    &
     &            BLK2LOC  ,NGX      ,NGY      ,NIBLO
      USE YOWNEMOFLDS , ONLY : NEMO2WAM
      USE YOWMESPAS, ONLY : LFDBIOOUT,LGRIBOUT
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG, NFRE, NFRE_RED,                       &
     &             SWAMPWIND,LL1D     ,LLUNSTR
      USE YOWPCONS , ONLY : G        ,RAD      ,DEG      ,ZMISS    ,    &
     &            ROAIR
      USE YOWSHAL  , ONLY : BATHY, WVENVI, BATHYMAX
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,CDATEA   ,    &
     &            CDATEE   ,CDATEF   ,CDTPRO   ,CDATER   ,CDATES   ,    &
     &            IDELPRO  ,IDELWI   ,IDELWO   ,                        &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,              &
     &            IREFDATE ,ISTREAM  ,NLOCGRB  ,IREFRA
      USE YOWSPEC  , ONLY : NSTART   ,NEND     ,FF_NOW   ,FL1      ,    &
     &            NBLKS    ,NBLKE
      USE YOWTABL  , ONLY :  FAC0     ,FAC1     ,FAC2     ,FAC3    ,    &
     &            FAK      ,FRHF      ,DFIMHF    , OMEGA   ,THH     ,   &
     &            DFDTH    ,IM_P      ,IM_M     ,TA       ,TB      ,    &
     &            TC_QL    ,TT_4M     ,TT_4P    ,TFAKH

      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : ICPLEN   ,USERID   ,RUNID    ,PATH     ,    &
     &            CPATH
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL ,ONLY : LPREPROC
      USE YOWPD, ONLY : MNP => npa
      USE UNSTRUCT_BOUND, ONLY : IOBPD
#endif
      USE YOWUNIT  , ONLY : IU12     ,IU14     ,IU15
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEWO  ,CDATEFL  ,    &
     &                      NXFFS    ,NXFFE    ,NYFFS    ,NYFFE    ,    &
     &                      LLNEWCURR, WSPMIN   ,FF_NEXT
      USE YOWABORT , ONLY : WAM_ABORT
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END

! -------------------------------------------------------------------

      IMPLICIT NONE 
#include "abort1.intfb.h"
#include "cigetdeac.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "iniwcst.intfb.h"
#include "init_fieldg.intfb.h"
#include "mchunk.intfb.h"
#include "mfredir.intfb.h"
#include "mstart.intfb.h"
#include "mswell.intfb.h"
#include "outspec.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "prewind.intfb.h"
#include "readmdlconf.intfb.h"
#include "savspec.intfb.h"
#include "savstress.intfb.h"

      INTEGER(KIND=JWIM), PARAMETER :: NC=1
      INTEGER(KIND=JWIM), PARAMETER :: NR=1
      INTEGER(KIND=JWIM), PARAMETER :: NGPTOTG=NC*NR
      INTEGER(KIND=JWIM), PARAMETER :: NFIELDS=1
      INTEGER(KIND=JWIM) :: ILEN, IREAD, IOPTI
      INTEGER(KIND=JWIM) :: IJ, K, M
      INTEGER(KIND=JWIM) :: IPRM, ICHNK
      INTEGER(KIND=JWIM) :: IU05, IU07 

      INTEGER(KIND=JWIM) :: I4(2)
      INTEGER(KIND=JWIM) :: MASK_IN(NGPTOTG)

      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JWRB) :: THETA, FETCH, FRMAX 
      REAL(KIND=JWRB) :: FM, ALFA, GAMMA, SA, SB, THETAQ
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: X4(2)
      REAL(KIND=JWRB) :: FIELDS(NGPTOTG,NFIELDS)
      TYPE(FORCING_FIELDS) :: FIELDG

      CHARACTER(LEN=1) :: CLTUNIT
      CHARACTER(LEN=70) :: HEADER
      CHARACTER(LEN=120) :: SFILENAME, ISFILENAME

      CHARACTER(LEN=14), PARAMETER :: ZERO='              '
      CHARACTER(LEN=14), PARAMETER :: CDUM='00000000000000'

      LOGICAL :: LLINIT
      LOGICAL :: LLINIT_FIELDG
      LOGICAL :: LWCUR
      LOGICAL :: LLINIALL, LLOCAL

#ifndef WAM_HAVE_UNWAM
      LOGICAL::LPREPROC
#endif

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ HEADER,                                         &
     &          CLDOMAIN,                                               &
     &          NANG, IFRE1, FR1, NFRE, NFRE_RED,                       &
     &          IOPTI, ITEST, ITESTB,                                   &
     &          ALFA, FM, GAMMA, SA, SB, THETA, FETCH, SWAMPWIND ,      &
     &          USERID, RUNID, PATH, CPATH,                             &
     &          CDATEA, IDELWI, CLTUNIT,                                &
     &          LLUNSTR, LPREPROC,                                      &
     &          LGRIBOUT, NGRIB_VERSION,                                &
     &          MARSTYPE, YCLASS, YEXPVER, NPROMA_WAM

!     IOPTI : IT SELECTS COLD START SPECTRAL FORM
!             (0, 1, or 2 see MSTART)
!             (3 see MSWELL)
!     ITEST, ITESTB : TEST OUTPUT LEVEL SEE MODULE YOWTEST
!     ALFA, FM, GAMMA, SA, SB, THETA, FETCH : JONSWAP PARAMETERS
!     SWAMPWIND : CONSTANT WIND SPEED USED BY SWAMP CASE.
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER 
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     CPATH : PATH FOR OUTPUT TO DISK
!     CDATEA : INPUT DATE.
!     IDELWI : INPUT WIND TIME STEP IN UNIT DEFINED BY CLTUNIT
!     CLTUNIT : INPUT WIND TIME STEP UNIT (S : seconds or H : hours)
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     MARSTYPE : DATA TYPE USED TO CODE DATA IN GRIB
!     YCLASS   : DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER : EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     NPROMA_WAM: NUMBER OF THE GRID POINTS PER LOOP WHEN CUT INTO CHUNKS


! ----------------------------------------------------------------------

      CALL MPL_INIT(KOUTPUT=1)

IF (LHOOK) CALL DR_HOOK('PRESET',0,ZHOOK_HANDLE)


      PRPLRADI=1.0_JWRB
      CALL INIWCST(PRPLRADI)

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      HEADER = ZERO
      CLDOMAIN  = 'g'
      NANG      = 0
      IFRE1     = 3
      FR1       = 4.177248E-02_JWRB
      NFRE      = 0
      NFRE_RED  = 0
      IOPTI  =    1
      ITEST  =   -9
      ITESTB =   -9
      ALFA   =    0.0_JWRB
      FM     =    0.0_JWRB
      GAMMA  =    0.0_JWRB
      SA     =    0.0_JWRB
      SB     =    0.0_JWRB
      THETA  =    0.0_JWRB
      FETCH  =    0.0_JWRB
      SWAMPWIND = 18.45_JWRB
      USERID = ZERO
      RUNID  = ZERO
      PATH   = ZERO
      CPATH  = ZERO
      CDATEA = ZERO
      CLTUNIT= 'H' 
      IDELWI =    0

      LLUNSTR  =.FALSE.
      LPREPROC =.FALSE.

      LGRIBOUT = .TRUE.
      NGRIB_VERSION = 1

      MARSTYPE = 'an'
      YCLASS   = 'rd'
      YEXPVER  = USERID//'a'

      NPROMA_WAM = 0

!*    1. DEFINE UNIT NAMES.
!        ------------------

      IU05 = 5
      IU06 = 6

      IU07 = IWAM_GET_UNIT(IU06, 'wam_grid_tables', 'r', 'u', 0, 'READWRITE')

      IU12 = 12
      IU14 = 14
      IU15 = 15

!     1.1 INITIALISATION OF MPP VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------

      LFDBIOOUT = .FALSE. 
      IRANK=1
      NPROC=1
      LL1D=.FALSE.
      X4(:)=1.0_JWRB
      NPRECR = KIND(X4)
      I4(:)=1
      NPRECI = KIND(I4)
      KTAG=100

      LLNEWCURR=.TRUE.

! ----------------------------------------------------------------------

! ALLOCATE NECESSARY ARRAYS

      ALLOCATE (NSTART(NPROC),NEND(NPROC))
      ALLOCATE (NBLKS(NPROC),NBLKE(NPROC))

! ----------------------------------------------------------------------

!*    2. READ NAMELIST NALINE.
!        ---------------------

      CDATEWO = ' '
      CDAWIFL = ' '
      CDATEFL = ' '

      READ (IU05, NALINE)

      IF( NANG <= 0 ) CALL WAM_ABORT( "Expected positive value for NANG", __FILENAME__, __LINE__ )
      IF( NFRE <= 0 ) CALL WAM_ABORT( "Expected positive value for NFRE", __FILENAME__, __LINE__ )
      IF( NFRE_RED <= 0 ) NFRE_RED = NFRE
      IF( IFRE1 <= 0 ) CALL WAM_ABORT( "Expected positive value for IFRE1",  __FILENAME__, __LINE__ )

      IF (NFRE_RED > NFRE ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE REDUCED NUMBER OF FREQUENCIES NFRE_RED *'
        WRITE (IU06,*) '* IS LARGER THAN THE TOTAL NUMNBER NFRE  !!  *'
        WRITE (IU06,*) '* NFRE_RED = ', NFRE_RED
        WRITE (IU06,*) '* NFRE     = ', NFRE
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF

!     INITIALISE THE FREQUENCY AND DIRECTION ARRAYS
      CALL MFREDIR

      IF (CLTUNIT == 'H') IDELWI = IDELWI*3600
      CDATEF = CDATEA
      CDATER = '000000000000'
      CDATES = '000000000000'
      ICPLEN=LEN_TRIM(CPATH)

      NTOTENS = 0
      NENSFNB = 0  
      ISTREAM =1045 !! if changed to an ifs stream also change LNEWLVTP
      NLOCGRB = 1
      NSYSNB  = -1
      NMETNB  = -1
      IREFDATE= 0

      LGRHDIFS=.FALSE.
      LNEWLVTP=.FALSE.

      LWCUR=.FALSE.

! ----------------------------------------------------------------------

!*    3. READ PREPROC OUTPUT.
!        --------------------

      CALL READMDLCONF (IU07)

      NINF=1
      NSUP=NIBLO

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IJS = 1
        IJL = NIBLO
        NTOTIJ = IJL-IJS+1
        IF(NPROMA_WAM == 0 ) NPROMA_WAM = NTOTIJ
        NPROMA_WAM = MIN(NPROMA_WAM, NIBLO)
        NXFFS=1
        NXFFE=MNP
        NYFFS=1
        NYFFE=1
        NSTART=1
        NEND=IJL
        IJSLOC=1
        IJLLOC=IJL
        IJGLOBAL_OFFSET=0
        NBLKS=NSTART
        NBLKE=NEND

        CALL MCHUNK

        IF (ALLOCATED(BLK2LOC%IFROMIJ)) CALL BLK2LOC%DEALLOC()
        CALL BLK2LOC%ALLOC(NPROMA_WAM,NCHNK)

        DO ICHNK = 1, NCHNK
          DO IPRM = 1, NPROMA_WAM 
            IJ = IJFROMCHNK(IPRM, ICHNK)
            IF (IJ > 0) THEN
              BLK2LOC%IFROMIJ(IPRM,ICHNK)=IJ
              BLK2LOC%KFROMIJ(IPRM,ICHNK)=1
              BLK2LOC%JFROMIJ(IPRM,ICHNK)=1
            ELSE
              BLK2LOC%IFROMIJ(IPRM,ICHNK)=BLK2LOC%IFROMIJ(1,ICHNK)
              BLK2LOC%KFROMIJ(IPRM,ICHNK)=BLK2LOC%KFROMIJ(1,ICHNK)
              BLK2LOC%JFROMIJ(IPRM,ICHNK)=BLK2LOC%JFROMIJ(1,ICHNK)
            ENDIF
          ENDDO
        ENDDO
#else
      CALL WAM_ABORT("ecwam not compiled with UNWAM support",__FILENAME__,__LINE__)
#endif

      ELSE
        NTOTIJ = NIBLO 
        IF(NPROMA_WAM == 0 ) NPROMA_WAM = NTOTIJ
        NPROMA_WAM = MIN(NPROMA_WAM, NIBLO)
        NXFFS=1
        NXFFE=NGX
        NYFFS=1
        NYFFE=NGY
        NSTART=1
        NEND=NIBLO
        IJSLOC=1
        IJLLOC=NIBLO
        IJGLOBAL_OFFSET=0
        NBLKS=NSTART
        NBLKE=NEND

        CALL MCHUNK

        IF (ALLOCATED(BLK2LOC%IFROMIJ)) CALL BLK2LOC%DEALLOC()
        CALL BLK2LOC%ALLOC(NPROMA_WAM,NCHNK)

        DO ICHNK = 1, NCHNK
          DO IPRM = 1, NPROMA_WAM 
            IJ = IJFROMCHNK(IPRM, ICHNK)
            IF (IJ > 0) THEN
              BLK2LOC%IFROMIJ(IPRM,ICHNK)=BLK2GLO%IXLG(IJ)
              BLK2LOC%KFROMIJ(IPRM,ICHNK)=BLK2GLO%KXLT(IJ)
              BLK2LOC%JFROMIJ(IPRM,ICHNK)=NGY-BLK2GLO%KXLT(IJ)+1
            ELSE
              BLK2LOC%IFROMIJ(IPRM,ICHNK)=BLK2LOC%IFROMIJ(1,ICHNK)
              BLK2LOC%KFROMIJ(IPRM,ICHNK)=BLK2LOC%KFROMIJ(1,ICHNK)
              BLK2LOC%JFROMIJ(IPRM,ICHNK)=BLK2LOC%JFROMIJ(1,ICHNK)
            ENDIF
          ENDDO
        ENDDO

      ENDIF

      IF (ALLOCATED(WVENVI%UCUR)) CALL WVENVI%DEALLOC()
      CALL WVENVI%ALLOC(NPROMA_WAM,NCHNK)

      DO ICHNK = 1, NCHNK
        DO IPRM = 1, NPROMA_WAM 
          IJ = IJFROMCHNK(IPRM,ICHNK)
          IF (IJ > 0 ) THEN
            WVENVI%DEPTH(IPRM,ICHNK) = BATHY( BLK2LOC%IFROMIJ(IPRM,ICHNK), BLK2LOC%KFROMIJ(IPRM,ICHNK) )
          ELSE
            WVENVI%DEPTH(IPRM,ICHNK) = BATHYMAX
          ENDIF
          WVENVI%UCUR(IPRM,ICHNK) = 0.0_JWRB
          WVENVI%VCUR(IPRM,ICHNK) = 0.0_JWRB
        ENDDO
      ENDDO

      DEALLOCATE(BATHY)


!!!   deallocate big arrays that were read in with READPRE

      IF (ALLOCATED(FAC0)) DEALLOCATE(FAC0)
      IF (ALLOCATED(FAC1)) DEALLOCATE(FAC1)
      IF (ALLOCATED(FAC2)) DEALLOCATE(FAC2)
      IF (ALLOCATED(FAC3)) DEALLOCATE(FAC3)
      IF (ALLOCATED(FAK)) DEALLOCATE(FAK)
      IF (ALLOCATED(FRHF)) DEALLOCATE(FRHF)
      IF (ALLOCATED(DFIMHF)) DEALLOCATE(DFIMHF)

      IF (ALLOCATED(OMEGA)) DEALLOCATE(OMEGA)
      IF (ALLOCATED(THH))   DEALLOCATE(THH)
      IF (ALLOCATED(DFDTH)) DEALLOCATE(DFDTH)
      IF (ALLOCATED(IM_P)) DEALLOCATE(IM_P)
      IF (ALLOCATED(IM_M)) DEALLOCATE(IM_M) 
      IF (ALLOCATED(TA)) DEALLOCATE(TA)
      IF (ALLOCATED(TB)) DEALLOCATE(TB)
      IF (ALLOCATED(TC_QL)) DEALLOCATE(TC_QL)
      IF (ALLOCATED(TT_4M)) DEALLOCATE(TT_4M)
      IF (ALLOCATED(TT_4P)) DEALLOCATE(TT_4P)
      IF (ALLOCATED(TFAKH)) DEALLOCATE(TFAKH)


!*    3.* SET GRIB HEADERS FOR INPUTS/OUTPUTS
!         -----------------------------------
      IF (.NOT. LGRHDIFS) THEN
!       FOR INTEGRATED PARAMETERS
        CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!       FOR SPECTRA 
        CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)
      ENDIF

! ----------------------------------------------------------------------

!*    4. PRINTER PROTOCOL OF INPUT.
!        --------------------------

      WRITE (IU06,'(A70)') HEADER
      WRITE (IU06,'('' MODEL OPTIONS  :'',/)')

      WRITE (IU06,*)'  '
      WRITE (IU06,*)' ******************************************'
      IF (LGRIBOUT) THEN
        WRITE (IU06,*)' THE OUTPUT SPECTRA WILL BE GRIBBED '
        WRITE (IU06,*)' INTO ',NFRE*NANG,' FIELDS'
      ELSE
        WRITE (IU06,*)' THE OUTPUT OF SPECTRA WILL BE BINARY  '
      ENDIF
      WRITE (IU06,*)' ******************************************'
      WRITE (IU06,*)'  '

      IF (IOPTI == 0) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' INPUT PARAMETERS.'')')
      ELSEIF (IOPTI == 1) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' WAVE ENERGY IS ZERO IN CALM WIND AREAS.'')')
      ELSEIF (IOPTI == 2) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' LOCAL WIND.'')')
        WRITE (IU06,'('' PARAMETERS USED IN CALM WIND AREAS.'')')
      ELSEIF (IOPTI == 3) THEN
        WRITE (IU06,'('' INITIAL VALUES ARE COMPUTED FROM'',            &
     &   '' IMPOSED SWELL SYSTEMS.'')')
      ELSE
        WRITE (IU06,'('' INVALID INPUT OPTION. PROGRAM WILL ABORT '')')
        CALL ABORT1
      ENDIF

      WRITE(IU06,*) '  '

      WRITE (IU06,*) ' TEST OUTPUT LEVEL IS .......... ITEST = ', ITEST

      WRITE (IU06,'('' JONSWAP PARAMETERS  :'',/)')
      WRITE (IU06,'('' ALFA : '',F10.5,'' FM : '',F10.5,'' GAMMA : '',  &
     &              F10.5,'' SA : '',F10.5,'' SB : '',F10.5)')          &
     &              ALFA, FM, GAMMA, SA, SB
      WRITE (IU06,'('' MEAN WAVE DIRECTION :  THETA = '',F10.5,         &
     &              '' DEGREE'')')  THETA
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' WIND INPUT TIMESTEP (SECONDS)      : ',IDELWI
      WRITE (IU06,*) '  '
      WRITE (IU06,*) ' END OF USER INPUT PROTOCOLL'
      WRITE (IU06,'(''  NUMBER OF DIRECTION BINS  NANG = '',I4)') NANG
      WRITE (IU06,'(''  NUMBER OF FREQUENCY BINS  NFRE = '',I4)') NFRE

! ----------------------------------------------------------------------

!*    5. PREPARE WINDFIELD.
!        ------------------

      CDATEWL = ZERO
      CDTPRO  = CDATEA
      CDATEE  = CDATEA
      CDAWIFL = CDATEA
      IDELPRO = IDELWI
      IDELWO  = IDELWI

      LWCOU=.FALSE.
      IREFRA=0
      ILEN=1

      IREAD=1

      IF (.NOT.ALLOCATED(FL1)) ALLOCATE (FL1(NPROMA_WAM,NANG,NFRE,NCHNK))

      IF (.NOT.ALLOCATED(FF_NOW%UWND)) CALL FF_NOW%ALLOC(NPROMA_WAM,NCHNK) 

      WSPMIN = 0.0_JWRB
      DO ICHNK=1,NCHNK
        FF_NOW%WSWAVE(:,ICHNK) = WSPMIN
        FF_NOW%WDWAVE(:,ICHNK) = 0.0_JWRB
        FF_NOW%USTRA(:,ICHNK) = 0.0_JWRB
        FF_NOW%VSTRA(:,ICHNK) = 0.0_JWRB
        FF_NOW%UFRIC(:,ICHNK) =  FF_NOW%WSWAVE(:,ICHNK)*0.035847_JWRB
        FF_NOW%TAUW(:,ICHNK) = 0.1_JWRB*FF_NOW%UFRIC(:,ICHNK)
        FF_NOW%TAUWDIR(:,ICHNK) = 0.0_JWRB
        FF_NOW%Z0M(:,ICHNK) = 0.00001_JWRB
        FF_NOW%CHRNCK(:,ICHNK) = 0.018_JWRB
        FF_NOW%AIRD(:,ICHNK) = ROAIR
        FF_NOW%WSTAR(:,ICHNK) = 0.0_JWRB
        FF_NOW%CICOVER(:,ICHNK) = 0.0_JWRB
        FF_NOW%CITHICK(:,ICHNK) = 0.0_JWRB
      ENDDO

      IF (.NOT.ALLOCATED(FF_NEXT%UWND)) CALL FF_NEXT%ALLOC(NPROMA_WAM,NCHNK)

      IF (.NOT.ALLOCATED(NEMO2WAM%NEMOSST)) CALL NEMO2WAM%ALLOC(NPROMA_WAM,NCHNK)

      IF (IOPTI > 0 .AND. IOPTI /= 3) THEN

!!!! might need to restict call when needed !!!
!!! remove that call in 40R3
      CALL CIGETDEAC

        LLINIT = .FALSE.
        LLINIT_FIELDG = .TRUE.

        CALL PREWIND (BLK2LOC, WVENVI, FF_NOW, FF_NEXT,            &
     &                NXFFS, NXFFE, NYFFS, NYFFE, LLINIT_FIELDG,   &
     &                LLINIT, IREAD,                               &
     &                NFIELDS, NGPTOTG, NC, NR,                    &
     &                FIELDS, LWCUR, MASK_IN,                      &
     &                NEMO2WAM)

      DO ICHNK=1,NCHNK
        FF_NOW%TAUW(:,ICHNK) = 0.1_JWRB * FF_NOW%UFRIC(:,ICHNK)**2
      ENDDO

      ENDIF

! ----------------------------------------------------------------------

!*    6. DEFINE FETCH AND MAXIMUM PEAK FREQUENCY.
!        ----------------------------------------


      IF (FETCH < 0.1E-5_JWRB) FETCH = 0.5_JWRB*DELPHI
      FRMAX = FM
      IF (IOPTI /= 0 .AND. IOPTI /= 3) THEN
        WRITE (IU06,*) ' FETCH USED (METRES)       : ', FETCH
        WRITE (IU06,*) ' MAXIMUM PEAK FREQUENCY IS : ', FRMAX
      ENDIF

! ----------------------------------------------------------------------

!*    7. GENERATE AND WRITE START FILES.
!        -------------------------------

      IF (IOPTI /= 3) THEN
        THETAQ = THETA * RAD
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          CALL MSTART (IOPTI, FETCH, FRMAX, THETAQ,                      &
     &                 FM, ALFA, GAMMA, SA, SB,                          &
     &                 1, NPROMA_WAM, FL1(:,:,:,ICHNK),                  &
     &                 FF_NOW%WSWAVE(:,ICHNK), FF_NOW%WDWAVE(:,ICHNK))
        ENDDO
!$OMP END PARALLEL DO

        CDTPRO  = ZERO
        CDATEWO = ZERO
        CDAWIFL = ZERO
        CDATEFL = ZERO

      ELSE

        LLINIALL=.FALSE.
        LLOCAL=.TRUE.
        CALL FIELDG%ALLOC(NXFFS, NYFFS, NXFFE, NYFFE)

        CALL INIT_FIELDG(BLK2LOC, LLINIALL, LLOCAL, &
     &                   NXFFS, NXFFE, NYFFS, NYFFE, FIELDG)

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ICHNK)
        DO ICHNK = 1, NCHNK
          CALL MSWELL (1, NPROMA_WAM, BLK2LOC%IFROMIJ(:,ICHNK), BLK2LOC%JFROMIJ(:,ICHNK), NXFFS, NXFFE, NYFFS, &
     &                 NYFFE, FIELDG, FL1(:,:,:,ICHNK) )
        ENDDO
!$OMP END PARALLEL DO

        CALL FIELDG%DEALLOC()


        IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
!         reset points with no flux out of the boundary to 0
          DO ICHNK = 1, NCHNK
            DO M=1,NFRE
              DO K=1,NANG
                DO IPRM = 1, KIJL4CHNK(ICHNK)
                  IJ = IJFROMCHNK(IPRM, ICHNK)
                  FL1(IPRM, K, M, ICHNK) = FL1(IPRM, K, M, ICHNK) * IOBPD(K,IJ)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
#endif
        ENDIF

      ENDIF

! ----------------------------------------------------------------------

!*    8. DISPOSE START FILES.
!        --------------------

!     STRESS RELATED FIELDS :
      IF (.NOT.LGRIBOUT) THEN
        CALL SAVSTRESS(WVENVI, FF_NOW, NBLKS, NBLKE, CDATEA, CDATEA) 
      ENDIF


!     SPECTRA :
      WRITE (IU06,*) ' SAVING SPECTRA '
      CALL FLUSH(IU06)
      IF (LGRIBOUT) THEN
!       THE COLD START SPECTRA WILL BE SAVED AS GRIB FILES.
        CDTPRO  = CDATEA
        CALL OUTSPEC(FL1, FF_NOW)

      ELSE
        CALL SAVSPEC(FL1, NBLKS, NBLKE, CDATEA, CDATEA, CDUM)
      ENDIF

! ----------------------------------------------------------------------

!*    9. END OF JOB: DELETE WORK FILES.
!        ------------------------------

      WRITE (IU06,*) ' '
      WRITE (IU06,*) ' PROGRAM PRESET: ALL DONE'

      CALL MPL_END()

IF (LHOOK) CALL DR_HOOK('PRESET',1,ZHOOK_HANDLE)

END PROGRAM preset
