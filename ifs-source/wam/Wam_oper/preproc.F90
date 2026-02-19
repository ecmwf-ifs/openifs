! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

PROGRAM preproc 

! ----------------------------------------------------------------------

!**** *PROGRAM PREPROC* - PREPARE DATA (BUT NOT WINDS) FOR INPUT
!                         TO WAM WAVE MODELS.

!     SUSANNE HASSELMANN  MPI     JUNE 1986.

!     ANNEGRET SPEIDEL    MPI  OCTOBER 1988. MODFIED FOR CYCLE_2.

!     K. HUBBERT          POL     JUNE 1989  DEPTH AND CURRENT
!                                            REFRACTION.

!     H. GUNTHER   ECMWF/GKSS    APRIL 1990  LAND POINTS ARE REMOVED
!                                            FROM BLOCKS AND THE CODE
!                                            HAS BEEN RESTRUCTURED.

!     R. PORTZ     MPI         JANUARY 1991  NESTED GRID OPTION.

!     H. GUNTHER   ECMWF/GKSS    APRIL 1991  CYCLE_4 MODIFICATIONS.
!                                            MULTI-PART REMOVED.
!                                            NEW SOURCE FUNCTIONS.
!                                            LOG. DEPTH TABLE.
!     J. BIDLOT    ECMWF         SEPTEMBER 1996 REDUCED GRID.

!     J. BIDLOT    ECMWF OCTOBER 1998 : INTRODUCE USE OF MODULES
!                                       IN PLACE OF COMMON BLOCKS 
!     P.A.E.M. JANSSEN  ECMWF JULY 2010 : INTRODUCE 2ND ORDER INTERACTION
!                                         COEFFICIENTS
!
!*    PURPOSE.
!     --------

!       TO ARRANGE A GRID FOR THE WAM WAVE MODEL AND COMPUTE
!       ALL FIXED MODEL PARAMETERS WHICH ARE STORED IN DIFFERENT
!       MODULE.

!     METHOD.
!     -------

!       A REPRESENTATIVE TOPOGRAPHIC DATA SET ON LAT-LONG
!       COORDINATES CONTAINING THE MODEL SQUARE BOX REGION IS
!       READ IN.THE MODEL REGION IS EXTRACTED AND INTERPOLATED
!       ONTO GIVEN LAT-LONG GRID INCREMENTS (SEE SUB TOPOAR).
!       THE PROGRAM CHECKS FOR A PERIODIC LATITUDE GRID. IF THE
!       GRID IS NOT PERIODIC A CLOSED BASIN IS ASSUMED.
!       THE PROGRAM DOES NOT DISTINGUISH BETWEEN DEEP AND SHALLOW
!       WATER.

!       -BLOCK STRUCTURE :
!        GRID POINTS ARE COLLECTED INTO A 1-DIMENSIONAL ARRAY,
!        BLOCKS OF MAXIMALLY NIBLO ELEMENTS,  GRID POINTS
!        (ONLY SEAPOINTS) ARE COUNTED ALONG LINES OF LATITUDES
!        FROM WEST TO EAST WORKING FROM SOUTH TO NORTH.
!        BLOCKS OVERLAP OVER TWO LATITUDE LINES,TO COMPUTE NORTH
!        -SOUTH ADVECTION TERMS.

!       -NESTED GRIDS: THE GRID GENERATED CAN BE A
!         - COARSE GRID WHICH MEANS OUTPUT OF SPECTRA
!                       FOR A FOLLOW UP FINE GRID RUN.
!         - FINE   GRID WHICH MEANS INPUT OF SPECTRA
!                       FROM  AN EARLIER COARSE GRID RUN.
!         - COARSE AND FINE GRID

!       - REFRACTION: CONTROLLED BY THE REFRACTION OPTION
!         A CURRENT FIELD IS READ, INTERPOLATED TO THE MODEL
!         GRID AND STORED IN THE GRID OUTPUT FILE.

!       - PARAMETERS FOR ARRAY DIMENSIONS: THE PRORAM CHECKS
!         ALL DIMENSIONS INTERNALLY. ONLY THE BLOCK LENGTH
!         (NIBLO) IS USED FOR THE SET UP OF THE GRID, ALL
!         THE OTHER PARAMETERS HAVE TO BE LARGE ENOUGH TO
!         GET A SUCCESFULL RUN OF THE JOB. AT THE END OF
!         THE OUTPUT PROTOCOLL A LIST IS PRINTED FOR THE
!         OPTIMAL SETTINGS OF THE DIMENSION.

!**   INTERFACE.
!     ----------

!       *PROGRAM* *PREPROC*

!       *IU01*   - LOGICAL UNIT FOR INPUT OF TOPOGRAPHIC DATA.
!                  (SEE SUB TOPOAR AND MUBUF).
!       *IU02*   - LOGICAL UNIT FOR INPUT OF CURRENTS.
!       *IU03*   - LOGICAL UNIT FOR INPUT OF COARSE GRID
!                  BOUNDARY ORGANISATION (MODULE CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  FORMATED IF IFORM = 2 OTHERWISE UNFORMATED.
!                  (SEE SUB MBOUNF).
!       *IU06*   - LOGICAL UNIT FOR PRINTER OUTPUT UNIT
!       *IU07*   - LOGICAL UNIT FOR OUTPUT OF GRID ORGANISATION
!                  AND COMPUTED CONSTANTS. (UNFORMATED)
!                  (SEE SUB OUTCOM).
!       *IU08*   - LOGICAL UNITS FOR OUTPUT OF MODULE UBUF.
!                  (UNFORMATED) (SEE SUB MUBUF).
!       *IU09*   - LOGICAL UNIT FOR UNFORMATED OUTPUT OF COARSE
!                  GRID BOUNDARY ORGANISATION (MODULE CBOUND),
!                  IF THIS IS A COARSE GRID PREPROC.
!                  (SEE SUB MBOUNC).
!       *IU10*   - LOGICAl UNIT FOR UNFORMATED OUTPUT OF FINE
!                  GRID BOUNDARY ORGANISATION (MODULE CBOUND).
!                  IF THIS IS A FINE GRID PREPROC.
!                  (SEE SUB MBOUNF).
!       *IU17*   - SAME AS IU07 BUT FORMATED.
!       *IU19*   - SAME AS IU09 BUT FORMATED.
!       *IU20*   - SAME AS IU10 BUT FORMATED.

!       ALL UNITS ARE DEFINE IN SECTION 1. OF THIS PROGRAM.

!       MODULES YOWPARAM, YOWCOUPL, YOWCURR, YOWFRED, YOWINDNL, YOWGRID,
!       YOWMAP, YOWCOUT, YOWTABL, AND YOWSHAL ARE WRITTEN TO UNIT
!       IU07 AND/OR IU17.
!       ALL FREQUENCY AND DIRECTION DEPENDENT ARRAYS ARE WRITTEN FROM
!       1 TO THE USED NUMBER OF FREQUENCIES AND THE USED NUMBER OF
!       DIRECTIONS.
!       OTHER ARRAYS ARE WRITTEN ACCORDING TO THEIR DIMENSIONS.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWPARAM , ONLY : NIBLO    ,NIBLD    ,NIBLC    ,              &
     &            NANG     ,NFRE     ,NFRE_RED
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF
      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,    &
     &            TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH
      USE YOWMAP   , ONLY : NX       ,NY       ,IPER     ,IRGG     ,    &
     &            KXLTMIN  ,KXLTMAX  ,                                  &
     &            AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,    &
     &            XDELLO   ,ZDELLO   ,NLONRGG  ,LAQUA
      USE YOWSHAL  , ONLY : BATHYMAX
      USE YOWTEST  , ONLY : IU06
      USE YOWPCONS , ONLY : OLDPI    ,CIRC     ,RAD
      USE YOWUNIT  , ONLY : NPROPAGS ,IU07     ,IU08
      USE YOWUNPOOL ,ONLY : LLUNSTR  ,LPREPROC
      USE UNWAM

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "check.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "iniwcst.intfb.h"
#include "mbounc.intfb.h"
#include "mbounf.intfb.h"
#include "mcout.intfb.h"
#include "mfredir.intfb.h"
#include "mgrid.intfb.h"
#include "mtabs.intfb.h"
#include "mubuf.intfb.h"
#include "nlweigt.intfb.h"
#include "outcom.intfb.h"
#include "secondhh.intfb.h"
#include "secondhh_gen.intfb.h"
#include "topoar.intfb.h"
#include "uiprep.intfb.h"

      INTEGER(KIND=JWIM) :: IU01, IU02, IU03, IU09, IU10, IU17,IU19,IU20
      INTEGER(KIND=JWIM) :: K, IX, ICL, IFORM, LNAME,IINPC,LFILE

      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JWRB) :: OLDRAD
      REAL(KIND=JWRB) :: XLAT, XLATD, XLATMAX,PLONS, COSPHMIN
      REAL(KIND=JWRB), ALLOCATABLE :: BATHY(:,:)

      CHARACTER(LEN=1) :: C1 
      CHARACTER(LEN=80) :: FILENAME

      LOGICAL :: LLEXIST
      LOGICAL :: LLGRID

! ----------------------------------------------------------------------

      PRPLRADI=1.0_JWRB

      CALL INIWCST(PRPLRADI)

!*    1. INITIALISATION OF INPUT/OUTPUT UNITS.
!        -------------------------------------

      IU02  = 2
      IU06  = 6

! ----------------------------------------------------------------------

!*    2. USER INPUT AND LINEPRINTER PROTOCOL.
!        ------------------------------------

      CALL UIPREP (IFORM, LLGRID)

! ----------------------------------------------------------------------

!*    2.1 ALLOCATE NECESSARY ARRAYS
!     -----------------------------

      ALLOCATE(BATHY(NX, NY))

!     FR IS ALLOCATED IN UIPREP 
      ALLOCATE(DFIM(NFRE)) 
      ALLOCATE(GOM(NFRE)) 
      ALLOCATE(C(NFRE)) 
      ALLOCATE(TH(NANG)) 
      ALLOCATE(COSTH(NANG)) 
      ALLOCATE(SINTH(NANG)) 

      ALLOCATE(DELLAM(NY)) 
      ALLOCATE(SINPH(NY)) 
      ALLOCATE(COSPH(NY)) 
!     NLONRGG IS ALLOCATED IN UIPREP
      ALLOCATE(ZDELLO(NY))

      IU01 = IWAM_GET_UNIT(IU06, 'wam_topo', 'r', 'f', 0, 'READWRITE')

      FILENAME='wam_topo'
      LLEXIST=.FALSE.
      LNAME = LEN_TRIM(FILENAME)
      INQUIRE(FILE=FILENAME(1:LNAME),EXIST=LLEXIST)
      IF (.NOT. LLEXIST) THEN
        WRITE (IU06,*) '*************************************'
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO INQUIRE  *'
        WRITE (IU06,*) '*  IN PREPROC:                     *'
        WRITE (IU06,*) '*  COULD NOT FIND FILE ',FILENAME
        WRITE (IU06,*) '*                                   *'
        WRITE (IU06,*) '*************************************'
        CALL ABORT1
      ENDIF
      IU01 = IWAM_GET_UNIT(IU06, FILENAME(1:LNAME), 'r', 'f', 0, 'READWRITE')

      IF (IFORM /= 2) THEN
        IU07 = IWAM_GET_UNIT(IU06, 'wam_grid_tables', 'w', 'u', 0, 'READWRITE')
      ELSE
        IU17 = IWAM_GET_UNIT(IU06, 'wam_grid_tables_form', 'w', 'f', 0, 'READWRITE')
      ENDIF

      DO ICL=0,NPROPAGS
        WRITE(C1,'(I1)') ICL
        FILENAME='wam_subgrid_'//C1
        LFILE=0
        IF (FILENAME /= ' ') LFILE=LEN_TRIM(FILENAME)
        IU08(ICL) = IWAM_GET_UNIT(IU06,FILENAME(1:LFILE) , 'w', 'u', 0, 'READWRITE')
      ENDDO

      IF (IBOUNC == 1) THEN
!       Information of the nested grid(s) that will be produce by a coarse grid run
        IF (IFORM /= 2) THEN
          IU09=IWAM_GET_UNIT(IU06,'wam_nested_grids_info','w', 'u', 0, 'READWRITE')
        ELSE
          IU19=IWAM_GET_UNIT(IU06,'wam_nested_grids_info_form','w','f',0,'READWRITE')
        ENDIF
      ENDIF

      IF (IBOUNF == 1) THEN
!       Information of the nested grid(s) that were produced by a coarse grid run
        IF (IFORM /= 2) THEN
          IU03=IWAM_GET_UNIT(IU06,'wam_nested_grids_from_coarse_info','r', 'u', 0, 'READWRITE')
        ELSE
          IU03=IWAM_GET_UNIT(IU06,'wam_nested_grids_info_from_coarse_form','r','f',0,'READWRITE')
        ENDIF

!       Information about the boundary points that will be needed for a fine
!       grid run.
        IF (IFORM /= 2) THEN
          IU10=IWAM_GET_UNIT(IU06,'wam_boundary_grid_info', 'w', 'u', 0, 'READWRITE')
        ELSE
          IU20=IWAM_GET_UNIT(IU06,'wam_boundary_grid_info_form','w','f',0,'READWRITE')
        ENDIF
      ENDIF

! ----------------------------------------------------------------------

!*    3. INITIALISE TOTAL NUMBER OF BLOCKS,
!*       AND GRID INCREMENTS IN RADIANS AND METRES.
!        ------------------------------------------

      DELPHI = XDELLA*CIRC/360.0_JWRB
      OLDRAD=OLDPI/180.0_JWRB
      DO K=1,NY
        XLAT       = (AMOSOP + REAL(K-1)*XDELLA)*RAD
        XLATD      = (AMOSOP + REAL(K-1)*XDELLA)
        SINPH(K)   = SIN(XLAT)
        COSPH(K)   = COS(XLAT)
        IF (.NOT.LLGRID) THEN
          IF (IRGG == 1) THEN
#ifdef _CRAYFTN
              NLONRGG(K) = MAX(NINT(NX*COSD(XLATD)),2)
#else
              NLONRGG(K) = MAX(NINT(NX*COS(RAD*XLATD)),2)
#endif
            IF (MOD(NLONRGG(K),2) == 1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NX
          ENDIF      
        ENDIF      

        PLONS=(AMOEAP-AMOWEP) + IPER*XDELLO
        IF (NX == 1 .AND. NY == 1) THEN
          NLONRGG(K) = NX
          ZDELLO(K)  = XDELLO
          DELLAM(K)  = ZDELLO(K)*CIRC/360.0_JWRB
          EXIT
        ENDIF
        IF (IPER == 1) THEN
          ZDELLO(K)  = 360.0_JWRB/REAL(NLONRGG(K))
        ELSE
          ZDELLO(K)  = PLONS/REAL(NLONRGG(K)-1)
        ENDIF
        DELLAM(K)  = ZDELLO(K)*CIRC/360.0_JWRB
      ENDDO

!     IF THE POLES ARE INCLUDED, ARTIFICIALLY REMOVE THE SINGULARITY
      XLATMAX=87.5_JWRB
      COSPHMIN=COS(XLATMAX*RAD)
      DO K=1,NY
        IF (COSPH(K) <= COSPHMIN) THEN
          COSPH(K)=COS(XLATMAX*RAD)
          SINPH(K)=SIN(XLATMAX*RAD)
        ENDIF
      ENDDO 

      IF (ALLOCATED(KXLTMIN)) DEALLOCATE(KXLTMIN)
      ALLOCATE(KXLTMIN(1))
      KXLTMIN(1) = 1
      IF (ALLOCATED(KXLTMAX)) DEALLOCATE(KXLTMAX)
      ALLOCATE(KXLTMAX(1))
      KXLTMAX(1) = NY

! ----------------------------------------------------------------------

!*    4. COMPUTE GRID INDEPENDENT MODULE.
!        -------------------------------


!*    4.1 MODULE FREDIR (FREQUENCY/DIRECTION CONST).
!         ------------------------------------------

      CALL MFREDIR

!*    4.2 MODULE INDNL (WEIGHT OF NON-LINEAR INTERACTION).
!         ------------------------------------------------

      CALL NLWEIGT

!*    4.3 MODULE SHALLOW (SHALLOW WATER TABLES).
!         --------------------------------------

      CALL MTABS

!*    4.5 HIGHER HARMONICS TABLES
!         -----------------------
!
!         FIRST ROUTINE IS FOR DEEP WATER ONLY AND IS APPLIED FOR 
!         ALTIMETRIC CORRECTIONS
!         SECOND ROUTINE IS MORE GENERAL, HENCE APPLIES TO DEEP AND 
!         SHALLOW WATER
!     

      CALL SECONDHH

      CALL SECONDHH_GEN

! ----------------------------------------------------------------------

!*    5. GENERATE OUTPUT GRID INFORMATION.
!        ---------------------------------


!*    5.1 READ IN TOPOGRAPHY AND ARRANGE ON REQUESTED MODEL AREA OF
!*        REQUESTED RESOLUTION.
!         ---------------------------------------------------------

      IF (.NOT.LAQUA) THEN
        CALL TOPOAR (IU01, BATHY)
      ELSE
!       AQUA PLANET SET TO DEEP EVERYWHERE
!       EXCEPT AT THE POLES THAT ARE EXCLUDED AS LAND.
        BATHY(:,:)=BATHYMAX
        DO IX=1,NX
          BATHY(IX,1)=-999.0_JWRB
          BATHY(IX,NY)=-999.0_JWRB
        ENDDO
      ENDIF

!*    5.2 COMPUTATION OF BLOCKS.
!         ----------------------

      CALL MGRID (BATHY)

      IF (LLUNSTR .AND. LPREPROC) THEN

        CALL INIT_UNWAM

      END IF ! LLUNSTR

      NIBLD=0
      NIBLC=0

!*    5.3 COMPUTE OUTPUT POINT INDICES (MODULE YOWCOUT).
!         ----------------------------------------------

      CALL MCOUT

! ----------------------------------------------------------------------

!*    6. COMPUTE NEST INFORMATION.
!        -------------------------


!*    6.1 COMPUTE FINE GRID NEST INFORMATION (MODULE YOWFPBO).
!         ---------------------------------------------------
      IF (.NOT. LLUNSTR) THEN

        IF (IBOUNF == 1) THEN
          CALL MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
        ELSE
          IINPC  = 0
          NBOUNF = 0
        ENDIF

!*    6.2 COMPUTE COARSE GRID NEST INFORMATION (MODULE YOWCPBO).
!         -----------------------------------------------------

        IF (IBOUNC == 1) THEN
          CALL MBOUNC (IU09, IU19, IFORM)
        ELSE
          NBOUNC = 0
        ENDIF

! ----------------------------------------------------------------------

!*    7. GENERATE CURRENT FIELD.
!        -----------------------

!     THE READING OF THE INPUT CURRENTS IN NOW DONE IN INITMDL 

! ----------------------------------------------------------------------

!*    8. GENERATE AND WRITE MODULE UBUF.
!        -------------------------------

          CALL MUBUF (IU01, BATHY, IU08, NPROPAGS)
 
      END IF ! .NOT. LLUNSTR


! ----------------------------------------------------------------------

!*    9. OUTPUT OF MODULES.
!        ------------------

      CALL OUTCOM (IU07, IU17, IFORM)

! ----------------------------------------------------------------------

!*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
!*        OUTPUT OF NECESSARY DIMENSIONS.
!         --------------------------------------------------

      CALL CHECK (IINPC)
 
END PROGRAM
