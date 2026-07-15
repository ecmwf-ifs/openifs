! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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

      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWCPBO  , ONLY : IBOUNC   ,NBOUNC
      USE YOWFPBO  , ONLY : IBOUNF   ,NBOUNF
      USE YOWMAP   , ONLY : NGX      ,NGY      ,IPER     ,IRGG     ,              &
     &                      KXLTMIN  ,KXLTMAX  ,                                  &
     &                      AMOWEP   ,AMOSOP   ,AMOEAP   ,AMONOP   ,XDELLA   ,    &
     &                      XDELLO   ,ZDELLO   ,NLONRGG  ,LAQUA
      USE YOWSHAL  , ONLY : BATHYMAX
      USE YOWTEST  , ONLY : IU06
      USE YOWPCONS , ONLY : OLDPI    ,CIRC     ,RAD
      USE YOWUNIT  , ONLY : NPROPAGS ,IU07     ,IU08
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL ,ONLY : LPREPROC
      USE UNWAM     ,ONLY : INIT_UNWAM
#endif
      USE YOWABORT , ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"
#include "check.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "iniwcst.intfb.h"
#include "mbounc.intfb.h"
#include "mbounf.intfb.h"
#include "mgrid.intfb.h"
#include "mubuf.intfb.h"
#include "outcom.intfb.h"
#include "topoar.intfb.h"
#include "uiprep.intfb.h"

      INTEGER(KIND=JWIM) :: IU01, IU02, IU03, IU09, IU10, IU17,IU19,IU20
      INTEGER(KIND=JWIM) :: K, IX, ICL, IFORM, LNAME,IINPC,LFILE

      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JWRB) :: XLAT
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

      ALLOCATE(BATHY(NGX, NGY))


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

!*    3. GRID DEFINITION
!        ---------------

      ALLOCATE(ZDELLO(NGY))
      DO K=1,NGY
        XLAT = (AMOSOP + REAL(K-1)*XDELLA)*RAD
        IF (.NOT.LLGRID) THEN
          IF (IRGG == 1) THEN
            NLONRGG(K) = MAX(NINT(NGX*COS(XLAT)),2)
            IF (MOD(NLONRGG(K),2) == 1) NLONRGG(K) = NLONRGG(K)+1
          ELSE
            NLONRGG(K) = NGX
          ENDIF      
        ENDIF      

        IF (NGX == 1 .AND. NGY == 1) THEN
          NLONRGG(K) = NGX
          ZDELLO(K) = XDELLO
          EXIT
        ENDIF

        IF (IPER == 1) THEN
          ZDELLO(K) = 360.0_JWRB/REAL(NLONRGG(K),JWRB)
        ELSE
          ZDELLO(K) = (AMOEAP-AMOWEP)/REAL(NLONRGG(K)-1,JWRB)
        ENDIF

      ENDDO


      IF (ALLOCATED(KXLTMIN)) DEALLOCATE(KXLTMIN)
      ALLOCATE(KXLTMIN(1))
      KXLTMIN(1) = 1
      IF (ALLOCATED(KXLTMAX)) DEALLOCATE(KXLTMAX)
      ALLOCATE(KXLTMAX(1))
      KXLTMAX(1) = NGY

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
        DO IX=1,NGX
          BATHY(IX,1)=-999.0_JWRB
          BATHY(IX,NGY)=-999.0_JWRB
        ENDDO
      ENDIF

!*    5.2 COMPUTATION OF BLOCKS.
!         ----------------------

      CALL MGRID (BATHY)

      IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
        IF (LPREPROC) THEN
          CALL INIT_UNWAM
        ENDIF
#else
        CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
      END IF ! LLUNSTR

!*    5.3 COMPUTE OUTPUT POINT INDICES (MODULE YOWCOUT).
!         ----------------------------------------------

! ----------------------------------------------------------------------

!*    6. COMPUTE NEST INFORMATION.
!        -------------------------


      IF (.NOT. LLUNSTR) THEN
!*      6.1 COMPUTE FINE GRID NEST INFORMATION (MODULE YOWFPBO).
!           ---------------------------------------------------

        IF (IBOUNF == 1) THEN
          CALL MBOUNF (IU03, IU10, IU20, IFORM, IINPC)
        ELSE
          IINPC  = 0
          NBOUNF = 0
        ENDIF

!*      6.2 COMPUTE COARSE GRID NEST INFORMATION (MODULE YOWCPBO).
!           -----------------------------------------------------

        IF (IBOUNC == 1) THEN
          CALL MBOUNC (IU09, IU19, IFORM)
        ELSE
          NBOUNC = 0
        ENDIF

! ----------------------------------------------------------------------

!*      8. GENERATE AND WRITE MODULE UBUF.
!          -------------------------------

          CALL MUBUF (IU01, IU08, NPROPAGS)
 
      END IF ! .NOT. LLUNSTR


! ----------------------------------------------------------------------

!*    9. OUTPUT OF MODULES.
!        ------------------

      CALL OUTCOM (IU07, BATHY)

! ----------------------------------------------------------------------

!*    10. CONSISTENCY CHECK OF COMPUTED BLOCK PARAMETERS AND
!*        OUTPUT OF NECESSARY DIMENSIONS.
!         --------------------------------------------------

      CALL CHECK (IINPC)
 
END PROGRAM
