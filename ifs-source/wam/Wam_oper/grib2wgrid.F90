! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GRIB2WGRID (IU06, KPROMA,                                &
     &                 KGRIB_HANDLE, KGRIB, ISIZE,                  &
     &                 LLUNSTR,                                     &
     &                 NGY, KRGG, KLONRGG_LOC,                      &
     &                 NXS, NXE, NYS, NYE,                          &
     &                 XLON, YLAT,                                  &
     &                 PMISS, PPREC, PPEPS,                         &
     &                 CDATE, IFORP, IPARAM, KZLEV, KKK, MMM, FIELD)
! ----------------------------------------------------------------------    

!***  *GRIB2WGRID* - UNPACKS GRIB DATA FIELD
!                    AND PUTS IT ON THE WAVE MODEL GRID USING BILINEAR
!                    INTERPOLATION IF NECESSARY.
!      J. BIDLOT  ECMWF  FEBRUARY 2000. 
!      J. BIDLOT  ECMWF  NOVEMBER 2002  : ADD INTERPOLATION OF OCEAN
!                                         MODEL DATA.
!      J. BIDLOT  ECMWF  MARCH 2010 : MOVE TO GRIBAPI

!     PURPOSE.                                                          
!     --------                                                          

!     IT UNPACKS A GRIB DATA FIELD INCLUDING THE INFORMATION ON
!     ITS GRID AND USE IT TO PRODUCE A FIELD ON THE WAVE MODEL GRID
!     USING BILINEAR INTERPOLATION IF THE 2 GRIDS DO NOT MATCH.
!     THE CLOSEST GRID POINT IS USED INSTEAD IF FIELDS ARE
!     ARE ENCOUNTERED WITH MISSING DATA POINTS.
!     FOR NON WAVE FIELDS IF THE CLOSEST VALUE IS STILL MISSING THEN
!     THE AVERAGE OF THE NEIGHBORING NON MISSING POINTS IS TAKEN.
!     IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!     (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!     XLON and YLAT WERE INTRODUCED TO LIMIT EXTRACTION FOR POINTS THAT ARE
!     SPECIFIED AS NON MISSING WITH THOSE TWO ARRAYS.

!**   INTERFACE.
!     ----------  

!      *CALL GRIB2WGRID* (IU06, KPROMA,
!    &                    KGRIB_HANDLE, KGRIB, ISIZE,
!    &                    LLUNSTR,
!    &                    NGY, KRGG, KLONRGG_LOC,
!    &                    NXS, NXE, NYS, NYE,
!    &                    XLON, YLAT,
!    &                    PMISS, PPREC, PPEPS,
!    &                    CDATE, IFORP, IPARAM, KZLEV, KKK, MMM, FIELD)

!        *IU06*   - OUTPUT UNIT.
!        *KPROMA* - IF OPENMP IS USED THEN IT WILL BE THE
!                   NUMBER OF GRID POINTS THREADS. OTHERWISE IT CAN BE USED
!                   TO SPLIT THE GRID POINTS INTO CHUNCKS
!        *KGRIB_HANDLE GRIB API HANDLE TO THE DATA
!        *KGRIB*  - GRIB CODED DATA ARRAY
!        *ISIZE*  - SIZE OF KGRIB
!        *LLUNSTR - FLAG SPECIFYING IF THE UNSTRUCTURED GRID OPTION USED FOR THE MODEL

!      WAVE MODEL GRID SPECIFICATION (ONLY MEANINGFUL IF STRUCTURED GRID):
!        *NGY*    - TOTAL NUMBER OF LATITUDES
!        *KRGG*   - GRID DEFINITION PARAMETER (O=REGULAR, 1=IRREGULAR)
!        *KLONRGG_LOC - IF STRUCTURED GRID :: NUMBER OF GRID POINTS FOR EACH LATITUDE
!                    UNSTRUCTURED GRID :: TOTAL NUMBER OF GRID POINTS.

!        *NXS:NXE* - FIRST DIMENSION OF FIELD, XLON, YLAT USED.
!        *NYS:NYE* - SECOND DIMENSION OF FIELD, XLON, YLAT USED.
!        *XLON*   - LONGITUDE OF POINTS FOR WHICH FIELDG MUST HAVE VALUES.
!        *YLAT    - LATITUDE  OF POINTS FOR WHICH FIELDG MUST HAVE VALUES.
!!!      other points in the list (if any) are specified as missing (=PMISS).

!        *PMISS*  - VALUE FOR MISSING DATA.
!        *PPREC*  - ONLY USED FOR WAVE SPECTRAL FIELD.
!                   SMALL NUMBER USED IN SPECTRAL PACKING OF 251.
!        *PPEPS*  - ONLY USED FOR WAVE SPECTRAL FIELD.
!                   REFERENCE VALUE FOR SPECTRAL PACKING OF 251.
!        OUTPUT:
!        *CDATE*  - DATE/TIME OF THE DATA READ.         
!        *IFORP*  - FORCAST PERIOD IN SECONDS.                          
!        *IPARAM* - DATA CODE
!        *KZLEV*  - REFERENCE LEVEL IN full METER
!                   OR FOR THE OCEAN DATA WHERE IT WILL BE THE 
!                   REFERENCE DEPTH. 
!        *KKK*    - ONLY USED FOR WAVE SPECTRAL FIELD: DIRECTION NUMBER.
!        *MMM*    - ONLY USED FOR WAVE SPECTRAL FIELD: FREQUENCY NUMBER.
!        *FIELD*  - UNPACKED DATA ON WAVE MODEL GRID.

!     EXTERNALS.                                                        
!     ----------                                                        

!     *GRIBAPI*         UNPACKS MARS DATA.                               

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB   ,ONLY : IGRIB_GET_VALUE, IGRIB_SET_VALUE, JPGRIB_SUCCESS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
                                                                        
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "adjust.intfb.h"
#include "incdate.intfb.h"
#include "wstream_strg.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, KPROMA, KGRIB_HANDLE, ISIZE
      INTEGER(KIND=JWIM), INTENT(IN) :: NGY, KRGG
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      INTEGER(KIND=JWIM), INTENT(IN) :: KGRIB(ISIZE)
      INTEGER(KIND=JWIM), INTENT(IN) :: KLONRGG_LOC(NGY)
      INTEGER(KIND=JWIM), INTENT(OUT) :: IFORP, IPARAM, KZLEV, KKK, MMM

      REAL(KIND=JWRB), INTENT(IN) :: PMISS, PPREC, PPEPS
      REAL(KIND=JWRB) ,DIMENSION(NXS:NXE, NYS:NYE), INTENT(IN) :: XLON, YLAT
      REAL(KIND=JWRB) ,DIMENSION(NXS:NXE, NYS:NYE), INTENT(OUT) :: FIELD

      CHARACTER(LEN=14), INTENT(OUT) :: CDATE

      LOGICAL, INTENT(IN) :: LLUNSTR


      INTEGER(KIND=JWIM) :: LL, LS, LE
      INTEGER(KIND=JWIM) :: ITABLE
      INTEGER(KIND=JWIM) :: NC, NR, I, J, JSN, K, L, JRGG, IREPR, IR, IVAL, IDUM
      INTEGER(KIND=JWIM) :: IRET, ILEN1, IEN
      INTEGER(KIND=JWIM) :: ISCAN, ILOC, ISTAG, ICRLST, ICFG3, ICFG4, KSKIP, NGCOR
      INTEGER(KIND=JWIM) :: ITOP, IBOT, NDIM
      INTEGER(KIND=JWIM) :: NRFULL, ISTART, ISTOP
      INTEGER(KIND=JWIM) :: NUMBEROFVALUES
      INTEGER(KIND=JWIM) :: IPLPRESENT, NB_PL
      INTEGER(KIND=JWIM) :: IPERIODIC
      INTEGER(KIND=JWIM) :: KKMIN, KKMAX, ID, KK, II
      INTEGER(KIND=JWIM) :: KSN, KK1, KSN1, II1
      INTEGER(KIND=JWIM) :: KSN1LIM, IIP, IIP1, KCL, ICL, NCOUNT
      INTEGER(KIND=JWIM) :: IYYYYMMDD, IHHMM
      INTEGER(KIND=JWIM) :: IDIRSCALING, IFRESCALING
      INTEGER(KIND=JWIM) :: ILEVTYPE, ISTREAM
      INTEGER(KIND=JWIM) :: KPMONOP, KPMOEAP, KPMOSOP, KPMOWEP
      INTEGER(KIND=JWIM) :: KRMONOP, KRMOEAP, KRMOSOP, KRMOWEP 
      INTEGER(KIND=JWIM) :: KSNLIM 
      INTEGER(KIND=JWIM) :: KPARFRSTT, KPARFRSTB
      INTEGER(KIND=JWIM), ALLOCATABLE :: RLONRGG(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IGRIDCOORDINATE(:)
      INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: PL

      REAL(KIND=JWRB) :: RAD, DEG
      REAL(KIND=JWRB) :: DELLA, DELLO
      REAL(KIND=JWRB) :: PMOWEP, PMOSOP, PMOEAP, PMONOP
      REAL(KIND=JWRB) :: RMOWEP, RMOSOP, RMOEAP, RMONOP
      REAL(KIND=JWRB) :: YFRST, YLAST
      REAL(KIND=JWRB) :: FCST, STEP, START_STEP, END_STEP
      REAL(KIND=JWRB) :: XK, DK1, DK2, XI, XII, RMOWEP_KK, DII1, DII2
      REAL(KIND=JWRB) :: XIIP, DIIP1, DIIP2
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: WK(2)
      REAL(KIND=JWRB), ALLOCATABLE :: RDELLO(:), WORK(:,:,:)
      REAL(KIND=JWRB), ALLOCATABLE :: RLAT(:)
      REAL(KIND=JWRB), ALLOCATABLE :: VALUES(:)

      CHARACTER(LEN=2)  :: CDUM 
      CHARACTER(LEN=4)  :: CSTREAM
      CHARACTER(LEN=8)  :: CSTEPUNITS
      CHARACTER(LEN=8)  :: CSTEPTYPE
      CHARACTER(LEN=12) :: CGRIDTYPE

      LOGICAL :: LLINTERPOL, LLNEAREST, LLNONWAVE, LLOCEAN, LASTREAM
      LOGICAL :: LLNEAREST_LOC
      LOGICAL :: LLSCANNS
      LOGICAL :: LLDIRFLD
      LOGICAL :: LLSKIP

      DATA KPARFRSTT, KPARFRSTB /0, 0/

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('GRIB2WGRID',0,ZHOOK_HANDLE)

      CALL GSTATS(1703,0)

      CSTREAM='****'
      ISTAG=0
      RAD=4.0_JWRB*ATAN(1.0_JWRB)/180.0_JWRB
      DEG=1.0_JWRB/RAD

      FIELD(:,:) = PMISS

      PMOWEP=0.0_JWRB
      PMOSOP=0.0_JWRB
      PMOEAP=0.0_JWRB
      PMONOP=0.0_JWRB

      OUTDM: DO K = NYS, NYE
        JSN = NGY-K+1
        DO I = NXS, MIN(KLONRGG_LOC(JSN), NXE)
          IF (YLAT(I,K) == PMISS .OR. XLON(I,K) == PMISS) CYCLE
          PMOWEP=XLON(I,K)
          PMOSOP=YLAT(I,K)
          PMOEAP=XLON(I,K)
          PMONOP=YLAT(I,K)
          EXIT OUTDM
        ENDDO
      ENDDO OUTDM
      DO K = NYS, NYE
        JSN = NGY-K+1
        DO I = NXS, MIN(KLONRGG_LOC(JSN), NXE)
          IF (YLAT(I,K) == PMISS .OR. XLON(I,K) == PMISS) CYCLE
          PMOWEP=MIN(PMOWEP,XLON(I,K))
          PMOSOP=MIN(PMOSOP,YLAT(I,K))
          PMOEAP=MAX(PMOEAP,XLON(I,K))
          PMONOP=MAX(PMONOP,YLAT(I,K))
        ENDDO
      ENDDO

!*    UNPACK MARS FIELDS.                                           
!     -------------------                                           

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Nj',NRFULL)
      NR=NRFULL

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'jScansPositively',ISCAN)
      IF (ISCAN == 0) THEN
        LLSCANNS=.TRUE.
      ELSEIF (ISCAN == 1) THEN
        LLSCANNS=.FALSE.
      ELSE
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*   ERROR IN SUB. GRIB2WGRID      *'
        WRITE(IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
        WRITE(IU06,*) '*  ISCAN = ', ISCAN
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'paramId',IVAL)
      ITABLE=IVAL/1000
      IPARAM=IVAL-ITABLE*1000

!     MAKE A DISTINCTION FOR OCEAN MODEL DATA AND TEST CONFIGURATION.
!     ??? The option argument KRET does not seem to work
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'editionNumber',IEN)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'section1Length',ILEN1)

      IF ((IEN == 1 .AND. ILEN1 > 28) .OR. IEN /= 1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'localDefinitionNumber',ILOC,KRET=IRET)
        IF (IRET /= JPGRIB_SUCCESS) THEN
          WRITE(IU06,*) '   Data do not contain localDefinitionNumber !'
          WRITE(IU06,*) '   The program will continue wihtout it.'
          ILOC=-1
        ENDIF
      ELSE
        ILOC=-1
      ENDIF
      IF (ILOC == 4) THEN
        WRITE(IU06,*) '   OCEAN MODEL DATA DECODED, PARAM= ',IPARAM

        LLOCEAN=.TRUE.
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'flagForNormalOrStaggeredGrid',ISTAG)
        IF (ISTAG /= 1 .AND. ISTAG /= 0) THEN 
          WRITE(IU06,*) '***************************************'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRIB2WGRID     *'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  IT WAS ASSUMED THAT THE OCEAN GRID *' 
          WRITE(IU06,*) '*  IS NORMAL OR STAGGERED             *' 
          WRITE(IU06,*) '*  THIS IS NOT THE CASE !!!           *'
          WRITE(IU06,*) '*  ISTAG SHOULD BE 1 or 2             * '
          WRITE(IU06,*) '*  BUT IT IS ',ISTAG
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
          WRITE(IU06,*) '***************************************'
          CALL ABORT1
        ENDIF

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'flagForIrregularGridCoordinateList',ICRLST)
        IF (ICRLST /= 2) THEN
          WRITE(IU06,*) '***************************************'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRIB2WGRID     *'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  IT WAS ASSUMED THAT THE OCEAN GRID *' 
          WRITE(IU06,*) '*  COULD ONLY BE IRREGULAR IN THE Y-  *' 
          WRITE(IU06,*) '*  DIRECTION.                         *'
          WRITE(IU06,*) '*  ICRLST = ',ICRLST
          WRITE(IU06,*) '*  IT SHOULD BE = 2                   *'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
          WRITE(IU06,*) '***************************************'
          CALL ABORT1
        ENDIF

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'coordinate3Flag',ICFG3)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'coordinate4Flag',ICFG4)
        IF (ICFG3 /= 3 .OR. ICFG4 /= 4) THEN
          WRITE(IU06,*) '***************************************'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRIB2WGRID     *'
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*  IT WAS ASSUMED THAT THE OCEAN GRID *' 
          WRITE(IU06,*) '*  WAS FOR AN HORIZONTAL SECTION ONLY *' 
          WRITE(IU06,*) '*  THIS IS NOT THE CASE !!!           *' 
          WRITE(IU06,*) '*  ICFG3 SHOULD BE 3 BUT IS ',ICFG3
          WRITE(IU06,*) '*  ICFG4 SHOULD BE 4 BUT IS ',ICFG4
          WRITE(IU06,*) '*                                     *'
          WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
          WRITE(IU06,*) '***************************************'
          CALL ABORT1
        ENDIF

      ELSE
        LLOCEAN=.FALSE.
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'gridType', CGRIDTYPE)
      IF (CGRIDTYPE(1:10) == 'regular_gg') THEN
        JRGG=0
        IREPR=4
      ELSEIF (CGRIDTYPE(1:10) == 'reduced_gg') THEN
        JRGG=1
        IREPR=4
      ELSEIF (CGRIDTYPE(1:7) == 'regular') THEN
        JRGG=0
        IREPR=0
      ELSEIF (CGRIDTYPE(1:7) == 'reduced') THEN
        JRGG=1
        IREPR=0
      ELSE
        WRITE(IU06,*) '*********************************'
        WRITE(IU06,*) '*  ERROR IN SUB. GRIB2WGRID     *'
        WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!! *'
        WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
        WRITE(IU06,*) '*********************************'
        CALL ABORT1
      ENDIF

      IF (JRGG == 1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'PLPresent',IPLPRESENT)
        IF (IPLPRESENT == 1) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberOfPointsAlongAMeridian',NB_PL)
          ALLOCATE(PL(NB_PL))
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'pl',PL)
        ELSE
          WRITE(IU06,*) 'NUMBER OF POINTS PER LATITUDE MISSING !!!'
          CALL ABORT1
        ENDIF
        NC=0
        DO J=1,NB_PL
          NC = MAX(NC,PL(J))
        ENDDO
        IR=0
        DO J=1,NB_PL
          IF (PL(J) /= 0) IR=IR+1
        ENDDO
        NR=IR

      ELSEIF (JRGG == 0) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Ni',IVAL)
        NC=IVAL
      ELSE
        WRITE(IU06,*) '  GRIB2WGRID:  STRUCTURE OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'levtype',ILEVTYPE, KRET=IRET)
      IF (IRET /= JPGRIB_SUCCESS) ILEVTYPE=0

      IF (ILOC >= 0) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stream',ISTREAM)
      ELSE
        ISTREAM=0
      ENDIF
      IDUM=0
      CALL WSTREAM_STRG(ISTREAM,CSTREAM,IDUM,IDUM,CDUM,IDUM,LASTREAM)

      IF (CSTREAM == '****' .OR.                                        &
     &   (LASTREAM .AND. ILEVTYPE /= 209 .AND. ILEVTYPE /= 212 .AND.    &
     &    .NOT.LLOCEAN) ) THEN 
        LLNONWAVE=.TRUE.
      ELSE
        LLNONWAVE=.FALSE.
      ENDIF

      IF (IREPR /= 0 .AND. IREPR /= 4) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRIB2WGRID     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
        WRITE(IU06,*) '*  GRIB2WGRID CAN ONLY DEAL WITH      *'
        WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
        WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)             *' 
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF

      IF (ILEVTYPE == 105 .OR. ILEVTYPE == 160) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'level',KZLEV)
      ELSE
        KZLEV=0
      ENDIF

!*    DETERMINE INFORMATION ABOUT THE DECODED DATA 
!     --------------------------------------------

!     START DATE. 
      IF (ILOC == 11) THEN
!       Supplementary data used by the analysis:
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'dateOfAnalysis',IYYYYMMDD)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'timeOfAnalysis',IHHMM)
      ELSE
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'dataDate',IYYYYMMDD)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'time',IHHMM)
      ENDIF

      WRITE(CDATE(1:12),'(I8.8,I4.4)') IYYYYMMDD,IHHMM 
      CDATE(13:14)='00'

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stepType',CSTEPTYPE)

!     FORECAST STEP (in seconds)
      IF (CSTEPTYPE(1:7) == 'instant') THEN
        CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','s')
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'step',STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'startStep',START_STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'endStep',END_STEP)
!       THE DATA ARE VALID BETWEEN TWO TIMES. TAKE THE MIDDLE POINT
        IF (START_STEP /= END_STEP) THEN
          STEP=(END_STEP-START_STEP)/2
        ENDIF
        IFORP=STEP
      ELSE
        WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP TYPE !!!'
        WRITE(*,*) 'stepType = ',CSTEPTYPE
        CALL ABORT1
      ENDIF

!     GET DATE OF THE FORECAST INSTEAD OF STARTING DATE
      CALL INCDATE (CDATE,IFORP)


!*    DETERMINE GRID PARAMETERS.                                    
!     --------------------------                                    

      ALLOCATE(RLONRGG(NR))
      RLONRGG(:)=0
      ALLOCATE(RDELLO(NR))
      ALLOCATE(RLAT(NR))

      IF (.NOT.LLOCEAN) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'latitudeOfFirstGridPointInDegrees',YFRST)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'latitudeOfLastGridPointInDegrees',YLAST)

        IF (LLSCANNS) THEN
          RMONOP = YFRST 
          RMOSOP = YLAST 
        ELSE
          RMONOP = YLAST 
          RMOSOP = YFRST 
        ENDIF

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'longitudeOfFirstGridPointInDegrees',RMOWEP)


!!!   THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!!   THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!!   A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!!   RMOEAP IS IMPLICITLY KNOWN.
        IF (IREPR == 4 .AND. LLNONWAVE) THEN
          DELLO = 360.0_JWRB/MAX(1,NC)
          RMOEAP = RMOWEP+360.0_JWRB - DELLO
          IPERIODIC = 1
        ELSE
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'longitudeOfLastGridPointInDegrees',RMOEAP)

          CALL ADJUST (RMOWEP, RMOEAP)
          IPERIODIC = 0
          DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
          IF (RMOEAP-RMOWEP+1.5_JWRB*DELLO >= 360.0_JWRB) IPERIODIC = 1
        ENDIF

      ELSE
!       the ocean data are implicitly global 
        IPERIODIC = 1
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'iIncrement',IVAL)
        DELLO = IVAL*1.E-6_JWRB

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberInTheGridCoordinateList',NGCOR)
        ALLOCATE(IGRIDCOORDINATE(NGCOR))
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'gridCoordinate',IGRIDCOORDINATE)

        KSKIP=ISTAG
        DO KK=1,NR
          RLAT(KK)=IGRIDCOORDINATE(KK+KSKIP)*1.E-6_JWRB
        ENDDO
        DEALLOCATE(IGRIDCOORDINATE)
        IF (LLSCANNS) THEN
          RMONOP = RLAT(1) 
          RMOSOP = RLAT(NR) 
        ELSE
          RMONOP = RLAT(NR) 
          RMOSOP = RLAT(1) 
        ENDIF
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                              &
     &                      'coordinate3OfFirstGridPoint',IVAL)
        RMOWEP = IVAL*1.E-6_JWRB
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                              &
     &                      'coordinate3OfLastGridPoint',IVAL)
        RMOEAP = IVAL*1.E-6_JWRB
      ENDIF

      IF (JRGG == 1) THEN
        ISTART=1
        DO WHILE(PL(ISTART) == 0 .AND. ISTART < NB_PL)
           ISTART=ISTART+1
        ENDDO
        ISTART=ISTART-1

        ISTOP=0
        DO WHILE(PL(NB_PL-ISTOP) == 0 .AND. ISTOP < NB_PL)
          ISTOP=ISTOP+1
        ENDDO

        DO J=1,NR
          IF (LLSCANNS) THEN
            JSN=NR-J+1
          ELSE
            JSN=J
          ENDIF
          RLONRGG(JSN) = PL(J+ISTART) 

        ENDDO

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                              &
     &                      'latitudeOfFirstGridPointInDegrees',YFRST)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                              &
     &                      'latitudeOfLastGridPointInDegrees',YLAST)

        IF (ISTART /= 0 .OR. ISTOP /= 0) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                            &
     &                       'jDirectionIncrementInDegrees',DELLA)

          YFRST = YFRST-ISTART*DELLA 
          YLAST = YLAST+ISTOP*DELLA 
        ENDIF

        IF (LLSCANNS) THEN
          RMONOP = YFRST 
          RMOSOP = YLAST 
        ELSE
          RMONOP = YLAST 
          RMOSOP = YFRST 
        ENDIF


      ELSEIF (JRGG == 0) THEN
        RLONRGG(:)=NC
      ELSE
        WRITE(IU06,*) ' SUB GRIB2WGRID : REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

!     FIND WHETHER INTERPOLATION IS NEEDED

      DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

      KPMOWEP=NINT(PMOWEP*100.0_JWRB)
      KPMOEAP=NINT(PMOEAP*100.0_JWRB)
      KPMONOP=NINT(PMONOP*100.0_JWRB)
      KPMOSOP=NINT(PMOSOP*100.0_JWRB)
      KRMOWEP=NINT(RMOWEP*100.0_JWRB)
      KRMOEAP=NINT(RMOEAP*100.0_JWRB)
      KRMONOP=NINT(RMONOP*100.0_JWRB)
      KRMOSOP=NINT(RMOSOP*100.0_JWRB)

      IF (IPERIODIC == 1) THEN
        DELLO=360.0_JWRB/MAX(1,NC)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = 360.0_JWRB/MAX(1,RLONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,RLONRGG(JSN)-1) 
        ENDDO
      ENDIF


      ITOP=0
      IF (KPMONOP > KRMONOP) THEN
        IF (KPARFRSTT /= IPARAM) THEN
          KPARFRSTT=IPARAM
        ENDIF
        ITOP=INT((PMONOP-RMONOP)/DELLA)+2
      ENDIF
      IBOT=0
      IF (KPMOSOP < KRMOSOP) THEN
        IF (KPARFRSTB /= IPARAM) THEN
          KPARFRSTB=IPARAM
        ENDIF
        IBOT=INT((RMOSOP-PMOSOP)/DELLA)+2
      ENDIF

      IF (KPMONOP < KRMOSOP .OR. KPMOSOP > KRMONOP) THEN
         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  FATAL ERROR IN SUB. GRIB2WGRID *'
         WRITE(IU06,*) '*  ============================   *'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  THE MODEL DOMAIN IS OUTSIDE    *' 
         WRITE(IU06,*) '*  THE INPUT DOMAIN FOR PARAMETER *'
         WRITE(IU06,*) '  ',IPARAM 
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  MODEL:             INPPUT:     *'
         WRITE(IU06,*) '*  PMOSOP: ', PMOSOP, 'RMOSOP : ',RMOSOP
         WRITE(IU06,*) '*  PMONOP: ', PMONOP, 'RMONOP : ',RMONOP
         WRITE(IU06,*) '*  PMOWEP: ', PMOWEP, 'RMOWEP : ',RMOWEP
         WRITE(IU06,*) '*  PMOEAP: ', PMOEAP, 'RMOEAP : ',RMOEAP
         WRITE(IU06,*) '*  DELLO: ', DELLO
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*     THE PROGRAM ABORTS          *'
         WRITE(IU06,*) '***********************************'
         CALL ABORT1
      ENDIF

      IF ((JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1              &
     &     .AND. KPMOWEP < KRMOWEP) .OR.                                &
     &    (JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1              &
     &     .AND. KPMOEAP > NINT((RMOEAP+DELLO)*100.0_JWRB))) THEN

         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  WARNING IN SUB. GRIB2WGRID     *'
         WRITE(IU06,*) '*  ============================   *'
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  THE MODEL DOMAIN IS PARTIALLY OUTSIDE *' 
         WRITE(IU06,*) '*  THE INPUT DOMAIN FOR PARAMETER *'
         WRITE(IU06,*) '  ',IPARAM 
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  PMOSOP: ', PMOSOP, 'RMOSOP : ',RMOSOP
         WRITE(IU06,*) '*  PMONOP: ', PMONOP, 'RMONOP : ',RMONOP
         WRITE(IU06,*) '*  PMOWEP: ', PMOWEP, 'RMOWEP : ',RMOWEP
         WRITE(IU06,*) '*  PMOEAP: ', PMOEAP, 'RMOEAP : ',RMOEAP
         WRITE(IU06,*) '*  DELLO: ', DELLO
         WRITE(IU06,*) '*                                 *'
         WRITE(IU06,*) '*  MISSING DATA WILL BE USED FOR  *'
         WRITE(IU06,*) '*  THE PART THAT IS NOT COMMON    *'
        WRITE(IU06,*) '***********************************'
      ENDIF


      IF (LLUNSTR) THEN
!!!!!! when we will start input on unstructered grid, we will need to adapt this bit
!!!!!! for now always interpolation because the input grid is structured
        LLINTERPOL=.TRUE.

      ELSE
        LLINTERPOL=.TRUE.

        IF (KPMONOP == KRMONOP .AND. KPMOSOP == KRMOSOP .AND.           &
     &      KPMOWEP == KRMOWEP .AND. KPMOEAP == KRMOEAP       ) THEN
           IF (JRGG == KRGG .AND. NC == NXE .AND. NR == NGY) THEN

              LLINTERPOL=.FALSE.

              IF (KRGG == 1) THEN
                DO J = 1, NGY 
                  IF (RLONRGG(J) /= KLONRGG_LOC(J)) THEN
                    LLINTERPOL=.TRUE.
                    EXIT
                  ENDIF
                ENDDO
              ENDIF

           ENDIF
        ENDIF
      ENDIF

      IF (.NOT.LLSCANNS) LLINTERPOL=.TRUE.

!     GET THE DATA
      CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'missingValue',PMISS)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                                &
     &                    'numberOfEffectiveValues',NUMBEROFVALUES)
!! for reason I do not understand, I had a user who could not read grib2 wind data with
!! numberOfEffectiveValues
!! Instead, it worked with the following:
!!     &                    'getNumberOfValues',NUMBEROFVALUES)

      ALLOCATE(VALUES(NUMBEROFVALUES))
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'values',VALUES)

!     TRANSFORM WAVE SPECTRAL VALUE TO THEIR ACTUAL SCALE
      KKK=0
      MMM=0
      IF (IPARAM == 251 .AND. CSTREAM /= '****') THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'directionNumber',KKK)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'frequencyNumber',MMM)

!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(LL, LS, LE, L)
        DO LL = 1, NUMBEROFVALUES, KPROMA
          LS=LL
          LE=MIN(LS+KPROMA-1, NUMBEROFVALUES)
          DO L=LS,LE
            IF (VALUES(L) /= PMISS) THEN
              VALUES(L) = 10.0_JWRB**(VALUES(L)-ABS(PPREC))- PPEPS
            ELSE
              VALUES(L) = PMISS
            ENDIF
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

      ENDIF


      IF (.NOT.LLINTERPOL) THEN

!       REARRANGE DATA FIELD.
!       --------------------

        L = 0                                                          
        DO K = NYS, NYE                                                
          JSN = NGY-K+1
          DO I = NXS, MIN(KLONRGG_LOC(JSN), NXE)
            L = L+1                                                     
            FIELD(I,K) = VALUES(L)
          ENDDO
        ENDDO
        DEALLOCATE(VALUES)

      ELSE

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        IF ((IPARAM == 230 .OR. IPARAM == 235 .OR. IPARAM == 238 .OR.   &
     &       IPARAM == 242 .OR. IPARAM == 249 .OR. IPARAM == 122 .OR.   &
     &       IPARAM == 125 .OR. IPARAM == 128 .OR. IPARAM == 113) .AND. &
     &       CSTREAM /= '****') THEN

!         WE HAVE ASSUMED THAT THE FIELDS ARE GIVEN IN DEGREES !!!
          LLDIRFLD=.TRUE.
          NDIM=2
        ELSE
          LLDIRFLD=.FALSE.
          NDIM=1
        ENDIF

        IF (.NOT. LLNONWAVE) THEN 
!         FOR WAVE DATA ALWAYS USE THE CLOSEST GRID POINT IF
!         CORNER DATA ARE MISSING.
          LLNEAREST=.TRUE.
        ELSEIF (ITOP /= 0 .OR. IBOT /= 0) THEN
!         MISSING DATA ARE USED TO PAD THE INPUT DATA
          LLNEAREST=.TRUE.
        ELSE
!         FOR NON WAVE DATA USE THE CLOSEST GRID POINT IF
!         CORNER DATA ARE FOUND TO BE MISSING.
          L=1
          IF (L <= NUMBEROFVALUES) THEN
            DO WHILE(VALUES(L) /= PMISS)
              L=L+1
              IF (L > NUMBEROFVALUES) EXIT
            ENDDO
          ENDIF
          IF (L > NUMBEROFVALUES) THEN
            LLNEAREST=.FALSE.
          ELSE
            LLNEAREST=.TRUE.
          ENDIF
        ENDIF

        KKMIN=1-ITOP
        KKMAX=NR+IBOT
        ALLOCATE(WORK(1-IPERIODIC:NC+IPERIODIC,KKMIN:KKMAX,NDIM))

        DO ID=1,NDIM
          DO KK=KKMIN,0
            DO II=1-IPERIODIC,NC+IPERIODIC
              WORK(II,KK,ID) = PMISS
            ENDDO
          ENDDO
        ENDDO

        DO ID=1,NDIM
          DO KK=NR+1,KKMAX
            DO II=1-IPERIODIC,NC+IPERIODIC
              WORK(II,KK,ID) = PMISS
            ENDDO
          ENDDO
        ENDDO

        IF (LLSCANNS) THEN
          L = 0
          IF (LLDIRFLD) THEN
            DO KK=1,NR
              JSN=NR-KK+1
              DO II=1,RLONRGG(JSN)
                L = L+1
                IF (VALUES(L) /= PMISS) THEN
                  WORK(II,KK,1) = SIN(RAD*VALUES(L))
                  WORK(II,KK,2) = COS(RAD*VALUES(L))
                ELSE
                  WORK(II,KK,1) = PMISS 
                  WORK(II,KK,2) = PMISS
                ENDIF
              ENDDO
            ENDDO
          ELSE
            DO KK=1,NR
              JSN=NR-KK+1
              DO II=1,RLONRGG(JSN)
                L = L+1
                WORK(II,KK,1) = VALUES(L)
              ENDDO
            ENDDO
          ENDIF
        ELSE
          L = 0
          IF (LLDIRFLD) THEN
            DO KK=NR,1,-1
              JSN=KK
              DO II=1,RLONRGG(JSN)
                L = L+1
                IF (VALUES(L) /= PMISS) THEN
                  WORK(II,KK,1) = SIN(RAD*VALUES(L))
                  WORK(II,KK,2) = COS(RAD*VALUES(L))
                ELSE
                  WORK(II,KK,1) = PMISS 
                  WORK(II,KK,2) = PMISS
                ENDIF
              ENDDO
            ENDDO
          ELSE
            DO KK=NR,1,-1
              JSN=KK
              DO II=1,RLONRGG(JSN)
                L = L+1
                WORK(II,KK,1) = VALUES(L)
              ENDDO
            ENDDO
          ENDIF
        ENDIF

        DEALLOCATE(VALUES)

        IF (IPERIODIC == 1) THEN
          DO ID=1,NDIM
            DO KK=1,NR                                                
            JSN=NR-KK+1
            WORK(0,KK,ID)= WORK(RLONRGG(JSN),KK,ID)
            WORK(RLONRGG(JSN)+1,KK,ID)= WORK(1,KK,ID)
            ENDDO
          ENDDO
        ENDIF

!       loop over all wave model latitudes.
        DO K = NYS, NYE                                                
          JSN = NGY-K+1
!         loop over all wave model grid points on each latitude.
          DO I = NXS, MIN(KLONRGG_LOC(JSN), NXE)

!           *** skip all missing points !
            IF (YLAT(I,K) == PMISS .OR. XLON(I,K) == PMISS) CYCLE 

            LLSKIP=.FALSE.
            IF (.NOT.LLOCEAN) THEN
              XK = RMONOP - YLAT(I,K) 
              XK = (XK/DELLA)+0.5_JWRB*(1.0_JWRB+SIGN(1.0_JWRB,XK))
              KK = MAX(KKMIN,INT(XK))

              KSN= NR-KK+1
              KK1=MIN(KK+1,KKMAX)
              KSN1=NR-KK1+1
              DK1=ABS(XK-REAL(KK))
              DK2=1.0_JWRB-DK1

              XI = XLON(I,K) - RMOWEP

            ELSE
              KK=1
              XK = YLAT(I,K)
              DO WHILE (RLAT(KK) > XK .AND. KK <= NR) 
                KK=KK+1
              ENDDO
              KK=MAX(KK-1,1)

              KSN= NR-KK+1
              KK1=MIN(KK+1,NR)
              KSN1=NR-KK1+1
              DK1=RLAT(KK)-XK
              DK2=1.0_JWRB-DK1

              RMOWEP_KK=RMOWEP-0.5_JWRB*DELLO*ISTAG*MOD(KK,2)
              XI = XLON(I,K) - RMOWEP_KK 

            ENDIF

            XI = MOD(XI+720.0_JWRB,360.0_JWRB)

            KSNLIM=MIN(MAX(KSN,1),NR)
            XII = XI/RDELLO(KSNLIM)+1.0_JWRB
            II = INT(XII)
            IF (II < 1-IPERIODIC .OR. II > RLONRGG(KSNLIM)) LLSKIP=.TRUE.
            II = MIN(MAX(1-IPERIODIC,II),RLONRGG(KSNLIM))
            II1= II+1
            IF (II1 > RLONRGG(KSNLIM)+IPERIODIC) THEN
               IF (.NOT.LLSKIP) LLNEAREST_LOC=.TRUE.
            ENDIF
            II1= MIN(II1,RLONRGG(KSNLIM)+IPERIODIC)
            DII1=XII-REAL(II)
            DII2=1.0_JWRB-DII1

            IF (LLOCEAN) THEN
!             when the grid is staggered RMOWEP_KK varies
              RMOWEP_KK=RMOWEP-0.5_JWRB*DELLO*ISTAG*MOD(KK1,2)
              XI = XLON(I,K) - RMOWEP_KK 
            ENDIF

            KSN1LIM=MIN(MAX(KSN1,1),NR)
            XIIP = XI/RDELLO(KSN1LIM)+1.0_JWRB
            IIP = INT(XIIP)
            IF (IIP < 1-IPERIODIC .OR. IIP > RLONRGG(KSN1LIM)) LLSKIP=.TRUE.
            IIP = MIN(MAX(1-IPERIODIC,IIP),RLONRGG(KSN1LIM))
            IIP1 = IIP+1
            IF (IIP1 > RLONRGG(KSN1LIM)+IPERIODIC) THEN
               IF (.NOT.LLSKIP) LLNEAREST_LOC=.TRUE.
            ENDIF
            IIP1 = MIN(IIP1,RLONRGG(KSN1LIM)+IPERIODIC) 
            DIIP1=XIIP-REAL(IIP)
            DIIP2=1.0_JWRB-DIIP1

            IF (LLSKIP) THEN
              FIELD(I,K)=PMISS
            ELSE IF (LLNEAREST .OR. LLNEAREST_LOC) THEN
!             DETERMINE WHETHER ANY OF THE 4 CONERS HAS A MISSING DATA
!             IF SO USE THE CLOSEST GRID POINT VALUE.
              IF (WORK(II,  KK, 1) == PMISS .OR.                        &
     &            WORK(II1, KK, 1) == PMISS .OR.                        &
     &            WORK(IIP, KK1,1) == PMISS .OR.                        &
     &            WORK(IIP1,KK1,1) == PMISS    ) THEN

                IF (DK1 <= 0.5_JWRB) THEN
                  KCL=KK
                  ICL=II1
                  IF (DII1 <= 0.5_JWRB) ICL=II
                ELSE
                  KCL=KK1
                  ICL=IIP1
                  IF (DIIP1 <= 0.5_JWRB) ICL=IIP
                ENDIF
                IF (LLDIRFLD) THEN
                  IF (WORK(ICL,KCL,2) /= PMISS) THEN
                    FIELD(I,K)=DEG*ATAN2(WORK(ICL,KCL,1),WORK(ICL,KCL,2))
                  ELSE
                    FIELD(I,K)=PMISS
                  ENDIF
                ELSE
                  FIELD(I,K)=WORK(ICL,KCL,1)
                ENDIF

!               FOR NON WAVE FIELD, NON MISSING VALUE OVER SEA IS NEEDED
                IF (LLNONWAVE) THEN
                  IF (FIELD(I,K) == PMISS) THEN
                    NCOUNT=0
                    FIELD(I,K)=0.0_JWRB
                    IF (WORK(II,KK,1) /= PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(II,KK,1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF (WORK(II1,KK,1) /= PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(II1,KK,1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF (WORK(IIP,KK1,1) /= PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(IIP,KK1,1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF (WORK(IIP1,KK1,1) /= PMISS) THEN
                      FIELD(I,K)=FIELD(I,K)+WORK(IIP1,KK1,1)
                      NCOUNT=NCOUNT+1
                    ENDIF
                    IF (NCOUNT > 0) THEN
                      FIELD(I,K)=FIELD(I,K)/NCOUNT
                    ENDIF
                  ENDIF
                ENDIF

              ELSE

                DO ID=1,NDIM
                  WK(ID)=DK2*(DII2 *WORK(II,  KK, ID)+                  &
     &                        DII1 *WORK(II1, KK, ID)) +                &
     &                   DK1*(DIIP2*WORK(IIP, KK1,ID)+                  &
     &                        DIIP1*WORK(IIP1,KK1,ID))
                ENDDO
                IF (LLDIRFLD) THEN
                  FIELD(I,K)=DEG*ATAN2(WK(1),WK(2))
                ELSE
                  FIELD(I,K)=WK(1)
                ENDIF
              ENDIF

            ELSE

              DO ID=1,NDIM
                WK(ID)=DK2*(DII2 *WORK(II,  KK, ID)+                    &
     &                      DII1 *WORK(II1, KK, ID)) +                  &
     &                 DK1*(DIIP2*WORK(IIP, KK1,ID)+                    &
     &                      DIIP1*WORK(IIP1,KK1,ID))
              ENDDO
              IF (LLDIRFLD) THEN
                FIELD(I,K)=DEG*ATAN2(WK(1),WK(2))
              ELSE
                FIELD(I,K)=WK(1)
              ENDIF
            ENDIF

          ENDDO
        ENDDO

        DEALLOCATE(WORK)
      ENDIF

      DEALLOCATE(RLONRGG)
      DEALLOCATE(RDELLO)
      DEALLOCATE(RLAT)

      CALL GSTATS(1703,1)

      IF (LHOOK) CALL DR_HOOK('GRIB2WGRID',1,ZHOOK_HANDLE)


END SUBROUTINE GRIB2WGRID
