! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE INTRPOLCHK (IU06, NCA, NRA,                            &
     &                       NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,   &
     &                       AMOWEP, AMOSOP, AMOEAP, AMONOP,            &
     &                       NCAD, NRAD,                                &
     &                       RMONOP, RMOSOP, RMOWEP, RMOEAP,            &
     &                       ILONRGG, IPERIODIC, LLINTERPOL)
! ----------------------------------------------------------------------    

!***  *INTRPOLCHK* - CHECK IF FIELDS PASSED FROM THE ATMOSPHERE NEED TO BE INTERPOLATED

!        *IU06*   - OUTPUT UNIT.
!       ATMOSPHERIC MODEL GRID:
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!       WAVE MODEL GRID SPECIFICATION:
!        *NGX  *  - NUMBER OF COLUMNS IN ARRAY FIELD USED.              
!        *NGY  *  - NUMBER OF ROWS    IN ARRAY FIELD USED.              
!        *KRGG*   - GRID DEFINITION PARAMETER (0=REGULAR, 1=IRREGULAR)
!        *KLONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE
!        *XDELLA* - GRID POINT SPACING BETWEEN LATITUDES.
!        *ZDELLO* - GRID POINT SPACING PER LATITUDES. 
!        *AMOWEP* - MOST WESTERN LONGITUDE IN GRID (  1, ? ).
!        *AMOSOP* - MOST SOUTHERN LATITUDE IN GRID.( ? ,NGY).
!        *AMOEAP* - MOST EASTERN LONGITUDE IN GRID (NGX, ? ).
!        *AMONOP* - MOST NORTHERN LATITUDE IN GRID ( ? , 1 ).
!       OUTPUT: 
!        *NCAD*   - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR DECODED
!        *NRAD*   - NUMBER OF ATM. ROWS OF LATITUDES DECODED
!        *RMOWEP* - MOST WESTERN LONGITUDE IN ATM GRID
!        *RMOSOP* - MOST SOUTHERN LATITUDE IN ATM GRID
!        *RMOEAP* - MOST EASTERN LONGITUDE IN ATM GRID
!        *RMONOP* - MOST NORTHERN LATITUDE IN ATM GRID
!        *ILONRGG*- NUMBER OF GRID POINTS FOR EACH LATITUDE (ATM. GRID)
!        *IPERIODIC* - SPECIFIES IF ATM. GRID IS PERIODIC (=1) OR NOT (=0)
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)

!     EXTERNALS.                                                        
!     ----------                                                        

!        none
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS2
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE YOWGRIB   ,ONLY : IGRIB_GET_VALUE

! ----------------------------------------------------------------------

      IMPLICIT NONE 
#include "abort1.intfb.h"
#include "adjust.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NCA, NRA
      INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG
      REAL(KIND=JWRB), INTENT(IN) :: XDELLA
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO
      REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP

      INTEGER(KIND=JWIM), INTENT(OUT) :: NCAD, NRAD
      REAL(KIND=JWRB), INTENT(OUT) :: RMONOP, RMOSOP, RMOWEP, RMOEAP
      INTEGER(KIND=JWIM), DIMENSION(NRA), INTENT(OUT) :: ILONRGG
      INTEGER(KIND=JWIM), INTENT(OUT) :: IPERIODIC
      LOGICAL, INTENT(OUT) :: LLINTERPOL


      INTEGER(KIND=JWIM) :: I, J, K, JSN, IR, ISTART, ISTOP
      INTEGER(KIND=JWIM) :: KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: IPLPRESENT, NB_PL, IVAL, ISCAN
      INTEGER(KIND=JWIM) :: JRGG, IREPR
      INTEGER(KIND=JWIM) :: KAMOWEP, KAMOEAP, KAMONOP, KAMOSOP
      INTEGER(KIND=JWIM) :: KRMOWEP, KRMOEAP, KRMONOP, KRMOSOP
      INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: PL

      REAL(KIND=JWRB) :: DELLA, DELLO
      REAL(KIND=JWRB) :: YFRST, YLAST
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=12) :: CGRIDTYPE

      LOGICAL :: LLSCANNS

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INTRPOLCHK',0,ZHOOK_HANDLE)

!!!   use grib handle passed down from IFS

      KGRIB_HANDLE=NGRIB_HANDLE_IFS2

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Nj',NRAD)

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
        WRITE(IU06,*) '*  ERROR IN SUB. INTRPOLCHK*'
        WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!! *'
        WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
        WRITE(IU06,*) '*********************************'
        CALL ABORT1
      ENDIF

      IF (JRGG == 1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'PLPresent',IPLPRESENT)
        IF (IPLPRESENT == 1) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'numberOfPointsAlongAMeridian',NB_PL)  
          ALLOCATE(PL(NB_PL))
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'pl',PL)
        ELSE
          WRITE(IU06,*) '*********************************'
          WRITE(IU06,*) '*  ERROR IN SUB. INTRPOLCHK*'
          WRITE(IU06,*) '*  NUMBER OF POINTS PER LATITUDE MISSING !!!'
          WRITE(IU06,*) '*********************************'
          CALL ABORT1
        ENDIF
        NCAD=0
        DO J=1,NB_PL
          NCAD = MAX(NCAD,PL(J))
        ENDDO
        IR=0
        DO J=1,NB_PL
          IF (PL(J) /= 0) IR=IR+1
        ENDDO
        NRAD=IR

      ELSEIF (JRGG == 0) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Ni',IVAL)
        NCAD=IVAL
      ELSE
        WRITE(IU06,*)' SUB INTRPOLCHK : REPRESENTATION OF THE FIELD NOT KNOWN'
        WRITE(IU06,*)'  JRGG= ',JRGG
        CALL ABORT1
      ENDIF

      IF (NCAD /= NCA .OR. NRAD /= NRA) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*    WARNING   IN SUB. INTRPOLCHK     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NCAD and NCA are NOT equal:         *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')                                     &
     &             '*  NCAD = ', NCAD, ' and  NCA = ',NCA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*              OR                     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NRAD and NRA are NOT equal:         *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')                                     &
     &             '*  NRAD = ', NRAD, ' and  NRA = ',NRA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '***************************************'
      ENDIF

      IF (IREPR /= 0 .AND. IREPR /= 4) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INTRPOLCHK     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  UNKNOWN GRID REPRESENTATION = ',IREPR
        WRITE(IU06,*) '*  IT CAN ONLY DEAL WITH        *'
        WRITE(IU06,*) '*  LATITUDE/LONGITUDE GRID (IREPR=0)  *'
        WRITE(IU06,*) '*   OR GAUSSIAN (IREPR=4)             *' 
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*     THE PROGRAM ABORTS              *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'jScansPositively',ISCAN)
      IF (ISCAN == 0) THEN
        LLSCANNS=.TRUE.
      ELSEIF (ISCAN == 1) THEN
        LLSCANNS=.FALSE.
      ELSE
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  ERROR IN SUB. INTRPOLCHK       *'
        WRITE(IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
        WRITE(IU06,*) '*  ISCAN = ', ISCAN
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      ILONRGG(:)=0

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfFirstGridPointInDegrees',YFRST) 
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfLastGridPointInDegrees',YLAST) 

      IF (LLSCANNS) THEN
        RMONOP = YFRST 
        RMOSOP = YLAST 
      ELSE
        RMONOP = YLAST 
        RMOSOP = YFRST 
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'longitudeOfFirstGridPointInDegrees', RMOWEP) 


!!!   THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!!   THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!!   A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!!   RMOEAP IS IMPLICITLY KNOWN.
      IF (IREPR == 4) THEN
        DELLO = 360.0_JWRB/MAX(1,NCAD)
        RMOEAP = RMOWEP+360.0_JWRB - DELLO
        IPERIODIC = 1
      ELSE
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'longitudeOfLastGridPointInDegrees', RMOEAP)

        CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NCAD-1)
        IF (RMOEAP-RMOWEP+1.5_JWRB*DELLO >= 360.0_JWRB) IPERIODIC = 1
      ENDIF

      IF (JRGG == 1) THEN
        ISTART=1
        DO WHILE(PL(ISTART) == 0 .AND. ISTART > NB_PL)
          ISTART=ISTART+1
        ENDDO
        ISTART=ISTART-1

        ISTOP=0
        DO WHILE(PL(NB_PL-ISTOP) == 0 .AND. ISTOP < NB_PL)
          ISTOP=ISTOP+1
        ENDDO

        DO J=1,NRAD-ISTART
          IF (LLSCANNS) THEN
            JSN=NRAD-J+1
          ELSE
            JSN=J
          ENDIF
          ILONRGG(JSN) = PL(J+ISTART) 
        ENDDO

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfFirstGridPointInDegrees',YFRST) 
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfLastGridPointInDegrees',YLAST)

        IF (ISTART /= 0 .OR. ISTOP /= 0) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'jDirectionIncrementInDegrees',DELLA) 

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
        ILONRGG(:)=NCAD
      ELSE
        WRITE(IU06,*) ' SUB INTRPOLCHK: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF


!     FIND WHETHER INTERPOLATION IS NEEDED

      DELLA=(RMONOP-RMOSOP)/MAX(1,NRAD-1)

      KAMOWEP=NINT(AMOWEP*100.0_JWRB)
      KAMOEAP=NINT(AMOEAP*100.0_JWRB)
      KAMONOP=NINT(AMONOP*100.0_JWRB)
      KAMOSOP=NINT(AMOSOP*100.0_JWRB)
      KRMOWEP=NINT(RMOWEP*100.0_JWRB)
      KRMOEAP=NINT(RMOEAP*100.0_JWRB)
      KRMONOP=NINT(RMONOP*100.0_JWRB)
      KRMOSOP=NINT(RMOSOP*100.0_JWRB)

      IF (IPERIODIC == 1) THEN
        DELLO=360.0_JWRB/MAX(1,NCAD)
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NCAD-1)
      ENDIF

      IF (KAMONOP > KRMONOP .OR. KAMONOP < KRMOSOP .OR.               &
     &    KAMOSOP < KRMOSOP .OR. KAMOSOP > KRMONOP .OR.               &
     &    (JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1            &
     &     .AND. KAMOWEP < KRMOWEP) .OR.                              &
     &    (JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1            &
     &     .AND. KAMOEAP > NINT((RMOEAP+DELLO)*100.0_JWRB)) ) THEN

         WRITE(IU06,*) '                               '
         WRITE(IU06,*) ' SUB. INTRPOLCHK :             '
         WRITE(IU06,*) ' THE MODEL DOMAIN IS OUTSIDE   ' 
         WRITE(IU06,*) ' THE INPUT DOMAIN FOR FIELDS   '
         WRITE(IU06,*) ' AMOSOP: ', AMOSOP, 'RMOSOP : ',RMOSOP
         WRITE(IU06,*) ' AMONOP: ', AMONOP, 'RMONOP : ',RMONOP
         WRITE(IU06,*) ' AMOWEP: ', AMOWEP, 'RMOWEP : ',RMOWEP
         WRITE(IU06,*) ' AMOEAP: ', AMOEAP, 'RMOEAP : ',RMOEAP
         WRITE(IU06,*) ' DELLO: ', DELLO
         WRITE(IU06,*) ' MODEL GRID POINTS OUTSIDE WILL HAVE' 
         WRITE(IU06,*) ' THE BOUNDARY VALUES OF THE INPUT DOMAIN'
      ENDIF

      LLINTERPOL=.TRUE.

      IF (KAMONOP == KRMONOP .AND. KAMOSOP == KRMOSOP .AND.             &
     &    KAMOWEP == KRMOWEP .AND. KAMOEAP == KRMOEAP      ) THEN
         IF (JRGG == KRGG .AND. NCAD == NGX .AND. NRAD == NGY) THEN

            LLINTERPOL=.FALSE.

            IF (KRGG == 1) THEN
              DO J=1,NGY
                IF (ILONRGG(J) /= KLONRGG(J)) THEN
                  LLINTERPOL=.TRUE.
                  EXIT
                ENDIF
              ENDDO
            ENDIF

         ENDIF
      ENDIF

      IF (LLINTERPOL) THEN

!       INTERPOLATE TO WAVE MODEL GRID
!       ------------------------------

        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE FIELDS FROM ATM. MODEL'
        IF (IREPR == 0 .AND. JRGG == 0) THEN
          WRITE(IU06,*) ' ON A REGULAR LATITUDE/LONGITUDE GRID '
        ELSEIF (IREPR == 0 .AND. JRGG == 1) THEN
          WRITE(IU06,*) ' ON A REDUCED LATITUDE/LONGITUDE GRID '
        ELSEIF (IREPR == 4 .AND. JRGG == 0) THEN
          WRITE(IU06,*) ' ON A REGULAR GAUSSIAN GRID '
        ELSE
          WRITE(IU06,*) ' ON A REDUCED GAUSSIAN GRID '
        ENDIF
        WRITE(IU06,*)' ARE TO BE INTERPOLATED ONTO WAVEMODEL GRID'
        WRITE(IU06,*) ' '

      ELSE
        WRITE(IU06,*) ' '
        WRITE(IU06,*) ' THE FIELDS FROM ATM. MODEL MATCH THE WAVEMODEL GRID '
        WRITE(IU06,*) ' '
      ENDIF

      IF (ALLOCATED(PL)) DEALLOCATE(PL)

      IF (LHOOK) CALL DR_HOOK('INTRPOLCHK',1,ZHOOK_HANDLE)

      END SUBROUTINE INTRPOLCHK
