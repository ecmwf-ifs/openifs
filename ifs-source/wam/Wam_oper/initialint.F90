! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE INITIALINT (IU06, NCA, NRA,                     &
     &                       NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,   &
     &                       AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC, &
     &                       ILONRGG,                                   &
     &                       LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)
! ----------------------------------------------------------------------    

!***  *INITIALINT* - INITIALIZES THE INTERPOLATION PARAMETERS TO BE USED
!                    IN INTERPOLATING FIELDS PASSED FROM THE ATMOSPHIRIC 
!                    SIDE TO THE WAVE SIDE OF THE MODEL.

!      S. ABDALLA  ECMWF  DECEMBER 2000. 

!     PURPOSE.                                                          
!     --------                                                          

!     IT COMPUTES THE WEIGHT VALUES AND THE 4 GRID POINTS TO BE USED
!     FOR INTERPOLATING THE FIELDS PASSED FROM THE ATMOSPHIRIC SIDE TO  
!     THE WAVE SIDE OF THE MODEL. "INITIALINT" MUST BE CALLED BEFORE
!     ANY CALL TO "FLDINTER".

!     LIMITATIONS.
!     ------------
!     1. IT IS ASSUMED THAT THE FIELDS ARE NOT ANGULAR (DIRECTIONAL).
!     2. IT IS ASSUMED THAT THERE IS NO MISSING DATA POINTS.
!     3. IT ONLY WORKS IF THE INPUT AND OUTPUT GRIDS ARE LATITUDE/LONGITUDE
!        (REGULAR OR IRREGULAR) OR GAUSSIAN GRIDS (FULL OR REDUCED) !!!!

!**   INTERFACE.                                                        
!     ----------                                                        

!      *CALL* *INITIALINT* (IU06, NCA, NRA,
!    &                      NGX, NGY, KRGG, KLONRGG, XDELLA, ZDELLO,
!    &                      AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
!    &                      ILONRGG,
!    &                      LLINTERPOL,DK1, DII1, DIIP1, KK, II, IIP)

!        *IU06*   - OUTPUT UNIT.
!       ATMOSPHERIC MODEL GRID:
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *ILONRGG - NUMBER OF GRID POINTS FOR EACH LATITUDE (ATM. GRID)
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
!       OUTPUT: INTERPOLATION COEFFICIENTS & INDECES
!        *IPERIODIC* - SPECIFIES IF ATM. GRID IS PERIODIC (=1) OR NOT (=0)
!        *LLINTERPOL* - FLAG (TRUE=DO INTERPOLATION, FALSE=ASSIGN ONLY)
!        *DK1M*   - COEFFICIENT
!        *DII1M*  - COEFFICIENT
!        *DIIP1M* - COEFFICIENT
!        *KKM*    - INDEX
!        *IIM*    - INDEX
!        *IIPM*   - INDEX

!     EXTERNALS.                                                        
!     ----------                                                        

!        none
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

      USE YOWGRIB   ,ONLY : IGRIB_GET_VALUE

! ----------------------------------------------------------------------

      IMPLICIT NONE 
#include "abort1.intfb.h"
#include "adjust.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, NCA, NRA
      INTEGER(KIND=JWIM), INTENT(IN) :: NGX, NGY, KRGG
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IPERIODIC
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(IN) :: KLONRGG
      INTEGER(KIND=JWIM), DIMENSION(NGY), INTENT(INOUT) :: KK
      INTEGER(KIND=JWIM), DIMENSION(NRA), INTENT(INOUT) :: ILONRGG
      INTEGER(KIND=JWIM), DIMENSION(NGX, NGY), INTENT(INOUT) :: II, IIP

      REAL(KIND=JWRB), INTENT(IN) :: XDELLA
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(IN) :: ZDELLO
      REAL(KIND=JWRB), INTENT(IN) :: AMOWEP, AMOSOP, AMOEAP, AMONOP
      REAL(KIND=JWRB), DIMENSION(NGY), INTENT(INOUT) :: DK1
      REAL(KIND=JWRB), DIMENSION(NGX, NGY), INTENT(INOUT) :: DII1, DIIP1

      LOGICAL, INTENT(INOUT) :: LLINTERPOL


      INTEGER(KIND=JWIM) :: I, J, K, JSN, IR, ISTART, ISTOP
      INTEGER(KIND=JWIM) :: KSN, KK1, KSN1, III, IIIP
      INTEGER(KIND=JWIM) :: KGRIB_HANDLE
      INTEGER(KIND=JWIM) :: NC, NR, IPLPRESENT, NB_PL, IVAL, ISCAN
      INTEGER(KIND=JWIM) :: JRGG, IREPR
      INTEGER(KIND=JWIM) :: KAMOWEP, KAMOEAP, KAMONOP, KAMOSOP
      INTEGER(KIND=JWIM) :: KRMOWEP, KRMOEAP, KRMONOP, KRMOSOP
      INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: PL

      REAL(KIND=JWRB) :: DELLA, DELLO
      REAL(KIND=JWRB) :: YFRST, YLAST
      REAL(KIND=JWRB) :: RMONOP, RMOSOP, RMOWEP, RMOEAP
      REAL(KIND=JWRB) :: XK, XLAT, XI, XII, XIIP
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION(:), ALLOCATABLE :: RDELLO

      CHARACTER(LEN=12) :: CGRIDTYPE

      LOGICAL :: LLSCANNS

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('INITIALINT',0,ZHOOK_HANDLE)

!!!   use grib handle passed down from IFS

      KGRIB_HANDLE=NGRIB_HANDLE_IFS

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Nj',NR)

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
        WRITE(IU06,*) '*  ERROR IN SUB. INITIALINT*'
        WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!! *'
        WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
        WRITE(IU06,*) '*********************************'
        CALL ABORT1
      ENDIF

      IF (JRGG == 1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'PLPresent',IPLPRESENT)
        IF (IPLPRESENT == 1) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                            &
     &                        'numberOfPointsAlongAMeridian',NB_PL)
          ALLOCATE(PL(NB_PL))
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'pl',PL)
        ELSE
          WRITE(IU06,*) '*********************************'
          WRITE(IU06,*) '*  ERROR IN SUB. INITIALINT*'
          WRITE(IU06,*) '*  NUMBER OF POINTS PER LATITUDE MISSING !!!'
          WRITE(IU06,*) '*********************************'
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
        WRITE(IU06,*)                                                   &
     &    ' SUB INITIALINT : REPRESENTATION OF THE FIELD NOT KNOWN'
        WRITE(IU06,*)'  JRGG= ',JRGG
        CALL ABORT1
      ENDIF

      IF (NC /= NCA .OR. NR /= NRA) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*    WARNING   IN SUB. INITIALINT     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NC and NCA are NOT equal:           *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')                                     &
     &             '*  NC = ', NC, '   and    NCA = ',NCA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*              OR                     *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '* NR and NRA are NOT equal:           *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,'(A,I6,A,I6,A)')                                     &
     &             '*  NR = ', NR, '   and    NRA = ',NRA,'  *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '***************************************'
      ENDIF

      IF (IREPR /= 0 .AND. IREPR /= 4) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. INITIALINT     *'
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
        WRITE(IU06,*) '*  ERROR IN SUB. INITIALINT       *'
        WRITE(IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
        WRITE(IU06,*) '*  ISCAN = ', ISCAN
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      ILONRGG(:)=0
      ALLOCATE(RDELLO(NR))

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                                &
     &                    'latitudeOfFirstGridPointInDegrees',YFRST)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                                &
     &                    'latitudeOfLastGridPointInDegrees',YLAST)

      IF (LLSCANNS) THEN
        RMONOP = YFRST 
        RMOSOP = YLAST 
      ELSE
        RMONOP = YLAST 
        RMOSOP = YFRST 
      ENDIF

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                                &
     &                    'longitudeOfFirstGridPointInDegrees',RMOWEP)


!!!   THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!!   THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!!   A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!!   RMOEAP IS IMPLICITLY KNOWN.
      IF (IREPR == 4) THEN
        DELLO = 360.0_JWRB/MAX(1,NC)
        RMOEAP = RMOWEP+360.0_JWRB - DELLO
        IPERIODIC = 1
      ELSE
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,                              &
     &                      'longitudeOfLastGridPointInDegrees',RMOEAP)

        CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
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

        DO J=1,NR-ISTART
          IF (LLSCANNS) THEN
            JSN=NR-J+1
          ELSE
            JSN=J
          ENDIF
          ILONRGG(JSN) = PL(J+ISTART) 
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
        ILONRGG(:)=NC
      ELSE
        WRITE(IU06,*)                                                   &
     &    ' SUB INITIALINT: REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF


!     FIND WHETHER INTERPOLATION IS NEEDED

      DELLA=(RMONOP-RMOSOP)/MAX(1,NR-1)

      KAMOWEP=NINT(AMOWEP*100.0_JWRB)
      KAMOEAP=NINT(AMOEAP*100.0_JWRB)
      KAMONOP=NINT(AMONOP*100.0_JWRB)
      KAMOSOP=NINT(AMOSOP*100.0_JWRB)
      KRMOWEP=NINT(RMOWEP*100.0_JWRB)
      KRMOEAP=NINT(RMOEAP*100.0_JWRB)
      KRMONOP=NINT(RMONOP*100.0_JWRB)
      KRMOSOP=NINT(RMOSOP*100.0_JWRB)

      IF (IPERIODIC == 1) THEN
        DELLO=360.0_JWRB/MAX(1,NC)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = 360.0_JWRB/MAX(1,ILONRGG(JSN)) 
        ENDDO
      ELSE
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        DO J=1,NR
          JSN=NR-J+1
          RDELLO(JSN) = (RMOEAP-RMOWEP)/MAX(1,ILONRGG(JSN)-1) 
        ENDDO
      ENDIF

      IF (KAMONOP > KRMONOP .OR. KAMONOP < KRMOSOP .OR.               &
     &    KAMOSOP < KRMOSOP .OR. KAMOSOP > KRMONOP .OR.               &
     &    (JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1            &
     &     .AND. KAMOWEP < KRMOWEP) .OR.                              &
     &    (JRGG == 0 .AND. IREPR /= 4 .AND. IPERIODIC /= 1            &
     &     .AND. KAMOEAP > NINT((RMOEAP+DELLO)*100.0_JWRB)) ) THEN

         WRITE(IU06,*) '                               '
         WRITE(IU06,*) ' SUB. INITIALINT :             '
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
         IF (JRGG == KRGG .AND. NC == NGX .AND. NR == NGY) THEN

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

        DO K=1,NGY                                                

          JSN=NGY-K+1

          XLAT = AMOSOP + (JSN-1)*XDELLA
          XLAT=MIN(MAX(XLAT,RMOSOP),RMONOP)
          XK = RMONOP - XLAT 
          XK = XK/DELLA+1.0_JWRB
          KK(K) = MAX(1,INT(XK))
          KSN = NR-KK(K)+1
          KK1 = MIN(KK(K)+1,NR)
          KSN1 = NR-KK1+1
          DK1(K)=XK-REAL(KK(K))

          DO I=1,KLONRGG(JSN)

            XI = AMOWEP+(I-1)*ZDELLO(JSN) - RMOWEP
            XI = MOD(XI+720.0_JWRB,360.0_JWRB)

            XII = XI/RDELLO(KSN)+1.0_JWRB
            III = MAX(1-IPERIODIC,INT(XII))
            II(I,K) = MIN(III,ILONRGG(KSN))
            DII1(I,K) = XII-REAL(II(I,K))

            XIIP = XI/RDELLO(KSN1)+1.0_JWRB
            IIIP = MAX(1-IPERIODIC,INT(XIIP))
            IIP(I,K) = MIN(IIIP,ILONRGG(KSN1)) 
            DIIP1(I,K)=XIIP-REAL(IIP(I,K))

          ENDDO
        ENDDO

      ELSE
        WRITE(IU06,*) ' '
        WRITE(IU06,*)                                                   &
     &   ' THE FIELDS FROM ATM. MODEL MATCH THE WAVEMODEL GRID '
        WRITE(IU06,*) ' '
      ENDIF

      DEALLOCATE (RDELLO)
      IF (ALLOCATED(PL)) DEALLOCATE(PL)

      IF (LHOOK) CALL DR_HOOK('INITIALINT',1,ZHOOK_HANDLE)

      END SUBROUTINE INITIALINT
