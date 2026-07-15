! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MICEP (IPARAM, KIJS, KIJL, IFROMIJ, JFROMIJ,    &
 &                NXS, NXE, NYS, NYE, FIELDG,              &
 &                CICVR, CITH, NEMOCICOVER, NEMOCITHICK)

!-------------------------------------------------------------------

!**** *MICEP* - CLEAN UP SEA FRACTION FOR ALL SEA POINTS.
!               DETERNINE THE SEA ICE THICKNESS IF NOT SUPPLIED
!               AND CHECk THE CONSISTENCY BETWEEN SEA ICE COVER AND
!               THICKNESS. IMPOSE MINIMUM THICKNESS.

!     R. PORTZ     MPI HAMBURG   OCTOBER 1992
!     J. BIDLOT    ECMWF         JUNE 1996    MESSAGE PASSING
!     J. BIDLOT    ECMWF         JANUARY 1998 : CORRECT THRESHOLD FOR
!                                               ICE TO BE LT 271.50
!     J. BIDLOT    ECMWF         FEBRUARY 2000 : USE OF SEA ICE FRACTION
!     J. BIDLOT    ECMWF         AUG 2006 :  


!     PURPOSE.
!     --------
!            SELECTS POINTS (TEMP < 271.5 KELVIN) IN THE TEMPERATURE
!            GRID AND DEFINE SEA COVER AS 1., OTHERWISE 0. 
!            BUT IF SEA ICE FRACTION IS GIVEN DEFINE MODEL SEA ICE
!            COVER MAKING SURE IT IS ALWAYS >=0. AND <=1. 

!**   INTERFACE
!     ---------
!             *CALL MICEP* *(IPARAM, KIJS, KIJL, IFROMIJ, JFROMIJ,
!                            NXS, NXE, NYS, NYE, FIELDG,
!                            CICVR, CITH, NEMOCICOVER, NEMOCITHICK)

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELDG*CICOVER
!      *KIJS*      INDEX OF FIRST GRIDPOINT
!      *KIJL*      INDEX OF LAST GRIDPOINT
!      *IFROMIJ*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *JFROMIJ*   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!      *NXS:NXE*   FIRST DIMENSION OF FIELDG
!      *NYS:NYE*   SECOND DIMENSION OF FIELDG
!      *FIELDG*    INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!      *CICVR*     SEA ICE COVER.
!      *CITH*      SEA ICE THICKNESS.
!      *NEMOCICOVER NEMO SEA ICE COVER (if used)
!      *NEMOCITHICK NEMO SEA ICE THICKNESS (if used)



!     METHOD.
!     -------
!             NONE

!    EXTERNAL.
!    ---------
!             *NONE* 

!---------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU, JWRO
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWICE   , ONLY : CITHRSH  ,LICERUN ,LMASKICE   ,LICETH     , &
     &               HICMIN
      USE YOWMAP   , ONLY : NGX      ,NGY     ,CLDOMAIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : SWAMPCITH
      USE YOWPCONS , ONLY : ZMISS
      USE YOWTEST  , ONLY : IU06
      USE YOWCOUP  , ONLY : LWCOU    ,LWNEMOCOUCIC, LWNEMOCOUCIT
      USE YOWNEMOFLDS, ONLY : LNEMOICEREST

      USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IPARAM, KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(INOUT) :: CICVR, CITH
      REAL(KIND=JWRO), DIMENSION (KIJS:KIJL), INTENT(IN) :: NEMOCICOVER, NEMOCITHICK


      INTEGER(KIND=JWIM) :: IJ, IX, IY
      INTEGER(KIND=JWIM) :: NICE

!     CONSTANTS FOR PARAMETRISATION OF SEA ICE THICKNESS:
      REAL(KIND=JWRB), PARAMETER :: C1=0.2_JWRB
      REAL(KIND=JWRB), PARAMETER :: C2=0.4_JWRB

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

!     1. INITIALIZATION
!     ------------------

      IF (LHOOK) CALL DR_HOOK('MICEP',0,ZHOOK_HANDLE)

!*    2. WE DEFINE ICE FOR THOSE POINTS WHERE SST < 271.5 KELVIN). 
!        OR CLEAN SEA ICE FRACTION TO DEFINE MODEL SEA ICE COVER.
!        -------------------------------------------------------

      IF (LWNEMOCOUCIC) THEN
        IF (LWCOU) THEN
          DO IJ=KIJS,KIJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            IF (FIELDG%LKFR(IX,IY) <= 0.0_JWRB ) THEN
!            if lake cover = 0, we assume open ocean point, then get sea ice directly from NEMO 
              CICVR(IJ) = NEMOCICOVER(IJ)
            ELSE
!            get ice information from atmopsheric model
              IF (FIELDG%CICOVER(IX,IY) == ZMISS .OR.                   &
     &            FIELDG%CICOVER(IX,IY) < 0.01_JWRB .OR.                &
     &            FIELDG%CICOVER(IX,IY) > 1.01_JWRB ) THEN
                CICVR(IJ) = 0.0_JWRB
              ELSEIF (FIELDG%CICOVER(IX,IY) > 0.95_JWRB) THEN
                CICVR(IJ) = 1.0_JWRB
              ELSE
                CICVR(IJ) = FIELDG%CICOVER(IX,IY)
              ENDIF
            ENDIF
          ENDDO
        ELSE
          DO IJ=KIJS,KIJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            CICVR(IJ) = NEMOCICOVER(IJ)
          ENDDO
        ENDIF
      ELSEIF (IPARAM == 31) THEN
        DO IJ=KIJS,KIJL
          IX = IFROMIJ(IJ)
          IY = JFROMIJ(IJ)
          IF (FIELDG%CICOVER(IX,IY) == ZMISS .OR.                       &
     &        FIELDG%CICOVER(IX,IY) < 0.01_JWRB .OR.                    &
     &        FIELDG%CICOVER(IX,IY) > 1.01_JWRB ) THEN
            CICVR(IJ) = 0.0_JWRB
          ELSEIF (FIELDG%CICOVER(IX,IY) .GT. 0.95_JWRB) THEN
            CICVR(IJ) = 1.0_JWRB
          ELSE
            CICVR(IJ) = FIELDG%CICOVER(IX,IY)
          ENDIF
        ENDDO
      ELSEIF (IPARAM == 139) THEN
        DO IJ=KIJS,KIJL
          IX = IFROMIJ(IJ)
          IY = JFROMIJ(IJ)
          IF (FIELDG%CICOVER(IX,IY) < 271.5_JWRB) THEN
            CICVR(IJ) = 1.0_JWRB
          ELSE
            CICVR(IJ) = 0.0_JWRB
          ENDIF
        ENDDO
      ENDIF 


      IF (.NOT. LICERUN .OR. LMASKICE) THEN
!       SEA ICE THICKNESS IN CASE IT IS NOT SUPPLIED AS INPUT:
        DO IJ=KIJS,KIJL
          CITH(IJ)=0.0_JWRB
        ENDDO
      ELSEIF (CLDOMAIN == 's') THEN
        DO IJ=KIJS,KIJL
          IF (CICVR(IJ) > 0.0_JWRB) THEN
            CITH(IJ)=SWAMPCITH
          ELSE
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO
      ELSEIF ((.NOT.LICETH) .AND. (.NOT.LWNEMOCOUCIT)) THEN
!       SEA ICE THICKNESS IS PARAMETERISED:
        DO IJ=KIJS,KIJL
          IF (CICVR(IJ) > 0.0_JWRB) THEN
            CITH(IJ)=MAX(C1+C2*CICVR(IJ),0.0_JWRB)
          ELSE
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ELSEIF (LWNEMOCOUCIT) THEN
        IF (LWCOU) THEN
          DO IJ=KIJS,KIJL
            IX = IFROMIJ(IJ)
            IY = JFROMIJ(IJ)
            IF (FIELDG%LKFR(IX,IY) <= 0.0_JWRB ) THEN
!             if lake cover = 0, we assume open ocean point, then get sea ice thickness directly from NEMO 
              IF (LNEMOICEREST) THEN
                CITH(IJ)=NEMOCITHICK(IJ)
              ELSE
                CITH(IJ)=CICVR(IJ)*NEMOCITHICK(IJ)
              ENDIF
            ELSE
!           We should get ice thickness information from atmopsheric model
!           but it is not yet coded. For now, parameterise it from the cover...
              IF (CICVR(IJ) > 0.0_JWRB) THEN
                CITH(IJ)=MAX(C1+C2*CICVR(IJ),0.0_JWRB)
              ELSE
               CITH(IJ)=0.0_JWRB
              ENDIF
            ENDIF
          ENDDO

        ELSE
          IF (LNEMOICEREST) THEN
            DO IJ=KIJS,KIJL
               CITH(IJ)=NEMOCITHICK(IJ)
            ENDDO
          ELSE
            DO IJ=KIJS,KIJL
               CITH(IJ)=CICVR(IJ)*NEMOCITHICK(IJ)
            ENDDO
          ENDIF
        ENDIF

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=KIJS,KIJL
          IF (CICVR(IJ) > 0.0_JWRB .AND. CITH(IJ) < 0.5_JWRB*HICMIN) THEN
            CICVR(IJ)=0.0_JWRB
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ELSE

!!!!   define a representative sea ice thickness that account for sea ice coverage
        DO IJ=KIJS,KIJL
          CITH(IJ)=CICVR(IJ)*CITH(IJ)
        ENDDO

!       CONSISTENCY CHECK:
!       no ice if thickness < 0.5*HICMIN
        DO IJ=KIJS,KIJL
          IF (CICVR(IJ) > 0.0_JWRB .AND. CITH(IJ) < 0.5_JWRB*HICMIN) THEN
            CICVR(IJ)=0.0_JWRB
            CITH(IJ)=0.0_JWRB
          ENDIF
        ENDDO

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MICEP',1,ZHOOK_HANDLE)

END SUBROUTINE MICEP
