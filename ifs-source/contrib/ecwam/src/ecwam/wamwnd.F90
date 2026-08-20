! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WAMWND (KIJS, KIJL,                 &
     &                   IFROMIJ, JFROMIJ,           &
     &                   NXS, NXE, NYS, NYE, FIELDG, &
     &                   UCUR, VCUR,                 &
     &                   U10, US,                    &
     &                   THW, ADS,                   &
     &                   WSTAR, CITH,                &
     &                   USTRA, VSTRA,               &
     &                   LWCUR, ICODE_WND)

! ----------------------------------------------------------------------

!     MODIFIED  S. ABDALLA OCTOBER 2001  INCLUSION OF AIR DENSITY & Zi/L

!     MODIFIED J. BIDLOT 2008 :: INCLUDE SURFACE CURRENTS FOR
!                                RELATIVE WIND CALCULATION. 
!                                USE OF NEUTRAL WIND SPEED FROM WAM

!**** *WAMWND* - TRANSFORMS INPUT FORCING TO BLOCKED WAM POINTS.


!*    PURPOSE.
!     --------

!       CONVERTS THE INTERPOLATED INPUT FIELDS TO WAM BLOCKS FOR ALL
!       POINTS IN THE GRID ON A PE, EXCEPT FOR U and V CURRENTS !!!

!**   INTERFACE.
!     ----------

!       *CALL WAMWND (KIJS, KIJL,
!                     IFROMIJ, JFROMIJ,
!                     NXS, NXE, NYS, NYE, FIELDG,
!                     UCUR, VCUR,
!                     U10, US,
!                     THW, ADS, WSTAR, CITH,
!                     USTRA, VSTRA,
!                     LWCUR, ICODE_WND)
!          *KIJS:KIJL* DIMENSION OF PASSED ARRAYS
!          *IFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *JFROMIJ*  POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          *NXS:NXE*  FIRST DIMENSION OF FIELDG
!          *NYS:NYE*  SECOND DIMENSION OF FIELDG
!          *FIELDG* - INPUT FORCING FIELDS ON THE WAVE MODEL GRID
!          *UCUR* - U-COMPONENT OF THE SURFACE CURRENT
!          *VCUR* - V-COMPONENT OF THE SURFACE CURRENT
!          *U10*  - INTERPOLATED WINDS AT ALL POINTS AND BLOCKS.
!          *US*   - INTERPOLATED FRICTION VELOCITY
!          *THW*  - INTERPOLATED WIND DIRECTION AT ALL POINTS.
!          *ADS*  - INTERPOLATED AIR DENSITY AT ALL POINTS.
!          *WSTAR*- INTERPOLATED CONVECTIVE VELOCITY AT ALL POINTS
!          *CITH* - SEA ICE THICKNESS. 
!          *USTRA*- U-COMPONENT OF THE ATMOSPHERIC SURFACE STRESS
!          *VSTRA*- V-COMPONENT OF THE ATMOSPHERIC SURFACE STRESS
!          *LWCUR*  - LOGICAL INDICATES THE PRESENCE OF SURFACE U AND V CURRENTS
!          *ICODE_WND* - INTEGER INDICATES WHAT IS SAVED IN FIELDG%UWND, and
!                        FIELDG%VWND
!                        AND
!                        WHAT IS UPDATED:
!                        ICODE_WND=3 THEN U10, otherwise US is UPDATED !!
!                        ++++++++++++++++++++++++++++++++++++++++++++++++


!     METHOD.
!     -------

!       THE INTERPOLATED VALUES ARE TRANSFORMED TO
!       MAGNITUDE AND DIRECTION. INPUT MAY BE WIND IN 10M HEIGHT ,
!       SURFACE WINDS OR FRICTION VELOCITIES. THE INPUT GRID HAS TO BE
!       A LATITUDE/LONGITUDE GRID EITHER PERIODIC OR NON PERIODIC.

!     EXTERNALS.
!     ----------

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FORCING_FIELDS

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWPCONS , ONLY : G        ,ZPI      ,ZMISS    ,EPSUS
      USE YOWPHYS  , ONLY : XKAPPA
      USE YOWSTAT  , ONLY : IREFRA   ,LRELWIND
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : WSPMIN   ,LLWSWAVE ,LLWDWAVE ,RWFAC

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: KIJS, KIJL
      INTEGER(KIND=JWIM), DIMENSION(KIJS:KIJL), INTENT(IN) :: IFROMIJ  ,JFROMIJ
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(IN) :: FIELDG
      INTEGER(KIND=JWIM), INTENT(IN) :: ICODE_WND
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(IN) :: UCUR, VCUR 
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(INOUT) :: U10, US
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(OUT) :: THW, ADS, WSTAR, CITH
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL), INTENT(OUT) :: USTRA, VSTRA 
      LOGICAL, INTENT(IN) :: LWCUR


      INTEGER(KIND=JWIM) :: IJ, IX, JY

      REAL(KIND=JWRB) :: RESCALE
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), DIMENSION (KIJS:KIJL) :: UU, VV, WSPEED

      LOGICAL :: LCORREL

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMWND',0,ZHOOK_HANDLE)

!     CORRECT FOR RELATIVE WINDS WITH RESPECT TO THE SURFACE CURRENTS.
      LCORREL=.FALSE.
      IF (LWCOU) THEN
!       In coupled experiments the winds are already relative to the currents
        IF (LWCUR .AND. (.NOT.LRELWIND) ) LCORREL=.TRUE.
      ELSE
        IF (LRELWIND .AND. (IREFRA == 2 .OR. IREFRA == 3)) LCORREL=.TRUE.
      ENDIF

!*    2. TRANSFORM GRIDDED WIND INPUT INTO BLOCK
!        ----------------------------------------

      DO IJ = KIJS, KIJL
        IX = IFROMIJ(IJ)
        JY = JFROMIJ(IJ)
        UU(IJ) = FIELDG%UWND(IX,JY)
        VV(IJ) = FIELDG%VWND(IX,JY)
        ADS(IJ) = FIELDG%AIRD(IX,JY)
        WSTAR(IJ)= FIELDG%WSTAR(IX,JY)
        CITH(IJ)= FIELDG%CITHICK(IX,JY)
        USTRA(IJ)= FIELDG%USTRA(IX,JY)
        VSTRA(IJ)= FIELDG%VSTRA(IX,JY)
      ENDDO


!*    3. PROCESS WINDS ACCORDING TO TYPE
!        ----------------------------------------

      SELECT CASE (ICODE_WND)
      CASE(3)
!       USE NEUTRAL WIND SPEED AND DIRECTION FROM A PREVIOUS WAM RUN.
        IF (LLWSWAVE .AND. LLWDWAVE) THEN
          DO IJ = KIJS, KIJL
            IX = IFROMIJ(IJ)
            JY = JFROMIJ(IJ)
            U10(IJ) = FIELDG%WSWAVE(IX,JY)
            THW(IJ) = FIELDG%WDWAVE(IX,JY)
          ENDDO

!         THERE MIGHT BE POINTS FOR WHICH NO WAM VALUES ARE AVAILABLE
!         USE U and V FROM ATMOSPHERE INSTEAD
          DO IJ = KIJS, KIJL
            IF (U10(IJ) <= 0.0_JWRB) THEN
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF (WSPEED(IJ) > 0.0_JWRB) THEN
                U10(IJ) = WSPEED(IJ)
                THW(IJ) = ATAN2(UU(IJ),VV(IJ))
              ELSE
                U10(IJ) = 0.0_JWRB
                THW(IJ) = 0.0_JWRB
              ENDIF
            ENDIF
          ENDDO

!         CORRECT FOR RELATIVE WINDS WITH RESPECT TO THE SURFACE CURRENTS.
          IF (LCORREL) THEN
            DO IJ = KIJS, KIJL
              UU(IJ) = U10(IJ)*SIN(THW(IJ))
              VV(IJ) = U10(IJ)*COS(THW(IJ))
              UU(IJ) = UU(IJ) - RWFAC*UCUR(IJ)
              VV(IJ) = VV(IJ) - RWFAC*VCUR(IJ)
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF (WSPEED(IJ) > 0.0_JWRB) THEN
                U10(IJ) = WSPEED(IJ)
                THW(IJ) = ATAN2(UU(IJ),VV(IJ))
              ELSE
                U10(IJ) = 0.0_JWRB
                THW(IJ) = 0.0_JWRB
              ENDIF
            ENDDO
          ENDIF

!       USE NEUTRAL WIND SPEED FROM A PREVIOUS WAM RUN.
!       -----------------------------------------------
        ELSEIF (LLWSWAVE) THEN

          DO IJ = KIJS, KIJL
            IX = IFROMIJ(IJ)
            JY = JFROMIJ(IJ)

            IF (FIELDG%WSWAVE(IX,JY) /= ZMISS .AND.                      &
     &          FIELDG%WSWAVE(IX,JY) > 0.0_JWRB ) THEN
              WSPEED(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
              IF (WSPEED(IJ) > 0.0_JWRB) THEN
                RESCALE=FIELDG%WSWAVE(IX,JY)/WSPEED(IJ)
                UU(IJ) = UU(IJ) * RESCALE 
                VV(IJ) = VV(IJ) * RESCALE
              ENDIF
            ENDIF
          ENDDO

          IF (LCORREL) THEN
            DO IJ = KIJS, KIJL
              UU(IJ) = UU(IJ) + RWFAC*UCUR(IJ)
              VV(IJ) = VV(IJ) + RWFAC*VCUR(IJ)
            ENDDO
          ENDIF

          DO IJ = KIJS, KIJL
            U10(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
            IF (U10(IJ) /= 0.0_JWRB) THEN 
              THW(IJ) = ATAN2(UU(IJ),VV(IJ))
            ELSE
              THW(IJ) = 0.0_JWRB
            ENDIF
          ENDDO


        ELSE

          IF (LCORREL) THEN
            DO IJ = KIJS, KIJL
              UU(IJ) = UU(IJ) + RWFAC*UCUR(IJ)
              VV(IJ) = VV(IJ) + RWFAC*VCUR(IJ)
            ENDDO
          ENDIF

          DO IJ = KIJS, KIJL
            U10(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
            IF (U10(IJ) /= 0.0_JWRB) THEN 
              THW(IJ) = ATAN2(UU(IJ),VV(IJ))
            ELSE
              THW(IJ) = 0.0_JWRB
            ENDIF
          ENDDO

        ENDIF

!       IMPOSE A MINIMUM WIND SPEED.
        DO IJ = KIJS, KIJL
          U10(IJ) = MAX(U10(IJ), WSPMIN)
        ENDDO

      CASE(1)

!*    3.2  INPUT IS FRICTION VELOCITY.
!          ---------------------------

        DO IJ = KIJS, KIJL
          US(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
          IF (US(IJ) /= 0.0_JWRB) THEN 
            THW(IJ) = ATAN2(UU(IJ),VV(IJ))
          ELSE
            THW(IJ) = 0.0_JWRB
          ENDIF
          US(IJ) = MAX(US(IJ),EPSUS)
        ENDDO

      CASE(2)

!*    3.3 INPUT WINDS ARE SURFACE STRESSES.
!         ---------------------------------

        DO IJ = KIJS, KIJL
          US(IJ) = SQRT(UU(IJ)**2 + VV(IJ)**2)
          IF (US(IJ) /= 0.0_JWRB) THEN 
            THW(IJ) = ATAN2(UU(IJ),VV(IJ))
          ELSE
            THW(IJ) = 0.0_JWRB
          ENDIF
          US(IJ) = SQRT(MAX(US(IJ),0.0_JWRB)/MAX(ADS(IJ),1.))
          US(IJ) = MAX(US(IJ),EPSUS)
        ENDDO

      END SELECT 

      DO IJ = KIJS, KIJL
        IF (THW(IJ) < 0.0_JWRB) THW(IJ) = THW(IJ) + ZPI
      ENDDO

! ----------------------------------------------------------------------

!*    4. TEST OUTPUT OF WAVE MODEL BLOCKS
!        ---------------------------------

      IF (LHOOK) CALL DR_HOOK('WAMWND',1,ZHOOK_HANDLE)

      END SUBROUTINE WAMWND
