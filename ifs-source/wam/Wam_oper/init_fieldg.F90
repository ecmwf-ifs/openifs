! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE INIT_FIELDG(BLK2LOC, LLINIALL, LLOCAL,     &
 &                     NXS, NXE, NYS, NYE, FIELDG) 

! ----------------------------------------------------------------------

!**** *INIT_FIELDG* - INITIALISES DATA STRUCTURE FIELDG 

!*    PURPOSE.
!     --------

!     INITIALISES THE DATA STRUCTURE FIELDG 

!**   INTERFACE.
!     ----------

!     *CALL* *INIT_FIELDG(BLK2LOC,LLINIALL,LLOCAL,
!                         NXS, NXE, NYS, NYE, FIELDG)

!          BLK2LOC   POINTERS FROM LOCAL GRID POINTS TO 2-D MAP
!          LLINIALL  LOGICAL : IF TRUE ALL STRUCTURE IS INITIALISED
!                              OTHERWISE ONLY XLON and YLAT
!          LLOCAL    LOGICAL : IF TRUE ONLY THE LOCAL SEA POINTS ARE GIVEN
!                              VALID COORDINATES (XLON, YLAT) 
!                              PROVIDED THE LOCAL INDEX ARRAYS EXIST.
!         *NXS:NXE*  FIRST DIMENSION OF FIELDG
!         *NYS:NYE*  SECOND DIMENSION OF FIELDG
!          FIELDG    INPUT FORCING FIELDS ON THE WAVE MODEL GRID

!     REFERENCE.
!     -----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDLOC, FORCING_FIELDS

      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK, KIJL4CHNK
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO, NLONRGG
      USE YOWPARAM , ONLY : NGY
      USE YOWPCONS , ONLY : ZMISS    ,ROAIR    ,WSTAR0
      USE YOWPHYS  , ONLY : PRCHAR
      USE YOWWIND  , ONLY : WSPMIN
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWPD    , ONLY : MNP=>npa , XP=>x, YP=>y

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

      TYPE(WVGRIDLOC), DIMENSION(NPROMA_WAM, NCHNK), INTENT(IN) :: BLK2LOC
      LOGICAL, INTENT(IN) :: LLINIALL, LLOCAL
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), DIMENSION(NXS:NXE, NYS:NYE), INTENT(INOUT) :: FIELDG


      INTEGER(KIND=JWIM) :: IJ, IX, JY, JSN
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',0,ZHOOK_HANDLE)


      CALL GSTATS(1501,0)
      IF (LLINIALL) THEN
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY, IX)
        DO JY = NYS, NYE
          DO IX = NXS, NXE
            FIELDG(IX,JY)%XLON    = ZMISS 
            FIELDG(IX,JY)%YLAT    = ZMISS
            FIELDG(IX,JY)%UWND    = 0.0_JWRB
            FIELDG(IX,JY)%VWND    = 0.0_JWRB
            FIELDG(IX,JY)%WSWAVE  = WSPMIN 
            FIELDG(IX,JY)%WDWAVE  = 0.0_JWRB
            FIELDG(IX,JY)%UFRIC   = 0.0_JWRB
            FIELDG(IX,JY)%CICOVER = 0.0_JWRB
            FIELDG(IX,JY)%CITHICK = 0.0_JWRB
            FIELDG(IX,JY)%LKFR    = 0.0_JWRB
            FIELDG(IX,JY)%AIRD    = ROAIR
            FIELDG(IX,JY)%WSTAR   = WSTAR0
            FIELDG(IX,JY)%UCUR    = 0.0_JWRB
            FIELDG(IX,JY)%VCUR    = 0.0_JWRB
            FIELDG(IX,JY)%TAUW    = 0.0_JWRB
            FIELDG(IX,JY)%TAUWDIR = 0.0_JWRB
            FIELDG(IX,JY)%Z0M     = 0.0_JWRB
            FIELDG(IX,JY)%Z0B     = 0.0_JWRB
            FIELDG(IX,JY)%CHRNCK  = PRCHAR 
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY, IX)
        DO JY = NYS, NYE
          DO IX = NXS, NXE
            FIELDG(IX,JY)%XLON   = ZMISS 
            FIELDG(IX,JY)%YLAT   = ZMISS
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF

!     COORDINATES OF THE POINTS IN FIELDG THAT ARE NEEDED

      IF (LLOCAL) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IJ, IX, JY, JSN)
        DO ICHNK = 1, NCHNK
          KIJS=1
          KIJL=KIJL4CHNK(ICHNK)
          DO IJ=KIJS,KIJL
            IX = BLK2LOC(IJ,ICHNK)%IFROMIJ
            JY = BLK2LOC(IJ,ICHNK)%JFROMIJ
            JSN= BLK2LOC(IJ,ICHNK)%KFROMIJ
            IF (LLUNSTR) THEN
              FIELDG(IX,JY)%XLON = XP(IJ)
              FIELDG(IX,JY)%YLAT = YP(IJ)
            ELSE
              FIELDG(IX,JY)%XLON = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG(IX,JY)%YLAT = AMOSOP + (JSN-1)*XDELLA
            ENDIF
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO

      ELSE
        IF (LLUNSTR) THEN
!!!!
          write(*,*) 'In INIT_FIELDG : not yet ready for unstructured grid '
          call abort1

        ELSE
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY, IX, JSN)
          DO JY = NYS, NYE
            JSN = NGY-JY+1
            DO IX = NXS, MIN(NLONRGG(JSN), NXE)
              FIELDG(IX,JY)%XLON = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG(IX,JY)%YLAT = AMOSOP + (JSN-1)*XDELLA
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDIF

      ENDIF
      CALL GSTATS(1501,1)


IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_FIELDG
