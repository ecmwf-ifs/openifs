! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

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
      USE YOWMAP   , ONLY : AMOWEP   ,AMOSOP   ,XDELLA   ,ZDELLO, NLONRGG, NGY
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPCONS , ONLY : ZMISS    ,ROAIR    ,WSTAR0
      USE YOWPHYS  , ONLY : PRCHAR
      USE YOWWIND  , ONLY : WSPMIN
#ifdef WAM_HAVE_UNWAM
      USE YOWPD    , ONLY : XP=>x, YP=>y
#endif
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"

      TYPE(WVGRIDLOC), INTENT(IN) :: BLK2LOC
      LOGICAL, INTENT(IN) :: LLINIALL, LLOCAL
      INTEGER(KIND=JWIM), INTENT(IN) :: NXS, NXE, NYS, NYE
      TYPE(FORCING_FIELDS), INTENT(INOUT) :: FIELDG


      INTEGER(KIND=JWIM) :: IJ, IX, JY, JSN
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',0,ZHOOK_HANDLE)

IF (LLUNSTR) THEN
#ifndef WAM_HAVE_UNWAM
  CALL WAM_ABORT("UNWAM support not available",__FILENAME__,__LINE__)
#endif
ENDIF


      CALL GSTATS(1501,0)
      IF (LLINIALL) THEN
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY, IX)
        DO JY = NYS, NYE
          DO IX = NXS, NXE
            FIELDG%XLON(IX,JY)    = ZMISS
            FIELDG%YLAT(IX,JY)    = ZMISS
            FIELDG%UWND(IX,JY)    = 0.0_JWRB
            FIELDG%VWND(IX,JY)    = 0.0_JWRB
            FIELDG%WSWAVE(IX,JY)  = WSPMIN
            FIELDG%WDWAVE(IX,JY)  = 0.0_JWRB
            FIELDG%UFRIC(IX,JY)   = 0.0_JWRB
            FIELDG%CICOVER(IX,JY) = 0.0_JWRB
            FIELDG%CITHICK(IX,JY) = 0.0_JWRB
            FIELDG%LKFR(IX,JY)    = 0.0_JWRB
            FIELDG%AIRD(IX,JY)    = ROAIR
            FIELDG%WSTAR(IX,JY)   = WSTAR0
            FIELDG%USTRA(IX,JY)   = 0.0_JWRB
            FIELDG%VSTRA(IX,JY)   = 0.0_JWRB
            FIELDG%UCUR(IX,JY)    = 0.0_JWRB
            FIELDG%VCUR(IX,JY)    = 0.0_JWRB
            FIELDG%TAUW(IX,JY)    = 0.0_JWRB
            FIELDG%TAUWDIR(IX,JY) = 0.0_JWRB
            FIELDG%Z0M(IX,JY)     = 0.0_JWRB
            FIELDG%Z0B(IX,JY)     = 0.0_JWRB
            FIELDG%CHRNCK(IX,JY)  = PRCHAR
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ELSE
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(JY, IX)
        DO JY = NYS, NYE
          DO IX = NXS, NXE
            FIELDG%XLON(IX,JY)   = ZMISS
            FIELDG%YLAT(IX,JY)   = ZMISS
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
            IX = BLK2LOC%IFROMIJ(IJ,ICHNK)
            JY = BLK2LOC%JFROMIJ(IJ,ICHNK)
            JSN= BLK2LOC%KFROMIJ(IJ,ICHNK)
            IF (LLUNSTR) THEN
#ifdef WAM_HAVE_UNWAM
              FIELDG%XLON(IX,JY) = XP(IJ)
              FIELDG%YLAT(IX,JY) = YP(IJ)
#endif
            ELSE
              FIELDG%XLON(IX,JY) = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG%YLAT(IX,JY) = AMOSOP + (JSN-1)*XDELLA
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
              FIELDG%XLON(IX,JY) = AMOWEP + (IX-1)*ZDELLO(JSN)
              FIELDG%YLAT(IX,JY) = AMOSOP + (JSN-1)*XDELLA
            ENDDO
          ENDDO
!$OMP     END PARALLEL DO
        ENDIF

      ENDIF
      CALL GSTATS(1501,1)


IF (LHOOK) CALL DR_HOOK('INIT_FIELDG',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_FIELDG
