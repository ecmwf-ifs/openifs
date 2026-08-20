! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE NEWWIND (CDATE, CDATEWH, LLNEWFILE,           &
 &                  WVPRPT, FF_NOW, FF_NEXT)
! ----------------------------------------------------------------------

!**** *NEWWIND* - HANDLING OF WINDS OUTSIDE MULTITASKED AREA    

!     P.A.E.M. JANSSEN  KNMI/ECMWF  SEPTEMBER 1994
!     J. BIDLOT         ECMWF       FEBRUARY 1996  MESSAGE PASSING 
!     S. ABDALLA        ECMWF       OCTOBER 2001   INCLUSION OF AIR
!                                                  DENSITY AND Zi/L
!     J. BIDLOT         ECMWF       AUGUST 2008: CLEAN-UP 

!*    PURPOSE.
!     --------

!       GETS FORCING FIELDS WHEN NEEDED.                                   

!**   INTERFACE.
!     ----------

!       *CALL* *NEWWIND (CDATE, CDATEWH, LLNEWFILE,
!                        WVPRPT, FF_NOW, FF_NEXT,
!      *CDATE*    - START DATE OF SOURCE FUNCTION INTEGRATION
!      *CDATEWH*  - DATE OF THE NEXT FORCING FIELDS
!      *LLNEWFILE*- TRUE IF NEW WIND FILE HAS BEEN OPENED
!      *WVPRPT*   - WAVE PROPERTIES FIELDS
!      *FF_NOW*   - DATA STRUCTURE WITH THE CURRENT FORCING FIELDS
!      *FF_NEXT*  - DATA STRUCTURE WITH THE NEXT FORCING FIELDS

!     METHOD.
!     -------

!       READ NEW WINDS AND ASSOCIATED FIELDS WHEN NEEDED.

!     EXTERNALS.
!     ---------

!       *INCDATE*   - UPDATE DATE TIME GROUP.

!     REFERENCE.
!     ----------

!       NONE

! ----------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : FREQUENCY, FORCING_FIELDS

      USE YOWPCONS , ONLY : ACD      ,BCD      ,EPSMIN
      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWPARAM , ONLY : NFRE
      USE YOWPHYS  , ONLY : ALPHA
      USE YOWSTAT  , ONLY : IDELWO
      USE YOWTEST  , ONLY : IU06
      USE YOWWIND  , ONLY : CDATEWL  ,CDAWIFL  ,CDATEFL  ,CDTNEXT  ,         &
     &                      NSTORE   ,WSPMIN_RESET_TAUW  ,USTMIN_RESET_TAUW
      USE YOWWNDG  , ONLY : ICODE    ,ICODE_CPL

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "cireduce.intfb.h"
#include "incdate.intfb.h"

      CHARACTER(LEN=14), INTENT(IN)    :: CDATE
      CHARACTER(LEN=14), INTENT(INOUT) :: CDATEWH
      LOGICAL, INTENT(INOUT)           :: LLNEWFILE
      TYPE(FREQUENCY), INTENT(INOUT) :: WVPRPT
      TYPE(FORCING_FIELDS), INTENT(INOUT)  :: FF_NOW
      TYPE(FORCING_FIELDS), INTENT(IN)     :: FF_NEXT


      INTEGER(KIND=JWIM) :: ICODE_WND
      INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJ

      REAL(KIND=JWRB) :: WGHT, TLWMAX
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('NEWWIND',0,ZHOOK_HANDLE)

      IF (LWCOU) THEN
        ICODE_WND = ICODE_CPL
      ELSE
        ICODE_WND = ICODE
      ENDIF

      WGHT = 1.0_JWRB/MAX(WSPMIN_RESET_TAUW,EPSMIN)

      IF (CDATE >= CDATEWH) THEN

!*    2. NEW WIND INPUT.
!        ---------------
        IF (CDATE >= CDATEFL) THEN
            LLNEWFILE = .TRUE.
        ENDIF

!*    2.2 NEW WINDS ARE READ IN.
!         ----------------------

        CDATEWL = CDTNEXT

        CALL GSTATS(1492,0)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(ICHNK, KIJS, KIJL, IJ, TLWMAX)
        DO ICHNK = 1, NCHNK
          KIJS = 1
          KIJL = NPROMA_WAM
          IF (ICODE_WND == 3 ) THEN
            DO IJ = KIJS, KIJL
              FF_NOW%WSWAVE(IJ,ICHNK) = FF_NEXT%WSWAVE(IJ,ICHNK)
! adapt first estimate of wave induced stress for low winds
! to a fraction of the simple relation u*^2 = Cd(U10) * U10^2
! where this fraction varies from 0 for U10=0 to 1 for U10=WSPMIN_RESET_TAUW
              IF (FF_NOW%WSWAVE(IJ,ICHNK) < WSPMIN_RESET_TAUW) THEN
                TLWMAX = WGHT * (ACD+BCD*FF_NOW%WSWAVE(IJ,ICHNK)) * FF_NOW%WSWAVE(IJ,ICHNK)**3
                FF_NOW%TAUW(IJ,ICHNK) = MIN(FF_NOW%TAUW(IJ,ICHNK), TLWMAX)
              ENDIF
            ENDDO
          ELSE
            DO IJ = KIJS, KIJL
              FF_NOW%UFRIC(IJ,ICHNK) = FF_NEXT%UFRIC(IJ,ICHNK)
! update the estimate of TAUW
              FF_NOW%TAUW(IJ,ICHNK) = FF_NOW%UFRIC(IJ,ICHNK)**2 * (1.0_JWRB-(ALPHA/FF_NOW%CHRNCK(IJ,ICHNK))**2)
! adapt first estimate of wave induced stress for low winds
              IF (FF_NOW%UFRIC(IJ,ICHNK) < USTMIN_RESET_TAUW) FF_NOW%TAUW(IJ,ICHNK) = 0.0_JWRB
            ENDDO
          ENDIF

          FF_NOW%WDWAVE(KIJS:KIJL,ICHNK)  = FF_NEXT%WDWAVE(KIJS:KIJL,ICHNK)
          FF_NOW%AIRD(KIJS:KIJL,ICHNK)    = FF_NEXT%AIRD(KIJS:KIJL,ICHNK)
          FF_NOW%WSTAR(KIJS:KIJL,ICHNK)   = FF_NEXT%WSTAR(KIJS:KIJL,ICHNK)
          FF_NOW%CICOVER(KIJS:KIJL,ICHNK) = FF_NEXT%CICOVER(KIJS:KIJL,ICHNK)
          FF_NOW%CITHICK(KIJS:KIJL,ICHNK) = FF_NEXT%CITHICK(KIJS:KIJL,ICHNK)
          FF_NOW%USTRA(KIJS:KIJL,ICHNK)   = FF_NEXT%USTRA(KIJS:KIJL,ICHNK)
          FF_NOW%VSTRA(KIJS:KIJL,ICHNK)   = FF_NEXT%VSTRA(KIJS:KIJL,ICHNK)

        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1492,1)


        CALL INCDATE(CDATEWH, IDELWO)

!       UPDATE THE SEA ICE REDUCTION FACTOR
        CALL CIREDUCE (WVPRPT, FF_NOW)

      ENDIF

IF (LHOOK) CALL DR_HOOK('NEWWIND',1,ZHOOK_HANDLE)

END SUBROUTINE NEWWIND
