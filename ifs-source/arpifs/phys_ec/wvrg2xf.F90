! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE WVRG2XF(YDSURF,P_ZFLOCAL,K_NWVFIELDS,KGPTOT,KPROMA)

!**** *WVRG2XF*  - stores wave fields into work array 

!     Purpose.
!     --------
!                 stores wave fields into work array

!**   Interface.
!     ----------
!        *CALL* *WVRG2XF(...)*

!        Explicit arguments : 
!        ------------------
!            ZFLOCAL   - The fields to be passed from WAVE to ATM
!                        It only contains the contribution for local
!                        grid points.
!            NWVFIELDS - Number of fields passed from WAVE to ATM

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!      LARS ISAKSEN *ECMWF*
!      JIM DOYLE

!     Modifications:
!     --------------
!      S. Abdalla 01-11-27: Generalised interface.
!      J. Bidlot  01-11-20: work only on local grid points.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!      M.Hamrud      01-Jul-2006  Revised surface fields
!     ------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TSURF), INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: K_NWVFIELDS 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZFLOCAL(KGPTOT,K_NWVFIELDS) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPTOT
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IFLD
INTEGER(KIND=JPIM) :: IEND, IBL, IST, JROF, JSTGLO

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WVRG2XF',0,ZHOOK_HANDLE)
ASSOCIATE(SD_WS=>YDSURF%SD_WS, YSD_WS=>YDSURF%YSD_WS)
!     ------------------------------------------------------------------

CALL GSTATS(1426,0)
!$OMP PARALLEL DO SCHEDULE(STATIC)&
!$OMP&PRIVATE(JSTGLO,IST,IEND,IBL,IFLD,JROF)
DO JSTGLO=1,KGPTOT,KPROMA
  IST    = 1
  IEND   = MIN(KPROMA,KGPTOT-JSTGLO+1)
  IBL    = (JSTGLO-1)/KPROMA+1

  DO IFLD=1, K_NWVFIELDS
    DO JROF =IST,IEND
      SD_WS(JROF,YSD_WS%YWS(IFLD)%MP,IBL) = P_ZFLOCAL(JSTGLO+JROF-1,IFLD)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1426,1)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WVRG2XF',1,ZHOOK_HANDLE)
END SUBROUTINE WVRG2XF
