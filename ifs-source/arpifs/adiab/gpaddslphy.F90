! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE GPADDSLPHY(YDGEOMETRY,YDSLPHY,YGFL,LDLFSTEP,KST,KEND,PTMST, &
 & PSAVTEND,PGFLSLP,&
 & PSLB1U0,PSLB1V0,PSLB1T0,PSLB1GFL0)  

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMSLPHY     , ONLY : TSLPHY
USE YOM_YGFL     , ONLY : TYPE_GFLD

!**** *GPADDSLPHY - add physics tendencies to RHS of semi-lagrangian buffers

!     Purpose. add physics tendencies to RHS 
!     -------- 

!**   Interface.
!     ----------
!        *CALL* *GPADDSLPHY(...)*

!        Explicit arguments : see CPG
!        -------------------- 

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        None

!     Author.
!     -------
!        Nils Wedi  *ECMWF*

!     Modifications.
!     --------------
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   F.Vana        11-Sep-2020 Cleaning
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSLPHY)      ,INTENT(IN)    :: YDSLPHY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
LOGICAL           ,INTENT(IN)    :: LDLFSTEP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSAVTEND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDSLPHY%NVTEND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFLSLP(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIMSLP) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLB1U0(YDGEOMETRY%YRDIM%NPROMA, &
 & YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLB1V0(YDGEOMETRY%YRDIM%NPROMA, &
 & YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLB1T0(YDGEOMETRY%YRDIM%NPROMA, &
 & YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSLB1GFL0(YDGEOMETRY%YRDIM%NPROMA, &
 & YDGEOMETRY%YRDIMV%NFLSA:YDGEOMETRY%YRDIMV%NFLEN, YGFL%NDIMSLP) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPADDSLPHY',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, NFLEVG=>YDDIMV%NFLEVG, &
 & MT_SAVTEND=>YDSLPHY%MT_SAVTEND, MU_SAVTEND=>YDSLPHY%MU_SAVTEND, &
 & MV_SAVTEND=>YDSLPHY%MV_SAVTEND, NVTEND=>YDSLPHY%NVTEND)
!     ------------------------------------------------------------------

IF(.NOT.LDLFSTEP) THEN
  DO JLEV=1,NFLEVG
    DO JROF=KST,KEND
      PSLB1U0(JROF,JLEV)=PSLB1U0(JROF,JLEV)+PSAVTEND(JROF,JLEV,MU_SAVTEND)*PTMST
      PSLB1V0(JROF,JLEV)=PSLB1V0(JROF,JLEV)+PSAVTEND(JROF,JLEV,MV_SAVTEND)*PTMST
      PSLB1T0(JROF,JLEV)=PSLB1T0(JROF,JLEV)+PSAVTEND(JROF,JLEV,MT_SAVTEND)*PTMST
    ENDDO
  ENDDO

  DO JGFL=1,NUMFLDS
    IF(YCOMP(JGFL)%LPHY) THEN
      DO JLEV=1,NFLEVG
        DO JROF=KST,KEND
          PSLB1GFL0(JROF,JLEV,YCOMP(JGFL)%MPSLP)=PSLB1GFL0(JROF,JLEV,YCOMP(JGFL)%MPSLP)+&
           & PGFLSLP(JROF,JLEV,YCOMP(JGFL)%MPSLP)*PTMST  
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPADDSLPHY',1,ZHOOK_HANDLE)
END SUBROUTINE GPADDSLPHY
