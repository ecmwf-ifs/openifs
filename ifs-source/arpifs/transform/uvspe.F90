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

SUBROUTINE UVSPE(YDGEOMETRY,PETPSI,PDIKHI,PU,PV,KFLSUR,KFIELD,KOUTP)

!**** *UVSPE - Compute spectral vorticity and divergence from physical
!                space u and v.

!     Purpose.
!     --------
!     Compute vorticity and divergence or psi and khi from u and v.

!**   Interface.
!     ----------
!        *CALL* *UVSPE(...)

!        Explicit arguments :
!        --------------------
!           PETPSI             : vorticity/stream function (output)
!           PDIKHI             : divergence/vel.potential  (output)
!           PU                 : u velocity (input)
!           PV                 : v velocity (input)
!           KFLSUR             : first dimension of spectral arrays (input)
!           KFIELD             : Number of fields (input)
!           KOUTP              : KOUTP=1 output vor/div (input)
!                              : KOUTP=2 output phi/kai

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation about spectral transforms.

!     Externals.  DIR_TRANS - direct transform (TRANS library)
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE about spectral transforms.

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 97-06-11 Based on UVSPE

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV, NPRTRV

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PETPSI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDIKHI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(YDGEOMETRY%YRGEM%NGPTOT,KFIELD) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PV(YDGEOMETRY%YRGEM%NGPTOT,KFIELD) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOUTP 

REAL(KIND=JPRB) :: ZUV(YDGEOMETRY%YRGEM%NGPTOT,2*KFIELD,1)

INTEGER(KIND=JPIM) :: ILF,ISTG,IFLDPP,IREST,IGF
INTEGER(KIND=JPIM) :: INUMLF(NPRTRV),IVSETUV(KFIELD)
INTEGER(KIND=JPIM) :: JROCB,JLEV,JSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "dir_trans.h"

!     ------------------------------------------------------------------

!*    INITIALISATION

IF (LHOOK) CALL DR_HOOK('UVSPE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NRESOL=>YDDIM%NRESOL, NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & NVALUE=>YDLAP%NVALUE, RLAPIN=>YDLAP%RLAPIN)
ILF = 0
IF( KFIELD == 1 )THEN
  IVSETUV(1) = NPRTRV
  IF(NPRTRV == MYSETV) THEN
    ILF = 1
  ENDIF
  IGF = 1
ELSE
  IFLDPP = KFIELD/NPRTRV
  IREST = KFIELD-IFLDPP*NPRTRV
  DO JROCB=1,NPRTRV
    IF( JROCB <= IREST )THEN
      INUMLF(JROCB) = IFLDPP+1
    ELSE
      INUMLF(JROCB) = IFLDPP
    ENDIF
  ENDDO
  ILF = INUMLF(MYSETV)
  ISTG = 1
  DO JROCB=1,NPRTRV
    IVSETUV(ISTG:ISTG+INUMLF(JROCB)-1) = JROCB
    ISTG = ISTG+INUMLF(JROCB)
  ENDDO
  IGF = ISTG-1
ENDIF

ZUV(:,1:KFIELD,1) = PU(:,:)
ZUV(:,KFIELD+1:2*KFIELD,1) = PV(:,:)
CALL DIR_TRANS(PSPVOR=PETPSI(1:ILF,:),PSPDIV=PDIKHI(1:ILF,:),&
 & KVSETUV=IVSETUV(1:IGF),KRESOL=NRESOL,KPROMA=NGPTOT,PGP=ZUV)  

IF (KOUTP == 2) THEN
  DO JLEV=1,ILF
    DO JSP=1,NSPEC2
      PETPSI(JLEV,JSP) = PETPSI(JLEV,JSP)*RLAPIN(NVALUE(JSP))
      PDIKHI(JLEV,JSP) = PDIKHI(JLEV,JSP)*RLAPIN(NVALUE(JSP))
    ENDDO
  ENDDO
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UVSPE',1,ZHOOK_HANDLE)
END SUBROUTINE UVSPE
