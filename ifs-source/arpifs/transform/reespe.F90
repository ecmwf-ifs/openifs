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

SUBROUTINE REESPE(YDGEOMETRY,KFLSUR,KFIELD,PSPEC,PREEL)

!**** *REESPE* -Direct transform to Spherical harmonics

!     Purpose.
!     --------
!     Direct transform from grid-point to spectral space.

!**   Interface.
!     ----------
!        *CALL* *REESPE(...)

!        Explicit arguments :
!        --------------------
!           KFLSUR   : first dimension of PSPEC
!           KFIELD   : number of fields to be transformed
!           PSPEC    : spectral array (output)
!           PREEL    : grid point array (input)

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.  DIR_TRANS - direct transform (TRANS lib)
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : MYSETV, NPRTRV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPEC(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREEL(YDGEOMETRY%YRGEM%NGPTOT,KFIELD,1)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ILF,ISTG,IFLDPP,IREST,IGF
INTEGER(KIND=JPIM) :: INUMLF(NPRTRV),IVSETSC(KFIELD)
INTEGER(KIND=JPIM) :: JROCB
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "dir_trans.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('REESPE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NRESOL=>YDDIM%NRESOL, NSPEC2=>YDDIM%NSPEC2, &
 & NGPTOT=>YDGEM%NGPTOT)
!     ------------------------------------------------------------------

!  Initialisation
ILF = 0
IF( KFIELD == 1 )THEN
  IVSETSC(1) = NPRTRV
  IF(NPRTRV == MYSETV) THEN
    ILF = 1
  ENDIF
  ISTG=2
ELSE
  IFLDPP=KFIELD/NPRTRV
  IREST=KFIELD-IFLDPP*NPRTRV
  DO JROCB=1,NPRTRV
    IF( JROCB <= IREST )THEN
      INUMLF(JROCB)=IFLDPP+1
    ELSE
      INUMLF(JROCB)=IFLDPP
    ENDIF
  ENDDO
  ILF = INUMLF(MYSETV)
  ISTG = 1
  DO JROCB=1,NPRTRV
    IVSETSC(ISTG:ISTG+INUMLF(JROCB)-1) = JROCB
    ISTG = ISTG+INUMLF(JROCB)
  ENDDO
ENDIF
IGF = ISTG-1

CALL DIR_TRANS(PSPSCALAR=PSPEC(1:ILF,:),KRESOL=NRESOL,KPROMA=NGPTOT,&
 & KVSETSC=IVSETSC(1:IGF),PGP=PREEL)  

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('REESPE',1,ZHOOK_HANDLE)
END SUBROUTINE REESPE

