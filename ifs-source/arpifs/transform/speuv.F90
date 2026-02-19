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

SUBROUTINE SPEUV(YDGEOMETRY,PETPSI,PDIKHI,PU,PV,KFLSUR,KFIELD,KINP)

!**** *SPEUV - Transform from spectral fields to winds in physical space - DM

!     Purpose.
!     --------

!     Compute wind from vorticity/divergence or streamfunction/
!     velocity potential.

!**   Interface.
!     ----------
!        *CALL* *SPEUV(...)

!        Explicit arguments :
!        --------------------
!           PETPSI                     : vorticity/stream function (input)
!           PDIKHI                     : divergence/vel.pot.       (input)
!           PU                         : u velocity                (output)
!           PV                         : v velocity                (output)
!           KFLSUR                     : first dimension of spectral arrays
!           KFIELD                     : Number of fields
!           KINP                       : KINP=1 - input vor/div
!                                      : KINP=2 - input psi/kai

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals. INV_TRANS - inverse transform (TRANS lib)
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE

!     Author.
!     -------
!      MPP Group *ECMWF*
!      Original : 95-10-01 Based on SPEUV 

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : MYSETV, NPRTRV

IMPLICIT NONE

TYPE(GEOMETRY),    INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETPSI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIKHI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PU(YDGEOMETRY%YRGEM%NGPTOT,KFIELD) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PV(YDGEOMETRY%YRGEM%NGPTOT,KFIELD) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KINP 
REAL(KIND=JPRB) :: ZETPSI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)         , ZDIKHI(KFLSUR,YDGEOMETRY%YRDIM%NSPEC2)
REAL(KIND=JPRB) :: ZUV(YDGEOMETRY%YRGEM%NGPTOT,2*KFIELD,1)

INTEGER(KIND=JPIM) :: ILF,ISTG,IFLDPP,IREST,IGF,IM,ISP
INTEGER(KIND=JPIM) :: INUMLF(NPRTRV),IVSETUV(KFIELD)
INTEGER(KIND=JPIM) :: JROCB,JMLOC,JLEV,JN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "inv_trans.h"

#include "abor1.intfb.h"

!  INITIALISATION

IF (LHOOK) CALL DR_HOOK('SPEUV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NRESOL=>YDDIM%NRESOL, NSMAX=>YDDIM%NSMAX, NSPEC2=>YDDIM%NSPEC2, &
 & NUMP=>YDDIM%NUMP, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & MYMS=>YDLAP%MYMS, NASM0=>YDLAP%NASM0, RLAPDI=>YDLAP%RLAPDI)
ILF = 0
IF( KFIELD == 1 )THEN
  IVSETUV(1) = NPRTRV
  IF(NPRTRV == MYSETV) THEN
    ILF = 1
  ENDIF
  IGF = 1
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
    IVSETUV(ISTG:ISTG+INUMLF(JROCB)-1) = JROCB
    ISTG = ISTG+INUMLF(JROCB)
  ENDDO
  IGF = ISTG-1
ENDIF

!     PREPARE THE FIELDS AND RESOLVE EQ OF POISSON

IF (KINP == 1) THEN

!    INPUT VORTICITY/DIVERGENCE.

  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    IF( IM == 0 )THEN
      DO JLEV=1,ILF
        DO JN=IM,NSMAX
          ISP=NASM0(IM)+(JN-IM)*2
          ZETPSI(JLEV,ISP  ) = PETPSI(JLEV,ISP)
          ZETPSI(JLEV,ISP+1) = 0.0_JPRB
          ZDIKHI(JLEV,ISP  ) = PDIKHI(JLEV,ISP)
          ZDIKHI(JLEV,ISP+1) = 0.0_JPRB
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,ILF
        DO JN=IM,NSMAX
          ISP=NASM0(IM)+(JN-IM)*2
          ZETPSI(JLEV,ISP  ) = PETPSI(JLEV,ISP  )
          ZETPSI(JLEV,ISP+1) = PETPSI(JLEV,ISP+1)
          ZDIKHI(JLEV,ISP  ) = PDIKHI(JLEV,ISP  )
          ZDIKHI(JLEV,ISP+1) = PDIKHI(JLEV,ISP+1)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
ELSEIF(KINP == 2) THEN

!    INPUT STREAMFUNCTION/VELOCITY POTENTIAL.

  DO JMLOC=1,NUMP
    IM = MYMS(JMLOC)
    IF( IM == 0 )THEN
      DO JLEV=1,ILF
        DO JN=IM,NSMAX
          ISP=NASM0(IM)+(JN-IM)*2
          ZETPSI(JLEV,ISP  ) = PETPSI(JLEV,ISP)*RLAPDI(JN)
          ZETPSI(JLEV,ISP+1) = 0.0_JPRB
          ZDIKHI(JLEV,ISP  ) = PDIKHI(JLEV,ISP)*RLAPDI(JN)
          ZDIKHI(JLEV,ISP+1) = 0.0_JPRB
        ENDDO
      ENDDO
    ELSE
      DO JLEV=1,ILF
        DO JN=IM,NSMAX
          ISP=NASM0(IM)+(JN-IM)*2
          ZETPSI(JLEV,ISP  ) = PETPSI(JLEV,ISP  )*RLAPDI(JN)
          ZETPSI(JLEV,ISP+1) = PETPSI(JLEV,ISP+1)*RLAPDI(JN)
          ZDIKHI(JLEV,ISP  ) = PDIKHI(JLEV,ISP  )*RLAPDI(JN)
          ZDIKHI(JLEV,ISP+1) = PDIKHI(JLEV,ISP+1)*RLAPDI(JN)
        ENDDO
      ENDDO
    ENDIF
  ENDDO

ELSE
  WRITE(NULOUT,*) ' PB KINP IN SPEUV ',KINP
  CALL ABOR1(' SPEUV : INVALID KINP ')
ENDIF

CALL INV_TRANS(PSPVOR=ZETPSI(1:ILF,:),PSPDIV=ZDIKHI(1:ILF,:),&
 & KVSETUV=IVSETUV(1:IGF),KRESOL=NRESOL,KPROMA=NGPTOT,PGP=ZUV)  

PU(:,:) = ZUV(:,1:KFIELD,1)
PV(:,:) = ZUV(:,KFIELD+1:2*KFIELD,1)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPEUV',1,ZHOOK_HANDLE)
END SUBROUTINE SPEUV

