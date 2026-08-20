! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LAITBS1C(YDDIMV, KFLEV, PVINTW, KLEV, PRDETAR, PXSL, PXF )


#ifdef DOC
!**** *LAITBS1C  -  semi-LAgrangian scheme:
!                   Vertical cubic interpolation with conservative
!                   Bermejo-Staniforth fixer.

!     Purpose.
!     --------
!       Vertical cubic interpolation with conservative Bermejo-Staniforth fixer.

!**   Interface.
!     ----------
!        *CALL* *LAITBS1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          PVINTW  - weights for cubic vertical interpolation
!          PRDETAR - attribute VRDETAR
!          PXSL    - semi-lagrangian variable.

!        OUTPUT:
!          PXF     - interpolated variable.

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        No external.
!        Called by LARTQM1C.

!     Reference.
!     ----------

!     Author.
!     -------
!        Filip Vana (ECMWF)
!           after laitqm1c

!     Modifications.
!     --------------
!        Original   18-Jan-2023

!     ------------------------------------------------------------------
#endif

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB

IMPLICIT NONE

!     DUMMY INTEGER SCALARS
TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
INTEGER(KIND=JPIM) :: KFLEV
INTEGER(KIND=JPIM) :: KLEV(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PDVER(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PRDETAR(YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: PVINTW(YDDIMV%NFLEVG,2:4)
REAL(KIND=JPRB) :: PXSL(YDDIMV%NFLEVG)

!     OUTPUT:
REAL(KIND=JPRB) :: PXF(YDDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ILEVM2, ILEVV, JLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZFMAX, ZFMIN, ZII, ZIS, ZSI, ZSS, ZSURPL, ZXF


!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

ILEVM2=KFLEV-2

ZSURPL=0._JPRB  ! surplus, to be redistributed...

DO JLEV=1,KFLEV

  ILEVV=KLEV(JLEV)

  IF (ILEVV == 0) THEN
    ZSS=PXSL(ILEVV+1)
  ELSE
    ZSS=PXSL(ILEVV)
  ENDIF

  ZSI=PXSL(ILEVV+1)
  ZIS=PXSL(ILEVV+2)

  IF (ILEVV == ILEVM2) THEN
    ZII=PXSL(ILEVV+2)
  ELSE
    ZII=PXSL(ILEVV+3)
  ENDIF

  ZFMAX=MAX(ZSI,ZIS)
  ZFMIN=MIN(ZSI,ZIS)

  ZXF=    ZSS &
   &+(ZSI-ZSS)*PVINTW(JLEV,2)&
   &+(ZIS-ZSS)*PVINTW(JLEV,3)&
   &+(ZII-ZSS)*PVINTW(JLEV,4)

  ZXF=ZXF+ZSURPL  ! conservation step 1

  PXF(JLEV)=MAX(ZFMIN,MIN(ZFMAX,ZXF))

  ZSURPL=(ZXF-PXF(JLEV))*PRDETAR(JLEV) ! conservation step 2

ENDDO


END SUBROUTINE LAITBS1C
