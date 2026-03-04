! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE LAITQM1C(YDDIMV, KFLEV, PVINTW, KLEV, PXSL, PXF )


#ifdef DOC
!**** *LAITQM1C  -  semi-LAgrangian scheme:
!                   Vertical cubic quasi-monotone interpolation.

!     Purpose.
!     --------
!       Vertical cubic quasi-monotone interpolation.

!**   Interface.
!     ----------
!        *CALL* *LAITQM1C

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KFLEV   - vertical dimension.
!          PVINTW  - weights for cubic vertical interpolation
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
!        Joao Teixeira  *ECMWF*

!     Modifications.
!     --------------
!        Original   04-Jun-1995
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 

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
REAL(KIND=JPRB) :: PVINTW(YDDIMV%NFLEVG,2:4)
REAL(KIND=JPRB) :: PXSL(YDDIMV%NFLEVG)

!     OUTPUT:
REAL(KIND=JPRB) :: PXF(YDDIMV%NFLEVG)

!     LOCAL INTEGER SCALARS
INTEGER(KIND=JPIM) :: ILEVM2, ILEVV, JLEV

!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: ZFMAX, ZFMIN, ZII, ZIS, ZSI, ZSS


!     ------------------------------------------------------------------

!*       1.    INTERPOLATIONS.
!              ---------------

ILEVM2=KFLEV-2

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

  PXF(JLEV)=      ZSS &
   &+(ZSI-ZSS)*PVINTW(JLEV,2)&
   &+(ZIS-ZSS)*PVINTW(JLEV,3)&
   &+(ZII-ZSS)*PVINTW(JLEV,4)

  ZFMAX=MAX(ZSS,ZSI,ZIS,ZII)
  ZFMIN=MIN(ZSS,ZSI,ZIS,ZII)
  PXF(JLEV)=MAX(ZFMIN,MIN(ZFMAX,PXF(JLEV)))

ENDDO


END SUBROUTINE LAITQM1C
