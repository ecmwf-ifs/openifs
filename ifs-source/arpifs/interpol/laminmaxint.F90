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

SUBROUTINE LAMINMAXINT(KPROMA,KPROMB,KST,KPROF,KFLEV, &
 &               KFLDN,KFLDX,KL0,PXSL,PMINVFLD,PMAXVFLD)

! Purpose :
! -------
!   LAMINMAXINT - FIND MINIMUM/MAXIMUM G.P. field VALUES SUROUNDING A departure point

! Interface :
! ---------
!   INPUT:
!     KPROMA  - horizontal dimension for grid-point quantities
!     KPROMB  - horizontal dimension for interpolation point quantities
!     KST     - first element of arrays where computations are performed
!     KPROF   - depth of work
!     KFLEV   - vertical dimension
!     KFLDN   - number of the first field
!     KFLDX   - number of the last field
!     KL0     - indices of the four western points of the 16 point
!               interpolation grid
!     PXSL    - quantity to be quasi-monotically limited 
!   OUTPUT:
!     PXF     - final (limited) interpolated variable

! Externals :
! ---------
!   None.

! Method :
! ------
!   See documentation.

! Reference :
! ---------

! Author :
! ------
!   -Mar-2012 M. Diamantakis

! Modifications :
! -------------

! End Modifications
!------------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB     ,JPIA
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK, JPHOOK

!------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KPROMB,KFLEV,0:3)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXSL(KPROMA*(KFLDX-KFLDN+1))
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PMINVFLD(KPROMB,KFLEV)
REAL(KIND=JPRB)   ,INTENT(OUT), OPTIONAL :: PMAXVFLD(KPROMB,KFLEV)

!------------------------------------------------------------------------------

INTEGER(KIND=JPIA) :: IV1L1, IV1L2, IV2L1, IV2L2
INTEGER(KIND=JPIM) :: JLEV, JROF

REAL(KIND=JPRB) :: ZF111, ZF112, ZF121, ZF122
REAL(KIND=JPRB) :: ZF211, ZF212, ZF221, ZF222
REAL(KIND=JPRB) :: ZMIN, ZMAX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAMINMAXINT',0,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

! Limit interpolated value by min/max of surrounding gp values

! offsets for getting values on interpolation stencil
IV1L1=1+KPROMA
IV1L2=2+KPROMA
IV2L1=IV1L1+KPROMA
IV2L2=IV1L2+KPROMA

IF (PRESENT(PMAXVFLD)) THEN
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1)
      ZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2)
      ZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1)
      ZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2)
      ZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1)
      ZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2)
      ZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1)
      ZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2)
      ZMIN=MIN(ZF111,ZF211,ZF121,ZF221,ZF112,ZF212,ZF122,ZF222)
      ZMAX=MAX(ZF111,ZF211,ZF121,ZF221,ZF112,ZF212,ZF122,ZF222)
      PMINVFLD(JROF,JLEV)=ZMIN
      PMAXVFLD(JROF,JLEV)=ZMAX
    ENDDO
  ENDDO
ELSE
  DO JLEV=1,KFLEV
    DO JROF=KST,KPROF
      ZF111=PXSL(KL0(JROF,JLEV,1)+IV1L1)
      ZF211=PXSL(KL0(JROF,JLEV,1)+IV1L2)
      ZF121=PXSL(KL0(JROF,JLEV,2)+IV1L1)
      ZF221=PXSL(KL0(JROF,JLEV,2)+IV1L2)
      ZF112=PXSL(KL0(JROF,JLEV,1)+IV2L1)
      ZF212=PXSL(KL0(JROF,JLEV,1)+IV2L2)
      ZF122=PXSL(KL0(JROF,JLEV,2)+IV2L1)
      ZF222=PXSL(KL0(JROF,JLEV,2)+IV2L2)
      ZMIN=MIN(ZF111,ZF211,ZF121,ZF221,ZF112,ZF212,ZF122,ZF222)
      PMINVFLD(JROF,JLEV)=ZMIN
    ENDDO
  ENDDO
ENDIF

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LAMINMAXINT',1,ZHOOK_HANDLE)
END SUBROUTINE LAMINMAXINT

