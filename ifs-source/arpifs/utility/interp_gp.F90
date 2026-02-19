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

SUBROUTINE INTERP_GP(KINTERP,KPARAM,KLEVEL,KDGLG_I,KLOENG_I,PLATIG_I,&
 & KDGLG_O,KLOENG_O,PLATIG_O,PVAH,PVBH,PGLFLDI,PGLFLDO,CDLEVTYPE)

!**** *INTERP_GP* - Interpolate grid-point fields

!     Purpose.
!     --------
!     Wrapper for different interpolations in grid-point space applied
!     to non-distributed fields

!**   Interface.
!     ----------
!        *CALL* *INTERP_GP*(...)

! Explicit arguments : KINTERP - interpolation type
! -------------------- KPARAM  - GRIB code of field to be interpolated
!                      KLEVEL  - level of field
!                      KDGLG_I - number of Gaussian lats, input field
!                      KLOENG_I - number of long. points, input field
!                      PLATIG_I - Gaussian latitudes, input field
!                      KDGLG_O - number of Gaussian lats, output field
!                      KLOENG_O - number of long. points, output field
!                      PLATIG_O - Gaussian latitudes, output field
!                      PVAH     - A values for vertical grid 
!                      PVBH     - B values for vertical grid 
!                      PGLFLDI  - input field (full field)
!                      PGLFLDO  - output field (full field)
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!       Original : 15-Dec-2005
!       E.Holm    04-Jan-2008 Conserving interpolation interface
!                             changed to use surface pressure
!       P. de Rosnay and E. Holm December 2017: account for CDLEVTYPE 
!                     for conservative interpolation of surface fields
!       G. Arduini Jun 2021: add possibility of using closest interp 
!                            for snow fields if LSNOWTRAJCONS=false
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_GRIB_CODES  , ONLY : NGRBV,NGRBU,NGRBSD,NGRBRSN,NGRBTSN,NGRBWSN,NGRBASN
USE YOM_GRID_BICONSERV, ONLY : RGPPRS_HR, RGPPRS_LR
USE YOMLUN    ,ONLY : NULOUT
USE YOMTRAJ  , ONLY : LSNOWTRAJCONS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KINTERP
INTEGER(KIND=JPIM),INTENT(IN) :: KPARAM
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVEL
INTEGER(KIND=JPIM),INTENT(IN) :: KDGLG_I
INTEGER(KIND=JPIM),INTENT(IN) :: KLOENG_I(:)
REAL(KIND=JPRB)   ,INTENT(IN) :: PLATIG_I(:)
INTEGER(KIND=JPIM),INTENT(IN) :: KDGLG_O
INTEGER(KIND=JPIM),INTENT(IN) :: KLOENG_O(:)
REAL(KIND=JPRB)   ,INTENT(IN) :: PLATIG_O(:)
REAL(KIND=JPRB)   ,INTENT(IN) :: PVAH(:)
REAL(KIND=JPRB)   ,INTENT(IN) :: PVBH(:)
REAL(KIND=JPRB)   ,INTENT(IN) :: PGLFLDI(:)
REAL(KIND=JPRB)   ,INTENT(OUT):: PGLFLDO(:)
CHARACTER(LEN=*)  ,OPTIONAL,INTENT(IN) :: CDLEVTYPE

INTEGER(KIND=JPIM) :: IGPTOTG_I,IGPTOTG_O
LOGICAL :: LLPARITY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
CHARACTER(LEN=3)   :: CLEVTYPE = 'SFC' 

#include "abor1.intfb.h"
#include "grid_bicubic.intfb.h"
#include "grid_bilinear.intfb.h"
#include "grid_closest.intfb.h"
#include "grid_biconserv.intfb.h"
#include "grid_fce.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTERP_GP',0,ZHOOK_HANDLE)

!Preset defaults
CLEVTYPE = 'ML'
WRITE(NULOUT,*) 'default: ',CLEVTYPE
! Set logical switches according to arguments
IF(PRESENT(CDLEVTYPE))   CLEVTYPE    = CDLEVTYPE
WRITE(NULOUT,*) 'effective: ',CLEVTYPE
!
IGPTOTG_I = SUM(KLOENG_I(1:KDGLG_I))
IGPTOTG_O = SUM(KLOENG_O(1:KDGLG_O))

SELECT CASE (KINTERP)
CASE(0)
  WRITE(NULOUT,*) 'CASE0: ',KPARAM,KLEVEL,CLEVTYPE
  CALL GRID_CLOSEST(KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
   & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
   & PGLFLDI,PGLFLDO)
CASE(1)
  WRITE(NULOUT,*) 'CASE1: ',KPARAM,KLEVEL,CLEVTYPE
  CALL GRID_BILINEAR(KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
   & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
   & PGLFLDI,PGLFLDO)
CASE(2)
  WRITE(NULOUT,*) 'CASE2: ',KPARAM,KLEVEL,CLEVTYPE
  CALL GRID_BICUBIC(KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
   & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
   & PGLFLDI,PGLFLDO)
CASE(3)
  IF (KLEVEL==0) THEN
      WRITE(NULOUT,*) 'CASE3 interpol',KPARAM,KLEVEL,CLEVTYPE
      CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
       & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
       & PGLFLDI,PGLFLDO)
  ELSEIF (KLEVEL>0 .AND. CLEVTYPE=='SFC') THEN
     IF (.NOT.LSNOWTRAJCONS.AND.&
        &(KPARAM==NGRBSD.OR.KPARAM==NGRBRSN.OR.KPARAM==NGRBTSN.OR.&
        & KPARAM==NGRBWSN.OR.KPARAM==NGRBASN)) THEN 
       WRITE(NULOUT,*) 'CASE3 interpol with closest ',KPARAM,KLEVEL,CLEVTYPE
       CALL GRID_CLOSEST(KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
        & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
        & PGLFLDI,PGLFLDO)
     ELSE
       WRITE(NULOUT,*) 'CASE3 interpol with conservative interp ',KPARAM,KLEVEL,CLEVTYPE
       CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
        & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
        & PGLFLDI,PGLFLDO)
     ENDIF
  ELSEIF (KLEVEL>0) THEN
    WRITE(NULOUT,*) 'CASE3 klevl>0 ? ',KPARAM,KLEVEL,CLEVTYPE
    IF (.NOT.ALLOCATED(RGPPRS_LR)) THEN 
      IF (.NOT.ALLOCATED(RGPPRS_HR)) THEN
        CALL ABOR1('INTER_GP:RGPPRS_HR NOT ALLOCATED')
      ENDIF
      IF (IGPTOTG_I>IGPTOTG_O) THEN
        ALLOCATE (RGPPRS_LR(IGPTOTG_O))
        CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
                          & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),&
                          & IGPTOTG_O,&
                          & RGPPRS_HR,RGPPRS_LR)
      ELSEIF (IGPTOTG_I<=IGPTOTG_O) THEN 
        ALLOCATE (RGPPRS_LR(IGPTOTG_I))
        CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),&
                          & IGPTOTG_O,&
                          & KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
                          & RGPPRS_HR,RGPPRS_LR)
      ENDIF 
    ENDIF
    LLPARITY=.NOT.(KPARAM==NGRBV.OR.KPARAM==NGRBU )
    IF (IGPTOTG_I>IGPTOTG_O) THEN
      CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
       & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
       & PGLFLDI,PGLFLDO,LLPARITY,KLEVEL,RGPPRS_HR,RGPPRS_LR)
    ELSEIF (IGPTOTG_I<=IGPTOTG_O) THEN 
      CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
       & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
       & PGLFLDI,PGLFLDO,LLPARITY,KLEVEL,RGPPRS_LR,RGPPRS_HR)
    ENDIF
  ENDIF 
CASE(4)
  WRITE(NULOUT,*) 'CASE4: ',KPARAM,KLEVEL,CLEVTYPE
  CALL GRID_BICONSERV(PVAH,PVBH, KDGLG_I,KLOENG_I,PLATIG_I,IGPTOTG_I,&
   & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
   & PGLFLDI,PGLFLDO)
CASE(5)
  WRITE(NULOUT,*) 'CASE5: ',KPARAM,KLEVEL,CLEVTYPE
  CALL GRID_FCE(KDGLG_I,KLOENG_I(1:KDGLG_I),PLATIG_I(1:KDGLG_I),IGPTOTG_I,&
   & KDGLG_O,KLOENG_O(1:KDGLG_O),PLATIG_O(1:KDGLG_O),IGPTOTG_O,&
   & PGLFLDI(1:IGPTOTG_I),PGLFLDO(1:IGPTOTG_O))
END SELECT
IF (LHOOK) CALL DR_HOOK('INTERP_GP',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE INTERP_GP
