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

MODULE YOM4FPOS

USE PARKIND1  ,ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

! === DYNAMIC POST-PROCESSING ===

TYPE TRQ3FP

! Number of fields
INTEGER(KIND=JPIM) :: NPPFIELDG = 0
! Maximum number of levels (0 may stand for surface fields ?? will see later)
INTEGER(KIND=JPIM) :: NPXLEV = 0
! Spectral fit to be perfomed on source geometry (1), target geometry (2), or not (0). What about 2D fields ?
! Does "2" implies spectral output ?
INTEGER(KIND=JPIM) :: NFIT = 0
! Kind of generic spectral filter : 
! no filter (1), gaussian filter (2) or low-pass filter on homogenous resolution spectrum for for stretched geometry (3)
INTEGER(KIND=JPIM) :: NSPFIL = 1
! Possible level values
REAL(KIND=JPRB), ALLOCATABLE :: ZPXLEV(:)
! Fields internal codes
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOD(:)
! actual number of levels for each field
INTEGER(KIND=JPIM), ALLOCATABLE :: NPPLEVG(:)
! level values for each field
INTEGER(KIND=JPIM), ALLOCATABLE :: ILEV(:,:)
! number of domains for each level of each field
INTEGER(KIND=JPIM), ALLOCATABLE :: IDOM(:,:)
! domain indexes for  each level of each field
INTEGER(KIND=JPIM), ALLOCATABLE :: IDMP(:,:,:)

END TYPE TRQ3FP

TYPE TRQ2FP

INTEGER(KIND=JPIM) :: NPPFIELDG = 0
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOD(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDOM(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IDMP(:,:)

END TYPE TRQ2FP


TYPE TRQFPDYN

! Pressure levels
TYPE(TRQ3FP) :: YPR
! isentropic (Theta) levels 
TYPE(TRQ3FP) :: YTH
! Potential Vorticity levels
TYPE(TRQ3FP) :: YPV
! Height (above orography) levels
TYPE(TRQ3FP) :: YHO
! Eta (s) levels
TYPE(TRQ3FP) :: YES
! Flight (altitude in feet) levels
TYPE(TRQ3FP) :: YFL
! isothermic (T) levels
TYPE(TRQ3FP) :: YIT
! 2D fields
TYPE(TRQ2FP) :: Y2D

END TYPE TRQFPDYN

CONTAINS

SUBROUTINE ALLOCATE_TRQ3FP(KFIELD,KLEV,KDOM,LDPRINT,KULOUT,YD)

! To allocate and initialize 3D dynamic fields pp request

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
INTEGER(KIND=JPIM), INTENT(IN) :: KLEV
INTEGER(KIND=JPIM), INTENT(IN) :: KDOM
LOGICAL,            INTENT(IN) :: LDPRINT
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT
TYPE(TRQ3FP),       INTENT(OUT) :: YD

INTEGER(KIND=JPIM) :: JDOM, JLEV, JFLD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOM4FPOS:ALLOCATE_TRQ3FP',0,ZHOOK_HANDLE)

YD%NPXLEV=KLEV

ALLOCATE(YD%ZPXLEV(KLEV))
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%ZPXLEV    ',SIZE(YD%ZPXLEV   ),SHAPE(YD%ZPXLEV   )

YD%NPPFIELDG=KFIELD

ALLOCATE(YD%ICOD(KFIELD))
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%ICOD    ',SIZE(YD%ICOD   ),SHAPE(YD%ICOD   )

ALLOCATE(YD%NPPLEVG(KFIELD))
YD%NPPLEVG(:)=KLEV
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%NPPLEVG    ',SIZE(YD%NPPLEVG   ),SHAPE(YD%NPPLEVG   )

ALLOCATE(YD%ILEV(KLEV,KFIELD))
DO JFLD=1,KFIELD
  DO JLEV=1,KLEV
    YD%ILEV(JLEV,JFLD)=JLEV
  ENDDO
ENDDO
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%ILEV   ',SIZE(YD%ILEV   ),SHAPE(YD%ILEV   )

ALLOCATE(YD%IDOM(KLEV,KFIELD))
YD%IDOM(:,:)=KDOM
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%IDOM   ',SIZE(YD%IDOM   ),SHAPE(YD%IDOM   )

ALLOCATE(YD%IDMP(KDOM,KLEV,KFIELD))
DO JFLD=1,KFIELD
  DO JLEV=1,KLEV
    DO JDOM=1,KDOM
      YD%IDMP(JDOM,JLEV,JFLD)=JDOM
    ENDDO
  ENDDO
ENDDO
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%IDMP    ',SIZE(YD%IDMP   ),SHAPE(YD%IDMP   )

IF (LHOOK) CALL DR_HOOK('YOM4FPOS:ALLOCATE_TRQ3FP',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOCATE_TRQ3FP


SUBROUTINE ALLOCATE_TRQ2FP(KFIELD,KDOM,LDPRINT,KULOUT,YD)

! To allocate and initialize 2D dynamic fields pp request

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD
INTEGER(KIND=JPIM), INTENT(IN) :: KDOM
LOGICAL,            INTENT(IN) :: LDPRINT
INTEGER(KIND=JPIM), INTENT(IN) :: KULOUT
TYPE(TRQ2FP),       INTENT(OUT) :: YD

INTEGER(KIND=JPIM) :: JDOM, JFLD

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('YOM4FPOS:ALLOCATE_TRQ2FP',0,ZHOOK_HANDLE)

YD%NPPFIELDG=KFIELD

ALLOCATE(YD%ICOD(KFIELD))
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%ICOD    ',SIZE(YD%ICOD   ),SHAPE(YD%ICOD   )

ALLOCATE(YD%IDOM(KFIELD))
YD%IDOM(:)=KDOM
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%IDOM   ',SIZE(YD%IDOM   ),SHAPE(YD%IDOM   )

ALLOCATE(YD%IDMP(KDOM,KFIELD))
DO JFLD=1,KFIELD
  DO JDOM=1,KDOM
    YD%IDMP(JDOM,JFLD)=JDOM
  ENDDO
ENDDO
IF(LDPRINT) WRITE(KULOUT,'(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)') 'YD%IDMP    ',SIZE(YD%IDMP   ),SHAPE(YD%IDMP   )

IF (LHOOK) CALL DR_HOOK('YOM4FPOS:ALLOCATE_TRQ2FP',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOCATE_TRQ2FP
!     ------------------------------------------------------------------
END MODULE YOM4FPOS

