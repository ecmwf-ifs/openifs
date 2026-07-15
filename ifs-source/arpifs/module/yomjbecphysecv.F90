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

MODULE YOMJBECPHYSECV

!     Purpose.
!     --------
!       Data and controls for extended control variable used
!       for the optimisation of physical paremetrisation

!     Author.
!     -------
!       S. Massart

!     Modifications.
!     --------------
!       Original    April-2020
!       S. Massart  July-2021   Addition of 2D/3D Spectral fields
! ------------------------------------------------------------------

USE PARKIND1, ONLY: JPIM, JPRB

IMPLICIT NONE

SAVE

  INTEGER(KIND=JPIM), PARAMETER :: MMISSING=-999_JPIM
  REAL(KIND=JPRB),    PARAMETER :: ZMISSING=-999._JPRB

  TYPE, PUBLIC :: TECVPHYS
    !   ------------------------------------------------------------------
    !   CNAME      Field name
    !   NINDEX     Index in the ECV 3D fields
    !   LLOGCV     True to have a log normal control variable
    !   LRESCALE   True to reascle the background error
    !   LLAND      True if over land only
    !   RBCAP      Bottom cap for first guess
    !   ------------------------------------------------------------------
    CHARACTER(LEN=16)    :: CNAME = ''
    INTEGER(KIND=JPIM)   :: NINDEX = MMISSING
    LOGICAL              :: LLOGCV = .FALSE.
    LOGICAL              :: LRESCALE = .TRUE.
    LOGICAL              :: LLAND = .FALSE.
    REAL(KIND=JPRB)      :: RBCAP = ZMISSING
  END TYPE TECVPHYS

  !   ------------------------------------------------------------------
  !   NPHYS_ECV     Number of physical parameters to optimise
  !   YPHYS_ECV     Structure for each parameter
  !   ------------------------------------------------------------------

  INTEGER(KIND=JPIM)              :: NPHYS_ECV
  INTEGER(KIND=JPIM), ALLOCATABLE :: NPHYS_ECV_VIDS(:)
  TYPE(TECVPHYS)                  :: YPHYS_ECV_SDFOR
  TYPE(TECVPHYS)                  :: YPHYS_ECV_GPPBFAS
  TYPE(TECVPHYS)                  :: YPHYS_ECV_RECBFAS


  !   ------------------------------------------------------------------
  !   LSOLARCST     True to add  the optimisation of solar constant
  !   LSIGFLTORO    True to add  the optimisation of std. dev. filtered orography
  !   LBFASECV      True for the optimisation of BFAS coeficients
  !   LBFAS_GPP     True for the optimisation of BFAS/GPP
  !   LBFAS_REC     True for the optimisation of BFAS/REC
  !   ------------------------------------------------------------------

  LOGICAL :: LSOLARCST
  LOGICAL :: LSIGFLTORO
  LOGICAL :: LBFASECV
  LOGICAL :: LBFAS_GPP
  LOGICAL :: LBFAS_REC

#include "abor1.intfb.h"
!-----------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------

SUBROUTINE GET_PHYS_ECV(YDGEOMETRY, YDECV, KID, KBLOCK, PPHYS)

!     Purpose.
!     --------
!       Get the physics parameters from the ECV structure

USE GEOMETRY_MOD,             ONLY : GEOMETRY
USE YOMHOOK,                  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR,                   ONLY : LECV
USE YOMLUN,                   ONLY : NULOUT
USE ECV_DEFINITIONS,          ONLY : ECVFIELD_ACCESS,VID
USE FIELD_CONTAINER_GP_MOD,   ONLY : FIELD_CONTAINER_GP
USE FIELD_CONTAINER_OPER_MOD, ONLY : FIELD_CONTAINER_GPNORM
USE YOMJBECV,                 ONLY : NDIAECV

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_CONTAINER_GP),INTENT(INOUT) :: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KID
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB)         ,INTENT(OUT)   :: PPHYS(YDGEOMETRY%YRDIM%NPROMA)

TYPE(ECVFIELD_ACCESS) :: YLFAC
INTEGER(KIND=JPIM)    :: JECV,  IPHYS, JINDEX, JROF
INTEGER(KIND=JPIM)    :: JKGLO, IST, ICEND, IBL, ILEV
LOGICAL               :: LLGET

REAL(KIND=JPHOOK)     :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:GET_PHYS_ECV',0,ZHOOK_HANDLE)

ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA,&
        & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT,&
        & NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS)

! 1. INITIALISATION


PPHYS(:) = 0.0_JPRB


! 2. GET PHYSICS FIELDS

IF (LECV) THEN

  DO WHILE (YDECV%FIELD_ITERATED(KFIDS=(/KID/),KBLOCK=KBLOCK,FAC=YLFAC) == 1)
    IF (.NOT. ASSOCIATED(YLFAC%RR1)) &
      &  CALL ABOR1("YOMJBECPHYSECV:GET_PHYS_ECV RR1 NOT ASSOCIATED")
    JKGLO = 1+(KBLOCK-1)*NPROMA
    IST = 1
    ICEND = MIN(NPROMA,NGPTOT-JKGLO+1)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JROF)
    DO JROF=IST,ICEND
      PPHYS(JROF) = YLFAC%RR1(JROF)
    ENDDO
!$OMP END PARALLEL DO
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:GET_PHYS_ECV',1,ZHOOK_HANDLE)

END SUBROUTINE GET_PHYS_ECV

!-----------------------------------------------------------------------

SUBROUTINE GET_PHYS_ECV_AD(YDGEOMETRY, YDECV, KID, KBLOCK, PPHYS)

!     Purpose.
!     --------
!       Get the physics parameters from the ECV structure (adjoint version)

USE GEOMETRY_MOD,           ONLY : GEOMETRY
USE YOMHOOK,                ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMVAR,                 ONLY : LECV
USE ECV_DEFINITIONS,        ONLY : ECVFIELD_ACCESS,VID
USE FIELD_CONTAINER_GP_MOD, ONLY : FIELD_CONTAINER_GP

IMPLICIT NONE

TYPE(GEOMETRY)          ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELD_CONTAINER_GP),INTENT(INOUT) :: YDECV
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KID
INTEGER(KIND=JPIM)      ,INTENT(IN)    :: KBLOCK
REAL(KIND=JPRB)         ,INTENT(INOUT) :: PPHYS(YDGEOMETRY%YRDIM%NPROMA)
call abor1("oifs/fc-only - GET_PHYS_ECV_AD should never be called")

END SUBROUTINE GET_PHYS_ECV_AD

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE ECPHYS_UPPER_AIR_GP(YDGEOMETRY,PECV,KECV)

!     Purpose.
!     --------
!       Allow to cap the min value while calling UPPER_AIR (grid point)


USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN) :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN) :: KECV
REAL(KIND=JPRB)   , INTENT(IN) :: PECV(YDGEOMETRY%YRGEM%NGPTOT,YDGEOMETRY%YRDIMV%NFLEVG)

INTEGER(KIND=JPIM):: JSTGLO, ICEND, ISTC, IBL, IPHYS, JINDEX, JROF, JK
REAL(KIND=JPRB)   :: ZBCAP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:ECPHYS_UPPER_AIR_GP',0,ZHOOK_HANDLE)

ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA, &
        & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT,&
        & NSPEC2=>YDGEOMETRY%YRDIM%NSPEC2,&
        & NFLEVL=>YDGEOMETRY%YRDIMV%NFLEVL,&
        & MYLEVS=>YDGEOMETRY%YRMP%MYLEVS)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSTGLO,ICEND,ISTC,IBL,&
!$OMP& IPHYS,JINDEX,ZBCAP,JROF)
DO JSTGLO=1,NGPTOT,NPROMA
  ICEND = MIN(NPROMA,NGPTOT-JSTGLO+1)
  ISTC  = 1
  IBL   = (JSTGLO-1)/NPROMA+1
!  DO IPHYS = 1, NPHYS_ECV
!    JINDEX = YPHYS_ECV(IPHYS)%NINDEX
!    ZBCAP = YPHYS_ECV(IPHYS)%RBCAP
!    DO JROF=ISTC,ICEND
!      YRECV5%RGPECV3D(JROF,JINDEX,IBL,KECV) = YRECV5%RGPECV3D(JROF,JINDEX,IBL,KECV) &
!      & + PECV(JSTGLO+JROF-1,JINDEX)
!      IF (ZBCAP /=  ZMISSING) THEN
!          YRECV5%RGPECV3D(JROF,JINDEX,IBL,KECV) = MAX(ZBCAP, &
!            &  YRECV5%RGPECV3D(JROF,JINDEX,IBL,KECV))
!      ENDIF
!    ENDDO
!  ENDDO
ENDDO
!$OMP END PARALLEL DO

END ASSOCIATE

IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:ECPHYS_UPPER_AIR_GP',1,ZHOOK_HANDLE)

END SUBROUTINE ECPHYS_UPPER_AIR_GP

!-----------------------------------------------------------------------

SUBROUTINE SUINFCE_ECPHYS(YDGEOMETRY,YDFIELDS,YD_JB_STRUCT)

!     Purpose.
!     --------
!       Change of background error standard deviation if log normal change of variable
!       and set it to 0 over ocean if needed

USE GEOMETRY_MOD   , ONLY : GEOMETRY
USE FIELDS_MOD     , ONLY : FIELDS
USE YOMJG          , ONLY : TYPE_JB_STRUCT
USE YOMJBECV       , ONLY : YRECVDATA
USE YOMHOOK        , ONLY : LHOOK, DR_HOOK, JPHOOK
USE ECV_DEFINITIONS, ONLY : ECVFIELD_ACCESS,VID

IMPLICIT NONE

TYPE(GEOMETRY)      ,INTENT(IN)    :: YDGEOMETRY
TYPE(FIELDS)        ,INTENT(INOUT) :: YDFIELDS
TYPE(TYPE_JB_STRUCT),INTENT(INOUT) :: YD_JB_STRUCT

TYPE(ECVFIELD_ACCESS)           :: YLFAC
REAL(KIND=JPRB),    ALLOCATABLE :: ZFCE(:)
INTEGER(KIND=JPIM)              :: JECV, JPHYS, JF, JK, JBLK
INTEGER(KIND=JPIM)              :: JKGLO, IST, ICEND
INTEGER(KIND=JPIM), ALLOCATABLE :: IFCEBALPHA(:)
LOGICAL                         :: LLCHANGE, LLOGCV, LLAND

REAL(KIND=JPHOOK)            :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:SUINFCE_ECPHYS',0,ZHOOK_HANDLE)

ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA, &
        & NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS,&
        & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, &
        & NOFCEF=>YD_JB_STRUCT%JB_DATA%NOFCEF, &
        & YDSURF=>YDFIELDS%YRSURF)

! 1. INITIALISATION

ALLOCATE(ZFCE(NPROMA))
ALLOCATE(IFCEBALPHA(NPHYS_ECV))

! 2. LOOP OVER THE ECV FIELDS

JECV = 0
JPHYS = 0
DO JF=1,YD_JB_STRUCT%CONFIG%N_SPJB_VARS
  IF (YD_JB_STRUCT%SPJB_VARS_INFO(JF)%L_IN_ECV) THEN
    JECV = JECV + 1
    IF (TRIM(YRECVDATA%CSETDESC(JECV)) ==  'EC_PHYS') THEN
      JPHYS = JPHYS + 1
      IF (JPHYS > NPHYS_ECV) CALL ABOR1('YOMJBECPHYSECV:SUINFCE_ECPHYS SOMETHING IS WRONG')
      IFCEBALPHA(JPHYS) = YD_JB_STRUCT%SPJB_VARS_INFO(JF)%IPTFCE
    ENDIF
  ENDIF
ENDDO

! 3. LOOP OVER THE ECV FIELDS TO CHANGE FCE

DO JPHYS = 1, NPHYS_ECV
  IF (NPHYS_ECV_VIDS(JPHYS) == VID%SDFOR) THEN
    LLOGCV = YPHYS_ECV_SDFOR%LLOGCV
    LLAND = YPHYS_ECV_SDFOR%LLAND
  ELSEIF (NPHYS_ECV_VIDS(JPHYS) == VID%GPPBFAS) THEN
    LLOGCV = YPHYS_ECV_GPPBFAS%LLOGCV
    LLAND = YPHYS_ECV_GPPBFAS%LLAND
  ELSEIF (NPHYS_ECV_VIDS(JPHYS) == VID%RECBFAS) THEN
    LLOGCV = YPHYS_ECV_RECBFAS%LLOGCV
    LLAND = YPHYS_ECV_RECBFAS%LLAND
  ELSE
    CALL ABOR1('YOMJBECPHYSECV:SUINFCE_ECPHYS TEST MISSING')
  ENDIF

  LLCHANGE = LLOGCV .OR. LLAND
  IF (LLCHANGE) THEN
    DO JBLK=1,NGPBLKS
      !
      ! Dimension
      !
      JKGLO = 1+(JBLK-1)*NPROMA
      IST = 1
      ICEND = MIN(NPROMA,NGPTOT-JKGLO+1)
      !
      ! Get FCE
      !
      ZFCE(IST:ICEND) = YD_JB_STRUCT%JB_DATA%FCE%BUF(IST:ICEND, IFCEBALPHA(JPHYS), JBLK)
      !
      ! Normalise with the first guess
      !
      IF (LLOGCV) THEN
        DO WHILE (YDFIELDS%FIELD5_ECV%FIELD_ITERATED(KFIDS=(/NPHYS_ECV_VIDS(JPHYS)/),KBLOCK=JBLK,FAC=YLFAC) == 1)
          IF (.NOT. ASSOCIATED(YLFAC%RR1)) &
            &  CALL ABOR1('YOMJBECPHYSECV:SUINFCE_ECPHYS RR1 NOT ASSOCIATED')
          DO JK = IST, ICEND
            IF (YLFAC%RR1(JK) > TINY(1._JPRB)) THEN
              ZFCE(JK) = ZFCE(JK)/YLFAC%RR1(JK)
            ELSE
              ZFCE(JK) = 0._JPRB
            ENDIF
           ENDDO
        ENDDO
      ENDIF
      !
      ! Add land/sea mask
      !
      IF (LLAND) THEN
        DO JK = IST, ICEND
          ZFCE(JK) = ZFCE(JK) &
            & * YDSURF%SD_VF(JK,YDSURF%YSD_VF%YLSM%MP,JBLK)
        ENDDO
      ENDIF
      !
      ! Put FCE
      !
      YD_JB_STRUCT%JB_DATA%FCE%BUF(IST:ICEND, IFCEBALPHA(JPHYS), JBLK) = ZFCE(IST:ICEND)
    ENDDO
  ENDIF
ENDDO

DEALLOCATE(ZFCE)
DEALLOCATE(IFCEBALPHA)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YOMJBECPHYSECV:SUINFCE_ECPHYS',1,ZHOOK_HANDLE)

END SUBROUTINE SUINFCE_ECPHYS

!     ------------------------------------------------------------------

END MODULE YOMJBECPHYSECV
