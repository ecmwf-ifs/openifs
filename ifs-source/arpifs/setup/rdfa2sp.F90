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

SUBROUTINE RDFA2SP(YDGEOMETRY,YDRIP,KFIELDG,KFIELDS,PSPBUFL,YDFLDSC,KUNIT,KFILE,YDML_LBC)

!**** *RDFA2SP*  - Read & distribute spectral fields

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014
!
!   Modifications.
!   --------------
!                 19th July 2018  Y. Michel: specify KRESOL in EDIST_SPEC
!     ------------------------------------------------------------------

USE YOMRIP       , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMCT0       , ONLY : LELAM
USE YOMMP0       , ONLY : NSTRIN, MYPROC, MYSETV, NPROC,&
                        & NPRINTLEV

USE YEMLBC_MODEL , ONLY : TELBC_MODEL
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE FA_MOD       , ONLY : FA_COM, FA_COM_DEFAULT
USE FADUP_MOD    , ONLY : FADUPN1, FADUPN2, FADUPN3, FADUPN4,&
                        & FADUP_PARAMS
USE MFIOOPTS_MOD , ONLY : MFIOOPTS_GETOPTS, MFIOOPTS
USE OML_MOD      , ONLY : OML_MAX_THREADS, OML_MY_THREAD

!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY)    , INTENT (IN)    :: YDGEOMETRY
TYPE (TRIP)        , INTENT (IN)    :: YDRIP
INTEGER (KIND=JPIM), INTENT (IN)    :: KFIELDG
INTEGER (KIND=JPIM), INTENT (IN)    :: KFIELDS
REAL (KIND=JPRB)   , INTENT (OUT)   :: PSPBUFL(YDGEOMETRY%YRDIM%NSPEC2,KFIELDS) 
TYPE (IOFLDDESC)   , INTENT (IN)    :: YDFLDSC (KFIELDG)
INTEGER (KIND=JPIM), INTENT (IN),  OPTIONAL :: KUNIT           ! Used if present
INTEGER (KIND=JPIM), INTENT (IN),  OPTIONAL :: KFILE           ! Used if present and KUNIT not present
TYPE (TELBC_MODEL) , INTENT (IN),  OPTIONAL :: YDML_LBC

!     -------------------------------------------------------------------------

!     ISTRIN : maximum number of I/O procs on a V-set
!     IFLDSPL: total number of fields to read on this processor
!     IFLDSCH: number of fields in each V-set
!     INFD   : number of fields in the current chunk for each I/O processor
!              (local variable with respect to the V-set)
!     IFLDOFF: fields offset after distribution among ISTRIN procs.
!              (local variable with respect to the V-set)
!     IFLDG  : field index in global arrays
!     ZSPBUFG: a chunk of of spectral buffer fields (global)

INTEGER (KIND=JPIM) :: JFLDL, JFLDL1, JFLDL2
INTEGER (KIND=JPIM) :: JFLDG, JFLDG1, JFLDG2
INTEGER (KIND=JPIM) :: IFLDSPL
INTEGER (KIND=JPIM) :: IPROC
INTEGER (KIND=JPIM) :: IFLDG, ISTRIN, IOLEV
INTEGER (KIND=JPIM) :: IREP, IUNIT
INTEGER (KIND=JPIM) :: ILNOMA
INTEGER (KIND=JPIM) :: ITID, JTID
INTEGER (KIND=JPIM) :: ILONG
INTEGER (KIND=JPIM) :: IPOSG (KFIELDG), ISORTG (KFIELDG)
INTEGER (KIND=JPIM) :: IPOSS (KFIELDS), ISORTS (KFIELDS)

TYPE (IOFLDDESC)    :: YLFLDSC (KFIELDG)
INTEGER (KIND=JPIM) :: INFD (NPROC), IFLDOFF (NPROC)
LOGICAL             :: LLOPENMP
CHARACTER (LEN=16)  :: CLNOMA
INTEGER (KIND=JPIM) :: IFROM (KFIELDG)
REAL(KIND=JPRB),     ALLOCATABLE :: ZSPBUFG (:,:)
TYPE (FA_COM),       POINTER :: YLFA
TYPE (FADUP_PARAMS), POINTER :: YLDFP (:)
TYPE (MFIOOPTS)    :: YLOPTS
LOGICAL            :: LLSORT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

#include "openfa.intfb.h"
#include "sumpioh.intfb.h"
#include "dist_spec.h"
#include "edist_spec.h"
#include "facilo_mt.h"
#include "abor1.intfb.h"

!     -------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('RDFA2SP',0,ZHOOK_HANDLE)

CALL MFIOOPTS_GETOPTS (YLOPTS)

IF (PRESENT (KUNIT)) THEN
  IUNIT = KUNIT
ELSEIF (PRESENT (KFILE)) THEN
  CALL OPENFA(YDGEOMETRY,YDRIP,KFILE,IUNIT,YDML_LBC)
ELSE
  CALL ABOR1 ('RDFA2SP: KUNIT OR KFILE IS EXPECTED')
ENDIF


LLSORT = YLOPTS%LROOO

ISTRIN=MIN(NSTRIN,NPROC)

CALL SUMPIOH (NPROC, ISTRIN, KFIELDG, INFD, IFLDOFF)

IFLDSPL = INFD (MYPROC)

ALLOCATE (ZSPBUFG (YDGEOMETRY%YRDIM%NSPEC2G, IFLDSPL))

IF (LLSORT) THEN

! Read FA index

  DO JFLDG = 1, KFIELDG
    CALL FANFAN (IREP, IUNIT, YDFLDSC (JFLDG)%CPREF, YDFLDSC (JFLDG)%IOLEV,&
               & YDFLDSC (JFLDG)%CSUFF, CLNOMA, ILNOMA)
    CALL LFINFO (IREP, IUNIT, CLNOMA, ILONG, IPOSG (JFLDG))
  ENDDO

  IPOSS = PACK (IPOSG, MASK=YDFLDSC%IVSET == MYSETV)

! Sort fields by position in file

  CALL QSORTI4 (KFIELDG, ISORTG, IPOSG)
  CALL QSORTI4 (KFIELDS, ISORTS, IPOSS)

ELSE

  ISORTG = (/ (JFLDG, JFLDG = 1, KFIELDG) /)
  ISORTS = (/ (JFLDL, JFLDL = 1, KFIELDS) /)

ENDIF

YLFLDSC = YDFLDSC (ISORTG)

IF (IFLDSPL > 0) THEN

  LLOPENMP = (IFLDSPL > 2) .AND. YLOPTS%LFADPT_OPENMP

  IF (LLOPENMP) THEN
    CALL FADUPN1 (YLDFP, IUNIT)
    IUNIT = - IUNIT
  ENDIF

! Read global fields

!$OMP PARALLEL PRIVATE (YLFA, JFLDL, IFLDG, IOLEV, IREP, CLNOMA, ILNOMA, &
!$OMP&                  ITID, JTID, JFLDL1, JFLDL2)                      &
!$OMP& IF (LLOPENMP)

  IF (LLOPENMP) THEN
    CALL FADUPN2 (YLFA, YLDFP, IUNIT)
    JTID = OML_MAX_THREADS ()
    ITID = OML_MY_THREAD ()
  ELSE
    YLFA => FA_COM_DEFAULT
    JTID = 1
    ITID = 1
  ENDIF

  JFLDL1 = 1 + ((ITID-1) * IFLDSPL) / JTID
  JFLDL2 =     ((ITID  ) * IFLDSPL) / JTID

  DO JFLDL = JFLDL1, JFLDL2

    IFLDG = IFLDOFF(MYPROC) + JFLDL
    IOLEV = YLFLDSC(IFLDG)%IOLEV

    IF (NPRINTLEV > 0) THEN
      CALL FANFAN_MT (YLFA, IREP, IUNIT, YLFLDSC(IFLDG)%CPREF, IOLEV,&
                    & YLFLDSC(IFLDG)%CSUFF, CLNOMA, ILNOMA)
      WRITE (NULOUT, '("RDFA2SP:  ",I4," ",A," READ FROM UNIT = ",I3)') JFLDL, CLNOMA, IUNIT
    ENDIF

    CALL FACILO_MT (YLFA, IREP, IUNIT, YLFLDSC(IFLDG)%CPREF, IOLEV,&
                  & YLFLDSC(IFLDG)%CSUFF, ZSPBUFG(1,JFLDL), .TRUE.) 

  ENDDO

  IF (LLOPENMP) THEN
    CALL FADUPN3 (YLFA, YLDFP, IUNIT)
  ENDIF

!$OMP END PARALLEL

  IF (LLOPENMP) THEN
    CALL FADUPN4 (YLDFP)
    IUNIT = - IUNIT
  ENDIF

ENDIF

CALL FAIRME (IREP, IUNIT, 'KEEP')

! Distribute fields

DO IPROC = 1, NPROC
  JFLDG1 = IFLDOFF (IPROC) + 1
  JFLDG2 = IFLDOFF (IPROC) + INFD (IPROC)
  IFROM (JFLDG1:JFLDG2) = IPROC
ENDDO

IF (LELAM) THEN
  CALL EDIST_SPEC (PSPECG=ZSPBUFG, KFDISTG=KFIELDG, KFROM=IFROM, KSORT=ISORTS,&
       & KVSET=YLFLDSC%IVSET, PSPEC=PSPBUFL, LDIM1_IS_FLD=.FALSE.,&
       & KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ELSE
  CALL DIST_SPEC (PSPECG=ZSPBUFG, KFDISTG=KFIELDG, KFROM=IFROM, KSORT=ISORTS,&
       & KVSET=YLFLDSC%IVSET, PSPEC=PSPBUFL, LDIM1_IS_FLD=.FALSE.,&
       & KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ENDIF

DEALLOCATE (ZSPBUFG)

IF (LHOOK) CALL DR_HOOK('RDFA2SP',1,ZHOOK_HANDLE)

END SUBROUTINE RDFA2SP

