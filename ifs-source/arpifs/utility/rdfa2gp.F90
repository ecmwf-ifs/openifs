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

SUBROUTINE RDFA2GP (YDGEOMETRY, YDRIP, KFIELDG, PBUF, YDFLDSC, KUNIT, KFILE, YDML_LBC)

!**** *RDFA2GP*  - Read fields from one or several FA files and 
!                  distribute

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE*
!      Original : 11-09-2012
!
!     Modifications.
!     --------------
!      P. Marguinaud : 10-10-2013 : Use DIST_GRID & EDIST_GRID
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      P. Marguinaud : 10-10-2014 : Read fields in file order
!      O. Marsden      April 2017 : Add YDML_LBC argument for accessing TEFRLC and NEN2
!------------------------------------------------------------------------------

USE YOMRIP       , ONLY : TRIP
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMCT0       , ONLY : LELAM
USE YOMMP0       , ONLY : NPROC, MYPROC, NSTRIN, NPRINTLEV
USE YOMLUN       , ONLY : NULOUT, NULERR

USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE FA_MOD       , ONLY : FA_COM, FA_COM_DEFAULT
USE FADUP_MOD    , ONLY : FADUPN1, FADUPN2, FADUPN3, FADUPN4,&
 &                        FADUP_PARAMS
USE MFIOOPTS_MOD , ONLY : MFIOOPTS_GETOPTS, MFIOOPTS
USE OML_MOD      , ONLY : OML_MAX_THREADS, OML_MY_THREAD
USE YEMLBC_MODEL   , ONLY : TELBC_MODEL

!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY)    , INTENT (IN)  :: YDGEOMETRY
TYPE (TRIP)        , INTENT (IN)  :: YDRIP
INTEGER (KIND=JPIM), INTENT (IN)  :: KFIELDG 
REAL (KIND=JPRB)   , INTENT (OUT) :: PBUF (YDGEOMETRY%YRGEM%NGPTOT, KFIELDG, 1)
TYPE (IOFLDDESC)   , INTENT (IN)  :: YDFLDSC (KFIELDG) 
INTEGER (KIND=JPIM), INTENT (IN),  OPTIONAL :: KUNIT           ! Used if present
INTEGER (KIND=JPIM), INTENT (IN),  OPTIONAL :: KFILE           ! Used if present and KUNIT not present
TYPE (TELBC_MODEL),  INTENT (IN),  OPTIONAL :: YDML_LBC 

!------------------------------------------------------------------------------

#include "openfa.intfb.h"
#include "sumpioh.intfb.h"
#include "abor1.intfb.h"
#include "dist_grid.h"
#include "edist_grid.h"
#include "facilo_mt.h"

!------------------------------------------------------------------------------

INTEGER (KIND=JPIM) :: IPFLD (NPROC), IPFLDOFF (NPROC)
LOGICAL             :: LLFULL, LLCHKF

REAL (KIND=JPRB), ALLOCATABLE :: ZFIELDBUF (:,:)

TYPE (IOFLDDESC)    :: YLFLDSC (KFIELDG)
INTEGER (KIND=JPIM) :: JFLDG, JFLDL, JFLDL1, JFLDL2
INTEGER (KIND=JPIM) :: IUNIT
INTEGER (KIND=JPIM) :: IFROM (KFIELDG)
INTEGER (KIND=JPIM) :: IREP, ILONG, IPOSEX
INTEGER (KIND=JPIM) :: ILEV
INTEGER (KIND=JPIM) :: ILNOMA
INTEGER (KIND=JPIM) :: IPROC
INTEGER (KIND=JPIM) :: ITID, JTID
INTEGER (KIND=JPIM) :: IERROR
CHARACTER (LEN=16)  :: CLNOMA
LOGICAL :: LLEXIS, LLCOSP
INTEGER (KIND=JPIM) :: INGRIB, IARG1, IARG2, IARG3
TYPE (FA_COM),       POINTER :: YLFA
TYPE (FADUP_PARAMS), POINTER :: YLDFP (:)
TYPE (MFIOOPTS) :: YLOPTS
INTEGER (KIND=JPIM) :: ILEN (KFIELDG), IPOS (KFIELDG), ISORT (KFIELDG)
LOGICAL :: LLOPENMP
LOGICAL :: LLABORT
LOGICAL :: LLSORT

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK ('RDFA2GP',0,ZHOOK_HANDLE)


!------------------------------------------------------------------------------
CALL MFIOOPTS_GETOPTS (YLOPTS)

LLSORT = YLOPTS%LROOO

IF (ANY (YDFLDSC%CSUFF == ""))&
  & CALL ABOR1 ('RDFA2GP: INCOMPLETE FIELD DESCRIPTION')

CALL SUMPIOH (NPROC, NSTRIN, KFIELDG, IPFLD, IPFLDOFF)

LLFULL = .TRUE.
LLCHKF = .TRUE.

IF (PRESENT (KUNIT)) THEN
  IUNIT = KUNIT
ELSEIF (PRESENT (KFILE)) THEN
  CALL OPENFA (YDGEOMETRY, YDRIP, KFILE, IUNIT, YDML_LBC=YDML_LBC)
  SELECT CASE (KFILE)
    CASE (9) ! Lateral coupling
     LLFULL = .FALSE.
     LLCHKF = .TRUE.
    CASE (13, 14) ! SURFEX initial conditions (YRML_LBC is not ready when these files are read)
     LLFULL = .TRUE.
     LLCHKF = .FALSE.
  END SELECT
ELSE
  CALL ABOR1 ('RDFA2GP: KUNIT OR KFILE IS EXPECTED')
ENDIF

IF (LLSORT) THEN

! Read FA index

  DO JFLDG = 1, KFIELDG
    CALL FANFAN (IREP, IUNIT, YDFLDSC (JFLDG)%CPREF, YDFLDSC (JFLDG)%IOLEV,&
               & YDFLDSC (JFLDG)%CSUFF, CLNOMA, ILNOMA)
    CALL LFINFO (IREP, IUNIT, CLNOMA, ILEN (JFLDG), IPOS (JFLDG))
  ENDDO

! Sort fields by position in file

  CALL QSORTI4 (KFIELDG, ISORT, IPOS)

ELSE

  ISORT = (/ (JFLDG, JFLDG = 1, KFIELDG) /)

ENDIF

YLFLDSC = YDFLDSC (ISORT)

ALLOCATE (ZFIELDBUF (YDGEOMETRY%YRGEM%NGPTOTG, IPFLD (MYPROC)))

! IO tasks read in fields & scatter

IF (IPFLD (MYPROC) > 0) THEN

  LLOPENMP = (IPFLD (MYPROC) > 1) .AND. YLOPTS%LFADPT_OPENMP

  IF (LLOPENMP) THEN
    CALL FADUPN1 (YLDFP, IUNIT)
    IUNIT = - IUNIT
  ENDIF

  IERROR = 0

!$OMP PARALLEL PRIVATE (YLFA, JFLDL, JFLDG, LLABORT, ILEV,   &
!$OMP&                  CLNOMA, ILNOMA, IREP, ILONG, IPOSEX, &
!$OMP&                  ITID, JTID, JFLDL1, JFLDL2,          &
!$OMP&                  INGRIB, IARG1, IARG2, IARG3,         &
!$OMP&                  LLEXIS, LLCOSP)                      &
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

  JFLDL1 = 1 + ((ITID-1) * IPFLD (MYPROC)) / JTID
  JFLDL2 =     ((ITID  ) * IPFLD (MYPROC)) / JTID

  DO JFLDL = JFLDL1, JFLDL2

    JFLDG = IPFLDOFF (MYPROC) + JFLDL

    LLABORT = YLFLDSC (JFLDG)%LREQD

    IF (YLFLDSC(JFLDG)%LIOLV) THEN
      ILEV = YLFLDSC(JFLDG)%IOLEV
    ELSE
      ILEV = YLFLDSC(JFLDG)%ILEVG
    ENDIF

    CALL FANFAN_MT (YLFA, IREP, IUNIT, YLFLDSC(JFLDG)%CPREF, ILEV,&
                  & YLFLDSC(JFLDG)%CSUFF, CLNOMA, ILNOMA)

    IF (NPRINTLEV > 0) THEN
      WRITE (NULOUT, '("RDFA2GP:  ",I4," ",A," READ FROM UNIT = ",I3)')&
           & JFLDL, CLNOMA, IUNIT
    ENDIF

    CALL LFINFO_MT (YLFA%LFI, IREP, IUNIT, CLNOMA (1:ILNOMA), ILONG, IPOSEX)

    IF (ILONG > 0) THEN

      IF (LLCHKF) THEN
        CALL FANION_MT (YLFA, IREP, IUNIT, YLFLDSC(JFLDG)%CPREF, ILEV,&
                      & YLFLDSC(JFLDG)%CSUFF, LLEXIS, LLCOSP, INGRIB, IARG1, IARG2, IARG3)
        
        IF (INGRIB == 4) THEN
          IF (ILEV <= YDML_LBC%NEN2 .OR. LLFULL) THEN
!$OMP CRITICAL
            WRITE (NULOUT, *) 'RDFA2GP: FIELD '//TRIM (CLNOMA (1:ILNOMA))//' SHOULD NOT BE HOLLOW'
            WRITE (NULERR, *) 'RDFA2GP: FIELD '//TRIM (CLNOMA (1:ILNOMA))//' SHOULD NOT BE HOLLOW'
            IERROR = IERROR + 1
!$OMP END CRITICAL
          ENDIF
        ENDIF
      ENDIF
     
! Field was found

      CALL FACILO_MT (YLFA, IREP, IUNIT, YLFLDSC(JFLDG)%CPREF, ILEV,     &
                    & YLFLDSC(JFLDG)%CSUFF, ZFIELDBUF(:,JFLDL), .FALSE., &
                    & YLFLDSC(JFLDG)%LUNDF, YLFLDSC(JFLDG)%XUNDF)

    ELSEIF (.NOT. LLABORT) THEN

! Field was not found, but was not required

      ZFIELDBUF(:,JFLDL) = 0._JPRB

      WRITE(NULOUT,'("RDFA2GP: ",A," IS MISSING; SET TO ZERO")') CLNOMA 

    ELSE

      WRITE (NULOUT, '("RDFA2GP: ",A,'' IS MISSING; UNIT = '',I5)') CLNOMA, IUNIT
      CALL ABOR1 ('RDFA2GP: FIELD IS MISSING :'//TRIM (CLNOMA (1:ILNOMA)))

    ENDIF

  ENDDO

  IF (LLOPENMP) THEN
    CALL FADUPN3 (YLFA, YLDFP, IUNIT)
  ENDIF

!$OMP END PARALLEL

  IF (LLOPENMP) THEN
    CALL FADUPN4 (YLDFP)
    IUNIT = - IUNIT
  ENDIF

  IF (IERROR > 0) THEN
    CALL ABOR1 ('RDFA2GP: AN ERROR OCCURRED')
  ENDIF

ENDIF

IF (PRESENT (KUNIT)) THEN
ELSEIF (PRESENT (KFILE)) THEN
  CALL FAIRME (IREP, IUNIT, 'KEEP')
ENDIF

DO IPROC = 1, NPROC
  IFROM (IPFLDOFF (IPROC)+1:IPFLDOFF (IPROC)+IPFLD (IPROC)) = IPROC
ENDDO

IF (LELAM) THEN
  CALL EDIST_GRID (PGPG=ZFIELDBUF,KFDISTG=KFIELDG,KFROM=IFROM,&
     & PGP=PBUF,KSORT=ISORT,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ELSE
  CALL DIST_GRID (PGPG=ZFIELDBUF,KFDISTG=KFIELDG,KFROM=IFROM,&
     & PGP=PBUF,KSORT=ISORT,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ENDIF

DEALLOCATE (ZFIELDBUF)

!------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK ('RDFA2GP',1,ZHOOK_HANDLE)
END SUBROUTINE RDFA2GP

