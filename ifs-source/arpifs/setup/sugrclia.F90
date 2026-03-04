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

!OPTION! -O nomove
SUBROUTINE SUGRCLIA(YDGEOMETRY,YDSURF,YDMCC,YDRIP,YDML_LBC,LDINT,KLI,KINITMONTH)

!**** *SUGRCLIA*  - Initialize the clim. gridpoint fields from *FA*

!     Purpose.
!     --------
!           Initialize the clim. gridpoint fields from *FA*

!**   Interface.
!     ----------
!        *CALL* *SUGRCLIA(...)*

!        Explicit arguments :
!        ------------------
!         LDINT   - .TRUE. to make time interpolated fields
!         KLI     - Number of climatology files actually read
!         KINITMONTH : month of the initial file ; used to control the climatology, if provided in the interface

!        Implicit arguments :
!        --------------------
!        See 'USE MODULE' above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        Calls 'FA' and LFI routines.
!        Is called by SUGRIDF

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!        R. EL KHATIB *METEO-FRANCE*
!        after subroutines SUGRCLIADM and SUCACLIA 

!     Modifications.
!     --------------
!        ORIGINAL : 00-02-16
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        F.Taillefer   20-Jan-2009 Use in CANARI (instead of SUCACLIA)
!        Apr 2008  K. Yessad: use DISGRID instead of DISGRID_C + cleanings
!        R. El Khatib : 23-Apr-2010 use disgrid_mod instead of disgrid
!        P. Marguinaud: 26-May-2010 fix bug in data distribution (NPROC==1)
!        P. Marguinaud: 09-Sep-2012 Refactor using IOFLDDESC_MOD and IOGRCLIA_MOD
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!        K. Yessad (July 2014): Move some variables.
!        P. Marguinaud: 10-Oct-2014 Cleaning
!     ------------------------------------------------------------------

USE YOMMCC             , ONLY : TMCC
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMRIP             , ONLY : TRIP
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMMP0             , ONLY : MYPROC
USE YOMRIP0            , ONLY : NINDAT
USE QACLIM             , ONLY : LCLIM
USE MPL_MODULE         , ONLY : MPL_BROADCAST

USE IOGRCLIA_MOD       , ONLY :&
                               & NIOGRCLIACT_READ,   &
                               & IOGRCLIA_COUNT,     &
                               & IOGRCLIA_SELECTD,   &
                               & IOGRCLIA_SELECTF
USE YEMLBC_MODEL         , ONLY : TELBC_MODEL
USE IOFLDDESC_MOD      , ONLY : IOFLDDESC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(TMCC)        ,INTENT(INOUT) :: YDMCC
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
TYPE(TELBC_MODEL) ,INTENT(IN)    :: YDML_LBC
LOGICAL           ,INTENT(IN)    :: LDINT 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLI 
INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: KINITMONTH

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: IREP, J
INTEGER(KIND=JPIM) :: IDIM                      ! Number of files to read
LOGICAL            :: LLFILE(2)                 ! .TRUE. if climatology file fits the model date
REAL(KIND=JPRB)    :: ZALFA                     ! coefficient for time interpolation

REAL (KIND=JPRB), ALLOCATABLE, TARGET :: ZGPBUFL (:,:,:)
TYPE (IOFLDDESC), ALLOCATABLE :: YLFLDSC (:) 
REAL (KIND=JPRB), POINTER     :: ZGPBUFL12 (:,:)
INTEGER (KIND=JPIM)           :: IFNUM

INTEGER(KIND=JPIM) :: IFILE (2)     ! Identificator of files to read
INTEGER(KIND=JPIM) :: IUNTIN (2)    ! unit number of files to read
INTEGER(KIND=JPIM) :: IDATEF (11,2) ! date of files to read
INTEGER(KIND=JPIM) :: IMMCLI (2)    ! month of the files

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "openfa.intfb.h"
#include "sutimincli.intfb.h"
#include "rdfa2gp.intfb.h"
#include "fcttim.func.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRCLIA',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

KLI = 0
LCLIM = .FALSE.

IFNUM = 0
CALL IOGRCLIA_COUNT(YDGEOMETRY,YDSURF,YDMCC,NIOGRCLIACT_READ,IFNUM)

IF (IFNUM > 0) THEN

  IF (LDINT) THEN
    IDIM = 2
    IFILE (1) = 17_JPIM
    IFILE (2) = 18_JPIM
    DO J = 1, IDIM
      CALL OPENFA (YDGEOMETRY, YDRIP, IFILE(J), IUNTIN(J), &
                 & KDATE=IDATEF(:,J), YDML_LBC=YDML_LBC)
    ENDDO
  ELSE
    IDIM = 1
    IFILE = 10
    CALL OPENFA (YDGEOMETRY, YDRIP, IFILE(1), IUNTIN(1), YDML_LBC=YDML_LBC, &
               & KDATE=IDATEF(:,1),KINITMONTH=KINITMONTH)
  ENDIF

  IF (MYPROC == 1) THEN

! Compute ZALFA

    DO J = 1, IDIM
      IF (IUNTIN(J) > 0) THEN
! File exists, set the month
        IMMCLI(J) = IDATEF(2,J)
      ELSE
! Set an impossible month
        IMMCLI(J) = -999
      ENDIF
    ENDDO

    IF (LDINT) THEN 
      CALL SUTIMINCLI (NDD(NINDAT), NMM(NINDAT), NCCAA(NINDAT),&
                     & IMMCLI, NULOUT, LLFILE, ZALFA)
    ELSE
      LLFILE(:) = .TRUE.
    ENDIF
    
  ENDIF

! Broadcast parameters computed by #1

  CALL MPL_BROADCAST (IDIM, KTAG=0_JPIM, KROOT=1_JPIM, CDSTRING='SUGRCLIA:')
  CALL MPL_BROADCAST (LLFILE, KTAG=0_JPIM, KROOT=1_JPIM, CDSTRING='SUGRCLIA:')
  CALL MPL_BROADCAST (ZALFA, KTAG=0_JPIM, KROOT=1_JPIM, CDSTRING='SUGRCLIA:')

  KLI = COUNT (LLFILE)

! Get fields ids

  ALLOCATE (YLFLDSC (IFNUM), ZGPBUFL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM, IDIM))

  CALL IOGRCLIA_SELECTD(YDGEOMETRY,YDSURF,YDMCC,NIOGRCLIACT_READ,YLFLDSC)

! Read & distribute
  DO J = 1, IDIM
    IF (LLFILE (J)) THEN
      CALL RDFA2GP (YDGEOMETRY, YDRIP, IFNUM, ZGPBUFL(:,:,J), &
                  & YLFLDSC, KUNIT=IUNTIN(J), YDML_LBC=YDML_LBC)
    ENDIF
  ENDDO

! Average local buffers with ZALPHA; result is ZGPBUFL12

  LCLIM = .TRUE.

  IF (ALL (LLFILE) .AND. LDINT) THEN
    ZGPBUFL (:,:,1) = ZALFA * ZGPBUFL (:,:,1) + (1.0_JPRB-ZALFA) * ZGPBUFL (:,:,2)
    ZGPBUFL12 => ZGPBUFL (:,:,1)
  ELSEIF (LLFILE(1)) THEN
    ZGPBUFL12 => ZGPBUFL (:,:,1)
  ELSEIF (LLFILE(2)) THEN
    ZGPBUFL12 => ZGPBUFL (:,:,2)
  ELSE
    ZGPBUFL (:,:,1) = 0._JPRB
    ZGPBUFL12 => ZGPBUFL (:,:,1)
    LCLIM = .FALSE.
  ENDIF

  CALL IOGRCLIA_SELECTF(YDGEOMETRY,YDSURF,YDMCC,NIOGRCLIACT_READ,ZGPBUFL12,YLFLDSC)

  DEALLOCATE (YLFLDSC, ZGPBUFL)

  DO J = 1, IDIM
    IF (LLFILE (J)) THEN
      CALL FAIRME (IREP, IUNTIN(J), 'KEEP')
    ENDIF
  ENDDO

ENDIF


IF (LHOOK) CALL DR_HOOK('SUGRCLIA',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRCLIA

