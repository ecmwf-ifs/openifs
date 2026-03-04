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

SUBROUTINE SUGRIDA(YDGEOMETRY,YDSURF,YDMODEL,KFILE)

!**** *SUGRIDA*  - Initialize the gridpoint fields from *FA*

!     Purpose.
!     --------
!           Initialize the gridpoint fields of the model from FA.

!**   Interface.
!     ----------
!        *CALL* *SUGRIDA(...)*

!        Explicit arguments : KFILE - indicator for which file is to be read
!        ------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.     RDGPFA - read ARPEGE field
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      MPP Group *METEO-FRANCE*
!      Original : 96-04-30 (From SUGRIDA)

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      P. Marquet   : 02-07-17 add VCLIA for aerosol files
!      M.Hamrud      01-Jul-2006 Revised surface fields
!      A. Alias     : 05-10-10 GMGEC/EAC Modifs
!      P. Marquet   : add *VCLIS for ozone forcing files 
!      M. Deque     : allows NVCLIS=1 with LOZONE=.F.
!      K. Yessad    : read surface trajectory files (LTRAJHR_SURF=T).
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      R. El Khatib : 26-Feb-2008 Bugfix
!      R. ElKhatib   24-Oct-2008 Merge YEMOP and YOMOPH, USE OPENFA
!      R. El Khatib : 23-Apr-2010 use disgrid_mod instead of disgrid
!        A. Voldoire  : Dec 2010 in LSFXLSM case, read ITM from GPARBUF
!      R. El Khatib 30-Mar-2012 fanmsg moved to the setup
!      P. Marguinaud: 11-Sep-2012 Refactor using IOGRIDA_MOD, IOFLDDESC_MOD, RDFA2GP
!                                 and SUGRIDA_FIX_TOZ
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      P. Marguinaud: 10-Oct-2014 : Use SUGRIDA_FIXUP
!      O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : NCONF
USE YOMMP0             , ONLY : MYPROC

USE IOGRIDA_MOD        , ONLY : NIOGRIDACT_READ,&
 &                              IOGRIDA_COUNT,& 
 &                              IOGRIDA_SELECTD
USE IOGRID_MOD         , ONLY : IOGRID_SELECTF,&
 &                              NIOGRIDCT_READ
USE IOFLDPTR_MOD       , ONLY : IOFLDPTR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       , INTENT(INOUT) :: YDSURF
TYPE(MODEL)       , INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN)    :: KFILE 
!     ------------------------------------------------------------------
#include "rdfa2gp.intfb.h"
#include "sugrida_fixup.intfb.h"

TYPE (IOFLDPTR),  ALLOCATABLE :: YLFLDGT (:)
REAL (KIND=JPRB), ALLOCATABLE :: ZBUFL (:,:)
INTEGER(KIND=JPIM) :: IFNUM
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGRIDA',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

IFNUM = 0
CALL IOGRIDA_COUNT(YDGEOMETRY,YDSURF,YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_PHY_G%YRDPHY,YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDMODEL%YRML_PHY_MF%YRPHYDS,NIOGRIDACT_READ,IFNUM)

IF (IFNUM > 0) THEN

  ALLOCATE (YLFLDGT (IFNUM), ZBUFL (YDGEOMETRY%YRGEM%NGPTOT, IFNUM))

  CALL IOGRIDA_SELECTD(YDGEOMETRY, YDSURF, YDMODEL%YRML_AOC%YRMCC, YDMODEL%YRML_PHY_G%YRDPHY, YDMODEL%YRML_PHY_MF%YRPHY, &
 & YDMODEL%YRML_PHY_MF%YRPHYDS, NIOGRIDACT_READ, YLFLDGT)

  IF (NCONF == 701) CALL SUGRIDA_CHECK_C701
  
! read & distribute in local buffers
  CALL RDFA2GP (YDGEOMETRY, YDMODEL%YRML_GCONF%YRRIP, IFNUM, ZBUFL, &
              & YLFLDGT%YFLDDSC, KFILE=KFILE, YDML_LBC=YDMODEL%YRML_LBC)

! fill model variables
  CALL IOGRID_SELECTF (YDGEOMETRY, NIOGRIDCT_READ, ZBUFL, YLFLDGT)
  
! fix irregular fields

  CALL SUGRIDA_FIXUP(YDGEOMETRY,YDSURF,YDMODEL,IFNUM,ZBUFL,YLFLDGT%YFLDDSC)

  DEALLOCATE (YLFLDGT, ZBUFL)

ENDIF 


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUGRIDA',1,ZHOOK_HANDLE)

CONTAINS

SUBROUTINE SUGRIDA_CHECK_C701

USE YOMDAG, ONLY: NDATGU, NRESGU, NECHGU, NTYPGU
USE YOMCT0, ONLY: LINFLAT
USE YOMMP0, ONLY: NPROC, NSTRIN

#include "abor1.intfb.h"
#include "openfa.intfb.h"
#include "sumpioh.intfb.h"

CHARACTER(LEN=16)   :: CLIDEN
INTEGER (KIND=JPIM) :: IREP
LOGICAL             :: LLISIO
INTEGER(KIND=JPIM)  :: IUNTIN
INTEGER(KIND=JPIM)  :: IDATEF(11)

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('SUGRIDA:SUGRIDA_CHECK_C701',0,ZHOOK_HANDLE)

CALL SUMPIOH (NPROC, NSTRIN, LDISIO=LLISIO, KMYPROC=MYPROC)
  
IF (LLISIO) THEN

! Open file (only for nstrin procs)
  CALL OPENFA (YDGEOMETRY, YDMODEL%YRML_GCONF%YRRIP, KFILE, IUNTIN, &
             & YDML_LBC=YDMODEL%YRML_LBC, KDATE=IDATEF)

  CALL FALSIF (IREP, IUNTIN, CLIDEN)
  
  IF (CLIDEN(1:3) == 'CLI') THEN
    NTYPGU=2
  ELSEIF (CLIDEN(1:3) == 'ARP') THEN
    NTYPGU=3
  ELSEIF (CLIDEN(1:3) == 'EME') THEN
    NTYPGU=4
  ELSEIF (CLIDEN(1:3) == 'CEP') THEN
    NTYPGU=5
  ELSE
    NTYPGU=1
    WRITE(UNIT=NULOUT,FMT='(A45,A16)')&
     & ' TYPE D IDENTIFICATEUR NON PREVU : CLIDEN  : ',CLIDEN  
  ENDIF
  
  NDATGU=IDATEF(1)*10000+IDATEF(2)*100+IDATEF(3)
  NRESGU=IDATEF(4)*100+IDATEF(5)
  
  IF (IDATEF(9) == 10.OR.IDATEF(9) == 0) THEN
    IF (IDATEF(7) == 0) THEN
      WRITE(UNIT=NULOUT,FMT='(A45)')&
       & ' LE GUESS EST UNE ANALYSE NON INITIALISEE    '  
      IF (.NOT.LINFLAT) THEN
        NECHGU=0
      ELSE
        NECHGU=6
      ENDIF
    ELSE
      NECHGU=IDATEF(7)
    ENDIF
  ELSEIF (IDATEF(9) == 1) THEN
    IF (IDATEF(7) /= 0) THEN
      WRITE(UNIT=NULOUT,FMT='(A45)')&
       & ' INCOHERENCE ON A UNE ANALYSE                '  
      CALL ABOR1('SUGRIDA: ABOR1 CALLED')
    ENDIF
    WRITE(UNIT=NULOUT,FMT='(A45)')&
     & ' LE GUESS EST UNE ANALYSE INITIALISEE        '  
    NECHGU=0
  ELSE
    CALL ABOR1('SUGRIDA: ABOR1 CALLED')
  ENDIF

! Close file
  CALL FAIRME (IREP, IUNTIN, 'UNKNOWN')

ENDIF

IF (LHOOK) CALL DR_HOOK ('SUGRIDA:SUGRIDA_CHECK_C701',1,ZHOOK_HANDLE)

END SUBROUTINE SUGRIDA_CHECK_C701
  
END SUBROUTINE SUGRIDA

