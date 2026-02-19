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

#ifdef VPP
!OCL SCALAR
#endif
SUBROUTINE SUALLO(YDGEOMETRY,YDMODEL)

!**** *SUALLO * - Routine to allocate space for global variables

!     Purpose.
!     --------
!           Allocate space for the global fields.

!**   Interface.
!     ----------
!        *CALL* *SUALLO*

!     Explicit arguments :  None
!     --------------------
!        Called by SU0YOMA.

!     Implicit arguments :
!     --------------------
!        Pointers of comdecks

!     Method.
!     -------
!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-26

!     Modifications.
!     --------------
!      Modified 04-09-15 by F.Bouyssel: allocation of STPREH
!      M.Hamrud      01-Jul-2006 Revised surface fields
!      K. Yessad 15-May-2006: memory optimisations for stretched geometry
!      M. Bellus 28-Sep-2006: ALARO-0 phasing (Smith RH for Lopez scheme)
!      N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!      30-Jun-2008 J. Masek   Auxiliary quantities for SLHD interpolators.
!      K. Yessad (Sep 2008): LREGETA -> LREGETA + LVFE_REGETA.
!      A. Fouilloux (sep 2009) : remove call to sualobs
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived types TVAB, TVETA, TVSLETA, TCSGLEG
!      R. El Khatib 10-Aug-2011 NIOLEVG management
!      K. Yessad (Apr 2012): allocation of NVAUTF and NVAUTH moved in SUVERT.
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July): call to SUALDYN moved in SUDYN.
!      JJMorcrette 20130805 MACC-derived aerosol climatology
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : NPRINTLEV
USE YOMCT0       , ONLY : LECMWF, LALLOPR
USE YOMLUN       , ONLY : NULOUT
USE OML_MOD      , ONLY : OML_INIT_LOCK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(MODEL)    ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM) :: IU

LOGICAL ::  LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
#ifdef __INTEL_COMPILER
INTRINSIC :: SHAPE ! Fails with Intel compiler as it thinks we use ATLAS shape function
#endif


!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUALLO',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, &
 & YDMP=>YDGEOMETRY%YRMP, YDPHY=>YDMODEL%YRML_PHY_MF%YRPHY,YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY, &
 & YDDPHY=>YDMODEL%YRML_PHY_G%YRDPHY, YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH, YDECLD=>YDMODEL%YRML_PHY_EC%YRECLD, &
 & YDEAERD=>YDMODEL%YRML_PHY_RAD%YREAERD, &
 & YDERAD=>YDMODEL%YRML_PHY_RAD%YRERAD, YDEOVLP=>YDMODEL%YRML_PHY_RAD%YREOVLP)
ASSOCIATE(NDGENG=>YDDIM%NDGENG, NDGSAG=>YDDIM%NDGSAG, &
 & NVCLIS=>YDDPHY%NVCLIS, NTOZ2D=>YDDPHY%NTOZ2D, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & LRRMES=>YDPHY%LRRMES, LAGPHY=>YDEPHY%LAGPHY)
!     ------------------------------------------------------------------

!*       1.    ALLOCATE SPACE FOR ARRAYS.
!              --------------------------

!*       1.1   GRID POINT FIELDS.

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
IF (LLP) WRITE(NULOUT,'('' SUALLO PRINTOUTS '')')

!*       1.16  PHYSICS

IF (LRRMES.AND..NOT.LAGPHY) THEN
  ALLOCATE(YDTOPH%RMESOU(NFLEVG))
  IF(LLP)WRITE(IU,9) 'RMESOU   ',SIZE(YDTOPH%RMESOU),SHAPE(YDTOPH%RMESOU)
  ALLOCATE(YDTOPH%RMESOT(NFLEVG))
  IF(LLP)WRITE(IU,9) 'RMESOT   ',SIZE(YDTOPH%RMESOT),SHAPE(YDTOPH%RMESOT)
  ALLOCATE(YDTOPH%RMESOQ(NFLEVG))
  IF(LLP)WRITE(IU,9) 'RMESOQ   ',SIZE(YDTOPH%RMESOQ),SHAPE(YDTOPH%RMESOQ)
ENDIF

!*       1.18  EC PHYSICS CONSTANTS.

ALLOCATE(YDECLD%CETA(NFLEVG))
IF(LLP)WRITE(IU,9) 'CETA     ',SIZE(YDECLD%CETA),SHAPE(YDECLD%CETA)
ALLOCATE(YDEAERD%CVDAES(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAES   ',SIZE(YDEAERD%CVDAES),SHAPE(YDEAERD%CVDAES)
ALLOCATE(YDEAERD%CVDAEL(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAEL   ',SIZE(YDEAERD%CVDAEL),SHAPE(YDEAERD%CVDAEL)
ALLOCATE(YDEAERD%CVDAEU(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAEU   ',SIZE(YDEAERD%CVDAEU),SHAPE(YDEAERD%CVDAEU)
ALLOCATE(YDEAERD%CVDAED(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAED   ',SIZE(YDEAERD%CVDAED),SHAPE(YDEAERD%CVDAED)

ALLOCATE(YDERAD%CVDAESS(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAESS  ',SIZE(YDERAD%CVDAESS),SHAPE(YDERAD%CVDAESS)
ALLOCATE(YDERAD%CVDAEDU(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAEDU  ',SIZE(YDERAD%CVDAEDU),SHAPE(YDERAD%CVDAEDU)
ALLOCATE(YDERAD%CVDAEOM(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAEOM  ',SIZE(YDERAD%CVDAEOM),SHAPE(YDERAD%CVDAEOM)
ALLOCATE(YDERAD%CVDAEBC(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAEBC  ',SIZE(YDERAD%CVDAEBC),SHAPE(YDERAD%CVDAEBC)
ALLOCATE(YDERAD%CVDAESU(NFLEVG+1))
IF(LLP)WRITE(IU,9) 'CVDAESU  ',SIZE(YDERAD%CVDAESU),SHAPE(YDERAD%CVDAESU)

ALLOCATE(YDEOVLP%RA1OVLP(NFLEVG))
IF(LLP)WRITE(IU,9) 'RA1OVLP  ',SIZE(YDEOVLP%RA1OVLP),SHAPE(YDEOVLP%RA1OVLP)

!*       1.20  CLIMATOLOGICAL PHYSICS FIELDS

IF (LECMWF) THEN
  IF(NVCLIS*NTOZ2D > 0) THEN
    ALLOCATE(YDMODEL%YRML_CHEM%YROZO%TOZ2DG(NFLEVG*NVCLIS*NTOZ2D,NDGSAG:NDGENG))
    IF(LLP)WRITE(IU,9) 'TOZ2DG   ',SIZE(YDMODEL%YRML_CHEM%YROZO%TOZ2DG),SHAPE(YDMODEL%YRML_CHEM%YROZO%TOZ2DG)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!*       2.    ALLOCATE SPACE FOR ANALYSIS JOBS: CONT. VAR., COST FUNCTION,
!              ------------------------------------------------------------
!              OBSERVATIONS,GUESS.
!              -------------------


! Initialise lock of OpenMP
CALL OML_INIT_LOCK()

!     ------------------------------------------------------------------
9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALLO',1,ZHOOK_HANDLE)
END SUBROUTINE SUALLO
