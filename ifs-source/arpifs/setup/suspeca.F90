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

SUBROUTINE SUSPECA(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSPEC,KFILE,LDATA,LDINOR,YDMCUF)

!**** *SUSPECA*  - Initialize the spectral fields from FA (Arpege or Aladin)

!     Purpose.
!     --------
!         Initialize the spectral fields of the model from ARPEGE ALADIN file.

!**   Interface.
!     ----------
!        *CALL* *SUSPECA(.....)*

!        Explicit arguments :
!        --------------------
!        KFILE : an indicator for which spectral file is to be read  (I)
!        LDATA : .TRUE. if data is returned from file                (O)
!        PSPVOR etc. - spectral fields                               (O)
!        LDINOR - switch for initializing the orography              (I)

!        Implicit arguments :
!        --------------------
!        See modules above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!        see calls below.

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      Ryad El Khatib *METEO-FRANCE*
!      Original     : 99-04-23 (replacement of suspecasm/dm)

!     Modifications.
!     --------------
!      R. El Khatib : 03-08-05 gfl
!      R. El Khatib : 03-04-17 Cleanups
!      R. El Khatib : 03-05-14 New interfaces euvgeovd
!      G. Hello     : supress the call to gnhpdvdconv at this level
!      M.Hamrud     : 01-Oct-2003 CY28 Cleaning
!      M.Hamrud     : 10-Jan-2004 CY28R1 Cleaning
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      R. Brozkova  :  05-07-23 Handle the divergence "X" term not to be read
!      R. El Ouaraini & R. El Khatib : 05-07-28 : New fields for the monitoring
!                     of update frequency of the coupling files for Aladin
!      Y. Bouteloup : 05-10-18 YCVGQ not read even it is spectral
!      K. Yessad    : 28 Feb 2006: cleanings (useless dummies, indentations).
!      Y. Seity     : 23-06 2006: remove YCVGQ specificities and replace by more
!                     general tests on NREQIN
!      R. El Khatib : 30-Mar-2010 reduce overhead of call to spreord
!      P.Marguinaud : 28-05-2010 Change SUMPIOH interface
!      R. El Khatib : 10-Aug-2011 NIOLEVG management
!      R. El Khatib : 01-Feb-2012 Extend I/O processors up to NPRTRV*NSTRIN
!      P.Marguinaud : 26-Apr-2012 Refactor using IOSPECA_MOD + fix bug in PSPGFL
!                                 initialization
!      P.Marguinaud : 11-Sep-2012 Refactor using IOFLDDESC_MOD, add call to SUSPECA_GP
!      P.Marguinaud : 10-Oct-2013 Use DIST_SPEC & EDIST_SPEC
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      P.Marguinaud : 10-Oct-2014 Read spectral fields using the IO server and use RDFA2SP
!     -------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0                 , ONLY : N3DINI, LSUSPECA_GP
USE YOMMP0                 , ONLY : NPRTRV, MYSETV
USE IOSPECA_MOD            , ONLY : IOSPECA_COUNT,&
                                  & IOSPECA_CTX,&
                                  & NSPECACT_READ,&
                                  & IOSPECA_START,&
                                  & IOSPECA_FINISH,&
                                  & IOSPECA_SELECTD,&
                                  & IOSPECA_VSETOFF,&
                                  & IOSPECA_SELECTF
USE ERLBC_MOD              , ONLY : YSUSPCTX
USE IOFLDDESC_MOD          , ONLY : IOFLDDESC
USE YEMLBC_MODEL             , ONLY : TELBC_MODEL
USE SPECTRAL_FIELDS_MOD
USE YOMMCUF                , ONLY : TMCUF


!     -------------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT (IN)        :: YDGEOMETRY
TYPE(TGFL)          , INTENT (INOUT)     :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN) :: YDML_GCONF
TYPE(TDYN)          , INTENT (IN)        :: YDDYN
TYPE(TDYNA)         , INTENT (IN)        :: YDDYNA
TYPE(TELBC_MODEL)       , INTENT (IN)        :: YDML_LBC
TYPE(SPECTRAL_FIELD), INTENT (INOUT)     :: YDSPEC
INTEGER(KIND=JPIM)  , INTENT (IN)        :: KFILE 
LOGICAL             , INTENT (OUT)       :: LDATA 
LOGICAL             , INTENT (IN)        :: LDINOR 
TYPE(TMCUF)         , INTENT(INOUT), OPTIONAL :: YDMCUF
!     -------------------------------------------------------------------------
!     IFLDSPG: total number of fields to read
!     IVSETOFF: V-set offset for re-ordered fields
!     ZSPBUFL: a chunk of of spectral buffer fields (local)
INTEGER (KIND=JPIM) :: IFLDSPG, IFLDSPS
INTEGER (KIND=JPIM) :: I3DINI
LOGICAL   :: LLINOR

INTEGER(KIND=JPIM) :: IVSETOFF (NPRTRV)
INTEGER(KIND=JPIM) :: IREP
TYPE (IOFLDDESC),    ALLOCATABLE :: YLFLDSC (:)
REAL(KIND=JPRB),     ALLOCATABLE :: ZSPBUFL (:,:)
TYPE (IOSPECA_CTX) :: YLIOSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "suspeca_gp.intfb.h"
#include "suspeca_fixup.intfb.h"
#include "rdfa2sp.intfb.h"
#include "suspeca_map_part2.intfb.h"

!     -------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSPECA',0,ZHOOK_HANDLE)
!     -------------------------------------------------------------------------

IF (LSUSPECA_GP) THEN
  
  CALL SUSPECA_GP(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSPEC,KFILE,LDATA,LDINOR)

  IF (LHOOK) CALL DR_HOOK('SUSPECA',1,ZHOOK_HANDLE)
  RETURN

ENDIF

CALL SUSPECA_MAP_PART2(YDGEOMETRY,YDML_GCONF%YGFL,YDDYN,YDDYNA,IREP,YSUSPCTX,YDSPEC)


IF (IREP == 0) THEN
  IF (LHOOK) CALL DR_HOOK('SUSPECA',1,ZHOOK_HANDLE)
  RETURN
ENDIF

I3DINI=N3DINI
LLINOR=LDINOR

IF (KFILE == 4) THEN
  CALL ABOR1('SUSPECA : BLENDING REMOVED')
ELSE
  LDATA = .TRUE.
ENDIF

! Count fields
CALL IOSPECA_COUNT(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_READ,IFLDSPG,LDINOR=LLINOR,K3DINI=I3DINI,YDMCUF=YDMCUF)

IF (IFLDSPG > 0) THEN

  CALL IOSPECA_START(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,NSPECACT_READ,YLIOSP,YDSPEC)

  ! Setup various local variables related to file
  ALLOCATE (YLFLDSC (IFLDSPG))
  CALL IOSPECA_SELECTD(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_READ,YLFLDSC,LDINOR=LLINOR,K3DINI=I3DINI,YDMCUF=YDMCUF)

  ! Setup chunks of fields
  CALL IOSPECA_VSETOFF (YLFLDSC%IVSET, KVSETOFF=IVSETOFF)

  IFLDSPS = COUNT (YLFLDSC%IVSET == MYSETV)


  ALLOCATE (ZSPBUFL (YDGEOMETRY%YRDIM%NSPEC2, IFLDSPS))
  CALL RDFA2SP (YDGEOMETRY, YDML_GCONF%YRRIP, IFLDSPG, IFLDSPS, &
              & ZSPBUFL, YLFLDSC, KFILE=KFILE, YDML_LBC=YDML_LBC)

  CALL IOSPECA_SELECTF(YDGEOMETRY,YDML_GCONF%YGFL,YDDYN, &
 &                     NSPECACT_READ, YLIOSP, IFLDSPS, IVSETOFF(MYSETV), 0, YLFLDSC, ZSPBUFL, YDSPEC, YDMCUF=YDMCUF)

  DEALLOCATE (ZSPBUFL)

! Compute Vorticity, Divergence & mean wind from geographical wind

  CALL IOSPECA_FINISH (YDGEOMETRY, YDDYNA,NSPECACT_READ, YLIOSP, YDSPEC)

! Release descriptors

  DEALLOCATE (YLFLDSC)

  CALL SUSPECA_FIXUP(YDGEOMETRY,YDML_GCONF%YGFL,KFILE,YDSPEC%GFL,YDSPEC%SP)

ENDIF ! (IFLDSPG > 0)

IF (LHOOK) CALL DR_HOOK('SUSPECA',1,ZHOOK_HANDLE)

END SUBROUTINE SUSPECA

