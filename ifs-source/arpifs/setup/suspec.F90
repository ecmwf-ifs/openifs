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

SUBROUTINE SUSPEC(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,&
           &YDSP,CDPATH,KFILE,LDRDGRIDSP,LDINOR,LDSPOR,YDMCUF)

!**** *SUSPEC*  - Routine to initialize the spectral fields of the model.

!     Purpose.
!     --------
!           Initialize the spectral fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUSPEC(.....)*

!        Explicit arguments :
!        --------------------
!        KFILE : an indicator for which spectral file is to be read
!               KFILE = 0 the CFNISH file is read
!               KFILE = 1 the CFNGSH file is read
!               KFILE = 2 the CFNRF file is read
!               KFILE = 3 the CFNBS file is read (bias, for recursive dfi)
!               KFILE =12 the CFNBGHRSH file is read
!               PSPVOR etc. - spectral fields
!               LDRDGRIDSP: T (resp. F): GMV, GMVS and spectral GFL read as
!                g.p. fields (resp. spectral fields) on files.
!               LDINOR - switch for initializing the orography
!               LDSPOR - switch indicating whether spectral orography
!                        field has been read

!        Implicit arguments :
!        --------------------
!        The spectral fields of the model.
!        The boundary condition fields of the model.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-12-14

!     Modifications.
!     --------------
!      E.Holm        13-Nov-2008 (KFILE=12)
!      K. Yessad (July 2014): Move some variables.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      P.Marguinaud  10-Oct-2014 Use SUSPECA_GP
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      K. Yessad (Dec 2016): Prune obsolete options.
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN                 , ONLY : TDYN
USE YOMDYNA                , ONLY : TDYNA
USE YEMLBC_MODEL             , ONLY : TELBC_MODEL
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE YOMGFL                 , ONLY : TGFL
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN                 , ONLY : NULOUT
USE YOMMP0                 , ONLY : NPRINTLEV
USE YOMCT0                 , ONLY : LR2D, N3DINI, LECMWF
USE YOMARG                 , ONLY : NGRIBFILE
USE SPECTRAL_FIELDS_MOD    , ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)
USE YOMMCUF                , ONLY : TMCUF

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               ,INTENT(INOUT) :: YDGEOMETRY
TYPE(TGFL)                   ,INTENT(INOUT) :: YDGFL
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN)    :: YDML_GCONF
TYPE(TDYN)                   ,INTENT(IN)    :: YDDYN
TYPE(TDYNA)                  ,INTENT(IN)    :: YDDYNA
TYPE(TELBC_MODEL)                ,INTENT(IN)    :: YDML_LBC
TYPE(SPECTRAL_FIELD)         ,INTENT(INOUT) :: YDSP
CHARACTER(LEN=*)    ,OPTIONAL,INTENT(IN)    :: CDPATH
INTEGER(KIND=JPIM)  ,OPTIONAL,INTENT(IN)    :: KFILE 
LOGICAL             ,OPTIONAL,INTENT(IN)    :: LDRDGRIDSP
LOGICAL             ,OPTIONAL,INTENT(IN)    :: LDINOR 
LOGICAL             ,OPTIONAL,INTENT(INOUT) :: LDSPOR 
TYPE(TMCUF)         ,OPTIONAL,INTENT(INOUT) :: YDMCUF

!     ------------------------------------------------------------------

LOGICAL :: LLDATA,LLINOR,LLSPOR,LLRDGRIDSP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "suspeca_gp.intfb.h"
#include "suorog.intfb.h"
#include "suspeca.intfb.h"
#include "suspecb.intfb.h"
#include "suspecg.intfb.h"
#include "suspecg1.intfb.h"
#include "suspecg2.intfb.h"

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSPEC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NFLEVL=>YDDIMV%NFLEVL)
!     ------------------------------------------------------------------

CALL GSTATS(19,0)

IF (NPRINTLEV >=2) THEN
  WRITE(NULOUT,'(A,I3)') ' ENTERING SUSPEC WITH KFILE= ',KFILE
  CALL FLUSH(NULOUT)
ENDIF

LLINOR=.FALSE.
IF(PRESENT(LDINOR)) LLINOR = LDINOR

IF (PRESENT(LDRDGRIDSP)) THEN
  LLRDGRIDSP=LDRDGRIDSP
ELSE
  LLRDGRIDSP=.FALSE.
ENDIF

!     ------------------------------------------------------------------

!*       1.    INITIALIZE SPECTRAL FIELDS FOR BAROTROPIC MODELS.
!              -------------------------------------------------

IF (LR2D) THEN

  CALL SUSPECB(YDGEOMETRY,YDML_GCONF,YDDYN,KFILE,YDSP)
  LLINOR=.TRUE.
  LLSPOR=.TRUE.

!      -----------------------------------------------------------

!*       2.    INITIALIZE SPECTRAL FIELDS (3D)
!              -------------------------------

ELSE

  IF ((NGRIBFILE /= 1).AND.(N3DINI == 0)) THEN

!*       2.1   INITIALIZE SPECTRAL FIELDS FROM *FA* FILE.

    IF(LLRDGRIDSP) THEN
      ! fields are stored in grid-point on files.
      CALL SUSPECA_GP(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSP,KFILE,LLDATA,LLINOR,YDMCUF=YDMCUF)
      IF (NPRINTLEV >=2) THEN
        WRITE(NULOUT,*) ' SUSPEC: SUSPECA_GP called '
        CALL FLUSH(NULOUT)
      ENDIF
    ELSE
      ! fields are stored in spectral on files.
      CALL SUSPECA(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYN,YDDYNA,YDML_LBC,YDSP,KFILE,LLDATA,LLINOR,YDMCUF=YDMCUF)
      IF (NPRINTLEV >=2) THEN
        WRITE(NULOUT,*) ' SUSPEC: SUSPECA called '
        CALL FLUSH(NULOUT)
      ENDIF
    ENDIF
    IF(PRESENT(LDSPOR))LDSPOR=.TRUE.

  ELSEIF ( N3DINI == 0 ) THEN

!*       2.2   INITIALIZE SPECTRAL FIELDS FROM GRIB FILES

    IF(LLRDGRIDSP) THEN
      CALL ABOR1('SUSPEC: option LLRDGRIDSP not yet implemented for GRIB')
      ! future CALL SUGRIDSPG.
    ELSE
      CALL SUSPECG(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA%LNHDYN,YDSP%SP3D,YDSP%SP2D,CDPATH,KFILE,LLINOR,LLSPOR)
    ENDIF

  ELSEIF( N3DINI == 1 ) THEN

!*       2.3   INITIALIZE SPECTRAL FIELDS FROM ARTIFICIAL DATA
    IF(LLRDGRIDSP) THEN
      CALL ABOR1('SUSPEC: option LLRDGRIDSP not yet implemented for GRIB')
      ! future CALL SUGRIDSPG1.
    ELSE
      CALL SUSPECG1(YDGEOMETRY,YDML_GCONF%YGFL,YDSP%VOR(:,:),YDSP%DIV(:,:),&
       & YDSP%T(:,:),YDSP%GFL(:,:,:),&
       & YDSP%SP(:),YDSP%OROG(:),&
       & LLINOR,LDSPOR)
    ENDIF

  ELSEIF( N3DINI == 2 ) THEN

    !*    2.4   INITIALIZE SPECTRAL FIELDS FOR ACADEMIC TEST CASES
    !              STARTING FROM ARTIFICIAL DATA (NSUPERSEDE=0 could work)
    IF(LLRDGRIDSP) THEN
      CALL ABOR1('SUSPEC: option LLRDGRIDSP not yet implemented for GRIB')
      ! future CALL SUGRIDSPG2 ???
    ELSE  
      CALL SUSPECG2(YDGEOMETRY,YDGFL,YDML_GCONF,YDSP%SP3D,YDSP%SP2D)
      LLSPOR = .TRUE.
    ENDIF

  ELSE
      CALL ABOR1('SUSPEC: ONLY N3DINI = 0, 1 or 2 are implemented currently')
  ENDIF

!*       2.5   INITIALIZE OVERDIMENSIONED SPECTRAL FIELDS

  YDSP%SP3D(NFLEVL+1:,:,:) = 0.0_JPRB
   
ENDIF

!     ------------------------------------------------------------------

!*       3.    INITIALIZE GRIDPOINT OROGRAPHY.
!              -------------------------------

!*       3.1   INITIALIZE GRIDPOINT OROGRAPHY

IF(LECMWF) THEN
  IF(LLINOR.AND.LLSPOR) THEN
    WRITE(NULOUT,'(" OROGRAPHY INITIALISED FROM SPECTRAL INPUT")')
    CALL SUOROG(YDGEOMETRY,YDSP%OROG)
  ENDIF
  IF(PRESENT(LDSPOR)) LDSPOR=LLSPOR
ELSE
  IF(LLINOR) THEN
    WRITE(NULOUT,'(" OROGRAPHY INITIALISED FROM SPECTRAL INPUT")')
    CALL SUOROG(YDGEOMETRY,YDSP%OROG)
  ENDIF
ENDIF

!     ------------------------------------------------------------------

CALL GSTATS(19,1)

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSPEC',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPEC
