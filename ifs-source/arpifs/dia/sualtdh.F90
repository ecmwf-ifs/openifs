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

SUBROUTINE SUALTDH(YDDIMV,YDML_DIAG,YDARPHY,YDPHY)
!****
!     ------------------------------------------------------------------

!     SUALTDH * DECLARATION DES TABLEAUX DE CUMULS
!     DIAGNOSTICS PAR DOMAINES HORIZONTAUX
!     ------------------------------------

!       -------------------------------------------------------
!     BUT
!     ---
!       COMPTE TENU DES INDICATEURS LOGIQUES, UTILISANT LE RESULTAT
!     DE *SUMDDH* ET LE NOMBRE DE DOMAINES INTERNES (NDHIDH),
!     CE PROGRAMME ALLOUE LA MEMOIRE DES TABLEAUX DE CUMULS

!     ARGUMENTS D ENTREE
!     ------------------

!     ARGUMENTS IMPLICITES
!     --------------------
!       UNITE NULOUT
!       INDICATEURS LOGIQUES, COMMON /YOMLDDH/
!         DANS LE COMDECK * YOMMDDH * POUR LES MASQUES
!       DIMENSIONS PROPRES AUX MASQUES DIAGNOSTIQUES DDH, COMMON /DDHDIM/
!       ET DIMENSIONS POUR LES TABLEAUX DE CUMULS
!         EN PARTICULIER NDHIDH
!       POINTEURS PROPRES AUX MASQUES DIAGNOSTIQUES DDH, COMMON /POMMDDH/
!       DECLARATION DES MASQUES ET LEURS AUXILLIAIRES
!         DANS LE COMDECK * YOMTDDH * POUR LES TABLEAUX DE CUMULS
!       POINTEURS PROPRES AUX TABLEAUX DE CUMULS DDH, COMMON /POMTDDH/
!       DECLARATION DES TABLEAUX DE CUMULS

!     AUTRES ENTREES
!     --------------

!     SORTIES
!     -------
!        INITIALISE LES POINTEURS /POMTDDH/
!        INITIALISE LES TABLEAUX DE CUMULS A ZERO

!     ECRIT PAR
!     --------- ALAIN JOLY

!     25/1/91

!     Modifications:
!     --------------
!        R. El Khatib : 01-08-07 Pruning options
!        M.Hamrud     : 20-Apr-2004 OpenMP improvements of DDH
!        T.Kovacic    : 16-Mar-2007 Allocation of APFT  
!        T.Kovacic    : 01-May-2007 INIAPFT_BP002
!        O.Riviere: 01-oct-2008 removal of iniapft_bp002
!        R.Brozkova   : 19-May-2011 temporary fix of ddh for ALARO
!        K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE MODEL_DIAGNOSTICS_MOD , ONLY : MODEL_DIAGNOSTICS_TYPE
USE YOMDIMV               , ONLY : TDIMV
USE PARKIND1              , ONLY : JPIM, JPRB
USE YOMHOOK               , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0                , ONLY : NPRINTLEV
USE YOMCT0                , ONLY : LALLOPR
USE YOMLUN                , ONLY : NULOUT
USE OML_MOD               , ONLY : OML_MAX_THREADS
USE YOMPHY                , ONLY : TPHY
USE YOMARPHY              , ONLY : TARPHY

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV), INTENT(IN)   :: YDDIMV
TYPE(TARPHY),INTENT(IN)   :: YDARPHY
TYPE(MODEL_DIAGNOSTICS_TYPE),INTENT(INOUT):: YDML_DIAG
TYPE(TPHY)  ,INTENT(IN)   :: YDPHY
INTEGER(KIND=JPIM) :: IU,ITHRMAX
LOGICAL :: LLP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "iniapft_bp002.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUALTDH',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & LSDDH=>YDML_DIAG%YRLDDH%LSDDH, &
 & NDHCS=>YDML_DIAG%YRMDDH%NDHCS, NDHCV=>YDML_DIAG%YRMDDH%NDHCV, NDHIDH=>YDML_DIAG%YRMDDH%NDHIDH, & 
 & LMPA=>YDARPHY%LMPA, L3MT=>YDPHY%L3MT, LSTRAPRO=>YDPHY%LSTRAPRO,YDTDDH=>YDML_DIAG%YRTDDH)
!     ------------------------------------------------------------------

LLP = NPRINTLEV >= 1.OR. LALLOPR
IU = NULOUT
ITHRMAX=OML_MAX_THREADS()
IF (LLP) THEN
  WRITE (NULOUT,*) ' DDH * DECLARATION DES TABLEAUX DE CUMULS'
  WRITE (NULOUT,*) ' ----------------------------------------'
  WRITE (NULOUT,*) ' '
  WRITE (NULOUT,*) ' DDH * NDHCV = ',NDHCV,' NDHCS = ',NDHCS,&
   & ' NDHIDH = ',NDHIDH  
ENDIF


!     DECLARATIONS PROPREMENT DITES
!     -----------------------------

IF(NDHCV*(NFLEVG+1)*NDHIDH >= 1) THEN
  ALLOCATE(YDTDDH%HDCVB0(NDHCV*(NFLEVG+1),NDHIDH,ITHRMAX))
  IF(LLP)WRITE(IU,9) 'HDCVB0   ',SIZE(YDTDDH%HDCVB0),SHAPE(YDTDDH%HDCVB0)
  ALLOCATE(YDTDDH%HDCVB1(NDHCV*(NFLEVG+1),NDHIDH,ITHRMAX))
  IF(LLP)WRITE(IU,9) 'HDCVB1   ',SIZE(YDTDDH%HDCVB1),SHAPE(YDTDDH%HDCVB1)
  ALLOCATE(YDTDDH%CFLDNAMES3D(NDHCV))
  YDTDDH%CFLDNAMES3D(:)='@@@@@@@@@@'
  ALLOCATE(YDTDDH%CFLDTYPES3D(NDHCV))
  YDTDDH%CFLDTYPES3D(:)=' '
ENDIF
IF(NDHCS*NDHIDH >= 1) THEN
  ALLOCATE(YDTDDH%HDCS0(NDHCS,NDHIDH,ITHRMAX))
  IF(LLP)WRITE(IU,9) 'HDCS0    ',SIZE(YDTDDH%HDCS0),SHAPE(YDTDDH%HDCS0)
  ALLOCATE(YDTDDH%HDCS1(NDHCS,NDHIDH,ITHRMAX))
  IF(LLP)WRITE(IU,9) 'HDCS1   ',SIZE(YDTDDH%HDCS1),SHAPE(YDTDDH%HDCS1)
  ALLOCATE(YDTDDH%CFLDNAMES2D(NDHCS))
  YDTDDH%CFLDNAMES2D(:)='@@@@@@@@@@'
  ALLOCATE(YDTDDH%CFLDTYPES2D(NDHCS))
  YDTDDH%CFLDTYPES2D(:)=' '
ENDIF

IF (.NOT.LMPA.AND.( L3MT.OR.LSTRAPRO) ) CALL INIAPFT_BP002(YDML_DIAG%YRLDDH,YDML_DIAG%YRSDDH,YDPHY,NULOUT)

!*    INITIALISATION DES TABLEAUX DE CUMULS
!     -------------------------------------

IF ( LSDDH ) THEN
  YDTDDH%HDCVB0(:,:,:)=0.0_JPRB
  YDTDDH%HDCVB1(:,:,:)=0.0_JPRB
  YDTDDH%HDCS0(:,:,:)=0.0_JPRB
  YDTDDH%HDCS1(:,:,:)=0.0_JPRB
ENDIF


9 FORMAT(1X,'ARRAY ',A10,' ALLOCATED ',8I8)

!     ------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUALTDH',1,ZHOOK_HANDLE)
END SUBROUTINE SUALTDH
