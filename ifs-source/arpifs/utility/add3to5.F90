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

SUBROUTINE ADD3TO5(YDMP,YDML_GCONF,YDSP3)

!**** *ADD5TO3* - ROUTINE TO ADD SP5 AND SP3 SPECTRAL ARRAYS

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        *CALL* *ADD3TO5*

!        IMPLICIT ARGUMENTS :  SP5Ax and SPAx
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      J.N. Thepaut based on ADD5TO3  *ECMWF*
!      ORIGINAL : 94-07-29

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!      O.Marsden     August 2016 Remove use of SPA3, replace by explicit argument YDSP3
!     ------------------------------------------------------------------
!        O.Riviere     23-Nov-10 Sum on SPGFL done only if pointer allocated

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMMP                  , ONLY : TMP
USE PARKIND1               , ONLY : JPRB
USE YOMHOOK                , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMSP5                 , ONLY : SPA5
USE YOMCT0                 , ONLY : LELAM
USE SPECTRAL_FIELDS_MOD    , ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

IMPLICIT NONE

TYPE(TMP)                    , INTENT(IN) :: YDMP
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(IN) :: YDML_GCONF
TYPE(SPECTRAL_FIELD)         , INTENT(IN) :: YDSP3

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

!*    1.    ADD SP5Ax to SPAx.
!           ------------------

IF (LHOOK) CALL DR_HOOK('ADD3TO5',0,ZHOOK_HANDLE)
ASSOCIATE(NUMSPFLDS=>YDML_GCONF%YGFL%NUMSPFLDS, &
 & NFD2D=>YDML_GCONF%YRDIMF%NFD2D, &
 & NPSP=>YDMP%NPSP)
SPA5%VOR(:,:) = YDSP3%VOR(:,:) + SPA5%VOR(:,:)
SPA5%DIV(:,:) = YDSP3%DIV(:,:) + SPA5%DIV(:,:)
SPA5%T(:,:)   = YDSP3%T(:,:)   + SPA5%T(:,:)
IF (NUMSPFLDS>0) THEN
  SPA5%GFL(:,:,:)   = YDSP3%GFL(:,:,:)   + SPA5%GFL(:,:,:)
ENDIF

IF(NPSP == 1) THEN
  SPA5%SP2D(:,1:NFD2D) = YDSP3%SP2D(:,1:NFD2D) + SPA5%SP2D(:,1:NFD2D)
ENDIF

IF (LELAM) SPA5%SP1D(:,:) = YDSP3%SP1D(:,:) + SPA5%SP1D(:,:)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('ADD3TO5',1,ZHOOK_HANDLE)
END SUBROUTINE ADD3TO5
