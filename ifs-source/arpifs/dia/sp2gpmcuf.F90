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

SUBROUTINE SP2GPMCUF(YDGEOMETRY,YDMCUF,PMCUFGP)

!**** *SP2GPMCUF*  - Convert spectral monitoring of coupling update frequencies into gridpoint fields

!     Purpose.
!     --------
!        To convert spectral monitoring of coupling update frequencies into gridpoint fields prior to interpolations

!**   Interface.
!     ----------
!        *CALL* *SP2GPMCUF(...)

!        Explicit arguments :
!        --------------------
!           YDGEOMETRY : input model geometry
!           PMCUFGP : monitoring of coupling-update frequency in gridpoint




!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------
!        See documentation

!     Externals.    
!     ----------    

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib  *METEO-FRANCE*
!        Original : 29-Sep-2017

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMCUF            , ONLY : TMCUF
USE YOMCT0             , ONLY : LELAM

!     ------------------------------------------------------------------

IMPLICIT NONE
TYPE(GEOMETRY)   ,INTENT(IN)  :: YDGEOMETRY
TYPE(TMCUF)      ,INTENT(IN)  :: YDMCUF
REAL(KIND=JPRB)  ,INTENT(OUT) :: PMCUFGP(YDGEOMETRY%YRDIM%NPROMA,YDMCUF%NCUFNR,YDGEOMETRY%YRDIM%NGPBLKS)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JCUFNR, JKGLO, IBL, IST, IEND, JROF
REAL(KIND=JPRB) :: ZMCUF(YDGEOMETRY%YRGEM%NGPTOT)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "esperee.intfb.h"
#include "speree.intfb.h"

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SP2GPMCUF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDGEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, NSPEC2=>YDDIM%NSPEC2, NPROMA=>YDDIM%NPROMA)
!      -----------------------------------------------------------

DO JCUFNR=1,YDMCUF%NCUFNR
  IF (LELAM) THEN
    CALL ESPEREE(YDGEOMETRY,1,1,YDMCUF%RMCUFFP(1:NSPEC2,0,JCUFNR),ZMCUF)
  ELSE
    CALL SPEREE(YDGEOMETRY,1,1,YDMCUF%RMCUFFP(1:NSPEC2,0,JCUFNR),ZMCUF)
  ENDIF
  DO JKGLO=1,NGPTOT,NPROMA
    IBL=(JKGLO-1)/NPROMA+1
    IST=1
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    DO JROF=IST,IEND
      PMCUFGP(JROF,JCUFNR,IBL)=ZMCUF(JKGLO+JROF-1)
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SP2GPMCUF',1,ZHOOK_HANDLE)
END SUBROUTINE SP2GPMCUF
