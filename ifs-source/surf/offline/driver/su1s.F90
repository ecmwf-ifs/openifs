SUBROUTINE SU1S
USE PARKIND1  ,ONLY : JPIM     ,JPRB,  JPRD
USE YOMHOOK, ONLY   : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLOG1S , ONLY : LDBGS1   ,LACCUMW , LSEMISS , CFFORC &
                     &,CMODID   ,CVERID  , LRESET  , NACCUR &
                     &,NDIMCDF  &
                     &,LWREFL   ,LWRWAT  , LWRSUS  , LWRSUB &
                     &,LWREVA   ,LWRCLD  , LWRGG   , LWRCLM &
		             &,LWRLKE   ,LWROCP  , LWROCD  , LWROCR &
                     &,LOFFL    ,CFOUT    ,CFSURF  , CFINIT &
                     &,LWRCO2,LWRVEG,LWRFRA,LWRFOR,LWRBIO,LWRTIL,LWRVTY,CSITE &
                     &,LWREXT,LWRD2M,LWRGGD,NDLEVEL,IDBGS1,NCDFTYPE,LNCSNC

USE YOMLUN1S , ONLY : NULNAM   ,NULOUT

#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SU1s* - Initialize common YOMLOG1S controlling multi-column surface model 

!     Purpose.
!     --------
!           Initialize YOMLOG1S, the common that includes 
!           the basic switches for the multi-column model.

!**   Interface.
!     ----------
!        *CALL* *SU1S from SU0YOM1S

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation 
!        of the one-column surface model

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo   *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-22
!        BART VD HURK (KNMI): ADJUSTED FOR MULTI-COLUMN MODE
!        E. DUTRA - NEW OUTPUT FILE FOR LAKES 04/07/2008

!     ------------------------------------------------------------------
#endif
IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "nam1s.h"

IF (LHOOK) CALL DR_HOOK('SU1S',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

!*       1.1   Set default values for the multi-column surface model switches.
!              ------------------------------------------------------------

IDBGS1=0
LDBGS1=.FALSE.
LACCUMW=.TRUE.
LRESET=.TRUE.
LSEMISS=.FALSE.
CFFORC='netcdf'
CFOUT='netcdf'
CFSURF='netcdf'
CFINIT='netcdf'
CMODID='tessel'
CVERID='1.1'
NACCUR=1
NDLEVEL=0
NDIMCDF=2
NCDFTYPE=4
LNCSNC=.FALSE.
LWREFL=.TRUE.
LWRWAT=.TRUE.
LWRSUS=.TRUE.
LWRSUB=.FALSE. !.TRUE.
LWREVA=.TRUE.
LWRCLD=.TRUE.
LWRGG=.TRUE.
LWRCLM=.TRUE.
LOFFL=.TRUE.
LWRLKE=.FALSE.
LWROCP=.FALSE.   !KPP
LWROCD=.FALSE.   !KPP
LWROCR=.FALSE.   !KPP

LWRCO2=.TRUE.
LWRVEG=.FALSE.
LWRBIO=.FALSE.
LWRTIL=.FALSE.
LWREXT=.FALSE.
LWRVTY=.FALSE.
LWRFRA=.FALSE.
LWRFOR=.FALSE.
CSITE='exp?'
LWRD2M=.TRUE.
LWRGGD=.TRUE.

!     ------------------------------------------------------------------

!*       2.    Modify default values.
!              ----------------------

REWIND(NULNAM)
READ(NULNAM,NAM1S)

!     ------------------------------------------------------------------

!*       3.    CHECK

IF(.NOT.LACCUMW .AND. LRESET) THEN
  LRESET = .FALSE.
  WRITE(NULOUT,*)'LRESET =.T. INCOMPATIBLE WITH LACCUM = .F.'
  WRITE(NULOUT,*)'LRESET SET TO .F.'
ENDIF

IF(NACCUR /= 1.AND. NACCUR /= 2)THEN
  NACCUR=1
  WRITE(NULOUT,*)'NACCUR SET TO ',NACCUR
ENDIF

IF ( NCDFTYPE /= 3 .AND. NCDFTYPE /= 4 .AND. NCDFTYPE /= 33 .AND. NCDFTYPE /= 44  ) THEN
  WRITE(NULOUT,*) ' NCDFTYPE must be 3,33 or 4 :',NCDFTYPE
  CALL ABORT()
ENDIF
!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SU1S
