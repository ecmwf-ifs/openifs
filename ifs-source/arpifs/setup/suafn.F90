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

SUBROUTINE SUAFN(YDNAMFPSCI,YDAFN,LDNEWCL,YGFL,LDNHDYN,KFPCONF,YDCOMPO,YDEAERATM,CDFPSFXFNAME)

!**** *SUAFN*  - INITIALIZE ARPEGE FIELD DESCRIPTORS

!     PURPOSE.
!     --------
!        SETS DEFAULT VALUES, THEN READS NAMELIST NAMAFN.

!**   INTERFACE.
!     ----------
!       *CALL* *SUAFN*

!        EXPLICIT ARGUMENTS
!        --------------------
!           KFPCONF : configuration of the post-processing :
!                     0 : vertical interpolation only (<CFPFMT='MODEL'>)
!                     1 : gridpoint post-processing, possibly with spectral filters (<NFPOS=1>)
!                     2 : gridpoint/spectral post-processing (spectral outputs possible) (<NFPOS=2>)
!           CDFPSFXFNAME : array of filename of surfex climatology file on target geometries

!        IMPLICIT ARGUMENTS
!        --------------------
!        See #include below.

!     METHOD.
!     -------
!        See documentation about FULL-POS. 

!     EXTERNALS.
!     ----------
!        Calls SUAFN1 and SUAFN2.
!        Is called by SU0YOMA.

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib 17-Aug-2016 No printouts if no Fullpos
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMAFN   , ONLY : TAFN
USE YOMFPC   , ONLY : TNAMFPSCI
USE YOMCOMPO , ONLY : TCOMPO
USE YOEAERATM, ONLY : TEAERATM

IMPLICIT NONE

TYPE(TNAMFPSCI),    INTENT(IN)    :: YDNAMFPSCI
TYPE(TAFN)         ,INTENT(OUT)   :: YDAFN
LOGICAL            ,INTENT(IN)    :: LDNEWCL(2)
TYPE(TYPE_GFLD)    ,INTENT(IN)    :: YGFL
LOGICAL            ,INTENT(IN)    :: LDNHDYN
INTEGER(KIND=JPIM), INTENT(IN)    :: KFPCONF
TYPE(TCOMPO)       ,INTENT(IN), OPTIONAL    :: YDCOMPO
TYPE(TEAERATM)     ,INTENT(IN), OPTIONAL    :: YDEAERATM
!Phasing - comment:
!YDEAERATM was made OPTIONAL for the phasing of the LAM ->
!Temporary solution - the correct way would be to introduce YDMODEL
!into ebicli.F90 or to make a new suafn.F90 for the LAM
!This probably applies to YDCOMPO too, which is now needed
!here for IFS with atmospheric composition enabled.
CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: CDFPSFXFNAME(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "suafn1.intfb.h"
#include "suafn2.intfb.h"
#include "suafn3.intfb.h"

IF (LHOOK) CALL DR_HOOK('SUAFN',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1. SET DEFAULT VALUES
!           -------------------

CALL SUAFN1(YDNAMFPSCI,YDAFN,LDNEWCL,YGFL,KFPCONF,YDCOMPO,YDEAERATM,CDFPSFXFNAME)

!*       3. SET EQUIVALENCES
!           ----------------

CALL SUAFN2(LDNHDYN,YDAFN,KFPCONF)

!*       4. PRINT OUT FINAL VALUES
!           ----------------------

CALL SUAFN3(YDAFN)

IF (LHOOK) CALL DR_HOOK('SUAFN',1,ZHOOK_HANDLE)
END SUBROUTINE SUAFN
