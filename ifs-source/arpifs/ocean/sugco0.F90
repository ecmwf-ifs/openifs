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

SUBROUTINE SUGCO0(YDMCC)

!**** *SUGCO0 * - ROUTINE TO RESET TO 0 THE COUPLED FIELDS

!     PURPOSE.
!     --------
!        RESET YOMGCO ARRAY TO 0

!**   INTERFACE.
!     ----------
!        *CALL* *SUGCO0*

!     EXPLICIT ARGUMENTS :  NONE
!     --------------------

!     IMPLICIT ARGUMENTS :
!     --------------------
!        COMMON  YOMGCO

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "IN CORE MODEL"

!     AUTHOR.
!     -------
!      C. DREVETON D'APRES SUGFU0 (25.05.92)

!     Modifications:
!     --------------
!      Modified 03-04-22 by P. Marquet (Cloud Forcing diag. in CUMCODM)
!     F. Vana  05-Mar-2015  Support for single precision
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB, JPIM
USE PARKIND_OCEAN, ONLY : JPRO
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMGCO   , ONLY : STRSU    ,STRSV    ,FCHAS    ,FRSOS    ,&
 & FHUMS    ,RUIST    ,FCHLL    ,FCHLN    ,FCHSS    ,&
 & FHUML    ,FHUMN    ,FRLDS    ,QWPRO    ,TMAX2    ,&
 & TMIN2    ,HSNOW    ,TSUR2    ,TSTS     ,TSFL     ,&
 & XPUQ     ,XPVQ     ,FCRFTH   ,FCRFSO

USE YOMMCC   , ONLY : TMCC

USE CPLNG    , ONLY : CPL_OUT

!      ----------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM) :: JF

!      ----------------------------------------------------------------

!*       1.    SET YOMGCO TO 0.
!              ---------------

TYPE(TMCC),INTENT(INOUT):: YDMCC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SUGCO0',0,ZHOOK_HANDLE)
ASSOCIATE(LMCCDYNSEAICE=>YDMCC%LMCCDYNSEAICE, &
 & LNEMOOCEICEMIX=>YDMCC%LNEMOOCEICEMIX, &
 & LNEMOATMFLDS=>YDMCC%LNEMOATMFLDS, &
 & CPLNG_FLD=>YDMCC%CPLNG_FLD, &
 & CPLNG_NUM_FIELDS=>YDMCC%CPLNG_NUM_FIELDS)

FCHAS=0._JPRB
FRSOS=0._JPRB
FHUMS=0._JPRB
FCHLL=0._JPRB
FCHLN=0._JPRB
FCHSS=0._JPRB

STRSU=0._JPRB
STRSV=0._JPRB
RUIST=0._JPRB
FHUML=0._JPRB
FHUMN=0._JPRB
FRLDS=0._JPRB
QWPRO=0._JPRB
TMAX2=-1000._JPRB
TMIN2=1000._JPRB
HSNOW=0._JPRB
TSUR2=0._JPRB
TSTS =0._JPRB
TSFL =0._JPRB
XPUQ =0._JPRB
XPVQ =0._JPRB
FCRFTH=0._JPRB
FCRFSO=0._JPRB

IF(LMCCDYNSEAICE) THEN
  ! Reset to zero
   DO JF=1,CPLNG_NUM_FIELDS
      IF (CPLNG_FLD(JF)%INOUT==CPL_OUT) THEN
         CPLNG_FLD(JF)%D(:,:,:)=0.0_JPRO
      ENDIF
   ENDDO
ENDIF

IF (LNEMOATMFLDS) THEN
   CPLNG_FLD(YDMCC%IP_A_SST_ATM)%D(:,:,:)=0.0_JPRO
   CPLNG_FLD(YDMCC%IP_A_TSK_ATM)%D(:,:,:)=0.0_JPRO
ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGCO0',1,ZHOOK_HANDLE)
END SUBROUTINE SUGCO0
