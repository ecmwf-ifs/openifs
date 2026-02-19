! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_STRATOLOSS(YGFL,PY0,PY,PDT)

!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!   Simple solver based on Euler Backward for trace gases defined in troposphere
!   but not available in stratospheric chemistry
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_stratoloss* IS CALLED FROM *CHEM_bascoetm5*.

! INPUTS:
! -------
! PDT   :  Time step in seconds 
! PY0(NCHEM)        : initial volume ratios OF TRACERS           (mol/mol)
!
!
! OUTPUTS:
! -------
! PY (NCHEM)        : final   volume ratios OF TRACERS           (mol/mol)
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-03-12



USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE BASCOETM5_MODULE, ONLY :  IPAR, IETH,  IOLE,  IALD2  , IPAN,&
   &   IORGNTR,  ISO2, ISO4, IMGLY, IRN222,  IPB210,   ICH3OH,&
   &   IHCOOH,  IMCOOH,   IC2H6,&
   &   IETHOH,   IC3H8,   IC3H6,  IACET


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD),INTENT(INOUT) :: YGFL
REAL(KIND=JPRB),INTENT(IN)    :: PDT
REAL(KIND=JPRB),INTENT(IN)    :: PY0(YGFL%NCHEM)   ! initial concentrations
REAL(KIND=JPRB),INTENT(INOUT) :: PY(YGFL%NCHEM)    ! final concentrations

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! Rate for radon: ~ 5.5 day
REAL(KIND=JPRB), PARAMETER ::  ZRR_RN222=2.10E-6_JPRB 
! Rate for CB05: ~ 10 day (relevant for lower stratosphere?)
REAL(KIND=JPRB), PARAMETER ::  ZRRSTRAT=1.157E-6_JPRB 
! Rate for SO4: ~ 31 day (relevant for SO4?)
REAL(KIND=JPRB), PARAMETER ::  ZRRSTRAT_MONTH=3.73E-7_JPRB 
REAL(KIND=JPRB)            :: ZTSCALI

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_STRATOLOSS',0,ZHOOK_HANDLE )

     !* Computation of denomenator, only once:
     ZTSCALI =  1./(1._JPRB+ZRRSTRAT*PDT)
     PY(ICH3OH)=PY0(ICH3OH)*ZTSCALI
     PY(IHCOOH)=PY0(IHCOOH)*ZTSCALI
     PY(IMCOOH)=PY0(IMCOOH)*ZTSCALI
     PY(IETHOH)=PY0(IETHOH)*ZTSCALI
     PY(IC2H6)=PY0(IC2H6)*ZTSCALI
     PY(IC3H8)=PY0(IC3H8)*ZTSCALI
     PY(IC3H6)=PY0(IC3H6)*ZTSCALI
     PY(IPAR)=PY0(IPAR)*ZTSCALI
     PY(IETH)=PY0(IETH)*ZTSCALI
     PY(IOLE)=PY0(IOLE)*ZTSCALI
     PY(IALD2)=PY0(IALD2)*ZTSCALI
     PY(IMGLY)=PY0(IMGLY)*ZTSCALI
     PY(IORGNTR)=PY0(IORGNTR)*ZTSCALI
     PY(IPAN)=PY0(IPAN)*ZTSCALI
!     PY(ISO2)=PY0(ISO2)*ZTSCALI
! Stratospheric SO4 blows up. Here introduce a slow decay rate (month-time scale) to prevent this.
! In reality there is only SO4 sedimentation     
!     PY(ISO4)=PY0(ISO4)/(1._JPRB+ZRRSTRAT_MONTH*PDT)
     PY(IACET)=PY0(IACET)*ZTSCALI

     
     PY(IRN222) = PY0(IRN222)/(1._JPRB+ZRR_RN222*PDT)
     PY(IPB210) = PY0(IPB210)+PY0(IRN222)-PY(IRN222)


IF (LHOOK) CALL DR_HOOK('TM5_STRATOLOSS',1,ZHOOK_HANDLE )
END SUBROUTINE TM5_STRATOLOSS
