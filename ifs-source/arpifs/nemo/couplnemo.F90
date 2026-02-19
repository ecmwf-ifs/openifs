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

SUBROUTINE COUPLNEMO(YDGEM,YDMCC,YDRIP,KSTEP)
!
!**** *COUPLNEMO*  - Run the NEMO model.
!
!     Purpose.
!     --------
!       Update the NEMO ocean state by calling the time
!       stepping routine of NEMO until it is uptodate.
!
!**   Interface.
!     ----------
!       *CALL*  *COUPLNEMO(KSTEP)*
!
!     Input:
!     -----
!
!     Output:
!     ------
!
!     Method:
!     ------
!       
!     Externals:ifs/nemo/couplnemo.F90
!     ---------
!
!     Reference:
!     ---------
!
!     Author:
!     -------
!       K. Mogensen, ECMWF
!
!     Modifications.
!     --------------
!       S. Keeley + K. Mogensen Januar 2012 update to support LIM2.
!       T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!       K. Yessad (July 2014): Move some variables.
!     F. Vana  05-Mar-2015  Support for single precision
!     -----------------------------------------------------------
   
USE YOMGEM   , ONLY : TGEM
USE PARKIND1 , ONLY : JPIM, JPRB
USE PARKIND_OCEAN, ONLY : JPRO
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE MPL_MODULE, ONLY : MPL_COMM
USE YOMCT0   , ONLY : NFRCO
USE YOMMP0   , ONLY : MYPROC, NPROC
USE YOMRIP   , ONLY : TRIP
USE YOMGCO   , ONLY : STRSU, STRSV, FCHAS, FRSOS, FHUMS
USE YOMMCC   , ONLY : TMCC
! NEMO time step control
USE YOMNEMO  , ONLY : NEMOCSTEP, NEMONSTEP
! Coupling interface
USE CPLNG

!     -----------------------------------------------------------

IMPLICIT NONE

TYPE(TGEM)        , INTENT(IN) :: YDGEM
TYPE(TMCC)         ,INTENT(INOUT):: YDMCC
TYPE(TRIP)         ,INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM), INTENT(IN) :: KSTEP 

!     -----------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRO) :: ZRUNOFF(YDGEM%NGPTOT)
REAL(KIND=JPRO) :: ZOCERUNOFF(YDGEM%NGPTOT)
REAL(KIND=JPRO) :: Z2DECV(YDGEM%NGPTOT)
REAL(KIND=JPRO) :: ZFACT
INTEGER(KIND=JPIM) :: JSTPNEMO
INTEGER(KIND=JPIM) :: JF
INTEGER(KIND=JPIM) :: IIDATE, IITIME
INTEGER(KIND=JPIM) :: IL_QS, IL_QNS
!     -----------------------------------------------------------

#include "sugco0.intfb.h"

!     -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPLNEMO',0,ZHOOK_HANDLE)
ASSOCIATE(NGPTOT=>YDGEM%NGPTOT, &
 & LMCCDYNSEAICE=>YDMCC%LMCCDYNSEAICE, &
 & LNEMOFLUXNC=>YDMCC%LNEMOFLUXNC, &
 & TSTEP=>YDRIP%TSTEP,LNEMOATMFLDS=>YDMCC%LNEMOATMFLDS,&
 & LNEMOOCEICEMIX=>YDMCC%LNEMOOCEICEMIX,CPLNG_FLD=>YDMCC%CPLNG_FLD,&
 & LNEMOQNSICEFILT=>YDMCC%LNEMOQNSICEFILT, &
 & LNEMOACCUMFLUX=>YDMCC%LNEMOACCUMFLUX, &
 & L2DECV2NEMO=>YDMCC%L2DECV2NEMO, &
 & CPLNG_NUM_FIELDS=>YDMCC%CPLNG_NUM_FIELDS )
!     -----------------------------------------------------------

! Normalize the data

IF (LMCCDYNSEAICE) THEN

   IF (LNEMOOCEICEMIX) THEN
      IL_QS=YDMCC%IP_A_QS_MIX
      IL_QNS=YDMCC%IP_A_QNS_MIX
   ELSE
      IL_QS=YDMCC%IP_A_QS_OCE
      IL_QNS=YDMCC%IP_A_QNS_OCE
   ENDIF

   IF (LNEMOACCUMFLUX) THEN
      ZFACT=1.0_JPRO/NFRCO
      DO JF=1,CPLNG_NUM_FIELDS
         IF (CPLNG_FLD(JF)%INOUT==CPL_OUT) THEN
            CPLNG_FLD(JF)%D(:,:,:)=&
               & CPLNG_FLD(JF)%D(:,:,:)*ZFACT
         ENDIF
      ENDDO
   ENDIF
   ! Something to look at later!!!
   ZRUNOFF(:)=0.0_JPRO
   ZOCERUNOFF(:)=0.0_JPRO
   ! We can't pass CPLNG_FLD(YDMCC%IP_A_2DECV_SKT)%D(:,1,1)
   ! directly to NEMO since it might not be alllocated
   IF (L2DECV2NEMO) THEN
      Z2DECV(:)=CPLNG_FLD(YDMCC%IP_A_2DECV_SKT)%D(:,1,1)
   ELSE
      Z2DECV(:)=0.0_JPRO
   ENDIF

   ! Update fluxes in the NEMO interface

#ifdef WITH_NEMO
   CALL NEMOGCMCOUP_LIM2_UPDATE( MYPROC-1, NPROC, MPL_COMM, NGPTOT, &
      &                          CPLNG_FLD(YDMCC%IP_A_TAUX_OCE)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_TAUY_OCE)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_TAUX_ICE)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_TAUY_ICE)%D(:,1,1), &
      &                          CPLNG_FLD(IL_QS)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_QS_ICE)%D(:,1,1), &
      &                          CPLNG_FLD(IL_QNS)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_QNS_ICE)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_DQNS_DT)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_EVAP_TOTAL)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_EVAP_ICE)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_PRECIP_LIQUID)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_PRECIP_SOLID)%D(:,1,1), &
      &                          ZRUNOFF, ZOCERUNOFF, & 
      &                          CPLNG_FLD(YDMCC%IP_A_TOTAL_CC)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_LOW_CC)%D(:,1,1), &
      &                          CPLNG_FLD(YDMCC%IP_A_IST_ATM)%D(:,1,1), &
      &                          Z2DECV, &
      &                          KSTEP, LNEMOFLUXNC, LNEMOOCEICEMIX, &
      &                          LNEMOQNSICEFILT, L2DECV2NEMO )
#else
   CALL ABOR1('COUPLNEMO: COMPILED WITHOUT WITH_NEMO')
#endif

ELSE

  ZFACT=1.0_JPRO/(NFRCO*REAL(TSTEP,JPRO))
  CPLNG_FLD(YDMCC%IP_A_TAUX)%D(:,1,1)=REAL(STRSU(:)*ZFACT,JPRO)
  CPLNG_FLD(YDMCC%IP_A_TAUY)%D(:,1,1)=REAL(STRSV(:)*ZFACT,JPRO)
  CPLNG_FLD(YDMCC%IP_A_QS)%D(:,1,1)=REAL(FRSOS(:,1)*ZFACT,JPRO)
  CPLNG_FLD(YDMCC%IP_A_QNS)%D(:,1,1)=REAL(FCHAS(:,1)*ZFACT,JPRO)
  CPLNG_FLD(YDMCC%IP_A_WATER)%D(:,1,1)=REAL(FHUMS(:,1)*ZFACT,JPRO)

  ! Update fluxes in the NEMO interface

#ifdef WITH_NEMO
  CALL NEMOGCMCOUP_UPDATE( MYPROC-1, NPROC, MPL_COMM, NGPTOT, &
     &                     CPLNG_FLD(YDMCC%IP_A_TAUX)%D(:,1,1), & 
     &                     CPLNG_FLD(YDMCC%IP_A_TAUY)%D(:,1,1), &
     &                     CPLNG_FLD(YDMCC%IP_A_QS)%D(:,1,1), &
     &                     CPLNG_FLD(YDMCC%IP_A_QNS)%D(:,1,1), &
     &                     CPLNG_FLD(YDMCC%IP_A_WATER)%D(:,1,1), &
     &                     KSTEP, LNEMOFLUXNC )
#else
  CALL ABOR1('COUPLNEMO: COMPILED WITHOUT WITH_NEMO')
#endif

ENDIF

#ifdef WITH_NEMO
IF(LNEMOATMFLDS) THEN
   CALL NEMOGCMCOUP_UPDATE_ADD( MYPROC-1, NPROC, MPL_COMM, NGPTOT, &
      & CPLNG_FLD(YDMCC%IP_A_SST_ATM)%D(:,1,1)-273.15_JPRO, &
      & CPLNG_FLD(YDMCC%IP_A_TSK_ATM)%D(:,1,1)-273.15_JPRO, &
      & KSTEP, LNEMOFLUXNC )
ENDIF
#endif

! Run the NEMO model for NEMONSTEP time steps

DO JSTPNEMO=NEMOCSTEP,NEMOCSTEP+NEMONSTEP-1

   WRITE(NULOUT,*)'CALLING NEMO FOR STEP : ', JSTPNEMO

   ! Advance the NEMO model 1 time step
#ifdef WITH_NEMO
   CALL NEMOGCMCOUP_STEP( JSTPNEMO, IIDATE, IITIME )
#else
  CALL ABOR1('COUPLNEMO: COMPILED WITHOUT WITH_NEMO')
#endif

ENDDO

NEMOCSTEP=NEMOCSTEP+NEMONSTEP

! Reset to zero
CALL SUGCO0(YDMCC)

!     -----------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('COUPLNEMO',1,ZHOOK_HANDLE)
END SUBROUTINE COUPLNEMO
