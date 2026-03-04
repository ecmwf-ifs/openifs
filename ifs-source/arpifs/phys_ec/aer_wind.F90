! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_WIND ( &
 & KIDIA , KFDIA , KLON  , KTILES, &
 & PUMLEV, PVMLEV, PTMLEV, PQMLEV, PAPHMS, PGEOMLEV, &
 & PUSTRTI,PVSTRTI,PAHFSTI,PEVAPTI,PFRTI , PZ0M, &
 & PWS    ,PGUST  ,PUST)  
!     ------------------------------------------------------------------

!**   *AER_WIND* - Prepares wind, gust and friction velocity for the 
!                 aerosol/surface interaction

!     Original  A. BELJAARS   ECMWF    24/08/2007.
!     Modified 

!     PURPOSE
!     -------

!     Compute velocity scales for aerosol modelling before VDF has been called

!     INTERFACE
!     ---------

!     *AER_WIND* IS CALLED BY ?

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        Start index
!     *KFDIA*        End index
!     *KLON*         Number of grid points per packet
!     *KTILES*       Number of tiles

!     INPUT PARAMETERS (REAL):

!     *PUMLEV*       X-velocity at t-1, lowest model level      (m/s)
!     *PVMLEV*       Y-velocity at t-1, lowest model level      (m/s)
!     *PTMLEV*       Geopotential at t-1                        (K)
!     *PQMLEV*       Geopotential at t-1                        (kg/kg)
!     *PAPHMS*       Surface pressure at t-1                    (Pa)
!     *PGEOMLEV*     Geopotentail at -1, lowest model level     (m2/s2)
!     *PUSTRTI*      X-stress                                   (N/m2)
!     *PVSTRTI*      Y-stress                                   (N/m2)
!     *PAHFSTI*      sensible heat flux                         (W/m2)
!     *PEVAPTI*      moisture flux                              (kg/m2s)
!     *PFRTI*        tile fractions                             (-)
!     *PZ0M*         roughness length for momentum              (m)

!     OUTPUT PARAMETERS (REAL):

!     *PWS*          Wind speed (average of horizontal wind speed)   (m/s)
!     *PGUST         Wind gust (maximum 3 second gust in the hour)   (m/s)
!     *PUST*         Friction velocity                               (m/s)       


!     METHOD
!     ------

!     SEE DOCUMENTATION

!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM,JPRB
USE YOMHOOK  , ONLY : LHOOK,DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG,RD,RCPD,RETV
USE PARPHY   , ONLY : REPDU2,RKAP
USE YOMCT3   , ONLY : NSTEP

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUMLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVMLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQMLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHMS(KLON) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGEOMLEV(KLON) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVSTRTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PAHFSTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVAPTI(KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI  (KLON,KTILES) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PZ0M(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PWS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGUST(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PUST(KLON)

!*            LOCAL STORAGE
!             ----- -------

REAL(KIND=JPRB) ::    ZAHFSM(KLON),ZEVAPM(KLON),ZUSTRM(KLON),ZVSTRM(KLON),&
                     &ZDU2(KLON),ZRHO(KLON),ZBUOM(KLON),ZUSTAR(KLON)

INTEGER(KIND=JPIM) :: JL,JTILE

REAL(KIND=JPRB) ::    ZDUA,ZZCDN,ZROWT,ZROWQ,ZWST2,ZCON2,ZIPBL,&
                     &ZDUA2,ZUGN,ZCZI,Z1D3,ZIDL,ZEPUST
REAL(KIND=JPHOOK) ::    ZHOOK_HANDLE

!     ------------------------------------------------------------------


!*       1.     Initialize constants and arrays
!               ---------- --------- --- ------

IF (LHOOK) CALL DR_HOOK('AER_WIND',0,ZHOOK_HANDLE)
ZIPBL = 1000._JPRB
ZCON2 = 2.0_JPRB/3._JPRB
ZUGN=7.71_JPRB
ZCZI=0.5_JPRB/12._JPRB
Z1D3=1._JPRB/3._JPRB
ZEPUST=0.0001_JPRB

DO JL=KIDIA,KFDIA
  ZRHO(JL)=PAPHMS(JL)/( RD*PTMLEV(JL)*(1.0_JPRB+RETV*PQMLEV(JL)) )
  ZDU2(JL)=MAX(REPDU2,PUMLEV(JL)**2+PVMLEV(JL)**2)
ENDDO


!*       2.     Time step 0 does not have fluxes yet (neutral estimate)
!               ---- ---- - ---- --- ---- ------ --- -------- ---------

IF (NSTEP == 0) THEN
  DO JTILE=1,KTILES
    DO JL=KIDIA,KFDIA
      ZDUA=SQRT(ZDU2(JL))
      ZZCDN=(RKAP/LOG(1.0_JPRB+PGEOMLEV(JL)/(RG*PZ0M(JL))))**2
      PUSTRTI(JL,JTILE)=ZRHO(JL)*PUMLEV(JL)*ZDUA*ZZCDN
      PVSTRTI(JL,JTILE)=ZRHO(JL)*PVMLEV(JL)*ZDUA*ZZCDN
      PAHFSTI(JL,JTILE)=0.0_JPRB
      PEVAPTI(JL,JTILE)=0.0_JPRB
    ENDDO
  ENDDO
ENDIF


!        2. Compute area averaged fluxes
!           ------- ---- -------- ------ 

ZAHFSM(KIDIA:KFDIA) = 0.0_JPRB
ZEVAPM(KIDIA:KFDIA) = 0.0_JPRB
ZUSTRM(KIDIA:KFDIA) = 0.0_JPRB
ZVSTRM(KIDIA:KFDIA) = 0.0_JPRB
DO JTILE=1,KTILES
  DO JL=KIDIA,KFDIA
    ZAHFSM(JL)=ZAHFSM(JL)+PFRTI(JL,JTILE)*PAHFSTI(JL,JTILE)
    ZEVAPM(JL)=ZEVAPM(JL)+PFRTI(JL,JTILE)*PEVAPTI(JL,JTILE)
    ZUSTRM(JL)=ZUSTRM(JL)+PFRTI(JL,JTILE)*PUSTRTI(JL,JTILE)
    ZVSTRM(JL)=ZVSTRM(JL)+PFRTI(JL,JTILE)*PVSTRTI(JL,JTILE)
  ENDDO
ENDDO


!        2. Mean horizontal wind speed, friction velocity (incl gustiness effects)
!           ---------------------------------------------------------------------- 

DO JL=KIDIA,KFDIA
  ZROWQ=ZEVAPM(JL)
  ZROWT=ZAHFSM(JL)/RCPD
  ZBUOM(JL)=RG*(-RETV*ZROWQ-ZROWT/PTMLEV(JL))/ZRHO(JL)
  ZUSTAR(JL)=MAX(ZEPUST,SQRT(SQRT(ZUSTRM(JL)**2+ZVSTRM(JL)**2)/ZRHO(JL)))
  IF (ZBUOM(JL) <= 0.0_JPRB) THEN
    PUST(JL)=ZUSTAR(JL)
    ZDUA2=ZDU2(JL)
  ELSE
    ZWST2=(ZBUOM(JL)*ZIPBL)**ZCON2
    ZDUA2=ZDU2(JL)+ZWST2
    PUST(JL)=ZUSTAR(JL)*SQRT(SQRT(ZDUA2/ZDU2(JL)))
  ENDIF
  PWS(JL) =SQRT(ZDUA2)
ENDDO


!        3. wind gusts (3s extreme in 1 hour)
!           ---- ----- ----------------------
DO JL=KIDIA,KFDIA

!     AREA AVERAGE OF ABSOLUTE 10 M (TO BE USED FOR GUSTS)

  ZIDL=-ZIPBL*RKAP*ZBUOM(JL)/PUST(JL)**3

  IF (ZIDL >= 0.) THEN
    PGUST(JL)=PWS(JL)+ZUSTAR(JL)*ZUGN
  ELSE
    PGUST(JL)=PWS(JL)+ZUSTAR(JL)*ZUGN*(1.0_JPRB-ZCZI*ZIDL)**Z1D3
  ENDIF
  
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_WIND',1,ZHOOK_HANDLE)
END SUBROUTINE AER_WIND
