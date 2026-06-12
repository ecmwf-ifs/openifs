MODULE SRFSN_ASN_MOD
CONTAINS
SUBROUTINE SRFSN_ASN(KIDIA,KFDIA,KLON,PTMST,LLNOSNOW,PASNM1M,&
  & PMSN,PTSNM1M,PSNOWF,YDSOIL,YDCST,PASN)


USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOS_SOIL , ONLY : TSOIL 
USE YOS_CST  , ONLY : TCST
! (C) Copyright 2015- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFSN_RSN* - Snow albedo
!     PURPOSE.
!     --------
!          THIS ROUTINE CONTROLS THE ALBEDO EVOLUTION

!**   INTERFACE.
!     ----------
!          *SRFSN_ASN* IS CALLED FROM *SRFSN_DRIVER*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET


!     INPUT PARAMETERS (REAL):
!     *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LLNOSNOW*   NO-SNOW/SNOW MASK (TRUE IF NO-SNOW)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PASNM1M*    SNOW ALBEDO                                      0-1
!    *PTSNM1M*    TEMPERATURE OF SNOW LAYER                        K
!    *PSNOWF*     TOTAL SNOW FLUX AT THE SURFACE                 KG/M**2/S
!    *PMSN*       DIFFERENCE BETWEEN LATENT HEAT OF MELTING      J/m**2
!                 AND LATENT HEAT OF FREEZING. PROXY FOR MELTING COND 

!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PASN*       SNOW ALBEDO                                      0-1
!     

!     METHOD.
!     -------
!          

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          

!     Modifications:
!     Original   E. Dutra      ECMWF     04/12/2015
!                G. Arduini    ECMWF     01/09/2021

!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments 
INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON
REAL(KIND=JPRB)   , INTENT(IN)   :: PTMST
LOGICAL           , INTENT(IN)   :: LLNOSNOW(:) 

REAL(KIND=JPRB)   , INTENT(IN)   :: PASNM1M(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PMSN(:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PTSNM1M(:,:)
REAL(KIND=JPRB)   , INTENT(IN)   :: PSNOWF(:)
TYPE(TSOIL)       , INTENT(IN)   :: YDSOIL
TYPE(TCST)        , INTENT(IN)   :: YDCST


REAL(KIND=JPRB)   , INTENT(OUT)  :: PASN(:)

! Local variables 
INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_ASN_MOD:SRFSN_ASN',0,ZHOOK_HANDLE)

!    -----------------------------------------------------------------

DO JL=KIDIA,KFDIA
  IF (LLNOSNOW(JL)) THEN
    PASN(JL)= YDSOIL%RALFMAXSN
  ELSE
!     IF (PMSN(JL) > 0.0_JPRB .OR. PTSNM1M(JL,1) > YDCST%RTT-2._JPRB ) THEN
    IF (PMSN(JL) > 0.0_JPRB  ) THEN
    ! MELTING CONDITIONS
      PASN(JL)=(PASNM1M(JL)-YDSOIL%RALFMINSN)*&
       & EXP(-YDSOIL%RTAUF*PTMST/YDCST%RDAY)+&
       & YDSOIL%RALFMINSN
    ELSE
    ! NORMAL CONDITIONS
      PASN(JL)=MAX(YDSOIL%RALFMINSN,PASNM1M(JL)-YDSOIL%RTAUA*PTMST/YDCST%RDAY)
    ENDIF
    ! UPDATE ALBEDO DUE TO SNOWFALL EVENTS
    PASN(JL)=PASN(JL)+ MIN(MAX(PSNOWF(JL)*PTMST , 0._JPRB)/(10._JPRB) , 1._JPRB) *&
     & (YDSOIL%RALFMAXSN-PASN(JL))
    PASN(JL)=MIN(YDSOIL%RALFMAXSN,MAX(PASN(JL),YDSOIL%RALFMINSN))
    
  ENDIF
ENDDO



!    -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SRFSN_ASN_MOD:SRFSN_ASN',1,ZHOOK_HANDLE)

END SUBROUTINE SRFSN_ASN
END MODULE SRFSN_ASN_MOD
