MODULE SRFWLS_MOD
CONTAINS
SUBROUTINE SRFWLS(KIDIA,KFDIA,KLON,&
 & PTMST,PWLM1M,PCVL,PCVH,PWLMX,&
 & PFRTI,PEVAPTI,PRSFC,PRSFL,&
 & LDLAND,&
 & YDSOIL,YDVEG,&
 & PTSFC,PTSFL)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_THF   , ONLY : RHOH2O
USE YOS_SOIL  , ONLY : TSOIL
USE YOS_VEG   , ONLY : TVEG

#ifdef DOC

! (C) Copyright 2011- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFWLS* - COMPUTES CHANGES IN THE SKIN RESERVOIR.
!     PURPOSE.
!     --------
!          THIS ROUTINE COMPUTES THE CHANGES IN THE SKIN RESERVOIR AND
!     THE RUN-OFF BEFORE THE SNOW MELTS.

!**   INTERFACE.
!     ----------
!          *SRFWLS* IS CALLED FROM *SURFS*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----

!     INPUT PARAMETERS (INTEGER):
!    *KIDIA*      START POINT
!    *KFDIA*      END POINT
!    *KLON*       NUMBER OF GRID POINTS PER PACKET

!     INPUT PARAMETERS (REAL):
!    *PTMST*      TIME STEP                                      S

!     INPUT PARAMETERS (LOGICAL):
!    *LDLAND*     LAND/SEA MASK (TRUE/FALSE)

!     INPUT PARAMETERS AT T-1 OR CONSTANT IN TIME (REAL):
!    *PWLM1M*     SKIN RESERVOIR WATER CONTENT                   kg/m**2
!    *PCVL*       LOW VEGETATION COVER  (CORRECTED)              (0-1)
!    *PCVH*       HIGH VEGETATION COVER (CORRECTED)              (0-1)
!    *PWLMX*      MAXIMUM SKIN RESERVOIR CAPACITY                kg/m**2
!    *PFRTI*      TILE FRACTIONS                                 (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
!    *PEVAPTI*      SURFACE MOISTURE FLUX, FOR EACH TILE       KG/M2/S
!    *PRSFC*      CONVECTIVE RAIN FLUX AT THE SURFACE          KG/M**2/S
!    *PRSFL*      LARGE SCALE RAIN FLUX AT THE SURFACE         KG/M**2/S

!     OUTPUT PARAMETERS AT T+1 (UNFILTERED,REAL):
!    *PTSFC*      CONVECTIVE THROUGHFALL AT THE SURFACE        KG/M**2/S
!    *PTSFL*      LARGE SCALE THROUGHFALL AT THE SURFACE       KG/M**2/S
!                  (NB: THROUGHFALL=RAINFALL-INTERCEPTION)

!     METHOD.
!     -------
!          STRAIGHTFORWARD ONCE THE DEFINITION OF THE CONSTANTS IS
!     UNDERSTOOD. FOR THIS REFER TO DOCUMENTATION.

!     EXTERNALS.
!     ----------
!          NONE.

!     REFERENCE.
!     ----------
!          SEE SOIL PROCESSES' PART OF THE MODEL'S DOCUMENTATION FOR
!     DETAILS ABOUT THE MATHEMATICS OF THIS ROUTINE.

!     Original   
!     --------
!          Simplified version based on SRFWL
!     M. Janiskova              E.C.M.W.F.     27-07-2011  

!     Modifications
!     -------------

!     ------------------------------------------------------------------
#endif

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)   :: KLON

REAL(KIND=JPRB),    INTENT(IN)   :: PTMST
REAL(KIND=JPRB),    INTENT(IN)   :: PWLM1M(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVL(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PCVH(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PWLMX(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PFRTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PEVAPTI(:,:)
REAL(KIND=JPRB),    INTENT(IN)   :: PRSFC(:)
REAL(KIND=JPRB),    INTENT(IN)   :: PRSFL(:)

LOGICAL,            INTENT(IN)   :: LDLAND(:)

TYPE(TSOIL),        INTENT(IN)   :: YDSOIL
TYPE(TVEG),         INTENT(IN)   :: YDVEG

REAL(KIND=JPRB),    INTENT(OUT)  :: PTSFC(:)
REAL(KIND=JPRB),    INTENT(OUT)  :: PTSFL(:)

INTEGER(KIND=JPIM) :: JL

REAL(KIND=JPRB) :: ZPWL(KLON)
REAL(KIND=JPRB) :: ZEPPRCP,&
 & ZEPTINY, ZIPRCP, ZMPRCP, ZPSFR, ZQHFLW, ZTMST, &
 & ZTPRCP, ZVINTER, ZWL


REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SRFWLS_MOD:SRFWLS',0,ZHOOK_HANDLE)
ASSOCIATE(RPSFR=>YDSOIL%RPSFR, &
 & RVINTER=>YDVEG%RVINTER)

!*         1.    SET UP SOME CONSTANTS.
!                --- -- ---- ----------
!*               PHYSICAL CONSTANTS.

ZVINTER=RVINTER
ZPSFR=1.0_JPRB/RPSFR

!*    SECURITY PARAMETERS
ZEPTINY=10._JPRB*TINY(RHOH2O)
ZEPPRCP=ZEPTINY

!*    COMPUTATIONAL CONSTANTS.
ZTMST=1.0_JPRB/PTMST

!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
!*          3.  UPWARDS EVAPORATION.
!               ------- ------------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN

!           INITIALISE PWL (TO MAKE THE CODE SIMPLER).
    ZPWL(JL)=PWLM1M(JL)

!           EVAPORATION OF THE SKIN RESERVOIR (EL < 0).
    IF (PEVAPTI(JL,3) < 0.0_JPRB) THEN
      ZWL=PWLM1M(JL)
      ZQHFLW=PTMST*PFRTI(JL,3)*PEVAPTI(JL,3)
      ZPWL(JL)=ZWL+ZQHFLW
      ZPWL(JL)=MAX(0.0_JPRB,ZPWL(JL))
    ENDIF
  ELSE
!    Sea points
    ZPWL(JL)=0.0_JPRB
  ENDIF
ENDDO

!           6.  INTERCEPTION OF PRECIPITATION BY THE VEGETATION.
!               ------------ -- ------------- -- --- -----------
DO JL=KIDIA,KFDIA
  IF (LDLAND(JL)) THEN
!          LARGE SCALE PRECIPITATION.
    IF (PRSFL(JL) > ZEPPRCP) THEN
      ZTPRCP=PRSFL(JL)
      ZIPRCP=ZTPRCP*(PCVL(JL)+PCVH(JL))*ZVINTER
      ZMPRCP=MIN(PWLMX(JL)-ZEPTINY-ZPWL(JL),PTMST*ZIPRCP)
      PTSFL(JL)=PRSFL(JL)-ZMPRCP*ZTMST
    ELSE
      PTSFL(JL)=0.0_JPRB
    ENDIF

!          CONVECTIVE PRECIPITATION.
    IF (PRSFC(JL) > ZEPPRCP) THEN
      ZTPRCP=PRSFC(JL)*RPSFR
      ZIPRCP=ZTPRCP*(PCVL(JL)+PCVH(JL))*ZVINTER
      ZMPRCP=MIN(PWLMX(JL)-ZEPTINY-ZPWL(JL),PTMST*ZIPRCP)*ZPSFR
      PTSFC(JL)=PRSFC(JL)-ZMPRCP*ZTMST
    ELSE
      PTSFC(JL)=0.0_JPRB
    ENDIF

!          SEA POINTS.
  ELSE
    PTSFC(JL)=PRSFC(JL)
    PTSFL(JL)=PRSFL(JL)
  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRFWLS_MOD:SRFWLS',1,ZHOOK_HANDLE)

END SUBROUTINE SRFWLS
END MODULE SRFWLS_MOD




