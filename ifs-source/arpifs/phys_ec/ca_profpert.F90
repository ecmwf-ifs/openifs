! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE CA_PROFPERT(YDECUCONVCA,KLON,KLEV,KIDIA,KFDIA,PTSPHY,PQP,PTENT,PTENQ,PAPRSF,PCUCONVCA,PNLCONVCA)

!**** *CA_PROFPERT* - Perturb T and q profiles where CA is active

!     Purpose.
!     --------
!        Perturb the T and q input profiles to the convection scheme 
!        The result then store in the specific tendency used in convection


!**   Interface.
!     ----------
!        *CALL* *CA_PROFPERT(KLON,KLEV,KIDIA,KFDIA,PTSPHY,PQP,PTENT,PTENQ,PAPRSF,PCUCONVCA,PNLCONVCA)

!        Explicit arguments :
!        --------------------

!        KIDIA     : START OF HORIZONTAL LOOP
!        KFDIA     : END   OF HORIZONTAL LOOP
!        KLON      : HORIZONTAL DIMENSION
!        KLEV      : END OF VERTICAL LOOP AND VERTICAL DIMENSION
!        PTSPHY    : time-step for physics
!        PQP       : model q profiles
!        PTENT     : model tendency of T to create updated state for convection
!        PTENQ     : model tendency of q to create updated state for convection
!        PAPRSF    : PRESSURE ON FULL LEVELS
!        PCUCONVCA : CA state
!        PNLCONVCA : indicates if CA was alive at a gridpoint before

!        Implicit arguments :  none.
!        --------------------

!     Method.
!     -------
!        
!     Externals.
!     ----------
!      none

!     Reference.
!     ----------
!        none

!     Author.
!     -------
!        Martin Steinheimer  *ECMWF*

!     Modifications.
!     --------------
!        Original : 30-04-2010
!        F. Vana  16-May-2012 rewritten to tendency fluctuation + optimization
!     ------------------------------------------------------------------

  USE PARKIND1  ,ONLY : JPIM     ,JPRB
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK 
  USE YOMCST   , ONLY : RPI
  USE YOE_CUCONVCA, ONLY : TECUCONVCA
  USE YOMCT3   , ONLY : NSTEP
  IMPLICIT NONE

!     ------------------------------------------------------------------

 TYPE(TECUCONVCA)   ,INTENT(INOUT) :: YDECUCONVCA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
  INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
  INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV
  REAL(KIND=JPRB),INTENT(IN)       :: PTSPHY
  REAL(KIND=JPRB),INTENT(IN)       :: PQP(KLON,KLEV)
  REAL(KIND=JPRB),INTENT(INOUT)    :: PTENT(KLON,KLEV)
  REAL(KIND=JPRB),INTENT(INOUT)    :: PTENQ(KLON,KLEV)
  REAL(KIND=JPRB),INTENT(IN)       :: PAPRSF(KLON,KLEV)
  REAL(KIND=JPRB),INTENT(IN)       :: PCUCONVCA(KLON)
  REAL(KIND=JPRB),INTENT(IN)       :: PNLCONVCA(KLON)

!     ------------------------------------------------------------------

  INTEGER(KIND=JPIM)               :: JK
  REAL(KIND=JPHOOK)                  :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('CA_PROFPERT',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

  IF ((NSTEP > 1) .AND. (YDECUCONVCA%CA_FORC == 'CONV_TQ' )) THEN
    !compute T and Q tendency update for convection
    DO JK=1,KLEV
      WHERE ((PNLCONVCA(KIDIA:KFDIA) == 1).AND.(PCUCONVCA(KIDIA:KFDIA) > 2).AND.&
         &   (PAPRSF(KIDIA:KFDIA,JK) > 20000._JPRB))
        PTENT(KIDIA:KFDIA,JK) = PTENT(KIDIA:KFDIA,JK) -0.2_JPRB*SIN(2*RPI/(PAPRSF(KIDIA:KFDIA,KLEV)-20000._JPRB)*&
         & (PAPRSF(KIDIA:KFDIA,JK)-20000._JPRB)) /  PTSPHY
        PTENQ(KIDIA:KFDIA,JK) = PTENQ(KIDIA:KFDIA,JK) + PQP(KIDIA:KFDIA,JK)*0.02_JPRB / PTSPHY
      ENDWHERE
      WHERE ((PNLCONVCA(KIDIA:KFDIA) == 1).AND.(PCUCONVCA(KIDIA:KFDIA) <= 2).AND.&
         &   (PAPRSF(KIDIA:KFDIA,JK) > 20000._JPRB))
        PTENT(KIDIA:KFDIA,JK) = PTENT(KIDIA:KFDIA,JK)+0.2_JPRB*SIN(2*RPI/(PAPRSF(KIDIA:KFDIA,KLEV)-20000._JPRB)*&
         & (PAPRSF(KIDIA:KFDIA,JK)-20000._JPRB)) / PTSPHY
        PTENQ(KIDIA:KFDIA,JK) = PTENQ(KIDIA:KFDIA,JK) - PQP(KIDIA:KFDIA,JK)*0.02_JPRB / PTSPHY
      ENDWHERE
    ENDDO
  ENDIF  
!     ------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('CA_PROFPERT',1,ZHOOK_HANDLE)
!     ------------------------------------------------------------------
END SUBROUTINE CA_PROFPERT
