! (C) Copyright 2003- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

SUBROUTINE SURDI(YDERDI)

!**** *SURDI*   - INITIALIZE COMMON YOERDI CONTROLLING RADINT

!     PURPOSE.
!     --------
!           INITIALIZE YOERDI, THE COMMON THAT CONTROLS THE
!           RADIATION INTERFACE

!**   INTERFACE.
!     ----------
!        CALL *SURDI* FROM *SURAD* and *SUECRAD*
!              -----        -----       -------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOERDI

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS MODEL

!     AUTHOR.
!     -------
!      Original  JEAN-JACQUES MORCRETTE  *ECMWF*
!      ORIGINAL : 88-12-15

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      Modified   P. Viterbo   24-05-2004  surf library
!      JJMorcrette   2004-10-07 Gas concentrations
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOERDI   , ONLY : TERDI
! Molar masses
USE YOMCST    ,ONLY : RMD, RMCO2, RMCH4, RMN2O, RMNO2, RMCFC11, RMCFC12, RMO3, RMHCFC22, RMCCL4

IMPLICIT NONE

TYPE(TERDI),INTENT(INOUT):: YDERDI
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

IF (LHOOK) CALL DR_HOOK('SURDI',0,ZHOOK_HANDLE)
ASSOCIATE(RCARDI=>YDERDI%RCARDI, RCCCL4=>YDERDI%RCCCL4, RCCFC11=>YDERDI%RCCFC11, &
 & RCCFC12=>YDERDI%RCCFC12, RCCFC22=>YDERDI%RCCFC22, RCCH4=>YDERDI%RCCH4, &
 & RCCL4=>YDERDI%RCCL4, RCCO2=>YDERDI%RCCO2, RCFC11=>YDERDI%RCFC11, &
 & RCFC12=>YDERDI%RCFC12, RCFC22=>YDERDI%RCFC22, RCH4=>YDERDI%RCH4, &
 & RCN2O=>YDERDI%RCN2O, RCNO2=>YDERDI%RCNO2, REPCLC=>YDERDI%REPCLC, &
 & REPH2O=>YDERDI%REPH2O, RN2O=>YDERDI%RN2O, RNO2=>YDERDI%RNO2, RO3=>YDERDI%RO3, &
 & RRAE=>YDERDI%RRAE, RSUNDUR=>YDERDI%RSUNDUR)
RRAE = 0.1277E-02_JPRB

!* Threshold for computing sunshine duration (W/m2)
RSUNDUR=120._JPRB

! Compute mass mixing ratios from reference volume mixing ratios
RCARDI  = RCCO2   * RMCO2   / RMD
RCH4    = RCCH4   * RMCH4   / RMD
RN2O    = RCN2O   * RMN2O   / RMD
RNO2    = RCNO2   * RMNO2   / RMD
RO3     = 1.E-06_JPRB*RMO3  / RMD
RCFC11  = RCCFC11 * RMCFC11 / RMD
RCFC12  = RCCFC12 * RMCFC12 / RMD
RCFC22  = RCCFC22 * RMHCFC22/ RMD
RCCL4   = RCCCL4  * RMCCL4  / RMD

REPCLC=1.E-12_JPRB
REPH2O=1.E-12_JPRB

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURDI',1,ZHOOK_HANDLE)
END SUBROUTINE SURDI
