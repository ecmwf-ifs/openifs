SUBROUTINE SURDI
! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURDI*   - INITIALIZE COMMON YOERDI CONTROLLING RADINT

!     PURPOSE.
!     --------
!           INITIALIZE YOERDI, THE COMMON THAT CONTROLS THE
!           RADIATION INTERFACE

!**   INTERFACE.
!     ----------
!        CALL *SURDI* FROM *SURAD*
!              -----        -----

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
!        Original  JEAN-JACQUES MORCRETTE  *ECMWF*
!        Modified   P. Viterbo   99-03-26    Tiling of the land surface
!        Modified   P. Viterbo   24-05-2004  surf library

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        JJMorcrette   2004-10-07 Gas concentrations
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERDI   , ONLY : RRAE     ,&
 & RCARDI   ,RCH4     ,RN2O     ,RO3      ,RCFC11   ,&
 & RCFC12   ,REPCLC   ,REPH2O   ,RSUNDUR  ,&
 & RCCO2    ,RCCH4    ,RCN2O    ,RCCFC11  ,RCCFC12

IMPLICIT NONE

REAL(KIND=JPRB) :: ZAIRMWG, ZC11MWG, ZC12MWG, ZCH4MWG, ZCO2MWG, ZN2OMWG, ZO3MWG
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      ----------------------------------------------------------------

!*       1.    SET DEFAULT VALUES.
!              -------------------

IF (LHOOK) CALL DR_HOOK('SURDI',0,ZHOOK_HANDLE)
RRAE = 0.1277E-02_JPRB

!* Threshold for computing sunshine duration (W/m2)
RSUNDUR=120._JPRB

!*  For sea ice, monthly values are based on Ebert and Curry, 1993, Table 2.
!   We take dry snow albedo as the representative value for non-summer
!   months, and bare sea-ice as the representative value for summer
!   months. The values for Antarctic are shifted six-months.
! All computations brought back to *SUSWN*

!*  Concentration of the various trace gases (IPCC/SACC values for 1990)
!        CO2         CH4        N2O        CFC11       CFC12
!      353ppmv     1.72ppmv   310ppbv     280pptv     484pptv

ZAIRMWG = 28.970_JPRB
ZCO2MWG = 44.011_JPRB
ZCH4MWG = 16.043_JPRB
ZN2OMWG = 44.013_JPRB
ZO3MWG  = 47.9982_JPRB
ZC11MWG = 137.3686_JPRB
ZC12MWG = 120.9140_JPRB

RCCO2   = 353.E-06_JPRB
RCCH4   = 1.72E-06_JPRB
RCN2O   = 310.E-09_JPRB
RCCFC11 = 280.E-12_JPRB
RCCFC12 = 484.E-12_JPRB
!RCARDI  = 353.E-06_JPRB*ZCO2MWG/ZAIRMWG
!RCH4    = 1.72E-06_JPRB*ZCH4MWG/ZAIRMWG
!RN2O    = 310.E-09_JPRB*ZN2OMWG/ZAIRMWG
!RO3     =   1.E-06_JPRB*ZO3MWG /ZAIRMWG
!RCFC11  = 280.E-12_JPRB*ZC11MWG/ZAIRMWG
!RCFC12  = 484.E-12_JPRB*ZC12MWG/ZAIRMWG

RCARDI  = RCCO2   * ZCO2MWG/ZAIRMWG
RCH4    = RCCH4   * ZCH4MWG/ZAIRMWG
RN2O    = RCN2O   * ZN2OMWG/ZAIRMWG
RO3     = 1.E-06_JPRB*ZO3MWG /ZAIRMWG
RCFC11  = RCCFC11 * ZC11MWG/ZAIRMWG
RCFC12  = RCCFC12 * ZC12MWG/ZAIRMWG

REPCLC=1.E-12_JPRB
REPH2O=1.E-12_JPRB

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURDI',1,ZHOOK_HANDLE)
END SUBROUTINE SURDI
