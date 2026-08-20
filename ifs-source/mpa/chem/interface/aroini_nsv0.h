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

INTERFACE
SUBROUTINE AROINI_NSV0(OUSECHEM,OORILAM,ODUST,OCO2, ODEPOS, KGFL,KGFL_EXT, &
                      KSV_CHEMBEG, KSV_CHEMEND, KSV_DSTBEG, KSV_DSTEND,    &
                      KSV_DSTDEPBEG, KSV_DSTDEPEND, KSV_AERBEG, KSV_AEREND,&
                      KSV_CO2, CLNAME)

USE PARKIND1  ,ONLY : JPIM

!**** *SUINI_NSV * - Setting up the boundaries of the scalar vector GFL_EXT

!     Purpose.
!     --------
!           Initialization of YOMNSV variables
!**   Interface.
!     ----------
!        *CALL* *AROINI_NSV0*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------
!           Common MODD_CH_M9
!           Common MODD_CH_AEROSOL

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        

!     Author.
!     -------
!        P. Tulet  *CNRM*

!     Modifications.
!     --------------

IMPLICIT NONE

LOGICAL,INTENT(IN):: OORILAM
LOGICAL,INTENT(IN):: ODUST
LOGICAL,INTENT(IN):: OCO2
LOGICAL,INTENT(IN):: OUSECHEM
LOGICAL,INTENT(IN):: ODEPOS
INTEGER(KIND=JPIM),INTENT(IN) ::KGFL 
INTEGER(KIND=JPIM),INTENT(INOUT) ::KGFL_EXT 
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_CHEMBEG
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_CHEMEND
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_DSTBEG
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_DSTEND
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_DSTDEPBEG
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_DSTDEPEND
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_AERBEG
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_AEREND
INTEGER(KIND=JPIM),INTENT(INOUT) ::KSV_CO2
CHARACTER(LEN=15), DIMENSION(KGFL),  INTENT(INOUT) :: CLNAME

END SUBROUTINE AROINI_NSV0
END INTERFACE
