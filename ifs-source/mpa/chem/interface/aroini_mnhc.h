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
SUBROUTINE AROINI_MNHC(OUSECHEM, OORILAM, ODUST, ODEPOS, &
                       OINITCHEM, OINITDUST, OINITORILAM,&
                       KDAY,KMONTH, KYEAR, KLUOUT, KPROC)

!!
!!*** *AROINI_MNHC*
!!
!!    PURPOSE
!!    -------
!        initialize chemical variables
!        initialize chemical core system
!!
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  *CNRM / GMEI* and contributors of MesoNH-C (K. Shure, C. Mari, V. Crassier)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/04/05
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE PARKIND1  ,ONLY : JPIM

!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL, INTENT(IN)  :: OUSECHEM
LOGICAL, INTENT(IN)  :: ODUST
LOGICAL, INTENT(IN)  :: OORILAM
LOGICAL, INTENT(IN)  :: ODEPOS
LOGICAL, INTENT(IN)  :: OINITCHEM
LOGICAL, INTENT(IN)  :: OINITDUST
LOGICAL, INTENT(IN)  :: OINITORILAM
INTEGER(KIND=JPIM), INTENT(IN)  :: KLUOUT   ! output listing channel
INTEGER(KIND=JPIM), INTENT(IN)  :: KYEAR    ! year of initialization
INTEGER(KIND=JPIM), INTENT(IN)  :: KMONTH   ! month of initialization
INTEGER(KIND=JPIM), INTENT(IN)  :: KDAY     ! day of initialization
INTEGER(KIND=JPIM), INTENT(IN)  :: KPROC    ! proc number
!
END SUBROUTINE AROINI_MNHC
END INTERFACE
