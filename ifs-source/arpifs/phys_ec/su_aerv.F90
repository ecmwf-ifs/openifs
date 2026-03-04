! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SU_AERV &
 & (  YDERAD,KLEV  , PETAH )  

!**** *SU_AERV* - PARAMETERS FOR THE VERTICAL DISTRIBUTIONS OF AEROSOLS.

!     PURPOSE.
!     --------

!          THIS ROUTINE COMPUTES THE VALUES *PVDAEN* (*N=*S,*L,*U OR *D
!     FOR SEA,LAND,URBAN OR DESERT) OF A SURFACE-NORMALISED VERTICAL
!     DISTRIBUTION OF AEROSOLS' OPTICAL DEPHTS FROM THE ARGUMENT *PETAH*
!     (VERTICAL COORDINATE) AT *KLEVP1* LEVELS. 


!**   INTERFACE.
!     ----------

!          *SU_AERV* IS CALLED FROM *SUECRAD*.
!          THERE ARE SIXTEEN DUMMY ARGUMENTS: *PETAH* IS THE VERTICAL
!     COORDINATE.
!                                            *PVDAEN* (*N=*SS,*DU,*OM, 
!     *BC OR *SU) ARE THE NORMALISED VERTICAL DISTRIBUTIONS.
!                                             *KLEVP1* IS THE NUMBER OF
!     LEVELS.

!     METHOD.
!     -------

!          STRAIGHTFORWARD, EQUIVALENT HEIGTHS ARE GIVEN IN METERS (8434
!     FOR THE ATMOSPHERE) AND TROPOSPHERIC AND STRATOSPHERIC PRESSURE
!     BOUNDARY VALUES ARE SET AT 101325 AND 19330 *PASCAL.

!     EXTERNALS.
!     ----------

!          NONE.

!     REFERENCE.
!     ----------

!          NONE.

!     AUTHOR
!     ------
!     JJMorcrette ECMWF, 20120815, heavily borrowed from JFGeleyn's *SUAERV*

!     MODIFICATIONS
!     -------------
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOERAD    ,ONLY : TERAD

IMPLICIT NONE

TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PETAH(KLEV+1) 

!     ------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JK

REAL(KIND=JPRB) :: ZHSBC, ZHSDU, ZHSOM, ZHSSS, ZHSSU
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*         1.     COMPUTATIONS.
!                 -------------

IF (LHOOK) CALL DR_HOOK('SU_AERV',0,ZHOOK_HANDLE)



ZHSSS=MAX(1.0_JPRB, 8434._JPRB / YDERAD%RAESHSS )
ZHSDU=MAX(1.0_JPRB, 8434._JPRB / YDERAD%RAESHDU )
ZHSOM=MAX(1.0_JPRB, 8434._JPRB / YDERAD%RAESHOM )
ZHSBC=MAX(1.0_JPRB, 8434._JPRB / YDERAD%RAESHBC )
ZHSSU=MAX(1.0_JPRB, 8434._JPRB / YDERAD%RAESHSU )

YDERAD%CVDAESS(1)=0._JPRB
YDERAD%CVDAEDU(1)=0._JPRB
YDERAD%CVDAEOM(1)=0._JPRB
YDERAD%CVDAEBC(1)=0._JPRB
YDERAD%CVDAESU(1)=0._JPRB
IF(PETAH(1) /= 0.0_JPRB) THEN
  YDERAD%CVDAESS(1)=PETAH(1)**ZHSSS
  YDERAD%CVDAEDU(1)=PETAH(1)**ZHSDU
  YDERAD%CVDAEOM(1)=PETAH(1)**ZHSOM
  YDERAD%CVDAEBC(1)=PETAH(1)**ZHSBC
  YDERAD%CVDAESU(1)=PETAH(1)**ZHSSU
ENDIF
DO JK=2,KLEV+1
  YDERAD%CVDAESS(JK)=PETAH(JK)**ZHSSS
  YDERAD%CVDAEDU(JK)=PETAH(JK)**ZHSDU
  YDERAD%CVDAEOM(JK)=PETAH(JK)**ZHSOM
  YDERAD%CVDAEBC(JK)=PETAH(JK)**ZHSBC
  YDERAD%CVDAESU(JK)=PETAH(JK)**ZHSSU
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SU_AERV',1,ZHOOK_HANDLE)
END SUBROUTINE SU_AERV
