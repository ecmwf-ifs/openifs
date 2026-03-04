! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWTT1 ( YDESWRT,KIDIA,KFDIA,KLON,KNU,KABS,KKIND, PU, PTR )

!**** *SWTT1* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWTT1* IS CALLED FROM *SW1S*.

!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! KKIND   : (KABS)              ; INDICES OF THE ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
!     AND HORNER'S ALGORITHM.

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 95-01-20
!        03-10-10 Deborah Salmond and Marta Janiskova Optimisation
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
   
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESW    , ONLY : TESWRT

IMPLICIT NONE

TYPE(TESWRT)      ,INTENT(IN)    :: YDESWRT
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KABS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
INTEGER(KIND=JPIM),INTENT(IN)    :: KKIND(KABS) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KABS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR(KLON,KABS) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZR1(KLON), ZR2(KLON), ZU(KLON)
REAL(KIND=JPRB) :: ZRR

INTEGER(KIND=JPIM) :: IA, JA, JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

!*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION

IF (LHOOK) CALL DR_HOOK('SWTT1',0,ZHOOK_HANDLE)
DO JA = 1,KABS
  IA=KKIND(JA)
  DO JL = KIDIA,KFDIA
    ZU(JL) = PU(JL,JA)
    ZR1(JL) = YDESWRT%APAD(KNU,IA,1) + ZU(JL) * (YDESWRT%APAD(KNU,IA,2) + ZU(JL)&
     & * ( YDESWRT%APAD(KNU,IA,3) + ZU(JL) * (YDESWRT%APAD(KNU,IA,4) + ZU(JL)&
     & * ( YDESWRT%APAD(KNU,IA,5) + ZU(JL) * (YDESWRT%APAD(KNU,IA,6) + ZU(JL)&
     & * ( YDESWRT%APAD(KNU,IA,7) ))))))  

    ZR2(JL) = YDESWRT%BPAD(KNU,IA,1) + ZU(JL) * (YDESWRT%BPAD(KNU,IA,2) + ZU(JL)&
     & * ( YDESWRT%BPAD(KNU,IA,3) + ZU(JL) * (YDESWRT%BPAD(KNU,IA,4) + ZU(JL)&
     & * ( YDESWRT%BPAD(KNU,IA,5) + ZU(JL) * (YDESWRT%BPAD(KNU,IA,6) + ZU(JL)&
     & * ( YDESWRT%BPAD(KNU,IA,7) ))))))  
    ZRR=1.0_JPRB/ZR2(JL)

!*         2.      ADD THE BACKGROUND TRANSMISSION

    PTR(JL,JA) = (ZR1(JL)*ZRR) * (1.0_JPRB-YDESWRT%D(KNU,IA)) + YDESWRT%D(KNU,IA)
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SWTT1',1,ZHOOK_HANDLE)
END SUBROUTINE SWTT1
