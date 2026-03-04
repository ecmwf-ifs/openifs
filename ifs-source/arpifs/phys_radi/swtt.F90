! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWTT ( YDESWRT, KIDIA, KFDIA, KLON, KNU, KA , PU, PTR)

!**** *SWTT* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
!     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
!     INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWTT* IS CALLED FROM *SW1S*, *SWNI*.

!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KA     :                     ; INDEX OF THE ABSORBER
! PU     : (KLON)             ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON)             ; TRANSMISSION FUNCTION

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
!        ORIGINAL : 88-12-15
!        M.Hamrud     01-Oct-2003 CY28 Cleaning
!        M.Janiskova  18-Mar-2008 modifications for optimizied TL/AD 
   
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOESW    , ONLY : TESWRT

IMPLICIT NONE

TYPE(TESWRT)      ,INTENT(IN)    :: YDESWRT
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
INTEGER(KIND=JPIM),INTENT(IN)    :: KA 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR(KLON) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZR1(KLON), ZR2(KLON)
REAL(KIND=JPRB) :: ZDIV

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------

!*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION

IF (LHOOK) CALL DR_HOOK('SWTT',0,ZHOOK_HANDLE)
DO JL = KIDIA,KFDIA
  ZR1(JL) = YDESWRT%APAD(KNU,KA,1) + PU(JL) * (YDESWRT%APAD(KNU,KA,2) + PU(JL)&
   & * ( YDESWRT%APAD(KNU,KA,3) + PU(JL) * (YDESWRT%APAD(KNU,KA,4) + PU(JL)&
   & * ( YDESWRT%APAD(KNU,KA,5) + PU(JL) * (YDESWRT%APAD(KNU,KA,6) + PU(JL)&
   & * ( YDESWRT%APAD(KNU,KA,7) ))))))  

  ZR2(JL) = YDESWRT%BPAD(KNU,KA,1) + PU(JL) * (YDESWRT%BPAD(KNU,KA,2) + PU(JL)&
   & * ( YDESWRT%BPAD(KNU,KA,3) + PU(JL) * (YDESWRT%BPAD(KNU,KA,4) + PU(JL)&
   & * ( YDESWRT%BPAD(KNU,KA,5) + PU(JL) * (YDESWRT%BPAD(KNU,KA,6) + PU(JL)&
   & * ( YDESWRT%BPAD(KNU,KA,7) ))))))  

!*         2.      ADD THE BACKGROUND TRANSMISSION

  ZDIV = 1.0_JPRB/ZR2(JL)
  PTR(JL) = (ZR1(JL) * ZDIV) * (1.0_JPRB - YDESWRT%D(KNU,KA)) + YDESWRT%D(KNU,KA)
ENDDO

IF (LHOOK) CALL DR_HOOK('SWTT',1,ZHOOK_HANDLE)
END SUBROUTINE SWTT
