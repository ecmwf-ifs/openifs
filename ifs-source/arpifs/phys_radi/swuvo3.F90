! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE SWUVO3 &
 & ( KIDIA,KFDIA,KLON,KNU,KABS,&
 & KEXPO3, PEXPO3,&
 & PU, PTR &
 & )  
  
!**** *SWUVO3* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

!     PURPOSE.
!     --------
!           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR OZONE
!     IN THE UV and VISIBLE SPECTRAL INTERVALS.

!**   INTERFACE.
!     ----------
!          *SWUVO3* IS CALLED FROM *SW1S*.

!        EXPLICIT ARGUMENTS :
!        --------------------
! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
! KABS   :                     ; NUMBER OF ABSORBERS
! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
!     ==== OUTPUTS ===
! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          TRANSMISSION FUNCTION ARE COMPUTED USING SUMS OF EXPONENTIALS

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
!        ORIGINAL : 00-12-18
!        Modified J. HAGUE          03-01-03 MASS Vector Functions       
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
   
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMJFH   , ONLY : N_VMASS

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KABS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KNU 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEXPO3(6)
REAL(KIND=JPRB),   INTENT(IN)    :: PEXPO3(6,2,7)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PU(KLON,KABS) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTR(KLON,KABS) 
!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZU(KLON)
REAL(KIND=JPRB) :: ZTMP1(KFDIA-KIDIA+1)
REAL(KIND=JPRB) :: ZTMP2(KFDIA-KIDIA+1)

INTEGER(KIND=JPIM) ::  JA, JL, IEXP, JX, JLEN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SWUVO3',0,ZHOOK_HANDLE)
IEXP=KEXPO3(KNU)

IF(N_VMASS > 0) THEN
  JLEN=KFDIA-KIDIA+1
ENDIF

DO JA = 1,KABS
  DO JL=KIDIA,KFDIA
    PTR(JL,JA)=0.0_JPRB
  ENDDO
  
  IF(N_VMASS <= 0) THEN ! Do not use Vector Mass

    DO JX=1,IEXP
      DO JL = KIDIA,KFDIA
        ZU(JL) = PU(JL,JA)
        PTR(JL,JA) = PTR(JL,JA)+PEXPO3(KNU,1,JX)*EXP(-PEXPO3(KNU,2,JX)*ZU(JL))
      ENDDO
    ENDDO

  ELSE  ! Use Vector MASS

    DO JX=1,IEXP
      DO JL = KIDIA,KFDIA
        ZTMP1(JL-KIDIA+1)=-PEXPO3(KNU,2,JX)*PU(JL,JA)
      ENDDO
  
      CALL VEXP(ZTMP2,ZTMP1,JLEN)
  
      DO JL = KIDIA,KFDIA
        PTR(JL,JA) = PTR(JL,JA)+PEXPO3(KNU,1,JX)*ZTMP2(JL-KIDIA+1)
      ENDDO
    ENDDO    

  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('SWUVO3',1,ZHOOK_HANDLE)
END SUBROUTINE SWUVO3
