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

SUBROUTINE SUVFE_BASIS(KORDER,PVAL,KBASIS,PKNOT,PSPLINE)

!**** *SUVFE_BASIS*  - define VFE vertical operator (B-splines): basis

!**   Interface.
!     ----------

!     *CALL* SUVFE_BASIS

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        KORDER    - order of spline
!                    ( KORDER =  polynomial degree + 1, for example linear functions KORDER=2 )
!        PVAL      - parameter value for which splines are evaluated
!        KBASIS    - number of basis functions we wants to compute:
!                    * this must be consistent with the choice of knot vector
!                    * be carefull, KBASIS may be not equal to NFLEVG, because
!                      every additional boundary condition increases KBASIS+1
!                      and sometimes we need to compute even more functions
!                      in order to compute boundary elements that satisfy
!                      specific boundary condition.
!        PKNOT     - knots sequence (first and last knots must have multiplicity k=KORDER).
!                    PKNOT(i+1) must be >= PKNOT(i) for every i;
!                    multiplicity of nodes k is knot bigger that KORDER => max(k) <= KORDER

!      * OUTPUT:
!        PSPLINE   - array of computed values of basis functions evaluated in PVAL.

!     Method.
!     -------
!        Method: Recursive algorithm of DeBoor

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN-NH documentation.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE 
!        Original : 2010-09

!     Modifications.
!     --------------
!   J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!-----------------------------------------------------------------------------------------

USE PARKIND1,ONLY : JPRB, JPIM
USE YOMHOOK ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!-----------------------------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KORDER
INTEGER(KIND=JPIM),INTENT(IN) :: KBASIS
REAL(KIND=JPRB),INTENT(IN)    :: PVAL
REAL(KIND=JPRB),INTENT(IN)    :: PKNOT(KBASIS+KORDER)
REAL(KIND=JPRB),INTENT(OUT)   :: PSPLINE(KBASIS)

!-----------------------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IKNOT,JI,JK
REAL(KIND=JPRB)    :: ZTEMP(KBASIS+KORDER),ZD,ZE
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!-----------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_BASIS',0,ZHOOK_HANDLE)
!-----------------------------------------------------------------------------------------

IKNOT = KBASIS + KORDER

! calculate the first order basis functions n(ji)(1)
DO JI=1,IKNOT-1
  IF ((PVAL >= PKNOT(JI)) .AND. (PVAL < PKNOT(JI+1)))THEN
    ZTEMP(JI) = 1.0_JPRB
  ELSE
    ZTEMP(JI) = 0.0_JPRB
  ENDIF
ENDDO

! calculate the higher order basis functions 
DO JK=2,KORDER
  DO JI=1,IKNOT-JK
    IF( ZTEMP(JI) /= 0.0_JPRB )THEN    ! if the lower order basis function is zero skip the calculation 
      ZD = ((PVAL-PKNOT(JI))*ZTEMP(JI))/(PKNOT(JI+JK-1)-PKNOT(JI))
    ELSE
      ZD = 0.0_JPRB
    ENDIF

    IF( ZTEMP(JI+1) /= 0.0_JPRB )THEN     ! if the lower order basis function is zero skip the calculation 
      ZE = ((PKNOT(JI+JK)-PVAL)*ZTEMP(JI+1))/(PKNOT(JI+JK)-PKNOT(JI+1))
    ELSE
      ZE = 0.0_JPRB
    ENDIF

    ZTEMP(JI) = ZD + ZE
  ENDDO
ENDDO

! last point of interval
IF( ABS(PVAL - 1.0_JPRB) < TINY(1.0_JPRB) )THEN
   ZTEMP         = 0.0_JPRB
   ZTEMP(KBASIS) = 1.0_JPRB
ENDIF

! put in an array
DO JI=1,KBASIS
  PSPLINE(JI) = ZTEMP(JI)
ENDDO

!-----------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_BASIS',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_BASIS
