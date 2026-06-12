MODULE KPP_TRIDCOF_MOD
CONTAINS
SUBROUTINE KPP_TRIDCOF &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,LDKPPCAL ,&
  &   PDIFF    ,PTRI0    ,PTRI1    ,PCU      ,PCC      ,&
  &   PCL      )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine sets coefficients for tridiagonal matrix. 

! Interface :
! ---------
!   Call *KPP_TRIDCOF* from *KPP_OCNINT*

! Method :
! ------
!   implicit integration scheme (Backward Euler Method) 

! Externals :
! ---------

! Reference :
! ---------
! Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.

! Modifications :
! -------------
!     06-Jun-1994  Bill Large                 
!            2002  Steve Woolnough, Reading Univ. 
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
REAL(KIND=JPRB),INTENT(IN)  :: PDIFF(KLON,0:KLEVO)    ! Diffusivity profile 
REAL(KIND=JPRB),INTENT(IN)  :: PTRI0(KLON,0:KLEVO)    ! Array for diff. eq.
REAL(KIND=JPRB),INTENT(IN)  :: PTRI1(KLON,0:KLEVO)    ! Array for diff. eq.
REAL(KIND=JPRB),INTENT(OUT) :: PCU(KIDIA:KFDIA,KLEVO) ! Upper coef. for (K-1)  
                                                      ! on k line of trid. mtx
REAL(KIND=JPRB),INTENT(OUT) :: PCC(KIDIA:KFDIA,KLEVO) ! central ...   (K  ) ..
REAL(KIND=JPRB),INTENT(OUT) :: PCL(KIDIA:KFDIA,KLEVO) ! lower .....   (K-1) ..

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)

INTEGER(KIND=JPIM) :: JZ                    ! Loop control variable
INTEGER(KIND=JPIM) :: JL                    ! Loop control variable

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_TRIDCOF_MOD:KPP_TRIDCOF',0,ZHOOK_HANDLE)

DO JL = KIDIA, KFDIA
  IF( LDKPPCAL(JL) ) THEN

!   coef. in the surface layer
!   note: cu(1) = 0. and cl(jpnz) = 0. are necessary conditions.
    PCU(JL,1) = 0.0_JPRB
    PCL(JL,1) =          - PTRI1(JL,1) * PDIFF(JL,1) 
              !   - dto/h(1)/dzb(1)*diff(1)
    PCC(JL,1) = 1.0_JPRB - PCL(JL,1)
              ! 1.+ dto/h(1)/dzb(1)*diff(1)

!   coef. in the bottom layer
    PCL(JL,KLEVO)= 0.0_JPRB

!   coef. inside the domain
    DO JZ = 2, KLEVO
      PCU(JL,JZ) = - PTRI0(JL,JZ) * PDIFF(JL,JZ-1)  
                 ! - dto/h(jz)/dzb(jz-1)*diff(jz-1)
      PCL(JL,JZ) = - PTRI1(JL,JZ) * PDIFF(JL,JZ)    
                 ! - dto/h(jz)/dzb(jz)*diff(jz)
      PCC(JL,JZ) = 1.0_JPRB - PCL(JL,JZ) - PCU(JL,JZ)
                 ! 1.+ dto/h(jz)/dzb(jz)*diff(jz)
    ENDDO   

  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('KPP_TRIDCOF_MOD:KPP_TRIDCOF',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_TRIDCOF
END MODULE KPP_TRIDCOF_MOD
