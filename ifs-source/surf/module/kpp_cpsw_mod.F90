MODULE KPP_CPSW_MOD
CONTAINS
SUBROUTINE KPP_CPSW &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEV     ,LDKPPCAL ,&
  &   PS       ,PT1      ,P0       ,PCPSW    )

! (C) Copyright 2008- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
! Compute Specific heat

! Interface :
! ---------
!       PRESSURE        P0       dbar
!       TEMPERATURE     PT1      deg.C (IPTS-68)
!       SALINITY        PS       (IPSS-78)
!       SPECIFIC HEAT   PCPSW    J/(KG DEG C)

! Method :
! ------

! Externals :
! ---------

! Reference :
! ---------
!   Millero et al., 1973, JGR, 78, 4499-4507
!   Millero et al., 1981, UNESCO Report No. 38 pp. 99-188.

! Modifications :
! -------------
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEV
REAL(KIND=JPRB),INTENT(IN)  :: PT1(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: PS(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB),INTENT(IN)  :: P0(KIDIA:KFDIA,KLEV)
REAL(KIND=JPRB),INTENT(OUT) :: PCPSW(KLON,KLEV)

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)

INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: JZ

REAL(KIND=JPRB) :: ZT
REAL(KIND=JPRB) :: ZP
REAL(KIND=JPRB) :: ZCP0
REAL(KIND=JPRB) :: ZCP1
REAL(KIND=JPRB) :: ZCP2
REAL(KIND=JPRB) :: ZSR
REAL(KIND=JPRB) :: ZA
REAL(KIND=JPRB) :: ZB
REAL(KIND=JPRB) :: ZC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_CPSW_MOD:KPP_CPSW',0,ZHOOK_HANDLE)

DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    DO JZ = 1, KLEV

      ZT = MAX( PT1(JL,JZ), -2.0_JPRB ) 

!     Scale pressure to bar
      ZP = P0(JL,JZ) * 0.1_JPRB

!     Sqrt salinity for fractional term
      ZSR = SQRT( ABS(PS(JL,JZ)) )

!     Specific heat ZCP0 for ZP=0 (Millero et al., UNESCO 1981)
      ZA = (-1.38385E-3_JPRB*ZT+0.1072763_JPRB)*ZT &
         &  -7.643575_JPRB
      ZB = (5.148E-5_JPRB*ZT-4.07718E-3_JPRB )*ZT &
         &  +0.1770383_JPRB
      ZC = (((2.093236E-5_JPRB*ZT-2.654387E-3_JPRB)*ZT &
         &   +0.1412855_JPRB)*ZT-3.720283_JPRB)*ZT    &
         &  + 4217.4_JPRB
      ZCP0 = (ZB*ZSR+ZA)*PS(JL,JZ)+ZC

!     ZCP1 pressure and temperature terms for PS = 0
      ZA = (((1.7168E-8_JPRB*ZT+2.0357E-6_JPRB)*ZT     &
         &    -3.13885E-4_JPRB)*ZT+1.45747E-2_JPRB)*ZT &
         &  - 0.49592_JPRB
      ZB = (((2.2956E-11_JPRB*ZT-4.0027E-9_JPRB)*ZT    &
         &    +2.87533E-7_JPRB)*ZT-1.08645E-5_JPRB)*ZT &
         &  + 2.4931E-4_JPRB
      ZC = ((6.136E-13_JPRB*ZT-6.5637E-11_JPRB)*ZT     &
         &  +2.6380E-9_JPRB )*ZT-5.422E-8_JPRB
      ZCP1 = ((ZC*ZP+ZB)*ZP+ZA)*ZP

!     ZCP2 pressure and temperature terms for PS > 0
      ZA = (((-2.9179E-10_JPRB*ZT+2.5941E-8_JPRB)*ZT   &
         &   +9.802E-7_JPRB)*ZT-1.28315E-4_JPRB)*ZT    &
         & +4.9247E-3
      ZB = (3.122E-8_JPRB*ZT-1.517E-6_JPRB)*ZT &
         & - 1.2331E-4_JPRB
      ZA = (ZA+ZB*ZSR)*PS(JL,JZ)
      ZB = ((1.8448E-11_JPRB*ZT-2.3905E-9_JPRB)*ZT &
         &  +1.17054E-7_JPRB)*ZT-2.9558E-6_JPRB
      ZB = (ZB+9.971E-8_JPRB*ZSR)* PS(JL,JZ)
      ZC = ( 3.513E-13_JPRB*ZT-1.7682E-11_JPRB)*ZT &
         & + 5.540E-10_JPRB
      ZC = (ZC-1.4300E-12_JPRB*ZT*ZSR)*PS(JL,JZ)
      ZCP2 = ((ZC*ZP+ZB)*ZP+ZA)*ZP

!     Specific heat 
      PCPSW(JL,JZ) = ZCP0 + ZCP1 + ZCP2

    ENDDO
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('KPP_CPSW_MOD:KPP_CPSW',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_CPSW
END MODULE KPP_CPSW_MOD
