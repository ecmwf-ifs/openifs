MODULE KPP_SWFRAC_MOD
CONTAINS
SUBROUTINE KPP_SWFRAC &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KZ       ,&
  &   LDINIT_KPP,LDKPPCAL,PFACT    ,PDO      ,PSWDK_SAVE,&
  &   PSWDK    ,YDOCEAN_ML)

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine computes fraction of solar short-wave flux penetrating 
!   to specified depth.

! Interface :
! ---------
!   Call *TRIDCOF* from *OCNINT*
!   KZ : >= 0 use PSWDK_SAVE to optimize expoential operation
!         < 0 calculate directly

! Method :
! ------
!     two band solar absorption model of Simpson and  Paulson (1977)

! Externals :
! ---------

! Reference :
! ---------
! Large, W. G., J. C. McWilliams, S. C. Doney (1994), Rev. Geophys.
! Simpson, J. J. and C. A. Paulson (1979), Quart. J. Roy. Meteorol. Soc.

! Modifications :
! -------------
!     06-Jun-1994  Bill Large
!            2002  Steve Woolnough, Reading Univ. Optimized with lookup table 
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1,     ONLY : JPIM, JPRB
USE YOMHOOK,      ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_OCEAN_ML, ONLY : TOCEAN_ML

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KZ            ! index of vertical grid point

REAL(KIND=JPRB),INTENT(IN)    :: PFACT(KLON)   ! scale factor applied to depth 
REAL(KIND=JPRB),INTENT(IN)    :: PDO(KLON)     ! depth ( <0.) for desired SW 
REAL(KIND=JPRB),INTENT(INOUT) :: PSWDK_SAVE(KLON,0:KLEVO) !save array fro PSWDK
REAL(KIND=JPRB),INTENT(OUT)   :: PSWDK(KLON)   ! Solar rad fractional decay

LOGICAL,INTENT(IN) :: LDKPPCAL(KLON)
LOGICAL,INTENT(IN) :: LDINIT_KPP

TYPE(TOCEAN_ML),INTENT(IN) :: YDOCEAN_ML

INTEGER :: JL

REAL(KIND=JPRB) :: ZR1(KIDIA:KFDIA)            ! fraction of a band
REAL(KIND=JPRB) :: ZR2(KIDIA:KFDIA)            ! fraction of the other band
REAL(KIND=JPRB) :: ZA1_INV
REAL(KIND=JPRB) :: ZA2_INV

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_SWFRAC_MOD:KPP_SWFRAC',0,ZHOOK_HANDLE)
ASSOCIATE(NJERLOV_KPP=>YDOCEAN_ML%NJERLOV_KPP, &
 & RA1_JERLOV=>YDOCEAN_ML%RA1_JERLOV, RA2_JERLOV=>YDOCEAN_ML%RA2_JERLOV, &
 & RFAC_JERLOV=>YDOCEAN_ML%RFAC_JERLOV, RMIN_JERLOV=>YDOCEAN_ML%RMIN_JERLOV)

!   Jerlov water type :      I         Ia         Ib          II        III
!              jwtype        1          2          3           4          5
!
!    RFAC_JERLOV   / 0.58_JPRB, 0.62_JPRB, 0.67_JPRB,  0.77_JPRB, 0.78_JPRB /
!    RA1_JERLOV    / 0.35_JPRB, 0.6_JPRB , 1.0_JPRB ,  1.5_JPRB , 1.4_JPRB  /
!    RA2_JERLOV    / 23.0_JPRB, 20.0_JPRB, 17.0_JPRB, 14.0_JPRB , 7.9_JPRB  /
!    RMIN_JERLOV   / -80.0_JPRB  /

ZA1_INV = 1.0_JPRB / RA1_JERLOV(NJERLOV_KPP)
ZA2_INV = 1.0_JPRB / RA2_JERLOV(NJERLOV_KPP)

! 1. Set table for radiation
!    -----------------------

IF ( LDINIT_KPP .AND. KZ >= 0 ) THEN 

  DO JL = KIDIA, KFDIA
    IF ( LDKPPCAL(JL) ) THEN
!      ZR1(JL) = MAX( PDO(JL)*PFACT(JL)/RA1_JERLOV(NJERLOV_KPP) , RMIN_JERLOV )
!      ZR2(JL) = MAX( PDO(JL)*PFACT(JL)/RA2_JERLOV(NJERLOV_KPP) , RMIN_JERLOV )
      ZR1(JL) = MAX( PDO(JL)*PFACT(JL)*ZA1_INV , RMIN_JERLOV )
      ZR2(JL) = MAX( PDO(JL)*PFACT(JL)*ZA2_INV , RMIN_JERLOV )
      PSWDK_SAVE(JL,KZ) = RFAC_JERLOV(NJERLOV_KPP) * EXP(ZR1(JL)) &
                   & + ( 1.0_JPRB - RFAC_JERLOV(NJERLOV_KPP) ) * EXP(ZR2(JL))
    ENDIF 
  ENDDO

ENDIF

! 2. Use the table
!    -------------

IF ( KZ >= 0 ) THEN

  DO JL = KIDIA, KFDIA
    IF ( LDKPPCAL(JL) ) THEN
      PSWDK(JL)=PSWDK_SAVE(JL,KZ)
    ENDIF
  ENDDO

ELSE ! KZ < 0, calc. explicitly

  DO JL = KIDIA, KFDIA
    IF ( LDKPPCAL(JL) ) THEN
!      ZR1(JL)   = MAX( PDO(JL)*PFACT(JL)/RA1_JERLOV(NJERLOV_KPP), RMIN_JERLOV )
!      ZR2(JL)   = MAX( PDO(JL)*PFACT(JL)/RA2_JERLOV(NJERLOV_KPP), RMIN_JERLOV )
      ZR1(JL)   = MAX( PDO(JL)*PFACT(JL)*ZA1_INV, RMIN_JERLOV )
      ZR2(JL)   = MAX( PDO(JL)*PFACT(JL)*ZA2_INV, RMIN_JERLOV )
      PSWDK(JL) = RFAC_JERLOV(NJERLOV_KPP) * EXP(ZR1(JL)) &
               & + ( 1.0_JPRB - RFAC_JERLOV(NJERLOV_KPP) ) * EXP(ZR2(JL))
    ENDIF
  ENDDO

ENDIF

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_SWFRAC_MOD:KPP_SWFRAC',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_SWFRAC
END MODULE KPP_SWFRAC_MOD
