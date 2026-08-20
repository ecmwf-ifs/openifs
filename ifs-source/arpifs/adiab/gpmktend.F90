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

SUBROUTINE GPMKTEND(YDGEOMETRY,YGFL,KST,KEND,PTSPHY,&
 & PUT1,PVT1 ,PTT1,&
 & PGFLT1,&
 & PUT9,PVT9 ,PTT9,&
 & PGFL,PGFL_DYN)  

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOETHF       , ONLY : YDTHF=>YRTHF ! allows use of included functions. REK.
USE YOMCST       , ONLY : YDCST=>YRCST ! allows use of included functions. REK.

!**** *GPMKTEND - Convert t+dt values to tendencies 

!     Purpose. Convert t+dt values to tendencies before calling
!     -------- ECMWF physics

!**   Interface.
!     ----------
!        *CALL* *GPMKTEND(...)*

!        Explicit arguments : see CPGLAG, all names the same
!        -------------------- except for DOCTOR norm modifications

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        None

!     Author.
!     -------
!      Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K. Yessad and D. Salmond (Feb 2006): adapt to LPC_FULL.
!      N. Wedi and D. Salmond 06-05-01 : liquid water and ice split moved from callpar
!      K. Yessad (Dec 2008): remove dummy CDLOCK
!      A. Tompkins 2007-15-06 : remove liquid water and ice split 
!      R. El Khatib 08-Jul-2022 Contribution to the encapsulation of YOMCST and YOETHF
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSPHY 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT9(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PGFL_DYN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM1)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL
REAL(KIND=JPRB) :: ZTS_R
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPMKTEND',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, NFLEVG=>YDDIMV%NFLEVG,&
 & RLVTT=>YDCST%RLVTT, RLSTT=>YDCST%RLSTT, RTT=>YDCST%RTT, &
 & R2ES=>YDTHF%R2ES, R3LES=>YDTHF%R3LES, R3IES=>YDTHF%R3IES, R4LES=>YDTHF%R4LES, &
 & R4IES=>YDTHF%R4IES, R5LES=>YDTHF%R5LES, R5IES=>YDTHF%R5IES, R5ALVCP=>YDTHF%R5ALVCP, &
 & R5ALSCP=>YDTHF%R5ALSCP, RALVDCP=>YDTHF%RALVDCP, RALSDCP=>YDTHF%RALSDCP, &
 & RTWAT=>YDTHF%RTWAT, RTICE=>YDTHF%RTICE, RTICECU=>YDTHF%RTICECU, RTWAT_RTICE_R=>YDTHF%RTWAT_RTICE_R, &
 & RTWAT_RTICECU_R=>YDTHF%RTWAT_RTICECU_R)
!     ------------------------------------------------------------------

ZTS_R=1.0_JPRB/PTSPHY
DO JLEV=1,NFLEVG
  DO JROF=KST,KEND
    PUT1(JROF,JLEV)=(PUT1(JROF,JLEV)-PUT9(JROF,JLEV))*ZTS_R
    PVT1(JROF,JLEV)=(PVT1(JROF,JLEV)-PVT9(JROF,JLEV))*ZTS_R
    PTT1(JROF,JLEV)=(PTT1(JROF,JLEV)-PTT9(JROF,JLEV))*ZTS_R
  ENDDO
ENDDO

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LT1) THEN
    DO JLEV=1,NFLEVG
      DO JROF=KST,KEND
        PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)=(PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)-&
         & PGFL(JROF,JLEV,YCOMP(JGFL)%MP9_PH))*ZTS_R  
!!! Save Dyn tendencies for CPQTUV
        PGFL_DYN(JROF,JLEV,YCOMP(JGFL)%MP1)=PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)
      ENDDO
    ENDDO
  ENDIF
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPMKTEND',1,ZHOOK_HANDLE)
END SUBROUTINE GPMKTEND
