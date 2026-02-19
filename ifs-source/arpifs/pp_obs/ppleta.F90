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

SUBROUTINE PPLETA(YDVAB,YDCVER,KPROMA,KST,KND,KOPLEV,PXLEV,PS,PRESF,PRESHB,PRESHT)  

!**** *PPLETA*  - FULL-POS upper air pressure post-processing from eta system

!     PURPOSE.
!     --------
!        To compute upper air pressures from a given vertical hybrid system 
!         and the surface pressure.

!**   INTERFACE.
!     ----------
!       *CALL* *PPLETA*

!        EXPLICIT ARGUMENTS
!        --------------------
!        * INPUT:
!        KPROMA  : horizontal dimension
!        KST     : start of work
!        KND     : depth of work
!        KOPLEV  : numer of levels among all eta system
!        PXLEV   : indexes of output levels (as a real array)
!        PS      : output surface pressure
!        * OUTPUT:
!        PRESF   : full level pressures
!        PRESHB  : true half level pressures (on bottom full levels) 
!        PRESHT  : pseudo - half level pressures (on top full levels) 

!        IMPLICIT ARGUMENTS
!        ------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMVERT  , ONLY : TVAB
USE YOMCVER , ONLY : TCVER
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TVAB)        ,INTENT(IN)    :: YDVAB
TYPE(TCVER)       ,INTENT(IN)    :: YDCVER
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOPLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PXLEV(KOPLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PS(KPROMA) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRESF(KPROMA,KOPLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRESHB(KPROMA,KOPLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRESHT(KPROMA,KOPLEV) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZPRESF(KPROMA,UBOUND(YDVAB%VBH,DIM=1)),ZPRESH(KPROMA,0:UBOUND(YDVAB%VBH,DIM=1))

INTEGER(KIND=JPIM) :: JL
INTEGER(KIND=JPIM) :: IPTRLEV(KOPLEV) ! IPTRLEV : absolute indexes of levels
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "gphpre.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPLETA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    COMPUTE OUTPUT PRESSURE LEVELS
!              ------------------------------

!        1.1   PRESSURES OF OUTPUT ETA SYSTEM

ZPRESH(KST:KND,UBOUND(YDVAB%VBH,DIM=1))=PS(KST:KND)

CALL GPHPRE(KPROMA,UBOUND(YDVAB%VBH,DIM=1),KST,KND,YDVAB,YDCVER,ZPRESH,PRESF=ZPRESF)

!        1.2   SELECT EFFECTIVE OUTPUT UPPER AIR PRESSURES

DO JL=1,KOPLEV
  IPTRLEV(JL)=INT(PXLEV(JL),JPIM)
  PRESF (KST:KND,JL)=ZPRESF(KST:KND,IPTRLEV(JL))
  PRESHB(KST:KND,JL)=ZPRESH(KST:KND,IPTRLEV(JL))
  PRESHT(KST:KND,JL)=ZPRESH(KST:KND,IPTRLEV(JL)-1)
ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPLETA',1,ZHOOK_HANDLE)
END SUBROUTINE PPLETA
