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

SUBROUTINE GPRT(LDSPRT,KPROMA,KSTART,KEND,KLEV,PRD,PRV,PR,&
 & PT,PTL,PTM,PQL,PQM,PRT,PRTL,PRTM,PRL,PRM)  

!*    *GPRT* 

!     Purpose
!     -------  To calculate RT and its derivates 
!     Interface
!     ---------

!     Explicit arguments
!     ------------------
!     Input:
!    -------
!              LDSPRT   : .TRUE. if PTL and PTM already contain
!                         the derivatives of 'TV')
!              KPROMA   : Horizontal dimension
!              KSTART   : Start index
!              KEND     : End index
!              KLEV     : number of levels
!              PRD      : Rd
!              PRV      : Rv
!              PR       : R
!              PT       : T
!              PTL ,PTM : Horizontal derivatives of T
!              PQL, PQM : Horizontal derivatives of q
!     Output:
!    --------
!              PRT       : RT
!              PRTL,PRTM : Horizontal derivatives of RT

!     Author
!     ------
!           J.Boutahar *MAROC-METEO*
!     Modifications
!     -------------
!           Original: 97/06/06
!        C. Fischer 02-06-27 : cdlock
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     K. Yessad (Dec 2008): remove dummy CDLOCK
!----------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

!----------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
LOGICAL           ,INTENT(IN)    :: LDSPRT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQM(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRT(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRTL(KPROMA,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRTM(KPROMA,KLEV) 
REAL(KIND=JPRB),OPTIONAL  ,INTENT(IN)    :: PRL(KPROMA,KLEV) 
REAL(KIND=JPRB),OPTIONAL  ,INTENT(IN)    :: PRM(KPROMA,KLEV) 

!----------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JLH
LOGICAL :: LLRDERS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPRT',0,ZHOOK_HANDLE)
!----------------------------------------------------------
!---------------------------------------------------------------
!*       1. Compute RT and its derivatives
LLRDERS = .FALSE.
IF(PRESENT(PRL)) LLRDERS = .TRUE.
DO JLEV=1,KLEV
  DO JLH=KSTART,KEND
    PRT(JLH,JLEV)=PR(JLH,JLEV)*PT(JLH,JLEV)
    IF (LDSPRT) THEN
      PRTL(JLH,JLEV)=PRD*PTL(JLH,JLEV)
      PRTM(JLH,JLEV)=PRD*PTM(JLH,JLEV)
    ELSEIF(LLRDERS) THEN
      PRTL(JLH,JLEV) = PRL(JLH,JLEV)*PT(JLH,JLEV) &
      & +PR(JLH,JLEV)*PTL(JLH,JLEV)  
      PRTM(JLH,JLEV) = PRM(JLH,JLEV)*PT(JLH,JLEV) &
      & +PR(JLH,JLEV)*PTM(JLH,JLEV)  
    ELSE
      PRTL(JLH,JLEV)=(PRV-PRD)*PT(JLH,JLEV)*PQL(JLH,JLEV)&
       & +PR(JLH,JLEV)*PTL(JLH,JLEV)  
      PRTM(JLH,JLEV)=(PRV-PRD)*PT(JLH,JLEV)*PQM(JLH,JLEV)&
       & +PR(JLH,JLEV)*PTM(JLH,JLEV)  
    ENDIF
  ENDDO
ENDDO

!----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPRT',1,ZHOOK_HANDLE)
END SUBROUTINE GPRT
