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

SUBROUTINE GL2LL(YDMP,KGLATNO,KGLONNO,KPE,KLOLATNO,KLOLONNO)

!****   gl2ll  ****  routine to convert from global lat/lon numbers
!                                          to local lat/lon numbers

!       purpose
!       -------

!**     Interface.
!       ----------
!         *call* *gl2ll(...) 

!        Explicit arguments : 
!        --------------------
!                kglatno  - global latitude number  (input)
!                kglonno  - global longitude number  (input)
!                kpe      - PE to whom the point belongs (output)
!                klolatno - local latitude number (output)
!                klolatno - local latitude number (output)

!        Implicit arguments :  None.
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!         Original   : 96-04-04
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Dec-2003 CY28R1 Cleaning
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib 16-May-2019 optimize memory access in NGPSET2PE
!----------------------------------------------------------------------
USE YOMMP    , ONLY : TMP
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NGPSET2PE, N_REGIONS, N_REGIONS_NS

IMPLICIT NONE

TYPE(TMP)         ,INTENT(IN)    :: YDMP
INTEGER(KIND=JPIM),INTENT(IN)    :: KGLATNO 
INTEGER(KIND=JPIM),INTENT(IN)    :: KGLONNO 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KPE 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLOLATNO 
INTEGER(KIND=JPIM),INTENT(OUT)   :: KLOLONNO 

INTEGER(KIND=JPIM) :: IASET, IBSET, IEND, IFIRST, IGL1, IGL2, IGLOFF, IST, JA, JB, JGL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

!----------------------------------------------------------------------
! find first a-set

IF (LHOOK) CALL DR_HOOK('GL2LL',0,ZHOOK_HANDLE)
ASSOCIATE(NFRSTLAT=>YDMP%NFRSTLAT, NLSTLAT=>YDMP%NLSTLAT, NONL=>YDMP%NONL, &
 & NPTRFRSTLAT=>YDMP%NPTRFRSTLAT, NSTA=>YDMP%NSTA)
IASET=1
DO JA=1,N_REGIONS_NS
  IGLOFF=NPTRFRSTLAT(JA)
  IF(JA > 1) THEN
    IF(NLSTLAT(JA-1) == NFRSTLAT(JA)) THEN
      IFIRST=NLSTLAT(JA-1)
      IF(KGLATNO > IFIRST) THEN
        IASET=JA
      ELSEIF (KGLATNO == IFIRST) THEN
        DO JB=1,N_REGIONS(JA)
          IF(NONL(IGLOFF,JB) > 0.AND.NSTA(IGLOFF,JB) <= KGLONNO) IASET=JA
        ENDDO
      ENDIF
    ELSE
      IFIRST=NLSTLAT(JA-1)+1
      IF(KGLATNO >= IFIRST) IASET=JA
    ENDIF
  ENDIF
ENDDO

!  find the b-set

IBSET=-999
IGLOFF=NPTRFRSTLAT(IASET)
IGL1  =NFRSTLAT(IASET)
IGL2  =NLSTLAT(IASET)
DO JGL=IGL1,IGL2
  IF(JGL == KGLATNO) THEN
    DO JB=1,N_REGIONS(IASET)
      IST=NSTA(IGLOFF+JGL-IGL1,JB)
      IEND=NSTA(IGLOFF+JGL-IGL1,JB)+NONL(IGLOFF+JGL-IGL1,JB)-1
      IF(KGLONNO >= IST.AND.KGLONNO <= IEND) THEN
        IBSET=JB
        KLOLONNO=KGLONNO-IST+1
      ENDIF
    ENDDO
  ENDIF
ENDDO

IF(IASET <= 0.OR.IASET > N_REGIONS_NS.OR.IBSET <= 0 &
 & .OR.IBSET > N_REGIONS(IASET)) CALL ABOR1('gl2ll : error')  
KLOLATNO=KGLATNO-NFRSTLAT(IASET)+1

!     find the PE that owns the grid point

KPE=NGPSET2PE(IBSET,IASET)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GL2LL',1,ZHOOK_HANDLE)
END SUBROUTINE GL2LL

