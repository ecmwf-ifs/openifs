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

SUBROUTINE GPNORM3(YDGEOMETRY,PGP,KDIM3,KRET,CDLABEL,KLEVS,LDLEVELS)

!**** *GPNORM3*  - ROUTINE FOR CALCULATING GRID POINT NORM OF
!                  ARRAYS PGP(NPROMA,NFLEVG,KDIM3,NGPBLKS)

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *GPNORM3(...)*

!        INPUT:

!        KRET      - field number whose norm is to be computed
!        PGP       - Input array

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------

!        M. Hortal  *ECMWF* after L. Isaksen's GPNORM2

!     MODIFICATIONS.
!     --------------
!       Original : 02-03-19  
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        G.Mozdzynski  19-Sept-2008 gpnorm optimisations

!     -----------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : NSPPR

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM3 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(:,:,:,:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KRET 
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: CDLABEL
INTEGER(KIND=JPIM),INTENT(IN),OPTIONAL :: KLEVS
LOGICAL           ,INTENT(IN),OPTIONAL :: LDLEVELS
 
REAL(KIND=JPRB),ALLOCATABLE :: ZGP(:,:,:)

INTEGER(KIND=JPIM) :: IEND, IST, JF, JL, JKGLO, IBL, ILEVS
LOGICAL            :: LLEVELS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
#include "gpnorm1.intfb.h"


!     -----------------------------------------------------------------

!--- LOOP OVER NPROMA PACKETS

IF (LHOOK) CALL DR_HOOK('GPNORM3',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT)
IF(KRET > KDIM3) THEN
  CALL ABOR1(' KRET > KDIM3 in GPNORM3')
ENDIF
IF(PRESENT(KLEVS)) THEN
  ILEVS = KLEVS
ELSE
  ILEVS = NFLEVG
ENDIF
IF(PRESENT(LDLEVELS)) THEN
  LLEVELS = LDLEVELS
ELSE
  LLEVELS = NSPPR>0
ENDIF
  
ALLOCATE(ZGP(NPROMA,ILEVS,NGPBLKS))
CALL GSTATS(1223,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IST,IEND,JL,JF,IBL)
DO JKGLO=1,NGPTOT,NPROMA
  IST=1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  DO JF=1,ILEVS
    DO JL=IST,IEND
      ZGP(JL,JF,IBL) = PGP(JL,JF,KRET,IBL)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1223,1)

CALL GPNORM1(YDGEOMETRY,ZGP,ILEVS,LLEVELS,CDLABEL)
DEALLOCATE(ZGP)

!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM3',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM3
