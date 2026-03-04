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

SUBROUTINE GPNORM2(YDGEOMETRY,KFIELDS,KSTAF,PGP,CDLABEL)

!**** *GPNORM2*  - ROUTINE FOR CALCULATING GRID POINT NORM OF 2-D FILEDS
!                  

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *GPNORM2(...)*

!        INPUT:

!        KFIELDS   - number of fields 
!        KSTAF     - first field 
!        PGP       - field

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      L. Isaksen *ECMWF*
!      Original : 96-07-07  Based on extgpf.F

!     MODIFICATIONS.
!     --------------
!      Y.Tremolet: 03-11-19 Special case when buffer not filled.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      M.Hamrud      01-Jul-2006 Revised surface fields
!      G.Mozdzynski  19-Sept-2008 gpnorm optimisations
!     -----------------------------------------------------------------

USE GEOMETRY_MOD ,ONLY : GEOMETRY
USE PARKIND1     ,ONLY : JPIM     ,JPRB
USE YOMHOOK      ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTAF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(:,:,:)
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: CDLABEL
REAL(KIND=JPRB) :: ZGP(YDGEOMETRY%YRDIM%NPROMA,KFIELDS,YDGEOMETRY%YRDIM%NGPBLKS)

INTEGER(KIND=JPIM) ::  IEND, IST, JF, JL, JKGLO, IBL
LOGICAL :: LLEVELS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "gpnorm1.intfb.h"

!     -----------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPNORM2',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT)
!     -----------------------------------------------------------------

!--- LOOP OVER NPROMA PACKETS

CALL GSTATS(1216,0)
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JKGLO,IST,IEND,IBL,JF,JL)
DO JKGLO=1,NGPTOT,NPROMA
  IST=1
  IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
  IBL=(JKGLO-1)/NPROMA+1
  DO JF=1,KFIELDS
    DO JL=IST,IEND
      ZGP(JL,JF,IBL) = PGP(JL,JF+KSTAF-1,IBL)
    ENDDO
  ENDDO
ENDDO
!$OMP END PARALLEL DO
CALL GSTATS(1216,1)

LLEVELS=.TRUE.
CALL GPNORM1(YDGEOMETRY,ZGP,KFIELDS,LLEVELS,CDLABEL)

!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM2',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM2
