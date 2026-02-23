! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPNORM_GFL(YDGEOMETRY,YDGFL,LDPRINT_TL,PGFLNORMS)

!**** *GPNORM_GFL*  - ROUTINE FOR CALCULATING GRID POINT NORM OF
!                     GFL VARIABLES

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *GPNORM_GFL*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!       F. Vana (after the IFS version)
!       Original : 31 Jan 2017

!     MODIFICATIONS.
!     --------------
!
!     -----------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOM_YGFL , ONLY : YGFL

IMPLICIT NONE

TYPE(GEOMETRY)  , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)      , INTENT(IN)    :: YDGFL
LOGICAL         , INTENT(IN) , OPTIONAL :: LDPRINT_TL
REAL(KIND=JPRB) , INTENT(OUT), OPTIONAL :: PGFLNORMS(3,YDGEOMETRY%YRDIMV%NFLEVG)


INTEGER(KIND=JPIM) :: JGFL
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPNORM_GFL',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS,YCOMP=>YGFL%YCOMP,NFLEVG=>YDDIMV%NFLEVG,GFL=>YDGFL%GFL)


IF (PRESENT(PGFLNORMS)) THEN
  PGFLNORMS (1,1:NFLEVG) = 0._JPRB
  PGFLNORMS (2,1:NFLEVG) = 99999._JPRB
  PGFLNORMS (3,1:NFLEVG) =-99999._JPRB
ENDIF

DO JGFL=1,NUMFLDS

  WRITE(NULOUT,*) 'GPNORM ',YCOMP(JGFL)%CNAME, &
    & ' : AVERAGE =  ',SUM(GFL(1,1:NFLEVG,JGFL,1))/NFLEVG, &
    & ' MINIMUM = ',MINVAL(GFL(1,1:NFLEVG,JGFL,1)), &
    & ' MAXIMUM = ',MAXVAL(GFL(1,1:NFLEVG,JGFL,1))

  IF (PRESENT(PGFLNORMS)) THEN
    PGFLNORMS (1,1:NFLEVG) = PGFLNORMS (1,1:NFLEVG) + GFL(1,1:NFLEVG,JGFL,1)/REAL(NUMFLDS,JPRB)
    PGFLNORMS (2,1:NFLEVG) = MIN(PGFLNORMS (2,1:NFLEVG),GFL(1,1:NFLEVG,JGFL,1))
    PGFLNORMS (3,1:NFLEVG) = MAX(PGFLNORMS (3,1:NFLEVG),GFL(1,1:NFLEVG,JGFL,1))
  ENDIF

ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM_GFL',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM_GFL
