! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPNORM_GMV(YDGEOMETRY,YDGMV,PSPNORMS,PGMVNORMS)

!**** *GPNORM_GMV*  - ROUTINE FOR CALCULATING GRID POINT NORM OF
!                     GMV VARIABLES

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!        *CALL* *GPNORM_GMV*

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
USE YOMGMV       , ONLY : TGMV
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT

IMPLICIT NONE

TYPE(GEOMETRY) , INTENT(IN)   :: YDGEOMETRY
TYPE(TGMV)     , INTENT(IN)   :: YDGMV
REAL(KIND=JPRB) , INTENT(OUT), OPTIONAL :: PSPNORMS(3,1)    !! for surface pressure
REAL(KIND=JPRB) , INTENT(OUT), OPTIONAL :: PGMVNORMS(3,YDGEOMETRY%YRDIMV%NFLEVG) 

INTEGER(KIND=JPIM), PARAMETER :: JPNFLDS = 3
INTEGER(KIND=JPIM) :: IMPT(JPNFLDS)
CHARACTER(LEN=16)  :: CLNAME(JPNFLDS)

INTEGER(KIND=JPIM) :: JF
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GPNORM_GMV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & GMV=>YDGMV%GMV, GMVS=>YDGMV%GMVS, YT0=>YDGMV%YT0)


WRITE(NULOUT,*) 'GPNORM   SURFACE PRESSURE',GMVS(1,YT0%MSP,1)
IF (PRESENT(PSPNORMS)) PSPNORMS(1:3,1)=GMVS(1,YT0%MSP,1)


IMPT(1)  = YT0%MU
IMPT(2)  = YT0%MV
IMPT(3)  = YT0%MT

CLNAME(1) = 'U VELOCITY'
CLNAME(2) = 'V VELOCITY'
CLNAME(3) = 'TEMPERATURE'

IF (PRESENT(PGMVNORMS)) THEN
  PGMVNORMS (1,1:NFLEVG) = 0._JPRB
  PGMVNORMS (2,1:NFLEVG) = 99999._JPRB
  PGMVNORMS (3,1:NFLEVG) =-99999._JPRB
ENDIF

DO JF = 1,JPNFLDS 
  WRITE(NULOUT,*) 'GPNORM ',CLNAME(JF),&
    & ' : AVERAGE =  ',SUM(GMV(1,1:NFLEVG,IMPT(JF),1))/NFLEVG, &
    & ' MINIMUM = ',MINVAL(GMV(1,1:NFLEVG,IMPT(JF),1)), &
    & ' MAXIMUM = ',MAXVAL(GMV(1,1:NFLEVG,IMPT(JF),1))
  IF (PRESENT(PGMVNORMS)) THEN
    PGMVNORMS (1,1:NFLEVG) = PGMVNORMS (1,1:NFLEVG) + GMV(1,1:NFLEVG,IMPT(JF),1)/REAL(JPNFLDS,JPRB)
    PGMVNORMS (2,1:NFLEVG) = MIN(PGMVNORMS (2,1:NFLEVG),GMV(1,1:NFLEVG,IMPT(JF),1))
    PGMVNORMS (3,1:NFLEVG) = MAX(PGMVNORMS (3,1:NFLEVG),GMV(1,1:NFLEVG,IMPT(JF),1))
  ENDIF
ENDDO

!     -----------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPNORM_GMV',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM_GMV
