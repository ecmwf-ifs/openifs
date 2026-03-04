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

SUBROUTINE READMASK(YDGEOMETRY,CDFILE,KGPTOT,PBUF)

!**** *READMASK*

!     PURPOSE.
!     --------

!     READ THE FIRST RECORD OF A FILE FOR VARIOUS MASKINGS

!**   INTERFACE.
!     ----------

!     CALL READMASK

!        EXPLICIT ARGUMENTS :
!        --------------------

!     METHOD.
!     -------

!     ECCODES/GRIB_API + DISGRID_SEND/DISGRID_RECV
   
!     EXTERNALS.
!     ----------


!     AUTHORS.
!     --------

!     K. MOGENSEN (September 2017): Intial version
       
!     MODIFICATIONS.
!     --------------
!     ------------------------------------------------------------------

USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1, ONLY : JPRD, JPIM, JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN  , ONLY : NULOUT
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE YOMMP0, ONLY : MYPROC
USE DISGRID_MOD, ONLY : DISGRID_SEND, DISGRID_RECV
USE GRIB_API_INTERFACE, ONLY : IGRIB_OPEN_FILE, IGRIB_NEW_FROM_FILE,&
 & IGRIB_GET_VALUE, IGRIB_RELEASE, IGRIB_CLOSE_FILE
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
CHARACTER(LEN=*), INTENT(IN) :: CDFILE
INTEGER(KIND=JPIM), INTENT(IN) :: KGPTOT
REAL(KIND=JPRB), DIMENSION(KGPTOT), INTENT(OUT) :: PBUF

!     ------------------------------------------------------------------
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB), ALLOCATABLE :: ZBUF(:)
INTEGER :: IUNITC,IGRIB,IVALS
!     ------------------------------------------------------------------

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READMASK',0,ZHOOK_HANDLE)

! READ FILE
! BEWARE THAT THIS IS NO ERROR CHECKING BESIDES THAT THE
! NUMBER OF POINT IS CORRECT.
IF(MYPROC==1) THEN
   WRITE(NULOUT,*)'READING FROM FILE :',TRIM(CDFILE)
   CALL IGRIB_OPEN_FILE(IUNITC,TRIM(CDFILE),'r')
   CALL IGRIB_NEW_FROM_FILE(IUNITC,IGRIB)
   CALL IGRIB_GET_VALUE(IGRIB,'numberOfCodedValues',IVALS)
   IF (IVALS/=YDGEOMETRY%YRGEM%NGPTOTG) THEN
      CALL ABOR1('WRONG NUMBER OF POINTS IN '//TRIM(CDFILE))
   ENDIF
   ALLOCATE(ZBUF(YDGEOMETRY%YRGEM%NGPTOTG))
   CALL IGRIB_GET_VALUE(IGRIB,'values',ZBUF)
   CALL IGRIB_RELEASE(IGRIB)
   CALL IGRIB_CLOSE_FILE(IUNITC)
   CALL DISGRID_SEND(YDGEOMETRY,1,ZBUF,1,PBUF)
   DEALLOCATE(ZBUF)
ELSE
   CALL DISGRID_RECV(YDGEOMETRY,1,1,PBUF,1)
ENDIF
!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('READMASK',1,ZHOOK_HANDLE)
END SUBROUTINE READMASK
