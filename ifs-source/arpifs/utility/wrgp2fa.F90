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

SUBROUTINE WRGP2FA(YDGEOMETRY,YDRIP,KFIELDS,PBUF,YDFLDSC,YDFACTX,CDFIC,KTAG)
!  ndistio     : if set to 1, then fields with be written in separate files
!                (1 file/nstrout proc).

!**** *WRGP2FA*  - Write gridpoint fields to *FA* file

!     Purpose.
!     --------
!        To write local gridpoint fields to *FA* file. The file is expected
!        to be opened and close out of this subroutine, which, in counterpart,
!        handles the communications.

!**   Interface.
!     ----------
!        *CALL* *WRGP2FA*

!        Explicit arguments :
!        --------------------

!            PBUF    : gridpoint fields
!            KFIELDS : number of gridpoint fields in array
!            KUNIT   : file logical unit number 
!            KTAG    : Tag used for gathering data on MPI #1 in random order

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           FA routines
!           Processor communication routines
!           DISGRID

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!         RYAD EL KHATIB *METEO-FRANCE*

!     Modifications.
!     --------------
!        ORIGINAL : 01-0130 from RDPXFA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F.Bouyssel    05-11-09 Modifications for cloud fields
!        R. El Khatib 12-Mar-2010 Distribute bunches of fields
!                                 + partial merge with wrmlppa
!        O. Vignes    08-Jun-2010 Merged send/recv loops, OpenMP,
!                                 GATHERV (after S. Saarinen)
!        P. Marguinaud 18-May-2010 Write fields in separate files (NDISTIO(1)==1)
!                                  Keep norms output in a single file
!        P. Marguinaud 01-Jan-2011 IO server & misc NDISTIO settings
!        P. Marguinaud 06-Jul-2011 If PUNDEF present, replace all PUNDEF values
!                                  with field average
!        R. El Khatib 06-Apr-2012 NECSX compiler workaround (NECSX is not thread-safe, anyway)
!        P.Marguinaud 26-Apr-2012 Remove unused code
!        P.Marguinaud 11-Sep-2012 Refactor using IOMULTIBUF_MOD, IOFLDDESC_MOD, MFIOOPTS_MOD
!                                 WRFLDCW_MOD, WRGATHFLNM + Compress data with OpenMP
!        P.Marguinaud 10-Oct-2013 Use YDFACTX, GATH_GRID & EGATH_GRID;
!                                 remove various options
!        T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK ,DR_HOOK, JPHOOK

USE YOMMP0             , ONLY : NPROC, MYPROC, NSTROUT
USE YOMCT0             , ONLY : LELAM
USE YOMOPH0            , ONLY : YMDLOPH
USE FPGPNORM_MOD       , ONLY : FPGPNORMX

USE IOFLDDESC_MOD      , ONLY : IOFLDDESC
USE IOCPTDESC_MOD      , ONLY : IOCPTDESC

USE FACTX_MOD          , ONLY : FACTX, COMPACTFLD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY),     INTENT (IN)           :: YDGEOMETRY
TYPE (TRIP),         INTENT (IN)           :: YDRIP
INTEGER (KIND=JPIM), INTENT (IN)           :: KFIELDS
REAL (KIND=JPRB),    INTENT (IN)           :: PBUF (YDGEOMETRY%YRGEM%NGPTOT,KFIELDS,1)
TYPE (IOFLDDESC),    INTENT (IN)           :: YDFLDSC (:) 
TYPE (FACTX),        INTENT (INOUT)        :: YDFACTX
CHARACTER(LEN=*),    INTENT (IN)           :: CDFIC
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KTAG

#include "wrgathflnm.intfb.h"
#include "sumpioh.intfb.h"

#include "gath_grid.h"
#include "egath_grid.h"


!     ------------------------------------------------------------------

!     INFD   : number of fields in the current chunk for each processor
!     INFL   : number of fields in the current chunk for myproc
!     IFLDOFF: fields offset after distribution among NSTROUT proc, for each processor


TYPE (IOCPTDESC), ALLOCATABLE :: YLCPDSC (:)
REAL(KIND=JPRB),  ALLOCATABLE :: ZGVALCO(:,:)
REAL(KIND=JPRB),  ALLOCATABLE :: ZGPG(:,:)

INTEGER(KIND=JPIM) :: JROC

INTEGER(KIND=JPIM) :: INFD (NPROC), IFLDOFF (NPROC)
INTEGER(KIND=JPIM) :: ITO (KFIELDS)
INTEGER(KIND=JPIM) :: INFL, IFLDG1, IFLDG2
INTEGER(KIND=JPIM) :: IDIMGVAL

CHARACTER(LEN=LEN(CDFIC)) :: CLFIC(1)
!     ------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRGP2FA',0,ZHOOK_HANDLE)




CALL SUMPIOH(NPROC,NSTROUT,KFIELDS,INFD,IFLDOFF)

INFL = INFD(MYPROC)

ALLOCATE (ZGPG (YDGEOMETRY%YRGEM%NGPTOTG, INFL), YLCPDSC (INFL))

ITO = 0
DO JROC = 1, NPROC
  ITO (IFLDOFF (JROC)+1:IFLDOFF (JROC)+INFD (JROC)) = JROC
ENDDO

IF (LELAM) THEN
  CALL EGATH_GRID (KFGATHG=KFIELDS,KTO=ITO,PGP=PBUF,PGPG=ZGPG,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ELSE
  CALL GATH_GRID (KFGATHG=KFIELDS,KTO=ITO,PGP=PBUF,PGPG=ZGPG,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ENDIF

IFLDG1=IFLDOFF(MYPROC)+1
IFLDG2=IFLDOFF(MYPROC)+INFL

CALL FPGPNORMX (ZGPG, YLCPDSC%ZAVE, YLCPDSC%ZMIN, YLCPDSC%ZMAX, &
              & LDUNDF=YDFLDSC (IFLDG1:IFLDG2)%LUNDF,           &
              & PUNDEF=YDFLDSC (IFLDG1:IFLDG2)%XUNDF)     

IF (YDFACTX%YIOOPTS%LDOCP) THEN

  IDIMGVAL=MAXVAL (YDFLDSC%NSIZPK)

  ALLOCATE (ZGVALCO (IDIMGVAL, INFL))

  CALL COMPACTFLD (YMDLOPH(1)%CFPCA, YDFACTX, YDFLDSC (IFLDG1:IFLDG2), ZGPG, YLCPDSC, ZGVALCO)

  CLFIC(1)=CDFIC
  CALL WRGATHFLNM (YDGEOMETRY%YVABIO, 1_JPIM, INFL, IDIMGVAL, INFD, IFLDOFF, ZGVALCO, YLCPDSC, &
                 & YMDLOPH, CLFIC, LDNORM=.TRUE., YDFACTX=YDFACTX, KTAG=KTAG, &
                 & CDNORMTITLE='GPNORMS OF FIELDS TO BE WRITTEN OUT ON FILE :',YDGEOMETRY=YDGEOMETRY)

  DEALLOCATE (ZGVALCO)

ENDIF

DEALLOCATE (ZGPG, YLCPDSC)



IF (LHOOK) CALL DR_HOOK('WRGP2FA',1,ZHOOK_HANDLE)

END SUBROUTINE WRGP2FA

