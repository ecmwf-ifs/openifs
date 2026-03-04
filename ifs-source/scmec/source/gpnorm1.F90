! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE GPNORM1(YDGEOMETRY,PGP,KFIELDS,LDLEVELS,CDLABEL,PNORMS)

!**** *gpnorm1 * - Routine to calculate reproducible grid point norm

!     Purpose.
!     --------
! Routine to calculate reproducible grid point norm.

!**   Interface.
!     ----------
!        *CALL* *GPNORM1(...) *

!        Explicit arguments :
!        --------------------
!        PGP     : if nproc is > 1 then PGP is the distributed field
!        KFIELDS  : number of levels
!        LDLEVELS : =F Only print average for all levels. Min and max are global

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      L. Isaksen *ECMWF*
!      Original : 96-07-07

!     Modifications.
!     --------------
!      Modified : 02-05-10 by S.Checker : changed MIMIMUM to MINIMUM
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      L.Isaksen : Support printing of average for all levels
!      G.Mozdzynski : 25-Oct-2005 use all available tasks and threads for norms
!      Y.Tremolet    03-Aug-2005 Optional output argument
!      N.Wedi    : Proper area-weighted mean
!      G.Mozdzynski  19-Sept-2008 gpnorm optimisations
!      A.Bogatchev: 12-Jun-2009 call to egpnorm_trans
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0       , ONLY : NPRINTLEV, MYPROC, NPROC
USE YOMCT0       , ONLY : LALLOPR, LELAM
USE YOMLUN       , ONLY : NULOUT
USE MPL_MODULE   , ONLY : MPL_BROADCAST

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGP(YDGEOMETRY%YRDIM%NPROMA,KFIELDS,YDGEOMETRY%YRDIM%NGPBLKS) 
LOGICAL           ,INTENT(IN)    :: LDLEVELS
CHARACTER(LEN=*)  ,INTENT(IN),OPTIONAL :: CDLABEL
REAL(KIND=JPRB), INTENT(OUT), OPTIONAL :: PNORMS(3,KFIELDS) 

!     ------------------------------------------------------------------

LOGICAL :: LLP
REAL(KIND=JPRB),ALLOCATABLE :: ZGP(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZAVE(:), ZMAX(:), ZMIN(:)
REAL(KIND=JPRB)   :: ZMAXBLK(YDGEOMETRY%YRDIM%NGPBLKS), ZMINBLK(YDGEOMETRY%YRDIM%NGPBLKS)
REAL(KIND=JPRB)   :: ZAVEALL, ZMAXALL, ZMINALL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

INTEGER(KIND=JPIM) :: IGPNMASTER,IFIELDS
INTEGER(KIND=JPIM) :: JF,JL,ITAG,JKGLO,IST,IEND,IBL

CHARACTER(LEN=16) :: CLABEL
LOGICAL :: LLAVE_ONLY

!     ------------------------------------------------------------------


IF (LHOOK) CALL DR_HOOK('GPNORM1',1,ZHOOK_HANDLE)
END SUBROUTINE GPNORM1
