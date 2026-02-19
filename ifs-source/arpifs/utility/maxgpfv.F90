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

SUBROUTINE MAXGPFV(YDGEOMETRY,YDSURF,PWSMAX,PWDMAX)

!**** *MAXGPFV*  - COMPUTE MAXIMUM VALUE OF A GRIDPOINT FIELD

!     PURPOSE.
!     --------
!        COMPUTE MAXIMUM VALUE.

!**   INTERFACE.
!     ----------
!       *CALL* *MAXGPFV*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         YDSURF : model surface fields structure

!        OUTPUT:
!         PWSMAX    : maximum surface soil moisture
!         PWDMAX    : maximum deep    soil moisture

!        IMPLICIT ARGUMENTS
!        --------------------
!        See #include below.

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------


!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      R. El Khatib & J-F Estrade : 04-08-03 No comminication if 1 proc only.
!      M.Hamrud      01-Jul-2006  Revised surface fields
!      K.Yessad (Oct 2009) Use DIWRGRID for DM communications.
!      R. El Khatib : 23-Apr-2010 use diwgrid_mod instead of diwgrid
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      R. El Khatib : 02-Mar-2015 Rewrite
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMMP0       , ONLY : NPROC, MYPROC
USE MPL_MODULE   , ONLY : MPL_BROADCAST

IMPLICIT NONE

TYPE(GEOMETRY)  ,INTENT(IN)  :: YDGEOMETRY
TYPE(TSURF)     ,INTENT(IN)  :: YDSURF
REAL(KIND=JPRB) ,INTENT(OUT) :: PWSMAX 
REAL(KIND=JPRB) ,INTENT(OUT) :: PWDMAX 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZMAXB(YDGEOMETRY%YRDIM%NGPBLKS), ZMAX(NPROC)
INTEGER(KIND=JPIM) :: IST, IEND, IBLK, JROC, JKGLO
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MAXGPFV',0,ZHOOK_HANDLE)
ASSOCIATE(NPROMA=>YDGEOMETRY%YRDIM%NPROMA, NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS, NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, &
 & SP_RR=>YDSURF%SP_RR, SP_SB=>YDSURF%SP_SB, YSP_RR=>YDSURF%YSP_RR, YSP_SB=>YDSURF%YSP_SB)
!     ------------------------------------------------------------------

IST=1

! 1. Surface moisture

IF (YSP_RR%YW%LSET) THEN

  ! Compute local max value
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IEND,IBLK)
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBLK=(JKGLO-1)/NPROMA+1
    ZMAXB(IBLK)=MAXVAL(SP_RR(IST:IEND,YSP_RR%YW%MP0,IBLK))
  ENDDO
!$OMP END PARALLEL DO
  ZMAX(MYPROC)=MAXVAL(ZMAXB(:))

  ! Broadcast max values on each task
  DO JROC=1,NPROC
    CALL MPL_BROADCAST(ZMAX(JROC),KTAG=JROC,KROOT=JROC,CDSTRING='MAXGPFV:')
  ENDDO

  PWSMAX=MAXVAL(ZMAX(:))
  WRITE (NULOUT,*) ' MAXGPFV : MAX. SURFACE SOIL MOISTURE = ', PWSMAX

ENDIF

! 2. Deep soil moisture

IF (YSP_SB%YQ%LSET.AND. SIZE(SP_SB,DIM=2) >= 1) THEN

  ! Compute local max value
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(IEND,IBLK)
  DO JKGLO=1,NGPTOT,NPROMA
    IEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBLK=(JKGLO-1)/NPROMA+1
    ZMAXB(IBLK)=MAXVAL(SP_SB(IST:IEND,1,YSP_SB%YQ%MP0,IBLK))
  ENDDO
!$OMP END PARALLEL DO
  ZMAX(MYPROC)=MAXVAL(ZMAXB(:))

  ! Broadcast max values on each task
  DO JROC=1,NPROC
    CALL MPL_BROADCAST(ZMAX(JROC),KTAG=JROC,KROOT=JROC,CDSTRING='MAXGPFV:')
  ENDDO

  PWDMAX=MAXVAL(ZMAX(:))
  WRITE (NULOUT,*) ' MAXGPFV : MAX. DEEP SOIL MOISTURE = ', PWDMAX

ENDIF
!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('MAXGPFV',1,ZHOOK_HANDLE)

END SUBROUTINE MAXGPFV
