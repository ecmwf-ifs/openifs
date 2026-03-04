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

!OCL NOEVAL
SUBROUTINE SURCORDI(YDGEOMETRY,YDDYN)

!**** *SURCORDI*   - Computation of RCORDIT, RCORDIH, RCORDIF

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SURCORDI

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ARPEGE/ALADIN DOCUMENTATION

!     Author.
!     -------
!        K. Yessad (CNRM/GMAP) after the old part 2.5 of SUHDF
!         (written by M. Hortal).
!        Original : August 2005

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMLUN   , ONLY : NULOUT
USE YOMSTA   , ONLY : RDTDZ1,RTSUR,RTTROP
USE YOMDYN   , ONLY : TDYN
USE YOMCST   , ONLY : RG, RD
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
INTEGER(KIND=JPIM) :: JLEV

REAL(KIND=JPRB) :: ZALPHA, ZPB, ZPBH, ZTB, ZTBH, ZZL  

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SURCORDI',0,ZHOOK_HANDLE)
ASSOCIATE(YDVAB=>YDGEOMETRY%YRVAB,YDVETA=>YDGEOMETRY%YRVETA,YDVFE=>YDGEOMETRY%YRVFE,YDSTA=>YDGEOMETRY%YRSTA, &
 & YDDIMV=>YDGEOMETRY%YRDIMV, &
 & YDLAP=>YDGEOMETRY%YRLAP,YDCSGLEG=>YDGEOMETRY%YRCSGLEG, &
 & YDCSGEOM=>YDGEOMETRY%YRCSGEOM,YDCSGEOM_NB=>YDGEOMETRY%YRCSGEOM_NB,YDGSGEOM=>YDGEOMETRY%YRGSGEOM, &
 & YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, YDSPGEOM=>YDGEOMETRY%YSPGEOM)
ASSOCIATE(NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & RCORDIF=>YDDYN%RCORDIF, &
 & RCORDIH=>YDDYN%RCORDIH, RCORDIT=>YDDYN%RCORDIT)
!     ------------------------------------------------------------------

!*       1.    Computation of RCORDIT, RCORDIH, RCORDIF:

! * Remark: this piece of code is fully consistent with
!   lvertfe=F, ndlnpr=0, laprxpk=T
!   where prehyd(l)=0.5*prehyd(lbar-1)+0.5*prehyd(lbar)
!   and not completely consistent with the other cases
!   for lvertfe,ndlnpr,laprxpk.

IF (NFLEVG > 1) THEN
  ZALPHA=-RDTDZ1*RD/RG
  RCORDIH(0)=0._JPRB
  DO JLEV=1,NFLEVG
    ZPB=0.5_JPRB*(YDVAB%VAH(JLEV)+YDVAB%VAH(JLEV-1)+(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))*VP00)
    ZPBH=YDVAB%VAH(JLEV)+YDVAB%VBH(JLEV)*VP00
    ZTB=RTSUR*(ZPB/VP00)**ZALPHA
    ZTBH=RTSUR*(ZPBH/VP00)**ZALPHA
    RCORDIF(JLEV)=0.5_JPRB*(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))*ZALPHA*ZTB*VP00/ZPB
    RCORDIH(JLEV)=YDVAB%VBH(JLEV)*ZALPHA*ZTBH*VP00/ZPBH
    IF (ZTB > RTTROP) THEN
      ZZL=0.5_JPRB*(YDVAB%VBH(JLEV)+YDVAB%VBH(JLEV-1))*ZALPHA*ZTB*VP00/ZPB
    ELSE
      ZZL=0._JPRB
    ENDIF
    RCORDIT(JLEV)=ZZL
  ENDDO
ELSE
  RCORDIT(1)=1._JPRB
ENDIF

! * Printings:
                                                                                
IF (NPRINTLEV > 0) THEN
  WRITE(NULOUT,'(A)')
  WRITE(NULOUT,'(A)') ' --- SURCORDI:'
  WRITE(NULOUT,'(A)') ' RCORDIT:'
  WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIT(JLEV), JLEV=1,NFLEVG)
  IF (NFLEVG > 1) THEN
    WRITE(NULOUT,'(A)') ' RCORDIH:'
    WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIH(JLEV), JLEV=0,NFLEVG)
    WRITE(NULOUT,'(A)') ' RCORDIF:'
    WRITE(NULOUT,'(4(1X,E19.13))') (RCORDIF(JLEV), JLEV=1,NFLEVG)
  ENDIF
  CALL FLUSH(NULOUT)
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURCORDI',1,ZHOOK_HANDLE)
END SUBROUTINE SURCORDI
