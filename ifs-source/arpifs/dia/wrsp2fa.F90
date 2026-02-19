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

SUBROUTINE WRSP2FA(YDGEOMETRY,YDXFU,YDRIP,KFIELDS,PBUF,YDFLDSC,YDFACTX,CDFIC,KTAG)

!**** *WRSP2FA*  - Write the spectral fields to FA (Arpege or Aladin)

!     Purpose.
!     --------
!         Write the spectral fields of the model to ARPEGE ALADIN file.

!**   Interface.
!     ----------
!        *CALL* *WRSP2FA(.....)*

!        Explicit arguments :
!        --------------------
!        SPVOR etc. - spectral fields

!        Implicit arguments :
!        --------------------
!        See modules above.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        see below

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Note de travail ARPEGE Nr 12 et 17

!     Author.
!     -------
!        Stjepan IVATEK-SAHDAN and Ryad El Khatib *METEO-FRANCE*

!     Modifications.
!     --------------
!   Original : 01-04-17 from suspeca
!   Modified 02-03-08 C. Fischer - rescaling of new NH variables
!   Modified 02-09-30 P. Smolikova - rescaling of d4 in NH
!   R. El Khatib : 03-08-05 gfl+remove dummies since it it called only here.
!   O.Spaniel    : 03-04-15 cleaning-see interface wrspeca.h
!   R. El Khatib : 03-04-17 Cleanups
!   Modified 03-05-01 A. Bogatchev - check sizes and gnhpdvdconv
!   G. Hello : supress the call to gnhpdvdconv
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   P.Termonia: 03-09-01 writing RMCUFFP
!   D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!   Y. Bouteloup : 05-10-18 YCVGQ not write even it is spectral
!   Y. Seity     : 23-06 2006: remove YCVGQ specificities and replace by more
!                  general tests on LREQOUT
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   R. El Khatib : 15-Sep-2008 I/O savings
!   R. El Khatib : 07-Aug-2009 Bugfix for GRIB 1 encoding
!   R. El Khatib : 30-Mar-2010 reduce overhead of call to spreord
!   P. Marguinaud: 18-May-2010 Use one file / NSTROUT proc (NDISTIO(1)==1)
!   R. El Khatib : 01-Feb-2012 Extend I/O processors up to NPRTRV*NSTROUT
!   R. El Khatib 10-Aug-2011 NIOLEVG management
!   P. Marguinaud : 26-Apr-2012 : Refactor using IOSPECA_MOD
!   P. Marguinaud : 11-Sep-2012 : Refactor using WRGATHFLNM, IOMULTIBUF_MOD, DIWRSPEC_MOD, 
!                                 IOFLDDESC_MOD, MFIOOPTS_MOD, WRFLDCW_MOD + Compress data with OpenMP
!   P. Marguinaud : 10-Oct-2013 : Use FACTX, GATH_SPEC, EGATH_SPEC; remove field
!                                 based IO server mode
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   O. Marsden: June 2015 CY42 YRGMV, YRGFL, YRSURF, YRGMV5, and YRGFL5 are now passed by argument
!     ----------------------------------------------------------------------------

USE YOMRIP      , ONLY : TRIP
USE YOMXFU      , ONLY : TXFU
USE GEOMETRY_MOD, ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LELAM
USE YOMMP0   , ONLY : MYSETV, MYSETW, NSTROUT, NPRTRV, NPRTRW, NPROC, NPRINTLEV
USE YOMLUN       , ONLY : NULOUT
USE YOMOPH0  , ONLY : YMDLOPH
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE IOCPTDESC_MOD, ONLY : IOCPTDESC

USE IOSPECA_MOD, ONLY : IOSPECA_VSETOFF

USE FACTX_MOD, ONLY : FACTX, COMPACTFLD

!      -----------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY),     INTENT (IN)           :: YDGEOMETRY
TYPE(TRIP)          ,INTENT(INOUT)         :: YDRIP
TYPE(TXFU)          ,INTENT(INOUT)         :: YDXFU
INTEGER (KIND=JPIM), INTENT (IN)           :: KFIELDS
REAL (KIND=JPRB),    INTENT (IN)           :: PBUF (:,:)
TYPE (IOFLDDESC),    INTENT (IN)           :: YDFLDSC (:) 
TYPE (FACTX),        INTENT (INOUT)        :: YDFACTX
CHARACTER(LEN=*),    INTENT (IN)           :: CDFIC
INTEGER (KIND=JPIM), INTENT (IN), OPTIONAL :: KTAG

#include "wrgathflnm.intfb.h"
#include "sumpioh.intfb.h"
#include "set2pe.intfb.h"
#include "egath_spec.h"
#include "gath_spec.h"
#include "abor1.intfb.h"

!     ISTROUT: maximum number of I/O procs on a V-set
!     IFLDSPG: total number of fields to read
!     IFLDSCH: number of fields in each V-set
!     IVSETOFF: V-set offset for re-ordered fields
!     INFD   : number of fields in the current chunk for each I/O processor
!              (local variable with respect to the V-set)
!     IFLDOFF: fields offset after distribution among ISTROUT procs.
!              (local variable with respect to the V-set)
!     ZDATA  : a spectral field written to file
!     ZSPBUFG: a chunk of of spectral buffer fields to distribute (global fields)
!     ZVALCO : packed fields to write out
!     ISIZPK : dimension of a field for packing

INTEGER(KIND=JPIM) :: IFLDSPG, JSETV
INTEGER(KIND=JPIM) :: JSETW, ISIZPK, JFLDG1, JFLDG2
INTEGER(KIND=JPIM) :: ISTROUT, IPROC, INFL

INTEGER(KIND=JPIM), POINTER  :: INFD(:), IFLDOFF(:)                             ! Current V-set
INTEGER(KIND=JPIM)           :: INFD_PROC (NPROC), IFLDOFF_PROC (NPROC)
INTEGER(KIND=JPIM), TARGET   :: INFD_ALLVSETS (NPRTRW,NPRTRV), IFLDOFF_ALLVSETS (NPRTRW,NPRTRV)

INTEGER(KIND=JPIM) :: IFLDSCH (NPRTRV), IVSETOFF (NPRTRV)
TYPE (IOCPTDESC),    ALLOCATABLE :: YLCPDSC (:)
REAL(KIND=JPRB),     ALLOCATABLE :: ZSPBUFG(:,:)
REAL(KIND=JPRB),     ALLOCATABLE :: ZVALCO(:,:)

INTEGER(KIND=JPIM),  ALLOCATABLE :: ITO (:)

CHARACTER(LEN=LEN(CDFIC)) :: CLFIC(1)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRSP2FA',0,ZHOOK_HANDLE)

IFLDSPG = SIZE (YDFLDSC)

IF (IFLDSPG == 0) THEN
  IF (LHOOK) CALL DR_HOOK('WRSP2FA',1,ZHOOK_HANDLE)
  RETURN
ENDIF

!IF (NPRINTLEV > 0) THEN
!  DO JFLDG1 = 1, IFLDSPG
!    IF (YDFLDSC(JFLDG1)%IOLEV > 0) THEN
!      WRITE (NULOUT, '("WRSP2FA: WRITE OUT ",A,I3.3,A)') TRIM(YDFLDSC(JFLDG1)%CPREF),YDFLDSC(JFLDG1)%IOLEV,YDFLDSC(JFLDG1)%CSUFF
!    ELSE
!      WRITE (NULOUT, '("WRSP2FA: WRITE OUT ",A,A)') TRIM(YDFLDSC(JFLDG1)%CPREF),YDFLDSC(JFLDG1)%CSUFF
!    ENDIF
!  ENDDO
!ENDIF

ISTROUT=MIN(NPRTRW,NSTROUT)

CALL IOSPECA_VSETOFF (YDFLDSC%IVSET, IFLDSCH, IVSETOFF)

ISIZPK = MAXVAL (YDFLDSC%NSIZPK) ! Overdimension enables to write non-packed data

INFD_PROC = 0_JPIM
IFLDOFF_PROC = 0_JPIM

DO JSETV = 1, NPRTRV

  INFD    => INFD_ALLVSETS (:,JSETV)
  IFLDOFF => IFLDOFF_ALLVSETS (:,JSETV)

  CALL SUMPIOH(NPRTRW,ISTROUT,IFLDSCH(JSETV),INFD,IFLDOFF)

  DO JSETW = 1, NPRTRW

    IF (INFD_ALLVSETS (JSETW,JSETV) == 0) CYCLE

    CALL SET2PE (IPROC, 0, 0, JSETW, JSETV)

    INFD_PROC (IPROC)    = INFD_ALLVSETS (JSETW,JSETV)
    IFLDOFF_PROC (IPROC) = IFLDOFF_ALLVSETS (JSETW,JSETV)

  ENDDO

ENDDO

INFD    => INFD_ALLVSETS (:,MYSETV)
IFLDOFF => IFLDOFF_ALLVSETS (:,MYSETV)

IF (SUM (INFD) /= SIZE (PBUF, 2)) CALL ABOR1 ('WRSP2FA: DIMENSION MISMATCH')
  
INFL = INFD (MYSETW)

ALLOCATE (ZSPBUFG (YDGEOMETRY%YRDIM%NSPEC2G,INFL))

ALLOCATE (ITO (IFLDSPG))

ITO = 0
DO JSETV = 1, NPRTRV
  DO JSETW = 1, NPRTRW
    JFLDG1 = IVSETOFF (JSETV) + IFLDOFF_ALLVSETS (JSETW, JSETV) + 1
    JFLDG2 = IVSETOFF (JSETV) + IFLDOFF_ALLVSETS (JSETW, JSETV) + INFD_ALLVSETS (JSETW, JSETV)
    CALL SET2PE (IPROC, 0, 0, JSETW, JSETV)
    ITO (JFLDG1:JFLDG2) = IPROC
  ENDDO
ENDDO

IF (LELAM) THEN
  CALL EGATH_SPEC (KFGATHG=IFLDSPG, KTO=ITO, LDIM1_IS_FLD=.FALSE.,&
                 & PSPEC=PBUF, KVSET=YDFLDSC%IVSET, PSPECG=ZSPBUFG,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ELSE
  CALL GATH_SPEC (KFGATHG=IFLDSPG, KTO=ITO, LDIM1_IS_FLD=.FALSE.,&
                & PSPEC=PBUF, KVSET=YDFLDSC%IVSET, PSPECG=ZSPBUFG,KRESOL=YDGEOMETRY%YRDIM%NRESOL)
ENDIF

IF (YDFACTX%YIOOPTS%LDOCP) THEN

  ALLOCATE (ZVALCO (ISIZPK, INFL), YLCPDSC (INFL))

  JFLDG1=IVSETOFF(MYSETV)+IFLDOFF(MYSETW)+1
  JFLDG2=IVSETOFF(MYSETV)+IFLDOFF(MYSETW)+INFL

  CALL COMPACTFLD (YMDLOPH(1)%CFPCA, YDFACTX, YDFLDSC (JFLDG1:JFLDG2), ZSPBUFG, YLCPDSC, ZVALCO)

  CLFIC(1)=CDFIC
  CALL WRGATHFLNM (YDGEOMETRY%YVABIO,1_JPIM, INFL, ISIZPK, INFD_PROC, IFLDOFF_PROC, ZVALCO, &
                 & YLCPDSC, YMDLOPH, CLFIC, YDFACTX=YDFACTX, KTAG=KTAG, YDGEOMETRY=YDGEOMETRY)

  DEALLOCATE (ZVALCO, YLCPDSC)

ENDIF

DEALLOCATE(ITO)
   
DEALLOCATE(ZSPBUFG)

IF (LHOOK) CALL DR_HOOK('WRSP2FA',1,ZHOOK_HANDLE)

END SUBROUTINE WRSP2FA

