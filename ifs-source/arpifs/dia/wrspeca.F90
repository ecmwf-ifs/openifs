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

SUBROUTINE WRSPECA(YDGEOMETRY,YDGFL,YDSURF,YDXFU,YDML_GCONF,YDDYN,YDDYNA,YDSPEC,YDFACTX,CDFIC,YDMCUF)

!**** *WRSPECA*  - Write the spectral fields to FA (Arpege or Aladin)

!     Purpose.
!     --------
!         Write the spectral fields of the model to ARPEGE ALADIN file.

!**   Interface.
!     ----------
!        *CALL* *WRSPECA(.....)*

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
!   O. Marsden: June 2015 CY42 YRGMV, YRGFL, and YRSURF are now passed by argument
!   O. Marsden: Sept 2016 Added spectral argument for call to WRSPECA_GP
!     ----------------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYN             , ONLY : TDYN
USE YOMDYNA            , ONLY : TDYNA
USE YOMXFU             , ONLY : TXFU
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL   , ONLY : TGFL
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT0   , ONLY : LWRSPECA_GP
USE YOMMP0   , ONLY : MYSETV, NPRTRV
USE YOMTAG, ONLY : MTAG_MFIO_WRSPECA
USE SPECTRAL_FIELDS_MOD
USE YOMMCUF                , ONLY : TMCUF
USE IOSPECA_MOD, ONLY : IOSPECA_SELECTD,&
                      & IOSPECA_SELECTF,&
                      & IOSPECA_START,&
                      & IOSPECA_COUNT,&
                      & IOSPECA_FINISH,&
                      & IOSPECA_VSETOFF,&
                      & IOSPECA_CTX,&
                      & NSPECACT_WRITE
USE IOFLDDESC_MOD, ONLY : IOFLDDESC

USE FACTX_MOD, ONLY : FACTX

USE YOMFP_SERV, ONLY : FP_SERV_C001

!      -----------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)               , INTENT(IN)    :: YDGEOMETRY
TYPE(TGFL)                   , INTENT(INOUT) :: YDGFL
TYPE(TSURF)                  , INTENT(INOUT) :: YDSURF
TYPE(TDYN)                   , INTENT(INOUT) :: YDDYN
TYPE(TDYNA)                  , INTENT(INOUT) :: YDDYNA
TYPE(MODEL_GENERAL_CONF_TYPE), INTENT(INOUT) :: YDML_GCONF
TYPE(TXFU)                   , INTENT(INOUT) :: YDXFU
TYPE(SPECTRAL_FIELD)         , INTENT(INOUT) :: YDSPEC
TYPE (FACTX)                 , INTENT(INOUT) :: YDFACTX ! logical unit for file
CHARACTER(LEN=*), INTENT(IN) :: CDFIC
TYPE(TMCUF)                  , INTENT(INOUT), OPTIONAL :: YDMCUF

#include "wrspeca_gp.intfb.h"
#include "wrspeca_fp.intfb.h"
#include "wrsp2fa.intfb.h"

!      -----------------------------------------------------------

!     IFLDSPG: total number of fields to read
!     IFLDSCH: number of fields in each V-set
!     IVSETOFF: V-set offset for re-ordered fields
!     ZDATA  : a spectral field written to file
!     ZSPBUFG: a chunk of of spectral buffer fields to distribute (global fields)
!     ZSPBUFL: local fields
!     ZVALCO : packed fields to write out

INTEGER(KIND=JPIM) :: IFLDSPG, IFLDSPL

INTEGER(KIND=JPIM) :: IFLDSCH (NPRTRV), IVSETOFF (NPRTRV)
TYPE (IOFLDDESC),    ALLOCATABLE :: YLFLDSC (:)
REAL(KIND=JPRB),     ALLOCATABLE :: ZSPBUFL(:,:)

TYPE (IOSPECA_CTX) :: YLIOSP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('WRSPECA',0,ZHOOK_HANDLE)

IF (FP_SERV_C001%LFP_CLIENT) THEN
  CALL WRSPECA_FP(YDML_GCONF%YRRIP,FP_SERV_C001,YDGEOMETRY)
  IF (.NOT. FP_SERV_C001%LFP_CLIENT_WRITE) THEN
    IF (LHOOK) CALL DR_HOOK('WRSPECA',1,ZHOOK_HANDLE)
    RETURN
  ENDIF
ENDIF

IF (LWRSPECA_GP) THEN
  CALL WRSPECA_GP(YDGEOMETRY,YDXFU,YDML_GCONF,YDDYN,YDDYNA,YDSPEC,YDGFL,YDSURF,YDFACTX,CDFIC, YDMCUF=YDMCUF)
  IF (LHOOK) CALL DR_HOOK('WRSPECA',1,ZHOOK_HANDLE)
  RETURN
ENDIF

CALL IOSPECA_COUNT(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_WRITE,IFLDSPG, YDMCUF=YDMCUF)

CALL IOSPECA_START(YDGEOMETRY,YDGFL,YDML_GCONF,YDDYNA,NSPECACT_WRITE,YLIOSP,YDSPEC)

ALLOCATE (YLFLDSC(IFLDSPG))

CALL IOSPECA_SELECTD(YDGEOMETRY,YDML_GCONF,YDDYNA,NSPECACT_WRITE,YLFLDSC,YDMCUF=YDMCUF)

CALL IOSPECA_VSETOFF (YLFLDSC%IVSET, IFLDSCH, IVSETOFF)

IFLDSPL = COUNT (YLFLDSC%IVSET == MYSETV)

ALLOCATE (ZSPBUFL (YDGEOMETRY%YRDIM%NSPEC2, IFLDSPL))

CALL IOSPECA_SELECTF(YDGEOMETRY,YDML_GCONF%YGFL,YDDYN, NSPECACT_WRITE, YLIOSP, IFLDSPL,&
                    & IVSETOFF (MYSETV), 0_JPIM, YLFLDSC, ZSPBUFL, YDSPEC=YDSPEC, YDMCUF=YDMCUF)

CALL WRSP2FA(YDGEOMETRY,YDXFU,YDML_GCONF%YRRIP,IFLDSPG,ZSPBUFL,YLFLDSC,YDFACTX,CDFIC,KTAG=MTAG_MFIO_WRSPECA)

DEALLOCATE (ZSPBUFL)
   
DEALLOCATE (YLFLDSC)

CALL IOSPECA_FINISH (YDGEOMETRY, YDDYNA,NSPECACT_WRITE, YLIOSP, YDSPEC=YDSPEC)

IF (LHOOK) CALL DR_HOOK('WRSPECA',1,ZHOOK_HANDLE)

END SUBROUTINE WRSPECA

