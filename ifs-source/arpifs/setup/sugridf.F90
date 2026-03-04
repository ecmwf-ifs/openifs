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

SUBROUTINE SUGRIDF(YDGEOMETRY,YDSURF,YDMODEL,KFILE,PSPOR,LDSPOR,PBLH,LDCLI,LDRESET,KINITMONTH)

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF, GPOPER
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : N3DINI, LCANARI
USE YOMARG             , ONLY : NGRIBFILE
USE YOMDYNCORE         , ONLY : LAQUA

!**** *SUGRIDF*  - Routine to initialize the gridpoint fields

!     Purpose.
!     --------
!           Initialize the gridpoint fields of the model.

!**   Interface.
!     ----------
!        *CALL* *SUGRIDF(...)*

!        Explicit arguments :
!        --------------------
!        KFILE : an indicator for which file is to be read
!               KFILE = 0 the CFNISH file is read
!               KFILE = 1 the CFNGSH file is read
!               KFILE = 2 the CFNRF  file is read
!               KFILE = 4 the CFNGG  file is read
!               KFILE = 6 a climat. file is read
!               KFILE = 8 the CFANS  file is read
!               KFILE =12 the CFNBGHRGG file is read
!        PSPOR
!        LDSPOR
!        PBLH
!        LDCLI : Monthly climatology required
!        LDRESET : (re)set fields to default values at start
!        KINITMONTH : month of the initial file ; used to control the climatology, if provided in the interface

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
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 91-11-15

!     Modifications.
!     --------------
!      02-03-29  Y. Bouteloup : Climatological ozone profiles
!      03-10-01  M.Hamrud : CY28 Cleaning
!      02-07-17  P. Marquet : add VCLIA for aerosol files
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      D.Salmond     22-Nov-2005 Mods for coarser/finer physics
!      06-01-06 A. Alias Sulfate and volcano aerosols added
!      M.Hamrud      01-Jul-2006 Revised surface fields
!      R. El Khatib  26-Nov-2008  Cleaning
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TGSGEOM
!      G.Balsamo     13-Jun-2013 Add lake model fields
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      R. El Khatib 25-Mar-2016 LDCLI, LDRESET
!      F. Taillefer April 2016  Add prints
!      R. El Khatib 09-Sep-2016 better interoperability GRIB2 vs FA
!      E. Dutra/Arduini  Jun 2017: snow multi-layer changes in SNOWG
!      M. Chrust    06-Apr-2021 Add high resolution background=12
!     ----------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)       ,INTENT(INOUT) :: YDSURF
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFILE
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSPOR(:)
LOGICAL           ,INTENT(OUT)   :: LDSPOR
REAL(KIND=JPRB)   ,OPTIONAL,INTENT(OUT)   :: PBLH(YDGEOMETRY%YRGEM%NGPTOT)
LOGICAL           ,OPTIONAL,INTENT(IN)   :: LDCLI
LOGICAL           ,OPTIONAL,INTENT(IN)   :: LDRESET
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN)   :: KINITMONTH
!     ----------------------------------------------------------------------
INTEGER(KIND=JPIM) :: JV, ICLI
INTEGER(KIND=JPIM) :: JBL

LOGICAL :: LLINT, LLCLI, LLRESET
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ----------------------------------------------------------------------

#include "abor1.intfb.h"
#include "gpnorm2.intfb.h"
#include "gpnorm3.intfb.h"
#include "sugrclia.intfb.h"
#include "sugrida.intfb.h"
#include "sugridg.intfb.h"
#include "sugridva.intfb.h"
#include "gp_sstaqua.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGRIDF',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB)
ASSOCIATE(NGPBLKS=>YDDIM%NGPBLKS, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & SD_VA=>YDSURF%SD_VA, SD_VC=>YDSURF%SD_VC, SD_VF=>YDSURF%SD_VF, &
 & SD_VV=>YDSURF%SD_VV, SD_WS=>YDSURF%SD_WS, SD_WW=>YDSURF%SD_WW, &
 & SP_RR=>YDSURF%SP_RR, SP_SB=>YDSURF%SP_SB, SP_SG=>YDSURF%SP_SG, &
 & SP_SL=>YDSURF%SP_SL, SP_CI=>YDSURF%SP_CI, SP_CL=>YDSURF%SP_CL, &
 & SP_X2=>YDSURF%SP_X2, SD_VX=>YDSURF%SD_VX, &
 & YSD_VAD=>YDSURF%YSD_VAD, YSD_VCD=>YDSURF%YSD_VCD, &
 & YSD_VFD=>YDSURF%YSD_VFD, YSD_VVD=>YDSURF%YSD_VVD, YSD_VXD=>YDSURF%YSD_VXD, &
 & YSD_WSD=>YDSURF%YSD_WSD, YSD_WWD=>YDSURF%YSD_WWD, YSP_RRD=>YDSURF%YSP_RRD, &
 & YSP_SBD=>YDSURF%YSP_SBD, YSP_SGD=>YDSURF%YSP_SGD, YSP_SLD=>YDSURF%YSP_SLD, &
 & YSP_CID=>YDSURF%YSP_CID, YSP_CLD=>YDSURF%YSP_CLD, YSP_X2D=>YDSURF%YSP_X2D )
!     ------------------------------------------------------------------
!*       1.    SET DEFAULT VALUES FOR FIELDS.
!              ------------------------------

CALL GSTATS(26,0)

IF (PRESENT(LDRESET)) THEN
  LLRESET=LDRESET
ELSE
  LLRESET=.TRUE.
ENDIF
IF (LLRESET) THEN
  DO JBL=1,NGPBLKS
    CALL GPOPER(YDGEOMETRY%YRDIM,YDMODEL%YRML_DYN%YRDYN,'SETDEFAULT',YDSURF,KBL=JBL)
  ENDDO
ENDIF

!     ------------------------------------------------------------------
!*       3.    READ IN FIELDS FROM ARPEGE FILE.
!              -------------------------------

IF ((NGRIBFILE /= 1).AND.(N3DINI == 0)) THEN
  IF (KFILE /= 6) THEN
    CALL SUGRIDA(YDGEOMETRY,YDSURF,YDMODEL,KFILE)
  ELSE
    CALL SUGRIDVA(YDGEOMETRY,YDSURF,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC)
  ENDIF
ENDIF

! WARNING : THIS GROUP IS RESERVED FOR ARPEGE, AND NOT ALLOWED TO IFS.
! THEREFORE NGRIBFILE IS NOT TESTED. REK
IF (PRESENT(LDCLI)) THEN
  LLCLI=LDCLI.AND.(YSD_VXD%NUMFLDS > 0)
ELSE
  LLCLI=(YSD_VXD%NUMFLDS > 0)
ENDIF
IF (LLCLI) THEN
  LLINT=LCANARI
  CALL SUGRCLIA(YDGEOMETRY,YDSURF,YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_GCONF%YRRIP,YDMODEL%YRML_LBC,LLINT,ICLI,KINITMONTH)
  IF ((LLINT .AND. ICLI /=2).AND..NOT.LCANARI) THEN
    WRITE(NULOUT,'('' LLINT = '',L2, '' ICLI = '',I2,'' LCANARI = '',L2)')  LLINT,ICLI,LCANARI
    CALL ABOR1('SUGRIDF : SUGRCLIA FAILED')
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*       4.    READ IN FIELDS FROM GRIB FILE.
!              -----------------------------

IF ((NGRIBFILE == 1).AND.(N3DINI == 0)) THEN
  IF (PRESENT(PBLH)) THEN
    CALL SUGRIDG(YDGEOMETRY,YDSURF,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_DYN%YRDYN,KFILE,LDSPOR,PBLH=PBLH)
  ELSE
    CALL SUGRIDG(YDGEOMETRY,YDSURF,YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_AOC%YRMCC,YDMODEL%YRML_DYN%YRDYN,KFILE,LDSPOR)
  ENDIF
ENDIF


IF (LAQUA) THEN
  CALL GP_SSTAQUA(YDGEOMETRY,YDSURF,YDGSGEOM_NB%GEMU,YDGSGEOM_NB%GELAM,SP_SB,SP_SG,SP_RR,SD_VF,SP_SL)
ENDIF

!*       7.    GRID-POINT NORMS
!              ----------------

WRITE(NULOUT,'(A)') ' SUGRIDF: STATISTICS FOR ALL SURFACE FIELDS'
IF (YSP_SBD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A,I1,A)') ' SOILB   ',YSP_SBD%NUMFLDS,' FIELDS '
  DO JV=1,YSP_SBD%NUMFLDS
    CALL GPNORM3(YDGEOMETRY,SP_SB,YSP_SBD%NUMFLDS,JV,KLEVS=YSP_SBD%NLEVS,LDLEVELS=.TRUE.)
  ENDDO
ENDIF
IF (YSP_SGD%NUMFLDS>0) THEN
  DO JV=1,YSP_SGD%NUMFLDS
    WRITE(NULOUT,'(A,I4)') ' SNOWG, FIELD ',JV
    CALL GPNORM3(YDGEOMETRY,SP_SG,YSP_SGD%NUMFLDS,JV,KLEVS=YSP_SGD%NLEVS,LDLEVELS=.TRUE.)
  ENDDO
ENDIF
IF (YSP_CID%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' CANRI '
  CALL GPNORM2(YDGEOMETRY,YSP_CID%NUMFLDS,1,SP_CI)
ENDIF
IF (YSP_CLD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' CLS '
  CALL GPNORM2(YDGEOMETRY,YSP_CLD%NUMFLDS,1,SP_CL)
ENDIF
IF (YSP_X2D%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' XTRP2 '
  CALL GPNORM2(YDGEOMETRY,YSP_X2D%NUMFLDS,1,SP_X2)
ENDIF
IF (YSP_SLD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' LAKEB '
  CALL GPNORM2(YDGEOMETRY,YSP_SLD%NUMFLDS,1,SP_SL)
ENDIF
IF (YSP_RRD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' RESVR '
  CALL GPNORM2(YDGEOMETRY,YSP_RRD%NUMFLDS,1,SP_RR)
ENDIF
IF (YSD_WSD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' FROM WAM '
  CALL GPNORM2(YDGEOMETRY,YSD_WSD%NUMFLDS,1,SD_WS)
ENDIF
IF (YSD_WWD%NUMFLDS>0) THEN
  WRITE(NULOUT,'(A)') ' TO WAM '
  CALL GPNORM2(YDGEOMETRY,YSD_WWD%NUMFLDS,1,SD_WW)
ENDIF


IF (YSD_VFD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VARSF'
  CALL GPNORM2(YDGEOMETRY,YSD_VFD%NUMFLDS,1,SD_VF)
ENDIF
IF (YSD_VVD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIV'
  CALL GPNORM2(YDGEOMETRY,YSD_VVD%NUMFLDS,1,SD_VV)
ENDIF
IF (YSD_VAD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIA'
  CALL GPNORM2(YDGEOMETRY,YSD_VAD%NUMFLDS,1,SD_VA)
ENDIF
IF (YSD_VCD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VO3ABC'
  CALL GPNORM2(YDGEOMETRY,YSD_VCD%NUMFLDS,1,SD_VC)
ENDIF
IF (YSD_VXD%NUMFLDS > 0) THEN
  WRITE(NULOUT,'(A)') ' VCLIX'
  CALL GPNORM2(YDGEOMETRY,YSD_VXD%NUMFLDS,1,SD_VX)
ENDIF

CALL GSTATS(26,1)

!     ------------------------------------------------------------------
END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUGRIDF',1,ZHOOK_HANDLE)
END SUBROUTINE SUGRIDF


