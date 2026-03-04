! (C) Copyright 1989- Meteo-France.

SUBROUTINE WRMLFP(YDRQSP,YDRQGP,KFPRESOL,KSPEC2,KSPEC2G,KFPRGPL,KFPRGPG, &
 & KFPROMA,YDGEOMETRY,YDSURF,YDECLDP,YDEPHY,YDRQPHY,CDCONF, &
 & LDLAG,LDOCE,YDOFN,KSTEP,PFPBUF,PTSTEP,PSPEC)

!**** *WRMLFP*  - writes out the model level and sfc fields in GRIB

!     Purpose.
!     --------
!     Write out the model level and sfc fields in GRIB

!**   Interface.
!     ----------
!        *CALL* *WRMLFP(...)

!        Explicit arguments :    
!        --------------------
!           KFPRESOL : pp resolution tag
!           KSPEC2   : pp local number of waves
!           KSPEC2G  : pp global number of waves
!           KFPRGPL  : pp local number of gridpoints
!           KFPROMA  : pp blocking factor
!           CDCONF - configuration of call
!           LDLAG - Post-processing physics fields if called in lagged mode
!           PTSTEP : model time step
!           LDOCE   : Post-processing ocean fields if called in ocean mode

!        Implicit arguments :      The state variables of the model
!        --------------------

!     Method.
!     -------
!        See documentation
!        - spectral part would only work if input-resol = output resol
!          for vertical and horizontal resolution
!        - surface fields pp. without any fluxes ...

!     Externals.   MODMSEC - modify GRIB headers, output resolution
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Nils Wedi *ECMWF*
!      Original : 97-02-20 adapted from WRMLPPG,WRHFPG

!     Modifications.
!     --------------

!      R. El Khatib : 01-04-10 Enable post-processing of filtered spectra
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      D.Dent       (02-05-15) : remove reference to extra field 099
!      R. El Khatib : 03-04-17 Fullpos improvments
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      A. Tompkins : 11-02-04: adding TCLW, TCIW
!      P. Viterbo  : 24-05-2004   Change surface units
!      S. Serrar   : 26-09-2005 total column CO2 and SF6 added for postprocessing
!      M.Hamrud      01-May-2006 Generalized IO scheme
!      R. Engelen  : 02-Jun-06  : CO2 and SF6 replaced by generic GHG
!      JJMorcrette : 20060630   : MODIS Albedo
!      S. Serrar   : 07-Sep-06  : post-process in lagged mode total columns for a few tracers
!      JJMorcrette  20060925    DU, BC, OM, SU, VOL, SOA climatological fields 
!      G. Balsamo  : 20070115   : Soil Type (SOTY)
!      P. Bechtold : 20070829   : CAPE
!      H. Hersbach : 04-Dec-2009: 10m-neutral wind and friction velocity
!      P. Bechtold : 20110809   : CIN, Convective Indices
!      R. El Khatib: 29-Feb-2012 simplified interface to norms computation
!      R. El Khatib: 10-Aug-2012 prepare for Fullpos-2
!      R. El Khatib 20-Aug-2012 GAUXBUF removed and replaced by HFPBUF
!      A. Inness   : 23-Mar-2012: Add total column CHEM fields
!      R. El Khatib 13-Dec-2012 Fullpos buffers reshaping
!      R. Forbes    01-Mar-2014 Add precip rates/type,TCSLW,I10FG
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      JJMorcrette   20130730 15-variable aerosol model + oceanic DMS
!      R. El Khatib 04-Aug-2014 Pruning of the conf. 927/928
!      R. Forbes    10-Jan-2015 Add freezing rain FZRA
!     ------------------------------------------------------------------

USE YOECLDP            , ONLY : TECLDP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMMP0             , ONLY : LSPLITOUT, NPROC, MYSETV
USE YOMCT0   , ONLY : NCONF    ,LOBSC1
USE YOEPHY             , ONLY : TEPHY
USE YOMVAR             , ONLY : LTWBGV   ,LTWCGL  
USE YOMGRIB            , ONLY : JPNUMLPP,NSTEPLPP
USE YOMFP4L, ONLY : TRQFP
USE TYPE_FPOFN, ONLY : TFPOFN
USE MPL_MODULE         , ONLY : MPL_BARRIER,MPL_ALLREDUCE
USE IOSPECE_MOD        , ONLY : IOSPECE_ML_SELECTF,  IOSPECE_ML_SELECTD,  IOSPECE_ML_COUNT
USE IOGRIDE_MOD        , ONLY : IOGRIDE_SELECTF,     IOGRIDE_SELECTD,     IOGRIDE_COUNT
USE IOGRIDUE_MOD       , ONLY : IOGRIDUE_ML_SELECTF, IOGRIDUE_ML_SELECTD, IOGRIDUE_ML_COUNT
USE IOGRIDOE_MOD       , ONLY : IOGRIDOE_SELECTF,    IOGRIDOE_SELECTD,    IOGRIDOE_COUNT
USE IOFLDDESC_MOD      , ONLY : IOFLDDESC
USE YOMIO_SERV         , ONLY : IO_SERV_C001

IMPLICIT NONE

TYPE (TRQFP),  INTENT(IN) :: YDRQSP
TYPE (TRQFP),  INTENT(IN) :: YDRQGP
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPRESOL
INTEGER(KIND=JPIM),INTENT(IN)  :: KSPEC2
INTEGER(KIND=JPIM),INTENT(IN)  :: KSPEC2G
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPRGPL
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPRGPG
INTEGER(KIND=JPIM),INTENT(IN)  :: KFPROMA
TYPE(GEOMETRY)  ,INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)     ,INTENT(IN)    :: YDSURF
TYPE(TECLDP)    ,INTENT(IN) :: YDECLDP
TYPE(TEPHY)     ,INTENT(IN) :: YDEPHY
TYPE (TRQFP),  INTENT(IN) :: YDRQPHY
CHARACTER(LEN=1),INTENT(IN)    :: CDCONF 
LOGICAL         ,INTENT(IN)    :: LDLAG 
LOGICAL         ,INTENT(IN)    :: LDOCE
TYPE (TFPOFN)   ,INTENT(IN)    :: YDOFN
INTEGER(KIND=JPIM),INTENT(IN)  :: KSTEP
REAL(KIND=JPRB) ,INTENT(IN)    :: PFPBUF(:,:,:)
REAL(KIND=JPRB) ,INTENT(IN), OPTIONAL :: PTSTEP
REAL(KIND=JPRB) ,INTENT(IN), OPTIONAL :: PSPEC(:,:)

INTEGER(KIND=JPIM), EXTERNAL :: ISRCHEQ
CHARACTER (LEN=1) ::  CLMODE

REAL(KIND=JPRB),    ALLOCATABLE :: ZREAL(:,:)
REAL(KIND=JPRB),    ALLOCATABLE :: ZPTRFLD(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IGRIBIOSH(:,:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IGRIBIOGG(:,:)
TYPE (IOFLDDESC),   ALLOCATABLE :: YLFLDSC (:)
INTEGER (KIND=JPIM) :: IFNUM_SP, IFNUM_GRE, IFNUM_GRUE, IFNUM_GROE,JF,ILST

LOGICAL :: LLDYN, LLPHY,LLUPD_PPT

#include "fcttim.func.h"

#include "wroutgpgb.intfb.h"
#include "wroutspgb.intfb.h"
#include "wrmlfp_io_serv.intfb.h"

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!*       1.    PREPARATIONS.
!              --------------

IF (LHOOK) CALL DR_HOOK('WRMLFP',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFIELDL=>YDRQSP%NFIELDL, NGPTOT=>YDGEM%NGPTOT)
IF (IO_SERV_C001%NPROC_IO > 0) THEN
  CALL WRMLFP_IO_SERV(YDRQSP,YDRQGP,YDGEOMETRY,YDSURF,YDECLDP,YDEPHY,YDRQPHY,KFPRGPL,KFPRGPG,KFPROMA,KSPEC2,KSPEC2G, &
   & CDCONF,LDLAG,LDOCE,PFPBUF,PTSTEP=PTSTEP,PSPEC=PSPEC)
ENDIF

!     *  COMPUTE FIELD POINTERS

LLDYN=CDCONF /= 'I'
LLPHY=CDCONF == 'I'

!      -----------------------------------------------------------

!*       2.    WRITE SPECTRAL DATA ON MODEL LEVELS.
!              ------------------------------------

!*       2.1   EXTRACT AND COPY SPECTRAL FIELDS
!              --------------------------------

IF (IO_SERV_C001%NPROC_IO == 0) THEN
  IF (YDRQSP%NFIELDG > 0.AND.LLDYN .AND. .NOT. (LDLAG.OR.LDOCE)) THEN

    IFNUM_SP = 0
    CALL IOSPECE_ML_COUNT (YDRQSP%NFIELDG,IFNUM_SP)

    ALLOCATE (YLFLDSC (IFNUM_SP))
    CALL IOSPECE_ML_SELECTD (YDRQSP,KSPEC2G,YLFLDSC)

    ALLOCATE (ZPTRFLD (KSPEC2, NFIELDL))
    CALL IOSPECE_ML_SELECTF (PSPEC, ZPTRFLD, PACK (YLFLDSC, MASK=YLFLDSC%IVSET == MYSETV))

! Initialise igribiosh data structure
!   1 - grib code
!   2 - level 
!   3 - bset to send this field
!   4 - local level for sending bset

    ALLOCATE (IGRIBIOSH (4, IFNUM_SP))
    IGRIBIOSH (1,1:IFNUM_SP) = YLFLDSC%IGRIB
    IGRIBIOSH (2,1:IFNUM_SP) = YLFLDSC%ILEVG
    IGRIBIOSH (3,1:IFNUM_SP) = YLFLDSC%IVSET
    IGRIBIOSH (4,1:IFNUM_SP) = YLFLDSC%ILEVL
    DEALLOCATE (YLFLDSC)

!  Gather global fields and write out

    CALL WROUTSPGB(YDGEOMETRY%YRDIM,ZPTRFLD,IFNUM_SP,IGRIBIOSH,'m',KFPRESOL,YDOFN%CSH,PTSTEP)

    DEALLOCATE (ZPTRFLD)
    DEALLOCATE (IGRIBIOSH)

  ENDIF
ENDIF
!      -----------------------------------------------------------
!*       4.0   WRITE OUT SURFACE GRID POINT FIELDS
!              -----------------------------------

IF ( LLPHY.AND.(.NOT.(LTWBGV.OR.LTWCGL).AND.((&
   & YDEPHY%LEPHYS).OR.LOBSC1.OR.NCONF == 131)) ) THEN  

  IFNUM_GRE = 0
  CALL IOGRIDE_COUNT(YDGEOMETRY%YRGEM,YDEPHY,YDRQPHY,KFPRGPL,KFPROMA,IFNUM_GRE)

  IF(IFNUM_GRE > 0) THEN

    ALLOCATE (YLFLDSC (IFNUM_GRE))
    CALL IOGRIDE_SELECTD(YDGEM,YDEPHY,YDRQPHY,KFPRGPL,KFPROMA,YLFLDSC)
   
    ALLOCATE (ZREAL (KFPRGPL, IFNUM_GRE))
    CALL IOGRIDE_SELECTF(YDGEM,YDEPHY,YDRQPHY,KFPRGPL,KFPROMA,ZREAL,PFPBUF,YLFLDSC)
   
    ALLOCATE (IGRIBIOGG (2, IFNUM_GRE))
    IGRIBIOGG (1,1:IFNUM_GRE) = YLFLDSC%IGRIB
    IGRIBIOGG (2,1:IFNUM_GRE) = YLFLDSC%ILEVG

    IF (IO_SERV_C001%NPROC_IO == 0) THEN 
      IF(LDLAG.OR.LDOCE) THEN
        CLMODE='a' 
      ELSE
        CLMODE='w'
      ENDIF
      CALL WROUTGPGB(ZREAL,KFPRGPL,IFNUM_GRE,IGRIBIOGG,'s',CLMODE,KFPRESOL,YDOFN%CGG,PTSTEP)
    ENDIF
    LLUPD_PPT = .FALSE.
    DO JF=1,IFNUM_GRE
      ILST = ISRCHEQ(JPNUMLPP,NSTEPLPP(:,1),1,YLFLDSC(JF)%IGRIB)
      IF(ILST <= JPNUMLPP) THEN
        NSTEPLPP(ILST,2) = KSTEP
        LLUPD_PPT = .TRUE.
      ENDIF
    ENDDO
    IF(LLUPD_PPT) CALL MPL_ALLREDUCE(NSTEPLPP(:,2),'MAX')
    DEALLOCATE (YLFLDSC)
    DEALLOCATE (IGRIBIOGG)
    DEALLOCATE (ZREAL)

  ENDIF

ENDIF

!      -----------------------------------------------------------

!*       5.0   WRITE OUT UPPER AIR GRID POINT FIELDS
!              --------------------------------------
IF (IO_SERV_C001%NPROC_IO == 0) THEN
  IF( LSPLITOUT )THEN
    CLMODE = 'w'
  ELSE
    CLMODE = 'a'
  ENDIF

  IF (LLDYN.AND.(YDRQGP%NFIELDG) > 0 .AND. .NOT. (LDLAG.OR.LDOCE)) THEN

    IFNUM_GRUE = 0
    CALL IOGRIDUE_ML_COUNT(YDRQGP,YDGEOMETRY%YRGEM,YDECLDP,KFPRGPL,KFPROMA,IFNUM_GRUE)

    IF(IFNUM_GRUE > 0) THEN

      ALLOCATE (YLFLDSC (IFNUM_GRUE))
      CALL IOGRIDUE_ML_SELECTD(YDRQGP,YDGEOMETRY%YRGEM,YDECLDP,KFPRGPL,KFPROMA,YLFLDSC)
      
      ALLOCATE (ZREAL (KFPRGPL, IFNUM_GRUE))
      CALL IOGRIDUE_ML_SELECTF(YDRQGP,YDGEM,YDECLDP,KFPRGPL,KFPROMA,ZREAL,PFPBUF,YLFLDSC)
      
      ALLOCATE (IGRIBIOGG (2, IFNUM_GRUE))
      IGRIBIOGG (1,1:IFNUM_GRUE) = YLFLDSC%IGRIB
      IGRIBIOGG (2,1:IFNUM_GRUE) = YLFLDSC%ILEVG
      DEALLOCATE (YLFLDSC)

      CALL WROUTGPGB(ZREAL,KFPRGPL,IFNUM_GRUE,IGRIBIOGG,'m',CLMODE,KFPRESOL,YDOFN%CUA,PTSTEP)
      DEALLOCATE (IGRIBIOGG)
      DEALLOCATE (ZREAL)

    ENDIF

  ENDIF

!*       6.0   WRITE OUT OCEAN MIXED LAYER FIELDS
!              ----------------------------------

!********* this is NOT Fullpos. REK

  IF(  .NOT. (LDLAG.OR.LDOCE) ) THEN
    IFNUM_GROE = 0
    CALL IOGRIDOE_COUNT(YDGEOMETRY,YDSURF,IFNUM_GROE)
    IF(IFNUM_GROE > 0) THEN
      
      ALLOCATE (YLFLDSC (IFNUM_GROE))
      CALL IOGRIDOE_SELECTD(YDGEOMETRY,YDSURF,YLFLDSC)
    
      ALLOCATE (ZREAL (NGPTOT,IFNUM_GROE))
      CALL IOGRIDOE_SELECTF(YDGEOMETRY,YDSURF,ZREAL, YLFLDSC)
      
      ALLOCATE (IGRIBIOGG (2, IFNUM_GROE))
      IGRIBIOGG (1,1:IFNUM_GROE) = YLFLDSC%IGRIB
      IGRIBIOGG (2,1:IFNUM_GROE) = YLFLDSC%ILEVG
      DEALLOCATE (YLFLDSC)
    
      CALL WROUTGPGB(ZREAL,NGPTOT,IFNUM_GROE,IGRIBIOGG,'m',CLMODE,KFPRESOL,YDOFN%CGG,PTSTEP)
      DEALLOCATE (IGRIBIOGG)
      DEALLOCATE (ZREAL)
    
    ENDIF
  ENDIF
ENDIF

IF( NPROC > 1 )THEN
  CALL GSTATS(792,0)
  CALL MPL_BARRIER(CDSTRING='WRMLFP:')
  CALL GSTATS(792,1)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRMLFP',1,ZHOOK_HANDLE)

END SUBROUTINE WRMLFP

