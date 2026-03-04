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

SUBROUTINE WRMLPPA(YDGEOMETRY,YDGFL,YDSURF,YDSPEC,YDCFU,YDXFU,YDMODEL,CDCONF,YDMCUF)

!     Purpose.
!     --------
!     Write out the model level fields to ARPEGE file

!**   Interface.
!     ----------
!        *CALL* *WRMLPPA(CDCONF)

!        Explicit arguments :     CDCONF - configuration of call
!        --------------------

!        Implicit arguments :      The state variables of the model
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Note de travail ARPEGE NR 17

!     Author.
!     -------
!      MPP Group *METEO-FRANCE*
!      Original 96-02-05 (From WRMLPP)

!     Modifications.
!     --------------
!   Modified 01-03-19 S. Ivatek-Sahdan call DIWRSPE -> call DIWRSPE0 
!   Modified 01-03-21 S. Ivatek-Sahdan Add call WRGRIDUA
!   Modified 01-03-22 S. Ivatek-Sahdan bf of mergeing
!   Modified 01-06-25 J. Masek - rescaling of new NH variables
!   R. El Khatib : 01-08-07 Pruning options + wrspeca
!   Modified 01-10-17 P. Marquet: INUMG instead of IFIELD in ZREALG
!   Modified 02-03-08 C. Fischer: rescaling of NH var. moved to wrspeca
!   Modified 02-09-30 P. Smolikova : interface to WRSPECA for d4 in NH
!   O.Spaniel    : 03-04-15 cleaning-added interface wrspeca.h
!   R. El Khatib : 03-04-17 Cleanups
!   K. YESSAD    : 03-07-04 Bugfix for lsprt=true.
!   R. El Khatib : 03-08-08 gfl+get rid of wrspeca dummies since it is called one here
!   R. El Khatib : 03-08-18 Roughness lengths not packed
!   G. Hello : 03-09-19 rescaling of NH var
!   M.Hamrud      01-Oct-2003 CY28 Cleaning
!   Modified 02-07-17 P. Marquet add VCLIA for aerosol files
!   M.Hamrud      01-Jul-2006 Revised surface fields
!   N. Wedi and K. Yessad (Jan 2008): different dev for NH model and PC scheme
!   Y.Seity        11-Jan-2008 Add CDCONF(1:1)='S' for surfex ouput files
!   R. El Khatib : 15-Sep-2008 I/O savings
!   A.Alias      : 07-Aug-2009 Add CDCONF(1:1)='Q' when surfex and hist. ouput files
!   R. El Khatib : 18-Mar-2010 WRGRIDA + barrier
!   P.Marguinaud : 18-May-2010 Use one file / NSTROUT proc (NDISTIO(1)==1)
!   P.Marguinaud : 14-Jan-2011 Use WRGRIDALL & IO server
!   R. El Khatib : 09-Aug-2011 LWRSPEC
!   R. El Khatib : 01-Feb-2012 LLISIO is now an argument of inifaout
!   D. Degrauwe  (Feb 2012): LARPEGEF_WRGP_HIST added
!   P.Marguinaud : Call WRSFX (SURFEX fields) + cleaning
!   P.Marguinaud : 26-Apr-2012 : Temperature & other parameters massage only for 
!                                WRGPA (no need to do that here for WRSPECA
!                                anymore)
!   P.Marguinaud : 11-Sep-2012 : Move all transforms into WRGPA + duplicate FA library :
!                                each thread gets its own copy (this will allow OpenMP on
!                                data compression)
!   P.Marguinaud : 15-May-2013 : Remove LARPEGEF_WRGP_HIST
!   P.Marguinaud : 10-Oct-2013 : Use FACTX
!   T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!   O. Marsden: June 2015 CY42 YRGMV, YRGFL, and YRSURF are now passed by argument
!   O. Marsden  Sept 2016 Added explicit spectral argument for passing to WRSPECA
!     -------------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE YOMGFL             , ONLY : TGFL
USE YOMCFU             , ONLY : TCFU
USE YOMXFU             , ONLY : TXFU
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCT3             , ONLY : NSTEP
USE YOMLUN             , ONLY : NULHWF, NULOUT
USE YOMMP0             , ONLY : MYPROC, LUSEWRGRIDALL, LOUTPUT
USE YOMOPH0            , ONLY : CFNHWF
USE YOMIO_SERV         , ONLY : IO_SERV_C001
USE FACTX_MOD          , ONLY : FACTX, NEWFACTX, FREEFACTX
USE MFIOOPTS_MOD       , ONLY : MFIOOPTS, MFIOOPTS_GETOPTS, MFIOFLAG, MFIOOPTS_GETFLAG
USE YOMMCUF            , ONLY : TMCUF
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (GEOMETRY),     INTENT(IN)    :: YDGEOMETRY
TYPE (TGFL),         INTENT(INOUT) :: YDGFL
TYPE (TSURF),        INTENT(INOUT) :: YDSURF
TYPE(SPECTRAL_FIELD),INTENT(INOUT) :: YDSPEC
TYPE(TCFU),          INTENT(INOUT) :: YDCFU
TYPE(TXFU),          INTENT(INOUT) :: YDXFU
TYPE(MODEL),         INTENT(INOUT) :: YDMODEL
CHARACTER (LEN=1),   INTENT(IN)    :: CDCONF
TYPE(TMCUF)        , INTENT(INOUT), OPTIONAL :: YDMCUF
!     ------------------------------------------------------------------
CHARACTER(LEN=16):: CLINC
CHARACTER(LEN=256) :: CLFIC

LOGICAL :: LLUSE_IOSERV

TYPE (MFIOFLAG)     :: YLFLAG
TYPE (MFIOOPTS)     :: YLOPTS
TYPE (FACTX)        :: YLFACTX

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "wrspeca.intfb.h"
#include "io_serv_log.intfb.h"
#include "wrgridall.intfb.h"
#include "wrgrida.intfb.h"
#include "wrgridua.intfb.h"
#include "wrxfu.intfb.h"
#include "wrfu.intfb.h"
#include "wrsfx.h"
#include "wrmlppa_io_serv.intfb.h"
#include "inifaoutinfo.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRMLPPA',0,ZHOOK_HANDLE)


!     ------------------------------------------------------------------

!*       1.   PREPARATIONS.
!             -------------

IF(LOUTPUT) WRITE(NULOUT,*) 'WRMLPPA NSTEP=',NSTEP,' CDCONF=',CDCONF

! currently used letters (in CY40): A,D,F,Q,x
IF (.NOT. ANY(SPREAD(CDCONF,1,6)==(/'A','D','F','Q','x','f'/))) &
 & CALL ABOR1(' WRMLPPA: wrong value for CDCONF ')

CALL MFIOOPTS_GETOPTS (YLOPTS, IO_SERV_C001)

LLUSE_IOSERV = IO_SERV_C001%LIO_SERV_WR

CALL NEWFACTX (YLFACTX, YLOPTS, CDCONF=CDCONF)
CLINC=''
CALL INIFAOUTINFO (YDMODEL%YRML_GCONF%YRRIP, YDXFU, CDCONF, CDINC=CLINC, CDFIC=CLFIC, KDATEF=YLFACTX%IDATEF)

IF (LLUSE_IOSERV) THEN

  IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'wrmlppa', 1_JPIM)

  IO_SERV_C001%IDATEF = YLFACTX%IDATEF

ENDIF

CALL MFIOOPTS_GETFLAG(YDCFU,YDXFU,YDMODEL%YRML_GCONF%YGFL,YDMODEL%YRML_PHY_MF%YRPHY,YLFLAG, &
 & CDCONF,NSTEP)

WRITE (NULOUT, '("WRMLPPA: CDCONF = ",A," LLWRSPECA = ",L2," LLWRGRIDA = ",L2,&
                & " LLWRGRIDUA = ",L2," LLWRXFU = ",L2," LLWRFU = ",L2)') &
                & CDCONF, YLFLAG%LLWRSPECA, YLFLAG%LLWRGRIDA, &
                & YLFLAG%LLWRGRIDUA, YLFLAG%LLWRXFU, YLFLAG%LLWRFU

IF (LLUSE_IOSERV) THEN

  IF (YLFLAG%LLWRSFX) THEN
    CALL WRSFX (YDMODEL%YRML_GCONF%YRRIP,YDGEOMETRY,YLFACTX,CLFIC)
  ELSE
    CALL WRMLPPA_IO_SERV(YDGEOMETRY, YDGFL, YDSURF, YDCFU, YDXFU, YDSPEC, YDMODEL, YLFLAG%LLWRGRIDA, YLFLAG%LLWRGRIDUA,&
 &                       YLFLAG%LLWRXFU,YLFLAG%LLWRFU,YLFLAG%LLWRSPECA,YDMCUF=YDMCUF)
  ENDIF

ELSE

  IF (YLFLAG%LLWRSFX) THEN
    CALL WRSFX (YDMODEL%YRML_GCONF%YRRIP,YDGEOMETRY,YLFACTX,CLFIC)
  ELSEIF (LUSEWRGRIDALL) THEN
    CALL WRGRIDALL(YDGEOMETRY,YDGFL,YDSURF,YDCFU,YDXFU,YDMODEL,YLFACTX,YLFLAG,CLFIC)
  ELSE
    IF (YLFLAG%LLWRGRIDA)  CALL WRGRIDA(YDGEOMETRY,YDSURF,YDXFU,YDMODEL,YLFACTX,CLFIC)
    IF (YLFLAG%LLWRGRIDUA) CALL WRGRIDUA(YDGEOMETRY,YDGFL,YDXFU,YDMODEL%YRML_GCONF,YLFACTX,CLFIC)
    IF (YLFLAG%LLWRXFU)    CALL WRXFU   (YDGEOMETRY,      YDXFU,YDMODEL%YRML_GCONF%YRRIP,YLFACTX,CLFIC)
    IF (YLFLAG%LLWRFU)     CALL WRFU    (YDGEOMETRY,      YDCFU,YDXFU,&
 &                                       YDMODEL%YRML_GCONF%YRRIP,YLFACTX,CLFIC)
  ENDIF
  
  IF (YLFLAG%LLWRSPECA ) THEN
    CALL WRSPECA(YDGEOMETRY, YDGFL, YDSURF, YDXFU, YDMODEL%YRML_GCONF, YDMODEL%YRML_DYN%YRDYN,YDMODEL%YRML_DYN%YRDYNA, &
     & YDSPEC, YLFACTX,CLFIC,YDMCUF=YDMCUF)
  ENDIF

ENDIF


!     ------------------------------------------------------------------

!*       4.    CLOSE FILES.
!              ------------

CALL FREEFACTX (YLFACTX)

IF ((MYPROC == 1.AND.NSTEP >= 0).AND.(CDCONF == 'A'.OR.CDCONF == 'F'.OR.CDCONF == 'Q'.OR.CDCONF=='f')) THEN  
  OPEN (UNIT=NULHWF,FILE=CFNHWF,FORM='FORMATTED')
  REWIND(NULHWF)
  WRITE(NULHWF,'(A)') TRIM (CLINC)
  CLOSE(UNIT=NULHWF,STATUS='KEEP')
ENDIF

IF (LLUSE_IOSERV) THEN
  IF (IO_SERV_C001%LDBUG) CALL IO_SERV_LOG (IO_SERV_C001, 'wrmlppa', 2_JPIM)
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('WRMLPPA',1,ZHOOK_HANDLE)

END SUBROUTINE WRMLPPA

