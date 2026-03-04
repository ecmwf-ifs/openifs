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

SUBROUTINE SUDIM(YDGEOMETRY,KSUPERSEDE)

!     ----------------------------------------------------------------------
!**** *SUDIM * - Setting up of the geometry dimensioning of the model.

!     Purpose.
!     --------
!           Initialization of YOMDIM and YOMDIMV, and some printings.

!**   Interface.
!     ----------
!        *CALL* *SUDIM*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 87-10-15

!     Modifications.
!     --------------
!      Modified by K. Yessad: 31-01-2007: modularisation and big cleanings.
!      J.Haseler     27-Feb-2007 Modify defaults for LTRAJIO and LGPINGP
!      JJMorcrette   20060925    DU, BC, OM, SO2, VOL, SOA climatological flds
!      JJMorcrette   20061020    diagnostic of aerosol physical fluxes
!      22-Mar-2007 S.Serrar  : defaults for ERA40 GFL fields
!      02-MAY-2008 S.Serrar  : check non consistent values for NERA40
!      18-Dec-2007   Y. Seity : add Hail
!      G.Balsamo     20081016    add lake and ocean dimensions
!      JJMorcrette   20090217    PP of aerosol and UV processor output fields
!      Y.Tremolet    27-Jan-2009 Model error for stratosphere
!      30-Jun-2008 J. Masek   New SLHD interpolators.
!      Modified by K Yessad : 15-Sep-08 Prune conf 951.
!      A. Fouilloux  09/09/2009: remove call to sudimo
!      K. Yessad (Aug 2009): prune conf 912, externalise conf 911.
!      15-Oct-2009 Y. Bouteloup : add radiative cloud water and ice (YIRAD and YLRAD)
!      Modified by A. Alias : 15-Oct-09 new keys in MCC : LCURR , LGELATO
!      K. Yessad (Jan 2010): remove useless variables.
!      F. Vana    22-Feb-2011 LUVDER = .true. when 3D turbulence
!      H. Hersbach   20110401: auxiliary diagnostic radiation fields
!      J.Hague     21-Mar-11: YDGOM Derived Type added
!      J. Hague May-2011:   NUSE_ECCI Setup moved to SUGOMS
!      R. El Khatib 10-Aug-2011 NIOLEVG management
!      K. YESSAD (Nov 2011): move GFL set-up in SUGFL1.
!      K. YESSAD (Nov 2011): move (NFLSUL,NFLSA,NFLEN) calc. from SUDIM2 to SUDIM1.
!      R. El Khatib 22-Mar-2012 Fix uninitialized variables + bounds checking issue
!      G. Mozdzynski Aug-2011: support higher order interpolation
!      E. Holm Jul-2012: Model error forcing setting for L137
!      K. Yessad (July 2013): minor modifications in printings (printing and setup in the same routine).
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring (merge former SUDIM1+SUDIM2).
!      K. Yessad (Dec 2016): Prune obsolete options.
!     ----------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE PARDIM   , ONLY : JPSLWIDE
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMARG   , ONLY : NULIM, NSUPERSEDE, NUCMAX, NUDGL, NUDLON, NUFLEV, NUSMAX, NUSTTYP, USTRET  
USE YOMMP0   , ONLY : LOPT_SCALAR
USE YOMLUN   , ONLY : NULOUT, NULNAM, NULERR
USE YOMCT0   , ONLY : LR2D, NCONF, LELAM, LIOLEVG, LECMWF
USE YOMVAR   , ONLY : LREPRO4DVAR
USE YOMMP0   , ONLY : MYPROC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT),TARGET :: YDGEOMETRY
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KSUPERSEDE

INTEGER(KIND=JPIM) :: IPROMA, IAL1, IAL2, IDLSUR
LOGICAL :: LLGRID
LOGICAL :: LLSUPERSEDE
REAL(KIND=JPRB) :: ZEPS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), POINTER :: NDGLG, NDLON, NDGUXG, NDLUXG,&
 & NFLEVG, NIOLEVG, NVARMAX, NMSMAX, NSMAX, NCMAX, NPROMA,&
 & NSTENCILWIDE

#include "namdim.nam.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "suedim.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUDIM',0,ZHOOK_HANDLE)

ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 & YDEDIM=>YDGEOMETRY%YREDIM)
ASSOCIATE(NNOEXTZG=>YDEDIM%NNOEXTZG, NNOEXTZL=>YDEDIM%NNOEXTZL, &
 & LOPTPROMA=>YDDIM%LOPTPROMA, NDGNH=>YDDIM%NDGNH, NDLSM=>YDDIM%NDLSM, &
 & NDLSUR=>YDDIM%NDLSUR, NDLUNG=>YDDIM%NDLUNG, NDGENG=>YDDIM%NDGENG, &
 & NDGSAG=>YDDIM%NDGSAG, NDGSUR=>YDDIM%NDGSUR, NDGUNG=>YDDIM%NDGUNG, &
 & NDSUR1=>YDDIM%NDSUR1, NSEFRE=>YDDIM%NSEFRE, NSPEC2G=>YDDIM%NSPEC2G, &
 & NSPECG=>YDDIM%NSPECG, &
 & NFLEN=>YDDIMV%NFLEN, NFLSA=>YDDIMV%NFLSA, NFLSUL=>YDDIMV%NFLSUL)

! Associate pointers for variables in namelist
NDGLG        => YDDIM%NDGLG
NDLON        => YDDIM%NDLON
NDGUXG       => YDDIM%NDGUXG
NDLUXG       => YDDIM%NDLUXG
NFLEVG       => YDDIMV%NFLEVG
NIOLEVG      => YDDIMV%NIOLEVG
NVARMAX      => YDDIM%NVARMAX
NMSMAX       => YDDIM%NMSMAX
NSMAX        => YDDIM%NSMAX
NCMAX        => YDDIM%NCMAX
NPROMA       => YDDIM%NPROMA
NSTENCILWIDE => YDDIM%NSTENCILWIDE

!     ------------------------------------------------------------------

!*       0. PRELIMINAR CALCULATIONS.
!        ---------------------------

ZEPS=1.E-12_JPRB

IF (PRESENT(KSUPERSEDE)) THEN
  LLSUPERSEDE=(KSUPERSEDE == 1)
ELSE
  LLSUPERSEDE=(NSUPERSEDE == 1)
ENDIF

!     ------------------------------------------------------------------
!*       1. SET DEFAULT VALUES AND MODIFY THEM.
!           -----------------------------------

!*       1.1 SET DEFAULT VALUES FOR HORIZONTAL GEOMETRY.

NDGLG =32
NDLON =64
NSMAX =21
NCMAX = 1

IF(LOPT_SCALAR) THEN
  NPROMA=16
ELSE
  NPROMA=2047
ENDIF
IPROMA=NPROMA

NDGUNG= 1
NDGUXG= NDGLG
NDLUNG= 1
NDLUXG= NDLON

NMSMAX = NSMAX

! Set default max stencil width / 2
NSTENCILWIDE=2

IF (LECMWF) THEN
  IF(LLSUPERSEDE) THEN
    NSMAX=NUSMAX
  ENDIF
ELSE
  IF (LLSUPERSEDE) THEN
    NDGLG=NUDGL
    NDLON=NUDLON
    NSMAX = NUSMAX
    NCMAX = NUCMAX
    IF (ANY(SPREAD(NCONF,1,4) == (/1,302,601,801/))) THEN
      IF (ABS(USTRET-1.0_JPRB) <= ZEPS) THEN
        NCMAX = NUSMAX
        NUCMAX = NUSMAX
      ENDIF
    ELSEIF (NCONF == 701) THEN
      NCMAX = NUSMAX
      NUCMAX = NUSMAX
    ENDIF
    IF (LELAM) THEN
      NMSMAX=NUSTTYP
      NVARMAX=NMSMAX
      NDLUNG=NULIM(3)
      NDLUXG=NULIM(4)
      NDGUNG=NULIM(5)
      NDGUXG=NULIM(6)
    ENDIF
  ENDIF
ENDIF

IF (.NOT. LELAM) THEN
  NMSMAX =NSMAX
  NVARMAX=NSMAX
ENDIF

!*       1.2 SET DEFAULT VALUES FOR VERTICAL GEOMETRY.

NFLEVG=19
NIOLEVG=19

IF (.NOT.LECMWF) THEN
  IF (LLSUPERSEDE) THEN
    NFLEVG = NUFLEV
    NIOLEVG = NUFLEV
  ENDIF
ENDIF

!     ------------------------------------------------------------------
!*       2. READ NAMELIST.
!           --------------

CALL POSNAM(NULNAM,'NAMDIM')
READ(NULNAM,NAMDIM)

!     ------------------------------------------------------------------
!*       3. DO CHECKINGS AND RESETTINGS FOR HORIZONTAL GEOMETRY.
!           -----------------------------------------------------

IF (NPROMA==0) NPROMA=IPROMA  ! restore default nproma

!        3.1    Determine whether NPROMA is negative and set LOPTPROMA

LOPTPROMA=(NPROMA > 0) .AND. (.NOT.LREPRO4DVAR)
NPROMA=IABS(NPROMA)

!        3.2    Some consistency checks and resettings.
!               Final values of NDLON, NDGLG, NSMAX, NMSMAX must be computed there, and not modified later.

IF (LLSUPERSEDE) THEN
  IF (NSMAX /= NUSMAX) THEN
    NSMAX=NUSMAX
    WRITE(NULOUT,*) 'NSMAX   OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (NDGLG /= NUDGL) THEN
    NDGLG=NUDGL
    WRITE(NULOUT,*) 'NDGLG  OVERWRITTEN BY FILE FRAME'
  ENDIF
  IF (NDLON /= NUDLON) THEN
    NDLON=NUDLON
    WRITE(NULOUT,*) 'NDLON   OVERWRITTEN BY FILE FRAME'
  ENDIF

  ! * LAM models settings.
  IF (LELAM) THEN
    IF (NMSMAX /= NUSTTYP) THEN
      NMSMAX=NUSTTYP
      WRITE(NULOUT,*) 'NMSMAX   OVERWRITTEN BY FILE FRAME'
    ENDIF
  ENDIF
ENDIF

!*       3.3    Check for aliasing and use of Eulerian scheme on linear grid.
!               Write warnings for not recommended settings, but no abort there!

IAL1=3*NSMAX+2-(2*NDGLG-1)
IAL2=3*NSMAX+2-(NDLON-1)
IF(IAL1 > 0.OR.IAL2 > 0)THEN
  WRITE(UNIT=NULOUT,FMT='('' ******* ALIASING IN THE STRETCHED ''&
   & ,'' COORDINATE MODEL ******** '')')  
ENDIF
IF(IAL1 > 2.OR.IAL2 > 2)THEN
  WRITE(UNIT=NULOUT,FMT='('' ******* ALIASING IN THE UNIFORM   ''&
   & ,'' COORDINATE MODEL ******** '')')  
ENDIF

IF( .NOT.LELAM )THEN
  LLGRID=NSMAX > (NDLON+3)/3
ELSE
  LLGRID=NSMAX > (NDGLG -1)/3 .OR. NMSMAX > (NDLON -1)/3
ENDIF

!        3.4    Final value of NCMAX.

IF (LLSUPERSEDE) THEN
  IF (NCMAX /= NUCMAX) THEN
    NCMAX=NUCMAX
    WRITE(NULOUT,*) 'NCMAX   OVERWRITTEN BY FILE FRAME'
  ENDIF
ENDIF
NCMAX = MAX(NCMAX,NSMAX)
IF (LELAM) NCMAX = MAX(NCMAX,NMSMAX)

!*       3.5    Initialize collocation grid dimensioning: NDGNH, NDGSUR, NDGSAG, NDGENG.

NDGNH =(NDGLG+1)/2
IF (LELAM) THEN
  NDGSUR=0
ELSE
  ! number of extra-latitudes needed for horizontal interpolations
  NDGSUR=NSTENCILWIDE
ENDIF
NDGSAG=1-NDGSUR
NDGENG=NDGLG+NDGSUR

IF (NDGSUR>JPSLWIDE .OR. NSTENCILWIDE>JPSLWIDE) THEN
  CALL ABOR1(' SUDIM: NDGSUR and NSTENCILWIDE should be <= JPSLWIDE ')
ENDIF

!*       3.6    Final value of NDGUNG, NDGUXG, NDLUNG, NDLUXG.

IF (LELAM) THEN
  IF (LLSUPERSEDE) THEN
    IF (NDLUNG /= NULIM(3)) THEN
      NDLUNG=NULIM(3)
      WRITE(NULOUT,*) 'NDLUNG  OVERWRITTEN BY FILE FRAME'
    ENDIF
    IF (NDLUXG /= NULIM(4)) THEN
      NDLUXG=NULIM(4)
      WRITE(NULOUT,*) 'NDLUXG  OVERWRITTEN BY FILE FRAME'
    ENDIF
    IF (NDGUNG /= NULIM(5)) THEN
      NDGUNG=NULIM(5)
      WRITE(NULOUT,*) 'NDGUNG  OVERWRITTEN BY FILE FRAME'
    ENDIF
    IF (NDGUXG /= NULIM(6)) THEN
      NDGUXG=NULIM(6)
      WRITE(NULOUT,*) 'NDGUXG  OVERWRITTEN BY FILE FRAME'
    ENDIF
  ENDIF
ELSE
  NDGUNG=NDGSAG
  NDGUXG=NDGENG
  NDLUNG=1
  NDLUXG=NDLON
ENDIF

!*       3.7    Compute NSEFRE, NSPECG, NSPEC2G.

IF (LELAM) THEN
  CALL SUEDIM(YDGEOMETRY%YRDIM,YDGEOMETRY%YREDIM,NSEFRE,NSPECG,KSUPERSEDE)
ELSE
  NSEFRE=(NSMAX+1)*(NSMAX+1)
  NSPECG=(NSMAX+1)*(NSMAX+2)/2
ENDIF
NSPEC2G=NSPECG*2

!*       3.8    Compute NDSUR1, NDLSUR, NDLSM (must be done after calling SUEDIM).

IF (LELAM) THEN
  IF (NNOEXTZG>0.OR.NNOEXTZL>0) THEN
    IDLSUR=NDLON
  ELSE
    IDLSUR=MAX(NDLON,2*NMSMAX+1)
  ENDIF
ELSE
  IDLSUR=MAX(NDLON,2*NSMAX+1)
ENDIF
NDSUR1=NSTENCILWIDE*2-1
NDLSUR=IDLSUR+NDSUR1
NDLSM=NDLSUR-1

!     ------------------------------------------------------------------
!*       4. DO CHECKINGS AND RESETTINGS FOR VERTICAL GEOMETRY.
!           --------------------------------------------------

!        4.1    Checkings and final calculation for NFLEVG and NIOLEVG.

IF (LLSUPERSEDE) THEN
  IF (LIOLEVG) THEN
    IF (NFLEVG /= NUFLEV) THEN
      NFLEVG=NUFLEV
      WRITE(NULOUT,*) 'NFLEVG   OVERWRITTEN BY FILE FRAME'
    ENDIF
    IF (NIOLEVG /= NUFLEV) THEN
      NIOLEVG=NUFLEV
      WRITE(NULOUT,*) 'NIOLEVG   OVERWRITTEN BY FILE FRAME'
    ENDIF
  ENDIF
ELSE
  IF (LIOLEVG) THEN
    NIOLEVG=NFLEVG
  ENDIF
ENDIF

! Reset NFLEVG and NIOLEVG to 1 in some cases.
IF(LR2D) THEN  
  NFLEVG=1
  NIOLEVG=1
ELSEIF (ANY(SPREAD(NCONF,1,4) == (/923,931,932,933/))) THEN
  NFLEVG=1
  NIOLEVG=1
ENDIF

! Check that NFLEVG <= NIOLEVG.
IF (NFLEVG > NIOLEVG) THEN
  WRITE(NULERR,'(''NFLEVG BIGGER THAN NIOLEVG'')')
  WRITE(NULOUT,'('' NFLEVG = '',I3,'' NIOLEVG = '',I3)') NFLEVG,NIOLEVG
  CALL ABOR1(' SUDIM: NFLEVG BIGGER THAN NIOLEVG ')
ENDIF

!        4.2    Compute NFLSUL, NFLSA and NFLEN from the final value of NFLEVG.

IF(LR2D) THEN
  NFLSUL=0
ELSE
!  IF(LSLAG) THEN ! Always assume LSLAG to avoid dependence for geaometry on
                  ! the choise of advection (which belongs to model object) MH               
    NFLSUL=1
!  ELSE
!    NFLSUL=0
!  ENDIF
ENDIF
NFLSA=1-NFLSUL
NFLEN=NFLEVG+NFLSUL

!     ------------------------------------------------------------------
!*       5. PRINTINGS.
!           ----------

!*       5.1  Print vertical dimensions

WRITE(NULOUT,*) ' PRINTINGS IN SUDIM: VERTICAL DIMENSIONING '
WRITE(UNIT=NULOUT,FMT='('' NFLEVG ='',I6,'' NIOLEVG ='',I6,&
 & '' NFLSUL='',I6,'' NFLSA ='',I6,'' NFLEN ='',I6)')&
 & NFLEVG,NIOLEVG,NFLSUL,NFLSA,NFLEN

!*       5.2  Print horizontal dimensions

WRITE(NULOUT,*) ' PRINTINGS IN SUDIM: HORIZONTAL DIMENSIONING '
WRITE(UNIT=NULOUT,FMT='('' NDGLG ='',I6)') NDGLG
WRITE(NULOUT,*) 'NSTENCILWIDE=',NSTENCILWIDE
WRITE(UNIT=NULOUT,FMT='('' NDGNH ='',I6,&
 & '' NDGSUR='',I6,'' NDGSAG ='',I6,'' NDGENG ='',I6,/,&
 & '' NDLON ='',I6,'' NDSUR1='',I6,'' NDLSUR='',I6,'' NDLSM ='',I6,/,&
 & '' NDGUNG ='',I6,'' NDGUXG ='',I6,'' NDLUNG ='',I6,'' NDLUXG ='',I6,/,&
 & '' NPROMA='',I6)')&
 & NDGNH ,NDGSUR,NDGSAG,NDGENG,NDLON ,NDSUR1,NDLSUR,NDLSM,&
 & NDGUNG,NDGUXG,NDLUNG,NDLUXG,NPROMA
WRITE(NULOUT,'('' NSMAX ='',I6,'' NMSMAX='',I6,&
 & '' NSEFRE='',I8,'' NSPECG ='',I8,'' NSPEC2G='',I8)')&
 & NSMAX,NMSMAX,NSEFRE,NSPECG,NSPEC2G
WRITE(NULOUT,'('' NVARMAX ='',I6)') NVARMAX
WRITE(NULOUT,'('' NCMAX ='',I6)') NCMAX

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUDIM',1,ZHOOK_HANDLE)
END SUBROUTINE SUDIM

