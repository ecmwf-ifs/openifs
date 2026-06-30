MODULE SUSRAD_MOD
CONTAINS
SUBROUTINE SUSRAD(KSW,KSIL,KTSW,KLWEMISS,&
 & LD_LLCCNL,LD_LLCCNO,KALBEDOSCHEME,KEMISSSCHEME,&
 & PTSTAND,PXP,PRCCNSEA,PRCCNLND,PRSUN,&
 & YDDIM,YDRAD,YDRDI,YDLW,YDSW)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_DIM   , ONLY : TDIM
USE YOS_RAD   , ONLY : TRAD
USE YOS_RDI   , ONLY : TRDI
USE YOS_LW    , ONLY : TLW
USE YOS_SW    , ONLY : TSW
USE ABORT_SURF_MOD

! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**   *SUSRAD* IS THE DRIVER TO ROUTINES SETTING-UP RADIATION CONSTANTS
!              This contains the fundamental model constants

!     INTERFACE.
!     ----------
!     CALL *SUSRAD* FROM *SUSURF*

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!     Original    P. Viterbo      May 2005

!     MODIFICATIONS
!     -------------
!     R. Hogan        23-01-2019  Longwave spectral emissivity in six intervals

USE SURWN_MOD

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN)    :: KSW ! Number of shortwave spectral albedo intervals
INTEGER(KIND=JPIM), INTENT(IN)    :: KSIL
INTEGER(KIND=JPIM), INTENT(IN)    :: KTSW 
INTEGER(KIND=JPIM), INTENT(IN)    :: KLWEMISS  ! Number of longwave spectral emissivity intervals
INTEGER(KIND=JPIM), INTENT(IN)    :: KALBEDOSCHEME
INTEGER(KIND=JPIM), INTENT(IN)    :: KEMISSSCHEME
LOGICAL           , INTENT(IN)    :: LD_LLCCNL 
LOGICAL           , INTENT(IN)    :: LD_LLCCNO  
REAL(KIND=JPRB)   , INTENT(IN)    :: PTSTAND 
REAL(KIND=JPRB)   , INTENT(IN)    :: PXP(6,6) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PRCCNSEA 
REAL(KIND=JPRB)   , INTENT(IN)    :: PRCCNLND 
REAL(KIND=JPRB)   , INTENT(IN)    :: PRSUN(:) 
TYPE(TDIM)        , INTENT(IN)    :: YDDIM
TYPE(TRAD)        , INTENT(INOUT) :: YDRAD
TYPE(TRDI)        , INTENT(INOUT) :: YDRDI
TYPE(TLW)         , INTENT(INOUT) :: YDLW
TYPE(TSW)         , INTENT(INOUT) :: YDSW

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUSRAD_MOD:SUSRAD',0,ZHOOK_HANDLE)
ASSOCIATE(LCCNL=>YDRAD%LCCNL, LCCNO=>YDRAD%LCCNO, NALBEDOSCHEME=>YDRAD%NALBEDOSCHEME, &
 & NEMISSSCHEME=>YDRAD%NEMISSSCHEME, &
 & NSW=>YDRAD%NSW, NTSW=>YDRAD%NTSW, NLWEMISS=>YDRAD%NLWEMISS, RCCNLND=>YDRAD%RCCNLND, &
 & RCCNSEA=>YDRAD%RCCNSEA, RALB_SNOW_FOREST=>YDRDI%RALB_SNOW_FOREST, &
 & RALBSEAD=>YDRDI%RALBSEAD, RALBSFO=>YDRDI%RALBSFO, REMISS_DESERT=>YDRDI%REMISS_DESERT, &
 & REMISS_LAND=>YDRDI%REMISS_LAND, REMISS_SNOW=>YDRDI%REMISS_SNOW, REMISS_SEA=>YDRDI%REMISS_SEA, &
 & REMISS_WEIGHT=>YDRDI%REMISS_WEIGHT, REMISS_OLD_WEIGHT=>YDRDI%REMISS_OLD_WEIGHT, &
 & REPALB=>YDRDI%REPALB, NUVVIS=>YDRAD%NUVVIS)

!  DEFINE RADIATION GENERAL SETTINGS

RCCNSEA=PRCCNSEA
RCCNLND=PRCCNLND

NTSW=KTSW
NSW=KSW

! Specify number of UV/Vis spectral intervals, the remainder being
! Near-IR. The boundary is at 0.69 microns; look at surwn_mod.F90 for
! more details.
IF (NSW == 1) THEN
  NUVVIS = 1
ELSEIF (NSW == 2) THEN
  NUVVIS = 1
ELSEIF (NSW == 4) THEN
  NUVVIS = 1 ! 1 UV/Vis band and 3 Near-IR bands
ELSEIF (NSW == 6) THEN
  NUVVIS = 3
ELSE
  CALL ABORT_SURF('SUSRAD_MOD: WRONG NUMBER OF SW INTERVALS')
ENDIF

NLWEMISS=KLWEMISS

LCCNL=LD_LLCCNL
LCCNO=LD_LLCCNO
NALBEDOSCHEME=KALBEDOSCHEME
NEMISSSCHEME =KEMISSSCHEME

! Define radiation lower boundary conditions settings
RALBSEAD = 0.06_JPRB ! Default diffuse sea albedo
REPALB=1.E-12_JPRB   ! Security
RALBSFO=0.15_JPRB    ! Default albedo of snow in presence of forest

! Default albedo for snow in the presence of forest: use previous
! spectral scaling
RALB_SNOW_FOREST(:,1) = YDRDI%RALBSFO*1.10598_JPRB ! Default UV/Vis albedo
RALB_SNOW_FOREST(:,2) = YDRDI%RALBSFO*0.90981_JPRB ! Default Near-IR albedo

! Moody et al. (RSE 2007) MODIS albedos for different forest types in
! the 0.3-0.7 micron band and 0.7 to 5.0 micron band.  In future it
! would be possible to define them separately in the three near-IR
! bands.
RALB_SNOW_FOREST(3,1:2) = [0.31_JPRB, 0.24_JPRB] ! Evergreen needle forest
RALB_SNOW_FOREST(4,1:2) = [0.39_JPRB, 0.27_JPRB] ! Deciduous needle forest
RALB_SNOW_FOREST(5,1:2) = [0.35_JPRB, 0.27_JPRB] ! Deciduous broad forest
RALB_SNOW_FOREST(6,1:2) = [0.44_JPRB, 0.33_JPRB] ! Evergreen broad forest
RALB_SNOW_FOREST(18,1:2)= [0.32_JPRB, 0.25_JPRB] ! Mixed forest
RALB_SNOW_FOREST(19,1:2)= [0.32_JPRB, 0.25_JPRB] ! Interupted forest (== mixed forest)

! Ensure sensible values outside NLWEMISS
REMISS_DESERT(:) = 0.93_JPRB
REMISS_LAND(:)   = 0.98_JPRB
REMISS_SNOW(:)   = 0.98_JPRB
REMISS_SEA(:)    = 0.94_JPRB
REMISS_WEIGHT(:) = 0.00_JPRB
REMISS_OLD_WEIGHT(:,:) = 0.0_JPRB

IF (NEMISSSCHEME == 0) THEN
  ! Old two-value scheme, where the first number is outside the IR
  ! window and the second is within, where using the old radiation
  ! scheme definitions, the window is 800-1250 cm-1.
  REMISS_DESERT(1:2) = [0.99_JPRB, 0.93_JPRB]
  REMISS_LAND  (1:2) = [0.99_JPRB, 0.96_JPRB]
  REMISS_SNOW  (1:2) = [0.99_JPRB, 0.98_JPRB]
  REMISS_SEA   (1:2) = [0.99_JPRB, 0.99_JPRB]

  ! Approximate weightings for broadband emissivity estimate using the
  ! Planck function at 15 degC.
  REMISS_WEIGHT(1:2) = [0.7173_JPRB, 0.2827_JPRB]

  ! 1-to-1 mapping from old to old scheme.
  REMISS_OLD_WEIGHT(1:2,1) = [1.0_JPRB, 0.0_JPRB]
  REMISS_OLD_WEIGHT(1:2,2) = [0.0_JPRB, 1.0_JPRB]

ELSEIF (NEMISSSCHEME == 1) THEN
  ! Emissivity values in six intervals corresponding to groups of RRTM
  ! bands, computed as a Planck-weighted average from the high
  ! resolution spectral data from Feldman et al., 2014: "Far-infrared
  ! surface emissivity and climate", Proc. Nat. Ac. Sci.

  ! Spectral interval number, wavelength bounds (um), RRTM bands
  ! covered, percentage of Planck function at 15 degC
  ! 1.       0-8.4746,  1-8,   14.9%
  ! 2.  8.4746-10.2041, 9-10,  11.1%
  ! 3. 10.2041-12.1951, 11,    12.6%
  ! 4. 12.1951-15.8730, 12-13, 19.0%
  ! 5. 15.8730-28.5714, 14-15, 29.0%
  ! 6. 28.5714-Inf,     16,    13.4%
  REMISS_DESERT(1:6) = [0.97395_JPRB, 0.89598_JPRB, 0.93831_JPRB, 0.96009_JPRB, 0.89529_JPRB, 0.91139_JPRB]
  REMISS_LAND  (1:6) = [0.98836_JPRB, 0.98179_JPRB, 0.98257_JPRB, 0.98685_JPRB, 0.98746_JPRB, 0.98746_JPRB]
  REMISS_SNOW  (1:6) = [0.97851_JPRB, 0.98896_JPRB, 0.98112_JPRB, 0.96727_JPRB, 0.98421_JPRB, 0.99240_JPRB]
  REMISS_SEA   (1:6) = [0.93780_JPRB, 0.94782_JPRB, 0.95285_JPRB, 0.91103_JPRB, 0.88590_JPRB, 0.87016_JPRB]

  ! Approximate weightings for broadband emissivity estimate; note
  ! that while this ought to change with temperature, uncertainties in
  ! these values are commensurate with the uncertainties in the
  ! emissivity values themselves, and a better diagnostic can be
  ! obtained from the radiation scheme itself
  REMISS_WEIGHT(1:6) = [0.149_JPRB, 0.111_JPRB, 0.126_JPRB, 0.190_JPRB, 0.290_JPRB, 0.134_JPRB]

  ! Approximate weightings for emissivity (column 1) outside
  ! atmospheric window and (column 2) inside atmospheric window
  ! (800-1250 cm-1), used to map back to old two-interval emissivity
  ! intervals still used by simplified-physics radiation scheme.
  REMISS_OLD_WEIGHT(1:6,1) = [0.18859_JPRB, 0.40422_JPRB, 0.23892_JPRB, 0.0_JPRB, 0.0_JPRB, 0.16827_JPRB]
  REMISS_OLD_WEIGHT(1:6,2) = [0.0_JPRB, 0.0_JPRB, 0.064137_JPRB, 0.44589_JPRB, 0.39119_JPRB, 0.098787_JPRB]

ELSE
  CALL ABORT_SURF('SUSRAD_MOD: NEMISSSCHEME must be 0 or 1')
END IF

! Initialise radiation SW and LW constants

CALL SURWN(KSIL,PTSTAND,PXP,PRSUN,YDDIM,YDRAD,YDLW,YDSW)

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUSRAD_MOD:SUSRAD',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE SUSRAD
END MODULE SUSRAD_MOD
