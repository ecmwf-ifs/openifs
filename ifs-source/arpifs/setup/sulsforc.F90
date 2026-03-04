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

SUBROUTINE SULSFORC(YDML_GCONF,KULOUT)

!**** *SULSFORC*   - Initialize forcing parameters

!     Purpose.
!     --------
!           Initialize YOMLSFORC and checkings

!**   Interface.
!     ----------
!        *CALL* *SULSFORC(KULOUT)

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation academic case and 1D model (SCUM)

!     Author.
!     -------
!        S. Malardel
!        Original : 06-01-04

!     Modifications.
!     --------------
!        E. Bazile and I. Beau : 2011-01-18 Nudging options.
!        K. Yessad (July 2014): Move some variables.
!        E. Bazile (March 2016): variable nudging
!        R. Roehrig (Sept 2018): add LSPS_FRC (imposing surface pressure evolution in MUSC)
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULNAM
USE YOMCT0   , ONLY : LSFORC
USE YOMLSFORC

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF
INTEGER(KIND=JPIM),INTENT(IN)   :: KULOUT ! logical unit of output file

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: ITOT_FORC, IFORCECH
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
!     ------------------------------------------------------------------

#include "posnam.intfb.h"
#include "abor1.intfb.h"

!     ------------------------------------------------------------------

#include "namlsforc.nam.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SULSFORC',0,ZHOOK_HANDLE)
ASSOCIATE(NGFL_FORC=>YDML_GCONF%YGFL%NGFL_FORC, &
 & NSTOP=>YDML_GCONF%YRRIP%NSTOP, TSTEP=>YDML_GCONF%YRRIP%TSTEP)
!     ------------------------------------------------------------------

!         1. Set implicit default values for YOMLSFORC

LGEOST_UV_FRC = .FALSE.
RCORIO_FORC = 0._JPRB
RZ0_FORC = 0._JPRB
RALB_FORC = 0.15
REMIS_FORC = 0.98
NL_GEOST_UV_TIME(:) =999
NGEOST_UV_FREQ = 999
NGEOST_U_DEB = 1
NGEOST_U_NUM = 0
NGEOST_V_DEB = 1
NGEOST_V_NUM = 0

LUV_ADV_FRC = .FALSE.
NL_UV_ADV_TIME(:) = 999
NUV_ADV_FREQ = 999
NU_ADV_DEB = 1
NU_ADV_NUM = 0
NV_ADV_DEB = 1
NV_ADV_NUM = 0

LT_ADV_FRC = .FALSE.
NL_T_ADV_TIME(:) = 999
NT_ADV_FREQ = 999
NT_ADV_DEB = 1
NT_ADV_NUM = 0

LQV_ADV_FRC = .FALSE.
NL_QV_ADV_TIME(:) = 999
NQV_ADV_FREQ = 999
NQV_ADV_DEB = 1
NQV_ADV_NUM = 0

LSW_FRC = .FALSE.
NL_LSW_TIME(:) = 999
NLSW_FREQ = 999
NLSW_DEB = 1
NLSW_NUM = 0

LSOMEGA_FRC = .FALSE.
NL_LSOMEGA_TIME(:) = 999
NLSOMEGA_FREQ = 999
NLSOMEGA_DEB = 1
NLSOMEGA_NUM = 0

! Nudging
LT_NUDG = .FALSE.
LQV_NUDG = .FALSE.
LUV_NUDG = .FALSE.
NL_T_NUDG_TIME(:) = 999
NL_QV_NUDG_TIME(:) = 999
NL_UV_NUDG_TIME(:) = 999
RELAX_TAUT = 0._JPRB
RELAX_TAUQ = 0._JPRB
RELAX_TAUU = 0._JPRB
NT_NUDG_FREQ=0
NQV_NUDG_FREQ=0
NUV_NUDG_FREQ=999
NT_NUDG_DEB=0
NQV_NUDG_DEB=0
NU_NUDG_DEB=0
NV_NUDG_DEB=0
NT_NUDG_NUM=0
NQV_NUDG_NUM=0
NUV_NUDG_NUM=0


!  Surface forcing
NSH_FORC_DEB = 1
NSH_FORC_NUM = 0
NL_SH_ADV_TIME(:) = 999
NLH_FORC_DEB = 1
NLH_FORC_NUM = 0
NL_LH_ADV_TIME(:) = 999
NTS_FORC_DEB = 1
NTS_FORC_NUM = 0
NL_TS_ADV_TIME(:) = 999
NUS_FORC_DEB = 1
NUS_FORC_NUM = 0
NL_US_ADV_TIME(:) = 999

LMUSCLFA=.FALSE.
NMUSCLFA=86

! Surface pressure
LSPS_FRC=.FALSE.

!         2. Read namelist for forcing

IF (LSFORC) THEN
  CALL POSNAM(NULNAM,'NAMLSFORC')
  READ (NULNAM,NAMLSFORC)
ENDIF

!         3. Check namelist values

ITOT_FORC = NGEOST_U_NUM + NGEOST_V_NUM&
 & + NU_ADV_NUM + NV_ADV_NUM + NT_ADV_NUM + NQV_ADV_NUM&
 & + MAX(NLSW_NUM,NLSOMEGA_NUM)
IF (ITOT_FORC > NGFL_FORC) CALL ABOR1('Size of NGFL_FORC too small')

IF (LSW_FRC .AND. LSOMEGA_FRC) THEN
  CALL ABOR1('You can not have a large scale vertical advection prescibed both in w and omega')
ENDIF

IF (LT_NUDG.AND.(RELAX_TAUT==0._JPRB)) THEN
   CALL ABOR1('Nudging case: you need to initialize the relaxation time RELAX_TAUT ')
ENDIF
IF (LQV_NUDG.AND.(RELAX_TAUQ==0._JPRB)) THEN
   CALL ABOR1('Nudging case: you need to initialize the relaxation time RELAX_TAUQ ')
ENDIF
IF (LUV_NUDG.AND.(RELAX_TAUU==0._JPRB)) THEN
   CALL ABOR1('Nudging case: you need to initialize the relaxation time RELAX_TAUU ')
ENDIF
IF (LT_NUDG.AND.(NT_NUDG_DEB == 0)) THEN 
   CALL ABOR1 ('NO TEMPERATURE PROFILE FOR NUDGING : NT_NUDG=0, SPECIFY THE RANK OF YOUR PROFILE IN THE GFL FORC')
ENDIF
IF (LQV_NUDG.AND.(NQV_NUDG_DEB == 0)) THEN 
   CALL ABOR1 ('NO HUMIDITY PROFILE FOR NUDGING : NQV_NUDG=0, SPECIFY THE RANK OF YOUR PROFILE IN THE GFL FORC')
ENDIF
IF (LUV_NUDG.AND.((NU_NUDG_DEB == 0).AND.(NV_NUDG_DEB==0))) THEN 
   CALL ABOR1 ('NO WIND PROFILE FOR NUDGING : NU/V_NUDG=0, SPECIFY THE RANK OF YOUR PROFILE IN THE GFL FORC')
ENDIF

IFORCECH= NGEOST_UV_FREQ*(NGEOST_U_NUM-1)
IF ( (NGEOST_UV_FREQ/=999) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NGEOST_UV_FREQ * NGEOST_U_NUM is lower than the forecast length')
ENDIF
IF ( (NGEOST_UV_FREQ/=999) .AND.(NL_GEOST_UV_TIME(1)/=999) ) THEN
  CALL ABOR1('NGEOST_UV : You can not prescribe both a frequency and a time table for the forcing')
ENDIF

IFORCECH= NUV_ADV_FREQ*(NU_ADV_NUM-1)
IF ( (NUV_ADV_FREQ/=999 ) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NUV_ADV_FREQ *  NU_ADV_NUM is lower than the forecast length')
ENDIF
IF ( (NUV_ADV_FREQ/=999 ) .AND. (NL_UV_ADV_TIME(1)/=999 ) ) THEN
  CALL ABOR1(' NUV_ADV: You can not prescribe both a frequency and a time table for the forcing')
ENDIF

IFORCECH= NT_ADV_FREQ*(NT_ADV_NUM-1)
IF ( (NT_ADV_FREQ/=999 ) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NT_ADV_FREQ * NT_ADV_NUM is lower than the forecast length')
ENDIF
IF ( (NT_ADV_FREQ/=999 ) .AND. (NL_T_ADV_TIME(1)/=999 ) ) THEN
  CALL ABOR1(' NT_ADV: You can not prescribe both a frequency and a time table for the forcing')
ENDIF

IFORCECH= NQV_ADV_FREQ*(NQV_ADV_NUM-1)
IF ( (NQV_ADV_FREQ/=999 ) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NQV_ADV_FREQ * NQV_ADV_NUM is lower than the forecast length')
ENDIF
IF ( (NQV_ADV_FREQ/=999 ) .AND. (NL_QV_ADV_TIME(1)/=999 ) ) THEN
  CALL ABOR1(' NQV_ADV: You can not prescribe both a frequency and a time table for the forcing')
ENDIF

IFORCECH= NLSW_FREQ*(NLSW_NUM-1)
IF ( (NLSW_FREQ/=999 ) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NLSW_FREQ * NLSW_NUM is lower than the forecast length')
ENDIF
IF ( (NLSW_FREQ/=999 ) .AND. (NL_LSW_TIME(1)/=999 ) ) THEN
  CALL ABOR1(' NLSW: You can not prescribe both a frequency and a time table for the forcing')
ENDIF

IFORCECH= NLSOMEGA_FREQ*(NLSOMEGA_NUM-1)
IF ( (NLSOMEGA_FREQ/=999 ) .AND. (IFORCECH.LT.(NSTOP*INT(TSTEP))) ) THEN
  CALL ABOR1('NLSOMEGA_FREQ * NLSOMEGA_NUM  is lower than the forecast length')
ENDIF
IF ( (NLSOMEGA_FREQ/=999 ) .AND. (NL_LSOMEGA_TIME(1)/=999 ) ) THEN
  CALL ABOR1(' NLSOMEGA: You can not prescribe both a frequency and a time table for the forcing')
ENDIF

!         4. Transfert local namelist array NL_xxx_TIME to module array  Nxxx_TIME

IF ( NL_GEOST_UV_TIME(1)/=999 ) THEN
  ALLOCATE(NGEOST_UV_TIME(NGEOST_U_NUM))
  NGEOST_UV_TIME(1:NGEOST_U_NUM)=NL_GEOST_UV_TIME(1:NGEOST_U_NUM)
ENDIF

IF (NL_UV_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NUV_ADV_TIME(NU_ADV_NUM))
  NUV_ADV_TIME(1:NU_ADV_NUM)=NL_UV_ADV_TIME(1:NU_ADV_NUM)
ENDIF

IF (NL_T_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NT_ADV_TIME(NT_ADV_NUM))
  NT_ADV_TIME(1:NT_ADV_NUM)=NL_T_ADV_TIME(1:NT_ADV_NUM)
ENDIF

IF (NL_QV_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NQV_ADV_TIME(NQV_ADV_NUM))
  NQV_ADV_TIME(1:NQV_ADV_NUM)=NL_QV_ADV_TIME(1:NQV_ADV_NUM)
ENDIF

IF (NL_LSW_TIME(1)/=999 ) THEN
  ALLOCATE(NLSW_TIME(NLSW_NUM))
  NLSW_TIME(1:NLSW_NUM)=NL_LSW_TIME(1:NLSW_NUM)
ENDIF

IF (NL_LSOMEGA_TIME(1)/=999 ) THEN
  ALLOCATE(NLSOMEGA_TIME(NLSOMEGA_NUM))
  NLSOMEGA_TIME(1:NLSOMEGA_NUM)=NL_LSOMEGA_TIME(1:NLSOMEGA_NUM)
ENDIF

IF (NL_T_NUDG_TIME(1)/=999 ) THEN
ALLOCATE(NT_NUDG_TIME(NT_NUDG_NUM))
NT_NUDG_TIME(1:NT_NUDG_NUM)=NL_T_NUDG_TIME(1:NT_NUDG_NUM)
ENDIF

IF (NL_QV_NUDG_TIME(1)/=999 ) THEN
ALLOCATE(NQV_NUDG_TIME(NQV_NUDG_NUM))
NQV_NUDG_TIME(1:NQV_NUDG_NUM)=NL_QV_NUDG_TIME(1:NQV_NUDG_NUM)
ENDIF

IF (NL_UV_NUDG_TIME(1)/=999 ) THEN
ALLOCATE(NUV_NUDG_TIME(NUV_NUDG_NUM))
NUV_NUDG_TIME(1:NUV_NUDG_NUM)=NL_UV_NUDG_TIME(1:NUV_NUDG_NUM)
ENDIF

! Surface ...

IF (NL_SH_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NT_SH_ADV_TIME(NSH_FORC_NUM))
  NT_SH_ADV_TIME(1:NSH_FORC_NUM)=NL_SH_ADV_TIME(1:NSH_FORC_NUM)
ENDIF

IF (NL_LH_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NT_LH_ADV_TIME(NLH_FORC_NUM))
  NT_LH_ADV_TIME(1:NLH_FORC_NUM)=NL_LH_ADV_TIME(1:NLH_FORC_NUM)
ENDIF

IF (NL_TS_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NT_TS_ADV_TIME(NTS_FORC_NUM))
  NT_TS_ADV_TIME(1:NTS_FORC_NUM)=NL_TS_ADV_TIME(1:NTS_FORC_NUM)
ENDIF

IF (NL_US_ADV_TIME(1)/=999 ) THEN
  ALLOCATE(NT_US_ADV_TIME(NUS_FORC_NUM))
  NT_US_ADV_TIME(1:NUS_FORC_NUM)=NL_US_ADV_TIME(1:NUS_FORC_NUM)
ENDIF


!         5. Check print

WRITE(UNIT=KULOUT,FMT='(&
 & ''LMUSCLFA = '',L5,'' NMUSCLFA = '',I3, '' LGEOST_UV_FRC = '',L5,'' RCORIO_FORC = '',E14.7,&
 & '' NGEOST_UV_FREQ = '',I5,'' NGEOST_UV_DEB = '',I3,'' NGEOST_UV_NUM = '',I3 )')&
 & LMUSCLFA, NMUSCLFA, LGEOST_UV_FRC, RCORIO_FORC, NGEOST_UV_FREQ, NGEOST_U_DEB, NGEOST_U_NUM

WRITE(UNIT=KULOUT,FMT='(''  LUV_ADV_FRC = '',L5,&
 & '' NUV_ADV_FREQ = '',I5,'' NU_ADV_DEB = '',I3,'' NU_ADV_NUM = '',I3,&
 & '' NV_ADV_DEB = '',I3,'' NV_ADV_NUM = '',I3 )')&
 & LUV_ADV_FRC, NUV_ADV_FREQ, NU_ADV_DEB, NU_ADV_NUM, NV_ADV_DEB, NV_ADV_NUM

WRITE(UNIT=KULOUT,FMT='(''  LT_ADV_FRC = '',L5,&
 & '' NT_ADV_FREQ = '',I5,'' NT_ADV_DEB = '',I3,'' NT_ADV_NUM = '',I3 )')&
 & LT_ADV_FRC, NT_ADV_FREQ, NT_ADV_DEB, NT_ADV_NUM

WRITE(UNIT=KULOUT,FMT='(''  LQV_ADV_FRC = '',L5,&
 & '' NQV_ADV_FREQ = '',I5,'' NQV_ADV_DEB = '',I3,'' NQV_ADV_NUM = '',I3 )')&
 & LQV_ADV_FRC, NQV_ADV_FREQ, NQV_ADV_DEB, NQV_ADV_NUM

WRITE(UNIT=KULOUT,FMT='(''  LSW_FRC = '',L5,&
 & '' NLSW_FREQ = '',I5,'' NLSW_DEB = '',I3,'' NLSW_NUM = '',I3 )')&
 & LSW_FRC, NLSW_FREQ, NLSW_DEB, NLSW_NUM

WRITE(UNIT=KULOUT,FMT='(''  LSOMEGA_FRC = '',L5,&
 & '' NLSOMEGA_FREQ = '',I5,'' NLSOMEGA_DEB = '',I3,'' NLSOMEGA_NUM = '',I3 )')&
 & LSOMEGA_FRC, NLSOMEGA_FREQ, NLSOMEGA_DEB, NLSOMEGA_NUM

WRITE(UNIT=KULOUT,FMT='(''  LT_NUDG = '',L5,&
 & '' NT_NUDG_DEB = '',I5,'' RELAX_TAUT = '',F8.1)')&
 & LT_NUDG, NT_NUDG_DEB, RELAX_TAUT

WRITE(UNIT=KULOUT,FMT='(''  LQV_NUDG = '',L5,&
 & '' NQV_NUDG_DEB = '',I5,'' RELAX_TAUQ = '',F8.1)')&
 & LQV_NUDG, NQV_NUDG_DEB, RELAX_TAUQ

WRITE(UNIT=KULOUT,FMT='(''  LUV_NUDG = '',L5,&
 & '' NU_NUDG_DEB = '',I5,'' NV_NUDG_DEB = '',I5,'' RELAX_TAUU = '',F8.1)')&
 & LUV_NUDG, NU_NUDG_DEB, NV_NUDG_DEB, RELAX_TAUU

WRITE(UNIT=KULOUT,FMT='(''  Forcage Surface Sens. Heat  '',&
 & '' NSH_FORC_DEB = '',I5,'' NSH_FORC_NUM = '',I5)')&
 & NSH_FORC_DEB, NSH_FORC_NUM

WRITE(UNIT=KULOUT,FMT='(''  Forcage Surface Lat. Heat  '',&
 & '' NLH_FORC_DEB = '',I5,'' NLH_FORC_NUM = '',I5)')&
 & NLH_FORC_DEB, NLH_FORC_NUM

WRITE(UNIT=KULOUT,FMT='(''  Ts  '',&
 & '' NTS_FORC_DEB = '',I5,'' NTS_FORC_NUM = '',I5)')&
 & NTS_FORC_DEB, NLH_FORC_NUM

WRITE(UNIT=KULOUT,FMT='(''  Ustar  '',&
 & '' NUS_FORC_DEB = '',I5,'' NUS_FORC_NUM = '',I5)')&
 & NUS_FORC_DEB, NUS_FORC_NUM

WRITE(UNIT=KULOUT,FMT='(''  RZ0_FORC  '', E14.7,&
 & '' RALB_FORC = '',E14.7,'' REMIS_FORC  = '',E14.7)')&
 & RZ0_FORC, RALB_FORC, REMIS_FORC

WRITE(UNIT=KULOUT,FMT='(''  LSPS_FRC = '',L5)')&
 & LSPS_FRC
! -----------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SULSFORC',1,ZHOOK_HANDLE)
END SUBROUTINE SULSFORC
