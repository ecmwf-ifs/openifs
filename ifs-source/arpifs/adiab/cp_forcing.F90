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

SUBROUTINE CP_FORCING(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY, YDMODEL,KST,KEND,PDT,PAPRS,PAPRSF,PFORC,&
 & PUT0,PVT0,PTT0,PQT0,&
 ! --- INOUT -----------------------------------------------------------------
 & PATND,PATND_Q,PQT1,YDDDH )

!****CP_FORCING** -  Computes large scale forcing at time t from PFORC
!                    and force u,v,T and/or q at t+1

!     Purpose.
!     --------
!    Add large scale tendencies computed from prescribed 
!    large scale forcing to variables at time t+dt

!    This routine is validated only for SCUM (Single Column Unified Model)

!    The forcing are always computed at time t
!    (i.e. de-centred scheme in 2-TL and centred scheme in 3-TL ;
!    if predictor/corrector is used, the type of scheme change 
!    between predictor and corrector )

!    For each type of forcing (type XXX),
!    if only one forcing is present in PFORC (NXXX_NUM=1),
!    this forcing is applied at each time step
!    If a different forcing is prescibed every NXXX_FREQ or at
!    individual time in NXXX_TIME(NXXX_NUM) (in second), 
!    the forcing at time t is computed
!    by a linear interpolation between the two prescibed forcings 
!    surrounding time t

! For the time being, 3 types of forcing are available 
! (but more may be included easily later ...)

!  Geostrophic forcing (the presciption of a geostrophic wind and
!    a Coriolis parameter is equivalent 
!    to the presciption of a large scale pressure gradient force)

!  Large scale advection of wind

!  Large scale advection of temperature

!  Large scale advection of specific humidity

!**   Interface.
!     ----------
!        *CALL* *CP_FORCING(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!        KST       : first element of work
!        KEND      : last element of work
!        PDT       : For a leap-frog scheme (three time level scheme):
!                     'dt' at the first time-step, '2 dt' otherwise.
!                    For a 2TL SL scheme: timestep 'dt'.
!        PAPRS     : half-level pressure.
!        PFORC     : forcings (the all part of GFL dedicated for forcing)
!        PUT0      : zonal wind at time t
!        PVT0      : meridional wind at time t

!         * INOUT:
!        PATND     : adiabatic Lagrangian tendencies for GMV.
!        PATND_Q   : adiabatic Lagrangian tendencies for "q".
!        PQT1      : specific humidity at time t+dt
!        YDDDH     : diagnostic superstructure

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See SCUM (Single Column Unified Model) documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See SCUM documentation

!     Author.
!     -------
!        Sylvie Malardel

!     Modifications.
!     --------------
!        Original : 06-06-06
!        E. Bazile and I. Beau 2011-01-18 : Nudging forcing.
!        2011-09-07, Karim Yessad: call to CP_FORCING moved from CPG_DYN to CPG_GP.
!        2011-09-07, J.M. Piriou: add DDH diagnostics of 1D model MUSC dynamical tendencies.
!        2012-10-15, E. Bazile: Bug correction (found by R. Roehrig) Qv tendency.
!        2016-03-16, E. Bazile: variable nudging
!        K. Yessad (Feb 2018): remove deep-layer formulations.
!        2018-09-19, R. Roehrig: add option to read forcing in ASCII files
!                               + add MUSC diagnostics (potential temperature tendencies)
!     ------------------------------------------------------------------

USE TYPE_MODEL   , ONLY : MODEL
USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM, JPRB
USE YOMHOOK      , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLSFORC
USE DDH_MIX      , ONLY : ADD_FIELD_3D, NEW_ADD_FIELD_3D, TYP_DDH
USE YOMCST       , ONLY : RG, RD, RCPD
USE YOMCT3       , ONLY : NSTEP

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL)       ,INTENT(INOUT) :: YDMODEL
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRS(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPRSF(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_GCONF%YGFL%NGFL_FORC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PATND(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDMODEL%YRML_DYN%YRDYNA%YYTTND%NDIM)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PATND_Q(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PQT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
TYPE(TYP_DDH)     ,INTENT(INOUT) :: YDDDH
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV, JT, JC
INTEGER(KIND=JPIM) :: IINDECH, IINDA, IINDB                              
REAL(KIND=JPRB)    :: ZA, ZB, ZINT, ZSTATI
REAL(KIND=JPRB)    :: ZFU(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZFV(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZFT(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZFQ(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZFTH(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZTMPVAR(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZDUDTDYN(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG),ZDVDTDYN(YDGEOMETRY%YRDIM%NPROMA, &
 & YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZZFU(YDGEOMETRY%YRDIMV%NFLEVG),ZZFV(YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZZFT(YDGEOMETRY%YRDIMV%NFLEVG),ZZFQ(YDGEOMETRY%YRDIMV%NFLEVG)
CHARACTER          :: CLSTEP*20
CHARACTER          :: CLFILE*200
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "wrscmr.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CP_FORCING',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, &
 &  YDTOPH=>YDMODEL%YRML_PHY_MF%YRTOPH, &
  & YDRIP=>YDMODEL%YRML_GCONF%YRRIP,YDLDDH=>YDMODEL%YRML_DIAG%YRLDDH,YDPHY2=>YDMODEL%YRML_PHY_MF%YRPHY2, &
  & YDDYNA=>YDMODEL%YRML_DYN%YRDYNA,YGFL=>YDMODEL%YRML_GCONF%YGFL)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & TSPHY=>YDPHY2%TSPHY, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NTRELAXQ=>YDTOPH%NTRELAXQ, NTRELAXT=>YDTOPH%NTRELAXT, &
 & NTRELAXU=>YDTOPH%NTRELAXU, &
 & LFLEXDIA=>YDLDDH%LFLEXDIA, LDDH_OMP=>YDLDDH%LDDH_OMP, &
 & RSTATI=>YDRIP%RSTATI,YYTTND=>YDDYNA%YYTTND,  &
 & NGFL_FORC=>YGFL%NGFL_FORC)
!     ------------------------------------------------------------------

!          0.   PRELIMINARY INITIALISATIONS
!               ---------------------------

ZSTATI=RSTATI-TSPHY*0.5_JPRB
PATND_Q(:,:)=0.0_JPRB

!     ------------------------------------------------------------------

!          1.   GEOSTROPHIC FORCING
!               -------------------

IF(LGEOST_UV_FRC) THEN

  IF (NGEOST_U_NUM == -1) THEN
    ! Reading forcings in the appropriate files        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/ug_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFU(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFU(JROF,JLEV) = ZZFU(JLEV)
        ZTMPVAR(JROF,JLEV)=-RCORIO_FORC*(PUT0(JROF,JLEV)-ZFU(JROF,JLEV))
        ZDVDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NGEOST_U_NUM == 1) THEN
    ! Constant forcing
    ! Zonal componant
    IINDA = NGEOST_U_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFU(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        ZTMPVAR(JROF,JLEV)=-RCORIO_FORC*(PUT0(JROF,JLEV)-ZFU(JROF,JLEV))
        ZDVDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NGEOST_U_NUM > 1) THEN

    IF (NGEOST_UV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDECH= INT(ZSTATI / NGEOST_UV_FREQ)
      IINDA = NGEOST_U_DEB + IINDECH
      IINDB = NGEOST_U_DEB + IINDECH+1

      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NGEOST_UV_FREQ,JPRB)
          ZB = PFORC(JROF,JLEV,IINDA)
          ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NGEOST_UV_FREQ,JPRB))+ZB
          ZTMPVAR(JROF,JLEV)=-RCORIO_FORC*(PUT0(JROF,JLEV)-ZFU(JROF,JLEV))
          ZDVDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NGEOST_UV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NGEOST_U_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NGEOST_UV_TIME(JT-1) and NGEOST_UV_TIME(JT)
        IF (  NGEOST_UV_TIME(JT)>=INT(ZSTATI)) THEN
          IINDA = NGEOST_U_DEB + JT-2
          IINDB = NGEOST_U_DEB + JT-1
          ZINT=REAL(NGEOST_UV_TIME(JT),JPRB)-REAL(NGEOST_UV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(NGEOST_UV_TIME(JT-1),JPRB))+ZB
              ZTMPVAR(JROF,JLEV)=-RCORIO_FORC*(PUT0(JROF,JLEV)-ZFU(JROF,JLEV))
              ZDVDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO

    ELSE
      CALL ABOR1('Problem in the time definition of the forcing GEOST_UV')
    ENDIF

  ENDIF

  IF (NGEOST_V_NUM == -1) THEN
    ! Reading forcings in the appropriate files        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/vg_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFV(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFV(JROF,JLEV) = ZZFV(JLEV)
        ZTMPVAR(JROF,JLEV)=RCORIO_FORC*(PVT0(JROF,JLEV)-ZFV(JROF,JLEV))
        ZDUDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NGEOST_V_NUM == 1) THEN
    ! Constant forcing
    ! Meridional componant
    IINDA = NGEOST_V_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFV(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        ZTMPVAR(JROF,JLEV)=RCORIO_FORC*(PVT0(JROF,JLEV)-ZFV(JROF,JLEV))
        ZDUDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NGEOST_V_NUM > 1) THEN

    IF (NGEOST_UV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDA = NGEOST_V_DEB + IINDECH
      IINDB = NGEOST_V_DEB + IINDECH+1
      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NGEOST_UV_FREQ,JPRB)
          ZB = PFORC(JROF,JLEV,IINDA)
          ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NGEOST_UV_FREQ,JPRB))+ZB
          ZTMPVAR(JROF,JLEV)=RCORIO_FORC*(PVT0(JROF,JLEV)-ZFV(JROF,JLEV))
          ZDUDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NGEOST_UV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NGEOST_U_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NGEOST_UV_TIME(JT-1) and NGEOST_UV_TIME(JT)
        IF (  NGEOST_UV_TIME(JT)>=INT(ZSTATI)) THEN
          !      PRINT*,'JT=',JT
          !      PRINT*,'NGEOST_UV_TIME(JT)=',NGEOST_UV_TIME(JT)
          IINDA = NGEOST_V_DEB + JT-2
          IINDB = NGEOST_V_DEB + JT-1
          ZINT=REAL(NGEOST_UV_TIME(JT),JPRB)-REAL(NGEOST_UV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(NGEOST_UV_TIME(JT-1),JPRB))+ZB
              ZTMPVAR(JROF,JLEV)=RCORIO_FORC*(PVT0(JROF,JLEV)-ZFV(JROF,JLEV))
              ZDUDTDYN(JROF,JLEV)=ZTMPVAR(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO

    ELSE
      CALL ABOR1('Problem in the time definition of the forcing GEOST_UV')
    ENDIF

  ENDIF

  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZFUGEO',ZFU,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZDUGEO',ZDUDTDYN,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZFVGEO',ZFV,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZDVGEO',ZDVDTDYN,NPROMA,NFLEVG)

ENDIF ! LGEOST_UV_FRC

!     ------------------------------------------------------------------

!          2.   LARGE SCALE ADVECTION
!               ---------------------

!*  Large scale advection of wind
!   -----------------------------

IF(LUV_ADV_FRC) THEN

  IF (NU_ADV_NUM == -1) THEN
    ! Reading forcings in the appropriate files        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/du_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFU(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFU(JROF,JLEV) = ZZFU(JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZFU(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZFU(JROF,JLEV) 
      ENDDO
    ENDDO      
  ELSEIF (NU_ADV_NUM == 1) THEN
    ! constant forcing
    IINDA = NU_ADV_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFU(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZFU(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZFU(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NU_ADV_NUM > 1) THEN 

    IF (NUV_ADV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDECH= INT(ZSTATI / NUV_ADV_FREQ)
      IINDA = NU_ADV_DEB + IINDECH
      IINDB = NU_ADV_DEB + IINDECH+1

      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NUV_ADV_FREQ,JPRB)
          ZB= PFORC(JROF,JLEV,IINDA)
          ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NUV_ADV_FREQ,JPRB))+ZB
          PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZFU(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZFU(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NUV_ADV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NU_ADV_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NUV_ADV_TIME(JT-1) and NUV_ADV_TIME(JT)
        IF (  NUV_ADV_TIME(JT)>=INT(ZSTATI)) THEN
          IINDA = NU_ADV_DEB + JT-2
          IINDB = NU_ADV_DEB + JT-1
          ZINT=REAL(NUV_ADV_TIME(JT),JPRB)-REAL(NUV_ADV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(NUV_ADV_TIME(JT-1),JPRB))+ZB
              PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZFU(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZFU(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO
    ELSE
      CALL ABOR1('Problem in the time definition of the forcing UV_ADV')
    ENDIF 

  ENDIF

  IF (NV_ADV_NUM == -1) THEN        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/dv_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFV(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFV(JROF,JLEV) = ZZFV(JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZFV(JROF,JLEV)    
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZFV(JROF,JLEV)
      ENDDO
    ENDDO      
  ELSEIF (NV_ADV_NUM == 1) THEN
    ! constant forcing
    IINDA = NV_ADV_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFV(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZFV(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZFV(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NV_ADV_NUM > 1) THEN 

    IF (NUV_ADV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDECH= INT(ZSTATI / NUV_ADV_FREQ)
      IINDA = NV_ADV_DEB + IINDECH
      IINDB = NV_ADV_DEB + IINDECH+1

      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NUV_ADV_FREQ,JPRB)
          ZB= PFORC(JROF,JLEV,IINDA)
          ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NUV_ADV_FREQ,JPRB))+ZB
          PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZFV(JROF,JLEV)
          PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZFV(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NUV_ADV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NU_ADV_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NUV_ADV_TIME(JT-1) and NUV_ADV_TIME(JT)
        IF (  NUV_ADV_TIME(JT)>=INT(ZSTATI)) THEN
          IINDA = NV_ADV_DEB + JT-2
          IINDB = NV_ADV_DEB + JT-1
          ZINT=REAL(NUV_ADV_TIME(JT),JPRB)-REAL(NUV_ADV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(NUV_ADV_TIME(JT-1),JPRB))+ZB
              PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZFV(JROF,JLEV)
              PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZFV(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO

    ELSE
      CALL ABOR1('Problem in the time definition of the forcing UV_ADV')
    ENDIF 

  ENDIF

  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZFU',ZFU,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZFV',ZFV,NPROMA,NFLEVG)

ENDIF ! LUV_ADV_FRC

!*  Large scale advection of temperature
!   ------------------------------------

IF(LT_ADV_FRC) THEN

  IF (NT_ADV_NUM == -1) THEN        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/dT_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFT(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFT(JROF,JLEV) = ZZFT(JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)+ZFT(JROF,JLEV)   
      ENDDO
    ENDDO      
  ELSEIF (NT_ADV_NUM == 1) THEN
    ! constant forcing
    IINDA = NT_ADV_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFT(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)+ZFT(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NT_ADV_NUM > 1) THEN 

    IF (NT_ADV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDECH= INT(ZSTATI / NT_ADV_FREQ)
      IINDA = NT_ADV_DEB + IINDECH
      IINDB = NT_ADV_DEB + IINDECH+1

      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NT_ADV_FREQ,JPRB)
          ZB= PFORC(JROF,JLEV,IINDA)
          ZFT(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NT_ADV_FREQ,JPRB))+ZB
          PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)+ZFT(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NT_ADV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NT_ADV_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NT_ADV_TIME(JT-1) and NT_ADV_TIME(JT)
        !      PRINT*,'JT=',JT
        !      PRINT*,'NT_ADV_TIME(JT)=',NT_ADV_TIME(JT)
        !      PRINT*,'ZSTATI=',ZSTATI
        IF (  NT_ADV_TIME(JT)>INT(ZSTATI)) THEN
          IINDA = NT_ADV_DEB + JT-2
          IINDB = NT_ADV_DEB + JT-1
          !      PRINT*,'IINDA=',IINDA 
          !      PRINT*,'IINDB=',IINDB
          ZINT=REAL(NT_ADV_TIME(JT),JPRB)-REAL(NT_ADV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFT(JROF,JLEV) = ZA*(ZSTATI-REAL(NT_ADV_TIME(JT-1),JPRB))+ZB
              PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)+ZFT(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO

    ELSE
      CALL ABOR1('Problem in the time definition of the forcing T_ADV')
    ENDIF 
  ENDIF

  IF (LMUSCLFA) THEN
    CALL WRSCMR(NMUSCLFA,'ZFT',ZFT,NPROMA,NFLEVG)
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFTH(JROF,JLEV) = ZFT(JROF,JLEV)*(PAPRSF(JROF,JLEV)/100000._JPRB)**(RD/RCPD)
      ENDDO
    ENDDO
    CALL WRSCMR(NMUSCLFA,'ZFTH',ZFTH,NPROMA,NFLEVG)
  ENDIF

ENDIF ! LT_ADV_FRC

!*  Large scale advection of specific humidity
!   ------------------------------------------

IF(LQV_ADV_FRC) THEN

  IF (NQV_ADV_NUM == -1) THEN        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/dq_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFQ(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFQ(JROF,JLEV) =  ZZFQ(JLEV)
        PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZFQ(JROF,JLEV)   
      ENDDO
    ENDDO      
  ELSEIF (NQV_ADV_NUM == 1) THEN
    ! Constant forcing
    IINDA = NQV_ADV_DEB
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFQ(JROF,JLEV) = PFORC(JROF,JLEV,IINDA)
        PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZFQ(JROF,JLEV)
      ENDDO
    ENDDO

  ELSEIF (NQV_ADV_NUM > 1) THEN

    IF (NQV_ADV_FREQ/=999) THEN
      ! forcing are known with a constant time interval
      ! from the beginning to the end of the simulation
      ! time interpolation of the two forcings surrounding time t
      IINDECH= INT(ZSTATI / NQV_ADV_FREQ)
      IINDA = NQV_ADV_DEB + IINDECH
      IINDB = NQV_ADV_DEB + IINDECH+1

      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / REAL(NQV_ADV_FREQ,JPRB)
          ZB= PFORC(JROF,JLEV,IINDA)
          ZFQ(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NQV_ADV_FREQ,JPRB))+ZB
          PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZFQ(JROF,JLEV)
        ENDDO
      ENDDO

    ELSEIF (NQV_ADV_TIME(1)/=999) THEN
      ! forcing are known at individual times
      ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
      DO JT=2,NQV_ADV_NUM
        !!! Compute tendency for ZSTATI (current time) in between 
        !!! NQV_ADV_TIME(JT-1) and NQV_ADV_TIME(JT)
        !      PRINT*,'JT=',JT
        !      PRINT*,'NQV_ADV_TIME(JT)=',NQV_ADV_TIME(JT)
        !      PRINT*,'ZSTATI=',ZSTATI
        IF (  NQV_ADV_TIME(JT)>=INT(ZSTATI)) THEN
          IINDA = NQV_ADV_DEB + JT-2
          IINDB = NQV_ADV_DEB + JT-1
          !      PRINT*,'IINDA=',IINDA 
          !      PRINT*,'IINDB=',IINDB
          ZINT=REAL(NQV_ADV_TIME(JT),JPRB)-REAL(NQV_ADV_TIME(JT-1),JPRB)
          DO JROF=KST,KEND
            DO JLEV=1,NFLEVG
              ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) ) / ZINT
              ZB = PFORC(JROF,JLEV,IINDA)
              ZFQ(JROF,JLEV) = ZA*(ZSTATI-REAL(NQV_ADV_TIME(JT-1),JPRB))+ZB
              PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZFQ(JROF,JLEV)
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO

    ELSE
      CALL ABOR1('Problem in the time definition of the forcing QV_ADV')
    ENDIF

  ENDIF

  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'ZFQ',ZFQ,NPROMA,NFLEVG)

ENDIF ! LQV_ADV_FRC

!     ------------------------------------------------------------------

!          3.   NUDGING PART
!               ------------

! GMV:

IF(LUV_NUDG) THEN
  IF ((NU_NUDG == -1).AND.(NV_NUDG == -1)) THEN        
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/u_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFU(JLEV)
    ENDDO
    CLOSE(82)

    WRITE(CLFILE,FMT='(3A)') 'files/v_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFV(JLEV)
    ENDDO
    CLOSE(82)

    ZTMPVAR(:,:)=0._JPRB
    ZFU(:,:)=0._JPRB
    DO JROF=KST,KEND
      DO JLEV=1,NTRELAXU
        ZFU(JROF,JLEV) = ZZFU(JLEV)
        ZTMPVAR(JROF,JLEV)=(ZFU(JROF,JLEV)-PUT0(JROF,JLEV))/RELAX_TAUU
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'U_NUDG',ZFU,NPROMA,NFLEVG)
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TU_NUDG',ZTMPVAR,NPROMA,NFLEVG)

    ZTMPVAR(:,:)=0._JPRB
    ZFV(:,:)=0._JPRB
    DO JROF=KST,KEND
      DO JLEV=1,NTRELAXU
        ZFV(JROF,JLEV) = ZZFV(JLEV)
        ZTMPVAR(JROF,JLEV)=(ZFV(JROF,JLEV)-PVT0(JROF,JLEV))/RELAX_TAUU
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO

    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'V_NUDG',ZFV,NPROMA,NFLEVG)
    IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TV_NUDG',ZTMPVAR,NPROMA,NFLEVG)

  ELSEIF ( (NUV_NUDG_NUM == 1 ) .AND. (NU_NUDG_DEB /= 0) ) THEN
      ZTMPVAR(:,:)=0._JPRB
      ZFU(:,:)=0._JPRB
      DO JROF=KST,KEND
      DO JLEV=1,NTRELAXU
        ZFU(JROF,JLEV)=PFORC(JROF,JLEV,NU_NUDG_DEB)
        ZTMPVAR(JROF,JLEV)=(ZFU(JROF,JLEV)-PUT0(JROF,JLEV))/RELAX_TAUU
        PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
      ENDDO
      IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'U_NUDG',ZFU,NPROMA,NFLEVG)
      IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TU_NUDG',ZTMPVAR,NPROMA,NFLEVG)
      ZTMPVAR(:,:)=0._JPRB
      ZFV(:,:)=0._JPRB
      DO JROF=KST,KEND
      DO JLEV=1,NTRELAXU
        ZFV(JROF,JLEV)=PFORC(JROF,JLEV,NV_NUDG_DEB)
        ZTMPVAR(JROF,JLEV)=(ZFV(JROF,JLEV)-PVT0(JROF,JLEV))/RELAX_TAUU
        PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
        PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
      ENDDO
      ENDDO
      IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'V_NUDG',ZFV,NPROMA,NFLEVG)
      IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TV_NUDG',ZTMPVAR,NPROMA,NFLEVG)
  ELSEIF (NUV_NUDG_NUM > 1) THEN
      IF (NUV_NUDG_FREQ/=999) THEN
          ! forcing are known with a constant time interval
          ! from the beginning to the end of the simulation
          ! time interpolation of the two forcings surrounding time t
          IINDECH= INT(ZSTATI / NUV_NUDG_FREQ)
          IINDA = NU_NUDG_DEB + IINDECH
          IINDB = MIN(NU_NUDG_DEB + IINDECH+1,NU_NUDG_DEB+NUV_NUDG_NUM-1)

          ZTMPVAR(:,:)=0._JPRB
          ZFU(:,:)=0._JPRB
          DO JROF=KST,KEND
          DO JLEV=1,NTRELAXU
             ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
               &   / REAL(NUV_NUDG_FREQ,JPRB)
             ZB= PFORC(JROF,JLEV,IINDA)
             ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NUV_NUDG_FREQ,JPRB))+ZB
             ZTMPVAR(JROF,JLEV)=(ZFU(JROF,JLEV)-PUT0(JROF,JLEV))/RELAX_TAUU
             PATND(JROF,JLEV,YYTTND%M_TNDU)=PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
             PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
          ENDDO
          ENDDO
          IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'U_NUDG',ZFU,NPROMA,NFLEVG)
          IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TU_NUDG',ZTMPVAR,NPROMA,NFLEVG)

          IINDA = NV_NUDG_DEB + IINDECH
          IINDB = MIN(NV_NUDG_DEB + IINDECH+1,NV_NUDG_DEB+NUV_NUDG_NUM-1)
          ZTMPVAR(:,:)=0._JPRB
          ZFV(:,:)=0._JPRB
          DO JROF=KST,KEND
          DO JLEV=1,NTRELAXU
             ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
               &   / REAL(NUV_NUDG_FREQ,JPRB)
             ZB= PFORC(JROF,JLEV,IINDA)
             ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NUV_NUDG_FREQ,JPRB))+ZB
             ZTMPVAR(JROF,JLEV)=(ZFV(JROF,JLEV)-PVT0(JROF,JLEV))/RELAX_TAUU
             PATND(JROF,JLEV,YYTTND%M_TNDV)=PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
             PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
          ENDDO
          ENDDO
          IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'V_NUDG',ZFV,NPROMA,NFLEVG)
          IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TV_NUDG',ZTMPVAR,NPROMA,NFLEVG)
      ELSEIF (NUV_NUDG_TIME(1)/=999) THEN
          IINDA = 1
          IINDB = 1
          DO JT=2,NUV_NUDG_NUM
             IF (  NUV_NUDG_TIME(JT)>INT(ZSTATI)) THEN
                IINDA = NU_NUDG_DEB + JT-2
                IINDB = NU_NUDG_DEB + JT-1
                ZINT=REAL(NUV_NUDG_TIME(JT),JPRB)-REAL(NUV_NUDG_TIME(JT-1),JPRB)
                ZTMPVAR(:,:)=0._JPRB
                ZFU(:,:)=0._JPRB
                DO JROF=KST,KEND
                DO JLEV=1,NTRELAXU
                   ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
                    &   / ZINT
                   ZB= PFORC(JROF,JLEV,IINDA)
                   ZFU(JROF,JLEV) = ZA*(ZSTATI-REAL(NUV_NUDG_TIME(JT-1),JPRB))+ZB
                   ZTMPVAR(JROF,JLEV)=(ZFU(JROF,JLEV)-PUT0(JROF,JLEV))/RELAX_TAUU
                   PATND(JROF,JLEV,YYTTND%M_TNDU)=&
                          & PATND(JROF,JLEV,YYTTND%M_TNDU)+ZTMPVAR(JROF,JLEV)
                   PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)=&
                          & PATND(JROF,JLEV,YYTTND%M_TNDU_NOC)+ZTMPVAR(JROF,JLEV)
                ENDDO
                ENDDO
                IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'U_NUDG',ZFU,NPROMA,NFLEVG)
                IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TU_NUDG',ZTMPVAR,NPROMA,NFLEVG)

                IINDA = NV_NUDG_DEB + JT-2
                IINDB = NV_NUDG_DEB + JT-1
                ZTMPVAR(:,:)=0._JPRB
                ZFV(:,:)=0._JPRB
                DO JROF=KST,KEND
                DO JLEV=1,NTRELAXU
                   ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
                     &   / ZINT
                   ZB= PFORC(JROF,JLEV,IINDA)
                   ZFV(JROF,JLEV) = ZA*(ZSTATI-REAL(NUV_NUDG_TIME(JT-1),JPRB))+ZB
                   ZTMPVAR(JROF,JLEV)=(ZFV(JROF,JLEV)-PVT0(JROF,JLEV))/RELAX_TAUU
                   PATND(JROF,JLEV,YYTTND%M_TNDV)=&
                          & PATND(JROF,JLEV,YYTTND%M_TNDV)+ZTMPVAR(JROF,JLEV)
                   PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)=&
                          & PATND(JROF,JLEV,YYTTND%M_TNDV_NOC)+ZTMPVAR(JROF,JLEV)
                ENDDO
                ENDDO
                IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'V_NUDG',ZFV,NPROMA,NFLEVG)
                IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TV_NUDG',ZTMPVAR,NPROMA,NFLEVG)
             ENDIF  
          ENDDO
      ENDIF
  ELSE 
      CALL ABOR1('Problem in the time definition of the forcing WIND_NUDG')
  ENDIF
ENDIF

! TEMPERATURE NUDGING
IF(LT_NUDG) THEN
  ZTMPVAR(:,:)=0._JPRB
  ZFT(:,:)=0._JPRB
  IF (NT_NUDG == -1) THEN 
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/T_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFT(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
      DO JLEV=1,NTRELAXT
        ZFT(JROF,JLEV) = ZZFT(JLEV)
        ZTMPVAR(JROF,JLEV)=(ZFT(JROF,JLEV)-PTT0(JROF,JLEV))/RELAX_TAUT
        PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)&
        & +ZTMPVAR(JROF,JLEV)
      ENDDO
    ENDDO      
  ELSEIF (NT_NUDG_NUM==1) THEN
     ! Constant profile
     DO JROF=KST,KEND
     DO JLEV=1,NTRELAXT
        ZFT(JROF,JLEV)=PFORC(JROF,JLEV,NT_NUDG_DEB)
        ZTMPVAR(JROF,JLEV)=(ZFT(JROF,JLEV)-PTT0(JROF,JLEV))/RELAX_TAUT
        PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)&
       & +ZTMPVAR(JROF,JLEV)
     ENDDO
     ENDDO
  ELSEIF (NT_NUDG_NUM > 1) THEN
     IF (NT_NUDG_FREQ/=999) THEN
          ! forcing are known with a constant time interval
          ! from the beginning to the end of the simulation
          ! time interpolation of the two forcings surrounding time t
          IINDECH= INT(ZSTATI / NT_NUDG_FREQ)
          IINDA = NT_NUDG_DEB + IINDECH
          IINDB = NT_NUDG_DEB + IINDECH+1
          IINDB = MIN(NT_NUDG_DEB + IINDECH+1,NT_NUDG_DEB+NT_NUDG_NUM-1)

          DO JROF=KST,KEND
          DO JLEV=1,NTRELAXT
             ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
               &   / REAL(NT_NUDG_FREQ,JPRB)
             ZB= PFORC(JROF,JLEV,IINDA)
             ZFT(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NT_NUDG_FREQ,JPRB))+ZB
             ZTMPVAR(JROF,JLEV)=(ZFT(JROF,JLEV)-PTT0(JROF,JLEV))/RELAX_TAUT
             PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)&
                & +ZTMPVAR(JROF,JLEV)
          ENDDO
          ENDDO
      ELSEIF (NT_NUDG_TIME(1)/=999) THEN
          IINDA = 1
          IINDB = 1
          DO JT=2,NT_NUDG_NUM
             IF (  NT_NUDG_TIME(JT)>INT(ZSTATI)) THEN
                IINDA = NT_NUDG_DEB + JT-2
                IINDB = NT_NUDG_DEB + JT-1
                ZINT=REAL(NT_NUDG_TIME(JT),JPRB)-REAL(NT_NUDG_TIME(JT-1),JPRB)
                DO JROF=KST,KEND
                DO JLEV=1,NTRELAXT
                   ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
                    &   / ZINT
                   ZB= PFORC(JROF,JLEV,IINDA)
                   ZFT(JROF,JLEV) = ZA*(ZSTATI-REAL(NT_NUDG_TIME(JT-1),JPRB))+ZB
                   ZTMPVAR(JROF,JLEV)=(ZFT(JROF,JLEV)-PTT0(JROF,JLEV))/RELAX_TAUT
                   PATND(JROF,JLEV,YYTTND%M_TNDT)=PATND(JROF,JLEV,YYTTND%M_TNDT)&
                     & +ZTMPVAR(JROF,JLEV)
                ENDDO
                ENDDO
             ENDIF  
          ENDDO
       ELSE
          CALL ABOR1('Problem in the time definition of the forcing NT_NUDG_FREQ')
       ENDIF
    ELSE 
       CALL ABOR1('Problem in the time definition of the forcing NT_NUDG_NUM')
  ENDIF
  IF (LMUSCLFA) THEN
    CALL WRSCMR(NMUSCLFA,'T_NUDG',ZFT,NPROMA,NFLEVG)
    CALL WRSCMR(NMUSCLFA,'TT_NUDG',ZTMPVAR,NPROMA,NFLEVG)
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFTH(JROF,JLEV) = ZFT(JROF,JLEV)*(PAPRSF(JROF,JLEV)/100000._JPRB)**(RD/RCPD)
      ENDDO
    ENDDO
    CALL WRSCMR(NMUSCLFA,'TH_NUDG',ZFTH,NPROMA,NFLEVG)
    DO JROF=KST,KEND
      DO JLEV=1,NFLEVG
        ZFTH(JROF,JLEV) = ZTMPVAR(JROF,JLEV)*(PAPRSF(JROF,JLEV)/100000._JPRB)**(RD/RCPD)
      ENDDO
    ENDDO
    CALL WRSCMR(NMUSCLFA,'TTH_NUDG',ZFTH,NPROMA,NFLEVG)
  ENDIF
ENDIF

! GFL:

IF(LQV_NUDG) THEN
  ZTMPVAR(:,:)=0._JPRB
  ZFQ(:,:)=0._JPRB
  IF (NQV_NUDG == -1) THEN 
    WRITE(CLSTEP,FMT='(I5)') NSTEP
    DO JC=1,5
      IF (CLSTEP(JC:JC) == ' ') CLSTEP(JC:JC)='0'
    ENDDO
    WRITE(CLFILE,FMT='(3A)') 'files/q_forcing_',CLSTEP(1:len_trim(CLSTEP)),'.txt'

    OPEN(UNIT=82,FILE=CLFILE,FORM='formatted')
    DO JLEV=1,NFLEVG
      read(82,*) ZZFQ(JLEV)
    ENDDO
    CLOSE(82)

    DO JROF=KST,KEND
    DO JLEV=1,NTRELAXQ
        ZFQ(JROF,JLEV) = ZZFQ(JLEV)
        ZTMPVAR(JROF,JLEV)=(ZFQ(JROF,JLEV)-PQT0(JROF,JLEV))/RELAX_TAUQ
        PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZTMPVAR(JROF,JLEV)
    ENDDO
    ENDDO      
  ELSEIF (NQV_NUDG_NUM==1) THEN
     ! Constant profile
     DO JROF=KST,KEND
     DO JLEV=1,NTRELAXQ
        ZFQ(JROF,JLEV)=PFORC(JROF,JLEV,NQV_NUDG_DEB)
        ZTMPVAR(JROF,JLEV)=(ZFQ(JROF,JLEV)-PQT0(JROF,JLEV))/RELAX_TAUQ
        PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZTMPVAR(JROF,JLEV)
     ENDDO
     ENDDO
  ELSEIF (NQV_NUDG_NUM > 1) THEN
     IF (NQV_NUDG_FREQ/=999) THEN
          ! forcing are known with a constant time interval
          ! from the beginning to the end of the simulation
          ! time interpolation of the two forcings surrounding time t
          IINDECH= INT(ZSTATI / NQV_NUDG_FREQ)
          IINDA = NQV_NUDG_DEB + IINDECH
          IINDB = NQV_NUDG_DEB + IINDECH+1
          IINDB = MIN(NQV_NUDG_DEB + IINDECH+1,NQV_NUDG_DEB+NQV_NUDG_NUM-1)

          DO JROF=KST,KEND
          DO JLEV=1,NTRELAXQ
             ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
               &   / REAL(NQV_NUDG_FREQ,JPRB)
             ZB= PFORC(JROF,JLEV,IINDA)
             ZFQ(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NQV_NUDG_FREQ,JPRB))+ZB
             ZTMPVAR(JROF,JLEV)=(ZFQ(JROF,JLEV)-PQT0(JROF,JLEV))/RELAX_TAUQ
             PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZTMPVAR(JROF,JLEV)
          ENDDO
          ENDDO
      ELSEIF (NQV_NUDG_TIME(1)/=999) THEN
          IINDA = 1
          IINDB = 1
          DO JT=2,NQV_NUDG_NUM
             IF (  NQV_NUDG_TIME(JT)>INT(ZSTATI)) THEN
                IINDA = NQV_NUDG_DEB + JT-2
                IINDB = NQV_NUDG_DEB + JT-1
                ZINT=REAL(NQV_NUDG_TIME(JT),JPRB)-REAL(NQV_NUDG_TIME(JT-1),JPRB)
                DO JROF=KST,KEND
                DO JLEV=1,NTRELAXQ
                   ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
                    &   / ZINT
                   ZB= PFORC(JROF,JLEV,IINDA)
                   ZFQ(JROF,JLEV) = ZA*(ZSTATI-REAL(NQV_NUDG_TIME(JT-1),JPRB))+ZB
                   ZTMPVAR(JROF,JLEV)=(ZFQ(JROF,JLEV)-PQT0(JROF,JLEV))/RELAX_TAUQ
                   PATND_Q(JROF,JLEV)=PATND_Q(JROF,JLEV)+ZTMPVAR(JROF,JLEV)
                ENDDO
                ENDDO
             ENDIF  
          ENDDO
       ELSE
          CALL ABOR1('Problem in the time definition of the forcing NQ_NUDG_FREQ')
       ENDIF
    ELSE 
       CALL ABOR1('Problem in the time definition of the forcing NQ_NUDG_NUM')
  ENDIF
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'Q_NUDG',ZFQ,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TQ_NUDG',ZTMPVAR,NPROMA,NFLEVG)
  IF (LMUSCLFA) CALL WRSCMR(NMUSCLFA,'TQVRELAX',PATND_Q,NPROMA,NFLEVG)
  IF (LFLEXDIA) THEN
    ! DDH diagnostics.
    ZTMPVAR(:,:)=0._JPRB
    DO JROF=KST,KEND
      DO JLEV=1,NTRELAXQ
        ZTMPVAR(JROF,JLEV)=PATND_Q(JROF,JLEV)*(PAPRS(JROF,JLEV)-PAPRS(JROF,JLEV-1))/RG
      ENDDO
    ENDDO
    IF (LDDH_OMP) THEN
      CALL NEW_ADD_FIELD_3D(YDMODEL%YRML_DIAG%YRMDDH,ZTMPVAR,'TQVRELAX',YDDDH)
    ELSE
      CALL ADD_FIELD_3D(YDLDDH,ZTMPVAR,'TQVRELAX','T','ARP',.TRUE.,.TRUE.)
    ENDIF
  ENDIF
ENDIF

!     ------------------------------------------------------------------

!          4.   INCREMENT P[X]T1 FOR GFL
!               ------------------------

! For GMV variables, this incrementation is done under CPEULDYN or LACDYN.
! For GFL variables, this incrementation must be done there because CPEULDYN and LACDYN
!  assume that the increment is zero, so nothing is done in these two routines.

DO JROF=KST,KEND
  DO JLEV=1,NFLEVG
    PQT1(JROF,JLEV)=PQT1(JROF,JLEV)+PATND_Q(JROF,JLEV)*PDT
  ENDDO
ENDDO

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CP_FORCING',1,ZHOOK_HANDLE)
END SUBROUTINE CP_FORCING
