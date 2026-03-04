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

SUBROUTINE GPCTY_FORC(&
 ! --- INPUT -----------------------------------------------------------------
 & YDGEOMETRY, YDML_GCONF,YDPHY2,KST,KEND,&
 & PFORC,&
 & PRE0F, PTT0, PR0,&
 ! --- OUTPUT ----------------------------------------------------------------            
 & PEVEL, PVVEL )

!****GPCTY_FORC** -  Computes large scale vertical velocities for hybrid coordinate
!                    from large scale forcing read in initial file and stored in PFORC
!                    (a forcing given in term of w=Dz/Dt or omega=DP/dt)

!     Purpose.
!     --------
!      Computes the vertical velocities
!       "etadot (d prehyd / d eta)" on interlayers and
!       "(omega / prehyd)" on layers
!      from a given large scale w=Dz/Dt or omega=Dp/Dt prescribed as a forcing

!      These two terms are necessary to compute the large scale vertical advection
!      in 1D model (SCUM)

!      The computations of PVVEL and PEVEL in this routine 
!      overwrite the computations of GPCTY (which give PEVEL=PVVEL=0 in 1D)

!      * The large scale vertical velocity w/omega is computed from the values stored in PFORC
!         If only one large scale vertical velocity field is present in PFORC (NLSW/OMEGA_NUM=1)
!         this forcing is applied at each time step
!         If a different forcing is prescibed every NLSW/OMEGA_FREQ (in second), w/OMEGA at time t
!         is computed by a linear interpolation between the two prescibed forcings 
!         surrounding time t

!      * PVVEL is computed at full level with the hydrostatic simplified formulation
!         omega/p = - (g w) / (R T) or directly as omega/p if LSOMEGA_FRC
!         ZEVEL is computed at full levels as PVVEL* pressure at full levels (PRE0F)
!         ZEVEL is interpolated by a simple averaging from full to half levels (==> PEVEL)

!      Restrictions of use: the code is currently not available for:
!      * LRUBC=T
!      * NDPSFI=1

!**   Interface.
!     ----------
!        *CALL* *GPCTY_FORC(...)

!        Explicit arguments :
!        --------------------
!         * INPUT:
!        KST       : first element of work
!        KEND      : last element of work
!        PFORC     : forcings (the all part of GFL dedicated for forcing)
!        PRE0F     : hydrostatic pressure "prehyd" on layers (full levels) at time t
!        PTT0      : temperature at time t (full levels)
!        PR0       : air constant "R" at t (full levels)
!         * OUTPUT:
!        PVVEL0    : "omega/prehyd" at t (full levels)
!        PEVEL0    : "etadot d prehyd / d eta" at t (half levels)
!                    not yet coded for LVERTFE=T (finite difference) 

!        Implicit arguments :   None.
!        --------------------

!     Method.
!     -------
!        See SCUM (Single Column Unified Model) documentation

!     Externals.    None.
!     ----------

!     Reference.
!     ----------
!        See  SCUM documentation

!     Author.
!     -------
!        Sylvie Malardel

!     Modifications.
!     --------------
!        Original : 06-06-06
!     ------------------------------------------------------------------

USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE PARKIND1               , ONLY : JPIM, JPRB
USE YOMHOOK                , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST                 , ONLY : RG
USE YOMPHY2                , ONLY : TPHY2
USE YOMLSFORC              , ONLY : LSW_FRC, NLSW_NUM, NLSW_DEB, NLSW_FREQ, NLSW_TIME,&
 &                                  LSOMEGA_FRC, NLSOMEGA_NUM, NLSOMEGA_DEB, NLSOMEGA_FREQ, NLSOMEGA_TIME
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(IN):: YDML_GCONF
TYPE(TPHY2)       ,INTENT(IN)    :: YDPHY2
INTEGER(KIND=JPIM),INTENT(IN)    :: KST 
INTEGER(KIND=JPIM),INTENT(IN)    :: KEND 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFORC(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YDML_GCONF%YGFL%NGFL_FORC)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRE0F(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR0(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PEVEL(YDGEOMETRY%YRDIM%NPROMA,0:YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVVEL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JROF, JLEV, JT     
INTEGER(KIND=JPIM) :: IINDECH, IINDA, IINDB                                     
REAL(KIND=JPRB)    :: ZEVEL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZFW(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG)
REAL(KIND=JPRB)    :: ZA, ZB, ZINT, ZSTATI

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GPCTY_FORC',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YGFL=>YDML_GCONF%YGFL,YDRIP=>YDML_GCONF%YRRIP)
ASSOCIATE(NGFL_FORC=>YGFL%NGFL_FORC, &
 & NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & RSTATI=>YDRIP%RSTATI, &
 & TSPHY=>YDPHY2%TSPHY)
!     ------------------------------------------------------------------

ZSTATI=RSTATI-TSPHY*0.5_JPRB
! Case of a vertical velocity prescribed in m/s
IF(LSW_FRC) THEN
!     ------------------------------------------------------------------

!*          Computation of large scale w (m/s) for time t
!              ---------------------------------------------

IF(NLSW_NUM == 1) THEN

  ! Case of constant forcing
  IINDA = NLSW_DEB
  DO JROF=KST,KEND
    DO JLEV=1,NFLEVG
      ZFW(JROF,JLEV) =  PFORC(JROF,JLEV,IINDA)
    ENDDO
  ENDDO

ELSEIF(NLSW_NUM > 1) THEN
   IF (NLSW_FREQ/=999) THEN
  ! forcing are known with a constant time interval
  ! from the beginning to the end of the simulation   
  ! linear interpolation between two forcing times
  IINDECH= INT(ZSTATI / NLSW_FREQ)
  IINDA = NLSW_DEB + IINDECH
  IINDB = NLSW_DEB + IINDECH+1
  DO JROF=KST,KEND
    DO JLEV=1,NFLEVG
      ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
       & / REAL(NLSW_FREQ,JPRB)
      ZB= PFORC(JROF,JLEV,IINDA)
      ZFW(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NLSW_FREQ,JPRB))+ZB
    ENDDO
  ENDDO

    ELSEIF (NLSW_TIME(1)/=999) THEN
    ! forcing are known at individual times
    ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
    DO JT=2,NLSW_NUM
!!! Compute tendency for ZSTATI (current time) in between 
!!! NLSW_TIME(JT-1) and NLSW_TIME(JT)
      IF (  NLSW_TIME(JT)>=INT(ZSTATI)) THEN
      IINDA = NLSW_DEB + JT-2
      IINDB = NLSW_DEB + JT-1
      ZINT=REAL(NLSW_TIME(JT),JPRB)-REAL(NLSW_TIME(JT-1),JPRB)
      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
           &   / ZINT
          ZB = PFORC(JROF,JLEV,IINDA)
          ZFW(JROF,JLEV) = ZA*(ZSTATI-REAL(NLSW_TIME(JT-1),JPRB))+ZB
        ENDDO
      ENDDO
      GOTO 9998
      ENDIF
    ENDDO
!!!
    ELSE
        CALL ABOR1('Problem in the time definition of the forcing LSW')
     ENDIF

  ENDIF

!     ------------------------------------------------------------------

!*           Computation of large scale omega/p and m*eta_dot

! transformation of the w(m/s) given by the large scale forcing 
! into omega/p and m*eta_dot at full levels (simplified hydrostatic formulation)
9998 CONTINUE
DO JROF=KST,KEND
  DO JLEV=1,NFLEVG
    PVVEL(JROF,JLEV) = - ZFW(JROF,JLEV) * RG&
     & / ( PR0(JROF,JLEV) * PTT0(JROF,JLEV) )
    ZEVEL(JROF,JLEV) = - ZFW(JROF,JLEV) * RG * PRE0F(JROF,JLEV)&
     & / ( PR0(JROF,JLEV) * PTT0(JROF,JLEV) )
  ENDDO
ENDDO

! imposed values of m*eta_dot (PEVEL) at the surface and at the model top
! (the surface condition may be questionned for deltam=1)
! Assumes that LRUBC=F and NDPSFI=0.

DO JROF=KST,KEND
  PEVEL(JROF,NFLEVG) = 0._JPRB
  PEVEL(JROF,0) = 0._JPRB
ENDDO

IF (YDGEOMETRY%YRCVER%LVERTFE) THEN
  ! * PEVEL at full levels, simply copy ZEVEL into PEVEL.
  PEVEL(KST:KEND,1:NFLEVG)=ZEVEL(KST:KEND,1:NFLEVG) 
ELSE
  ! * PEVEL at half levels, for the time being a simple average between
  !   ZEVEL at full levels is done but this is probably not fully consistent
  !   with some other vertical discretisation, a cleaner solution must be
  !   studied in the future. 
  DO JROF=KST,KEND
    DO JLEV=2,NFLEVG
      PEVEL(JROF,JLEV-1) = 0.5_JPRB*( ZEVEL(JROF,JLEV)+ZEVEL(JROF,JLEV-1) )
    ENDDO
  ENDDO
ENDIF 

ENDIF

! Case of a vertical velocity prescribed in Pa/s
IF(LSOMEGA_FRC) THEN
!     ------------------------------------------------------------------

!*          Computation of large scale w (m/s) for time t
!              ---------------------------------------------

IF(NLSOMEGA_NUM == 1) THEN

  ! Case of constant forcing
  IINDA = NLSOMEGA_DEB
  DO JROF=KST,KEND
    DO JLEV=1,NFLEVG
      ZFW(JROF,JLEV) =  PFORC(JROF,JLEV,IINDA)
    ENDDO
  ENDDO

ELSEIF(NLSOMEGA_NUM > 1) THEN

   IF (NLSOMEGA_FREQ/=999) THEN
  ! forcing are known with a constant time interval
  ! from the beginning to the end of the simulation   
  ! linear interpolation between two forcing times
  IINDECH= INT(ZSTATI / NLSOMEGA_FREQ)
  IINDA = NLSOMEGA_DEB + IINDECH
  IINDB = NLSOMEGA_DEB + IINDECH+1
  DO JROF=KST,KEND
    DO JLEV=1,NFLEVG
      ZA= ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
       & / REAL(NLSOMEGA_FREQ,JPRB)
      ZB= PFORC(JROF,JLEV,IINDA)
      ZFW(JROF,JLEV) = ZA*(ZSTATI-REAL(IINDECH,JPRB)*REAL(NLSOMEGA_FREQ,JPRB))+ZB
    ENDDO
  ENDDO

    ELSEIF (NLSOMEGA_TIME(1)/=999) THEN
    ! forcing are known at individual times
    ! time interpolation of the two forcings surrounding time t
      IINDA = 1
      IINDB = 1
    DO JT=2,NLSOMEGA_NUM
!!! Compute tendency for ZSTATI (current time) in between 
!!! NLSOMEGA_TIME(JT-1) and NLSOMEGA_TIME(JT)
      IF (  NLSOMEGA_TIME(JT)>=INT(ZSTATI)) THEN
      IINDA = NLSOMEGA_DEB + JT-2
      IINDB = NLSOMEGA_DEB + JT-1
      ZINT=REAL(NLSOMEGA_TIME(JT),JPRB)-REAL(NLSOMEGA_TIME(JT-1),JPRB)
      DO JROF=KST,KEND
        DO JLEV=1,NFLEVG
          ZA = ( PFORC(JROF,JLEV,IINDB)-PFORC(JROF,JLEV,IINDA) )&
           &   / ZINT
          ZB = PFORC(JROF,JLEV,IINDA)
          ZFW(JROF,JLEV) = ZA*(ZSTATI-REAL(NLSOMEGA_TIME(JT-1),JPRB))+ZB
        ENDDO
      ENDDO
      GOTO 9999
      ENDIF
    ENDDO
!!!
    ELSE
        CALL ABOR1('Problem in the time definition of the forcing LSOMEGA')
    ENDIF 

ENDIF

!     ------------------------------------------------------------------

!*           Computation of large scale omega/p and m*eta_dot

! computation of omega/p and m*eta_dot=omega at full levels (simplified hydrostatic formulation)

9999 CONTINUE

DO JROF=KST,KEND
  DO JLEV=1,NFLEVG
    PVVEL(JROF,JLEV) = ZFW(JROF,JLEV)/ PRE0F(JROF,JLEV)
    ZEVEL(JROF,JLEV) = ZFW(JROF,JLEV)
  ENDDO
ENDDO

! imposed values of m*eta_dot (PEVEL) at the surface and at the model top
! (the surface condition may be questionned for deltam=1)
! Assumes that LRUBC=F and NDPSFI=0.

DO JROF=KST,KEND
  PEVEL(JROF,NFLEVG) = 0._JPRB
  PEVEL(JROF,0) = 0._JPRB
ENDDO

IF (YDGEOMETRY%YRCVER%LVERTFE) THEN
  ! * PEVEL at full levels, simply copy ZEVEL into PEVEL.
  PEVEL(KST:KEND,1:NFLEVG)=ZEVEL(KST:KEND,1:NFLEVG) 
ELSE
  ! * PEVEL at half levels, for the time being a simple average between
  !   ZEVEL at full levels is done but this is probably not fully consistent
  !   with some other vertical discretisation, a cleaner solution must be
  !   studied in the future. 
  DO JROF=KST,KEND
    DO JLEV=2,NFLEVG
      PEVEL(JROF,JLEV-1) = 0.5_JPRB*( ZEVEL(JROF,JLEV)+ZEVEL(JROF,JLEV-1) )
    ENDDO
  ENDDO
ENDIF 

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('GPCTY_FORC',1,ZHOOK_HANDLE)
END SUBROUTINE GPCTY_FORC
