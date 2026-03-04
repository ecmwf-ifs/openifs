! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_INTERP8( KN_OUT, PX_OUT, PY_OUT, KN, PX_IN, PY_IN, &
!VH     &                   LD_OK,  CD_EXT_VALS, LD_SPLINE, LD_ON_LOG)
     &                   LD_OK,  CD_EXT_VALS)
!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!-----------------------------------------------------------------------
!  General interpolation routine.
!  BEWARE! This code is not optimized and should NOT be used inside
!          big/nested loops
!
!  kn_out is the dimension of x_out and y_out
!  px_out is the target array of x values
!  py_out is the output array of corresponding y values
!  kn  is the number of input (x_in,y_in) data points
!  px_in  is an array of input x values (e.g., altitudes)
!  py_in  is an array of input y values (e.g., temperatures or number densities)
!       corresponding to the values of x
!  ld_ok is return .true. if all ok
!  cd_ext_vals :     (to deal with values of target grid which are out of
!                  input grid)
!             'notouch' to not touch these values
!             'allzero' to set these output values to zero (DEFAULT)
!             'asbndry' to set these values to the boundary input values
!             'extrapo' to extrapolate linearly outside of the target grid.
!                   DANGEROUS! At least the extrapolated values won't
!                   change sign (floor/ceiling set to zero)
!-----------------------------------------------------------------------

USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULERR

IMPLICIT NONE

!-----------------------------------------------------------------------
!   ... Dummy args
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM), INTENT(IN) :: KN_OUT, KN
REAL(KIND=JPRB), INTENT(IN)    :: PX_IN(KN), PY_IN(KN), PX_OUT(KN_OUT)
REAL(KIND=JPRB), INTENT(INOUT) :: PY_OUT(KN_OUT)
!VH LOGICAL, optional, INTENT(IN)  :: LD_SPLINE, LD_ON_LOG
CHARACTER(LEN=7), INTENT(IN)   :: CD_EXT_VALS
LOGICAL,           INTENT(OUT) :: LD_OK

!-----------------------------------------------------------------------
!   ... Local variables
!-----------------------------------------------------------------------
LOGICAL :: LL_LOG_INTERP, LL_INTERP_SPLINE
INTEGER(KIND=JPIM) :: JI, ILO, IHI, JK, JKLO, JKHI
REAL(KIND=JPRB) :: ZP,  ZSIG, ZA, ZB, ZH
REAL(KIND=JPRB), DIMENSION(KN)     :: ZX, ZY, ZY2, ZU
REAL(KIND=JPRB), DIMENSION(KN_OUT) :: ZX_OUT1
CHARACTER(LEN=7) :: CL_EXT
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',0,ZHOOK_HANDLE )

      LD_OK = .FALSE.
      LL_LOG_INTERP = .FALSE.
!VH      if( PRESENT( LD_ON_LOG ) ) LL_LOG_INTERP = LD_ON_LOG
      LL_INTERP_SPLINE = .FALSE.
!VH      if( PRESENT( LD_SPLINE ) ) LL_INTERP_SPLINE = LD_SPLINE
      CL_EXT = 'allzero'
      !VH if( PRESENT( CD_EXT_VALS ) ) CL_EXT = CD_EXT_VALS
      CL_EXT = CD_EXT_VALS
      IF( CL_EXT/='notouch' .AND. CL_EXT/='allzero' .AND. CL_EXT/='asbndry' .AND. &
     &    CL_EXT/='extrapo' ) THEN
         WRITE(NULERR,*) 'INTERP: Error: "CD_EXT_VALS" keyword not understood'
         IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
         RETURN
      ENDIF


!-----------------------------------------------------------------------
!  Make sure vectors conform to increasing x.
!-----------------------------------------------------------------------
      IF( PX_IN(2) < PX_IN(1) ) THEN
         ZX(KN:1:-1) = PX_IN(1:KN)
         ZY(KN:1:-1) = PY_IN(1:KN)
       ELSE
         ZX = PX_IN
         ZY = PY_IN
      ENDIF
      ZX_OUT1 = PX_OUT
      IF( KN_OUT > 1 ) THEN
         IF( PX_OUT(2) < PX_OUT(1) ) ZX_OUT1(KN_OUT:1:-1) = PX_OUT(1:KN_OUT)
      ENDIF

!-----------------------------------------------------------------------
!  Take log of values IF necessary
!-----------------------------------------------------------------------
      IF( LL_LOG_INTERP ) THEN
         IF ( ANY( ZY(1:KN) <= 0. ) ) THEN
            WRITE(NULERR,*) ' INTERP: Error: some input data is neg/null,'
            WRITE(NULERR,*) '   but log interpolation was specified.'
           IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
           RETURN
         ENDIF
         ZY(1:KN) = LOG( ZY(1:KN) )
      ENDIF

!-----------------------------------------------------------------------
!  test if data is sorted
!-----------------------------------------------------------------------
      IF( ANY( ZX(2:KN) < ZX(1:KN-1) ) ) THEN
         WRITE(NULERR,*)' INTERP: Error: input grid not sorted'
         IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
         RETURN
      ENDIF
      IF( ANY( ZX_OUT1(2:KN_OUT) <= ZX_OUT1(1:KN_OUT-1) ) ) THEN
         WRITE(NULERR,*)' INTERP: Error: target grid not sorted'
         IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
         RETURN
      ENDIF
      IF( (ZX_OUT1(1) >= ZX(KN)) .OR. (ZX_OUT1(KN_OUT) <= ZX(1)) ) THEN
          WRITE(NULERR,*)' INTERP: Error: target grid outside of input grid'
          IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
        RETURN
      ENDIF

      IF( LL_INTERP_SPLINE ) THEN
!-----------------------------------------------------------------------
!  Calculate ZY2, the 2d derivative of interpolating function at input ZX
!-----------------------------------------------------------------------
         ZY2(1) = 0._JPRB
         ZY2(KN) = 0._JPRB
         ZU(1) = 0._JPRB
         DO JI= 2, KN-1
            ZSIG = ( ZX(JI) - ZX(JI-1) ) / ( ZX(JI+1) - ZX(JI-1) )
            ZP = ZSIG * ZY2(JI-1) + 2._JPRB
            ZY2(JI) = ( ZSIG - 1. ) / ZP
            ZU(JI) = ( 6. * ( (ZY(JI+1)-ZY(JI)) / (ZX(JI+1)-ZX(JI)) &
     &                    - (ZY(JI)-ZY(JI-1)) / (ZX(JI)-ZX(JI-1)) ) &
     &                  / (ZX(JI+1)-ZX(JI-1)) - ZSIG * ZU(JI-1) ) &
     &           / ZP
         ENDDO
         DO JI = KN-1, 1, -1
           ZY2(JI) = ZY2(JI)*ZY2(JI+1) + ZU(JI)
         ENDDO
      ENDIF

!-----------------------------------------------------------------------
!  Find interval (ILO,IHI) where the target grid is inside the input grid
!-----------------------------------------------------------------------
      DO ILO = 1, KN_OUT
         IF( ZX_OUT1(ILO) >= ZX(1) ) exit
      ENDDO
      DO IHI = KN_OUT, 1, -1
         IF( ZX_OUT1(IHI) < ZX(KN) ) exit
      ENDDO

!-----------------------------------------------------------------------
!  Interpolate linearily or on cubic spline between ILO and IHI
!-----------------------------------------------------------------------
      JKLO = 1
      JKHI = 2
      DO JI = ILO, IHI
         IF( (ZX(JKLO) > ZX_OUT1(JI)) .OR. (ZX(JKHI) <= ZX_OUT1(JI)) ) THEN
            DO JK = JKLO, KN-1
               IF( (ZX(JK) <= ZX_OUT1(JI)) .AND. (ZX(JK+1) > ZX_OUT1(JI)) ) &
     &            EXIT
            ENDDO
            JKLO = JK
            JKHI = JK + 1
         ENDIF
         ZH = ZX(JKHI) - ZX(JKLO)
         ZA = ( ZX(JKHI) - ZX_OUT1(JI) ) / ZH
         ZB = ( ZX_OUT1(JI) - ZX(JKLO) ) / ZH
         PY_OUT(JI) = ZA*ZY(JKLO) + ZB*ZY(JKHI)
         IF( LL_INTERP_SPLINE ) PY_OUT(JI) = PY_OUT(JI) &
     &            + ( (ZA*ZA*ZA-ZA)*ZY2(JKLO) + (ZB*ZB*ZB-ZB)*ZY2(JKHI) )*ZH*ZH / 6.
      ENDDO

!-----------------------------------------------------------------------
!  if ext == 'extrapo', extrapolate linearly below and/or above x(n)
!-----------------------------------------------------------------------
      IF( CL_EXT == 'extrapo' .AND. (IHI < KN_OUT) ) THEN
         ZA = (ZY(KN) - ZY(KN-1)) / (ZX(KN) - ZX(KN-1))
         PY_OUT(IHI+1:KN_OUT) = PY_OUT(IHI) &
     &               + ZA * (ZX_OUT1(IHI+1:KN_OUT) - ZX_OUT1(IHI))
         IF( .not. LL_LOG_INTERP ) THEN
            IF( PY_OUT(IHI) > 0. ) THEN
               PY_OUT(IHI+1:KN_OUT) = MAX( PY_OUT(IHI+1:KN_OUT), 0._JPRB )
             ELSE
               PY_OUT(IHI+1:KN_OUT) = MIN( PY_OUT(IHI+1:KN_OUT), 0._JPRB )
            ENDIF
         ENDIF
      ENDIF
      IF( CL_EXT == 'extrapo' .AND. (ILO > 1) ) THEN
         ZA = (ZY(2) - ZY(1)) / (ZX(2) - ZX(1))
         PY_OUT(1:ILO-1) = PY_OUT(ILO) &
     &               + ZA * (ZX_OUT1(1:ILO-1) - ZX_OUT1(ILO))
         IF( .not. LL_LOG_INTERP ) THEN
            IF( PY_OUT(ILO) > 0. ) THEN
               PY_OUT(1:ILO-1) = MAX( PY_OUT(1:ILO-1), 0._JPRB )
             ELSE
               PY_OUT(1:ILO-1) = MIN( PY_OUT(1:ILO-1), 0._JPRB )
            ENDIF
         ENDIF
      ENDIF

!-----------------------------------------------------------------------
!  Get back to original order, make a last test
!-----------------------------------------------------------------------
      IF( LL_LOG_INTERP ) THEN
         IF( CL_EXT == 'extrapo' ) THEN
            PY_OUT(1:KN_OUT) = EXP( PY_OUT(1:KN_OUT) )
          ELSE
            PY_OUT(ILO:IHI) = EXP( PY_OUT(ILO:IHI) )
         ENDIF
      ENDIF
      IF( CL_EXT == 'allzero' ) THEN
         IF( ILO > 1 ) PY_OUT(1:ILO-1) = 0.
         IF( IHI < KN_OUT ) PY_OUT(IHI+1:KN_OUT) = 0.
       ELSEIF( CL_EXT == 'asbndry' ) THEN
         IF( ILO > 1 ) PY_OUT(1:ILO-1) = PY_OUT(ILO)
         IF( IHI < KN_OUT ) PY_OUT(IHI+1:KN_OUT) = PY_OUT(IHI)
      ENDIF
      IF( KN_OUT > 1 ) THEN
         IF( PX_OUT(2) < PX_OUT(1) ) PY_OUT(1:KN_OUT) = PY_OUT(KN_OUT:1:-1)
      ENDIF
      IF( PX_OUT(KN_OUT) == PX_IN(KN) ) PY_OUT(KN_OUT) = PY_IN(KN)
      IF( ALL( PY_IN > 0. ) .AND. ANY( PY_OUT < 0. ) ) THEN
         WRITE(NULERR,*)' INTERP: Error: all input > 0., some output < 0.'
         IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
         RETURN
      ENDIF

      LD_OK = .TRUE.

IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP8',1,ZHOOK_HANDLE )
END SUBROUTINE BASCOE_INTERP8

