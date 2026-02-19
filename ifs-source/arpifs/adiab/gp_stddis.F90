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

SUBROUTINE GP_STDDIS( &
 ! --- INPUT --------------------------------------------------
 & YDDYNA,KPROMA,KSTART,KPROF,KFLEV,PDT,&
 & PUT0L,PDIVT0,&
 ! --- OUTPUT -------------------------------------------------
 & PSTDDISU,PSTDDISV,PSTDDISW)

!**** *GP_STDDIS*

!     Purpose.
!     --------

!**   Interface.
!     ----------
!        *CALL* *GP_STDDIS(..)

!        Explicit arguments :
!        --------------------

!        INPUT:
!          KPROMA      - horizontal dimension.
!          KSTART      - first element of work.
!          KPROF       - depth of work.
!          KFLEV       - number of layers.
!          PDT         - time step
!          PUT0L       - zonal derivative of U-wind at time t.
!          PDIVT0      - HOR. DIV. of wind at time t.

!        OUTPUT:
!          PSTDDISU   - stretching/shrinking deformation along zonal direction
!          PSTDDISV   - stretching/shrinking deformation along meridional direction
!          PSTDDISW   - stretching/shrinking deformation along vertical direction

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           none
!           Called by LATTE_STDDIS.

!     Reference.
!     ----------
!             Arpege documentation about semi-Lagrangian scheme.

!  Author.
!  -------
!    Original (S.Malardel): November 2013.

!  Modifications.
!  --------------
!  F. Vana 26-Nov-2921  Better memory handling

! End Modifications
! ------------------------------------------------------------------

USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYNA,  ONLY : TDYNA

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDYNA)       ,INTENT(IN)    :: YDDYNA
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT0L(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDIVT0(KPROMA,KFLEV) 

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTDDISU(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTDDISV(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTDDISW(KPROMA,KFLEV) 

!     -------------------------------------------------------

REAL(KIND=JPRB)    :: ZVT0M
REAL(KIND=JPRB)    :: ZALPHA
REAL(KIND=JPRB)    :: ZBND1, ZBND2
INTEGER(KIND=JPIM) :: JLEV,JROF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GP_STDDIS',0,ZHOOK_HANDLE)
!    --------------------------------------------------------

!*       0.   Initialisation

ZBND1=1.0_JPRB
ZBND2=0.0001_JPRB

!    --------------------------------------------------------
!*       1.   Streching of distance in the zonal direction
!                             and the meridional direction

IF (YDDYNA%LCOMADH) THEN
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ! zonal direction dx(t)/dx(t+dt)
      ZALPHA=1.0_JPRB+PUT0L(JROF,JLEV)*PDT
      PSTDDISU(JROF,JLEV)=MIN(ZBND1,MAX(ZBND2,ZALPHA))
      ! meridional direction dy(t)/dy(t+dt)
      ZVT0M=PDIVT0(JROF,JLEV)-PUT0L(JROF,JLEV)
      ZALPHA=1.0_JPRB+ZVT0M*PDT
      PSTDDISV(JROF,JLEV)=MIN(ZBND1,MAX(ZBND2,ZALPHA))
    ENDDO
  ENDDO
ELSE
  ! with stretching factor=1, COMAD correction are neutralized
  PSTDDISU(KSTART:KPROF,1:KFLEV)=1._JPRB
  PSTDDISV(KSTART:KPROF,1:KFLEV)=1._JPRB
ENDIF

!    --------------------------------------------------------

!*       2.   Streching of distance in the vertical direction

! Compute vertical weight correction as if 3D non divergent flow (D3=0)
IF (YDDYNA%LCOMADV) THEN
  IF (YDDYNA%LNHDYN) CALL ABOR1(" Vertical direction in COMAD not coded for NH.")
  DO JLEV=1,KFLEV
    DO JROF=KSTART,KPROF
      ZALPHA=1.0_JPRB-PDIVT0(JROF,JLEV)*PDT
      PSTDDISW(JROF,JLEV)=MIN(ZBND1,MAX(ZBND2,ZALPHA))
    ENDDO
  ENDDO
ELSE
  ! with stretching factor=1, COMAD correction are neutralized
  PSTDDISW(KSTART:KPROF,1:KFLEV)=1._JPRB
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GP_STDDIS',1,ZHOOK_HANDLE)
END SUBROUTINE GP_STDDIS
