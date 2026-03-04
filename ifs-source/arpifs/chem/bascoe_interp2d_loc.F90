! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_INTERP2D_LOC(KX, KY, PX, PY, LDXLON, PXO, PYO, KDX, PRD, k_ok  )


!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_INTERP2D_LOC* IS CALLED FROM *BASCOE_J_INTERP*.

! INPUTS:
! -------
!
!
! OUTPUTS:
! -------
!
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2014-02-01
!
! LOCAL:
! -------
USE PARKIND1 , ONLY : JPIM ,   JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
!-----------------------------------------------------------------------
!   ... Dummy args.
!  If ldxlon==.true., x is longitudes from 0 to (360-dx) degrees by dx
!-----------------------------------------------------------------------
INTEGER(KIND=JPIM), INTENT(IN)  :: KX, KY          ! dims of 2D grid
REAL(KIND=JPRB)   , INTENT(IN)  :: PX(KX), PY(KY)  ! vectors defining 2D grid - MUST increase
REAL(KIND=JPRB)   , INTENT(IN)  :: PXO, PYO        ! coords of point to interpolate
LOGICAL,            INTENT(IN)  :: LDXLON           ! .true. if x=longitudes (0 to 360-dx degrees)

INTEGER(KIND=JPIM), INTENT(OUT) :: KDX(4)         ! vector of 4 indexes of neighbour gridpoints
REAL(KIND=JPRB)   , INTENT(OUT) :: PRD(2)         ! ratios of x distance and y distance
INTEGER(KIND=JPIM), INTENT(OUT) :: k_ok            ! 0 if everyhting OK


! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

REAL(KIND=JPRB), parameter :: Zxmissval = 1.E+29_JPRB
REAL(KIND=JPRB) :: ZX_IN(KX), ZXTMP(KX), ZYTMP(KY), ZXOUT, ZYOUT
INTEGER(KIND=JPIM) :: id(1), ixm, iym, ixp, iyp

#include "abor1.intfb.h"

!################  subroutine start ################

IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',0,ZHOOK_HANDLE )

!-----------------------------------------------------------------------
!  Check input
!-----------------------------------------------------------------------
k_ok = 0
IF( PX(1) > PX(2) .or. PY(1) > PY(2) ) THEN
   CALL ABOR1(  'INTERP2D_LOC error: grid vectors (x,y) MUST increase' ) 
ENDIF
IF( LDXLON .and. PX(1) /= 0. ) THEN
   CALL ABOR1(  'INTERP2D_LOC error:1st longitude MUST be zero degrees ' ) 
ENDIF
IF( LDXLON .and. PX(KX) >= 360._JPRB ) THEN
   CALL ABOR1( 'INTERP2D_LOC error: last longitude MUST be < 360 degrees' )
ENDIF
ZX_IN = PX
ZXOUT = PXO
IF( LDXLON .and. ZXOUT < 0. ) ZXOUT = ZXOUT + 360.
IF( ZXOUT < ZX_IN(1) ) THEN
   k_ok = 1
   IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',1,ZHOOK_HANDLE )
   RETURN
 ELSEIF( ZXOUT > ZX_IN(KX) ) THEN
   IF( LDXLON .and. ZXOUT <= 3.6E2_JPRB ) THEN
      ZX_IN(1:KX-1) = PX(2:KX)
      ZX_IN(KX) = 3.6E2_JPRB
    ELSE
      k_ok = 2
      IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',1,ZHOOK_HANDLE )
      RETURN
   ENDIF
ENDIF
ZYOUT = PYO
IF( ZYOUT < PY(1) ) THEN
   k_ok = 3
   IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',1,ZHOOK_HANDLE )
   RETURN
ELSEIF( ZYOUT > PY(KY) ) THEN
   k_ok = 4
   IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',1,ZHOOK_HANDLE )
   RETURN
ENDIF

!-----------------------------------------------------------------------
!  Get x indexes and ratios of distances
!-----------------------------------------------------------------------
ZXTMP = ZX_IN
where( ZXTMP(:KX) < ZXOUT ) ZXTMP(:KX) = Zxmissval
id = MINLOC( ZXTMP )
ixp = id(1)
ixm = ixp - 1
IF( ixm == 0 ) THEN
   IF( ZXOUT /= ZX_IN(1) ) THEN
     CALL ABOR1( 'INTERP2D_LOC internal error: ixm=0' )
   ENDIF
   ixm = ixp
   PRD(1) = 0.
ELSE
   PRD(1) = ( ZXOUT - ZX_IN(ixm) ) / ( ZX_IN(ixp) - ZX_IN(ixm) )
ENDIF

!-----------------------------------------------------------------------
!  Get y indexes and ratios of distances
!-----------------------------------------------------------------------
ZYTMP = PY
where( ZYTMP(:KY) < ZYOUT ) ZYTMP(:KY) = Zxmissval
id = MINLOC( ZYTMP )
iyp = id(1)
iym = iyp - 1
IF( iym == 0 ) THEN
   IF( ZYOUT /= PY(1) ) THEN
     CALL ABOR1('TERP2D_LOC internal error: iym=0' )
   ENDIF
   iym = iyp
   PRD(2) = 0.
 ELSE
   PRD(2) = ( ZYOUT - PY(iym) ) / ( PY(iyp) - PY(iym) )
ENDIF

KDX = (/ ixm, iym, ixp, iyp /)


IF (LHOOK) CALL DR_HOOK('BASCOE_INTERP2D_LOC',1,ZHOOK_HANDLE )


END SUBROUTINE BASCOE_INTERP2D_LOC
