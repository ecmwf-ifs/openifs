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

SUBROUTINE SUVFE_KNOT(YDVFE,YDCVER, KTBC, KBBC, KBASIS, KORDER, &
 & KFLEV, PETA, PKNOT)

!**** *SUVFE_KNOT*  - Routine to Set Up Vertical Finite Element scheme:
!                     compute KNOTs used for basis definition.

!     Purpose.
!     --------
!           Calculates the sequence of knots used for basis definition in FE.

!**   Interface.
!     ----------
!     *CALL* SUVFE_KNOT

!     Explicit arguments :
!     --------------------
!   * INPUT:
!     KTBC/KBBC    : type of top/bottom boundary conditions
!                    (=0 -> f=0; =n>0 -> all derivative up to nth order are 0)
!     KBASIS       : number of basis functions to compute
!     KORDER       : order of spline (KORDER=2 is linear basis)
!     KFLEV        : number of vertical levels on the input
!     PETA         : full level eta coordinates on the input
!   * OUTPUT: 
!     PKNOT        : computed knots (regular for now)basis

!     Method.
!     -------
!        - traditional addition of multiple knots at the boundary
!        - uniform resolution
!        - domain size <0,1>

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN/LACE documentation on NH dynamics.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE/ALADIN 
!        Original : 2010-09

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
!      P.Smolikova (Sep 2020): VFE pruning.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMVERT  , ONLY : TVFE
USE YOMCVER  , ONLY : TCVER
USE YOMCT0   , ONLY : LECMWF

!-------------------------------------------------

IMPLICIT NONE

TYPE(TVFE)        , INTENT(INOUT) :: YDVFE
TYPE(TCVER)       , INTENT(IN)    :: YDCVER
INTEGER(KIND=JPIM), INTENT(IN)    :: KTBC(2), KBBC(2)
INTEGER(KIND=JPIM), INTENT(IN)    :: KBASIS
INTEGER(KIND=JPIM), INTENT(IN)    :: KORDER
INTEGER(KIND=JPIM), INTENT(IN)    :: KFLEV
REAL   (KIND=JPRB), INTENT(IN)    :: PETA(KFLEV)
REAL   (KIND=JPRB), INTENT(OUT)   :: PKNOT(KBASIS+KORDER)

!-------------------------------------------------

INTEGER(KIND=JPIM) :: IKNOTS_BC, INTERNALS_BC
INTEGER(KIND=JPIM) :: ITBC, IBBC, IOFF, ISETBC
INTEGER(KIND=JPIM) :: II, IJ, IK
INTEGER(KIND=JPIM) :: IKNOTS, INTERNALS, IFRST

REAL(KIND=JPRB)    :: ZW, ZINDX, ZPERC, ZSTRETCH, ZLEVS, ZNOTS
REAL(KIND=JPRB)    :: ZK(KFLEV+KORDER), ZDK
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE

!-------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_KNOT',0,ZHOOK_HANDLE)
!-------------------------------------------------

! other cases kept for research purpose
ISETBC = 0

! number of knots
IKNOTS_BC = KBASIS + KORDER

! number of internal nodes
INTERNALS_BC = KBASIS - KORDER

!-----------------------------------
! KNOTS WITHOUT IMPLICIT CONDITIONS
!-----------------------------------

! number of knots
IKNOTS = KFLEV + KORDER

! number of internal nodes
INTERNALS = KFLEV - KORDER

IF( INTERNALS < 0 )THEN
  CALL ABOR1("(KNOTS) ERROR IN KNOTS. DECREASE ORDER OF SPLINES NVFE_ORDER.")
ENDIF

IF(YDCVER%LPERCENTILS)THEN

  ! multiplicity knots
  DO II = 1, KORDER
    PKNOT(II) = 0.0_JPRB
    PKNOT(KBASIS + KORDER - II + 1) = 1.0_JPRB
  ENDDO

  ! percentils
  DO II = 1, KBASIS - KORDER
     ZNOTS    = REAL(KBASIS - KORDER, JPRB)
     ZLEVS    = REAL(KFLEV, JPRB)

     ZSTRETCH = 0.5_JPRB * (ZLEVS - 2.0_JPRB - ZNOTS)
     ZPERC =  (REAL(II , JPRB) + ZSTRETCH) / &
      &       (REAL(KBASIS - KORDER + 1, JPRB) + 2.0_JPRB * ZSTRETCH)
     ZINDX =  ZPERC * REAL(KFLEV - 1, JPRB) + 1.0_JPRB

     ! percentil is located in interval <PETA(IJ), PETA(IJ + 1)>
     IJ = INT(ZINDX, JPIM)

     IK = KORDER + II

     IF(PETA(IJ) == 0.0_JPRB)THEN
       PKNOT(IK) = PETA(IJ + 1)
     ELSEIF(PETA(IJ) == 1.0_JPRB)THEN
       PKNOT(IK) = PETA(IJ - 1)
     ELSE
       ZW = ZINDX - INT(ZINDX,JPIM)
       PKNOT(IK) = (1.0_JPRB - ZW) * PETA(IJ) + ZW * PETA(IJ + 1)
     ENDIF
  ENDDO
ELSE
  ! first full level to be used as a knot
  IFRST = MAX(INT((KFLEV - INTERNALS) / 2), 1)

  ! multiple knots at material boundaries
  DO II = 1, KORDER
    ZK(II) = 0.0_JPRB
    ZK(KFLEV+II) = 1.0_JPRB
  ENDDO

  DO II=1,INTERNALS
    IK = II + KORDER
    IJ = II + IFRST
    IF( PETA(IJ) == 1.0_JPRB )THEN
      ZK(IK) = (1.0_JPRB + PETA(IJ - 1)) * 0.5_JPRB
    ELSEIF( PETA(IJ) == 0.0_JPRB )THEN
      ZK(IK) = (0.0_JPRB + PETA(IJ + 1)) * 0.5_JPRB
    ELSE
      IF(MOD(KORDER,2)==0)THEN
        ZK(IK) = PETA(IJ)
      ELSE
        ZK(IK) = (PETA(IJ) + PETA(IJ + 1)) * 0.5_JPRB
      ENDIF
    ENDIF
  ENDDO

  !-----------------------------------
  ! INJECT BCs KNOTS
  !-----------------------------------
  ITBC = KTBC(1) + KTBC(2)
  IBBC = KBBC(1) + KBBC(2)

  IF(ITBC + IBBC == 0 )THEN
    PKNOT = ZK
  ELSE
    PKNOT(1:KORDER) = ZK(1:KORDER)
    IOFF = KORDER
    DO II = 1, ITBC
      IK = II + IOFF
      IJ = II + IFRST - ITBC
      IF( ISETBC == 0 )THEN
        ! use input levels as boundary knots
        IF(MOD(KORDER,2)==0)THEN
          PKNOT(IK) = PETA(IJ)
        ELSE
          PKNOT(IK) = (PETA(IJ) + PETA(IJ + 1)) * 0.5_JPRB
        ENDIF
      ELSEIF( ISETBC == 1 )THEN
        ! regular distribution of boundary knots
        ZDK  = (ZK(KORDER + 1) - ZK(KORDER)) / REAL(ITBC + 1, JPRB)
        PKNOT(IK) = ZDK * REAL(II, JPRB)
      ELSEIF( ISETBC == 2 )THEN
        PKNOT(IK) = ZK(KORDER + 1)
      ELSE
        CALL ABOR1("SUVFE_KNOT: unknown ISETBC")
      ENDIF
    ENDDO

    IOFF = IOFF + ITBC
    PKNOT(IOFF + 1 : IOFF + INTERNALS) = ZK(KORDER + 1 : KORDER + INTERNALS)
    IOFF = IOFF + INTERNALS

    DO II = 1, IBBC
      IK = II + IOFF
      IJ = II + IFRST + INTERNALS
      IF( ISETBC == 0 )THEN
        IF(PETA(IJ) == 1.0_JPRB)THEN
          ZDK = (1.0_JPRB - PKNOT(IJ - 1)) / 2.0_JPRB
          PKNOT(IK) = PKNOT(IJ - 1) + ZDK
        ELSE
          IF(MOD(KORDER,2)==0)THEN
            PKNOT(IK) = PETA(IJ)
          ELSE
            PKNOT(IK) = (PETA(IJ) + PETA(IJ + 1)) * 0.5_JPRB
          ENDIF
        ENDIF
      ELSEIF( ISETBC == 1 )THEN
        ZDK      = (ZK(KFLEV + 1) - ZK(KFLEV)) / REAL(IBBC + 1, JPRB)
        PKNOT(IK) = ZK(KFLEV) + ZDK * REAL(II, JPRB)
      ELSEIF( ISETBC == 2 )THEN
        PKNOT(IK) = ZK(KORDER+INTERNALS)
      ELSE
        CALL ABOR1("SUVFE_KNOT: unknown ISETBC")
      ENDIF
    ENDDO

    IOFF = IOFF + IBBC
    PKNOT(IOFF + 1 : IOFF + KORDER) = ZK(KFLEV + 1 : KFLEV + KORDER)
  ENDIF
ENDIF

IF (YDCVER%LVFE_VERBOSE) THEN
  WRITE(NULOUT,*) "(KNOT) SEQUENCE"
  DO II = 1, IKNOTS_BC - 1
    WRITE(NULOUT,'(I3.3,1X,2F10.5)') II, PKNOT(II), PKNOT(II + 1) - PKNOT(II)
  ENDDO
  WRITE(NULOUT,'(I3.3,1X,2F10.5)') IKNOTS_BC, PKNOT(IKNOTS_BC)
ENDIF

IF (LHOOK) CALL DR_HOOK('SUVFE_KNOT',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_KNOT
