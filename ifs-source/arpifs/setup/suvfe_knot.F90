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

SUBROUTINE SUVFE_KNOT(YDVFE,YDCVER,LDFIX_ORDER, KTBC, KBBC, KBASIS, KORDER, &
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
!     LDFIX_ORDER  : T/F = fixed spline order/fixed knots
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
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPRB, JPIM
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMVERT  , ONLY : TVFE
USE YOMCVER  , ONLY : TCVER

!-------------------------------------------------

IMPLICIT NONE

TYPE(TVFE)        , INTENT(INOUT) :: YDVFE
TYPE(TCVER)       , INTENT(IN)  :: YDCVER
LOGICAL           , INTENT(IN)  :: LDFIX_ORDER
INTEGER(KIND=JPIM), INTENT(IN)  :: KTBC(2), KBBC(2)
INTEGER(KIND=JPIM), INTENT(IN)  :: KBASIS
INTEGER(KIND=JPIM), INTENT(IN)  :: KORDER
INTEGER(KIND=JPIM), INTENT(IN)  :: KFLEV
REAL   (KIND=JPRB), INTENT(IN)  :: PETA(KFLEV)
REAL   (KIND=JPRB), INTENT(OUT) :: PKNOT(KBASIS+KORDER)

!-------------------------------------------------

REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
REAL(KIND=JPRB)    :: ZK(KFLEV+KORDER), ZDK, ZTK, ZBK

INTEGER(KIND=JPIM) :: IKNOTS_BC, INTERNALS_BC
INTEGER(KIND=JPIM) :: ITBC, IBBC, IOFF
INTEGER(KIND=JPIM) :: II, IJ, IK, IMUL
INTEGER(KIND=JPIM) :: IKNOTS, INTERNALS, IFRST
LOGICAL            :: LLPERCENTILS
REAL(KIND=JPRB)    :: ZW, ZINDX, ZPERC, ZSTRETCH, ZLEVS, ZNOTS

!-------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_KNOT',0,ZHOOK_HANDLE)
!-------------------------------------------------

! number of knots
IKNOTS_BC = KBASIS + KORDER

! number of internal nodes
INTERNALS_BC = KBASIS - KORDER

IF( LDFIX_ORDER )THEN

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

  LLPERCENTILS = .TRUE.

  IF(LLPERCENTILS)THEN

    ! multiplicity knots
    DO II = 1, KORDER
      PKNOT(II) = 0.0_JPRB
      PKNOT(KBASIS + KORDER - II + 1) = 1.0_JPRB
    ENDDO

    ! percentils
    DO II = 1, KBASIS - KORDER

       ZNOTS    = REAL(KBASIS - KORDER, JPRB)
       ZLEVS    = REAL(KFLEV, JPRB)
       ! ZSTRETCH = (YDCVER%RVFE_KNOT_STRETCH - ZNOTS + YDCVER%RVFE_KNOT_STRETCH * ZNOTS -  ZLEVS) &
       !         & / (1.0_JPRB + ZLEVS - 2.0_JPRB * YDCVER%RVFE_KNOT_STRETCH)

       ZSTRETCH = 0.5_JPRB * (ZLEVS - 2.0_JPRB - ZNOTS)
       ZPERC =  (REAL(II , JPRB) + ZSTRETCH) / (REAL(KBASIS - KORDER + 1, JPRB) + 2.0_JPRB * ZSTRETCH)
       ZINDX =  ZPERC * REAL(KFLEV - 1, JPRB) + 1.0_JPRB

       ! percentil is located in interval <PETA(IJ), PETA(IJ + 1)>
       IJ = INT(ZINDX, JPIM)

       IK = KORDER + II

       IF(PETA(IJ) == 0.0_JPRB)THEN
         PKNOT(IK) = PETA(IJ + 1)
       ELSEIF(PETA(IJ) == 1.0_JPRB)THEN
         PKNOT(IK) = PETA(IJ - 1)
       ELSE

         ! IF(LVFE_REGETA)THEN
         ! IJ = NINT(ZINDX,JPIM)
         ! PKNOT(IK) = PETA(IJ)
         ! ELSE
         ZW = ZINDX - INT(ZINDX,JPIM)
         PKNOT(IK) = (1.0_JPRB - ZW) * PETA(IJ) + ZW * PETA(IJ + 1)
         ! ENDIF
          
       ENDIF

       ! IF(YDCVER%LVFE_VERBOSE)THEN
       !   WRITE(NULOUT,'("DBG KNOT :: ",I4,2(1X,F10.5),1X,I4,3(1X,F10.7))') II, ZPERC*100.0_JPRB, ZINDX, IJ, PETA(IJ), PKNOT(IK), ZSTRETCH
       ! ENDIF

    ENDDO

    ! IMUL = KORDER - 2
    ! ZTK  = PKNOT(KORDER + 4)
    ! ZBK  = PKNOT(KBASIS - 3)
    DO II = 1, IMUL
     !  PKNOT(KORDER + II    ) = ZTK
    !   PKNOT(KBASIS - II + 1) = ZBK
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

    IF( ( ITBC + IBBC ) > 0 )THEN
      PKNOT(1:KORDER) = ZK(1:KORDER)
      IOFF = KORDER
      DO II = 1, ITBC
        IK = II + IOFF
        IJ = II + IFRST - ITBC
        IF( YDCVER%NVFE_BC == 0 )THEN
          ! use input levels as boundary knots
          IF(MOD(KORDER,2)==0)THEN
            PKNOT(IK) = PETA(IJ)
          ELSE
            PKNOT(IK) = (PETA(IJ) + PETA(IJ + 1)) * 0.5_JPRB
          ENDIF
        ELSEIF( YDCVER%NVFE_BC == 1 )THEN
          ! regular distribution of boundary knots
          ZDK  = (ZK(KORDER + 1) - ZK(KORDER)) / REAL(ITBC + 1, JPRB)
          PKNOT(IK) = ZDK * REAL(II, JPRB)
        ELSEIF( YDCVER%NVFE_BC == 2 )THEN
          ! increased multiplicity of knots (decrease order of spline)
          ZDK  = 0.5_JPRB * (ZK(KORDER + 1) + ZK(KORDER))
          PKNOT(IK) = ZDK
        ELSE
          CALL ABOR1("SUVFE_KNOT: unknown NVFE_BC")
        ENDIF
      ENDDO

      IOFF = IOFF + ITBC
      PKNOT(IOFF + 1 : IOFF + INTERNALS) = ZK(KORDER + 1 : KORDER + INTERNALS)
      IOFF = IOFF + INTERNALS

      DO II = 1, IBBC
        IK = II + IOFF
        IJ = II + IFRST + INTERNALS
        IF( YDCVER%NVFE_BC == 0 )THEN
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
        ELSEIF( YDCVER%NVFE_BC == 1 )THEN
          ZDK      = (ZK(KFLEV + 1) - ZK(KFLEV)) / REAL(IBBC + 1, JPRB)
          PKNOT(IK) = ZK(KFLEV) + ZDK * REAL(II, JPRB)
        ELSEIF( YDCVER%NVFE_BC == 2 )THEN
          PKNOT(IK) = ZK(KORDER+INTERNALS)
        ELSE
          CALL ABOR1("SUVFE_KNOT: unknown NVFE_BC")
        ENDIF
      ENDDO

      IOFF = IOFF + IBBC
      PKNOT(IOFF + 1 : IOFF + KORDER) = ZK(KFLEV + 1 : KFLEV + KORDER)


    ELSE
      PKNOT = ZK 
    ENDIF

  ENDIF

ELSE ! LDFIX_ORDER

  IF( INTERNALS_BC /= YDCVER%NVFE_INTERNALS )THEN
    CALL ABOR1("(KNOT) INCONSISTENT DIMENSIONS IN KNOT FOR  &
     & LVFE_FIX_ORDER=.FALSE - FIXED KNOTS")
  ENDIF
  PKNOT(1       :KORDER) = 0.0_JPRB
  PKNOT(KORDER+1:KBASIS) = YDVFE%VFE_KNOT
  PKNOT(KBASIS+1:IKNOTS_BC) = 1.0_JPRB

ENDIF ! LDFIX_ORDER

IF (YDCVER%LVFE_VERBOSE) THEN
  WRITE(NULOUT,*) "(KNOT) SEQUENCE"
  DO II = 1, IKNOTS_BC - 1
    WRITE(NULOUT,'(I3.3,1X,2F10.5)') II, PKNOT(II), PKNOT(II + 1) - PKNOT(II)
  ENDDO
  WRITE(NULOUT,'(I3.3,1X,2F10.5)') IKNOTS_BC, PKNOT(IKNOTS_BC)
ENDIF

IF (LHOOK) CALL DR_HOOK('SUVFE_KNOT',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_KNOT
