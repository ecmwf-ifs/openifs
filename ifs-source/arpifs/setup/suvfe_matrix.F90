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

SUBROUTINE SUVFE_MATRIX(YDCVER,LDINT_FROM_SURF,KTYPE, &
   & KKNOT_IN,PKNOT_IN,KORDER_IN,PTM_IN,KBASIS_IN,KOFF_IN,PBAF_IN,  &
   & KKNOT_W,PKNOT_W,KORDER_W,PTM_W,KBASIS_W,KOFF_W,PWEI,  &
   & KFLEV,PVFE   )

!**** *SUVFE_MATRIX*  - Set Up VFE vertical operator (B-splines): compute MATRIX

!**   Interface.
!     ----------
!     *CALL* SUVFE_MATRIX

!     Explicit arguments :
!     --------------------
!      * INPUT:
!        LDINT_FROM_SURF ! T/F integral from surface/top; only for KTYPE=-1
!        KTYPE        ! control what this routine computes
!                     * KTYPE =-1 - compute stiffness matrix for integral
!                         M(j,i)   - Int_(0,1) ( ( Int_(0,eta) e(i) deta` )
!                                    *w(j) ) deta, i=1,KFLEV  j=1,KFLEV
!                     * KTYPE = 0 - compute mass matrix
!                         M(j,i)   - Int_(0,1) ( e(i)*w(j) ) deta, i=1,KFLEV
!                                    j=1,KFLEV
!                     * KTYPE = 1 - compute stiffness matrix for first derivative
!                         M(j,i)   - Int_(0,1) ( de(i)/deta *w(j) ) deta,
!                                    i=1,KFLEV  j=1,KFLEV
!                     * KTYPE = 2 - compute stiffness matrix for second derivative
!                         M(j,i)   - Int_(0,1) ( d2 e(i)/deta^2 *w(j) ) deta,
!                                    i=1,KFLEV  j=1,KFLEV
!        KKNOT_IN     ! number of knots of e(i)
!        PKNOT_IN     ! knots of e(i)
!        KORDER_IN    ! order of B-spline e(i)
!        PTM_IN       ! transform matrix for order KORDER_IN spline
!        KBASIS_IN    ! number of e(i) functions
!        KOFF_IN      ! index of e(i), i=1+KOFF_IN,KFLEV+KOFF_IN in PBAF_IN field
!        PBAF_IN      ! basis functions e(i)
!        KKNOT_W      ! number of knots of w(i)
!        PKNOT_W      ! knots of w(i)
!        KORDER_W     ! order of B-spline w(i)
!        PTM_W        ! transform matrix for order KORDER_W spline
!        KBASIS_W     ! number of w(i) functions
!        KOFF_W       ! index of w(i), i=1+KOFF_W,KFLEV+KOFF_W in PWEI field
!        PWEI         ! weighting functions w(i)
!        KFLEV        ! number of model full levels

!      * OUTPUT:
!        PVFE         ! mass/stiffness matrix

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ALADIN/LACE documentation for NH dynamics.

!     Author.
!     -------
!        Jozef Vivoda, SHMU/LACE 
!        Original : 2010-09

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      P. Marguinaud (Oct-2016): Port to single precision
!      J. Vivoda and P. Smolikova (Sep 2017): new options for VFE-NH
! --------------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB      ,JPRD
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMMP0    ,ONLY : LOUTPUT, NPRINTLEV
USE YOMCVER   ,ONLY : TCVER
USE YOMLUN    ,ONLY : NULOUT
USE SUVFE_HLP ,ONLY : MULPOL, DPOL, IPOL, EVPOL, GLOBAL2LOCAL, FX2T, RTMIN, RTMAX

! --------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TCVER)       ,INTENT(IN)  :: YDCVER
LOGICAL           ,INTENT(IN)  :: LDINT_FROM_SURF
INTEGER(KIND=JPIM),INTENT(IN)  :: KTYPE
INTEGER(KIND=JPIM),INTENT(IN)  :: KKNOT_IN
REAL(KIND=JPRB)   ,INTENT(IN)  :: PKNOT_IN(KKNOT_IN)
INTEGER(KIND=JPIM),INTENT(IN)  :: KORDER_IN
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTM_IN(KORDER_IN,KORDER_IN)
INTEGER(KIND=JPIM),INTENT(IN)  :: KBASIS_IN
INTEGER(KIND=JPIM),INTENT(IN)  :: KOFF_IN
REAL(KIND=JPRB)   ,INTENT(IN)  :: PBAF_IN(KBASIS_IN,KORDER_IN,KORDER_IN)
INTEGER(KIND=JPIM),INTENT(IN)  :: KKNOT_W
REAL(KIND=JPRB)   ,INTENT(IN)  :: PKNOT_W(KKNOT_W)
INTEGER(KIND=JPIM),INTENT(IN)  :: KORDER_W
REAL(KIND=JPRB)   ,INTENT(IN)  :: PTM_W(KORDER_W,KORDER_W)
INTEGER(KIND=JPIM),INTENT(IN)  :: KBASIS_W
INTEGER(KIND=JPIM),INTENT(IN)  :: KOFF_W
REAL(KIND=JPRB)   ,INTENT(IN)  :: PWEI(KBASIS_W,KORDER_W,KORDER_W)
INTEGER(KIND=JPIM),INTENT(IN)  :: KFLEV

REAL(KIND=JPRB)   ,INTENT(OUT) :: PVFE(KFLEV,KFLEV)

! --------------------------------------------------------------------------

INTEGER(KIND=JPIM)              :: IUN, IOFF, ISEG, JSEG
INTEGER(KIND=JPIM)              :: JI,JJ,JK,IMIN,IMAX,IJMAX,IJMIN
REAL(KIND=JPRB)   , ALLOCATABLE :: ZKNT(:)
REAL(KIND=JPRB)   , ALLOCATABLE :: ZINT(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IW(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: IWF(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZPT_IN(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZPT_W(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZPOLY(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZDPF (:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZIPW(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZPA(:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZIPA(:)
REAL(KIND=JPRB)                 :: ZX1, ZX2 
REAL(KIND=JPRB)                 :: ZDETADT_IN, ZDETADT_W, ZSZ, ZWGH
REAL(KIND=JPRB)                 :: ZDTDLOC_IN, ZDTDLOC_W
REAL(KIND=JPRB)                 :: ZTMIN_IN, ZTMAX_IN, ZTMIN_W, ZTMAX_W, ZTW1, ZTW2
INTEGER(KIND=JPIM)              :: IFRST_IN,ILAST_IN 
INTEGER(KIND=JPIM)              :: IFRST_W,ILAST_W 
INTEGER(KIND=JPIM)              :: IDIM1, IDIM2, IDIM3, IVFE, IJVFE

! integral operator only
REAL(KIND=JPRB), ALLOCATABLE    :: ZCUM_IN (:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: ZCUM_W (:,:)
REAL(KIND=JPRB), ALLOCATABLE    :: Z_ACCUM_IN (:)
REAL(KIND=JPRB), ALLOCATABLE    :: Z_ACCUM_W (:)

!REAL(KIND=JPRB)                 :: EVPOL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! --------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_MATRIX',0,ZHOOK_HANDLE)
! --------------------------------------------------------------------------

!-----------------
! union of knots 
!-----------------
IUN = KKNOT_W + KKNOT_IN
ALLOCATE( ZKNT ( IUN ) )
ALLOCATE( ZINT ( IUN ) )
ALLOCATE( IW   ( 0:IUN-1 ) )
ALLOCATE( IWF  ( 0:IUN-1 ) )

IOFF = 0
ZKNT(IOFF+1:IOFF+KKNOT_IN) = PKNOT_IN
ZINT(IOFF+1:IOFF+KKNOT_IN) = 1.0_JPRB
IOFF = KKNOT_IN
ZKNT(IOFF+1:IOFF+KKNOT_W) = PKNOT_W
ZINT(IOFF+1:IOFF+KKNOT_W) = 2.0_JPRB

IF (JPRB == JPRD) THEN
  CALL DSORT(ZKNT,ZINT,IUN,2)
ELSE
  CALL SSORT(ZKNT,ZINT,IUN,2)
ENDIF

IWF(0) = 0
IW(0) = 0
DO JI=1,IUN-1
  IF( NINT(ZINT(JI)) == 1 )THEN
    IWF(JI) = IWF(JI-1) + 1  
  ELSE
    IWF(JI) = IWF(JI-1) 
  ENDIF
  IF( NINT(ZINT(JI)) == 2 )THEN
    IW(JI) = IW(JI-1) + 1  
  ELSE
    IW(JI) = IW(JI-1) 
  ENDIF
ENDDO

!------------
! dimensions and allocation
!------------
ALLOCATE( ZPT_IN  (KORDER_IN) )
ALLOCATE( ZPT_W   (KORDER_W) )

IF( KTYPE == 0 )THEN

  IDIM1 = KORDER_IN
  IDIM2 = KORDER_W

ELSEIF( KTYPE == 1 )THEN

  IDIM1 = KORDER_IN-1
  IDIM2 = KORDER_W

ELSEIF( KTYPE == 2 )THEN

  IDIM1 = KORDER_IN-2
  IDIM2 = KORDER_W

  ALLOCATE( ZDPF(KORDER_IN-1) )

ELSEIF( KTYPE == -1 )THEN

  IDIM1 = KORDER_IN+1
  IDIM2 = KORDER_W

  ALLOCATE( ZCUM_IN(KBASIS_IN,IUN) )
  ALLOCATE( Z_ACCUM_IN(KBASIS_IN) )
  
  ALLOCATE( ZCUM_W(KBASIS_W,IUN) )
  ALLOCATE( Z_ACCUM_W(KBASIS_W) )
  
  ALLOCATE( ZIPW(KORDER_W+1) )

ELSE

  CALL ABOR1("SUVFE_MATRIX: (ERROR) VFE setup, wrong KTYPE ...")

ENDIF

IDIM3 = IDIM1 + IDIM2 - 1
ALLOCATE( ZPOLY( IDIM1 ) )
ALLOCATE( ZPA  ( IDIM3 ) )
ALLOCATE( ZIPA ( IDIM3 + 1 ) )

IFRST_W = 1     + KOFF_W
ILAST_W = KFLEV + KOFF_W
IFRST_IN = 1     + KOFF_IN
ILAST_IN = KFLEV + KOFF_IN

!------------
! initial values
!------------
PVFE(:,:) = 0.0_JPRB

!-------------------------------
! accumulate integral of e(i) 
! - needed only for integral operator
!-------------------------------
IF( KTYPE == -1 )THEN
  ZCUM_IN = 0.0_JPRB
  Z_ACCUM_IN = 0.0_JPRB

  DO JK=1,IUN-1

    ZSZ  = ZKNT(JK+1) - ZKNT(JK)
    Z_ACCUM_IN = 0.0_JPRB

    IF( ABS(ZSZ) > YDCVER%RMINDETA )THEN
      IMIN = MAX(IFRST_IN,IWF(JK) - KORDER_IN + 1)
      IMAX = MIN(ILAST_IN,IWF(JK)               ) 

      ZTMIN_IN = FX2T(PKNOT_IN(IWF(JK)), PKNOT_IN(IWF(JK)+1), ZKNT(JK  ))
      ZTMAX_IN = FX2T(PKNOT_IN(IWF(JK)), PKNOT_IN(IWF(JK)+1), ZKNT(JK+1))

      ZDETADT_IN = PKNOT_IN(IWF(JK)+1) - PKNOT_IN(IWF(JK))

      DO JI=IMIN,IMAX
        ISEG = IWF(JK) - JI + 1

        CALL IPOL(KORDER_IN,PBAF_IN(JI,ISEG,:),0.0_JPRB,IDIM1,ZPOLY) 
        ZX1     = EVPOL(IDIM1,ZPOLY,ZTMIN_IN)
        ZX2     = EVPOL(IDIM1,ZPOLY,ZTMAX_IN)
        Z_ACCUM_IN(JI)= ( ZX2  - ZX1 ) *  ZDETADT_IN

      ENDDO
    ENDIF

    ! cumulate partial integrals
    DO JI=1,KBASIS_IN
      ZCUM_IN(JI,JK+1) = ZCUM_IN(JI,JK) + Z_ACCUM_IN(JI)
    ENDDO
  ENDDO
ENDIF

!------------
! integration 
!------------
IF( KTYPE == -1 )THEN
  Z_ACCUM_W = 0.0_JPRB
  ZCUM_W      = 0.0_JPRB
ENDIF

DO JK=1,IUN-1

  ! size of interval
  ZSZ = ZKNT(JK+1) - ZKNT(JK)

  IF( KTYPE == -1 )THEN
    Z_ACCUM_W = 0.0_JPRB
  ENDIF

  ! integrate on non-zero intervals only
  IF( ABS(ZSZ) > YDCVER%RMINDETA )THEN

    IMIN = MAX(IFRST_IN,IWF(JK) - KORDER_IN + 1)
    IMAX = MIN(ILAST_IN,IWF(JK)               ) 
    IJMIN = MAX(IFRST_W,IW(JK) - KORDER_W + 1)
    IJMAX = MIN(ILAST_W,IW(JK)               )

    ZTMIN_W = FX2T(PKNOT_W(IW(JK)), PKNOT_W(IW(JK)+1), ZKNT(JK  ))
    ZTMAX_W = FX2T(PKNOT_W(IW(JK)), PKNOT_W(IW(JK)+1), ZKNT(JK+1))

    ZTMIN_IN = FX2T(PKNOT_IN(IWF(JK)), PKNOT_IN(IWF(JK)+1), ZKNT(JK  ))
    ZTMAX_IN = FX2T(PKNOT_IN(IWF(JK)), PKNOT_IN(IWF(JK)+1), ZKNT(JK+1))

    ZDETADT_IN = PKNOT_IN(IWF(JK)+1) - PKNOT_IN(IWF(JK))
    ZDETADT_W  = PKNOT_W(IW(JK)+1)   - PKNOT_W(IW(JK))

    ZDTDLOC_IN = ZTMAX_IN - ZTMIN_IN
    ZDTDLOC_W  = ZTMAX_W  - ZTMIN_W

    ! zdruzenie intervalov 
    ! aku polohu maju KNOTS_W v lokalnej suradnici t funkcie IN
    ! ZTW1 = FX2T(PKNOT_W(IW(JK)), PKNOT_W(IW(JK)+1), PKNOT_IN(IWF(JK)))
    ! ZTW2 = FX2T(PKNOT_W(IW(JK)), PKNOT_W(IW(JK)+1), PKNOT_IN(IWF(JK+1)))

    DO JJ=IJMIN,IJMAX

      IJVFE = JJ-KOFF_W
      JSEG = IW(JK) - JJ + 1

      IF( JJ< IFRST_W .OR. JJ> ILAST_W )THEN
        WRITE(NULOUT,*) "(ERR) OVERDIMENSION OF DIM JJ"
        CALL EXIT(-1)
      ENDIF
      IF( JSEG<1 .OR. JSEG>KORDER_W )THEN
        WRITE(NULOUT,*) "(ERR) WRONG JSEG DIM"
        CALL EXIT(-1)
      ENDIF

      IF( KTYPE == -1 )THEN
        ! int_{x,0,eta} w(i) deta
        CALL IPOL(KORDER_W,PWEI(JJ,JSEG,:),0.0_JPRB,KORDER_W+1,ZIPW) 
        ZX1     = EVPOL(KORDER_W+1,ZIPW,ZTMIN_W)
        ZX2     = EVPOL(KORDER_W+1,ZIPW,ZTMAX_W)
        Z_ACCUM_W(JJ)= (ZX2 - ZX1) * ZDETADT_W
      ENDIF
 
      ! prevedme W funciu do rovnakych koordinat ako je IN funkcia
      ! CALL GLOBAL2LOCAL(KORDER_W,PTM_W,ZTW1,ZTW2,PWEI(JJ,JSEG,:),ZPT_W)
      CALL GLOBAL2LOCAL(KORDER_W,PTM_W,ZTMIN_W,ZTMAX_W,PWEI(JJ,JSEG,:),ZPT_W)

      DO JI=IMIN,IMAX

        IVFE = JI-KOFF_IN
        ISEG = IWF(JK) - JI + 1

        ! check bounds
        IF( JI< IFRST_IN .OR. JI> ILAST_IN )THEN
          WRITE(NULOUT,*) "(ERR) OVERDIMENSION OF DIM JI"
          CALL EXIT(-1)
        ENDIF
        IF( ISEG<1 .OR. ISEG>KORDER_IN )THEN
          WRITE(NULOUT,*) "(ERR) WRONG ISEG DIM"
          CALL EXIT(-1)
        ENDIF

        ! change basis into the coordinate 
        ! t on <0,1) on local interval <zknt(k),zknt(k+1))
        CALL GLOBAL2LOCAL(KORDER_IN,PTM_IN,ZTMIN_IN,ZTMAX_IN,PBAF_IN(JI,ISEG,:),ZPT_IN)

        IF( KTYPE == 0 )THEN
        
          ZPOLY = ZPT_IN
          ZWGH  = ZSZ

        ELSEIF( KTYPE == 1 )THEN

          ! de(i)/deta
          CALL DPOL(KORDER_IN,ZPT_IN,IDIM1,ZPOLY)
          ZWGH = 1.0_JPRB

        ELSEIF( KTYPE == 2 )THEN

          ! d2 e(i)/deta^2
          CALL DPOL(KORDER_IN  ,ZPT_IN,KORDER_IN-1,ZDPF)
          CALL DPOL(KORDER_IN-1,ZDPF ,IDIM1     ,ZPOLY)
          ZWGH = 1.0_JPRB/ZSZ

        ELSEIF( KTYPE == -1 )THEN

          ! int_{i,1,eta} e(i) deta 
          CALL IPOL(KORDER_IN,ZPT_IN,0.0_JPRB,IDIM1,ZPOLY)
          ! ZWGH = ZDETADT_W * ZDETADT_IN * ZDTDLOC_W * ZDTDLOC_IN
          ZWGH = ZSZ * ZSZ

        ENDIF

        CALL MULPOL(IDIM1,ZPOLY,IDIM2,ZPT_W,IDIM3,ZPA)
        CALL IPOL  (IDIM3,ZPA,0.0_JPRB,IDIM3+1,ZIPA) 
        ZX1 = EVPOL(IDIM3+1,ZIPA,RTMIN)
        ZX2 = EVPOL(IDIM3+1,ZIPA,RTMAX)
        ! ZX1 = EVPOL(IDIM3+1,ZIPA,ZTMIN_IN)
        ! ZX2 = EVPOL(IDIM3+1,ZIPA,ZTMAX_IN)

        PVFE(IJVFE,IVFE) = PVFE(IJVFE,IVFE) + ZWGH * (ZX2 - ZX1)

        IF( KTYPE == -1 )THEN
          PVFE(IJVFE,IVFE) = PVFE(IJVFE,IVFE) + Z_ACCUM_W(JJ) * ZCUM_IN(JI,JK)
        ENDIF

      ENDDO ! JI

      IF( KTYPE == -1 )THEN

        ! interaction of w(j) function with e(i) function in integral
        ! ( int_{eta,0,1} e(i) deta ) * ( int_{eta,k,k+1} w(j) deta }
        DO JI=1,IMIN-1
          IVFE = JI-KOFF_IN
          IF( IVFE > 0 )THEN
            PVFE(IJVFE,IVFE) = PVFE(IJVFE,IVFE) + Z_ACCUM_W(JJ) * ZCUM_IN(JI,JK)

            ! IF(YDCVER%LVFE_VERBOSE)THEN
            !   WRITE(NULOUT,*) "DBG PVFE 01 ::", IJVFE, IVFE, PVFE(IJVFE,IVFE), Z_ACCUM_W(JJ) * ZCUM_IN(JI,JK)
            ! ENDIF
          ENDIF
        ENDDO

      ENDIF

    ENDDO ! JJ
  ENDIF

  IF( KTYPE == -1 )THEN
    ! cumulate partial integrals of w
    DO JJ=1,KBASIS_W
      ZCUM_W(JJ,JK+1) = ZCUM_W(JJ,JK) + Z_ACCUM_W(JJ)
    ENDDO
  ENDIF

ENDDO ! JK

IF(LDINT_FROM_SURF .AND. KTYPE == -1)THEN
  ! stiffness matrix for integral from surface:
  ! int_0^1 (int_0^eta f deta) w deta = (int_0^1 f deta) * (int_0^1 w deta) -
  !                                     int_0^1 (int_eta^1 f deta) w deta 
  DO JJ=IFRST_W, ILAST_W
    IJVFE = JJ-KOFF_W
    DO JI=IFRST_IN,ILAST_IN
      IVFE = JI-KOFF_IN
      PVFE(IJVFE,IVFE) = ZCUM_IN(JI,IUN)*ZCUM_W(JJ,IUN) - PVFE(IJVFE,IVFE)
    ENDDO
  ENDDO
ENDIF


IF (LOUTPUT.AND.NPRINTLEV >= 1) THEN
  IF( KTYPE == -1 )THEN
    WRITE(NULOUT,*) "ZCUM_IN:: total aux matrix (integral)"
    DO JI=1,KBASIS_IN
      WRITE(NULOUT,'(I5,1X,300F13.8)') JI,ZCUM_IN(JI,IUN)
    ENDDO
  ENDIF
ENDIF

!------------------
! Deallocate arrays
!------------------

IF (ALLOCATED(ZKNT)) DEALLOCATE(ZKNT)
IF (ALLOCATED(ZINT)) DEALLOCATE(ZINT)
IF (ALLOCATED(ZPT_IN)) DEALLOCATE(ZPT_IN)
IF (ALLOCATED(ZPT_W)) DEALLOCATE(ZPT_W)
IF (ALLOCATED(ZPOLY)) DEALLOCATE(ZPOLY)
IF (ALLOCATED(ZDPF)) DEALLOCATE(ZDPF)
IF (ALLOCATED(ZIPW)) DEALLOCATE(ZIPW)
IF (ALLOCATED(ZPA)) DEALLOCATE(ZPA)
IF (ALLOCATED(ZIPA)) DEALLOCATE(ZIPA)
IF (ALLOCATED(ZCUM_IN)) DEALLOCATE(ZCUM_IN)
IF (ALLOCATED(ZCUM_W)) DEALLOCATE(ZCUM_W)
IF (ALLOCATED(Z_ACCUM_IN)) DEALLOCATE(Z_ACCUM_IN)
IF (ALLOCATED(Z_ACCUM_W)) DEALLOCATE(Z_ACCUM_W)
IF (ALLOCATED(IW)) DEALLOCATE(IW)
IF (ALLOCATED(IWF)) DEALLOCATE(IWF)

! --------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUVFE_MATRIX',1,ZHOOK_HANDLE)
END SUBROUTINE SUVFE_MATRIX
