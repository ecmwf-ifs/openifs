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

SUBROUTINE CPQTUV ( YDECLDP,YGFL,KPROMA, KSTART, KPROF, KFLEV, &
 & PDT, LDNODRYFLX,&
 !  Variables 2D Input
 & PUT9 ,PVT9 ,PTT9 ,PGFL, PGFL_DYN, PCONVCTY,&
 !  Variables 2D In/Out
 & PUT1 ,PVT1 ,PTT1 ,PGFLT1 )  

!     ------------------------------------------------------------------
!     Update wind, temperature and GFL variables.
!     ------------------------------------------------------------------

!     INPUTS
!     ------
!       KPROMA : Horizontal dimension
!       KSTART : Start point
!       KPROF  : End point
!       KFLEV  : Number of levels
!       PDT    : Timestep for physics
!       LDNODRYFLX : true: save total water change rate for CTY next time step
!       PUT9   : t-dt U-wind.
!       PVT9   : t-dt V-wind.
!       PTT9   : t-dt temperature.
!       PGFL   : t-dt GFL variables. 
!       PCONVCTY: vert. MF div. from conv. scheme (unit: s-1)
!     INPUT/OUTPUT
!     ------------
!       PUT1   : t+dt U-wind.
!       PVT1   : t+dt V-wind.
!       PTT1   : t+dt temperature.
!       PGFLT1 : t+dt GFL variables.

!     OUTPUTS
!     -------
!  ------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOECLDP  , ONLY : TECLDP
USE YOETHF   , ONLY : RALVDCP, RALSDCP

!  ------------------------------------------------------------

IMPLICIT NONE

TYPE(TECLDP)      ,INTENT(INOUT) :: YDECLDP
TYPE(TYPE_GFLD)   ,INTENT(IN)    :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDT 
LOGICAL           ,INTENT(IN)    :: LDNODRYFLX
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUT9(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PVT9(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTT9(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(KPROMA,KFLEV,YGFL%NDIM)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGFL_DYN(KPROMA,KFLEV,YGFL%NDIM1)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCONVCTY(KPROMA,KFLEV)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PUT1(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PVT1(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTT1(KPROMA,KFLEV) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(KPROMA,KFLEV,YGFL%NDIM1) 

!  ------------------------------------------------------------

INTEGER(KIND=JPIM) :: JLEV, JROF, JGFL
REAL(KIND=JPRB) :: ZLITOT, ZEPS
REAL(KIND=JPRB) :: ZDWATER_PHY(KPROMA,KFLEV)
REAL(KIND=JPRB) :: ZQWT1(KPROMA,KFLEV)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!  ------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('CPQTUV',0,ZHOOK_HANDLE)
ASSOCIATE(NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, NUMFLDS=>YGFL%NUMFLDS, &
 & YA=>YGFL%YA, YCOMP=>YGFL%YCOMP, YI=>YGFL%YI, YL=>YGFL%YL, YR=>YGFL%YR, &
 & YS=>YGFL%YS,YQ=>YGFL%YQ, &
 & RAMIN=>YDECLDP%RAMIN, RLMIN=>YDECLDP%RLMIN)
!  ------------------------------------------------------------

!*     1.    TIME STEPPING
!            -------------

ZEPS = RLMIN

DO JLEV = 1, KFLEV
  DO JROF = KSTART, KPROF
    PUT1(JROF,JLEV) = PUT9(JROF,JLEV) + PDT * PUT1(JROF,JLEV)
    PVT1(JROF,JLEV) = PVT9(JROF,JLEV) + PDT * PVT1(JROF,JLEV)
    PTT1(JROF,JLEV) = PTT9(JROF,JLEV) + PDT * PTT1(JROF,JLEV)
  ENDDO
ENDDO

DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%LT1) THEN
    DO JLEV = 1, KFLEV
      DO JROF = KSTART, KPROF
        PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1) = PGFL(JROF,JLEV,YCOMP(JGFL)%MP9_PH) +&
         & PDT * PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)  
      ENDDO
    ENDDO
  ENDIF
ENDDO

!---------------------------------------------------------------------------
! Checking and evaporation of small values for the cloud/precip variables
! Corrections have been moved from SLTEND to CPQTUV
!---------------------------------------------------------------------------
! The LNEG attribute could be used to have a more general loop for all GFL.
! But we would need one more attribute to fix the value of RMIN
!---------------------------------------------------------------------------

  DO JLEV=1,KFLEV
    DO JROF=KSTART, KPROF
      IF (YL%LT1) THEN
        IF (PGFLT1(JROF,JLEV,YL%MP1)<ZEPS) THEN
          PGFLT1(JROF,JLEV,YQ%MP1)=PGFLT1(JROF,JLEV,YQ%MP1)+PGFLT1(JROF,JLEV,YL%MP1)
          PTT1(JROF,JLEV)=PTT1(JROF,JLEV)-RALVDCP*PGFLT1(JROF,JLEV,YL%MP1)
          PGFLT1(JROF,JLEV,YL%MP1)=0.0_JPRB
        ENDIF
      ENDIF
      IF (YI%LT1) THEN
        IF (PGFLT1(JROF,JLEV,YI%MP1)<ZEPS) THEN
          PGFLT1(JROF,JLEV,YQ%MP1)=PGFLT1(JROF,JLEV,YQ%MP1)+PGFLT1(JROF,JLEV,YI%MP1)
          PTT1(JROF,JLEV)=PTT1(JROF,JLEV)-RALSDCP*PGFLT1(JROF,JLEV,YI%MP1)
          PGFLT1(JROF,JLEV,YI%MP1)=0.0_JPRB
        ENDIF
      ENDIF
      IF (YR%LT1) THEN
        IF (PGFLT1(JROF,JLEV,YR%MP1)<ZEPS) THEN
          PGFLT1(JROF,JLEV,YQ%MP1)=PGFLT1(JROF,JLEV,YQ%MP1)+PGFLT1(JROF,JLEV,YR%MP1)
          PTT1(JROF,JLEV)=PTT1(JROF,JLEV)-RALVDCP*PGFLT1(JROF,JLEV,YR%MP1)
          PGFLT1(JROF,JLEV,YR%MP1)=0.0_JPRB
        ENDIF
      ENDIF
      IF (YS%LT1) THEN
        IF (PGFLT1(JROF,JLEV,YS%MP1)<ZEPS) THEN
          PGFLT1(JROF,JLEV,YQ%MP1)=PGFLT1(JROF,JLEV,YQ%MP1)+PGFLT1(JROF,JLEV,YS%MP1)
          PTT1(JROF,JLEV)=PTT1(JROF,JLEV)-RALSDCP*PGFLT1(JROF,JLEV,YS%MP1)
          PGFLT1(JROF,JLEV,YS%MP1)=0.0_JPRB
        ENDIF
      ENDIF
      IF (YL%LT1.AND.YI%LT1.AND.YA%LT1) THEN
        ZLITOT= PGFLT1(JROF,JLEV,YL%MP1)+PGFLT1(JROF,JLEV,YI%MP1)
        IF (PGFLT1(JROF,JLEV,YA%MP1)<RAMIN.OR.ZLITOT<ZEPS) THEN 
          PGFLT1(JROF,JLEV,YQ%MP1)=PGFLT1(JROF,JLEV,YQ%MP1)+ZLITOT
          PTT1(JROF,JLEV)=PTT1(JROF,JLEV)-RALVDCP*PGFLT1(JROF,JLEV,YL%MP1) &
         &                               -RALSDCP*PGFLT1(JROF,JLEV,YI%MP1)
          PGFLT1(JROF,JLEV,YL%MP1)=0.0_JPRB
          PGFLT1(JROF,JLEV,YI%MP1)=0.0_JPRB
          PGFLT1(JROF,JLEV,YA%MP1)=0.0_JPRB
        ENDIF
      ENDIF
      IF (YA%LT1) THEN
        IF( PGFLT1(JROF,JLEV,YA%MP1) > 1.0_JPRB ) THEN
          PGFLT1(JROF,JLEV,YA%MP1)=1.0_JPRB
        ENDIF
      ENDIF
    ENDDO
  ENDDO

!!! Include physics net mass changes in CTY equation
! PCONVCTY has the dimension of DIV (s-1)
! PCONVCTY  = 1/rho dM/dz = g dM/dp (M = rho w=massflux)
! PCONVCTY is saved in GFL YPHYCTY to be used in GPCTY at next time step

IF (YGFL%YPHYCTY%LACTIVE) THEN

! If the coupling PHS/DYN for the convection scheme is not activated (RMFADVW=0)
! PCONVCTY = 0

  DO JLEV=1,KFLEV
    DO JROF=KSTART, KPROF
      PGFL(JROF,JLEV,YGFL%YPHYCTY%MP)=PCONVCTY(JROF,JLEV)
    ENDDO
  ENDDO
  ! Continuity equation for specific variables
  ! The mixing ratios have to be modifed to take care of the total mass change
  ! This term is already taken into account for the variables mixed by convection (s,u/v,qv,o3)

!What about ql,qi,qr,qs ?
!  DO JLEV=1,KFLEV
!    DO JROF=KSTART, KPROF
!      IF (YL%LT1) THEN
!      PGFLT1(JROF,JLEV,YL%MP1)= &
! &         PGFLT1(JROF,JLEV,YL%MP1)* &
! &         (1._JPRB-PCONVCTY(JROF,JLEV)*PDT)
!      ENDIF
!    ENDDO
!  ENDDO

!!!! Include change of total water in continuity equation
! total water change rate saved in GFL PHYCTY
  IF (LDNODRYFLX) THEN
  ! Compute total water change due to physics (Delta qt)_phys
  !This option is coded only if all water variable are active
  
    IF (YQ%LACTIVE .AND. YL%LACTIVE .AND. YI%LACTIVE &
 &    .AND. YR%LACTIVE .AND. YS%LACTIVE ) THEN
      DO JLEV=1,KFLEV
        DO JROF=KSTART, KPROF
          ZDWATER_PHY(JROF,JLEV) = &
 &          ( (PGFLT1(JROF,JLEV,YQ%MP1)-PGFL(JROF,JLEV,YQ%MP9_PH)) &
 &             -PGFL_DYN(JROF,JLEV,YQ%MP1)*PDT ) + &
 &          ( ( PGFLT1(JROF,JLEV,YL%MP1)-PGFL(JROF,JLEV,YL%MP9_PH)) &
 &             -PGFL_DYN(JROF,JLEV,YL%MP1)*PDT ) + &
 &          ( ( PGFLT1(JROF,JLEV,YI%MP1) -PGFL(JROF,JLEV,YI%MP9_PH)) &
 &             -PGFL_DYN(JROF,JLEV,YI%MP1)*PDT ) + &
 &          ( ( PGFLT1(JROF,JLEV,YR%MP1)-PGFL(JROF,JLEV,YR%MP9_PH)) &
 &             -PGFL_DYN(JROF,JLEV,YR%MP1)*PDT ) + &
 &          ( ( PGFLT1(JROF,JLEV,YS%MP1)-PGFL(JROF,JLEV,YS%MP9_PH)) &
 &             -PGFL_DYN(JROF,JLEV,YS%MP1)*PDT )

! all prognostic GFL are updated (not sure about other types)
          DO JGFL=1,NUMFLDS
            IF(YCOMP(JGFL)%LT1) THEN
              PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)= &
 &              PGFLT1(JROF,JLEV,YCOMP(JGFL)%MP1)* &
 &              (1._JPRB-ZDWATER_PHY(JROF,JLEV) )
            ENDIF
          ENDDO 
! new total water (shall we use the new or the current one here?)
        ZQWT1(JROF,JLEV) =   PGFLT1(JROF,JLEV,YQ%MP1) + &
 &                         PGFLT1(JROF,JLEV,YL%MP1) + &
 &                         PGFLT1(JROF,JLEV,YI%MP1) + &
 &                         PGFLT1(JROF,JLEV,YR%MP1) + &
 &                         PGFLT1(JROF,JLEV,YS%MP1)
  ! Correction will be done at the next time step
  ! information about total water change at each level saved in PHYCTY
        PGFL(JROF,JLEV,YGFL%YPHYCTY%MP)=PGFL(JROF,JLEV,YGFL%YPHYCTY%MP)+&
&                   ZDWATER_PHY(JROF,JLEV)/(1._JPRB-ZQWT1(JROF,JLEV))/PDT

        ENDDO
      ENDDO
    ENDIF

  ENDIF

ENDIF
!  ------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CPQTUV',1,ZHOOK_HANDLE)
END SUBROUTINE CPQTUV

