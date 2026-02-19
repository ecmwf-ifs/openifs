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

MODULE INTDYNSL_MOD

! Purpose :
! -------
!    To define and compute pointers and logical conditions used when
!    computing local quantities in the dynamics: SL scheme.
!    Allows to use some global structures under CALL_SL
!    (and also their TL and AD).

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    K. YESSAD (CNRM/GMAP)
!    Original : January 2011

! Modifications :
! -------------
!  K. YESSAD (Feb 2014): split into INTDYNSL_MOD
!  S. Malardel Nov 2013: pointers for COMAD weights
!  F. Vana  21-Nov-2017: Option LHOISLT
!  F. Vana  20-Feb-2019: Quintic vertical interpolation
!  F. Vana  18-Jul-2019: SLVF
!-----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDYN   , ONLY : JPSLDIMK

IMPLICIT NONE
SAVE

!=============================================================================

!      1.    TYPE DEFINITIONS
!            ----------------

!      1.11  Types TLSCAW and TRSCAW: pointers for interpolation weights computed in the SL scheme.

! Linear weights:
TYPE TLSCAW
INTEGER(KIND=JPIM) :: M_WDLO        ! distances for horizontal linear interpolations in longitude
INTEGER(KIND=JPIM) :: M_WDLAT       ! distance for horizontal linear interpolations in latitude
INTEGER(KIND=JPIM) :: M_WDVER       ! distance for vertical linear interpolation
INTEGER(KIND=JPIM) :: M_WDLOMAD     ! WDLO for LCOMADH
INTEGER(KIND=JPIM) :: M_WDLAMAD     ! WDLAT for LCOMADH
INTEGER(KIND=JPIM) :: M_WDVERMAD    ! WDVER for LCOMADV
INTEGER(KIND=JPIM) :: NDIM          ! total number of fields allocated

CONTAINS
  
PROCEDURE, PASS :: PRINT => PRINT_TLSCAW_CONFIGURATION 
    
END TYPE TLSCAW

! Other weights:
TYPE TRSCAW
INTEGER(KIND=JPIM) :: M_WCLO(JPSLDIMK)     ! weights for horizontal cubic interpolations in longitude
INTEGER(KIND=JPIM) :: M_WCLA(JPSLDIMK)     ! weights for horizontal cubic interpolations in latitude
INTEGER(KIND=JPIM) :: M_WVINTW             ! vertical cubic interpolation weights
INTEGER(KIND=JPIM) :: M_WCLOSLD(JPSLDIMK)  ! cf. WCLO but for SLHD case
INTEGER(KIND=JPIM) :: M_WCLASLD(JPSLDIMK)  ! cf. WCLA but for SLHD case
INTEGER(KIND=JPIM) :: M_WCLOSLT            ! cf. WCLO 
INTEGER(KIND=JPIM) :: M_WCLASLT            ! cf. WCLA 
INTEGER(KIND=JPIM) :: M_WVINTWSLD          ! cf. WVINTW but for SLHD case
INTEGER(KIND=JPIM) :: M_WVINTWSLT          ! cf. WVINTW 
INTEGER(KIND=JPIM) :: M_WVINTWSLVF         ! cf. WVINTW but for SLVF case 
INTEGER(KIND=JPIM) :: M_WCLOMAD(JPSLDIMK)  ! cf. WCLO but for COMAD case
INTEGER(KIND=JPIM) :: M_WCLAMAD(JPSLDIMK)  ! cf. WCLA but for COMAD case
INTEGER(KIND=JPIM) :: M_WVINTWMAD          ! cf. WVINTW but for COMAD case
INTEGER(KIND=JPIM) :: M_WVINTWS            ! vertical spline interpolation weights
INTEGER(KIND=JPIM) :: M_WVDERW             ! weights to compute vertical derivatives (Hermite cubic vertical interpolation)
INTEGER(KIND=JPIM) :: M_WHVW               ! Hermite vertical cubic interpolation weights
INTEGER(KIND=JPIM) :: M_CW                 ! C_k weights for the vertical WENO interpolation
INTEGER(KIND=JPIM) :: NDIM                 ! total number of fields allocated

CONTAINS
  
PROCEDURE, PASS :: PRINT => PRINT_TRSCAW_CONFIGURATION 

END TYPE TRSCAW

!      1.12  Types TSCO and TCCO: pointers for coordinates computed in the SL scheme.

TYPE TSCO
! spherical geometry:
!   cos(Longitude-Longitude(grid-point))*cos(Latitude) of the interpolation
!   point (geographical longitude and latitude).
! plane projection: x - coordinate (fractional system).
INTEGER(KIND=JPIM) :: M_COSCO
! spherical geometry:
!   sin(Longitude-Longitude(grid-point))*cos(Latitude) of the interpolation
!   point (geographical longitude and latitude).
! plane projection: y - coordinate (fractional system).
INTEGER(KIND=JPIM) :: M_SINCO
! sine of the interpolation point geographical latitude.
INTEGER(KIND=JPIM) :: M_SINLA
! cosine of the geographical angle between the interpolation point
! and the grid-point.
INTEGER(KIND=JPIM) :: M_COPHI
! total number of fields allocated.
INTEGER(KIND=JPIM) :: NDIM

CONTAINS
  
PROCEDURE, PASS :: PRINT => PRINT_TSCO_CONFIGURATION 

END TYPE TSCO

TYPE TCCO
INTEGER(KIND=JPIM) :: M_RLON ! computational sphere longitude of interpolation point
INTEGER(KIND=JPIM) :: M_RLAT ! computational sphere latitude of interpolation point
INTEGER(KIND=JPIM) :: M_RQX  ! first element of the wind displacement matrix (p,q)
INTEGER(KIND=JPIM) :: M_RQY  ! second element of the wind displacement matrix (p,q)
INTEGER(KIND=JPIM) :: NDIM   ! total number of fields allocated

CONTAINS
  
PROCEDURE, PASS :: PRINT => PRINT_TCCO_CONFIGURATION 

END TYPE TCCO

!=============================================================================

!      2.    DECLARATIONS
!            ------------

!      2.11  Types TLSCAW and TRSCAW.

!TYPE(TLSCAW), POINTER :: YYTLSCAW  => NULL()    ! at full levels
!TYPE(TLSCAW), POINTER :: YYTLSCAWH => NULL()    ! at half levels
!TYPE(TRSCAW), POINTER :: YYTRSCAW  => NULL()    ! at full levels
!TYPE(TRSCAW), POINTER :: YYTRSCAWH => NULL()    ! at half levels

!      2.12  Types TSCO and TCCO.

!TYPE(TSCO), POINTER :: YYTSCO => NULL()
!TYPE(TCCO), POINTER :: YYTCCO => NULL()

!=============================================================================

CONTAINS

!      3.    SET-UP

!      3.00  General set-up.

SUBROUTINE SUINTDYNSL(YDDYN,YDDYNA,YGFL,YDTLSCAW,YDTLSCAWH,YDTRSCAW,YDTRSCAWH,YDTSCO,YDTCCO)

!--------------------------------------------------------------------------
! Sets-up internal dynamics structures used under CALL_SL
! Must be called after SUGFL3.
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

USE YOMDYN   , ONLY : TDYN
USE YOMDYNA  , ONLY : TDYNA
USE YOM_YGFL , ONLY : TYPE_GFLD

TYPE(TDYN)      , INTENT(INOUT) :: YDDYN
TYPE(TDYNA)     , INTENT(INOUT) :: YDDYNA
TYPE(TYPE_GFLD) , INTENT(INOUT) :: YGFL
TYPE(TLSCAW)    , INTENT(INOUT) :: YDTLSCAW   ! at full levels
TYPE(TLSCAW)    , INTENT(INOUT) :: YDTLSCAWH  ! at half levels
TYPE(TRSCAW)    , INTENT(INOUT) :: YDTRSCAW   ! at full levels
TYPE(TRSCAW)    , INTENT(INOUT) :: YDTRSCAWH  ! at half levels
TYPE(TSCO)      , INTENT(INOUT) :: YDTSCO 
TYPE(TCCO)      , INTENT(INOUT) :: YDTCCO 

LOGICAL :: LLHVI,LLVINTWS
INTEGER(KIND=JPIM) :: JGFL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUINTDYNSL',0,ZHOOK_HANDLE)
ASSOCIATE(NUMFLDS=>YGFL%NUMFLDS, YCOMP=>YGFL%YCOMP, &
 & NSLDIMK=>YDDYN%NSLDIMK)
!--------------------------------------------------------------------------

! Types TLSCAW and TRSCAW:
CALL SUPTR_TLSCAW(YDTLSCAW)
CALL SUPTR_TLSCAW(YDTLSCAWH)
LLHVI=.FALSE.
DO JGFL=1,NUMFLDS
  IF(YCOMP(JGFL)%CSLINT == 'LAIHVT      ' .OR. &
     & YCOMP(JGFL)%CSLINT == 'LAIHVTQM    ' .OR. &
     & YCOMP(JGFL)%CSLINT == 'LAIHVTQMH   ' ) THEN
    LLHVI=.TRUE.
  ENDIF
ENDDO
LLVINTWS=YDDYNA%LVSPLIP
CALL SUPTR_TRSCAW(YDDYN,YDDYNA,NSLDIMK,LLVINTWS,LLHVI,YDTRSCAW)
CALL SUPTR_TRSCAW(YDDYN,YDDYNA,NSLDIMK,.FALSE.,.FALSE.,YDTRSCAWH)

! Types TSCO and TCCO:
CALL SUPTR_TSCO(YDTSCO)
CALL SUPTR_TCCO(YDTCCO)

!--------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUINTDYNSL',1,ZHOOK_HANDLE)
END SUBROUTINE SUINTDYNSL


!      3.11  Set-up for types TLSCAW and TRSCAW.

SUBROUTINE SUPTR_TLSCAW(YD)

!--------------------------------------------------------------------------
! Sets-up pointers for TLSCAW structure

! YD       : contains pointers.
!--------------------------------------------------------------------------

TYPE(TLSCAW),    INTENT(OUT)   :: YD

!--------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INCR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TLSCAW',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

INCR=1
YD%M_WDLO=INCR
INCR=INCR+4
YD%M_WDLAT=INCR
INCR=INCR+1
YD%M_WDVER=INCR
INCR=INCR+1
YD%M_WDLOMAD=INCR
INCR=INCR+4
YD%M_WDLAMAD=INCR
INCR=INCR+1
YD%M_WDVERMAD=INCR
INCR=INCR+1

YD%NDIM=INCR-1

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TLSCAW',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTR_TLSCAW

SUBROUTINE SUPTR_TRSCAW(YDDYN,YDDYNA,KSLDIMK,LDVINTWS,LDHVI,YD)

!--------------------------------------------------------------------------
! Sets-up pointers for TRSCAW structure

! KSLDIMK  : number of sets of horizontal non-linear weights
! LDVINTWS : T if 'VINTWS' active.
! LDHVI    : T if 'VDERW' and 'HVW' active.
! YD       : contains pointers.
!--------------------------------------------------------------------------

USE YOMDYN , ONLY : TDYN
USE YOMDYNA, ONLY : TDYNA
TYPE(TDYN)        , INTENT(INOUT) :: YDDYN
TYPE(TDYNA)       , INTENT(INOUT) :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)    :: KSLDIMK
LOGICAL,            INTENT(IN)    :: LDVINTWS
LOGICAL,            INTENT(IN)    :: LDHVI
TYPE(TRSCAW),       INTENT(OUT)   :: YD

!--------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: INCR,JS
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TRSCAW',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

INCR=1
DO JS=1,KSLDIMK
  YD%M_WCLO(JS)=INCR
  INCR=INCR+6
ENDDO
DO JS=1,KSLDIMK
  YD%M_WCLA(JS)=INCR
  INCR=INCR+3
ENDDO
YD%M_WVINTW=INCR
INCR=INCR+3
IF (YDDYNA%LSLTVWENO.OR.YDDYNA%LRHSVWENO) INCR=INCR+6
IF (YDDYNA%LSLVF) THEN
  YD%M_WVINTWSLVF=INCR
  INCR=INCR+3
  ! SLVF is only applied to temperature
  IF (YDDYNA%LRHSVWENO.AND.YDDYN%LVWENO_T) INCR=INCR+6
ELSE
  YD%M_WVINTWSLVF=YD%M_WVINTW
ENDIF
DO JS=1,KSLDIMK
  YD%M_WCLOSLD(JS)=INCR
  INCR=INCR+6
ENDDO
IF (YDDYN%LSLHDHEAT) THEN
  YD%M_WCLOSLT=INCR
  INCR=INCR+6
ELSE
  YD%M_WCLOSLT=YD%M_WCLOSLD(1)
ENDIF
DO JS=1,KSLDIMK
  YD%M_WCLASLD(JS)=INCR
  INCR=INCR+3
ENDDO
IF (YDDYN%LSLHDHEAT) THEN
  YD%M_WCLASLT=INCR
  INCR=INCR+3
ELSE
  YD%M_WCLASLT=YD%M_WCLASLD(1)
ENDIF
YD%M_WVINTWSLD=INCR
INCR=INCR+3
IF(YDDYN%LSLHDHEAT) THEN
  YD%M_WVINTWSLT=INCR
  INCR=INCR+3
ELSE
  YD%M_WVINTWSLT=YD%M_WVINTWSLD
ENDIF
DO JS=1,KSLDIMK
  YD%M_WCLOMAD(JS)=INCR
  INCR=INCR+6
ENDDO
DO JS=1,KSLDIMK
  YD%M_WCLAMAD(JS)=INCR
  INCR=INCR+3
ENDDO
YD%M_WVINTWMAD=INCR
INCR=INCR+3
IF (LDVINTWS) THEN
  YD%M_WVINTWS=INCR
  INCR=INCR+4
ELSE
  YD%M_WVINTWS=1
ENDIF
IF (LDHVI) THEN
  YD%M_WVDERW=INCR
  INCR=INCR+4
ELSE
  YD%M_WVDERW=1
ENDIF
IF (LDHVI) THEN
  YD%M_WHVW=INCR
  INCR=INCR+4
ELSE
  YD%M_WHVW=1
ENDIF
IF (YDDYNA%LSLTVWENO.OR.YDDYNA%LRHSVWENO) THEN
  YD%M_CW=INCR
  INCR=INCR+3
ELSE
  YD%M_CW=1
ENDIF
YD%NDIM=INCR-1

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TRSCAW',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTR_TRSCAW

!      3.12  Set-up for types TSCO and TCCO.

SUBROUTINE SUPTR_TSCO(YD)

!--------------------------------------------------------------------------
! Sets-up pointers for TSCO structure

! YD       : contains pointers.
!--------------------------------------------------------------------------

TYPE(TSCO),  INTENT(OUT)   :: YD

!--------------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TSCO',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

YD%M_COSCO=1
YD%M_SINCO=2
YD%M_SINLA=3
YD%M_COPHI=4
YD%NDIM=4

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TSCO',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTR_TSCO

SUBROUTINE SUPTR_TCCO(YD)

!--------------------------------------------------------------------------
! Sets-up pointers for TCCO structure

! YD       : contains pointers.
!--------------------------------------------------------------------------

TYPE(TCCO), INTENT(OUT)   :: YD

!--------------------------------------------------------------------------

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TCCO',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

YD%M_RLON=1
YD%M_RLAT=2
YD%M_RQX=3
YD%M_RQY=4
YD%NDIM=4

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:SUPTR_TCCO',1,ZHOOK_HANDLE)
END SUBROUTINE SUPTR_TCCO

!=============================================================================
SUBROUTINE PRINT_TLSCAW_CONFIGURATION(SELF, KDEPTH, KOUTNO, CNAME)
  IMPLICIT NONE
  CLASS(TLSCAW)   , INTENT(IN) :: SELF
  INTEGER(KIND=JPIM)         , INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM)         , INTENT(IN) :: KOUTNO
  CHARACTER(LEN=*), INTENT(IN) :: CNAME

  INTEGER(KIND=JPIM) :: IDEPTHLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

  IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TLSCAW_CONFIGURATION',0,ZHOOK_HANDLE)

  IDEPTHLOC = KDEPTH+2

  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%' // CNAME // ' : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDLO = ', SELF%M_WDLO    
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDLAT = ', SELF%M_WDLAT   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDVER = ', SELF%M_WDVER   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDLOMAD = ', SELF%M_WDLOMAD 
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDLAMAD = ', SELF%M_WDLAMAD 
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WDVERMAD = ', SELF%M_WDVERMAD
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIM = ', SELF%NDIM      

  IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TLSCAW_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_TLSCAW_CONFIGURATION
!=============================================================================
SUBROUTINE PRINT_TRSCAW_CONFIGURATION(SELF, KDEPTH, KOUTNO, CNAME)
  IMPLICIT NONE
  CLASS(TRSCAW)   , INTENT(IN) :: SELF
  INTEGER(KIND=JPIM)         , INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM)         , INTENT(IN) :: KOUTNO
  CHARACTER(LEN=*), INTENT(IN) :: CNAME

  INTEGER(KIND=JPIM) :: IDEPTHLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TRSCAW_CONFIGURATION',0,ZHOOK_HANDLE)
  IDEPTHLOC = KDEPTH+2
  
!!  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIM = ', SELF%NDIM      
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) //'SUM(M_WCLO) = ', SUM(SELF%M_WCLO)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SUM(M_WCLA) = ', SUM(SELF%M_WCLA)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTW = ',  SELF%M_WVINTW
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SUM(M_WCLOSLD) = ', SUM(SELF%M_WCLOSLD)  
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SUM(M_WCLASLD) = ', SUM(SELF%M_WCLASLD)  
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WCLOSLT = ', SELF%M_WCLOSLT            
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WCLASLT = ', SELF%M_WCLASLT           
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTWSLD = ', SELF%M_WVINTWSLD
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTWSLT = ', SELF%M_WVINTWSLT   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTWSLVF = ', SELF%M_WVINTWSLVF   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SUM(M_WCLOMAD) = ', SUM(SELF%M_WCLOMAD)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'SUM(M_WCLAMAD) = ', SUM(SELF%M_WCLAMAD)
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTWMAD = ', SELF%M_WVINTWMAD
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVINTWS = ', SELF%M_WVINTWS   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WVDERW = ', SELF%M_WVDERW     
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_WHVW = ', SELF%M_WHVW
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_CW = ', SELF%M_CW
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIM = ', SELF%NDIM
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TRSCAW_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_TRSCAW_CONFIGURATION
!=============================================================================
SUBROUTINE PRINT_TSCO_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TSCO) , INTENT(IN) :: SELF
  INTEGER(KIND=JPIM)     , INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM)     , INTENT(IN) :: KOUTNO

  INTEGER(KIND=JPIM) :: IDEPTHLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TSCO_CONFIGURATION',0,ZHOOK_HANDLE)
  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%yytsco : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_COSCO = ', SELF%M_COSCO
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_SINCO = ', SELF%M_SINCO   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_SINLA = ', SELF%M_SINLA   
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_COPHI = ', SELF%M_COPHI
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIM = ', SELF%NDIM      
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TSCO_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_TSCO_CONFIGURATION
!=============================================================================
SUBROUTINE PRINT_TCCO_CONFIGURATION(SELF, KDEPTH, KOUTNO)
  IMPLICIT NONE
  CLASS(TCCO) , INTENT(IN) :: SELF
  INTEGER(KIND=JPIM)     , INTENT(IN) :: KDEPTH
  INTEGER(KIND=JPIM)     , INTENT(IN) :: KOUTNO

  INTEGER(KIND=JPIM) :: IDEPTHLOC
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TCCO_CONFIGURATION',0,ZHOOK_HANDLE)
  IDEPTHLOC = KDEPTH+2
  
  WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%yytcco : '
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RLON = ',SELF%M_RLON
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RLAT = ',SELF%M_RLAT
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RQX = ',SELF%M_RQX
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'M_RQY = ',SELF%M_RQY
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'NDIM = ',SELF%NDIM
IF (LHOOK) CALL DR_HOOK('INTDYNSL_MOD:PRINT_TCCO_CONFIGURATION',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_TCCO_CONFIGURATION
!=============================================================================

END MODULE INTDYNSL_MOD
