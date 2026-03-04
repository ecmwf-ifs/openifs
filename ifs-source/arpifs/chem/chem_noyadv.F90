! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOCHECK
#endif

SUBROUTINE CHEM_NOYADV(YDGEOMETRY,YGFL,YDCHEM, KMODE,KCHEM_NOXADV,PGFL,PGFLT1)
!
!     Purpose.   fill family array before advection (KMODE=1) 
!                fill individual components after advection (KMODE=2)                
!     --------   
!
!*   Interface.
!     ----------
!
!        *CALL* *CHEM_NOYADV* (KMODE, PGFL, PGFLT1 )
!
!     Explicit arguments :
!        --------------------
!
!
!!     INPUT/OUTPUT:
!     -------------
!        KCHEM_NOXADV : Switch defining which family to fill: NOx or NOy and/or Clx 
!        PGFL         : GFL variables buffer t0 t1
!        PGFLT1       : GFL variables buffer t + dt
!
!        Implicit arguments :  None.
!        --------------------
!
!     Method.
!     -------
!     - Fill family tracer for Stratosphere before advection 
!     - (apply advection of individual family components and family tracer)
!     - overwrite family components with family tracer above specified pressure level
!     - Make sure to switch off global mass conservation for family tracer
!
!     Externals.   See includes below.
!     ----------
!
!
!     Author.
!     -------
!        Johannes Flemming   *ECMWF*
!        Vincent Huijnen  *KNMI* (Taken from chem_noxadv.F90)
!
! Modifications
! -------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCHEM  , ONLY : TCHEM

IMPLICIT NONE

! arguments

TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
TYPE(TCHEM)       ,INTENT(INOUT) :: YDCHEM
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFL(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM,YDGEOMETRY%YRDIM%NGPBLKS) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PGFLT1(YDGEOMETRY%YRDIM%NPROMA,YDGEOMETRY%YRDIMV%NFLEVG,YGFL%NDIM1, &
 & YDGEOMETRY%YRDIM%NGPBLKS) 
INTEGER(KIND=JPIM), INTENT(IN) :: KMODE 
INTEGER(KIND=JPIM), INTENT(IN) :: KCHEM_NOXADV 

! local arrays
INTEGER(KIND=JPIM) :: JGFL ,IST, IBL, JROF, JLEV, ICEND, IEND, JKGLO, JLEV_TROP,JTRAC
INTEGER(KIND=JPIM) :: ICLXA,IBRXA,INOXA,INOYA
REAL(KIND=JPRB)    :: ZFAMRAT, ZFAMA 
REAL(KIND=JPRB)    :: ZPTROP, ZPLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Local definition of generic family tracer
INTEGER(KIND=JPIM),DIMENSION(20) :: ICOMP
REAL(KIND=JPRB),   DIMENSION(20) :: ZRMCOMP
INTEGER(KIND=JPIM) :: INJ, IFAM

! Local definition of Clx 
 INTEGER(KIND=JPIM):: INJ_CLX 
 INTEGER(KIND=JPIM), PARAMETER :: INJ_CLXMAX = 12
 CHARACTER(LEN=6), DIMENSION(INJ_CLXMAX):: CLCOMP_CLX
! With corresponding number of Cl-atoms:
REAL(KIND=JPRB),DIMENSION(INJ_CLXMAX)::ZNCL_COMP_CLX 
REAL(KIND=JPRB),DIMENSION(INJ_CLXMAX) :: ZRMCOMP_CLX
INTEGER(KIND=JPIM),DIMENSION(INJ_CLXMAX) :: ICOMP_CLX


! Local definition of Brx 
 INTEGER(KIND=JPIM):: INJ_BRX 
 INTEGER(KIND=JPIM), PARAMETER :: INJ_BRXMAX = 8
 CHARACTER(LEN=6), DIMENSION(INJ_BRXMAX):: CLCOMP_BRX
! With corresponding number of BR-atoms:
REAL(KIND=JPRB),DIMENSION(INJ_BRXMAX)::ZNCL_COMP_BRX 
REAL(KIND=JPRB),DIMENSION(INJ_BRXMAX) :: ZRMCOMP_BRX
INTEGER(KIND=JPIM),DIMENSION(INJ_BRXMAX) :: ICOMP_BRX





! Local definition of NOy 
INTEGER(KIND=JPIM) :: INJ_NOY
INTEGER(KIND=JPIM), PARAMETER :: INJ_NOYMAX = 10
CHARACTER(LEN=6), DIMENSION(INJ_NOYMAX):: CLCOMP_NOY
! With corresponding number of N-atoms:
REAL(KIND=JPRB),DIMENSION(INJ_NOYMAX) :: ZNCL_COMP_NOY 
REAL(KIND=JPRB),DIMENSION(INJ_NOYMAX) :: ZRMCOMP_NOY
INTEGER(KIND=JPIM),DIMENSION(INJ_NOYMAX) :: ICOMP_NOY


! Local simple definition of NOx  - use is depreciated
INTEGER(KIND=JPIM), PARAMETER :: INJ_NOX = 2
CHARACTER(LEN=6), DIMENSION(INJ_NOX), PARAMETER :: &
& CLCOMP_NOX=(/'NO    ','NO2   ' /)
! With corresponding number of N-atoms:
REAL(KIND=JPRB),DIMENSION(INJ_NOX),PARAMETER  :: &
& ZNCL_COMP_NOX =(/1.,1./)
REAL(KIND=JPRB),DIMENSION(INJ_NOX) :: ZRMCOMP_NOX
INTEGER(KIND=JPIM),DIMENSION(INJ_NOX) :: ICOMP_NOX


IF (LHOOK) CALL DR_HOOK('CHEM_NOYADV',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP,YDVAB=>YDGEOMETRY%YRVAB)
ASSOCIATE(NCHEM=>YGFL%NCHEM, NDIM=>YGFL%NDIM, NDIM1=>YGFL%NDIM1, &
 & YCHEM=>YGFL%YCHEM, YDGSGEOM=>YDGEOMETRY%YRGSGEOM,&
 & NGPBLKS=>YDDIM%NGPBLKS, NPROMA=>YDDIM%NPROMA, &
 & NFLEVG=>YDDIMV%NFLEVG, &
 & NGPTOT=>YDGEM%NGPTOT)

!* find index for (stratospheric!) CLX components and CLX tracer, put into 'CLXA'
ICOMP(:) = -9999
ICOMP_CLX(:) = -9999
ICOMP_BRX(:) = -9999
ICOMP_NOY(:) = -9999
ICOMP_NOX(:) = -9999
ICLXA =-9999
IBRXA =-9999
INOYA =-9999
INOXA =-9999


SELECT CASE (TRIM(YDCHEM%CHEM_SCHEME))

  CASE ("mocage")         

    ! Initialize Chlorine family
    ! INJ_CLX should be lower than INJ_CLXMAX !
    INJ_CLX = 10
    CLCOMP_CLX(1:INJ_CLX)=(/ 'CL2O2 ','OCLO  ','BRCL  ','HOCL  ','CLONO2','CL    ','HCL   ','CLO   ','CLNO2 ','CL2   ' /)
    ZNCL_COMP_CLX(1:INJ_CLX) =(/2.,1.,1.,1.,1.,1.,1.,1.,1.,2./)

    ! initialize bromine family    
    INJ_BRX=6
    CLCOMP_BRX(1:INJ_BRX)=(/ 'BRCL  ','HOBR  ','BRONO2','BR    ','HBR   ','BRO   ' /)
    ! With corresponding number of Cl-atoms:
    ZNCL_COMP_BRx(1:INJ_BRX) =(/1.,1.,1.,1.,1.,1./)

    ! Initialize NOy family
    INJ_NOY = 10
    CLCOMP_NOY(1:INJ_NOY)=(/ 'N     ','NO    ','NO2   ','NO3   ','HO2NO2','N2O5  ','HNO3  ','CLNO2 ','CLONO2','BRONO2' /)
    ZNCL_COMP_NOY(1:INJ_NOY) =(/1.,1.,1.,1.,1.,2.,1.,1.,1.,1./)

  CASE ("mozart")         

    ! Initialize Chlorine family
    ! INJ_CLX should be lower than INJ_CLXMAX !
    INJ_CLX = 9
    CLCOMP_CLX(1:INJ_CLX)=(/ 'CL2O2 ','OCLO  ','BRCL  ','HOCL  ','CLONO2','CL    ','HCL   ','CLO   ','CL2   ' /)
    ! With corresponding number of Cl-atoms:
    ZNCL_COMP_CLX(1:INJ_CLX) =(/2.,1.,1.,1.,1.,1.,1.,1.,2./)

    ! initialize bromine family    
    INJ_BRX=6
    CLCOMP_BRX(1:INJ_BRX)=(/ 'BRCL  ','HOBR  ','BRONO2','BR    ','HBR   ','BRO   ' /)
    ! With corresponding number of Br-atoms:
    ZNCL_COMP_BRx(1:INJ_BRX) =(/1.,1.,1.,1.,1.,1./)

    INJ_NOY = 9
    CLCOMP_NOY(1:INJ_NOY)=(/ 'N     ','NO    ','NO2   ','NO3   ','HO2NO2','N2O5  ','HNO3  ','CLONO2','BRONO2' /)
    ! With corresponding number of N-atoms:
    ZNCL_COMP_NOY(1:INJ_NOY) =(/1.,1.,1.,1.,1.,2.,1.,1.,1./)

  CASE ("bascoe" , "bascoetm5")         

    ! Initialize Chlorine family
    ! INJ_CLX should be lower than INJ_CLXMAX !
    INJ_CLX = 11
    CLCOMP_CLX(1:INJ_CLX)=(/ 'CL2O2 ','OCLO  ','BRCL  ','HOCL  ','CLONO2','CL    ','HCL   ','CLO   ','CLNO2 ','CL2   ','CLOO  ' /)
    ! With corresponding number of Cl-atoms:
    ZNCL_COMP_CLX(1:INJ_CLX) =(/2.,1.,1.,1.,1.,1.,1.,1.,1.,2.,1./)

    ! initialize bromine family    
    INJ_BRX=7
    CLCOMP_BRX(1:INJ_BRX)=(/ 'BRCL  ','HOBR  ','BRONO2','BR    ','HBR   ','BRO   ','BR2   ' /)
    ! With corresponding number of Br-atoms:
    ZNCL_COMP_BRx(1:INJ_BRX) =(/1.,1.,1.,1.,1.,1.,2./)

    ! Initialize NOy family
    INJ_NOY = 10
    CLCOMP_NOY(1:INJ_NOY)=(/ 'N     ','NO    ','NO2   ','NO3   ','HO2NO2','N2O5  ','HNO3  ','CLNO2 ','CLONO2','BRONO2' /)
    ZNCL_COMP_NOY(1:INJ_NOY) =(/1.,1.,1.,1.,1.,2.,1.,1.,1.,1./)

  CASE DEFAULT
      
    CALL ABOR1(" CHEMISTRY SCHEME "//TRIM(YDCHEM%CHEM_SCHEME)//" NOT SUPPORTED" )   
        
END SELECT 



IF (KCHEM_NOXADV == 1_JPIM) THEN
! Initialization 'classic' NOx-advection option - depreciated
  DO JGFL = 1, NCHEM
    DO JTRAC=1,INJ_NOX
      IF (TRIM(YCHEM(JGFL)%CNAME) == TRIM(CLCOMP_NOX(JTRAC)) )    ICOMP_NOX(JTRAC)=JGFL
    ENDDO
    ! Use NOXA-tracer to fill the stratospheric NOX field
    IF (TRIM(YCHEM(JGFL)%CNAME) == 'NOXA' )    INOXA=JGFL
  ENDDO
  
  IF (ANY( ICOMP_NOX    == -9999).OR. INOXA == -9999) THEN
     CALL ABOR1(' NOX indiv comps and/or NOXA not defined ')
  ENDIF 
  
  DO JTRAC=1,INJ_NOX
    ZRMCOMP_NOX(JTRAC) = ZNCL_COMP_NOX(JTRAC) / YCHEM(ICOMP_NOX(JTRAC))%RMOLMASS
  ENDDO

ELSEIF (KCHEM_NOXADV > 1_JPIM) THEN
! Alternatively Initialization 'revised' NOy-advection option 
  DO JGFL = 1, NCHEM
    DO JTRAC=1,INJ_NOY
      IF (TRIM(YCHEM(JGFL)%CNAME) == TRIM(CLCOMP_NOY(JTRAC)) )    ICOMP_NOY(JTRAC)=JGFL
    ENDDO
    ! Use NOXA-tracer to fill the stratospheric NOY field
    IF (TRIM(YCHEM(JGFL)%CNAME) == 'NOXA' )    INOYA=JGFL
  ENDDO
  
  IF (ANY( ICOMP_NOY(1:INJ_NOY)    == -9999).OR. INOYA == -9999) THEN
     CALL ABOR1(' NOY indiv comps and/or NOXA not defined ')
  ENDIF 
  
  DO JTRAC=1,INJ_NOY
    ZRMCOMP_NOY(JTRAC) = ZNCL_COMP_NOY(JTRAC) / YCHEM(ICOMP_NOY(JTRAC))%RMOLMASS
  ENDDO

ENDIF

IF (KCHEM_NOXADV == 3_JPIM .OR. KCHEM_NOXADV==4) THEN
! Additionally fill the CLx-family tracer (if available)

  DO JGFL = 1, NCHEM
    DO JTRAC=1,INJ_CLX
      IF (TRIM(YCHEM(JGFL)%CNAME) == TRIM(CLCOMP_CLX(JTRAC)) )    ICOMP_CLX(JTRAC)=JGFL
    ENDDO
    ! Use CLXA-tracer to fill the stratospheric CLX field
    IF (TRIM(YCHEM(JGFL)%CNAME) == 'CLXA' .OR. TRIM(YCHEM(JGFL)%CNAME) == 'CLY')    ICLXA=JGFL
  ENDDO
  
  IF (ANY( ICOMP_CLX(1:INJ_CLX)    == -9999).OR. ICLXA == -9999) THEN
     CALL ABOR1(' CLX indiv comps and/or CLXA not defined ')
  ENDIF 
  
  DO JTRAC=1,INJ_CLX
    ZRMCOMP_CLX(JTRAC) = ZNCL_COMP_CLX(JTRAC) / YCHEM(ICOMP_CLX(JTRAC))%RMOLMASS
  ENDDO

ENDIF


IF (KCHEM_NOXADV == 4_JPIM) THEN
! Additionally fill the BRx-family tracer (if available)

  DO JGFL = 1, NCHEM
    DO JTRAC=1,INJ_BRX
      IF (TRIM(YCHEM(JGFL)%CNAME) == TRIM(CLCOMP_BRX(JTRAC)) )    ICOMP_BRX(JTRAC)=JGFL
    ENDDO
    ! Use BRXA-tracer to fill the stratospheric BRX field. Use 'BRXA' or 'BRY?'
    IF (TRIM(YCHEM(JGFL)%CNAME) == 'BRXA' .OR. TRIM(YCHEM(JGFL)%CNAME) == 'BRY')    IBRXA=JGFL
  ENDDO
  
  IF (ANY( ICOMP_BRX(1:INJ_BRX) == -9999).OR. IBRXA == -9999) THEN
     CALL ABOR1(' BRX indiv comps and/or BRXA not defined ')
  ENDIF 
  
  DO JTRAC=1,INJ_BRX
    ZRMCOMP_BRX(JTRAC) = ZNCL_COMP_BRX(JTRAC) / YCHEM(ICOMP_BRX(JTRAC))%RMOLMASS
  ENDDO

ENDIF
! ZPTROP: lower pressure level where NOXA tracer should be filled
! ZPTROP = 32000_JPRB 

! fill CLXS array
IF (KMODE == 1 ) THEN 


  IF (KCHEM_NOXADV == 1_JPIM) THEN
    ICOMP(:) = -9999
    INJ=INJ_NOX
    IFAM=INOXA
    ICOMP(1:INJ_NOX)=ICOMP_NOX(1:INJ_NOX)
    ZRMCOMP(1:INJ_NOX)=ZRMCOMP_NOX(1:INJ_NOX)
  ELSE ! KCHEM_NOXADV > 1: use NOy
    ICOMP(:) = -9999
    INJ=INJ_NOY
    IFAM=INOYA
    ICOMP(1:INJ_NOY)=ICOMP_NOY(1:INJ_NOY)
    ZRMCOMP(1:INJ_NOY)=ZRMCOMP_NOY(1:INJ_NOY)
  ENDIF

  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    IST=1
    IEND=ICEND

    ! Fill FAM mass globally
    DO JROF=IST,IEND
      DO JLEV=1,NFLEVG
        PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  = 0._JPRB
        DO JTRAC=1,INJ
          PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  =  PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)+ &
           &                    ZRMCOMP(JTRAC) * PGFL(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP,IBL)   
        ENDDO
        PGFL(JROF,JLEV,YCHEM(IFAM)%MP9,IBL) = PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL) 
      ENDDO

    ENDDO
  ENDDO

  IF (KCHEM_NOXADV == 3_JPIM .OR. KCHEM_NOXADV==4_JPIM ) THEN
    ICOMP(:) = -9999
    INJ=INJ_CLX
    IFAM=ICLXA
    ICOMP(1:INJ_CLX)=ICOMP_CLX(1:INJ_CLX)
    ZRMCOMP(1:INJ_CLX)=ZRMCOMP_CLX(1:INJ_CLX)

    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
  
      ! Fill FAM mass globally for CLx
      DO JROF=IST,IEND
        DO JLEV=1,NFLEVG
          PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  = 0._JPRB
          DO JTRAC=1,INJ
            PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  =  PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)+ &
             &                    ZRMCOMP(JTRAC) * PGFL(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP,IBL)   
          ENDDO
          PGFL(JROF,JLEV,YCHEM(IFAM)%MP9,IBL) = PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL) 
        ENDDO
  
      ENDDO
    ENDDO

  ENDIF

  IF (KCHEM_NOXADV == 4_JPIM) THEN
    ICOMP(:) = -9999
    INJ=INJ_BRX
    IFAM=IBRXA
    ICOMP(1:INJ_BRX)=ICOMP_BRX(1:INJ_BRX)
    ZRMCOMP(1:INJ_BRX)=ZRMCOMP_BRX(1:INJ_BRX)

    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
  
      ! Fill FAM mass globally for BRX
      DO JROF=IST,IEND
        DO JLEV=1,NFLEVG
          PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  = 0._JPRB
          DO JTRAC=1,INJ
            PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)  =  PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL)+ &
             &                    ZRMCOMP(JTRAC) * PGFL(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP,IBL)   
          ENDDO
          PGFL(JROF,JLEV,YCHEM(IFAM)%MP9,IBL) = PGFL(JROF,JLEV,YCHEM(IFAM)%MP,IBL) 
        ENDDO
  
      ENDDO
    ENDDO

  ENDIF


! KMODE = 2
ELSE

  IF (KCHEM_NOXADV == 1_JPIM) THEN
    ICOMP(:) = -9999
    INJ=INJ_NOX
    IFAM=INOXA
    ICOMP(1:INJ_NOX)=ICOMP_NOX(1:INJ_NOX)
    ZRMCOMP(1:INJ_NOX)=ZRMCOMP_NOX(1:INJ_NOX)

  ELSE ! KCHEM_NOXADV > 1: use NOy
    ICOMP(:) = -9999
    INJ=INJ_NOY
    IFAM=INOYA
    ICOMP(1:INJ_NOY)=ICOMP_NOY(1:INJ_NOY)
    ZRMCOMP(1:INJ_NOY)=ZRMCOMP_NOY(1:INJ_NOY)
  ENDIF

  DO JKGLO=1,NGPTOT,NPROMA
    ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
    IBL=(JKGLO-1)/NPROMA+1
    IST=1
    IEND=ICEND

    DO JROF=IST,IEND

      ! Find index for simple tropopause level definition...
      ! Latitude: YDGSGEOM(IBL)%GELAT(JROF)
      ZPTROP=24000_JPRB - 14800_JPRB * (COS(YDGSGEOM(IBL)%GELAT(JROF)))**4_JPIM
      DO JLEV = 1,NFLEVG
        ZPLEV = YDVAB%VAH(JLEV) + YDVAB%VBH(JLEV) * 101300.0_JPRB 
        IF (ZPLEV< ZPTROP) THEN
          JLEV_TROP = JLEV
        ELSE
          ! level found; exit this loop
          CYCLE
        ENDIF
      ENDDO
!Only stratospheric levels to be updated - remove two levels to be sure.
      DO JLEV=1,JLEV_TROP-2
        ZFAMA   =0._JPRB
        ! new Family tracer
        DO JTRAC=1,INJ
          ZFAMA = ZFAMA + ZRMCOMP(JTRAC) * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
        ENDDO
        ZFAMRAT = PGFLT1(JROF,JLEV,YCHEM(IFAM)%MP1,IBL) / ZFAMA 

! make sure that individual NOX/NOY components after SL equals family after SL   
        DO JTRAC=1,INJ
          PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)   = ZFAMRAT * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  IF (KCHEM_NOXADV == 3_JPIM .OR. KCHEM_NOXADV == 4_JPIM ) THEN
  ! Additionally the Same procedure for CLx-family advection 
    ICOMP(:) = -9999
    INJ=INJ_CLX
    IFAM=ICLXA
    ICOMP(1:INJ_CLX)=ICOMP_CLX(1:INJ_CLX)
    ZRMCOMP(1:INJ_CLX)=ZRMCOMP_CLX(1:INJ_CLX)
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
  
      DO JROF=IST,IEND
  
        ! Find index for simple tropopause level definition...
        ! Latitude: YDGSGEOM(IBL)%GELAT(JROF)
        ZPTROP=24000_JPRB - 14800_JPRB * (COS(YDGSGEOM(IBL)%GELAT(JROF)))**4_JPIM
        DO JLEV = 1,NFLEVG
          ZPLEV = YDVAB%VAH(JLEV) + YDVAB%VBH(JLEV) * 101300.0_JPRB 
          IF (ZPLEV< ZPTROP) THEN
            JLEV_TROP = JLEV
          ELSE
            ! level found; exit this loop
            CYCLE
          ENDIF
        ENDDO
!Only stratospheric levels to be updated - remove two levels to be sure.
        DO JLEV=1,JLEV_TROP-2
          ZFAMA   =0._JPRB
          !new Family tracer
          DO JTRAC=1,INJ
            ZFAMA = ZFAMA + ZRMCOMP(JTRAC) * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
          ENDDO
          ZFAMRAT = PGFLT1(JROF,JLEV,YCHEM(IFAM)%MP1,IBL) / ZFAMA 
  
! make sure that individual CLx components after SL equals family after SL   
          DO JTRAC=1,INJ
            PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)   = ZFAMRAT * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDIF ! KCHEM_NOXADV=3
  

  IF ( KCHEM_NOXADV == 4_JPIM) THEN
  ! Additionally the Same procedure for BRX-family advection 
    ICOMP(:) = -9999
    INJ=INJ_BRX
    IFAM=IBRXA
    ICOMP(1:INJ_BRX)=ICOMP_BRX(1:INJ_BRX)
    ZRMCOMP(1:INJ_BRX)=ZRMCOMP_BRX(1:INJ_BRX)
    DO JKGLO=1,NGPTOT,NPROMA
      ICEND=MIN(NPROMA,NGPTOT-JKGLO+1)
      IBL=(JKGLO-1)/NPROMA+1
      IST=1
      IEND=ICEND
  
      DO JROF=IST,IEND
  
        ! Find index for simple tropopause level definition...
        ! Latitude: YDGSGEOM(IBL)%GELAT(JROF)
        ZPTROP=24000_JPRB - 14800_JPRB * (COS(YDGSGEOM(IBL)%GELAT(JROF)))**4_JPIM
        DO JLEV = 1,NFLEVG
          ZPLEV = YDVAB%VAH(JLEV) + YDVAB%VBH(JLEV) * 101300.0_JPRB 
          IF (ZPLEV< ZPTROP) THEN
            JLEV_TROP = JLEV
          ELSE
            ! level found; exit this loop
            CYCLE
          ENDIF
        ENDDO
!Only stratospheric levels to be updated - remove two levels to be sure.
        DO JLEV=1,JLEV_TROP-2
          ZFAMA   =0._JPRB
          !new Family tracer
          DO JTRAC=1,INJ
            ZFAMA = ZFAMA + ZRMCOMP(JTRAC) * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
          ENDDO
          ZFAMRAT = PGFLT1(JROF,JLEV,YCHEM(IFAM)%MP1,IBL) / ZFAMA 
  
! make sure that individual BRX components after SL equals family after SL   
          DO JTRAC=1,INJ
            PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)   = ZFAMRAT * PGFLT1(JROF,JLEV,YCHEM(ICOMP(JTRAC))%MP1,IBL)
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  ENDIF ! KCHEM_NOXADV=4  
  

ENDIF !KMODE = 2


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CHEM_NOYADV',1,ZHOOK_HANDLE)
END SUBROUTINE CHEM_NOYADV 
