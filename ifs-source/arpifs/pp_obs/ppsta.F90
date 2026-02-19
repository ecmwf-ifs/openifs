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

SUBROUTINE PPSTA(CDATM,KPROMA,KSTART,KPROF,KLEVP,KLOLEV,&
 & PRPRES,PLNPRES,PSTT,PSTFI)  

!**** *PPSTA*  - Integrate the STANDARD ATMOSPHERE

!     Purpose.
!     --------
!           Integrate STANDARD ATMOSPHERE for geopotental.

!**   Interface.
!     ----------
!        *CALL* *PPSTA(KPROMA,KSTART,KPROF,KLEVP,KLOLEV
!                     ,PRPRES,PLNPRES,PSTT,PSTFI)

!        Explicit arguments :
!        --------------------
!        YDSTA                     - DUMMY TSTA ARGUMENT.              (INPUT-C)
!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT-C)
!        KSTART                    - START OF WORK.                    (INPUT-C)
!        KPROF                     - DEPTH OF WORK.                    (INPUT-C)
!        KLEVP                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT-C)
!        KLOLEV                    - BEGINING FOR THE INTERPOLATION    (INPUT-C)
!        PRPRES(KPROMA,KLEVP)      - LIST OF PRESSURES                 (INPUT)
!        PLNPRES(KPROMA,KLEVP)     - LOG(PRPRES)                      (INPUT)
!        PSTT(KPROMA,KLEVP)        - ICAO TEMPERATURE                  (OUTPUT)
!        PSTFI(KPROMA,KLEVP)       - ICAO GEOPOTENTIAL                 (OUTPUT)

!        Implicit arguments :
!        --------------------
!        Common YOMSTA

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Erik Andersson ECMWF

!     Modifications.
!     --------------
!        Original : 82-04-10
!        Modified : ??-??-?? J.Latour ?? tidied up by M.Hamrud
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RG, RD
USE YOMSTA   , ONLY : RDTDZ1, RDTDZ2, RDTDZ3, RDTDZ4, &
 & RDTDZ5, RDTDZ6, RDTDZ7, RDTDZ8, RDTDZ9, RPABOV, RPMEPO, &
 & RPMES2, RPMESO, RPSTPO, RPSTR2, RPSTRA, RPTROP, RTABOV, &
 & RTMEPO, RTMES2, RTMESO, RTSTPO, RTSTR2, RTSTRA, RTSUR,  &
 & RTTROP, RZABOV, RZMEPO, RZMES2, RZMESO, RZSTPO, RZSTR2, &
 & RZSTRA, RZTROP, VDTDZ1, VDTDZ2, VDTDZ3, VDTDZ4, VDTDZ5, &
 & VDTDZ6, VDTDZ7, VDTDZ8, VDTDZ9, VPABOV, VPMEPO, VPMES2, &
 & VPMESO, VPSTPO, VPSTR2, VPSTRA, VPTROP, VTABOV, VTMEPO, &
 & VTMES2, VTMESO, VTSTPO, VTSTR2, VTSTRA, VTSUR , VTTROP, &
 & VZABOV, VZMEPO, VZMES2, VZMESO, VZSTPO, VZSTR2, VZSTRA, VZTROP
USE YOMVERT  , ONLY : VP00

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
CHARACTER(LEN=*)  ,INTENT(IN)    :: CDATM
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLOLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRES(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLNPRES(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTT(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTFI(KPROMA,KLEVP) 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IPRO(KPROMA,KLEVP)

INTEGER(KIND=JPIM) :: JPPRO

PARAMETER (JPPRO=9)
REAL(KIND=JPRB) :: ZZ(JPPRO),ZDTDZ(JPPRO),ZT(JPPRO),ZP(JPPRO)
REAL(KIND=JPRB) :: ZIP(JPPRO),ZLNP(JPPRO),ZIDTDZ(JPPRO)
LOGICAL :: LLGRZER

INTEGER(KIND=JPIM) :: IPR, JJPR, JL, JLEV

LOGICAL :: LLPREF

REAL(KIND=JPRB) :: ZROG, ZSTZ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PPSTA',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.    Initialize
!              ----------
!      ZZ : Altitude du changement de gradient vertical de temperarure
!      ZP : Pression a cette altitude
!      ZT : Temperature a cette meme altitude
!      ZDTDZ : Gradient vertical de temperature au-dessus de cette altitude

!*    1.0 INITIALISATION DE ZZ ET DE ZDTDZ

LLPREF=CDATM == 'PPREF'
IF(.NOT.LLPREF) THEN
  ZZ(1)=0.0_JPRB
  ZT(1)=RTSUR
  ZP(1)=VP00
  ZDTDZ(1)=RDTDZ1
  ZZ(2)=RZTROP
  ZT(2)=RTTROP
  ZP(2)=RPTROP
  ZDTDZ(2)=RDTDZ2
  ZZ(3)=RZSTRA
  ZT(3)=RTSTRA
  ZP(3)=RPSTRA
  ZDTDZ(3)=RDTDZ3
  ZZ(4)=RZSTR2
  ZT(4)=RTSTR2
  ZP(4)=RPSTR2
  ZDTDZ(4)=RDTDZ4
  ZZ(5)=RZSTPO
  ZT(5)=RTSTPO
  ZP(5)=RPSTPO
  ZDTDZ(5)=RDTDZ5
  ZZ(6)=RZMESO
  ZT(6)=RTMESO
  ZP(6)=RPMESO
  ZDTDZ(6)=RDTDZ6
  ZZ(7)=RZMES2
  ZT(7)=RTMES2
  ZP(7)=RPMES2
  ZDTDZ(7)=RDTDZ7
  ZZ(8)=RZMEPO
  ZT(8)=RTMEPO
  ZP(8)=RPMEPO
  ZDTDZ(8)=RDTDZ8
  ZZ(9)=RZABOV
  ZT(9)=RTABOV
  ZP(9)=RPABOV
  ZDTDZ(9)=RDTDZ9
ELSE
  ZZ(1)=0.0_JPRB
  ZT(1)=VTSUR
  ZP(1)=VP00
  ZDTDZ(1)=VDTDZ1
  ZZ(2)=VZTROP
  ZT(2)=VTTROP
  ZP(2)=VPTROP
  ZDTDZ(2)=VDTDZ2
  ZZ(3)=VZSTRA
  ZT(3)=VTSTRA
  ZP(3)=VPSTRA
  ZDTDZ(3)=VDTDZ3
  ZZ(4)=VZSTR2
  ZT(4)=VTSTR2
  ZP(4)=VPSTR2
  ZDTDZ(4)=VDTDZ4
  ZZ(5)=VZSTPO
  ZT(5)=VTSTPO
  ZP(5)=VPSTPO
  ZDTDZ(5)=VDTDZ5
  ZZ(6)=VZMESO
  ZT(6)=VTMESO
  ZP(6)=VPMESO
  ZDTDZ(6)=VDTDZ6
  ZZ(7)=VZMES2
  ZT(7)=VTMES2
  ZP(7)=VPMES2
  ZDTDZ(7)=VDTDZ7
  ZZ(8)=VZMEPO
  ZT(8)=VTMEPO
  ZP(8)=VPMEPO
  ZDTDZ(8)=VDTDZ8
  ZZ(9)=VZABOV
  ZT(9)=VTABOV
  ZP(9)=VPABOV
  ZDTDZ(9)=VDTDZ9
ENDIF

!  LLGRZER : .true. if one or more values of ZDTDZ are zero, .false. otherwise
LLGRZER=.FALSE.
!CP Re-written because may crash with some compiling options
IF (ANY(ZDTDZ(1:JPPRO) == 0.0_JPRB)) LLGRZER=.TRUE.

ZIDTDZ(1:JPPRO)=0.0_JPRB
DO JJPR=1,JPPRO
  ZLNP(JJPR)=LOG(ZP(JJPR))
  ZIP(JJPR)=1.0_JPRB/ZP(JJPR)
  IF (ZDTDZ(JJPR) /= 0.0_JPRB) ZIDTDZ(JJPR)=1.0_JPRB/ZDTDZ(JJPR)
ENDDO

!*    1.2 Calcul de l'Altitude (ZSTZ), la Temperature (PSTT),
!         le Geopotentiel (PSTFI).

DO JLEV=KLOLEV,KLEVP
  DO JL=KSTART,KPROF
    IF(PRPRES(JL,JLEV) >= ZP(2)) THEN
      IPRO(JL,JLEV)=1
    ELSEIF(PRPRES(JL,JLEV) >= ZP(3)) THEN
      IPRO(JL,JLEV)=2
    ELSEIF(PRPRES(JL,JLEV) >= ZP(4)) THEN
      IPRO(JL,JLEV)=3
    ELSEIF(PRPRES(JL,JLEV) >= ZP(5)) THEN
      IPRO(JL,JLEV)=4
    ELSEIF(PRPRES(JL,JLEV) >= ZP(6)) THEN
      IPRO(JL,JLEV)=5
    ELSEIF(PRPRES(JL,JLEV) >= ZP(7)) THEN
      IPRO(JL,JLEV)=6
    ELSEIF(PRPRES(JL,JLEV) >= ZP(8)) THEN
      IPRO(JL,JLEV)=7
    ELSEIF(PRPRES(JL,JLEV) >= ZP(9)) THEN
      IPRO(JL,JLEV)=8
    ELSE
      IPRO(JL,JLEV)=9
    ENDIF
  ENDDO
ENDDO

ZROG=RD/RG

!     if LLGRZER=.false. the loop does NOT need the test on ZDTDZ

IF(.NOT.LLGRZER) THEN

!ocl fusion
  DO JLEV=KLOLEV,KLEVP
    DO JL=KSTART,KPROF
      IPR=IPRO(JL,JLEV)
      ZSTZ=ZZ(IPR)+ZT(IPR)*ZIDTDZ(IPR)&
           & *(REAL(PRPRES(JL,JLEV)*ZIP(IPR),JPRD)**REAL(-ZDTDZ(IPR)*ZROG, JPRD)-1.0_JPRD)
      PSTT(JL,JLEV)=ZT(IPR)+ZDTDZ(IPR)*(ZSTZ-ZZ(IPR))
      PSTFI(JL,JLEV)=ZSTZ*RG
    ENDDO
  ENDDO
!ocl endfusion

ELSE

  DO JLEV=KLOLEV,KLEVP
    DO JL=KSTART,KPROF
      IPR=IPRO(JL,JLEV)
      IF(ZDTDZ(IPR) /= 0.0_JPRB) THEN
        ZSTZ=ZZ(IPR)+ZT(IPR)*ZIDTDZ(IPR)&
             & *(REAL(PRPRES(JL,JLEV)*ZIP(IPR),JPRD)**REAL(-ZDTDZ(IPR)*ZROG,JPRD)-1.0_JPRD)
        PSTT(JL,JLEV)=ZT(IPR)+ZDTDZ(IPR)*(ZSTZ-ZZ(IPR))
      ELSE
        ZSTZ=ZZ(IPR)-ZT(IPR)*ZROG*(PLNPRES(JL,JLEV)-ZLNP(IPR))
        PSTT(JL,JLEV)=ZT(IPR)
      ENDIF
      PSTFI(JL,JLEV)=ZSTZ*RG
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPSTA',1,ZHOOK_HANDLE)
END SUBROUTINE PPSTA
