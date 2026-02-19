!!!#ifdef RS6K
!!!@PROCESS HOT NOSTRICT
!!!#endif
SUBROUTINE SRTM_VRTQDR&
 &(YDDIMV,KIDIA, KFDIA, KLEV , KW,&
 &PREF , PREFD, PTRA , PTRAD,&
 &PDBT , PRDND, PRUP , PRUPD , PTDBT,&
 &PFD  , PFU  , PRMU0&
 & )  

!**** *SRTM_VRTQDR* - VERTICAL QUADRATURE

!     PURPOSE.
!     --------

!          THIS ROUTINE PERFORMS THE VERTICAL INTEGRATION

!**   INTERFACE.
!     ----------

!          *SRTM_VRTQDR* IS CALLED FROM *SRTM_SPCVRT*

!        IMPLICIT ARGUMENTS :
!        --------------------

!     ==== INPUTS ===
!     ==== OUTPUTS ===

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
!        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

!     AUTHOR.
!     -------
!        from Howard Barker
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 02-10-04
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D.Salmond  31-Oct-2007 Vector version in the style of RRTM from Meteo France & NEC
!        J.Hague  : 21-10-2014 Vector Optimisation for Cray
!        F. Vana  : 05-Mar-2015 Support for single precision
!     ------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOESRTM  , ONLY : JPGPT 
!USE YOESWN   , ONLY : NDBUG

IMPLICIT NONE

!     ------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA, KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KW(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREF(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PREFD(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRA(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTRAD(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBT(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PRDND(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRUP(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PRUPD(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTDBT(KIDIA:KFDIA,YDDIMV%NFLEVG+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PFD(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT) 
REAL(KIND=JPRB)   ,INTENT(INOUT)   :: PFU(KIDIA:KFDIA,YDDIMV%NFLEVG+1,JPGPT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRMU0(KIDIA:KFDIA)
!     ------------------------------------------------------------------

!              ------------

REAL(KIND=JPRB) :: ZTDN(KIDIA:KFDIA,YDDIMV%NFLEVG+1)  
INTEGER(KIND=JPIM) :: INDEX(KIDIA:KFDIA),ICOUNT,IC

INTEGER(KIND=JPIM) :: IKP, IKX, JK, JL

REAL(KIND=JPRB) :: ZREFLECT, ZEPSILON
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

! PREF(JK)   direct reflectance
! PREFD(JK)  diffuse reflectance
! PTRA(JK)   direct transmittance
! PTRAD(JK)  diffuse transmittance

! PDBT(JK)   layer mean direct beam transmittance
! PTDBT(JK)  total direct beam transmittance at levels

IF (LHOOK) CALL DR_HOOK('SRTM_VRTQDR',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

ZEPSILON = 100._JPRB*EPSILON(ZEPSILON)

!-- link lowest layer with surface

IC=0
DO JL = KIDIA, KFDIA
  IF (PRMU0(JL) > 0.0_JPRB) THEN
    IC=IC+1
    INDEX(IC)=JL
  ENDIF
ENDDO
ICOUNT=IC

IF(ICOUNT==0)THEN
  IF (LHOOK) CALL DR_HOOK('SRTM_VRTQDR',1,ZHOOK_HANDLE)
  RETURN
ENDIF

! IF ICOUNT=KFDIA-KIDIA+1 remove indexed addressing to enable vectorisation

IF(ICOUNT == KFDIA-KIDIA+1) THEN

  !-- link lowest layer with surface

  DO JL = KIDIA, KFDIA
    ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PREFD(JL,KLEV+1)*PREFD(JL,KLEV),ZEPSILON)
    PRUP(JL,KLEV)=PREF(JL,KLEV)+(PTRAD(JL,KLEV)* &
     & ((PTRA(JL,KLEV)-PDBT(JL,KLEV))*PREFD(JL,KLEV+1)+ &
     & PDBT(JL,KLEV)*PREF(JL,KLEV+1)))*ZREFLECT
    PRUPD(JL,KLEV)=PREFD(JL,KLEV)+PTRAD(JL,KLEV)* &
     & PTRAD(JL,KLEV)*PREFD(JL,KLEV+1)*ZREFLECT
  ENDDO

  !-- pass from bottom to top 

  DO JK=1,KLEV-1
    IKP=KLEV+1-JK                       
    IKX=IKP-1
    DO JL = KIDIA, KFDIA
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PRUPD(JL,IKP)*PREFD(JL,IKX),ZEPSILON)
      PRUP(JL,IKX)=PREF(JL,IKX)+(PTRAD(JL,IKX)* &
       & ((PTRA(JL,IKX)-PDBT(JL,IKX))*PRUPD(JL,IKP)+ &
       & PDBT(JL,IKX)*PRUP(JL,IKP)))*ZREFLECT
      PRUPD(JL,IKX)=PREFD(JL,IKX)+PTRAD(JL,IKX)* &
       & PTRAD(JL,IKX)*PRUPD(JL,IKP)*ZREFLECT
    ENDDO
  ENDDO

  !-- upper boundary conditions

  DO JL = KIDIA, KFDIA
    ZTDN(JL,1)=1.0_JPRB
    PRDND(JL,1)=0.0_JPRB
    ZTDN(JL,2)=PTRA(JL,1)
    PRDND(JL,2)=PREFD(JL,1)
  ENDDO

  !-- pass from top to bottom

  DO JK=2,KLEV
    IKP=JK+1
    DO JL = KIDIA, KFDIA
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PREFD(JL,JK)*PRDND(JL,JK),ZEPSILON)
      ZTDN(JL,IKP)=PTDBT(JL,JK)*PTRA(JL,JK)+ &
       & (PTRAD(JL,JK)*((ZTDN(JL,JK)-PTDBT(JL,JK))+ &
       & PTDBT(JL,JK)*PREF(JL,JK)*PRDND(JL,JK))) * ZREFLECT
      PRDND(JL,IKP)=PREFD(JL,JK)+PTRAD(JL,JK)*PTRAD(JL,JK) &
       & *PRDND(JL,JK)*ZREFLECT
    ENDDO
  ENDDO

  !-- up and down-welling fluxes at levels
  
  DO JK=1,KLEV+1
    DO JL = KIDIA, KFDIA
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB - PRDND(JL,JK)*PRUPD(JL,JK),ZEPSILON)
      PFU(JL,JK,KW(JL))=(PTDBT(JL,JK)*PRUP(JL,JK) + &
       & (ZTDN(JL,JK)-PTDBT(JL,JK))*PRUPD(JL,JK))*ZREFLECT
      PFD(JL,JK,KW(JL))=PTDBT(JL,JK) + (ZTDN(JL,JK)-PTDBT(JL,JK)+ &
       & PTDBT(JL,JK)*PRUP(JL,JK)*PRDND(JL,JK))*ZREFLECT
    ENDDO
  ENDDO

!     ------------------------------------------------------------------

ELSE ! ICOUNT /= KFDIA-KIDIA+1 Non-Vector (Original)

  !-- link lowest layer with surface

  DO IC=1,ICOUNT
    JL=INDEX(IC)
    ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PREFD(JL,KLEV+1)*PREFD(JL,KLEV),ZEPSILON)
    PRUP(JL,KLEV)=PREF(JL,KLEV)+(PTRAD(JL,KLEV)* &
     & ((PTRA(JL,KLEV)-PDBT(JL,KLEV))*PREFD(JL,KLEV+1)+ &
     & PDBT(JL,KLEV)*PREF(JL,KLEV+1)))*ZREFLECT
    PRUPD(JL,KLEV)=PREFD(JL,KLEV)+PTRAD(JL,KLEV)* &
     & PTRAD(JL,KLEV)*PREFD(JL,KLEV+1)*ZREFLECT
  ENDDO
  
  !-- pass from bottom to top 
  
  DO JK=1,KLEV-1
    IKP=KLEV+1-JK                       
    IKX=IKP-1
    DO IC=1,ICOUNT
      JL=INDEX(IC)
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PRUPD(JL,IKP)*PREFD(JL,IKX),ZEPSILON)
      PRUP(JL,IKX)=PREF(JL,IKX)+(PTRAD(JL,IKX)* &
       & ((PTRA(JL,IKX)-PDBT(JL,IKX))*PRUPD(JL,IKP)+ &
       & PDBT(JL,IKX)*PRUP(JL,IKP)))*ZREFLECT
      PRUPD(JL,IKX)=PREFD(JL,IKX)+PTRAD(JL,IKX)* &
       & PTRAD(JL,IKX)*PRUPD(JL,IKP)*ZREFLECT
    ENDDO
  ENDDO

  !-- upper boundary conditions
  
  DO IC=1,ICOUNT
    JL=INDEX(IC)
    ZTDN(JL,1)=1.0_JPRB
    PRDND(JL,1)=0.0_JPRB
    ZTDN(JL,2)=PTRA(JL,1)
    PRDND(JL,2)=PREFD(JL,1)
  ENDDO
  
  !-- pass from top to bottom
  
  DO JK=2,KLEV
    IKP=JK+1
    DO IC=1,ICOUNT
      JL=INDEX(IC)
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB -PREFD(JL,JK)*PRDND(JL,JK),ZEPSILON)
      ZTDN(JL,IKP)=PTDBT(JL,JK)*PTRA(JL,JK)+ &
       & (PTRAD(JL,JK)*((ZTDN(JL,JK)-PTDBT(JL,JK))+ &
       & PTDBT(JL,JK)*PREF(JL,JK)*PRDND(JL,JK))) * ZREFLECT
      PRDND(JL,IKP)=PREFD(JL,JK)+PTRAD(JL,JK)*PTRAD(JL,JK) &
       & *PRDND(JL,JK)*ZREFLECT
    ENDDO
  ENDDO
  
  !-- up and down-welling fluxes at levels
  
  DO JK=1,KLEV+1
    DO IC=1,ICOUNT
      JL=INDEX(IC)
      ZREFLECT=1.0_JPRB / MAX(1.0_JPRB - PRDND(JL,JK)*PRUPD(JL,JK),ZEPSILON)
      PFU(JL,JK,KW(JL))=(PTDBT(JL,JK)*PRUP(JL,JK) + &
       & (ZTDN(JL,JK)-PTDBT(JL,JK))*PRUPD(JL,JK))*ZREFLECT
      PFD(JL,JK,KW(JL))=PTDBT(JL,JK) + (ZTDN(JL,JK)-PTDBT(JL,JK)+ &
       & PTDBT(JL,JK)*PRUP(JL,JK)*PRDND(JL,JK))*ZREFLECT
    ENDDO
  ENDDO
 
ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SRTM_VRTQDR',1,ZHOOK_HANDLE)
END SUBROUTINE SRTM_VRTQDR
