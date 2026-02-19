! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS HOT  
#endif
SUBROUTINE LWVDR &
 & ( YDELWRAD,KIDIA,  KFDIA, KLON , KLEV  , KTRAER,&
 & PABCU,  PDBDT,&
 & PGA  ,  PGB,&
 & PCNTRB, PDISD, PDISU, PDWFSU, &
! for adjoint computation
 & PCOND  , PTTPA  , PTTPB , PTTA , PTTB  ,&
 & PZZA   , PZZB   , PXNA  , PXNB , PXDIVA,&
 & PXDIVB , PZZA1  , PZZB1 , PXNA1, PXNB1 ,&
 & PXDIVA1, PXDIVB1, PTTPA1, PTTPB1 &
 & )  

!**** *LWVDR*   - L.W., VERTICAL INTEGRATION, DISTANT LAYERS

!     PURPOSE.
!     --------
!           CARRIES OUT THE VERTICAL INTEGRATION ON THE DISTANT LAYERS

!**   INTERFACE.
!     ----------

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! PABCU  : (KLON,NUA,3*KLEV+1) ; ABSORBER AMOUNTS
! PDBDT  : (KLON,KLEV)         ; LAYER PLANCK FUNCTION GRADIENT
! PGA, PGB                     ; PADE APPROXIMANTS
!     ==== OUTPUTS ===
! PCNTRB : (KLON,KLEV+1,KLEV+1); ENERGY EXCHANGE MATRIX
! PDIS.. : (KLON,KLEV+1)       ; CONTRIBUTION BY DISTANT LAYERS
! PDWFSU : (KLON,NSIL)         ; SPECTRAL DOWNWARD FLUX AT SURFACE

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!          1. PERFORMS THE VERTICAL INTEGRATION CORRESPONDING TO THE
!     CONTRIBUTIONS OF THE DISTANT LAYERS USING TRAPEZOIDAL RULE

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 89-07-14
!------
!        M. Janiskova  16-Jan-2006 reduced routine based on lwvd 
!xxxxx
!        M.Janiskova   23-Nov-2006 reduced version for H2O and CO2
!                                  based on lwvd
!        M.Janiskova   23-Jan-2008 condition PCOND for reduced loops
!        J.Hague       08-Feb-2008 code optimization
!        M.Janiskova   20-May-2008 using shorter NPROMALW instead  
!                                  of KLON
!        J.Hague : 21-10-2014 Vector Optimisation for Cray
!-----------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOELW    , ONLY : NSIL     ,NIPD     ,NUA      ,NG1P1
USE YOELWRAD , ONLY : TELWRAD

IMPLICIT NONE

TYPE(TELWRAD)     ,INTENT(IN)    :: YDELWRAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAER 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PABCU(KLON,NUA,3*KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDBDT(KLON,NSIL,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGA(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PGB(KLON,NIPD,2,KLEV) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCNTRB(KLON,KLEV+1,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISD(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PDISU(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PDWFSU(KLON,NSIL) 
! for adjoint computation
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCOND(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPA(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPB(YDELWRAD%NPROMALW,KLEV+1,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVA(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVB(YDELWRAD%NPROMALW,KTRAER,KLEV,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PZZB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXNB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVA1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PXDIVB1(YDELWRAD%NPROMALW,KTRAER,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPA1(KLON,KLEV+1)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTTPB1(KLON,KLEV+1)

!-----------------------------------------------------------------------

!*       0.1   ARGUMENTS
!              ---------

!-----------------------------------------------------------------------

!              ------------

INTEGER(KIND=JPIM) :: IJKL, IKD1, IKD2, IKJ, IKJP1, IKM1, IKN,&
 & IKP1, IKU1, IKU2, JA, JK, JKJ, JL, JLK, IJA, IKJ1
LOGICAL :: LL_DO(KLON)
INTEGER(KIND=JPIM) :: IC, ILOOP, IJL

INTEGER(KIND=JPIM) :: JJ, JN, JX(KLON)

REAL(KIND=JPRB) :: ZTTP(KIDIA:KFDIA,KTRAER,KLEV+1)
REAL(KIND=JPRB) :: ZWW1, ZWW2, ZWW3, ZWW4, ZWW5, ZWW6
REAL(KIND=JPRB) :: ZWD1, ZWD2, ZWD3, ZWD4, ZWD5, ZWD6
REAL(KIND=JPRB) :: ZXD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!#include "lwttm.intfb.h"

!-----------------------------------------------------------------------

!*         1.    INITIALIZATION
!                --------------

!*         1.1     INITIALIZE LAYER CONTRIBUTIONS
!                  ------------------------------

IF (LHOOK) CALL DR_HOOK('LWVDR',0,ZHOOK_HANDLE)

!ILOOP = (KFDIA/(KFDIA-KIDIA+1)-1)*NPROMALW
ILOOP = KIDIA-1

DO JK = 1, KLEV+1
  DO JL = KIDIA,KFDIA
    PDISD(JL,JK) = 0.0_JPRB
    PDISU(JL,JK) = 0.0_JPRB
  ENDDO
ENDDO

!*         1.2     INITIALIZE TRANSMISSION FUNCTIONS
!                  ---------------------------------

DO JK = 1, KLEV+1
  DO JA = 1, KTRAER
    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

      ZTTP(JL,JA,JK) = 1.0_JPRB
!for adjoint computation
      PTTA(IJL,JA,1,JK) = 1.0_JPRB
      PTTB(IJL,JA,1,JK) = 1.0_JPRB
    ENDDO
  ENDDO
ENDDO

!     ------------------------------------------------------------------

!*         2.      VERTICAL INTEGRATION
!                  --------------------

!*         2.2     CONTRIBUTION FROM DISTANT LAYERS
!                  ---------------------------------

!*         2.2.1   DISTANT AND ABOVE LAYERS
!                  ------------------------

!*         2.2.2   FIRST UPPER LEVEL
!                  -----------------

DO JK = KLEV-1,1,-1
  IKP1=JK+1
  IKN=(JK-1)*NG1P1+1
  IKD1= JK  *NG1P1+1

  DO JA = 1 , 8
    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

      PZZA1(IJL,JA,JK) = SQRT(PABCU(JL,JA,IKN)-PABCU(JL,JA,IKD1)+1.D-25)
      ZXD = PGB( JL,JA,1,JK)&
       & + PZZA1(IJL,JA,JK)* (PGB( JL,JA,2,JK) + PZZA1(IJL,JA,JK) )
      PXNA1(IJL,JA,JK) = PGA( JL,JA,1,JK)&
       & + PZZA1(IJL,JA,JK)*PGA( JL,JA,2,JK)
      PXDIVA1(IJL,JA,JK) = 1.0_JPRB/ZXD
      ZTTP(JL,JA,1) = PXNA1(IJL,JA,JK)*PXDIVA1(IJL,JA,JK)
    ENDDO
  ENDDO

  DO JL = KIDIA,KFDIA
    PTTPA1(JL,JK) = ZTTP(JL,3,1)  ! for adjoint computation
    IF (PTTPA1(JL,JK) < 0.0_JPRB) THEN
      ZTTP(JL,3,1) = 0.0_JPRB  ! for adjoint computation
    ENDIF
    ZTTP(JL,9,1) = ZTTP(JL,8,1) ! for adjoint computation
  ENDDO

!*         2.2.3   HIGHER UP
!                  ---------

  DO JKJ=IKP1,KLEV
    IKJP1=JKJ+1
    IKD2= JKJ  *NG1P1+1
    IKJ1=JKJ-JK

    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP
!      IF (JK > 1 .AND. JK < KLEV-2 .AND. (IKJP1-JK) > 3) THEN
      IF (JK > 1 .AND. (IKJP1-JK) > 3) THEN
        PCOND(IJL,JKJ,JK) = ABS(PCNTRB(JL,IKJP1,JK+2)-PCNTRB(JL,IKJP1,JK+1))
      ELSE
        PCOND(IJL,JKJ,JK) = 9999._JPRB
      ENDIF
    ENDDO

    IC=0
    LL_DO(:)=.TRUE.
    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

!exp_A   
!      IF (PCOND(IJL,JKJ,JK) < 0.0001_JPRB) THEN
!exp_B      
!      IF (PCOND(IJL,JKJ,JK) < 0.00001_JPRB) THEN
!exp_C
      IF (PCOND(IJL,JKJ,JK) < 0.001_JPRB) THEN
        LL_DO(JL)=.FALSE.
        IC=IC+1
        PCNTRB(JL,IKJP1,JK) = PCNTRB(JL,IKJP1,JK+1)
        PDISD(JL,JK)=PDISD(JL,JK)+PCNTRB(JL,IKJP1,JK)
      ENDIF
    ENDDO

    JJ=0
    DO JL = KIDIA,KFDIA
      IF (LL_DO(JL)) THEN
        JJ=JJ+1
        JX(JJ)=JL
      ENDIF
    ENDDO
    JN=JJ

! Vector version without indirect addressing if IC=0
    IF (IC == 0) THEN

      DO JA = 1 , 8
        IJA = (JA-1)*KLEV+JKJ
        DO JL = KIDIA,KFDIA
            IJL = JL - ILOOP
            PZZA(IJL,JA,JKJ,JK) = SQRT(PABCU(JL,JA,IKN)-PABCU(JL,JA,IKD2)+1.D-25)
            ZXD = PGB( JL,JA,1,JKJ)&
             &+ PZZA(IJL,JA,JKJ,JK)* (PGB( JL,JA,2,JKJ) + PZZA(IJL,JA,JKJ,JK) )
            PXNA(IJL,JA,JKJ,JK) = PGA( JL,JA,1,JKJ)&
             & + PZZA(IJL,JA,JKJ,JK)*PGA( JL,JA,2,JKJ)
            PXDIVA(IJL,JA,JKJ,JK) = 1.0_JPRB/ZXD
            ZTTP(JL,JA,IKJ1+1) = PXNA(IJL,JA,JKJ,JK)*PXDIVA(IJL,JA,JKJ,JK)
        ENDDO
      ENDDO

      DO JL = KIDIA,KFDIA
          IJL = JL - ILOOP
          PTTPA(IJL,IKJ1+1,JK) = ZTTP(JL,3,IKJ1+1) ! for adjoint computation
          IF (PTTPA(IJL,IKJ1+1,JK) < 0.0_JPRB) THEN
            ZTTP(JL,3,IKJ1+1) = 0.0_JPRB ! for adjoint computation
          ENDIF
          ZTTP(JL,9,IKJ1+1) = ZTTP(JL,8,IKJ1+1) ! for adjoint computation
      ENDDO

      DO JA = 1, KTRAER
         IJA = (JA-1)*KLEV+JKJ
         DO JL = KIDIA,KFDIA
             IJL = JL - ILOOP
             PTTA(IJL,JA,JKJ,JK) = (ZTTP(JL,JA,IKJ1)+ZTTP(JL,JA,IKJ1+1))*0.5_JPRB
         ENDDO
       ENDDO
        
      DO JL = KIDIA,KFDIA
          IJL = JL - ILOOP

          ZWW1=PDBDT(JL,1,JKJ)*PTTA(IJL,1,JKJ,JK)
          ZWW2=PDBDT(JL,2,JKJ)*PTTA(IJL,2,JKJ,JK)*PTTA(IJL,7,JKJ,JK)
          ZWW3=PDBDT(JL,3,JKJ)*PTTA(IJL,4,JKJ,JK)*PTTA(IJL,8,JKJ,JK)
          ZWW4=PDBDT(JL,4,JKJ)*PTTA(IJL,5,JKJ,JK)*PTTA(IJL,9,JKJ,JK)
          ZWW5=PDBDT(JL,5,JKJ)*PTTA(IJL,3,JKJ,JK)
          ZWW6=PDBDT(JL,6,JKJ)*PTTA(IJL,6,JKJ,JK)
          PCNTRB(JL,IKJP1,JK)=ZWW1+ZWW2+ZWW3+ZWW4+ZWW5+ZWW6
          PDISD(JL,JK)=PDISD(JL,JK)+PCNTRB(JL,IKJP1,JK)
          IF (JK == 1) THEN
            PDWFSU(JL,1)=PDWFSU(JL,1)+ZWW1
            PDWFSU(JL,2)=PDWFSU(JL,2)+ZWW2
            PDWFSU(JL,3)=PDWFSU(JL,3)+ZWW3
            PDWFSU(JL,4)=PDWFSU(JL,4)+ZWW4
            PDWFSU(JL,5)=PDWFSU(JL,5)+ZWW5
            PDWFSU(JL,6)=PDWFSU(JL,6)+ZWW6
          ENDIF
      ENDDO

    ELSEIF(IC < KFDIA-KIDIA+1) THEN

      DO JA = 1 , 8
        IJA = (JA-1)*KLEV+JKJ
        DO JL = KIDIA,KFDIA
          IF (LL_DO(JL)) THEN
            IJL = JL - ILOOP
            PZZA(IJL,JA,JKJ,JK) = SQRT(PABCU(JL,JA,IKN)-PABCU(JL,JA,IKD2)+1.D-25)
            ZXD = PGB( JL,JA,1,JKJ)&
             &+ PZZA(IJL,JA,JKJ,JK)* (PGB( JL,JA,2,JKJ) + PZZA(IJL,JA,JKJ,JK) )
            PXNA(IJL,JA,JKJ,JK) = PGA( JL,JA,1,JKJ)&
             & + PZZA(IJL,JA,JKJ,JK)*PGA( JL,JA,2,JKJ)
            PXDIVA(IJL,JA,JKJ,JK) = 1.0_JPRB/ZXD
            ZTTP(JL,JA,IKJ1+1) = PXNA(IJL,JA,JKJ,JK)*PXDIVA(IJL,JA,JKJ,JK)
          ENDIF
        ENDDO
      ENDDO

      DO JL = KIDIA,KFDIA
        IF (LL_DO(JL)) THEN
          IJL = JL - ILOOP
          PTTPA(IJL,IKJ1+1,JK) = ZTTP(JL,3,IKJ1+1) ! for adjoint computation
          IF (PTTPA(IJL,IKJ1+1,JK) < 0.0_JPRB) THEN
            ZTTP(JL,3,IKJ1+1) = 0.0_JPRB ! for adjoint computation
          ENDIF
          ZTTP(JL,9,IKJ1+1) = ZTTP(JL,8,IKJ1+1) ! for adjoint computation
        ENDIF
      ENDDO

      DO JA = 1, KTRAER
         IJA = (JA-1)*KLEV+JKJ
         DO JL = KIDIA,KFDIA
           IF (LL_DO(JL)) THEN
             IJL = JL - ILOOP
             PTTA(IJL,JA,JKJ,JK) = (ZTTP(JL,JA,IKJ1)+ZTTP(JL,JA,IKJ1+1))*0.5_JPRB
           ENDIF
         ENDDO
       ENDDO
        
      DO JL = KIDIA,KFDIA
        IF (LL_DO(JL)) THEN
          IJL = JL - ILOOP

          ZWW1=PDBDT(JL,1,JKJ)*PTTA(IJL,1,JKJ,JK)
          ZWW2=PDBDT(JL,2,JKJ)*PTTA(IJL,2,JKJ,JK)*PTTA(IJL,7,JKJ,JK)
          ZWW3=PDBDT(JL,3,JKJ)*PTTA(IJL,4,JKJ,JK)*PTTA(IJL,8,JKJ,JK)
          ZWW4=PDBDT(JL,4,JKJ)*PTTA(IJL,5,JKJ,JK)*PTTA(IJL,9,JKJ,JK)
          ZWW5=PDBDT(JL,5,JKJ)*PTTA(IJL,3,JKJ,JK)
          ZWW6=PDBDT(JL,6,JKJ)*PTTA(IJL,6,JKJ,JK)
          PCNTRB(JL,IKJP1,JK)=ZWW1+ZWW2+ZWW3+ZWW4+ZWW5+ZWW6
          PDISD(JL,JK)=PDISD(JL,JK)+PCNTRB(JL,IKJP1,JK)
          IF (JK == 1) THEN
            PDWFSU(JL,1)=PDWFSU(JL,1)+ZWW1
            PDWFSU(JL,2)=PDWFSU(JL,2)+ZWW2
            PDWFSU(JL,3)=PDWFSU(JL,3)+ZWW3
            PDWFSU(JL,4)=PDWFSU(JL,4)+ZWW4
            PDWFSU(JL,5)=PDWFSU(JL,5)+ZWW5
            PDWFSU(JL,6)=PDWFSU(JL,6)+ZWW6
          ENDIF
        ENDIF
      ENDDO

    ENDIF

  ENDDO
ENDDO

!*         2.2.4   DISTANT AND BELOW LAYERS
!                  ------------------------

!*         2.2.5   FIRST LOWER LEVEL
!                  -----------------

DO JK=3,KLEV+1
  IKN=(JK-1)*NG1P1+1
  IKM1=JK-1
  IKJ=JK-2
  IKU1= IKJ  *NG1P1+1

  DO JA = 1 , 8
    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

      PZZB1(IJL,JA,JK) = SQRT(PABCU(JL,JA,IKU1)-PABCU(JL,JA,IKN)+1.D-25)
      ZXD = PGB( JL,JA,1,IKJ)&
       & + PZZB1(IJL,JA,JK)* (PGB( JL,JA,2,IKJ) + PZZB1(IJL,JA,JK) )
      PXNB1(IJL,JA,JK) = PGA( JL,JA,1,IKJ)&
       & + PZZB1(IJL,JA,JK)*PGA( JL,JA,2,IKJ)
      PXDIVB1(IJL,JA,JK) = 1.0_JPRB/ZXD
      ZTTP(JL,JA,1) = PXNB1(IJL,JA,JK)*PXDIVB1(IJL,JA,JK)
    ENDDO
  ENDDO

  DO JL = KIDIA,KFDIA
    PTTPB1(JL,JK) = ZTTP(JL,3,1) ! for adjoint computation
    IF (PTTPB1(JL,JK) < 0.0_JPRB) THEN
      ZTTP(JL,3,1) = 0.0_JPRB ! for adjoint computation
    ENDIF
    ZTTP(JL,9,1) = ZTTP(JL,8,1) ! for adjoint computation
  ENDDO

!*         2.2.6   DOWN BELOW
!                  ----------

  DO JLK=1,IKJ
    IJKL=IKM1-JLK
    IKU2=(IJKL-1)*NG1P1+1

    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

!      IF (JK > 4 .AND. (JK-IJKL) > 3) THEN
      IF ((JK-IJKL) > 3) THEN
        PCOND(IJL,JLK,JK) = ABS(PCNTRB(JL,IJKL,JK-1)-PCNTRB(JL,IJKL,JK-2))
      ELSE
        PCOND(IJL,JLK,JK) = 9999._JPRB
      ENDIF
    ENDDO

    IC=0
    LL_DO(:)=.TRUE.
    DO JL = KIDIA,KFDIA
      IJL = JL - ILOOP

!exp_A     
!      IF (PCOND(IJL,JLK,JK) < 0.0001_JPRB) THEN
!exp_B
!      IF (PCOND(IJL,JLK,JK) < 0.00001_JPRB) THEN
!exp_C
      IF (PCOND(IJL,JLK,JK) < 0.001_JPRB) THEN
        IC=IC+1
        LL_DO(JL)=.FALSE.
        PCNTRB(JL,IJKL,JK) = PCNTRB(JL,IJKL,JK-1)
        PDISU(JL,JK)=PDISU(JL,JK)+PCNTRB(JL,IJKL,JK)
      ENDIF
    ENDDO

    JJ=0
    DO JL = KIDIA,KFDIA
      IF (LL_DO(JL)) THEN
        JJ=JJ+1
        JX(JJ)=JL
      ENDIF
    ENDDO
    JN=JJ

! Vector version without indirect addressing if IC=0
    IF( IC == 0 ) THEN
    
      DO JA = 1 , 8
        IJA = (JA-1)*KLEV+JLK
        DO JL = KIDIA,KFDIA
            IJL = JL - ILOOP
            PZZB(IJL,JA,JLK,JK) = SQRT(PABCU(JL,JA,IKU2)-PABCU(JL,JA,IKN)+1.D-25)
            ZXD = PGB( JL,JA,1,IJKL)&
             & + PZZB(IJL,JA,JLK,JK)*(PGB( JL,JA,2,IJKL) + PZZB(IJL,JA,JLK,JK) )
            PXNB(IJL,JA,JLK,JK) = PGA( JL,JA,1,IJKL)&
             & + PZZB(IJL,JA,JLK,JK)*PGA( JL,JA,2,IJKL)
            PXDIVB(IJL,JA,JLK,JK) = 1.0_JPRB/ZXD
            ZTTP(JL,JA,JLK+1) = PXNB(IJL,JA,JLK,JK)*PXDIVB(IJL,JA,JLK,JK)
        ENDDO
      ENDDO

      DO JL = KIDIA,KFDIA
          IJL = JL - ILOOP
          PTTPB(IJL,JLK+1,JK) = ZTTP(JL,3,JLK+1) ! for adjoint computation
          IF (PTTPB(IJL,JLK+1,JK) < 0.0_JPRB) THEN
            ZTTP(JL,3,JLK+1) = 0.0_JPRB ! for adjoint computation
          ENDIF
          ZTTP(JL,9,JLK+1) = ZTTP(JL,8,JLK+1) ! for adjoint computation
      ENDDO

       DO JA = 1, KTRAER
         DO JL = KIDIA,KFDIA
             IJL = JL - ILOOP
             IJA = (JA-1)*KLEV+JLK
             PTTB(IJL,JA,JLK,JK) = (ZTTP(JL,JA,JLK)+ZTTP(JL,JA,JLK+1))*0.5_JPRB
         ENDDO
       ENDDO

      DO JL = KIDIA,KFDIA
          IJL = JL - ILOOP
          ZWD1 = PDBDT(JL,1,IJKL)*PTTB(IJL,1,JLK,JK)
          ZWD2 = PDBDT(JL,2,IJKL)*PTTB(IJL,2,JLK,JK)*PTTB(IJL,7,JLK,JK)
          ZWD3 = PDBDT(JL,3,IJKL)*PTTB(IJL,4,JLK,JK)*PTTB(IJL,8,JLK,JK)
          ZWD4 = PDBDT(JL,4,IJKL)*PTTB(IJL,5,JLK,JK)*PTTB(IJL,9,JLK,JK)
          ZWD5 = PDBDT(JL,5,IJKL)*PTTB(IJL,3,JLK,JK)
          ZWD6 = PDBDT(JL,6,IJKL)*PTTB(IJL,6,JLK,JK)
          PCNTRB(JL,IJKL,JK) = ZWD1+ZWD2+ZWD3+ZWD4+ZWD5+ZWD6
          PDISU(JL,JK)=PDISU(JL,JK)+PCNTRB(JL,IJKL,JK)
      ENDDO

    ELSEIF(IC < KFDIA-KIDIA+1) THEN

      DO JA = 1 , 8
        IJA = (JA-1)*KLEV+JLK
        DO JL = KIDIA,KFDIA
          IF (LL_DO(JL)) THEN
            IJL = JL - ILOOP

            PZZB(IJL,JA,JLK,JK) = SQRT(PABCU(JL,JA,IKU2)-PABCU(JL,JA,IKN)+1.D-25)
            ZXD = PGB( JL,JA,1,IJKL)&
             & + PZZB(IJL,JA,JLK,JK)*(PGB( JL,JA,2,IJKL) + PZZB(IJL,JA,JLK,JK) )
            PXNB(IJL,JA,JLK,JK) = PGA( JL,JA,1,IJKL)&
             & + PZZB(IJL,JA,JLK,JK)*PGA( JL,JA,2,IJKL)
            PXDIVB(IJL,JA,JLK,JK) = 1.0_JPRB/ZXD
            ZTTP(JL,JA,JLK+1) = PXNB(IJL,JA,JLK,JK)*PXDIVB(IJL,JA,JLK,JK)
          ENDIF
        ENDDO
      ENDDO

      DO JL = KIDIA,KFDIA
        IF (LL_DO(JL)) THEN
          IJL = JL - ILOOP

          PTTPB(IJL,JLK+1,JK) = ZTTP(JL,3,JLK+1) ! for adjoint computation
          IF (PTTPB(IJL,JLK+1,JK) < 0.0_JPRB) THEN
            ZTTP(JL,3,JLK+1) = 0.0_JPRB ! for adjoint computation
          ENDIF
          ZTTP(JL,9,JLK+1) = ZTTP(JL,8,JLK+1) ! for adjoint computation
        ENDIF
      ENDDO

       DO JA = 1, KTRAER
         DO JL = KIDIA,KFDIA
           IF(LL_DO(JL))THEN
             IJL = JL - ILOOP
             IJA = (JA-1)*KLEV+JLK
             PTTB(IJL,JA,JLK,JK) = (ZTTP(JL,JA,JLK)+ZTTP(JL,JA,JLK+1))*0.5_JPRB
           ENDIF
         ENDDO
       ENDDO

      DO JL = KIDIA,KFDIA
        IF (LL_DO(JL)) THEN
          IJL = JL - ILOOP
          ZWD1 = PDBDT(JL,1,IJKL)*PTTB(IJL,1,JLK,JK)
          ZWD2 = PDBDT(JL,2,IJKL)*PTTB(IJL,2,JLK,JK)*PTTB(IJL,7,JLK,JK)
          ZWD3 = PDBDT(JL,3,IJKL)*PTTB(IJL,4,JLK,JK)*PTTB(IJL,8,JLK,JK)
          ZWD4 = PDBDT(JL,4,IJKL)*PTTB(IJL,5,JLK,JK)*PTTB(IJL,9,JLK,JK)
          ZWD5 = PDBDT(JL,5,IJKL)*PTTB(IJL,3,JLK,JK)
          ZWD6 = PDBDT(JL,6,IJKL)*PTTB(IJL,6,JLK,JK)
          PCNTRB(JL,IJKL,JK) = ZWD1+ZWD2+ZWD3+ZWD4+ZWD5+ZWD6
          PDISU(JL,JK)=PDISU(JL,JK)+PCNTRB(JL,IJKL,JK)
        ENDIF
      ENDDO

    ENDIF

  ENDDO
ENDDO

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('LWVDR',1,ZHOOK_HANDLE)
END SUBROUTINE LWVDR
