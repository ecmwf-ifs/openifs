! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE OUT_ONEGRDPT(IU06)
!
!--------------------------------------------------------------------
!
!*****OUT_ONEGRDPT** OUTPUT OF ENERGY, PEAK FREQUENCY, DIMENSIONLESS
!                    ENERGY, ROUGHNESS, FETCH  FOR ONE GRID POINT.
!
!     P.JANSSEN JUNE 2005
!
!     PURPOSE
!     -------
!             TO OUTPUT INFORMATIONS AT A FIXED GRID POINT.
!
!     INTERFACE
!     ---------
!             *CALL* *OUT_ONEGRDPT*
!
!
!     METHOD
!     ------
!             NONE
!
!     EXTERNALS
!     ---------
!             NONE
!
!     REFERENCES
!     ----------
!             NONE
!
!--------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,FFLAG    ,GFLAG    ,NFLAG     ,   &
     &            IPFGTBL  ,NIPRMOUT ,ITOBOUT  ,IRCD     ,IRU10     ,   &
     &            IRHS     ,IRTP     ,IRT1     ,IRPHIAW  ,IRPHIOC   ,   &
     &            IRTAUOC  , IRHSWS  ,IRT1WS   ,IRBATHY 
      USE YOWGRID  , ONLY : DELPHI
      USE YOWINTP  , ONLY : GOUT
      USE YOWMAP   , ONLY : NGX      ,NGY
      USE YOWPARAM , ONLY : LLUNSTR
      USE YOWPCONS , ONLY : G        ,DEG      ,ZMISS    ,EPSUS    ,    &
     &            EPSU10   ,ZPI
      USE YOWPHYS  , ONLY : XKAPPA   ,XNLEV    ,RNUM     ,ALPHAMIN
      USE YOWSTAT  , ONLY : CDATEA   ,CDTPRO
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "difdate.intfb.h"

      INTEGER(KIND=JWIM) :: IU06, IU_INTP
      INTEGER(KIND=JWIM) :: ITIME, I, J
      INTEGER(KIND=JWIM) :: IPHS, IPCD, IPU10, IPTP, IPT1, IPPHIAW, IPPHIOC, IPTAUOC
      INTEGER(KIND=JWIM) :: IPHSWS, IPT1WS, IPBATHY

      REAL(KIND=JWRB) :: CD, U10, HS, HSWS, USTAR2, USTAR, TSTAR, DSTAR, WAGEP
      REAL(KIND=JWRB) :: E, ESTAR, FMSTAR, TSTAR_0, XP, BETA_K, ALPHA_K
      REAL(KIND=JWRB) :: E_LIM, E_STAR_OBS, FP, XNUSTAR, XNU_OBS
      REAL(KIND=JWRB) :: CDSQRTINV, Z0, BETA, DFETCH, FETCHSTAR  
      REAL(KIND=JWRB) :: T10, E10, FP10, FETCH10, T_0, E_OBS 
      REAL(KIND=JWRB) :: DEPTH, PHIAW, PHIOC, TAUOC
      REAL(KIND=JWRB) :: Z0VIS
      REAL(KIND=JWRB) :: XLOGE_YV, XLOGF_YV
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: Tws, Ews, Fws  

      CHARACTER(LEN=72) :: CPATH
      CHARACTER(LEN=100) :: TITL

      LOGICAL, SAVE :: FRSTIME
      LOGICAL :: LPARAM
      LOGICAL :: LDEPTH, LPHIAW, LPHIOC, LTAUOC 

      DATA FRSTIME/.TRUE./   
!
! -------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUT_ONEGRDPT',0,ZHOOK_HANDLE)
!
!       THE SWAMP CASE HAS BEEN REDUCED TO ONE LINE OF GRID POINTS
!       NECESSARY PARAMETERS ARE OUTPUT HERE
!

      IU_INTP = 999

      IF (FRSTIME) THEN
        CPATH = 'out_intp'
        OPEN (IU_INTP,FILE=CPATH,STATUS='UNKNOWN', FORM='FORMATTED')
        FRSTIME = .FALSE.
      ENDIF
         
      CALL DIFDATE(CDATEA,CDTPRO,ITIME)

      IPHS=ITOBOUT(IRHS)
      IF(IPHS.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPCD=ITOBOUT(IRCD)
      IF(IPCD.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPU10=ITOBOUT(IRU10)
      IF(IPU10.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPTP=ITOBOUT(IRTP)
      IF(IPTP.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPT1=ITOBOUT(IRT1)
      IF(IPT1.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPHSWS=ITOBOUT(IRHSWS)
      IF(IPHSWS.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPT1WS=ITOBOUT(IRT1WS)
      IF(IPT1WS.GT.0) THEN
        LPARAM = .TRUE.
      ELSE
        LPARAM = .FALSE.
      ENDIF

      IPBATHY=ITOBOUT(IRBATHY)
      IF(IPBATHY.GT.0) THEN
        LDEPTH = .TRUE.
      ELSE
        LDEPTH = .FALSE.
      ENDIF

      IPPHIAW=ITOBOUT(IRPHIAW)
      IF(IPPHIAW.GT.0) THEN
        LPHIAW = .TRUE.
      ELSE
        LPHIAW = .FALSE.
      ENDIF

      IPPHIOC=ITOBOUT(IRPHIOC)
      IF(IPPHIOC.GT.0) THEN
        LPHIOC = .TRUE.
      ELSE
        LPHIOC = .FALSE.
      ENDIF

      IPTAUOC=ITOBOUT(IRTAUOC)
      IF(IPTAUOC.GT.0) THEN
        LTAUOC = .TRUE.
      ELSE
        LTAUOC = .FALSE.
      ENDIF

      IF (LPARAM .AND. (.NOT. LLUNSTR)) THEN
        I = MAX(1,NGX/2)
        DO J = 1,NGY
          IF (GOUT(IPHS,I,J).NE.ZMISS) THEN
            IF(LDEPTH) THEN
              DEPTH=GOUT(IPBATHY,I,J)
            ELSE
              DEPTH=999.0_JWRB
            ENDIF
!
!            TAKE RELATION BETWEEN USTAR AND U10G FROM BUILDSTRESS
!
            CD     = GOUT(IPCD,I,J)
            U10    = GOUT(IPU10,I,J) 
            HS     = GOUT(IPHS,I,J)
            USTAR2 = MAX(CD*MAX(U10**2,EPSU10**2),EPSUS)
            USTAR  = SQRT(USTAR2)

            TSTAR = G*ITIME/USTAR
            DSTAR = G*DEPTH/(USTAR**2)
            E     = HS**2/16.0_JWRB
            ESTAR = G**2*E/(USTAR**4)
            FMSTAR=USTAR/(GOUT(IPT1,I,J)*G)

            HSWS  = GOUT(IPHSWS,I,J) 
            Tws   = G*ITIME/USTAR
            E     = HSWS**2/16.0_JWRB
            Ews   = G**2*E/(USTAR**4)
            Fws   = USTAR/(GOUT(IPT1WS,I,J)*G)

            TSTAR_0 = 4.26_JWRB*10.0_JWRB**5
            XP      = 1.5_JWRB
            BETA_K  = 100.0_JWRB
            ALPHA_K = 2.25_JWRB/10000.0_JWRB
              
            E_LIM = BETA_K**2/16.0_JWRB
            E_STAR_OBS = E_LIM/(1.+TSTAR_0/TSTAR)**XP 

            FP     = 1.0_JWRB/GOUT(IPTP,I,J)
            XNUSTAR = USTAR*FP/G
            XNU_OBS = (ALPHA_K/E_STAR_OBS)**(1.0_JWRB/3.0_JWRB)

            WAGEP = G / (ZPI*FP*USTAR)

            CDSQRTINV = MIN(1./SQRT(CD),50.0_JWRB)
            Z0        = XNLEV/(EXP(XKAPPA*CDSQRTINV)-1.0_JWRB)
            Z0VIS     = RNUM/USTAR 
            BETA      = MAX(G*(Z0-Z0VIS)/USTAR2,ALPHAMIN)

            DFETCH = (NGY-J+1)*DELPHI
            FETCHSTAR = G*DFETCH/USTAR2

            T10 = G*ITIME/U10
            E10 = G**2*E/U10**4
            FP10 = U10*FP/G
            FETCH10 = G*DFETCH/U10**2

            T_0 = 1.93_JWRB*10.0_JWRB**4
            XP      = 1.5_JWRB
            BETA_K  = 0.22_JWRB
            E_LIM = BETA_K**2/16.0_JWRB
            E_OBS = E_LIM/(1.+T_0/T10)**XP
 
            IF(LPHIAW) THEN
              PHIAW=GOUT(IPPHIAW,I,J)
            ELSE
              PHIAW=3.5_JWRB
            ENDIF

            IF(LPHIOC) THEN
!             make it positive for comparison with PHIAW
              PHIOC=-GOUT(IPPHIOC,I,J)
            ELSE
              PHIOC=3.5_JWRB
            ENDIF

            IF(LTAUOC) THEN
              TAUOC=GOUT(IPTAUOC,I,J)
            ELSE
              TAUOC=1.0_JWRB
            ENDIF

            WRITE(IU_INTP,60) NGY-J+1,DEPTH,ITIME/3600.0_JWRB,          &
     &                       LOG10(TSTAR),LOG10(FETCHSTAR),HS,          &
     &                       FP,LOG10(ESTAR),LOG10(E_STAR_OBS),         &
     &                       LOG10(XNUSTAR),LOG10(XNU_OBS),U10,         &
     &                       USTAR,CD,BETA,T10,FETCH10,E10,             &
     &                       TSTAR, FMSTAR, DSTAR, ESTAR,               &
     &                       PHIAW, PHIOC, TAUOC, Tws, Fws, Ews, HSWS,  &
     &                       WAGEP 

            WRITE(IU06,61) NGY-J+1,DEPTH,ITIME/3600.0_JWRB,             &
     &                     TSTAR,FETCHSTAR,HS,FP,ESTAR,                 &
     &                     E_STAR_OBS,XNUSTAR,XNU_OBS,U10,              &
     &                     USTAR,CD,BETA,T10,FETCH10,E10,               &
     &                     TSTAR, FMSTAR, DSTAR, ESTAR,                 &
     &                     PHIAW, PHIOC, TAUOC, Tws, Fws, Ews, HSWS,    &
     &                     WAGEP 

          ENDIF
        ENDDO
        CALL FLUSH(IU_INTP)
        CALL FLUSH(IU06)
      ELSE IF (.NOT. LLUNSTR) THEN
        WRITE(IU06,*) '***************************************'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. OUT_ONEGRDPT   *'
        WRITE(IU06,*) '*                                     *'
        WRITE(IU06,*) '*  MISSING SOME OUTPUT VARIABLES !    *'
        WRITE(IU06,*) '***************************************'
        CALL ABORT1
      ENDIF
   60 FORMAT(I4,F7.2,F7.2,2E10.3,2F8.3,4E12.3,F6.1,F7.3,                &
     &       2F12.5,2E10.3,15(1x,F15.5))
   61 FORMAT(I4,F7.2,F7.2,2E10.3,2F8.3,2F10.3,2F8.5,F6.1,F7.3,          &
     &       2F12.5,2E10.3,15(1x,F15.5))
!!
!     YOUNG-VERHAGEN LIMITS
! 
!     XLOGE_YV = 1.3*XD+LOG10(1.06/1000.)
!     XLOGF_YV = -0.375*XD+LOG10(0.2)
!
!     USTAR-SCALING
!
!     XLOGE_YV = 1.3*XD+LOG10(111.*1.06/1000.)
!     XLOGF_YV = -0.375*XD+LOG10(0.2*0.43)
!
!     CERC
!
!      XLOGE_YV = 1.5*XD+LOG10(0.04)
!      XLOGF_YV = -0.375*XD+LOG10(0.2*0.43)
!
! -------------------------------------------------------------------
 
      IF (LHOOK) CALL DR_HOOK('OUT_ONEGRDPT',1,ZHOOK_HANDLE)

      END SUBROUTINE OUT_ONEGRDPT
