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

SUBROUTINE SIMPLICO(KM,KSMAX,KFLEV,KFLSUR,PALPHA,PDENIM,&
 & PFPLUS,PFMINUS,PSIVP,PRLAPDI,PBDT2,PY,PX)  

!**** *SIMPLICO* - Resolution of the set of complex pentadiagonal systems
!                  generated in the case of implicit Coriolis terms

!     Purpose.    Resolution of a set of complex pentadiagonal systems.
!     --------

!**   Interface.
!     ----------
!        *CALL* *SIMPLICO(KM,KSMAX,KFLEV,KFLSUR,PALPHA,PDENIM,
!                         PFPLUS,PFMINUS,PSIVP,PRLAPDI,PBDT2,PY,PX)

!        Explicit arguments :
!        --------------------
!         KM :        - Zonal wavenumber                            (input)
!         KSMAX :     - Truncation limit                            (input)
!         KFLEV :     - Number of levels                            (input)
!         KFLSUR :    - Surdimension corresponding to KFLEV         (input)
!         PALPHA :    - Help array                                  (input)
!         PDENIM :    - Help array                                  (input)
!         PFPLUS :    - Help array                                  (input)
!         PFMINUS :   - Help array                                  (input)
!         PSIVP :     - Eigenvalues of vertical structure matrix    (input)
!         PRLAPDI :   - Eigenvalues of Laplacian                    (input)
!         PBDT2 :     - (BETA*DT)**2                                (input)
!         PY:         - known vector.                               (input)
!         PX:         - unknown vector.                             (output)

!        Implicit arguments :
!        --------------------
!        none.

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None.
!        Called by SPC.  

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        C. TEMPERTON 94/06/28

!     Modifications.
!     --------------
!       L. Isaksen 96-12-12 : Spectral data layout changed
!       M.Hamrud      01-Oct-2003 CY28 Cleaning
!       N. Wedi    2008-01-04 : Add KM=0 case.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAX 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLSUR 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHA(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDENIM(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFPLUS(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFMINUS(KM:KSMAX+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSIVP(KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRLAPDI(0:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBDT2 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PY(KFLSUR,2,KM:KSMAX) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PX(KFLSUR,2,KM:KSMAX) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZWORK(KFLSUR,2,KM:KSMAX)

INTEGER(KIND=JPIM) :: IHI, ILO, ITIMES, JLEV, JLPAR, JN

REAL(KIND=JPRB) :: ZDELT, ZDELTA, ZGAM, ZGAMMA, ZGI, ZGR, ZMU,&
 & ZNU, ZSUBI, ZSUBR, ZSUPI, ZSUPR, ZUNDER, ZWI, ZWR  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIMPLICO',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    FILL IN ZERO PARTS IF NEEDED
!              ----------------------------

IF (KM == 0) THEN
  PY(:,2,KM:KSMAX)=0.0_JPRB
  PX(:,2,KM:KSMAX)=0.0_JPRB
ENDIF

!     Solve "even" tridiagonal systems, then "odd" ones

IF (KM == KSMAX) THEN
  ITIMES=1
ELSE
  ITIMES=2
ENDIF

!     ------------------------------------------------------------------

!*       2.    INVERT PENTADIAGONAL SYSTEM 
!              ---------------------------

! operators are real only !
IF( KM == 0 ) THEN

  ODD_EVEN_LOOP_KM0 : DO JLPAR=1,ITIMES

    !*   2.1   FORWARD SWEEP FOR KM=0
    !          ----------------------

    ILO=JLPAR-1
    ZGAM=1.0_JPRB+PFPLUS(ILO)*PFMINUS(ILO+1)
    IF (JLPAR == 2) THEN
      ZGAM=ZGAM+PFMINUS(ILO)*PFPLUS(ILO-1)
    ENDIF
    ZSUPR=PFPLUS(ILO)*PFPLUS(ILO+1)
    
    DO JLEV=1,KFLEV
      ZGAMMA=ZGAM-PBDT2*PSIVP(JLEV)*PRLAPDI(ILO)
      ZUNDER=1.0_JPRB/(ZGAMMA*ZGAMMA)
      ZMU=-ZGAMMA*ZUNDER
      PX(JLEV,1,ILO)=ZMU*PY(JLEV,1,ILO)
      ZWORK(JLEV,1,ILO)=ZMU*ZSUPR
    ENDDO
    
    DO JN=ILO+2,KSMAX,2
      ZGAM=1.0_JPRB+PFMINUS(JN)*PFPLUS(JN-1)&
       & +PFPLUS(JN)*PFMINUS(JN+1)  
      ZSUPR=PFPLUS(JN)*PFPLUS(JN+1)
      ZSUBR=PFMINUS(JN-1)*PFMINUS(JN)
      DO JLEV=1,KFLEV
        ZGAMMA=ZGAM-PBDT2*PSIVP(JLEV)*PRLAPDI(JN)&
         & +ZSUBR*ZWORK(JLEV,1,JN-2)
        ZUNDER=1.0_JPRB/(ZGAMMA*ZGAMMA)
        ZMU=-ZGAMMA*ZUNDER
        ZGR=ZSUBR*PX(JLEV,1,JN-2)+PY(JLEV,1,JN)
        PX(JLEV,1,JN)=ZMU*ZGR
        ZWORK(JLEV,1,JN)=ZMU*ZSUPR
      ENDDO
    ENDDO
    
    !*   2.2   BACKWARD SWEEP FOR KM=0
    !          -----------------------

    IHI=KSMAX-MOD(KSMAX-KM+JLPAR-1,2)
    DO JLEV=1,KFLEV
      PX(JLEV,1,IHI)=-PX(JLEV,1,IHI)
    ENDDO
    
    DO JN=IHI-2,ILO,-2
      DO JLEV=1,KFLEV
        ZWR=ZWORK(JLEV,1,JN)
        PX(JLEV,1,JN)=ZWR*PX(JLEV,1,JN+2)-PX(JLEV,1,JN)
      ENDDO
    ENDDO
    
  ENDDO ODD_EVEN_LOOP_KM0

ELSE

  ODD_EVEN_LOOP : DO JLPAR=1,ITIMES

    !*   2.3   FORWARD SWEEP FOR KM>0
    !          ----------------------

    ILO=KM+JLPAR-1
    ZGAM=1.0_JPRB+PDENIM(ILO+1)*PFPLUS(ILO)*PFMINUS(ILO+1)
    ZDELTA=-PALPHA(ILO)+PDENIM(ILO+1)*PALPHA(ILO+1)*PFPLUS(ILO)*PFMINUS(ILO+1)
    IF (JLPAR == 2) THEN
      ZGAM=ZGAM+PDENIM(ILO-1)*PFMINUS(ILO)*PFPLUS(ILO-1)
      ZDELTA=ZDELTA+PDENIM(ILO-1)*PALPHA(ILO-1)*PFMINUS(ILO)*PFPLUS(ILO-1)
    ENDIF
    ZSUPR=PDENIM(ILO+1)*PFPLUS(ILO)*PFPLUS(ILO+1)
    ZSUPI=PALPHA(ILO+1)*ZSUPR
    
    DO JLEV=1,KFLEV
      ZGAMMA=ZGAM-PBDT2*PSIVP(JLEV)*PRLAPDI(ILO)
      ZUNDER=1.0_JPRB/(ZGAMMA*ZGAMMA+ZDELTA*ZDELTA)
      ZMU=-ZGAMMA*ZUNDER
      ZNU=ZDELTA*ZUNDER
      PX(JLEV,1,ILO)=ZMU*PY(JLEV,1,ILO)-ZNU*PY(JLEV,2,ILO)
      PX(JLEV,2,ILO)=ZNU*PY(JLEV,1,ILO)+ZMU*PY(JLEV,2,ILO)
      ZWORK(JLEV,1,ILO)=ZMU*ZSUPR-ZNU*ZSUPI
      ZWORK(JLEV,2,ILO)=ZNU*ZSUPR+ZMU*ZSUPI
    ENDDO

    DO JN=ILO+2,KSMAX,2
      ZGAM=1.0_JPRB+PDENIM(JN-1)*PFMINUS(JN)*PFPLUS(JN-1)&
       & +PDENIM(JN+1)*PFPLUS(JN)*PFMINUS(JN+1)  
      ZDELT=-PALPHA(JN)&
       & +PDENIM(JN-1)*PALPHA(JN-1)*PFMINUS(JN)*PFPLUS(JN-1)&
       & +PDENIM(JN+1)*PALPHA(JN+1)*PFPLUS(JN)*PFMINUS(JN+1)  
      ZSUPR=PDENIM(JN+1)*PFPLUS(JN)*PFPLUS(JN+1)
      ZSUPI=PALPHA(JN+1)*ZSUPR
      ZSUBR=PDENIM(JN-1)*PFMINUS(JN-1)*PFMINUS(JN)
      ZSUBI=PALPHA(JN-1)*ZSUBR
      DO JLEV=1,KFLEV
        ZGAMMA=ZGAM-PBDT2*PSIVP(JLEV)*PRLAPDI(JN)&
         & +ZSUBR*ZWORK(JLEV,1,JN-2)-ZSUBI*ZWORK(JLEV,2,JN-2)  
        ZDELTA=ZDELT+ZSUBI*ZWORK(JLEV,1,JN-2)+ZSUBR*ZWORK(JLEV,2,JN-2)
        ZUNDER=1.0_JPRB/(ZGAMMA*ZGAMMA+ZDELTA*ZDELTA)
        ZMU=-ZGAMMA*ZUNDER
        ZNU=ZDELTA*ZUNDER
        ZGR=ZSUBR*PX(JLEV,1,JN-2)-ZSUBI*PX(JLEV,2,JN-2)+PY(JLEV,1,JN)
        ZGI=ZSUBI*PX(JLEV,1,JN-2)+ZSUBR*PX(JLEV,2,JN-2)+PY(JLEV,2,JN)
        PX(JLEV,1,JN)=ZMU*ZGR-ZNU*ZGI
        PX(JLEV,2,JN)=ZNU*ZGR+ZMU*ZGI
        ZWORK(JLEV,1,JN)=ZMU*ZSUPR-ZNU*ZSUPI
        ZWORK(JLEV,2,JN)=ZNU*ZSUPR+ZMU*ZSUPI
      ENDDO
    ENDDO

    !*   2.4   BACKWARD SWEEP FOR KM>0
    !          -----------------------

    IHI=KSMAX-MOD(KSMAX-KM+JLPAR-1,2)
    DO JLEV=1,KFLEV
      PX(JLEV,1,IHI)=-PX(JLEV,1,IHI)
      PX(JLEV,2,IHI)=-PX(JLEV,2,IHI)
    ENDDO
    
    DO JN=IHI-2,ILO,-2
      DO JLEV=1,KFLEV
        ZWR=ZWORK(JLEV,1,JN)
        ZWI=ZWORK(JLEV,2,JN)
        PX(JLEV,1,JN)=ZWR*PX(JLEV,1,JN+2)-ZWI*PX(JLEV,2,JN+2)-PX(JLEV,1,JN)
        PX(JLEV,2,JN)=ZWI*PX(JLEV,1,JN+2)+ZWR*PX(JLEV,2,JN+2)-PX(JLEV,2,JN)
      ENDDO
    ENDDO
    
  ENDDO ODD_EVEN_LOOP
  
ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SIMPLICO',1,ZHOOK_HANDLE)
END SUBROUTINE SIMPLICO

