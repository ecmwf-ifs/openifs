! (C) Copyright 2004- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
SUBROUTINE VDFDIFC(KIDIA,KFDIA,KLON,KLEV,KTOP,KTRAC,&
                 & PTMST,PCM1,PTENC,PAPHM1,PCFH,PCFLX,PWD,PCFLXT)  
!     ------------------------------------------------------------------

!**   *VDFDIFC* - Does the implicit computation for diffusion of tracers
!     A. Beljaars       ECMWF    18-03-2004

!     Revisions: 
!     A. Beljaars    Jan-2014  Cleanup of confusing K-scaling + deposition vel. 

!     PURPOSE
!     -------

!     Solve the tridiagonal equations for the diffusion of tracers. 

!     INTERFACE
!     ---------

!     *VDFDIFC* is called by *VDFMAIN*

!     INPUT PARAMETERS (INTEGER):

!     *KIDIA*        Start point
!     *KFDIA*        End point
!     *KLEV*         Number of levels
!     *KLON*         Number of grid points per packet
!     *KTOP*         Index for boundary layer top
!     *KTRAC*        Number of tracers
!     *PCFLX*        Surface flux boundary condition for tracers    (kg/m2s)
!     *PWD*          Rho*Vd (density * deposition velocity)         (kg/m2s)                           (kg/m2s)

!     INPUT PARAMETERS (REAL):

!     *PTMST*        Time step                                      (s)
!     *PCM1*         Tracer concentration at T-1                    (kg/kg)
!     *PAPHM1*       Pressure AT T-1                                (Pa)
!     *PCFH*         Rho*K/dz (K-star in doc.)                      (kg/m2s)

!     UPDATED PARAMETERS (REAL):

!     *PTENC*        Tendency of tracer concentration               (1/s)

!     OUTPUT PARAMETERS (REAL):

!     *PCFLXT*       Total tracer flux (PCFLX+PWD*tracer concentr.) (kg/m2s)

!     METHOD
!     ------

!     *LU*-decomposition (downward scan), followed by 
!     back substitution (upward scan).

!     EXTERNALS.
!     ----------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,  JPRB
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOMCST   , ONLY : RG
IMPLICIT NONE


INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTOP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTRAC 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCM1(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PTENC(KLON,KLEV,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPHM1(KLON,KLEV+1) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFH(KLON,KLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCFLX(KLON,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWD  (KLON,KTRAC) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCFLXT(KLON,KTRAC) 

!*    LOCAL STORAGE
!     ----- -------

REAL(KIND=JPRB) ::  ZTCOE(KLON),ZEBSH(KLON,KLEV), ZQDP1(KLON)
REAL(KIND=JPRB) ::  ZDISC(KLON),ZFAC(KLON),ZCDIF(KLON,KLEV,KTRAC)

INTEGER(KIND=JPIM) :: ILEVM1,ITOPP1,JK,JL,JTR

REAL(KIND=JPRB) :: ZACL,ZBCL,ZALF,ZQDP,ZCONS13,ZCONS1,ZTMST
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('VDFDIFC',0,ZHOOK_HANDLE)

ZALF=1.0_JPRB     ! Implicitness coefficient for tracer diffusion 
ZTMST=ZALF*PTMST  ! Pseudo-time step (e.g. extrapolated in time with ZALF=1.5)
ZCONS1=RG*ZTMST

ZCONS13=1.0_JPRB/ZTMST
ILEVM1=KLEV-1
ITOPP1=KTOP+1

!*         1.1      Setting of right hand sides.

DO JTR=1,KTRAC
  DO JK=KTOP,KLEV
    DO JL=KIDIA,KFDIA
      ZCDIF(JL,JK,JTR)= PCM1(JL,JK,JTR)+ZTMST*PTENC(JL,JK,JTR)
    ENDDO
  ENDDO
ENDDO


!*         2.1     Top layer elimination.

DO JL=KIDIA,KFDIA
  ZTCOE(JL)=PCFH(JL,KTOP)
  ZQDP=ZCONS1/(PAPHM1(JL,ITOPP1)-PAPHM1(JL,KTOP))
  ZDISC(JL)=1.0_JPRB/(1.0_JPRB+PCFH(JL,KTOP)*ZQDP)
  ZEBSH(JL,KTOP)=ZDISC(JL)*(PCFH(JL,KTOP)*ZQDP)
ENDDO
DO JTR=1,KTRAC
  DO JL=KIDIA,KFDIA
    ZCDIF(JL,KTOP,JTR)=ZDISC(JL)*ZCDIF(JL,KTOP,JTR)
  ENDDO
ENDDO

!*         2.2     Elimination for middle layers.

DO JK=ITOPP1,ILEVM1
  DO JL=KIDIA,KFDIA
    ZQDP=ZCONS1/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
    ZFAC(JL)=ZTCOE(JL)*ZQDP
    ZTCOE(JL)=PCFH(JL,JK)
    ZDISC(JL)=1.0_JPRB/(1.0_JPRB+ZFAC(JL)*(1.0_JPRB-ZEBSH(JL,JK-1))+PCFH(JL,JK)*ZQDP)
    ZEBSH(JL,JK)=ZDISC(JL)*(PCFH(JL,JK)*ZQDP)
  ENDDO
  DO JTR=1,KTRAC
    DO JL=KIDIA,KFDIA
      ZCDIF(JL,JK,JTR)=ZDISC(JL)*(ZCDIF(JL,JK,JTR)+ZFAC(JL)*ZCDIF(JL,JK-1,JTR))
    ENDDO
  ENDDO
ENDDO

!*         2.3     Bottom layer elimination (preparation) 

JK=KLEV
DO JL=KIDIA,KFDIA
  ZQDP1(JL)=ZCONS1/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
  ZFAC(JL)=ZTCOE(JL)*ZQDP1(JL)
ENDDO

!*         3.     Linear relation between lowest 
!                 model level concentrations and fluxes.
!                 ZACL : Acl in Cl=Acl*Jc+Bcl                (m2s/kg)
!                 ZBCL : Bcl in Cl=Acl*Jc+Bcl                (kg/kg)
!                 Compute lowest model level concentrations and fluxes

DO JTR=1,KTRAC
  DO JL=KIDIA,KFDIA
    ZDISC(JL)=1.0_JPRB/(1.0+ZFAC(JL)*(1.0_JPRB-ZEBSH(JL,JK-1))+PWD(JL,JTR)*ZQDP1(JL))
    ZACL=-ZDISC(JL)*ZQDP1(JL)
    ZBCL=ZDISC(JL)*(ZCDIF(JL,JK,JTR)+ZFAC(JL)*ZCDIF(JL,JK-1,JTR))
    ZCDIF(JL,KLEV,JTR)=ZACL*PCFLX(JL,JTR)+ZBCL
    PCFLXT(JL,JTR)=PCFLX(JL,JTR)+PWD(JL,JTR)*ZCDIF(JL,KLEV,JTR)
  ENDDO
ENDDO

!*         4.     Back-substitution.

DO JTR=1,KTRAC
  DO JK=ILEVM1,KTOP,-1
    DO JL=KIDIA,KFDIA
      ZCDIF(JL,JK,JTR)=ZCDIF(JL,JK,JTR)+ZEBSH(JL,JK)*ZCDIF(JL,JK+1,JTR)
    ENDDO
  ENDDO
ENDDO

!*         5.    Update tendency (ZCDIF is C-hat i.e. at the end of ZTMST)

DO JTR=1,KTRAC
  DO JK=KTOP,KLEV
    DO JL=KIDIA,KFDIA
      PTENC(JL,JK,JTR)=(ZCDIF(JL,JK,JTR)-PCM1(JL,JK,JTR))*ZCONS13
    ENDDO
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('VDFDIFC',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFC





