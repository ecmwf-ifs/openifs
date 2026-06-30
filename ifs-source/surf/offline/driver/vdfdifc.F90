SUBROUTINE VDFDIFC(KIDIA,KFDIA,KLON,KLEV,KTOP,KTRAC,&
                 & PTMST,PCM1,PTENC,PAPHM1,PCFH,PCFLX)  
! (C) Copyright 2004- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

!**   *VDFDIFC* - Does the implicit computation for diffusion of tracers

!     A. Beljaars       ECMWF    18-03-2004

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
!     *PCFLX*        Surface flux boundary condition for tracers       (kg/m2s)

!     INPUT PARAMETERS (REAL):

!     *PTMST*        Time step                                         (s)
!     *PCM1*         Tracer concentration at T-1                       (kg/kg)
!     *PAPHM1*       Pressure AT T-1                                   (Pa)
!     *PCFH*         Prop. to exch. coeff. (K-star in doc.)

!     UPDATED PARAMETERS (REAL):

!     *PTENC*        Tendency of tracer concentration                  (1/s)

!     METHOD
!     ------

!     *LU*-decomposition (downward scan), followed by 
!     back substitution (upward scan).

!     EXTERNALS.
!     ----------

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM,  JPRB   ,JPRD
USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK

USE YOEVDF   , ONLY : RVDIFTS
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

!*    LOCAL STORAGE
!     ----- -------

REAL(KIND=JPRB) ::  ZTCOE(KLON),ZEBSH(KLON,KLEV)
REAL(KIND=JPRB) ::   ZACL(KLON),ZBCL(KLON,KTRAC),ZDISC(KLON),ZFAC(KLON),ZCDIF(KLON,KLEV,KTRAC)

INTEGER(KIND=JPIM) :: ILEVM1,ITOPP1,JK,JL,JTR

REAL(KIND=JPRB) :: ZALF,ZQDP,ZTPFAC2,ZCOEF2,ZCONS13 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('VDFDIFC',0,ZHOOK_HANDLE)
ZALF=RVDIFTS
ZTPFAC2=1.0_JPRB/RVDIFTS
ZCONS13=1.0_JPRB/PTMST
ILEVM1=KLEV-1
ITOPP1=KTOP+1

!*         1.1      Setting of right hand sides.

DO JTR=1,KTRAC
  DO JK=KTOP,KLEV
    DO JL=KIDIA,KFDIA
      ZCDIF(JL,JK,JTR)= ZTPFAC2*PCM1(JL,JK,JTR)+PTMST*PTENC(JL,JK,JTR)
    ENDDO
  ENDDO
ENDDO


!*         2.1     Top layer elimination.

DO JL=KIDIA,KFDIA
  ZTCOE(JL)=PCFH(JL,KTOP)
  ZQDP=1.0_JPRB/(PAPHM1(JL,ITOPP1)-PAPHM1(JL,KTOP))
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
    ZQDP=1.0_JPRB/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
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

!*         2.3     Bottom layer, linear relation between lowest  
!                  model level concentrations and fluxes.
!                  ZACL : Acl in Cl=Acl*Jc+Bcl                (m2s/kg)
!                  ZBCL : Bcl in Cl=Acl*Jc+Bcl                (kg/kg)


JK=KLEV
DO JL=KIDIA,KFDIA
  ZQDP=1.0_JPRB/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
  ZFAC(JL)=ZTCOE(JL)*ZQDP
  ZDISC(JL)=1.0_JPRB/(1.0_JPRB+ZFAC(JL)*(1.0_JPRB-ZEBSH(JL,JK-1)))
  ZACL(JL)=-ZDISC(JL)*RG*PTMST*ZQDP*ZALF
ENDDO
DO JTR=1,KTRAC
  DO JL=KIDIA,KFDIA
    ZBCL(JL,JTR)=ZDISC(JL)*(ZCDIF(JL,JK,JTR)+ZFAC(JL)*ZCDIF(JL,JK-1,JTR))*ZALF
  ENDDO
ENDDO



!*         3.     Compute lowest model level concentrations 

ZCOEF2=1.0_JPRB/ZALF
DO JTR=1,KTRAC
  DO JL=KIDIA,KFDIA
    ZCDIF(JL,KLEV,JTR)=ZCOEF2*(ZACL(JL)*PCFLX(JL,JTR)+ZBCL(JL,JTR))
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

!*         5.    Update tendency (ZCDIF is C-hat divided by alpha in kg/kg)


DO JTR=1,KTRAC
  DO JK=KTOP,KLEV
    DO JL=KIDIA,KFDIA
      PTENC(JL,JK,JTR)=(ZCDIF(JL,JK,JTR)-ZTPFAC2*PCM1(JL,JK,JTR))*ZCONS13
    ENDDO
  ENDDO
ENDDO


IF (LHOOK) CALL DR_HOOK('VDFDIFC',1,ZHOOK_HANDLE)
END SUBROUTINE VDFDIFC





