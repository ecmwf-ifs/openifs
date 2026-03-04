! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE N2O_SOLVER(YGFL,KIDIA,KFDIA,KLON,KN2O,PDT,PJ_VAL,PY0,PO1D,PY,PRSF1,PTEMP)

!**   DESCRIPTION 
!     ----------
!
!   Simple solver for N2O chemistry  
!--------------------------------------------------------------------------
!   Eulerian backward Iteration
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *N2O_SOLVER* IS CALLED FROM *CHEM_N2O*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KN2O  : Index of N2O tracer
! PDT   :  Time step in seconds 
! PJ_VAL  (KLON,2)         : N2O and O3_O1D photolysis rate
! PY0(KLON,NCHEM)        : initial concentrations of tracers         
! PO1D(KLON)             : est. background concentration of O1D         
! PY(KLON,NCHEM)         : final concentrations of tracers         
! PRSF1 (KLON)            : FULL-LEVEL PRESSURE           (Pa)
! PTEMP (KLON)            : Temperature
!
!
! OUTPUTS:
! -------
! PY (KLON,NCHEM)        : final   concentrations of tracers          
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-05-06



USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KN2O
REAL(KIND=JPRB),INTENT(IN)    :: PDT
REAL(KIND=JPRB),INTENT(IN)    :: PJ_VAL(KLON,2)   
REAL(KIND=JPRB),INTENT(IN)    :: PY0(KLON,YGFL%NCHEM)   ! initial concentrations
REAL(KIND=JPRB),INTENT(IN)    :: PO1D(KLON)             ! Background O1D field
REAL(KIND=JPRB),INTENT(OUT)   :: PY(KLON,YGFL%NCHEM)    ! final concentrations
REAL(KIND=JPRB),INTENT(IN)    :: PRSF1(KLON)            ! full level pressure
REAL(KIND=JPRB),INTENT(IN)    :: PTEMP(KLON)            ! temperature

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


REAL(KIND=JPRB)       :: ZRR_O1D_N2,ZRR_O1D_N2O_1,ZRR_O1D_N2O_2, &
    & ZTEMP, ZCFACTOR, ZN2, ZPN2O,ZXLN2O
! * counters
INTEGER(KIND=JPIM) :: JL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('N2O_SOLVER',0,ZHOOK_HANDLE )




    ! --- Long living compounds
    DO JL=KIDIA,KFDIA
       ! Get temperature
       ZTEMP=PTEMP(JL)

       ! Air density (molec/cm3) 
       ZCFACTOR = 7.24291E16_JPRB*PRSF1(JL)/ZTEMP 
       ZN2 = 0.781*ZCFACTOR    ! N2 number density

       ! Compute reaction rate
       ZRR_O1D_N2   = TBODY(ZCFACTOR,2.8E-36_JPRB,0.9E0_JPRB,0.E0_JPRB,0.E0_JPRB) 
       ZRR_O1D_N2O_1= DARR(4.63E-11_JPRB, 20.E0_JPRB)
       ZRR_O1D_N2O_2= DARR(7.25E-11_JPRB, 20.E0_JPRB)

       ! Compute updated N2O concentration
       ZPN2O  =  PO1D(JL)*ZN2*ZRR_O1D_N2
       ZXLN2O = (ZRR_O1D_N2O_1+ZRR_O1D_N2O_2)*PO1D(JL) + PJ_VAL(JL,1)
       PY(JL,KN2O) = (PY0(JL,KN2O)+ZPN2O*PDT)/(1._JPRB+ZXLN2O*PDT)

    ENDDO   !JL



IF (LHOOK) CALL DR_HOOK('N2O_SOLVER',1,ZHOOK_HANDLE )

CONTAINS



! Begin INLINED Rate Law Functions


    FUNCTION DARR(PA,PB)
      IMPLICIT NONE

      REAL(KIND=JPRB),INTENT(IN) :: PA, PB
      REAL(KIND=JPRB) :: DARR
      
      DARR = PA * EXP(PB/ZTEMP)
    END FUNCTION DARR
!
!
    FUNCTION TBODY( PSN, PK00, PN, PK0INF, PM )
      IMPLICIT NONE

      REAL(KIND=JPRB),INTENT(IN) ::PSN, PK00, PK0INF, PN, PM
      REAL(KIND=JPRB) :: TBODY,ZK0,ZKINF,ZX1
      
      ZK0=PK00*PSN*(300./ZTEMP)**PN
      IF (PK0INF > 0.) THEN
        ZKINF=PK0INF*(300./ZTEMP)**PM
        ZX1=1.+(LOG10(ZK0/ZKINF))**2
        ZX1=1./ZX1
        TBODY=(ZK0/(1.+ZK0/ZKINF))*0.6**ZX1
      ELSE
        TBODY=ZK0
      ENDIF
    END FUNCTION TBODY
    
    
END SUBROUTINE N2O_SOLVER
