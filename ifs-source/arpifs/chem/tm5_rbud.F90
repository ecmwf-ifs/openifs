! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_RBUD(YGFL,KIDIA, KFDIA, KLON, KLEV, KL, KOH, KNSTEP, PAIRD, PCR2, PCR3, &
  &  PBUDJ, PBUDR, PBUDX )


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!
!*** CALCULATION OF REACTION BUDGETS due to gas-phase chemistry and photolysis
!
!--------------------------------------------------------------------------
!
!
!
!**   INTERFACE.
!     ----------
!          *TM5_bdug* IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KL    :  Current level  
! KNSTEP:  number of time steps for which individual budget terms are accumulated 
! PCR2 (KLON,NPHOTO)    : budget contribution due to photolysis           
! PCR3 (KLON,NCHEM)     : budget contribution due to chem          !
!
!
! OUTPUTS:
! -------
! PBUDJ, PBUDR, PBUDX (KLON,KLEV,NBUD_EXTRA)  : splitted tendencies due to photolysis/oxidation (kg/kg/s)
!
! LOCAL:
! -------
!
! ZFAC
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2010-10-08



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE TM5_CHEM_MODULE , ONLY : NREAC, NRR, NRJ, NSOA_BUDG, NBUD_EXTRA, NBUD_EXTRA_CHEM, CHEM_BUDG ! Define CHEM_BUDG locally?!
USE TM5_PHOTOLYSIS , ONLY : NPHOTO,JO2,JACH2O

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON, KLEV, KL, KOH
INTEGER(KIND=JPIM),INTENT(IN) :: KNSTEP
REAL(KIND=JPRB),INTENT(IN)    :: PAIRD(KLON)
REAL(KIND=JPRB),INTENT(IN)    :: PCR2(KLON,NPHOTO)
REAL(KIND=JPRB),INTENT(IN)    :: PCR3(KLON,NREAC)
REAL(KIND=JPRB),INTENT(INOUT) :: PBUDJ(KLON,KLEV,NPHOTO+NSOA_BUDG)
REAL(KIND=JPRB),INTENT(INOUT) :: PBUDR(KLON,KLEV,YGFL%NCHEM)
REAL(KIND=JPRB),INTENT(INOUT) :: PBUDX(KLON,KLEV,NBUD_EXTRA)

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

! * counters
INTEGER(KIND=JPIM) :: JL, JR,JT, JH, IOH
REAL(KIND=JPRB)    :: ZFAC


IF (LHOOK) CALL DR_HOOK('TM5_RBUD',0,ZHOOK_HANDLE )
ASSOCIATE(NCHEM=>YGFL%NCHEM, YCHEM=>YGFL%YCHEM)

IOH = KOH
 
      ! Fill reaction tendencies for reaction with OH and photolysis 
      ! arrays nrr and nrj are used to determine which species are
      ! involved in a reaction 
      !
  !    DO JL=KIDIA,KFDIA
  !      DO JT=1,NCHEM
  !        ZFAC=YCHEM(JT)%RMOLMASS/(PAIRD(JL)*KNSTEP)
  !        ! search for all corresponding photolysis rates that contribute to loss of tracer JT
  !        DO JR=1_JPIM,NPHOTO
  !          IF (NRJ(JR) == JT) THEN
  !            PBUDJ(JL,KL,JR)=PBUDJ(JL,KL,JR)+ PCR2(JL,JR)*ZFAC  !units kg/kg/s  
  !          ENDIF
  !        ENDDO
  !        ! search for reactions of species JT with OH ...
  !        DO JR=1_JPIM,NREAC
  !          IF   ((NRR(JR,1) == JT .AND. NRR(JR,2) == IOH) & 
  !       &  .OR.  (NRR(JR,2) == JT .AND. NRR(JR,1) == IOH)) THEN
  !            PBUDR(JL,KL,JT)=PBUDR(JL,KL,JT)+ PCR3(JL,JR)*ZFAC  !units  kg/kg/s 
  !         ENDIF
  !        ENDDO !nj
  !      ENDDO
  !    ENDDO
  ! alternative solution, which is supposed to be faster:
      DO JL=KIDIA,KFDIA
        DO JR = 1,NPHOTO
          ZFAC=YCHEM(NRJ(JR))%RMOLMASS/(PAIRD(JL)*KNSTEP)
          ! search for all corresponding photolysis rates that contribute to loss of tracer JT
          PBUDJ(JL,KL,JR)=PBUDJ(JL,KL,JR)+ PCR2(JL,JR)*ZFAC  !units kg/kg/s  
        ENDDO
        DO JR = 1,NREAC
          ! search for reactions of species JT with OH ...
          IF ( NRR(JR,2) == IOH  ) THEN 
           JT = NRR(JR,1)
           ! Make sure that we are only considering oxidation reactions with OH
           IF ( JT > 0 ) THEN
             ZFAC=YCHEM(JT)%RMOLMASS/(PAIRD(JL)*KNSTEP)
             PBUDR(JL,KL,JT)=PBUDR(JL,KL,JT)+ PCR3(JL,JR)*ZFAC  !units  kg/kg/s 
           ENDIF
          ENDIF
        ENDDO
      ENDDO

      !
      ! fill the BUD2D array with reaction budgets for O3 prod / loss 
      ! (terms that are not yet provided above), or other interesting 
      ! chemistry budget terms. 
      !
      DO JL=KIDIA,KFDIA
        ZFAC=1._JPRB/(PAIRD(JL)*KNSTEP)
        DO JH = 1,NBUD_EXTRA_CHEM
          JR=CHEM_BUDG(JH)
          ! * Compute budget in terms of O3 mass?
          ! ZFAC=YCHEM(IO3)%RMOLMASS/(PAIRD(JL)*KNSTEP)
          ! * More general: compute just in terms of (arbitrary) unity molar mass. Modify 
          ! * script analysing these budgets to get appropriate mass for specific reaction budget. 
          ! ZFAC=1._JPRB/(PAIRD(JL)*KNSTEP)
          PBUDX(JL,KL,JH)=PBUDX(JL,KL,JH)+ PCR3(JL,JR)*ZFAC  !units kg/kg/s 
        ENDDO
      ! ENDDO
      ! Also include some specific photolysis budgets that are not nicely included above
      ! DO JL=KIDIA,KFDIA
          ! JO2
          JH = NBUD_EXTRA-1
          JR=JO2
          ! ZFAC=1._JPRB/(PAIRD(JL)*KNSTEP)
          PBUDX(JL,KL,JH)=PBUDX(JL,KL,JH)+ PCR2(JL,JR)*ZFAC  !units kg/kg/s 
          ! JACH2O
          JH = NBUD_EXTRA
          JR=JACH2O
          ! ZFAC=1._JPRB/(PAIRD(JL)*KNSTEP)
          PBUDX(JL,KL,JH)=PBUDX(JL,KL,JH)+ PCR2(JL,JR)*ZFAC  !units kg/kg/s 
      ENDDO
      


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5_RBUD',1,ZHOOK_HANDLE )

END SUBROUTINE TM5_RBUD

