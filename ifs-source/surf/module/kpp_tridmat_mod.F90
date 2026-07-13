MODULE KPP_TRIDMAT_MOD
CONTAINS
SUBROUTINE KPP_TRIDMAT &
  & ( KIDIA    ,KFDIA    ,KLON     ,KLEVO    ,KLEVOP1  ,&
  &   LDKPPCAL ,PCU      ,PCC      ,PCL      ,PRHS     ,&
  &   PYO      ,YDOCEAN_ML,&
  &   PYN )

! (C) Copyright 1994- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
! Purpose :
! -------
!   This routine solves tridiagonal matrix with Thomas algorithm.

! Interface :
! ---------
!   Call *KPP_TRIDMAT* from *KPP_OCNINT*

! Method :
! ------

! Externals :
! ---------

! Reference :
! ---------

! Modifications :
! -------------
!     06-Jun-1994  Bill Large
!            2002  Steve Woolnough, Reading Univ.
!     07-Oct-2008  Yuhei Takaya,    E.C.M.W.F.    Implemented to IFS.
! End Modifications :
!---------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
!USE YOMLUN1S , ONLY : NULOUT   
USE YOS_OCEAN_ML , ONLY : TOCEAN_ML
USE ABORT_SURF_MOD
IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN) :: KLON
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVO
INTEGER(KIND=JPIM),INTENT(IN) :: KLEVOP1
REAL(KIND=JPRB),INTENT(IN)    :: PCU (KIDIA:KFDIA,KLEVO)  ! upper coef for k-1 
                                                          ! on k of tridmtx 
REAL(KIND=JPRB),INTENT(IN)    :: PCC (KIDIA:KFDIA,KLEVO)  ! central .. (k  ) ..
REAL(KIND=JPRB),INTENT(IN)    :: PCL (KIDIA:KFDIA,KLEVO)  ! lower ..   (k-1) ..
REAL(KIND=JPRB),INTENT(IN)    :: PRHS (KIDIA:KFDIA,KLEVO) ! right hand side
REAL(KIND=JPRB),INTENT(IN)    :: PYO(KIDIA:KFDIA,KLEVOP1) ! old field
REAL(KIND=JPRB),INTENT(OUT)   :: PYN(KIDIA:KFDIA,KLEVOP1)! new field

LOGICAL,INTENT(IN)    :: LDKPPCAL(KLON)

TYPE(TOCEAN_ML),INTENT(IN) :: YDOCEAN_ML

INTEGER(KIND=JPIM)    :: JZ                ! loop control
INTEGER(KIND=JPIM)    :: JL                ! loop control
REAL(KIND=JPRB)       :: ZGAM(KLEVO)       ! array for tri-diagonal solver
REAL(KIND=JPRB)       :: ZBET              ! dummy variable for ...
REAL(KIND=JPHOOK)       :: ZHOOK_HANDLE

!---------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('KPP_TRIDMAT_MOD:KPP_TRIDMAT',0,ZHOOK_HANDLE)
ASSOCIATE(REPS_KPP=>YDOCEAN_ML%REPS_KPP)


DO JL = KIDIA, KFDIA
  IF ( LDKPPCAL(JL) ) THEN
    ZBET      = PCC(JL,1)
    PYN(JL,1) = PRHS(JL,1) / ZBET    

    DO JZ = 2, KLEVO

      ZGAM(JZ) = PCL(JL,JZ-1) / ZBET
      ZBET     = PCC(JL,JZ) - PCU(JL,JZ)*ZGAM(JZ)

      IF(ABS(ZBET) < REPS_KPP ) THEN
!        WRITE(NULOUT,*)'* TRIDMAT: ALGORITHM FOR SOLVING TRIDIAG MATRIX FAILS'
!        WRITE(NULOUT,*)'* ZBET=',ZBET
!        WRITE(NULOUT,*)'* JL=',JL
!        WRITE(NULOUT,*)'* JZ-1=',JZ-1,' PCC=',PCC(JL,JZ-1),'PCL=',PCL(JL,JZ-1)
!        WRITE(NULOUT,*)'* JZ=',JZ,' PCC=',PCC(JL,JZ),' PCU=',PCU(JL,JZ),&
!                      &' ZGAM=',ZGAM(JZ)
        CALL ABORT_SURF('TRIDMAT: ALGORITHM FOR SOLVING TRIDIAG MATRIX FAILS')
        ZBET = REPS_KPP
      ENDIF
      PYN(JL,JZ) = (PRHS(JL,JZ) - PCU(JL,JZ)*PYN(JL,JZ-1) ) / ZBET
    ENDDO

    DO JZ = KLEVO-1, 1, -1
      PYN(JL,JZ) = PYN(JL,JZ) - ZGAM(JZ+1)*PYN(JL,JZ+1)
    ENDDO   
    PYN(JL,KLEVO+1) = PYO(JL,KLEVO+1)

  ENDIF
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('KPP_TRIDMAT_MOD:KPP_TRIDMAT',1,ZHOOK_HANDLE)

END SUBROUTINE KPP_TRIDMAT
END MODULE KPP_TRIDMAT_MOD

