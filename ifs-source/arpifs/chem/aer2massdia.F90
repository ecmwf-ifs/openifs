! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

 SUBROUTINE AER2MASSDIA  &
 &    (YGFL,KIDIA  , KFDIA , KLON, KLEV , KVCLIS, KTRAC, KCHEM , KAERO,  &
 &     PTSTEP, KFLDX, KFLDX2 , KLEVX, &
 &     PEXTRA ,  &
 &     PAERDDP, PAERSDM, PAERSRC)

!**   DESCRIPTION 
!     ----------
!
!   Fill aer diagnostic in extra fields for massdia chem  
!
!
!
!**   INTERFACE.
!     ----------
!          *AER2MASSDIA  IS CALLED FROM *CALLPAR*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! KLEV  :  Number of Levels
! KTRAC :  Number tracers 
! KCHEM :  Array maps chemical species into tracer array  (NCHEM <= KTRAC)
! KAERO :  Array maps chemical aerosol into tracer array  
! PTSTEP:  Time step in seconds
! PAERDDP(KLON,NACTAERO) : aerosol dry deposition 
! PAERSDM(KLON,NACTAERO)  : aerosol sedimentation 
! PAERSRC(KLON,NACTAERO)  : aerosol source term 

!
! INOUPUTS:
!--------
!
! PEXTRA(KLON,KLEVX,KFLDX)     : extra 3d field for diagnostics
!
! LOCAL:
! -------
!
! 
!     Externals.  
!     ---------                  
!      CHEM_INI_EXT              
!
!     AUTHOR.
!     -------
!        JOHANNES FLEMMING  *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2013-11-05




USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMCHEM  , ONLY : IEXTR_DD, IEXTR_EM


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON , KLEV, KVCLIS
INTEGER(KIND=JPIM),INTENT(IN) :: KTRAC, KLEVX, KFLDX, KFLDX2
INTEGER(KIND=JPIM),INTENT(IN) :: KCHEM(YGFL%NCHEM) , KAERO(YGFL%NAERO)


REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(INOUT) :: PEXTRA(KLON,KLEVX,KFLDX)
REAL(KIND=JPRB),INTENT(IN)    :: PAERDDP(KLON, YGFL%NACTAERO)  
REAL(KIND=JPRB),INTENT(IN)    :: PAERSDM(KLON,YGFL%NACTAERO)   
REAL(KIND=JPRB),INTENT(IN)    :: PAERSRC(KLON,YGFL%NACTAERO)   



!
!*       0.5   LOCAL VARIABLES
!              ---------------
INTEGER(KIND=JPIM) ::  JL, JT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('AER2MASSDIA',0,ZHOOK_HANDLE)
ASSOCIATE(NACTAERO=>YGFL%NACTAERO, NAERO=>YGFL%NAERO, NCHEM=>YGFL%NCHEM)
DO JT=1,NACTAERO 
  DO JL=KIDIA,KFDIA
    PEXTRA(JL,JT+NCHEM, IEXTR_EM) = PEXTRA(JL,JT + NCHEM, IEXTR_EM ) - PAERSRC(JL,JT) * PTSTEP   
    PEXTRA(JL,JT+NCHEM, IEXTR_DD) = PEXTRA(JL,JT + NCHEM, IEXTR_DD ) - PAERDDP(JL,JT) * PTSTEP   
    PEXTRA(JL,JT+NCHEM, IEXTR_DD) = PEXTRA(JL,JT + NCHEM, IEXTR_DD ) - PAERSDM(JL,JT) * PTSTEP   
  ENDDO
ENDDO  
   
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER2MASSDIA',1,ZHOOK_HANDLE )
END SUBROUTINE AER2MASSDIA  


