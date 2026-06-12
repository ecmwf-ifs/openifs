MODULE CCETR_MOD
CONTAINS
SUBROUTINE CCETR(KIDIA,KFDIA,KLON,KVTYPE,LDLAND,PIA,PMU0,PABC,PLAI,YDAGS,PXIA)

! (C) Copyright 1998- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**  *CCETR* 

!     PURPOSE
!     -------
!     Calculates radiative transfer within the canopy

!     PARAMETER     DESCRIPTION                                   UNITS
!     ---------     -----------                                   -----
!     INPUT PARAMETERS (REAL):

!     *PIA*          ABSORBED PAR AT THE TOP OF THE CANOPY        W M-2   
!     *PMU0*        LOCAL COSINE OF INSTANTANEOUS MEAN SOLAR ZENITH ANGLE -
!     *PABC*         ABSCISSA NEEDED FOR INTEGRATION              -
!     *PLAI*         LEAF AREA INDEX                              M2/M2

!     OUTPUT PARAMETERS (REAL):
!     *PXIA*         INCIDENT RADIATION AFTER DIFFUSION           W M-2   


!     METHOD
!     ------
!     Calvet et al. 1998 Forr. Agri. Met. 
!     [from model of Jacobs(1994) and Roujean(1996)]
!
!     EXTERNALS
!     --------
!     none

!     REFERENCE
!     ---------
!     Calvet et al. 1998 Forr. Agri. Met. 
!      
!     AUTHOR
!     ------
!	  A. Boone           * Meteo-France *
!      (following Belair)

!     MODIFICATIONS
!     -------------
!      Original    27/10/97 
!      M.H. Voogt (KNMI) "C-Tessel"  09/2005
!      S. Lafont (ECMWF) externalised CTESSEL 04/2006     
!      G. Balsamo (ECMWF) 24/3/2014 cleaning and LDLAND protection
!     ---------------------------------------------------------------------
!
USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_AGS , ONLY : TAGS

IMPLICIT NONE
INTEGER(KIND=JPIM), INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM), INTENT(IN)    :: KLON
INTEGER(KIND=JPIM), INTENT(IN)    :: KVTYPE(:)
LOGICAL,            INTENT(IN)    :: LDLAND(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PIA(:)
REAL(KIND=JPRB),    INTENT(IN)    :: PMU0(:) 
REAL(KIND=JPRB),    INTENT(IN)    :: PABC
REAL(KIND=JPRB),    INTENT(IN)    :: PLAI(:)
TYPE(TAGS),         INTENT(IN)    :: YDAGS
REAL(KIND=JPRB),    INTENT(OUT)   :: PXIA(:)

!*         0.     LOCAL VARIABLES.
!                 ----- ----------

REAL(KIND=JPRB)  :: ZXFD(KLON),ZXSLAI(KLON),ZXIDF(KLON),ZXIDR(KLON)
INTEGER(KIND=JPIM) :: JL
! ZXFD   = fraction of diffusion
! ZXSLAI = LAI of upper layer 
! ZXIDF  = interception of diffusion
! ZXIDR  = direct interception
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!     ---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CCETR_MOD:CCETR',0,ZHOOK_HANDLE)
ASSOCIATE(RDIFRACF=>YDAGS%RDIFRACF, RXBOMEGA=>YDAGS%RXBOMEGA, &
 & RXBOMEGAM=>YDAGS%RXBOMEGAM, RXGT=>YDAGS%RXGT)

! initialization of local variable
DO JL=KIDIA,KFDIA
  ZXFD(JL)=0._JPRB
  ZXSLAI(JL)=0._JPRB
  ZXIDF(JL)=0._JPRB
  ZXIDR(JL)=0._JPRB
END DO


DO JL=KIDIA,KFDIA
  PXIA(JL)=0._JPRB
  IF (LDLAND(JL) .AND. PIA(JL) .GT. 0._JPRB) THEN

! fraction of diffusion
    ZXFD(JL)= RDIFRACF/(RDIFRACF+PMU0(JL))                 

! LAI of upper layer
    ZXSLAI(JL)=PLAI(JL)*(1.0_JPRB-PABC)                             

! interception of diffusion
!   ZXIDF(JL)=ZXFD(JL)*(1.0_JPRB-EXP(-0.8_JPRB*ZXSLAI(JL)*RXBOMEGA))   
    ZXIDF(JL)=ZXFD(JL)*(1.0_JPRB-EXP(-0.8_JPRB*ZXSLAI(JL)*RXBOMEGAM(KVTYPE(JL))))
! direct interception
!   ZXIDR(JL)=(1.0_JPRB-ZXFD(JL))*(1.0_JPRB-EXP(-RXGT*ZXSLAI(JL)*RXBOMEGA/PMU0(JL)))   
    ZXIDR(JL)=(1.0_JPRB-ZXFD(JL))*(1.0_JPRB-EXP(-RXGT*ZXSLAI(JL)*RXBOMEGAM(KVTYPE(JL))/PMU0(JL)))
! Adjusted radiation:
    PXIA(JL)=PIA(JL)*(1.0_JPRB-ZXIDF(JL)-ZXIDR(JL)) 
  ENDIF
END DO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('CCETR_MOD:CCETR',1,ZHOOK_HANDLE)
END SUBROUTINE CCETR
END MODULE  CCETR_MOD
 
