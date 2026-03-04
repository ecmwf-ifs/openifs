! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE TM5_STRATBC_CH4(YGFL,KIDIA,KFDIA,KLON,PTSTEP,KMONTH,PGELAT,KMODE,PCH4,PTENCH4)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!
!*** Stratospheric boundary conditions for CH4 as a function of latitude and month
!*** Based on HALOE observations
!
!--------------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *TM5_stratbc_ch4  IS CALLED FROM *CHEM_tm5*.

! INPUTS:
! -------
! KIDIA :  Start of Array  
! KFDIA :  End  of Array 
! KLON  :  Length of Arrays 
! PTSTEP:  Time step in seconds 
! KMODE : Selection of altitude level 
!
! PCH4(KLON)        : initial mass ratios OF CH4 TRACER  at two pressure levels    (kg/kg)
!
!
! OUTPUTS:
! -------
! PTENCH4 (KLON)   : tendency  change (kg/kg/sec) at two pressure levels
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        TM5-community    
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2013-04-23



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI, RMD   ! dry air molar mass
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMLUN   , ONLY : NULERR
USE TM5_CHEM_MODULE, ONLY : ICH4,HALOE_CH4

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD)   ,INTENT(INOUT):: YGFL
INTEGER(KIND=JPIM),INTENT(IN) :: KIDIA , KFDIA , KLON,KMONTH,KMODE
REAL(KIND=JPRB),INTENT(IN)    :: PTSTEP
REAL(KIND=JPRB),INTENT(IN)  :: PGELAT(KLON)   ! Latitude ( radians)
REAL(KIND=JPRB),INTENT(IN)  :: PCH4(KLON)    ! initial concentrations
REAL(KIND=JPRB),INTENT(OUT) :: PTENCH4(KLON) ! tendencies

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
REAL(KIND=JPRB)    :: ZLAT

!* nudging coefs...
REAL(KIND=JPRB),PARAMETER       ::  GNUD_HILAT = 2.9E-6   !  1 / (4   days)
REAL(KIND=JPRB),PARAMETER       ::  GNUD_TROP  = 4.8E-6   !  1 / (2.5 days)
REAL(KIND=JPRB)                 ::  ZNUD, ZCH4_CLIM,ZCH4_NEW

! * counters
INTEGER(KIND=JPIM) :: JL,J_LAT

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('TM5_STRATBC_CH4',0,ZHOOK_HANDLE )
ASSOCIATE(YCHEM=>YGFL%YCHEM)
IF (KMODE <= 2) THEN 
  DO JL=KIDIA,KFDIA
    ZLAT=(180.0_JPRB/RPI)*PGELAT(JL)
    J_LAT = NINT(ZLAT + 91.0_JPRB )
    J_LAT = MAX(1_JPIM,MIN(J_LAT,180_JPIM))
    ZNUD=GNUD_HILAT
    IF (ABS(ZLAT) < 30._JPRB) ZNUD=GNUD_TROP
!*  climatological and updated concentrations [kg/kg]
    ZCH4_CLIM = HALOE_CH4(KMONTH,J_LAT,KMODE)*YCHEM(ICH4)%RMOLMASS/RMD
    ZCH4_NEW  = (PCH4(JL) + PTSTEP*ZNUD*ZCH4_CLIM)/(1._JPRB+PTSTEP*ZNUD)
!*  tendency (kg/kg/sec)
    PTENCH4(JL) = (ZCH4_NEW - PCH4(JL))  / PTSTEP
  ENDDO
ELSE
  WRITE(NULERR,*)'This altitude for constraining strat. CH4 is not available: ',KMODE
ENDIF


END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('TM5_STRATBC_CH4',1,ZHOOK_HANDLE )

END SUBROUTINE TM5_STRATBC_CH4
