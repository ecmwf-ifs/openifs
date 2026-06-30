MODULE SURWN_MOD
CONTAINS
SUBROUTINE SURWN (KSIL,PTSTAND,PXP,PRSUN,YDDIM,YDRAD,YDLW,YDSW)

USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_DIM   , ONLY : TDIM
USE YOS_RAD   , ONLY : TRAD
USE YOS_LW    , ONLY : TLW
USE YOS_SW    , ONLY : TSW
USE ABORT_SURF_MOD
! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURWN*   - INITIALIZE COMMON YOESW

!     PURPOSE.
!     --------
!           INITIALIZE YOESW, THE COMMON THAT CONTAINS COEFFICIENTS
!           NEEDED TO RUN THE SHORTWAVE RADIATION SUBROUTINES

!**   INTERFACE.
!     ----------
!        *CALL* *SUSW

!        EXPLICIT ARGUMENTS :
!        --------------------

!      Input arguments
!    KSIL          :
!    PTSTAND       :
!    PXP           :
!    PRSUN         :

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 88-12-15
!        97-04-16 JJ Morcrette  2 and 4 interval spectral resolution
!        00-10-24 JJ Morcrette  sea-ice albedo revisited
!        00-12-14 JJ Morcrette 
!               and Ph.Dubuisson B.Bonnel 6 spectral interval resolution
!        01-04-17 Ph.Dubuisson, B.Bonnel, JJ.Morcrette 6 sp.int.resolu'n
!        01-06-28 B.Bonnel, JJ.Morcrette, Ph.Dubuisson  Rayleigh (2/4/6)
!        01-11-05 Ph.Dubuisson, JJMorcrette (new 2 intervals for TL/AD)
!         03-10-01 J.F. Estrade *ECMWF*  move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        JJMorcrette       20060530 MODIS Albedo
!        E.Dutra/G.Balsamo 20070501 add lake tile
!        JJMorcrette       20120723 fix spectral definition of snow albedo
!        R. El Khatib 27-Apr-2015 dummy init. if NWS==1
!        R. Hogan     14-Jan-2019 Replace LE4ALB with NALBEDOSCHEME
!        R. Hogan     26-Feb-2019 Removed RWEIGHT
!     ------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments
INTEGER(KIND=JPIM), INTENT(IN)    :: KSIL
REAL(KIND=JPRB)   , INTENT(IN)    :: PTSTAND 
REAL(KIND=JPRB)   , INTENT(IN)    :: PXP(6,6) 
REAL(KIND=JPRB)   , INTENT(IN)    :: PRSUN(:)
TYPE(TDIM)        , INTENT(IN)    :: YDDIM
TYPE(TRAD)        , INTENT(IN)    :: YDRAD
TYPE(TLW)         , INTENT(INOUT) :: YDLW
TYPE(TSW)         , INTENT(INOUT) :: YDSW

INTEGER(KIND=JPIM) :: IM,JM,JW,JT,I,J

REAL(KIND=JPRB) :: ZSUN2(2)
REAL(KIND=JPRB) :: ZSUN4(4),ZSUN6(6)
REAL(KIND=JPRB), ALLOCATABLE :: ZALBICE2(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZALBICE4(:,:)
REAL(KIND=JPRB), ALLOCATABLE :: ZALBICE6(:,:)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Are we using a MODIS land surface albedo climatology?
LOGICAL :: LLMODISALBEDO

IF (LHOOK) CALL DR_HOOK('SURWN_MOD:SURWN',0,ZHOOK_HANDLE)
ASSOCIATE(NMONTH=>YDDIM%NMONTH, NTILES=>YDDIM%NTILES, &
 & NSIL=>YDLW%NSIL, TSTAND=>YDLW%TSTAND, XP=>YDLW%XP, &
 & NSW=>YDRAD%NSW, NTSW=>YDRAD%NTSW)

LLMODISALBEDO = (YDRAD%NALBEDOSCHEME > 0)

! ALLOCATE 

IF(.NOT. ALLOCATED(YDSW%RALBICE_AR)) ALLOCATE(YDSW%RALBICE_AR(NMONTH,NTSW))
IF(.NOT. ALLOCATED(YDSW%RALBICE_AN)) ALLOCATE(YDSW%RALBICE_AN(NMONTH,NTSW))
ALLOCATE(ZALBICE2(NMONTH,2))
ALLOCATE(ZALBICE4(NMONTH,4))
ALLOCATE(ZALBICE6(NMONTH,6))

! LW constants

NSIL=KSIL
TSTAND=PTSTAND
DO I=1,6
  DO J=1,6
    XP(J,I)=PXP(J,I)
  ENDDO
ENDDO

!=====================================================================
!*  For sea ice, monthly values are based on Ebert and Curry, 1993, Table 2.

!*  For sea ice, monthly values are based on Ebert and Curry, 1993, Table 2.
!   We take dry snow albedo as the representative value for non-summer
!   months, and bare sea-ice as the representative value for summer
!   months. The values for Antarctic are shifted six-months.

!- 2-spectral intervals

!*  Sea ice surf. albedo for 0.25-0.69 micron (snow covered; Ebert and Curry, 1993)
ZALBICE2(1:NMONTH,1) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
 & 0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
 & 0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)  
!*  Sea ice surf. albedo for 0.69-4.00 microns (snow covered; Ebert and Curry, 1993)
!ZALBICE2(1:NMONTH,2) = (/0.664_JPRB,0.664_JPRB,0.664_JPRB,0.664_JPRB,&
!                    &0.664_JPRB,0.476_JPRB,0.288_JPRB,0.288_JPRB,&
!                    &0.664_JPRB,0.664_JPRB,0.664_JPRB,0.664_JPRB/)
ZALBICE2(1:NMONTH,2) = (/0.587_JPRB,0.587_JPRB,0.587_JPRB,0.587_JPRB,&
 & 0.587_JPRB,0.438_JPRB,0.288_JPRB,0.288_JPRB,&
 & 0.587_JPRB,0.587_JPRB,0.587_JPRB,0.587_JPRB/)  
                      
!- 4-spectral intervals

!*  Sea ice surf. albedo for 0.25-0.69 micron (snow covered; Ebert and Curry, 1993)
ZALBICE4(1:NMONTH,1) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
 & 0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
 & 0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)  
!*  Sea ice surf. albedo for 0.69-1.19 micron (snow covered; Ebert and Curry, 1993)
ZALBICE4(1:NMONTH,2) = (/0.832_JPRB,0.832_JPRB,0.832_JPRB,0.832_JPRB,&
 & 0.832_JPRB,0.638_JPRB,0.443_JPRB,0.443_JPRB,&
 & 0.832_JPRB,0.832_JPRB,0.832_JPRB,0.832_JPRB/)  
!*  Sea ice surf. albedo for 1.19-2.38 micron (snow covered; Ebert and Curry, 1993)
ZALBICE4(1:NMONTH,3) = (/0.250_JPRB,0.250_JPRB,0.250_JPRB,0.250_JPRB,&
 & 0.250_JPRB,0.153_JPRB,0.055_JPRB,0.055_JPRB,&
 & 0.250_JPRB,0.250_JPRB,0.250_JPRB,0.250_JPRB/)  
!*  Sea ice surf. albedo for 2.38-4.00 microns (snow covered; Ebert and Curry, 1993)
ZALBICE4(1:NMONTH,4) = (/0.025_JPRB,0.025_JPRB,0.025_JPRB,0.025_JPRB,&
 & 0.025_JPRB,0.030_JPRB,0.036_JPRB,0.036_JPRB,&
 & 0.025_JPRB,0.025_JPRB,0.025_JPRB,0.025_JPRB/)  
                      
!- 6-spectral intervals

!*  Sea ice surf. albedo for 0.185-0.25 micron (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,1) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
 & 0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
 & 0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)  
!*  Sea ice surf. albedo for 0.25-0.44 micron (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,2) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
 & 0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
 & 0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)  
!*  Sea ice surf. albedo for 0.44-0.69 micron (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,3) = (/0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB,&
 & 0.975_JPRB,0.876_JPRB,0.778_JPRB,0.778_JPRB,&
 & 0.975_JPRB,0.975_JPRB,0.975_JPRB,0.975_JPRB/)  
!*  Sea ice surf. albedo for 0.69-1.19 micron (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,4) = (/0.832_JPRB,0.832_JPRB,0.832_JPRB,0.832_JPRB,&
 & 0.832_JPRB,0.638_JPRB,0.443_JPRB,0.443_JPRB,&
 & 0.832_JPRB,0.832_JPRB,0.832_JPRB,0.832_JPRB/)  
!*  Sea ice surf. albedo for 1.19-2.38 micron (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,5) = (/0.250_JPRB,0.250_JPRB,0.250_JPRB,0.250_JPRB,&
 & 0.250_JPRB,0.153_JPRB,0.055_JPRB,0.055_JPRB,&
 & 0.250_JPRB,0.250_JPRB,0.250_JPRB,0.250_JPRB/)  
!*  Sea ice surf. albedo for 2.38-4.00 microns (snow covered; Ebert and Curry, 1993)
ZALBICE6(1:NMONTH,6) = (/0.025_JPRB,0.025_JPRB,0.025_JPRB,0.025_JPRB,&
 & 0.025_JPRB,0.030_JPRB,0.036_JPRB,0.036_JPRB,&
 & 0.025_JPRB,0.025_JPRB,0.025_JPRB,0.025_JPRB/)  

!     ----------------------------------------------------------------

IF(.NOT.ALLOCATED(YDSW%RSUN)) ALLOCATE(YDSW%RSUN(NTSW))

DO JW=1,NSW
  YDSW%RSUN(JW)=PRSUN(JW)
ENDDO

IF (NSW == 2) THEN

  DO JW=1,NSW
    DO JM=1,NMONTH
      IM=MOD(JM+5,NMONTH)+1
      YDSW%RALBICE_AR(IM,JW)=ZALBICE2(IM,JW)
      YDSW%RALBICE_AN(JM,JW)=YDSW%RALBICE_AR(IM,JW)
    ENDDO
  ENDDO

ELSEIF (NSW == 4) THEN

  DO JW=1,NSW
    DO JM=1,NMONTH
      IM=MOD(JM+5,NMONTH)+1
      YDSW%RALBICE_AR(IM,JW)=ZALBICE4(IM,JW)
      YDSW%RALBICE_AN(JM,JW)=YDSW%RALBICE_AR(IM,JW)
    ENDDO
  ENDDO

ELSEIF (NSW == 6) THEN

  DO JW=1,NSW
    DO JM=1,NMONTH
      IM=MOD(JM+5,NMONTH)+1
      YDSW%RALBICE_AR(IM,JW)=ZALBICE6(IM,JW)
      YDSW%RALBICE_AN(JM,JW)=YDSW%RALBICE_AR(IM,JW)
    ENDDO
  ENDDO

ELSEIF (NSW == 1) THEN

! Dummy initialization for a single SW interval
  YDSW%RALBICE_AR(:,:)=HUGE(1._JPRB)
  YDSW%RALBICE_AN(:,:)=HUGE(1._JPRB)

ELSE

  CALL ABORT_SURF('SURWN: WRONG NUMBER OF SW INTERVALS')

ENDIF

DEALLOCATE(ZALBICE2)
DEALLOCATE(ZALBICE4)
DEALLOCATE(ZALBICE6)
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SURWN_MOD:SURWN',1,ZHOOK_HANDLE)

END SUBROUTINE SURWN
END MODULE SURWN_MOD
