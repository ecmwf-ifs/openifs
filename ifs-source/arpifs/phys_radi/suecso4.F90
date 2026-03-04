! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

#ifdef RS6K
@PROCESS NOOPTIMIZE
#endif
!pgi$r opt=0 
SUBROUTINE SUECSO4 ( YD_RAERSO4, KCMIPFIXYR, KINDAT, KMINUT )

!**** *SUECSO4* - DEFINES HISTORICAL CLIMATOLOGICAL DISTRIBUTION OF AEROSOLS

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!      J.-J. MORCRETTE  E.C.M.W.F.    98/12/21

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. Elkhatib 12-10-2005 Split for faster and more robust compilation.
!      HHersbach 01-04-2011 Replace code by input from CMIP5-recommended CAM3.5 data set
!      F. Vana  05-Mar-2015  Support for single precision
!      C. Roberts 2017-05-30 Added NCMIPFIXYR.
!-----------------------------------------------------------------------

  USE PARKIND1           , ONLY : JPRD, JPIM     ,JPRB
  USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE YOEAERC            , ONLY : NLDECSO4DIM, NFDECSO4DIM, CLISTSO4, FILESO4
  USE YOMLUN             , ONLY : NULOUT, FOPEN
  USE REGLATLON_FIELD_MIX, ONLY : REGLATLON_FIELD, CREATE_REGLATLON_FIELD, STATS_REGLATLON
  USE EC_DATETIME_MOD,   ONLY : HOURDIFF, MININCR

  IMPLICIT NONE

  TYPE(REGLATLON_FIELD), INTENT(INOUT) :: YD_RAERSO4
  INTEGER(KIND=JPIM)   , INTENT(IN)    :: KCMIPFIXYR
  INTEGER(KIND=JPIM)   , INTENT(IN)    :: KINDAT  ! CCYYMMDD of start of forecast
  INTEGER(KIND=JPIM)   , INTENT(IN)    :: KMINUT  ! forecast step in minuts

!     -----------------------------------------------------------------

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!! Following values with SAVE attribute mean that all simultaneous OOPS MODEL objects will use same values
!! Will there be problems when assimilation windows cross a decade divide?
  INTEGER(KIND=JPIM)         ,SAVE :: INLATSO4,INLONSO4
  REAL(KIND=JPRB),ALLOCATABLE,SAVE :: ZSO4(:,:,:,:)
  INTEGER(KIND=JPIM)               :: IMONTH,ILAT,IRET,I,J
  INTEGER(KIND=JPIM)               :: IYR,IMM,IDD,IHH,IMI
  INTEGER(KIND=JPIM)               :: IYR1,IYR2,IYR5,IMON(2),IH21,IH2
  INTEGER(KIND=JPIM)               :: IDEC1,IDEC2
  REAL(KIND=JPRB)                  :: ZW,ZWDEC(2),ZWMON(2)

  REAL(KIND=JPRD),ALLOCATABLE      :: ZSO4_D(:)
  INTEGER(KIND=JPIM)         ,SAVE :: IDECH(2),INFDECSO4,INLDECSO4
  REAL(KIND=JPRB)            ,SAVE :: ZLAT1,ZLATN,ZLON1
  LOGICAL                    ,SAVE :: LLEW

#include "abor1.intfb.h"
#include "fcttim.func.h"


!     ------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('SUECSO4',0,ZHOOK_HANDLE)

!*         0.     INITIALIZATION (ONLY ONCE)
!                 --------------------------

  !! on first call to routine, YD_RAERSO4 has not yet been initialized
  IF (YD_RAERSO4%NLAT==-1) THEN
    CALL LINK_DATAFILES
  END IF

!*         1.0    FIND ACTUAL DATE AND TIME (UP TO THE MINUTE)
!                 --------------------------------------------

  CALL MININCR (NCCAA(KINDAT),NMM(KINDAT),NDD(KINDAT),0  ,0  , KMINUT, &
       &        IYR          ,IMM        ,IDD        ,IHH,IMI, IRET   )


!*         2.0    FIND INTERPOLATION DECADES AND WEIGHTS
!                 --------------------------------------
IF (KCMIPFIXYR>0) IYR=KCMIPFIXYR ! If using perpetual CMIP forcing

!-Assume that data is valid for the year 05 within a decade

  IYR5 =IYR-5
  IDEC1=IYR5/10 ; ZWDEC(2)=0.1*IYR5-IDEC1
  IDEC2=IDEC1+1 ; ZWDEC(1)=1.0_JPRB-ZWDEC(2)

  IDEC1=MAX(INFDECSO4,MIN(IDEC1,INLDECSO4))
  IDEC2=MAX(INFDECSO4,MIN(IDEC2,INLDECSO4))

!-If necessary, update input data
  IF(IDEC1/=IDECH(1) .OR. IDEC2/=IDECH(2)) THEN
     IDECH(1)=IDEC1 ; CALL READ_SO4_BURDEN(IDEC1,1)
     IDECH(2)=IDEC2 ; CALL READ_SO4_BURDEN(IDEC2,2)
     WRITE(NULOUT,'(A,2I5)')"SUECSO4: UPDATE HISTORY FILES FOR DECADES: ",10*IDECH
  ELSE
     WRITE(NULOUT,'(A,2I5)')"SUECSO4: USE AVAILABLE FILES FOR DECADES:  ",10*IDECH
  ENDIF

!-Allocate
  IF (YD_RAERSO4%NLAT==-1) THEN
    CALL CREATE_REGLATLON_FIELD(YD_RAERSO4,KNLAT=INLATSO4,KNLON=INLONSO4,&
         & PLAT1=ZLAT1,PLATN=ZLATN,PLON1=ZLON1,LDEW=LLEW,CDNAME="RAERSO4")
  ENDIF


!*         2.0    FIND INTERPOLATION MONTHS AND WEIGHTS
!                 -------------------------------------
!-Assume that data is valid for day 15 00UTC within a month
 
  IYR1   =IYR
  IMON(1)=IMM
  IF(IDD<15) IMON(1)=IMON(1)-1
  IF (IMON(1)==0) THEN
     IYR1   =IYR1-1
     IMON(1)=12 
  ENDIF

  IYR2=IYR1
  IMON(2)=IMON(1)+1
  IF (IMON(2)==13) THEN
     IYR2   =IYR2+1
     IMON(2)=1 
  ENDIF

!-Find interpolation weights
  CALL HOURDIFF(IYR2,IMON(2),15,00, IYR1,IMON(1),15  ,00,IH21,IRET)
  CALL HOURDIFF(IYR2,IMON(2),15,00, IYR ,IMM    ,IDD,IHH,IH2 ,IRET)
  ZWMON(1)=REAL(IH2)/REAL(IH21)
  ZWMON(2)=1.0_JPRB-ZWMON(1)

!-Write out weigths
! write(NULOUT,'(2I10,2X,I4,2(2I2.2,X),2(2I5,2F8.4)')KINDAT,KMINUT,IYR,IMM,IDD,IHH,IMI, &
!     &  10*IDECH,ZWDEC,IMON,ZWMON

  

!*         3.0    INTERPOLATE AND TRANSFORM FROM KG/M^2 TO G/M^2
!                 ----------------------------------------------

  YD_RAERSO4%PFLD(:,:)=0.0_JPRB
  DO I=1,2 
     DO J=1,2
        ZW=(1.E+3_JPRB)*ZWDEC(I)*ZWMON(J)
        YD_RAERSO4%PFLD(:,:)=YD_RAERSO4%PFLD(:,:)+ZW*ZSO4(:,:,IMON(J),I)
     ENDDO
  ENDDO
  CALL STATS_REGLATLON("SUBROUTINE SUECSO4:",YD_RAERSO4)

!     ------------------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('SUECSO4',1,ZHOOK_HANDLE)

!-----------------------------------------------------------------------------

CONTAINS

SUBROUTINE LINK_DATAFILES
  IMPLICIT NONE

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM)  :: IU, ILIM, KDEC, ISPACE
  CHARACTER(LEN=1000) :: CLINE

  IF (LHOOK) CALL DR_HOOK('SUECSO4:LINK_DATAFILES',0,ZHOOK_HANDLE)

  CALL FOPEN(KUNIT=IU,CDFILE=TRIM(CLISTSO4),CDFORM='FORMATTED',CDTEXT="SUECSO4: LIST OF SO4 HISTORY FILES")

!-Some indicators that things have not yet been set up properly
  IDECH(:)=0

  INFDECSO4=NLDECSO4DIM
  INLDECSO4=NFDECSO4DIM
  FILESO4(:)=""

  DO ILIM=1,10000
     READ(IU,'(A1000)',END=999)CLINE
     READ(CLINE,*)KDEC 
     KDEC=KDEC/10

     IF (KDEC<NFDECSO4DIM .OR. KDEC>NLDECSO4DIM)CYCLE
     INFDECSO4=MIN(INFDECSO4,KDEC)
     INLDECSO4=MAX(INLDECSO4,KDEC)

     ISPACE=SCAN(CLINE," ")
     IF (ISPACE<1) &
        &CALL ABOR1("LINK_DATAFILES: INPUT LIST HAS INCORRECT FORMAT")

     FILESO4(KDEC)=TRIM(ADJUSTL(CLINE(ISPACE:1000)))    

  ENDDO
  CALL ABOR1("LINK_DATAFILES: TOO MANY INPUT FILES")
 999 CONTINUE
  CLOSE(IU)

! SOME QC
  WRITE(NULOUT,'(2(A,I4))')"SUECSO4: SETUP SO4 HISTORY FILES FOR DECADES: ",10*INFDECSO4," - ",10*INLDECSO4
  DO KDEC=INFDECSO4,INLDECSO4
     IF (FILESO4(KDEC)=="") THEN
        CALL ABOR1("SUECSO4: SOME INTERMEDIATE DECADES ARE UNDEFINED: ")
     ELSE
        WRITE(NULOUT,'(I4,": ",A)')10*KDEC,TRIM(FILESO4(KDEC))
     ENDIF
  ENDDO

  IF (LHOOK) CALL DR_HOOK('SUECSO4:LINK_DATAFILES',1,ZHOOK_HANDLE)

END SUBROUTINE LINK_DATAFILES

!-----------------------------------------------------------------------------

SUBROUTINE READ_SO4_BURDEN(KDEC,KND)

  IMPLICIT NONE

  INTEGER(KIND=JPIM),INTENT(IN ) :: KDEC,KND

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  INTEGER(KIND=JPIM)             :: INLAT,INLON,ISPACE
  INTEGER(KIND=JPIM),PARAMETER   :: INC=10000
  CHARACTER(LEN=10000)           :: CLINE

  INTEGER(KIND=JPIM)             :: IU,ISKIP
  REAL(KIND=JPRD),ALLOCATABLE    :: ZDUM(:)

  IF (LHOOK) CALL DR_HOOK('SUECSO4:READ_SO4_BURDEN',0,ZHOOK_HANDLE)

  IF(KND<1 .OR. KND>2) CALL ABOR1("READ_SO4_BURDEN: KND SHOULD BE 1 OR 2")
  CALL FOPEN(KUNIT=IU,CDFILE=TRIM(FILESO4(KDEC)),CDFORM='FORMATTED',CDTEXT="READ_SO4_BURDEN")

!-Latitude
  READ(IU,'(A10000)')CLINE
  ISPACE=1+SCAN(CLINE," ")
  READ(CLINE(ISPACE:INC),*)INLAT
  IF (INLATSO4==0) THEN
      INLATSO4=INLAT
      ALLOCATE(ZDUM(INLAT))
      READ(CLINE(ISPACE:INC),*)INLAT,ZDUM
      ZLAT1=REAL(ZDUM(1),JPRB)
      ZLATN=REAL(ZDUM(INLAT),JPRB)
      DEALLOCATE(ZDUM)
  ELSEIF (INLATSO4/=INLAT) THEN
      CALL ABOR1("INCONSISTENCY IN # OF LATITUDES")
  ENDIF 

!-Longitude
  READ(IU,'(A10000)')CLINE
  ISPACE=1+SCAN(CLINE," ")
  READ(CLINE(ISPACE:INC),*)INLON
  IF (INLONSO4==0) THEN
      INLONSO4=INLON
      ALLOCATE(ZDUM(INLON))
      READ(CLINE(ISPACE:INC),*)INLON,ZDUM
      ZLON1=REAL(ZDUM(1),JPRB)
      LLEW=.FALSE. 
      IF (REAL(ZDUM(INLON),JPRB)<ZLON1) LLEW=.TRUE.
      DEALLOCATE(ZDUM)
  ELSEIF (INLONSO4/=INLON) THEN
      CALL ABOR1("INCONSISTENT # OF LONGITUDES")
  ENDIF

!-Data
  IF (.NOT.ALLOCATED(ZSO4)) ALLOCATE(ZSO4(INLAT,INLON,12,2))
  ISKIP=2 ! The first two lines state year and press level
  ALLOCATE(ZSO4_D(INLON))
  DO IMONTH=1,12
     DO ILAT=1-ISKIP,INLAT
        READ(IU,'(A10000)')CLINE
        IF (ILAT>0) THEN
           READ(CLINE,*) ZSO4_D(1:INLON)
           ZSO4(ILAT,1:INLON,IMONTH,KND) = REAL(ZSO4_D(1:INLON),JPRB)
!       ELSE
!          WRITE(NULOUT,'(A)')TRIM(CLINE)
        ENDIF
     ENDDO
  ENDDO
  CLOSE(IU)
  DEALLOCATE(ZSO4_D)
  WRITE(NULOUT,'(A)')"READ_SO4_BURDEN: CLOSE DATA FILE"
  WRITE(NULOUT,*)"READ_SO4_BURDEN: AVERAGE ZSO4: ",&
  &  KND,SUM(ZSO4(:,:,:,KND))/(INLAT*INLON*12)

  IF (LHOOK) CALL DR_HOOK('SUECSO4:READ_SO4_BURDEN',1,ZHOOK_HANDLE)

END SUBROUTINE READ_SO4_BURDEN

!-----------------------------------------------------------------------------
END SUBROUTINE SUECSO4
