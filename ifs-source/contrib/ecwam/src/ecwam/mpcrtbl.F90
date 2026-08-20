! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPCRTBL 

! ----------------------------------------------------------------------

!**** *MPCRTBL* -

!     J. BIDLOT     ECMWF    MARCH 1996    MESSAGE PASSING


!*    PURPOSE.
!     --------
!     CREATE A TABLE THAT ASSOCIATES EACH INTEGRATED OUTPUT PARAMETER
!     WITH A SPECIFIC OUTPUT PROCESSOR. THIS TABLE IS USED TO SPREAD THE
!     RETRIEVAL OF THOSE PARAMETER FIELD IN GRID FORMAT ON AS MANY 
!     PROCESSOR AS POSSIBLE

!     IT ALSO SETS THE RELEVANT INFORMATION FOR GRIB ENCODING.

!**   INTERFACE.
!     ----------


!     METHOD.
!     -------
!     SYSTEMATICALY ASSIGN A PROCESS TO A PARAMETER AND START AGAIN WHEN
!     RUN OUT OF PROCESS.

!     EXTERNALS.
!     ----------
!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCOUT  , ONLY : JPPFLAG  ,FFLAG    ,GFLAG    ,NFLAG     ,   &
     &            LFDB     ,                                            &
     &            IPFGTBL  ,NWRTOUTWAM, COUTNAME, NIPRMOUT,ITOBOUT  ,   &
     &            NTRAIN   ,LLPARTITION,NIPRMINFO,IPRMINFO          ,   &
     &            IRWDIR, IRCD ,IRU10  , IRALTHS ,IRALTHSC ,IRALTRC ,   &
     &            IRHS     ,IRTP     ,IRT1       ,IRPHIAW  ,IRPHIOC ,   &
     &            IRTAUOC   , IRHSWS   ,IRT1WS   ,IRBATHY  ,            &
     &            IFRSTPARTI, NINFOBOUT,INFOBOUT ,COUTDESCRIPTION
!      *IPRMINFO* INTEGER    AUXILIARY INFORMATION FOR OUTPUT OF INTEGRATED PARAMETERS
!                            IPRMINFO(:,1)  : GRIB TABLE NUMBER.
!                            IPRMINFO(:,2)  : GRIB PARAMETER IDENTIFIER.
!                            IPRMINFO(:,3)  : GRIB REFERENCE LEVEL IN FULL METER.
!                            IPRMINFO(:,4)  : 1 IF SEA ICE MASK IS IMPOSED ON OUTPUT FIELD.
!                            IPRMINFO(:,5)  : 1 IF TOO SHALLOW POINTS ARE SET TO MISSING.
      USE YOWGRID  , ONLY : IJSLOC   ,IJLLOC
      USE YOWMPP   , ONLY : NPROC
      USE YOWPHYS  , ONLY : XNLEV
      USE YOWTEST  , ONLY : IU06
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpabort.intfb.h"

      INTEGER(KIND=JWIM) :: IR, IFLAG, IT, IC, ITG, ITP, ITT, IZLEV

!     1. CREATE TABLE MAPPING OUTPUT INTEGRATED PARAMETER WITH PE RANK
!        -------------------------------------------------------------

      IZLEV=NINT(XNLEV)

      COUTDESCRIPTION(:)='VARIABLE NOT DEFINED, see MPCRTBL'
      COUTNAME(:)='undef'
      IPRMINFO(:,:)=0

!     PARAMETER 001
      IRHS = DEFINE_PARAMETER( 1, 'swh', 140229, 0, .True., .True., &
                             & 'SIGNIFICANT WAVE HEIGHT' )

!     PARAMETER 002
      IR = DEFINE_PARAMETER( 2, 'mwd', 140230, 0, .True., .True., &
                           & 'MEAN WAVE DIRECTION' )

!     PARAMETER 003
      IR = DEFINE_PARAMETER( 3, 'mwp', 140232, 0, .True., .True., &
                           & 'WAVE MEAN PERIOD (-1)' )


!     PARAMETER 004
      ! Use a spare extra grib parameter number
      IR = DEFINE_PARAMETER( 4, '004', 140084, 0, .False., .True., &
                           & 'FRICTION VELOCITY' )
      IF(GFLAG(IR) ) THEN
        WRITE(IU06,*) ' ******************* NOTE ********************'
        WRITE(IU06,*) ' GRIB OUTPUT POSSIBLE FOR ', COUTDESCRIPTION(IR)
        WRITE(IU06,*) ' USING SPARE PARAMETER NUMBER ', IPRMINFO(IR,2)
        WRITE(IU06,*) ' '
      ENDIF

!     PARAMETER 005
      IRWDIR = DEFINE_PARAMETER( 5, 'dwi', 140249, IZLEV, .False., .False., &
                               & 'WAVE MODEL WIND DIRECTION' )

!     PARAMETER 006
      IRTP = DEFINE_PARAMETER( 6, 'pp1d', 140231, 0, .True., .True., &
                             & 'WAVE PEAK PERIOD' )

!     PARAMETER 007
      IRCD = DEFINE_PARAMETER( 7, 'cdww', 140233, IZLEV, .False., .False., &
                             & 'DRAG COEFFICIENT' )

!     PARAMETER 008
      ! Use a spare extra grib parameter number
      IR = DEFINE_PARAMETER( 8, '008', 140083, 0, .True., .True., &
                           & 'NORMALISED WAVE STRESS' )
      IF(GFLAG(IR) ) THEN
        WRITE(IU06,*) ' ******************* NOTE ********************'
        WRITE(IU06,*) ' GRIB OUTPUT POSSIBLE FOR ', COUTDESCRIPTION(IR)
        WRITE(IU06,*) ' USING SPARE PARAMETER NUMBER ', IPRMINFO(IR,2)
        WRITE(IU06,*) ' '
      ENDIF

!     PARAMETER 009
      IR = DEFINE_PARAMETER( 9, 'msqs', 140244, 0, .True., .True., &
                           & 'MEAN SQUARE SLOPE' )

!     PARAMETER 010
      IRU10 = DEFINE_PARAMETER( 10, 'wind', 140245, IZLEV, .False., .False., &
                              & 'WAVE MODEL WIND SPEED' )

!     PARAMETER 011
      IRHSWS = DEFINE_PARAMETER( 11, 'shww', 140234, 0, .True., .True., &
                               & 'WIND SEA WAVE HEIGHT' )

!     PARAMETER 012
      IR = DEFINE_PARAMETER( 12, 'shts', 140237, 0, .True., .True., &
                           & 'TOTAL SWELL WAVE HEIGHT' )

!     PARAMETER 013
      IR = DEFINE_PARAMETER( 13, 'mdww', 140235, 0, .True., .True., &
                           & 'WIND SEA MEAN DIRECTION' )

!     PARAMETER 014
      IR = DEFINE_PARAMETER( 14, 'mdts', 140238, 0, .True., .True., &
                           & 'TOTAL SWELL WAVE MEAN DIRECTION' )

!     PARAMETER 015
      IR = DEFINE_PARAMETER( 15, 'mpww', 140236, 0, .True., .True., &
                           & 'WIND SEA MEAN PERIOD (-1)' )

!     PARAMETER 016
      IR = DEFINE_PARAMETER( 16, 'mpts', 140239, 0, .True., .True., &
                           & 'TOTAL SWELL MEAN PERIOD (-1)' )

!     PARAMETER 017
      IRALTHS = DEFINE_PARAMETER( 17, '017', 140246, 0, .True., .True., &
                                & 'ALTIMETER WAVE HEIGHT' )

!     PARAMETER 018
      IRALTHSC = DEFINE_PARAMETER( 18, '018', 140247, 0, .True., .True., &
                                 & 'CORRECTED ALT WAVE HEIGHT' )

!     PARAMETER 019
      IRALTRC = DEFINE_PARAMETER( 19, '019', 140248, 0, .True., .True., &
                                & 'ALTIMETER RANGE CORRECTION' )

!     PARAMETER 020
      IRT1 = DEFINE_PARAMETER( 20, 'mp1', 140220, 0, .True., .True., &
                             & 'WAVE MEAN PERIOD (1)' )

!     PARAMETER 021
      IR = DEFINE_PARAMETER( 21, 'mp2', 140221, 0, .True., .True., &
                           & 'WAVE MEAN PERIOD (2)' )

!     PARAMETER 022
      IR = DEFINE_PARAMETER( 22, 'wdw', 140222, 0, .True., .True., &
                           & 'WAVE DIRECTIONAL SPREAD' )

!     PARAMETER 023
      IRT1WS = DEFINE_PARAMETER( 23, 'p1ww', 140223, 0, .True., .True., &
                               & 'WIND SEA MEAN PERIOD (1)' )

!     PARAMETER 024
      IR = DEFINE_PARAMETER( 24, 'p1ps', 140226, 0, .True., .True., &
                           & 'TOTAL SWELL MEAN PERIOD (1)' )

!     PARAMETER 025
      IR = DEFINE_PARAMETER( 25, 'p2ww', 140224, 0, .True., .True., &
                           & 'WIND SEA MEAN PERIOD (2)' )

!     PARAMETER 026
      IR = DEFINE_PARAMETER( 26, 'p2ps', 140227, 0, .True., .True., &
                           & 'TOTAL SWELL MEAN PERIOD (2)' )

!     PARAMETER 027
      IR = DEFINE_PARAMETER( 27, 'dwww', 140225, 0, .True., .True., &
                           & 'WIND SEA DIRECTIONAL SPREAD' )

!     PARAMETER 028
      IR = DEFINE_PARAMETER( 28, 'dwps', 140228, 0, .True., .True., &
                            & 'TOTAL SWELL DIRECTIONAL SPREAD' )

!     PARAMETER 029
      IR = DEFINE_PARAMETER( 29, 'wsk', 140252, 0, .True., .True., &
                           & 'WAVE SPECTRAL KURTOSIS' )

!     PARAMETER 030
      IR = DEFINE_PARAMETER( 30, 'bfi', 140253, 0, .True., .True., &
                           & 'BENJAMIN-FEIR INDEX' )

!     PARAMETER 031
      IR = DEFINE_PARAMETER( 31, 'wsp', 140254, 0, .True., .True., &
                           & 'WAVE SPECTRAL PEAKEDNESS' )

!     PARAMETER 032
      IRBATHY = DEFINE_PARAMETER( 32, 'wmb', 140219, 0, .False., .True., &
                                & 'BATHYMETRY' )

!     PARAMETER 033
      IR = DEFINE_PARAMETER( 33, 'hmax', 140218, 0, .True., .True., &
                           & 'ENVELOP MAXIMUM WAVE HEIGHT' )

!     PARAMETER 034
      IR = DEFINE_PARAMETER( 34, 'tmax', 140217, 0, .True., .True., &
                           & 'MAXIMUM WAVE PERIOD' )

!     PARAMETER 035
      IR = DEFINE_PARAMETER( 35, 'ust', 140215, 0, .True., .True., &
                           & 'U-COMP SURFACE STOKES DRIFT' )

!     PARAMETER 036
      IR = DEFINE_PARAMETER( 36, 'vst', 140216, 0, .True., .True., &
                           & 'V-COMP SURFACE STOKES DRIFT' )

!     PARAMETER 037
      IR = DEFINE_PARAMETER( 37, 'ocu', 151131, 0, .False., .True., &
                           & 'U-COMP SURFACE CURRENT' )

!     PARAMETER 038
      IR = DEFINE_PARAMETER( 38, 'vcu', 151132, 0, .False., .True., &
                           & 'V-COMP SURFACE CURRENT' )

!     PARAMETER 039
      IRPHIOC = DEFINE_PARAMETER( 39, '039', 140212, 0, .False., .True., &
                                & 'NORMALISED ENERGY FLUX TO OCEAN' )

!     PARAMETER 040
      IRPHIAW = DEFINE_PARAMETER( 40, '040', 140211, 0, .False., .True., &
                                & 'NORMALISED ENERGY FLUX TO WAVES' )

!     PARAMETER 041
      IRTAUOC = DEFINE_PARAMETER( 41, '041', 140214, 0, .False., .True., &
                                & 'NORMALISED MOMENTUM FLUX TO OCEAN' )

!     PARAMETER 042
      ITP=0 ! total number of partitions parameters
      IFRSTPARTI = DEFINE_PARAMETER( 42, '042', 140121, 0, .True., .True., &
                                   & 'SWELL PARTITION 1 WAVE HEIGHT' )
      ITP=ITP+1

!     PARAMETER 043
      IR = DEFINE_PARAMETER( 43, '043', 140122, 0, .True., .True., &
                           & 'SWELL PARTITION 1 DIRECTION' )
      ITP=ITP+1

!     PARAMETER 044
      IR = DEFINE_PARAMETER( 44, '044', 140123, 0, .True., .True., &
                           & 'SWELL PARTITION 1 MEAN PERIOD' )
      ITP=ITP+1

!     PARAMETER 045
      IR = DEFINE_PARAMETER( 45, '045', 140124, 0, .True., .True., &
                           & 'SWELL PARTITION 2 WAVE HEIGHT' )
      ITP=ITP+1

!     PARAMETER 046
      IR = DEFINE_PARAMETER( 46, '046', 140125, 0, .True., .True., &
                           & 'SWELL PARTITION 2 DIRECTION' )
      ITP=ITP+1


!     PARAMETER 047
      IR = DEFINE_PARAMETER( 47, '047', 140126, 0, .True., .True., &
                           & 'SWELL PARTITION 2 MEAN PERIOD' )
      ITP=ITP+1

!     PARAMETER 048
      IR = DEFINE_PARAMETER( 48, '048', 140127, 0, .True., .True., &
                           & 'SWELL PARTITION 3 WAVE HEIGHT' )
      ITP=ITP+1

!     PARAMETER 049
      IR = DEFINE_PARAMETER( 49, '049', 140128, 0, .True., .True., &
                           & 'SWELL PARTITION 3 DIRECTION' )
      ITP=ITP+1

!     PARAMETER 050
      IR = DEFINE_PARAMETER( 50, '050', 140129, 0, .True., .True., &
                           & 'SWELL PARTITION 3 MEAN PERIOD' )
      ITP=ITP+1

!     PARAMETER 051
      IR = DEFINE_PARAMETER( 51, '051', 140210, 0, .False., .True., &
                           & 'MEAN SQUARE STRAIN IN ICE' )

!     PARAMETER 052
      IR = DEFINE_PARAMETER( 52, '052', 140120, 0, .True., .True., &
                           & 'WAVE HEIGHT WITH PERIOD > 10s' )

!     PARAMETER 053
      IR = DEFINE_PARAMETER( 53, '053', 140209, 0, .False., .False., &
                           & 'SURFACE AIR DENSITY' )

!     PARAMETER 054
      IR = DEFINE_PARAMETER( 54, '054', 140208, 0, .False., .False., &
                           & 'CONVECTIVE VELOCITY SCALE' )

!     PARAMETER 055
      IR = DEFINE_PARAMETER( 55, 'ci', 128031, 0, .False., .True., &
                           & 'SEA ICE COVER' )

!     PARAMETER 056
      IR = DEFINE_PARAMETER( 56, '056', 003092, 0, .False., .True., &
                           & 'SEA ICE THICKNESS' )

!     PARAMETER 057
      IR = DEFINE_PARAMETER( 57, '057', 140207, 0, .True., .True., &
                           & 'SPECTRAL SKWENESS' )

!     PARAMETER 058
      IR = DEFINE_PARAMETER( 58, 'sst', 151159, 0, .False., .False., &
                           & 'NEMO SST' )

!     PARAMETER 059
      IR = DEFINE_PARAMETER( 59, 'sst', 003091, 0, .False., .False., &
                           & 'NEMO SEA ICE COVER' )

!     PARAMETER 060
      IR = DEFINE_PARAMETER( 60, '060', 003092, 0, .False., .False., &
                           & 'NEMO SEA ICE THICKNESS' )

!     PARAMETER 061
      IR = DEFINE_PARAMETER( 61, 'ucurr', 003049, 0, .False., .False., &
                           & 'NEMO ZONAL CURRENT' )

!     PARAMETER 062
      IR = DEFINE_PARAMETER( 62, 'vcurr', 003050, 0, .False., .False., &
                           & 'NEMO MERIDIONAL CURRENT' )

      IR=IR+1
!     PARAMETER 063
      IR = DEFINE_PARAMETER( 63, '063', 140112, 0, .True., .True., &
                           & 'WAVE ENERGY FLUX MAGNITUDE' )

!     PARAMETER 064
      IR = DEFINE_PARAMETER( 64, '064', 140113, 0, .True., .True., &
                           & 'WAVE ENERGY FLUX DIRECTION' )

!     PARAMETER 065
      IR = DEFINE_PARAMETER( 65, '065', 140114, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 10<=T<=12' )

!     PARAMETER 066
      IR = DEFINE_PARAMETER( 66, '066', 140115, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 12<=T<=14' )

!     PARAMETER 067
      IR = DEFINE_PARAMETER( 67, '067', 140116, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 14<=T<=17' )

!     PARAMETER 068
      IR = DEFINE_PARAMETER( 68, '068', 140117, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 17<=T<=21' )

!     PARAMETER 069
      IR = DEFINE_PARAMETER( 69, '069', 140118, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 21<=T<=25' )

!     PARAMETER 070
      IR = DEFINE_PARAMETER( 70, '070', 140119, 0, .True., .True., &
                           & 'SIG. WAVE HEIGHT 25<=T<=30' )

!     PARAMETER 071
      IR = DEFINE_PARAMETER( 71, '071', 140098, 0, .True., .True., &
                           & 'WAVE INDUCED SEA LEVEL CORRECTION' )

!     PARAMETER 072
      IR = DEFINE_PARAMETER( 72, '072', 140099, 0, .True., .True., &
                           & 'SPECTRAL WIDTH INDEX' )

!     PARAMETER 073
      IR = DEFINE_PARAMETER( 73, '073', 140100, 0, .True., .True., &
                           & 'NUMBER OF FREAK WAVES EVENT' )

!     PARAMETER 074
      IR = DEFINE_PARAMETER( 74, '074', 140101, 0, .False., .True., &
                           & 'U-COMP ATMOSPHERIC STRESS' )

!     PARAMETER 075
      IR = DEFINE_PARAMETER( 75, '075', 140102, 0, .False., .True., &
                           & 'V-COMP ATMOSPHERIC STRESS' )

!     PARAMETER 076
      IR = DEFINE_PARAMETER( 76, '076', 140103, 0, .False., .True., &
                           & 'U-COMP STRESS INTO OCEANS' )

!     PARAMETER 077
      IR = DEFINE_PARAMETER( 77, '077', 140104, 0, .False., .True., &
                           & 'V-COMP STRESS INTO OCEANS' )

!     PARAMETER 078
      IR = DEFINE_PARAMETER( 78, '078', 140105, 0, .False., .True., &
                           & 'TURB ENERGY FLUX INTO OCEANS' )

!     PARAMETER 079
      IR = DEFINE_PARAMETER( 79, 'tdcmax', 140131, 0, .True., .True., &
                           & 'TIME DOMAIN MAXIMUM CREST HEIGHT' )

!     PARAMETER 080
      IR = DEFINE_PARAMETER( 80, 'tdhmax', 140132, 0, .True., .True., &
                           & 'TIME DOMAIN MAXIMUM WAVE HEIGHT' )

!     PARAMETER 081
      IR = DEFINE_PARAMETER( 81, 'stcmax', 140133, 0, .True., .True., &
                           & 'SPACE TIME MAXIMUM CREST HEIGHT' )

!     PARAMETER 082
      IR = DEFINE_PARAMETER( 82, 'sthmax', 140134, 0, .True., .True., &
                           & 'SPACE TIME MAXIMUM WAVE HEIGHT' )


!     add new definition here:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      DO IC=1,5
        ITG=JPPFLAG-5+IC
        WRITE(COUTNAME(ITG),'(I0.3)') ITG
        ITG = DEFINE_PARAMETER( ITG, 'COUTNAME(ITG)', 140079+IC, 0, .False., .False., &
                              & 'EXTRA_FIELD '//TRIM(COUTNAME(ITG)) )
      ENDDO


!     DETERMINE THE PROCERSSORS THAT WILL CARRY OUT THE OUTPUT
      NIPRMOUT=0
      IR=1
      DO IFLAG=1,JPPFLAG
        IF(FFLAG(IFLAG).OR.GFLAG(IFLAG).OR.NFLAG(IFLAG)) THEN

          IF(FFLAG(IFLAG)) THEN
!!!!        IN CASE OF NON GRIB OUTPUT, REDIRECT TO PE 1 
            IPFGTBL(IFLAG)=1
          ELSE IF (GFLAG(IFLAG)) THEN
            IF(LFDB) THEN
              IPFGTBL(IFLAG)=IR
              IR=IR+NWRTOUTWAM
              IF(IR.GT.NPROC) IR=1
            ELSE
!!!!        IN CASE OF NO FDB REDIRECT TO PE 1 
              IPFGTBL(IFLAG)=1
            ENDIF
          ELSE
!!!       IF NFLAG but not FFLAG and not GFLAG, no output, but still need to be on the list
            IPFGTBL(IFLAG)=-1
          ENDIF

          NIPRMOUT=NIPRMOUT+1
          ITOBOUT(IFLAG)=NIPRMOUT

        ELSE
          IPFGTBL(IFLAG)=0
          ITOBOUT(IFLAG)=0
        ENDIF
      ENDDO

      IF (NIPRMOUT > 0) THEN
        IF(ALLOCATED(INFOBOUT)) DEALLOCATE(INFOBOUT)
        ALLOCATE(INFOBOUT(NIPRMOUT,NINFOBOUT))
        DO IFLAG=1,JPPFLAG
          IT=ITOBOUT(IFLAG)
          IF(IT.GT.0) THEN
            INFOBOUT(IT,1)=IPRMINFO(IFLAG,1)
            INFOBOUT(IT,2)=IPRMINFO(IFLAG,2)
            INFOBOUT(IT,3)=IPRMINFO(IFLAG,3)
          ENDIF
        ENDDO
      ENDIF

!     WILL THERE BE OUTPUT OF PARTITONED PARAMETERS
      IF(NTRAIN*(ITP/NTRAIN).NE.ITP) THEN
        WRITE(NULERR,*) '******************************************'
        WRITE(NULERR,*) '*  FATAL ERROR IN SUB. MPCRTBL           *'
        WRITE(NULERR,*) '*  THE NUMBER OF PARTITONED PARAMETERS   *'
        WRITE(NULERR,*) '*  DOES NOT MATCH !!!!                   *'
        WRITE(NULERR,*) '*  NTRAIN = ',NTRAIN
        WRITE(NULERR,*) '*  ITP = ',ITP
        WRITE(NULERR,*) '******************************************'
        WRITE(IU06,*) '******************************************'
        WRITE(IU06,*) '*  FATAL ERROR IN SUB. MPCRTBL           *'
        WRITE(IU06,*) '*  THE NUMBER OF PARTITONED PARAMETERS   *'
        WRITE(IU06,*) '*  DOES NOT MATCH !!!!                   *'
        WRITE(IU06,*) '*  NTRAIN = ',NTRAIN
        WRITE(IU06,*) '*  ITP = ',ITP
        WRITE(IU06,*) '******************************************'
        CALL ABORT1
      ENDIF
      LLPARTITION=.FALSE.
      IFLAG=IFRSTPARTI
      DO ITT=1,ITP
        IF(IPFGTBL(IFLAG).NE.0) THEN
          LLPARTITION=.TRUE.
          EXIT
        ENDIF
        IFLAG=IFLAG+1
      ENDDO


! FIND P.E. FOR OUTPUT OF RESTART FILES
!c
!cc set to 1 in order to simplify ouput to FDB
      IPFGTBL(JPPFLAG+1)=1

!!!      IFLAG=JPPFLAG
!!!     DO WHILE ((IPFGTBL(IFLAG).EQ.0).AND.(IFLAG.GT.1))
!!!         IFLAG=IFLAG-1
!!!      END DO
!!!      IPFGTBL(JPPFLAG+1)=IPFGTBL(IFLAG)+1
!!!      IF(IPFGTBL(JPPFLAG+1).GT.NPROC)IPFGTBL(JPPFLAG+1)=1


      CONTAINS

            INTEGER FUNCTION DEFINE_PARAMETER( KPARAMETER, CNAME, KGRIB_PARAMID, KGRIB_REFLEVEL, LSEA_ICE_MASK, &
                                             & LSHALLOW_TO_MISSING, CDESCRIPTION )
                  INTEGER(KIND=JWIM), INTENT(IN) :: KPARAMETER          !  PARAMETER INDEX
                  CHARACTER(LEN=*),   INTENT(IN) :: CNAME               !  GRIB PARAMETER NAME
                  CHARACTER(LEN=*),   INTENT(IN) :: CDESCRIPTION        !  PARAMETER DESCRIPTION
                  INTEGER(KIND=JWIM), INTENT(IN) :: KGRIB_PARAMID       !  GRIB PARAMETER ID (6 digits : 3 for table, 3 for index)
                  INTEGER(KIND=JWIM), INTENT(IN) :: KGRIB_REFLEVEL      !  GRIB REFERENCE LEVEL IN FULL METER
                  LOGICAL,            INTENT(IN) :: LSEA_ICE_MASK       !  TRUE IF SEA ICE MASK IS IMPOSED ON OUTPUT FIELD.
                  LOGICAL,            INTENT(IN) :: LSHALLOW_TO_MISSING !  TRUE IF TOO SHALLOW POINTS ARE SET TO MISSING

                  IF(KPARAMETER.GT.JPPFLAG) CALL MPABORT('KPARAMETER > JPPFLAG IN MPCRTBL')

                  COUTNAME(KPARAMETER)   = CNAME
                  IPRMINFO(KPARAMETER,1) = KGRIB_PARAMID / 1000
                  IPRMINFO(KPARAMETER,2) = KGRIB_PARAMID - 1000 * IPRMINFO(KPARAMETER,1)
                  IPRMINFO(KPARAMETER,3) = KGRIB_REFLEVEL
                  IF(LSEA_ICE_MASK)       IPRMINFO(KPARAMETER,4) = 1
                  IF(LSHALLOW_TO_MISSING) IPRMINFO(KPARAMETER,5) = 1

                  COUTDESCRIPTION(KPARAMETER) = CDESCRIPTION
                  DEFINE_PARAMETER = KPARAMETER
            END FUNCTION DEFINE_PARAMETER

      END SUBROUTINE MPCRTBL
