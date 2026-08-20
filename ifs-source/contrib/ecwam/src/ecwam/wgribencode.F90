! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE WGRIBENCODE ( IU06, ITEST, &
&                        I1, I2, &
&                        FIELD, &
&                        ITABLE, IPARAM, &
&                        KLEV, &
&                        IK, IM, &
&                        CDATE, IFCST, MARSTYPE, &
&                        PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, &
&                        LGRHDIFS, &
&                        NDATE_TIME_WINDOW_END, &
&                        NGRBRESS, LNEWLVTP, LPADPOLES, &
&                        NLONRGG_SIZE, NLONRGG, IRGG, &
&                        AMONOP, AMOSOP, XDELLA, CLDOMAIN, &
&                        KCOUSTEP, LRSTST0, &
&                        ZMISS, &
&                        IGRIB_HANDLE)

! ----------------------------------------------------------------------

!****  *WGRIBENCODE*  ENCODES WAM MODEL FIELD INTO GRIB CODE AND OUTPUT

!       PURPOSE.
!       --------
!         SUBROUTINE PACKS WAVE FIELDS INTO THE GRIB CODE

!**    INTERFACE.
!      ----------

!          *IU06*       LOGFILE OUTPUT UNIT.
!          *ITEST*      TEST OUTPUT GIVEN IF ITEST GT 2.
!          *I1*         FIRST DIMENSION OF FIELD.
!          *I2*         SECOND DIMENSION OF FIELD.
!          *FIELD*      FIELD TO BE PACKED.
!          *ITABLE*     GRIB TABLE NUMBER.
!          *IPARAM*     PARAMETER IDENTIFIER.
!          *KLEV*       REFERENCE LEVEL IN full METER
!                       (SHOULD BE 0 EXCEPT FOR 233 AND 245)
!          *IK*         DIRECTION INDEX,
!                       ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *IM*         FREQUENCY INDEX,
!                       ONLY MEANINGFUL FOR SPECTRAL PARAMETERS.
!          *CDATE*      ACTUAL DATE AND TIME. IF NOT A FORECAST (IFCST<=0), OTHERWISE
!                       DATE AND TIME OF THE START OF THE FOERCAST.
!          *IFCST*      FORECAST STEP IN HOURS.
!          *MARSTYPE*   TYPE OF CURRENT FIELD
!          **           VARIOUS PARAMETERS, SEE BELOW
!          *IGRIB_HANDLE  GRIB HANDLE TO BE POPULATED WITH THE DATA, THE HANDLE SHOULD BE
!                         CLONED BEFORE SENDING TO THIS SUBROUTINE

!      METHOD.
!      -------

!      EXTERNALS.
!      ----------

!      REFERENCES.
!      -----------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWGRIB  , ONLY : IGRIB_GET_VALUE, IGRIB_SET_VALUE
      USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
      USE EC_LUN   , ONLY : NULERR

! ----------------------------------------------------------------------
      IMPLICIT NONE
#include "abort1.intfb.h"
#include "difdate.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ITEST, I1, I2
      INTEGER(KIND=JWIM), INTENT(IN) :: ITABLE, IPARAM, KLEV, IK, IM, IFCST
      INTEGER(KIND=JWIM), INTENT(INOUT) :: IGRIB_HANDLE

      REAL(KIND=JWRB), INTENT(INOUT) :: FIELD(I1,I2)

      CHARACTER(LEN=2), INTENT(IN) :: MARSTYPE
      CHARACTER(LEN=14), INTENT(IN) :: CDATE

      ! From yowgribhd
      LOGICAL, INTENT(IN)   :: LGRHDIFS         ! If true then grib header will use information as provided by the ifs.
      REAL(KIND=JWRB), INTENT(IN)      :: PPMISS  ! every spectral values less or equal ppmiss are replaced by the missing data indicator
      REAL(KIND=JWRB), INTENT(IN)      :: PPEPS   ! Small number used in spectral packing of 251
      REAL(KIND=JWRB), INTENT(IN)      :: PPREC   ! Reference value for spectral packing of 251
      REAL(KIND=JWRB), INTENT(IN)      :: PPRESOL ! Maximun resolution possible when encoding spectra (parameter 251).
      REAL(KIND=JWRB), INTENT(IN)      :: PPMIN_RESET      ! Can be used to set the minimum of ppmin in wgribout to a lower value.
      INTEGER(KIND=JWIM), INTENT(IN)   :: NTENCODE         ! Total number of grid points for encoding
      INTEGER(KIND=JWIM), INTENT(IN)   :: NDATE_TIME_WINDOW_END
      INTEGER(KIND=JWIM), INTENT(IN)   :: NGRBRESS         ! Number of bits used to encode spectra
      LOGICAL, INTENT(IN)   :: LNEWLVTP         ! If true the new levtype definition will be used.
      LOGICAL, INTENT(IN)   :: LPADPOLES        ! True if poles are padded when savind to grib.

      ! From yowgrid
      INTEGER(KIND=JWIM), INTENT(IN)   :: NLONRGG_SIZE 
      INTEGER(KIND=JWIM), INTENT(IN)   :: NLONRGG(NLONRGG_SIZE) 

      ! From yowmap
      INTEGER(KIND=JWIM), INTENT(IN) :: IRGG                ! Grid code: 0 = regular, 1 = irregular.
      REAL(KIND=JWRB),    INTENT(IN) :: AMONOP              ! Most Northern latitude in grid (degree).
      REAL(KIND=JWRB),    INTENT(IN) :: AMOSOP              ! Most Southern latitude in grid (degree).
      REAL(KIND=JWRB),    INTENT(IN) :: XDELLA              ! Grid increment for latitude (degree).

      ! From yowparam
      CHARACTER(LEN=1), INTENT(IN) :: CLDOMAIN   ! Defines the domain of the model (for the fdb and for selection of some variables)

      ! From yowcoup
      INTEGER(KIND=JWIM), INTENT(IN)     :: KCOUSTEP        ! Coupling time to the IFS (in seconds).

      ! From yowcout
      LOGICAL, INTENT(IN)     :: LRSTST0      ! True if GRIB header have to be reset such that the forecast step points to the start of the run.

      ! From yowpcons
      REAL(KIND=JWRB), INTENT(IN)        :: ZMISS           ! Missing data indicator (set in chief or via the ifs).


      INTEGER(KIND=JWIM) :: ICLASS, ISTEP, ISTEP_HRS 
      INTEGER(KIND=JWIM) :: IC, ITABPAR, IDATE, ITIME, IGRIB_VERSION, ILEVTYPE
      INTEGER(KIND=JWIM) :: ICOUNT, NN, I, J, JSN, KK, MM
      INTEGER(KIND=JWIM) :: IY1,IM1,ID1,IH1,IMN1,ISS1,IDATERES
      INTEGER(KIND=JWIM) :: IY2,IM2,ID2,IH2,IMN2,ISS2
      INTEGER(KIND=JWIM) :: IERR
      INTEGER(KIND=JWIM) :: NWINOFF
      INTEGER(KIND=JWIM) :: NPROMA, MTHREADS, JC, JCS, JCL, JJ, ITHRS
!$    INTEGER,EXTERNAL :: OMP_GET_MAX_THREADS

      REAL(KIND=JWRB) :: TEMP
      REAL(KIND=JWRB) :: ZMINSPEC, PPMAX, PPMIN, DELTAPP, ABSPPREC
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB), ALLOCATABLE :: VALUES(:)
      REAL(KIND=JWRB), ALLOCATABLE :: VALM(:)

      CHARACTER(LEN=12) :: C12
      CHARACTER(LEN=14) :: CDATE1, CDATE2

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE',0,ZHOOK_HANDLE)

      !CALL GSTATS(1709,0)

      IF (ITEST > 0) THEN
        WRITE(IU06,*) '   SUB. WGRIBENCODE CALLED FOR ',IPARAM
        CALL FLUSH(IU06)
      ENDIF

      ALLOCATE(VALUES(NTENCODE))

      MTHREADS=1
!$    MTHREADS=OMP_GET_MAX_THREADS()
      NPROMA=NTENCODE/MTHREADS + 1


!*    0. PUT FIELD INTO GLOBAL MATRIX VALUES.
!        -----------------------------------

      IF (IRGG == 1 .OR. CLDOMAIN == 'm' ) THEN
        ICOUNT=1
      ELSEIF (CLDOMAIN == 's' ) THEN
        ICOUNT=1
      ELSE
        ICOUNT = (NINT((90._JWRB - AMONOP ) / XDELLA))*I1 + 1
        VALUES(:)=ZMISS
      ENDIF

!     PAD THE POLES IF INCLUDED IN THE GRID
      IF (LPADPOLES) THEN
        IF ((NINT((90._JWRB - AMONOP ) / XDELLA)) == 0) THEN
          TEMP=0._JWRB
          NN=0
          J=2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF (FIELD(I,J) /= ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF (NN > 0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
        IF ((NINT((-90._JWRB - AMOSOP ) / XDELLA)) == 0) THEN
          TEMP=0._JWRB
          NN=0
          J=I2-1
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            IF (FIELD(I,J) /= ZMISS) THEN
              TEMP=TEMP+FIELD(I,J)
              NN=NN+1
            ENDIF
          ENDDO
          IF (NN > 0) THEN
            TEMP=TEMP/NN
          ELSE
            TEMP=ZMISS
          ENDIF
          J=I2
          JSN=I2-J+1
          DO I=1,NLONRGG(JSN)
            FIELD(I,J)=TEMP
          ENDDO
        ENDIF
      ENDIF

!     FILL THE ENCODING ARRAY
!     -----------------------
      DO J = 1, I2
        JSN = I2-J+1
        VALUES(ICOUNT:ICOUNT+NLONRGG(JSN)-1)=FIELD(1:NLONRGG(JSN), J)
        ICOUNT=ICOUNT+NLONRGG(JSN)
      ENDDO

      IF (IPARAM == 251) THEN

!       spectra are encoded on a log10 scale and small below a certain threshold are set to missing
        DELTAPP=(2**NGRBRESS-1)*PPRESOL
        ABSPPREC=ABS(PPREC)
        ZMINSPEC = LOG10(PPEPS)+ABSPPREC

        ALLOCATE(VALM(MTHREADS))
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC, JCS, JCL, JJ, ITHRS)
        DO JC= 1, NTENCODE, NPROMA
          JCS = JC
          JCL = MIN(JCS+NPROMA-1, NTENCODE)
          ITHRS = JC/NPROMA + 1
          VALM(ITHRS) = ZMINSPEC 
          DO JJ = JCS, JCL 
            IF(VALUES(JJ) /= ZMISS ) THEN
              VALUES(JJ) = LOG10(VALUES(JJ)+PPEPS)+ABSPPREC
              VALM(ITHRS) = MAX(VALM(ITHRS), VALUES(JJ))
            ELSE
              VALUES(JJ) = ZMINSPEC
            ENDIF
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

        PPMAX=MAXVAL(VALM(:))
        DEALLOCATE(VALM)

        PPMIN=MIN(PPMISS, PPMAX-DELTAPP)
        PPMIN=MIN(PPMIN, PPMIN_RESET)

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JC, JCS, JCL, JJ)
        DO JC= 1, NTENCODE, NPROMA
          JCS = JC
          JCL = MIN(JCS+NPROMA-1, NTENCODE)
          DO JJ = JCS, JCL 
            IF ( VALUES(JJ) <= PPMIN ) VALUES(JJ)=ZMISS
          ENDDO
        ENDDO
!$OMP END PARALLEL DO

      ENDIF

!*    1. FIX PARAMETERS AND PACK DATA.
!        -----------------------------

      CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'editionNumber',IGRIB_VERSION )

!     SPECIFIC VALUES :

!     MISSING DATA:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'missingValue',ZMISS)

!     GRIB TABLE AND PARAMETER NUMBER
      IF (ITABLE == 128) THEN
!       it seems that then default table is not used when defining paramId !
        ITABPAR=IPARAM
      ELSE
        ITABPAR=ITABLE*1000+IPARAM
      ENDIF
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'paramId',ITABPAR,IERR)
      IF (IERR /= 0) THEN
        WRITE(NULERR,*) ' *********************************************'
        WRITE(IU06,*) ' *********************************************'
        WRITE(IU06,*) ' ECCODES ERROR WHILE SETTING paramId ',ITABPAR
        WRITE(NULERR,*) ' ECCODES ERROR WHILE SETTING paramId ',ITABPAR
        WRITE(NULERR,*) ' ECCODES ERROR CODE ', IERR
        WRITE(IU06,*) ' ECCODES ERROR CODE ', IERR
        CALl FLUSH(IU06)
        IF (IERR == -36) THEN
          ITABPAR = 212*1000+IPARAM
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'paramId',ITABPAR)

          WRITE(NULERR,*) ' THE PARAMETER SHOULD BE ADDED TO THE LIST OF'
          WRITE(NULERR,*) ' PARAMETERS KNOWN BY ECCODES !!!' 
          WRITE(NULERR,*) ' IN THE MEAN TIME, THE PROGRAM WILL CONTINUE.' 
          WRITE(NULERR,*) ' USING EXPERIMENTAL PARAMETER TABLE 212' 
          WRITE(NULERR,*) ' WITH paramId= ', ITABPAR
          WRITE(NULERR,*) ' *********************************************'
          WRITE(IU06,*) ' THE PARAMETER SHOULD BE ADDED TO THE LIST OF'
          WRITE(IU06,*) ' PARAMETERS KNOWN BY ECCODES !!!' 
          WRITE(IU06,*) ' USING EXPERIMENTAL PARAMETER TABLE 212' 
          WRITE(IU06,*) ' WITH paramId= ', ITABPAR
          WRITE(IU06,*) ' *********************************************'

        ELSE
          WRITE(NULERR,*) ' *********************************************'
          CALL ABORT1
        ENDIF
      ENDIF

!     LEVEL DEFINITION
      IF (.NOT.LNEWLVTP) THEN
        IF ( IGRIB_VERSION == 1) THEN
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'levtype',ILEVTYPE)
          IF (KLEV /= 0) THEN
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',105)
          ELSE
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'levtype',102)
          ENDIF
        ENDIF
      ENDIF
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'level',KLEV)

      IF (.NOT. LGRHDIFS .OR. & 
     &   (MARSTYPE == 'an' .AND. IFCST == 0) .OR. &
     &   (MARSTYPE == '4v') .OR. &
         (MARSTYPE == 'fg' .AND. IFCST == 0) ) THEN

        READ(CDATE(1:8),'(I8)') IDATE
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'date',IDATE)
        READ(CDATE(9:12),'(I4)') ITIME

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'time',ITIME)

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepType','instant')
!!!   for compatibility with previous coding, impose:
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'timeRangeIndicator',10)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'endStep',IFCST)
        
        IF ( MARSTYPE == 'fg' .AND. IFCST == 0 ) THEN
          ICLASS = 1
        ELSEIF ( MARSTYPE == '4v' ) THEN
          ICLASS = 6
        ELSEIF ( MARSTYPE == 'an' .AND. IFCST == 0 ) THEN
          ICLASS = 2
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'stream',C12)
          IF (C12(1:4) == 'ewla' .OR. C12(1:4) == 'lwwv') THEN
            IDATERES=IDATE*100+ITIME/100
            IY2=IDATERES/1000000
            IDATERES=IDATERES-IY2*1000000
            IM2=IDATERES/10000
            IDATERES=IDATERES-IM2*10000
            ID2=IDATERES/100
            IDATERES=IDATERES-ID2*100
            IH2=IDATERES
            IMN2=0
            ISS2=0
            WRITE(CDATE2,'(I4,5I2)')IY2,IM2,ID2,IH2,IMN2,ISS2

            IDATERES=NDATE_TIME_WINDOW_END
            IF (IDATERES /= 0) THEN
              IY1=IDATERES/1000000
              IDATERES=IDATERES-IY1*1000000
              IM1=IDATERES/10000
              IDATERES=IDATERES-IM1*10000
              ID1=IDATERES/100
              IDATERES=IDATERES-ID1*100
              IH1=IDATERES
              IMN1=0
              ISS1=0
              WRITE(CDATE1,'(I4,5I2)')IY1,IM1,ID1,IH1,IMN1,ISS1
              CALL DIFDATE(CDATE2,CDATE1,NWINOFF)
              NWINOFF=NWINOFF/3600
            ELSE
!             this only works for 12 hour analysis windows !!!1
              IF (IH2+3 == 12 ) THEN
                NWINOFF=0
              ELSE 
                NWINOFF=12-MOD(IH2+3,12)
              ENDIF
              IF (IPARAM == 251) THEN
                CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'localFlag',4)
              ENDIF
            ENDIF
!           in hours
            CALL IGRIB_SET_VALUE(IGRIB_HANDLE, &
             'offsetToEndOf4DvarWindow',NWINOFF)
          ENDIF 
        ELSEIF ( MARSTYPE == 'cf' ) THEN
          ICLASS = 10 
        ELSEIF ( MARSTYPE == 'pf' ) THEN
          ICLASS = 11 
        ELSE
          ICLASS = 9
        ENDIF

        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'type',ICLASS)

        IF (ICLASS /= 9 .AND. ICLASS /= 10 .AND. ICLASS /= 11 &
           .AND. ICLASS /= 6 .AND. IFCST > 0) THEN
          WRITE(IU06,*)' SUB: WGRIBENCODE: THIS IS A FORECAST'
          WRITE(IU06,*)' BUT MARSTYPE DOES NOT KNOW ABOUT IT'
          WRITE(IU06,*)'  '
          WRITE(IU06,*)' CALL ABORT1 '
          WRITE(IU06,*)'  '
          CALL ABORT1
        ENDIF

      ELSE
!       TAKE THE COMMON INFORMATION FROM THE IFS
!       MOST OF IT WAS ALREADY OBTAINED WHEN CLONE OF IFS HANDLE WAS TAKEN
!       IN *PRESET_WGRIB_TEMPLATE*

        IF (LRSTST0) THEN
!         NEED TO TEMPORARLY RESET THE IFS FORECAST STEP
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','s')
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP)
          ISTEP=ISTEP-KCOUSTEP
          ISTEP_HRS=ISTEP/3600
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
          CALL IGRIB_GET_VALUE(IGRIB_HANDLE,'endStep',ISTEP_HRS)
        ENDIF

!!!   for compatibility with previous coding, impose:
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'timeRangeIndicator',10)

        IF ( IGRIB_VERSION == 1 ) THEN
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'unitOfTimeRange',1)
        ELSE
          CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'indicatorOfUnitOfTimeRange',1)
        ENDIF
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'stepUnits','h')
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'endStep',ISTEP_HRS)

      ENDIF

      IF (IPARAM == 251) THEN
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'directionNumber',IK)
        CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'frequencyNumber',IM)
      ENDIF

!     ENCODE DATA:
      CALL IGRIB_SET_VALUE(IGRIB_HANDLE,'values',VALUES)

      DEALLOCATE(VALUES)

      !CALL GSTATS(1709,1)

      IF (LHOOK) CALL DR_HOOK('WGRIBENCODE',1,ZHOOK_HANDLE)

END SUBROUTINE WGRIBENCODE
