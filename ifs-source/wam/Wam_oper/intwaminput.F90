! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

PROGRAM intwaminput

! ----------------------------------------------------------------------

!**** *INTWAMINPUT* -

!     J. BIDLOT    ECMWF  JULY 2000 
!     J. BIDLOT    ECMWF  MARCH 2010 : adapt to use grib_api 

!*    PURPOSE.
!     --------
!     TO TRANSFORM A GRIB INPUT WAVE FIELD (SPECTRA OR OTHERS) INTO
!     A CORRESPONDING FIELD WITH GRID, DIRECTION AND FREQUENCY AS
!     SPECIFIED BY THE MODEL INPUT GRID FILE (gridglou).
!
!     IF THERE IS NOT NEED FOR INTERPOLATION, THE INPUT FILE WILL
!     SIMPLY BE RECOPIED TO THE OUTPUT DESTINATION. NO CODING/DECODING
!     SHOULD OCCUR.

!     IN CASE THE INPUT FIELD IS PARAMETER 251 (WAVE SPECTRA) A FILE
!     NAMED wave_spectral_resolution WILL BE OUTPUT CONTAINING THE
!     WAVE SPECTRAL RESOLUTION (number of directions and frequencies)
!     AND WHETHER OR NOT THE FIELD HAS TO BE INTERPOLATED.

!**   INTERFACE.
!     ----------
!     gridglou must be present and connected to unit 7 (i.e fort.7)
!     The input should be in input_field
!     The output will be placed in output_field

!     METHOD.
!     -------
!     EXTERNALS.
!     ----------
!     REFERENCES
!     ----------
!       NONE.
! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,TH
      USE YOWGRIB_HANDLES , ONLY :NGRIB_HANDLE_WAM_I,NGRIB_HANDLE_WAM_S
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC    ,NGRBRESI ,PPMIN_RESET,  &
     &            NGRBRESS ,HOPERI   ,HOPERS   ,LGRHDIFS ,LNEWLVTP ,    &
     &            LPADPOLES,NGRIB_VERSION
      USE YOWGRID  , ONLY : DELPHI   ,IJS, IJL, NTOTIJ, NPROMA_WAM, NCHNK    
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,ZDELLO   ,NLONRGG  ,IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,    &
     &            KTAG     ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NGX      ,NGY      ,    &
     &            NIBLO    ,CLDOMAIN
      USE YOWPCONS , ONLY : ZMISS    ,EPSMIN
      USE YOWSTAT  , ONLY : MARSTYPE ,YCLASS   ,YEXPVER  ,ISTREAM  ,    &
     &            NLOCGRB , IREFRA   ,NENSFNB  ,NTOTENS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWUNPOOL ,ONLY : LLUNSTR
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE MPL_MODULE,ONLY : MPL_INIT, MPL_END

      USE YOWGRIB


! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "adjust.intfb.h"
#include "grib2wgrid.intfb.h"
#include "iwam_get_unit.intfb.h"
#include "iniwcst.intfb.h"
#include "kgribsize.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "readpre.intfb.h"
#include "wgribenout.intfb.h"
#include "wstream_strg.intfb.h"

      INTEGER(KIND=JWIM) :: NBIT = 1600000

      INTEGER(KIND=JWIM) :: I4(2), KRET, KPLENG, ISIZE, KLEN
      INTEGER(KIND=JWIM) :: ISYSTEM_STAT, ISTAT
      INTEGER(KIND=JWIM) :: KLENG, KLENP, KWORD, NC, NR, NTOT, I, J,    &
     &                      JRGG, IR, JQGAUSS
      INTEGER(KIND=JWIM) :: IREPR, IPARAM, ITABLE, KANG, KFRE, IALLFLD, &
     &                      JSN
      INTEGER(KIND=JWIM) :: NXS, NXE, NYS, NYE
      INTEGER(KIND=JWIM) :: IFR1, KFR1, IFR, IC, JC, MMSHIFT, MLAST
      INTEGER(KIND=JWIM) :: IFORP, KZLEV, K, M, KKK, MMM, JPARAM,       &
     &                      KKKPR, NGYFULL
      INTEGER(KIND=JWIM) :: IPERIODIC, NRFULL, JSTREAM, ISTART, ISTOP,  &
     &                      IDUM
      INTEGER(KIND=JWIM) :: IU07, ILEVTYPE, KSTREAM
      INTEGER(KIND=JWIM) :: KAMONOP, KAMOEAP, KAMOSOP, KAMOWEP
      INTEGER(KIND=JWIM) :: KRMONOP, KRMOEAP, KRMOSOP, KRMOWEP
      INTEGER(KIND=JWIM) :: LFILE, KFILE_HANDLE, KGRIB_HANDLE, IGRIB_LEN
      INTEGER(KIND=JWIM) :: IYYYYMMDD, IHHMM
      INTEGER(KIND=JWIM) :: IVAL
      INTEGER(KIND=JWIM) :: IRET
      INTEGER(KIND=JWIM) :: IPLPRESENT, NB_PL
      INTEGER(KIND=JWIM) :: IDIRSCALING, IFRESCALING
      INTEGER(KIND=JWIM) :: IUOUT
      INTEGER(KIND=JWIM), DIMENSION(:), ALLOCATABLE :: PL
      INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB_BUFR(:)
      INTEGER(KIND=JWIM), ALLOCATABLE ::  NLONRGG_LOC(:)
      INTEGER(KIND=JWIM), ALLOCATABLE ::  KLONRGG(:)

      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL(KIND=JWRB) :: X4(2)
      REAL(KIND=JWRB) :: STEP, START_STEP, END_STEP
      REAL(KIND=JWRB) :: RMONOP, RMOSOP, RMOEAP, RMOWEP
      REAL(KIND=JWRB) :: PRPLRADI
      REAL(KIND=JWRB) ::  DELLO, DELLA
      REAL(KIND=JWRB) :: ONETHIRD, TWOTHIRD
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: FIELD, TEMP
      REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:,:) :: XLON, YLAT
      REAL(KIND=JWRB), ALLOCATABLE :: SCFR(:), SCTH(:)


      CHARACTER(LEN=2) :: MARSFCTYPE
      CHARACTER(LEN=4) :: CSTREAM
      CHARACTER(LEN=8) :: CSTEPTYPE
      CHARACTER(LEN=14) :: CDATE
      CHARACTER(LEN=11) :: IFILENAME
      CHARACTER(LEN=12) :: OFILENAME
      CHARACTER(LEN=12) :: CGRIDTYPE
      CHARACTER(LEN=150):: CMDMSG

      LOGICAL :: LLINTERPOL, LLNONWAVE, LASTREAM
      LOGICAL :: LLFR1OK
      LOGICAL :: LLEXIST
      LOGICAL :: LFDB
      LOGICAL :: LLWAIT

! ----------------------------------------------------------------------

      CALL MPL_INIT()
      IU06 = 6
      IU07 = IWAM_GET_UNIT(IU06, 'wam_grid_tables', 'r', 'u', 0, 'READWRITE')

      WRITE(IU06,*) ' INTWAMINPUT STARTED '

!!!!!  because this program will write its output in grib, it is not yet ready for unsctructured grid
      LLUNSTR=.FALSE.

      PRPLRADI=1.0_JWRB
      CALL INIWCST(PRPLRADI)

!     1.1 INITIALISATION OF VARIABLES WITH DEFAULT VALUES
!         ---------------------------------------------------

      ITEST = 1 
      MARSTYPE = 'an'
      YCLASS   = 'od'
      YEXPVER  = '0001' 
      NENSFNB = 0
      NTOTENS = 0
      ISTREAM = 1045 !!! is changed to an ifs stream also change LNEWLVTP 
      LNEWLVTP=.FALSE.

      IRANK=1
      NPROC=1
      X4=1.0_JWRB
      NPRECR = KIND(X4)
      I4=1
      NPRECI = KIND(I4)
      KTAG=100

      ONETHIRD = 1.0_JWRB/3.0_JWRB
      TWOTHIRD = 2.0_JWRB/3.0_JWRB
      ! reset PPMIN to avoid imposing a minmimum value in *WGRIBENCODE*
      PPMIN_RESET=LOG10(PPEPS)+ABS(PPREC)+EPSMIN 
      LPADPOLES=.FALSE. ! do not pad poles in *WGRIBENCODE*

      LFDB=.FALSE.
      ALLOCATE (NSTART(NPROC),NEND(NPROC))

      KGRIB_HANDLE=0
! ----------------------------------------------------------------------

!*    2. READ PREPROC OUTPUT.
!        --------------------

      CALL READPRE (IU07)

      NINF=1
      NSUP=NIBLO

      NSTART=1
      NEND=IJL
      NTOTIJ = IJL-IJS+1

      NPROMA_WAM = NTOTIJ 
      NCHNK = 1


      IF (LLUNSTR) THEN
        WRITE(*,*) 'NOT YET READY FOR UNSTRUCTURED GRID '
        CALL ABORT1
      ELSE
        NXS = 1
        NXE = NGX
        NYS = 1
        NYE = NGY
        ALLOCATE(NLONRGG_LOC(NGY))
      ENDIF

      IF (.NOT.ALLOCATED(XLON)) ALLOCATE(XLON(NXS:NXE, NYS:NYE))
      XLON(:,:)=ZMISS
      IF (.NOT.ALLOCATED(YLAT)) ALLOCATE(YLAT(NXS:NXE, NYS:NYE))
      YLAT(:,:)=ZMISS

      IF (LLUNSTR) THEN
        WRITE(*,*) 'NOT YET READY FOR UNSTRUCTURED GRID '
!!!!    still too many NGX, NGY
        CALL ABORT1
      ELSE
        NLONRGG_LOC(:)=NLONRGG(:)
!$OMP   PARALLEL DO SCHEDULE(STATIC) PRIVATE(J,I,JSN)
        DO J = NYS, NYE
          JSN = NGY-J+1
          DO I = MAX(1,NXS), MIN(NLONRGG(JSN),NXE)
            XLON(I,J) = AMOWEP + (I-1)*ZDELLO(JSN)
            YLAT(I,J) = AMOSOP + (JSN-1)*XDELLA
          ENDDO
        ENDDO
!$OMP   END PARALLEL DO
      ENDIF


      IFILENAME='input_field'
      OFILENAME='output_field'
  
!     CONNECT TO INPUT FILE
      LFILE=0
      LLEXIST=.FALSE.
      IF (IFILENAME /= ' ') LFILE=LEN_TRIM(IFILENAME)
      INQUIRE(FILE=IFILENAME(1:LFILE),EXIST=LLEXIST)
      IF (LLEXIST) THEN
        CALL IGRIB_OPEN_FILE(KFILE_HANDLE,IFILENAME(1:LFILE),'r')
      ELSE
        WRITE(*,*)'****************************'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'*GRIB DATA NOT FOUND IN *'
        WRITE(*,*)  IFILENAME 
        WRITE(*,*)'*PROGRAM WILL ABORT        *'
        WRITE(*,*)'*                          *'
        WRITE(*,*)'****************************'
        CALL ABORT1
      ENDIF

!     LOAD THE DATA
1021  ISIZE=NBIT
      KBYTES=ISIZE*NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))

      CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,KGRIB_BUFR,KBYTES,IRET)

      IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!     *IGRIB_READ_FROM_FILE* does not read through the file if
!!!     the size is too small, so figure out the size and read again.
        CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INTWAMINPUT')
        DEALLOCATE(KGRIB_BUFR)
        GOTO 1021
      ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) '* INTWAMINPUT: END OF FILE ENCOUNTED'
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
        WRITE(IU06,*) '*************************************'
        WRITE(IU06,*) '* INTWAMINPUT: FILE HANDLING ERROR'
        WRITE(IU06,*) '*************************************'
        CALL ABORT1
      ENDIF

      KGRIB_HANDLE=-99
      CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB_BUFR)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Nj',NR)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'gridType', CGRIDTYPE)
      IF (CGRIDTYPE(1:10) == 'regular_gg') THEN
        JRGG=1
        IREPR=4
        JQGAUSS=1
      ELSEIF (CGRIDTYPE(1:10) == 'reduced_gg') THEN
        JRGG=1
        IREPR=4
        JQGAUSS=1
      ELSEIF (CGRIDTYPE(1:7) == 'regular') THEN
        JRGG=0
        IREPR=0
        JQGAUSS=0
      ELSEIF (CGRIDTYPE(1:7) == 'reduced') THEN
        JRGG=1
        IREPR=0
        JQGAUSS=0
      ELSE
        WRITE(IU06,*) '***********************************'
        WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!!'
        WRITE(IU06,*) '   gridType = ', CGRIDTYPE 
        WRITE(IU06,*) '***********************************'
        CALL ABORT1
      ENDIF

      IF (JRGG == 1) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'PLPresent',IPLPRESENT)
        IF (IPLPRESENT == 1) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'numberOfPointsAlongAMeridian',NB_PL)
          ALLOCATE(PL(NB_PL))
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'pl',PL)
        ELSE
          WRITE(IU06,*) 'NUMBER OF POINTS PER LATITUDE MISSING !!!'
          CALL ABORT1
        ENDIF
        NC=0
        DO J=1,NB_PL
          NC = MAX(NC,PL(J))
        ENDDO
        IR=0
        DO J=1,NB_PL
          IF (PL(J) /= 0) IR=IR+1
        ENDDO
        NR=IR
      ELSEIF (JRGG == 0) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'Ni',IVAL)
        NC=IVAL
      ELSE
        WRITE(IU06,*) '   STRUCTURE OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

!*    DETERMINE CODE FOR DATA FIELD TYPE.
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'paramId',IVAL)
      ITABLE=IVAL/1000
      IPARAM=IVAL-ITABLE*1000

!     DATE. 
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'dataDate',IYYYYMMDD)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'time',IHHMM)
      WRITE(CDATE(1:12),'(I8.8,I4.4)') IYYYYMMDD,IHHMM 
      CDATE(13:14)='00'

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stepType',CSTEPTYPE)
!     FORECAST STEP (in seconds)
      IF (CSTEPTYPE(1:7) == 'instant') THEN
        CALL IGRIB_SET_VALUE(KGRIB_HANDLE,'stepUnits','s')
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'step',STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'startStep',START_STEP)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'endStep',END_STEP)
!       THE DATA ARE VALID BETWEEN TWO TIMES. TAKE THE MIDDLE POINT
        IF (START_STEP /= END_STEP) THEN
          STEP=(END_STEP-START_STEP)/2
        ENDIF
        IFORP=STEP
      ELSE
        WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP TYPE !!!'
        WRITE(*,*) 'stepType = ',CSTEPTYPE
        CALL ABORT1
      ENDIF


!     DETERMINE GRID PARAMETERS.

      IF (.NOT.ALLOCATED(KLONRGG)) ALLOCATE(KLONRGG(NR))
      KLONRGG=0

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfFirstGridPointInDegrees',RMONOP)
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'latitudeOfLastGridPointInDegrees',RMOSOP)

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'longitudeOfFirstGridPointInDegrees',RMOWEP)

!!!   THERE IS A DANGER THAT THE DEFINITON FOR RMOEAP MIGHT VARY DUE TO
!!!   THE AMBIGOUS DEFINITION FOR IRREGULAR GRIDS. FOR NON WAVE FIELDS,
!!!   A GAUSSIAN GRID IMPLIES THAT THE GRID IS GLOBAL, THEREFORE
!!!   RMOEAP IS IMPLICITLY KNOWN.
      CSTREAM='****'
      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'levtype',ILEVTYPE, KRET=IRET)
      IF (IRET /= JPGRIB_SUCCESS) ILEVTYPE=0

      CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'stream',JSTREAM)
      CALL WSTREAM_STRG(JSTREAM,CSTREAM,NENSFNB,NTOTENS,MARSFCTYPE,     &
     &                  KSTREAM,LASTREAM)

      IF (CSTREAM == '****' .OR.                                        &
     &    (LASTREAM .AND. ILEVTYPE /= 209 .AND. ILEVTYPE /= 212 )) THEN 
        LLNONWAVE=.TRUE.
      ELSE
        LLNONWAVE=.FALSE.
      ENDIF

      IF (IREPR == 4 .AND. LLNONWAVE) THEN
        DELLO = 360.0_JWRB/MAX(1,NC)
        RMOEAP = RMOWEP+360.0_JWRB - DELLO
        IPERIODIC = 1
      ELSE
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'longitudeOfLastGridPointInDegrees',RMOEAP)
        IF (JRGG == 1) CALL ADJUST (RMOWEP, RMOEAP)
        IPERIODIC = 0
        DELLO=(RMOEAP-RMOWEP)/MAX(1,NC-1)
        IF (RMOEAP-RMOWEP+1.5_JWRB*DELLO >= 360.0_JWRB) IPERIODIC = 1
      ENDIF

      IF (JRGG == 1) THEN
        ISTART=1
        DO WHILE(PL(ISTART).EQ.0 .AND. ISTART.LT.NB_PL)
          ISTART=ISTART+1
        ENDDO
        ISTART=ISTART-1

        ISTOP=0
        DO WHILE(PL(NB_PL-ISTOP).EQ.0 .AND. ISTOP.LT.NB_PL)
          ISTOP=ISTOP+1
        ENDDO

        DO J=1,NR-ISTART
          JSN=NR-J+1
          KLONRGG(JSN) = PL(J+ISTART)
        ENDDO
        DEALLOCATE(PL)

        IF( JQGAUSS.NE.1 ) THEN
          CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'jDirectionIncrementInDegrees',DELLA)
          RMONOP =  RMONOP-ISTART*DELLA
          RMOSOP =  RMOSOP+ISTOP*DELLA
        ELSE
          DELLA = 0.0_JWRB
        ENDIF

      ELSEIF (JRGG == 0) THEN
        KLONRGG=NC
      ELSE
        WRITE(IU06,*) ' REPRESENTATION OF THE FIELD NOT KNOWN'
        CALL ABORT1
      ENDIF

      IF (IPARAM == 251) THEN
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberOfDirections',KANG)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'numberOfFrequencies',KFRE)

!       input data first frequency
        IF (.NOT.ALLOCATED(SCFR)) ALLOCATE(SCFR(KFRE))
        IF (.NOT.ALLOCATED(SCTH)) ALLOCATE(SCTH(KANG))

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'directionScalingFactor',IDIRSCALING)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'scaledDirections',SCTH)

        CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'frequencyScalingFactor',IFRESCALING)
        CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'scaledFrequencies',SCFR)

        KFR1=NINT(SCFR(1))
!       model first frequency
        IFR1=NINT(IFRESCALING*FR(1))

        IF (IFR1 /= KFR1 ) THEN
          LLFR1OK=.FALSE.
          IF (IFR1 < KFR1) THEN
            IC=1
            IFR=IFR1
            DO WHILE (IFR.LT.KFR1 .AND. IC.LT.NFRE) 
              IC=IC+1
              IFR=NINT(IFRESCALING*FR(IC))
            ENDDO
            MMSHIFT=IC-1
          ELSE
            IC=1
            IFR=KFR1
            DO WHILE (IFR.LT.IFR1 .AND. IC.LT.NFRE) 
              IC=IC+1
              IFR=NINT(SCFR(IC))
            ENDDO
            MMSHIFT=-(IC-1)
          ENDIF
        ELSE
          LLFR1OK=.TRUE.
          MMSHIFT=0
        ENDIF

      ENDIF

!      FIND WHETHER INTERPOLATION IS NEEDED
!      ------------------------------------

       LLINTERPOL=.TRUE.

       NGYFULL=INT(180.0_JWRB/XDELLA)+1

       KAMONOP=NINT(AMONOP*100.0_JWRB)
       KRMONOP=NINT(RMONOP*100.0_JWRB)
       KAMOSOP=NINT(AMOSOP*100.0_JWRB)
       KRMOSOP=NINT(RMOSOP*100.0_JWRB)
       KAMOWEP=NINT(AMOWEP*100.0_JWRB)
       KRMOWEP=NINT(RMOWEP*100.0_JWRB)
       KAMOEAP=NINT(AMOEAP*100.0_JWRB)
       KRMOEAP=NINT(RMOEAP*100.0_JWRB)
       IF ((KAMONOP == KRMONOP .OR. KRMONOP == 9000)  .AND.             &
     &     (KAMOSOP == KRMOSOP .OR. KRMOSOP == -9000) .AND.             &
     &     KAMOWEP == KRMOWEP .AND. KAMOEAP == KRMOEAP     ) THEN
         IF (JRGG == IRGG .AND. (NR == NGY .OR. NR == NGYFULL)) THEN
           IF (IRGG == 1) THEN
             LLINTERPOL=.FALSE.
             DO J=1,NGY
               IF (KLONRGG(J) /= NLONRGG(J)) THEN
                 LLINTERPOL=.TRUE.
                 EXIT
               ENDIF
             ENDDO
           ELSE
             IF (NC == NGX) LLINTERPOL=.FALSE.
           ENDIF
         ENDIF
       ENDIF
       IF (LLINTERPOL) THEN
         WRITE(IU06,*) ' SPATIAL INTERPOLATION REQUIRED '
         WRITE(IU06,*) ' RMONOP = ',RMONOP,' -> AMONOP = ',AMONOP
         WRITE(IU06,111) RMONOP,AMONOP
         WRITE(IU06,*) ' RMOSOP = ',RMOSOP,' -> AMOSOP = ',AMOSOP
         WRITE(IU06,111) RMOSOP,AMOSOP
         WRITE(IU06,*) ' RMOWEP = ',RMOWEP,' -> AMOWEP = ',AMOWEP
         WRITE(IU06,111) RMOWEP,AMOWEP
         WRITE(IU06,*) ' RMOEAP = ',RMOEAP,' -> AMOEAP = ',AMOEAP
         WRITE(IU06,111) RMOEAP,AMOEAP
         WRITE(IU06,*) ' JRGG = ',JRGG,' -> IRGG = ',IRGG
         WRITE(IU06,*) ' NR = ',NR,' -> NGY = ',NGY
         WRITE(IU06,*) ' NC = ',NC,' -> NGX = ',NGX
       ENDIF
111    FORMAT(4x,'HEX: ',3(Z16.16,2x))

       IF (IPARAM == 251) THEN
         IF (KANG /= NANG .OR. KFRE /= NFRE .OR. IFR1 /= KFR1 ) THEN
!        it's assumed here that we will not change the initial direction
!        without changing the total number of directions
           LLINTERPOL=.TRUE.
         ENDIF

         open(17,FILE='wave_spectral_resolution')
         WRITE(17,*) 'NANG ', KANG
         WRITE(17,*) 'NFRE ', KFRE
         WRITE(17,*) 'INTERPOL ', LLINTERPOL
         close(17)
       ENDIF

       IF (.NOT.LLINTERPOL) THEN
!      THERE IS NO NEED FOR INTERPOLATION, THE INPUT FILE
!      WILL SIMPLY BE RECOPIED TO THE OUTPUT DESTINATION.
         WRITE(IU06,*) ' THERE IS NO NEED FOR INTERPOLATION'
         WRITE(IU06,*) ' THE INPUT FILE WILL SIMPLY BE RECOPIED'
         WRITE(IU06,*) ' TO THE OUTPUT DESTINATION.'

         LLWAIT=.TRUE.
         ISYSTEM_STAT=-1
         ISTAT=-1
         CMDMSG="Unchanged CMDMSG"
#if defined(__PGI)
! Not supported by the PGI compiler yet
#else
         CALL EXECUTE_COMMAND_LINE('cp '//IFILENAME//' '//OFILENAME,    &
     &        WAIT=LLWAIT, EXITSTAT=ISYSTEM_STAT,CMDSTAT=ISTAT,CMDMSG=CMDMSG)
#endif
         IF (ISTAT /= 0 .OR. ISYSTEM_STAT /= 0) THEN
           WRITE(IU06,*) '*************************************'
           WRITE(IU06,*) '* INTWAMINPUT: FILE COPYING ERROR !'
           WRITE(IU06,*) '* ISTAT = ',ISTAT
           WRITE(IU06,*) '* ISYSTEM_STAT = ',ISYSTEM_STAT
           WRITE(IU06,*) '* CMDMSG = ',TRIM(CMDMSG)
           WRITE(IU06,*) '*************************************'
           CALL ABORT1
         ENDIF

         GOTO 8888 
       ELSE
!        PREPARE OUTPUT
         CALL IGRIB_GET_VALUE(KGRIB_HANDLE, 'localDefinitionNumber',NLOCGRB)

         CALL IGRIB_GET_VALUE(KGRIB_HANDLE,'editionNumber', NGRIB_VERSION)
!!!!! We might need to impose girb2 at later stage

!        GRIB HANDLES FOR OUTPUT
         LGRHDIFS=.FALSE.
!        FOR INTEGRATED PARAMETERS
         CALL PRESET_WGRIB_TEMPLATE("I",NGRIB_HANDLE_WAM_I)
!        FOR SPECTRA
         CALL PRESET_WGRIB_TEMPLATE("S",NGRIB_HANDLE_WAM_S)

         CALL IGRIB_OPEN_FILE(IUOUT,OFILENAME,'w')
       ENDIF

       IF (IPARAM == 250) THEN
         WRITE(IU06,*) '***************************************'
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '* THE FIELDS OF PARAMETER 250 MUST BE *'
         WRITE(IU06,*) '* CONVERTED TO PARAMETER 251 USING    *'
         WRITE(IU06,*) '* convert_grbspec                     *' 
         WRITE(IU06,*) '*                                     *'
         WRITE(IU06,*) '***************************************'
         CALL ABORT1

       ELSEIF (IPARAM == 251) THEN

         WRITE(IU06,*) ' '
         IF (KANG == NANG .AND. KFRE == NFRE .AND. LLFR1OK ) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL OR FREQUENCY INTERPOLATION',  &
     &                   ' REQUIRED. '
         ELSEIF (KANG == NANG .AND. .NOT.LLFR1OK ) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL INTERPOLATION REQUIRED'
           WRITE(IU06,*) ' BUT FREQUENCY INTERPOLATION IS REQUIRED'
           WRITE(IU06,*) ' BECAUSE OF DIFFERENT FIRST FREQUENCY'
           WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
           WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
         ELSEIF (KANG == 2*NANG .OR. 2*KANG == NANG .OR.                &
     &           KANG == 3*NANG .OR. 3*KANG == NANG .OR.                &
     &           4*KANG == 3*NANG .OR.                                  &
     &         2*KANG == 3*NANG .OR. 3*KANG == 2*NANG) THEN
           WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION REQUIRED.'
           WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
           WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
           WRITE(IU06,*) ' '
           IF (KFRE == NFRE .AND. LLFR1OK ) THEN
             WRITE(IU06,*) ' BUT NO FREQUENCY INTERPOLATION REQUIRED.'
           ELSEIF (.NOT.LLFR1OK .OR. KFRE /= NFRE ) THEN
             WRITE(IU06,*) ' AND FREQUENCY INTERPOLATION ALSO REQUIRED'
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
             WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
           ELSE
             WRITE(IU06,*) '*******************************************'
             WRITE(IU06,*) ' FREQUENCY REQUIRED.'
             WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW WHAT TO DO !'
             WRITE(IU06,*) ' IT ABORTS '
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) '*******************************************'
             CALL ABORT1
           ENDIF
         ELSEIF (KANG == NANG .AND. KFRE /= NFRE) THEN
           WRITE(IU06,*) ' NO DIRECTIONAL INTERPOLATION REQUIRED.'
           IF (.NOT.LLFR1OK .OR. KFRE /= NFRE ) THEN
             WRITE(IU06,*) ' BUT FREQUENCY INTERPOLATION IS REQUIRED.'
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) ' INPUT FIRST FREQUENCY: ',KFR1
             WRITE(IU06,*) ' MODEL FIRST FREQUENCY: ',IFR1
           ELSE
             WRITE(IU06,*) '*******************************************'
             WRITE(IU06,*) ' FREQUENCY INTERPOLATION IS REQUIRED.'
             WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW WHAT TO DO !'
             WRITE(IU06,*) ' IT ABORTS '
             WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
             WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
             WRITE(IU06,*) '*******************************************'
             CALL ABORT1
           ENDIF
         ELSE
           WRITE(IU06,*) '*********************************************'
           WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION IS REQUIRED.'
           WRITE(IU06,*) ' BUT THE PROGRAM DOES NOT KNOW HOW TO PROCEED'
           WRITE(IU06,*) ' IT ABORTS ' 
           WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
           WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
           WRITE(IU06,*) ' INPUT NUMBER OF FREQUENCIES: ',KFRE
           WRITE(IU06,*) ' MODEL NUMBER OF FREQUENCIES: ',NFRE
           WRITE(IU06,*) '*********************************************'
           CALL ABORT1
         ENDIF
         WRITE(IU06,*) ' '

         IF (.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NXS:NXE, NYS:NYE))
         IF (.NOT.ALLOCATED(TEMP)) ALLOCATE(TEMP(NGX,NGY))

         KKKPR=0

         IF (.NOT. LLFR1OK) THEN
!        SET MODEL VALUES TO MISSING FOR THE FIRST MMSHIFT FREQUENCY BINS
!        (WHEN MMSHIFT>0) SINCE NO INPUT DATA AVAILABLE FOR THOSE BINS.
           FIELD=ZMISS
           
           DO M=1,MMSHIFT
             DO K=1,NANG
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ENDDO
           ENDDO
         ENDIF

         DO IALLFLD=1,KANG*KFRE

!          GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

           CALL GRIB2WGRID (IU06, NPROMA_WAM,                           &
     &                      KGRIB_HANDLE, KGRIB_BUFR, ISIZE,            &
     &                      LLUNSTR,                                    &
     &                      NGY, IRGG, NLONRGG_LOC,                     &
     &                      NXS, NXE, NYS, NYE,                         &
     &                      XLON, YLAT,                                 &
     &                      ZMISS, PPREC, PPEPS,                        &
     &                      CDATE, IFORP, JPARAM, KZLEV,KKK, MMM, FIELD)


           DEALLOCATE(KGRIB_BUFR)

!          FREQUENCY INTERPOLATION
!          -----------------------
           M=MMM+MMSHIFT
           MLAST=M
           IF (M <= 0 .OR. M > NFRE) GOTO 2000

!          DIRECTION INTERPOLATION (piecewise interpolation)
!          -----------------------
           IF (KANG == NANG ) THEN
!            NO SPECTRAL INTERPOLATION REQUIRED
             K=KKK
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

           ELSEIF (2*KANG == NANG) THEN
!            DOUBLING THE NUMBER OF DIRECTIONS
             K=2*KKK-1
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             K=2*KKK
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

           ELSEIF (3*KANG == NANG) THEN
!            TRIPLING THE NUMBER OF DIRECTIONS
             K=3*KKK-2
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             K=3*KKK-1
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             K=3*KKK
             CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,              &
     &                     ITABLE, IPARAM, 0, K , M,                    &
     &                     CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

           ELSEIF (3*KANG == 2*NANG) THEN
!            3/2 OF THE NUMBER OF DIRECTIONS
             IF (MOD(KKK,2) == 1) THEN

!              KEEP INPUT (!! TO AVERAGE OUT FOR NEXT OUTPUT BIN)
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     TEMP(I,J)=0.5_JWRB*FIELD(I,J)
                   ELSE
                     TEMP(I,J)=ZMISS
                   ENDIF
                 ENDDO
               ENDDO

               K=3*KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ELSE
!              AVERAGE OUT THE INPUT WITH PREVIOUSLY SAVED ONE 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=TEMP(I,J)+0.5_JWRB*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=0.5_JWRB*FIELD(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO

               K=3*KKK/2-1
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,             &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
               K=3*KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ENDIF

           ELSEIF (4*KANG == 3*NANG) THEN
!            4/3 OF THE NUMBER OF DIRECTIONS
             IF (MOD(KKK,3) == 1) THEN
!              KEEP INPUT (!! TO AVERAGE OUT FOR NEXT OUTPUT BIN)
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     TEMP(I,J) = ONETHIRD*FIELD(I,J)
                   ELSE
                     TEMP(I,J) = ZMISS
                   ENDIF
                 ENDDO
               ENDDO

               K=4*KKK/3
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

             ELSEIF (MOD(KKK,3) == 2) THEN
!              AVERAGE OUT THE INPUT WITH PREVIOUSLY SAVED ONE 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=TEMP(I,J)+TWOTHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=TWOTHIRD*FIELD(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO

               K=4*KKK/3
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,             &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

!              KEEP INPUT (!! TO AVERAGE OUT FOR NEXT OUTPUT BIN)
!              (reuse TEMP)
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     TEMP(I,J) = TWOTHIRD*FIELD(I,J)
                   ELSE
                     TEMP(I,J) = ZMISS
                   ENDIF
                 ENDDO
               ENDDO

             ELSE
!              AVERAGE OUT THE INPUT WITH PREVIOUSLY SAVED ONE 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=TEMP(I,J)+ONETHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=ONETHIRD*FIELD(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=4*KKK/3-1
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,             &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
               K=4*KKK/3
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ENDIF

           ELSEIF (KANG == 2*NANG) THEN
!            HALVING THE NUMBER OF DIRECTIONS
             IF (MOD(KKK,2) == 1) THEN
!              SAVE DIRECTIONAL CONTRIBUTION TO BE COMBINED
!              WITH THE NEXT ONE.
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF (FIELD(I,J) /= ZMISS) THEN
                    TEMP(I,J)=FIELD(I,J)
                  ELSE
                    TEMP(I,J)=0.00_JWRB
                  ENDIF
                ENDDO
               ENDDO

               KKKPR=KKK
             ELSEIF (KKK == KKKPR+1) THEN
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF (FIELD(I,J) /= ZMISS) THEN
                    FIELD(I,J)=0.5_JWRB*(FIELD(I,J)+TEMP(I,J))
                  ELSE
                    IF (TEMP(I,J) /= 0.0_JWRB) THEN
                      FIELD(I,J)=0.5_JWRB*TEMP(I,J) 
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/2
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)


             ELSE
               WRITE(IU06,*) '************************************'
               WRITE(IU06,*) ' FOR THIS TYPE OF'
               WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION,'
               WRITE(IU06,*) ' THE DIRECTIONS IN THE INPUT FIELD ' 
               WRITE(IU06,*) ' MUST BE CONCECUTIVE, STARTING WITH'
               WRITE(IU06,*) ' THE FIRST ONE !!!!' 
               WRITE(IU06,*) ' THE PROGRAM ABORTS ' 
               WRITE(IU06,*) '*************************************'
               CALL ABORT1
             ENDIF
           ELSEIF (KANG == 3*NANG) THEN
!            1/3 OF THE NUMBER OF DIRECTIONS
             IF (MOD(KKK,3) == 1) THEN
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=0.0_JWRB
                 ENDDO
               ENDDO
             ENDIF
             IF (MOD(KKK,3) /= 0) THEN
!              SAVE DIRECTIONAL CONTRIBUTION TO BE COMBINED
!              WITH THE NEXT ONE.
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF (FIELD(I,J) /= ZMISS) THEN
                    TEMP(I,J)=TEMP(I,J)+ONETHIRD*FIELD(I,J)
                  ENDIF
                ENDDO
               ENDDO

             ELSE
               DO J=1,NGY
                JSN=NGY-J+1
                DO I=1,NLONRGG(JSN)
                  IF (FIELD(I,J) /= ZMISS) THEN
                    FIELD(I,J)=ONETHIRD*FIELD(I,J)+TEMP(I,J)
                  ELSE
                    IF (TEMP(I,J) /= 0.0_JWRB) THEN
                      FIELD(I,J)=TEMP(I,J)
                    ENDIF
                  ENDIF
                ENDDO
               ENDDO

               K=KKK/3
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ENDIF

           ELSEIF (2*KANG == 3*NANG) THEN
!            2/3 OF THE NUMBER OF DIRECTIONS
             IF (MOD(KKK,3) == 1) THEN
!              KEEP INPUT
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=FIELD(I,J)
                 ENDDO
               ENDDO

             ELSEIF (MOD(KKK,3) == 2) THEN
!              ADD CONTRIBUTION OF PREVIOUSLY SAVED TO INPUT 
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=TWOTHIRD*TEMP(I,J)+ONETHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=ONETHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=TWOTHIRD*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)+1
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,             &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

!              KEEP INPUT
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   TEMP(I,J)=FIELD(I,J)
                 ENDDO
               ENDDO

             ELSE
!              ADD INPUT TO PREVIOUSLY SAVED
               DO J=1,NGY
                 JSN=NGY-J+1
                 DO I=1,NLONRGG(JSN)
                   IF (FIELD(I,J) /= ZMISS) THEN
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)+TWOTHIRD*FIELD(I,J)
                     ELSE
                       TEMP(I,J)=TWOTHIRD*FIELD(I,J)
                     ENDIF
                   ELSE
                     IF (TEMP(I,J) /= ZMISS) THEN
                       TEMP(I,J)=ONETHIRD*TEMP(I,J)
                     ENDIF
                   ENDIF
                 ENDDO
               ENDDO
               K=2*(KKK/3)
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, TEMP,             &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

             ENDIF

           ELSE
             WRITE(IU06,*) '************************************'
             WRITE(IU06,*) ' NOT READY FOR THIS TYPE OF'
             WRITE(IU06,*) ' DIRECTIONAL INTERPOLATION,'
             WRITE(IU06,*) ' INPUT NUMBER OF DIRECTIONS: ',KANG 
             WRITE(IU06,*) ' MODEL NUMBER OF DIRECTIONS: ',NANG 
             WRITE(IU06,*) '*************************************'
             CALL ABORT1
           ENDIF



2000       CONTINUE
!          CONTINUE READING AND DECODING UNTIL LAST FIELD IS READ.
!          -------------------------------------------------------
           IF (IALLFLD == KANG*KFRE) EXIT

           CALL IGRIB_RELEASE(KGRIB_HANDLE)
           KGRIB_HANDLE=0

2021       ISIZE=NBIT
           KBYTES=ISIZE*NPRECI
           ALLOCATE(KGRIB_BUFR(ISIZE))
           CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE,KGRIB_BUFR,           &
     &                               KBYTES,IRET)
           IF (IRET == JPGRIB_BUFFER_TOO_SMALL) THEN
!!!          *IGRIB_READ_FROM_FILE* does not read through the file if
!!!          the size is too small, so figure out the size and read again.
             CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'INTWAMINPUT 2')
             DEALLOCATE(KGRIB_BUFR)
             GOTO 2021
           ELSEIF (IRET == JPGRIB_END_OF_FILE) THEN
             IF (IALLFLD == 1) THEN
               WRITE(IU06,*) ' ONLY ONE SPECTRAL FIELD HAS BEEN READ!'
               EXIT
             ELSE
               WRITE(IU06,*) '**********************************'
               WRITE(IU06,*) '*   END OF FILE ENCOUNTED'
               WRITE(IU06,*) '**********************************'
               CALL ABORT1
             ENDIF
           ELSEIF (IRET /= JPGRIB_SUCCESS) THEN
             WRITE(IU06,*) '*************************************'
             WRITE(IU06,*) '* INTWAMINPUT: FILE HANDLING ERROR'
             WRITE(IU06,*) '*************************************'
             CALL ABORT1
           ENDIF

           KGRIB_HANDLE=-99
           CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,KGRIB_BUFR)

         ENDDO


!        SET MODEL VALUES TO MISSING IF NO VALUES WERE PROVIDED 
!        BY INPUT DATA 
         IF (IALLFLD > 1) THEN
           FIELD=ZMISS
           DO M=MLAST+1,NFRE
             DO K=1,NANG
               CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,            &
     &                       ITABLE, IPARAM, 0, K , M,                  &
     &                       CDATE, IFORP, MARSTYPE, LFDB, IUOUT)
             ENDDO
           ENDDO
         ENDIF

       ELSE
!      ALL THE OTHER WAVE PARAMETERS

!        GET THE DATA AND PERFORM INTERPOLATION TO WAVE MODEL GRID.

         IF (.NOT.ALLOCATED(FIELD)) ALLOCATE(FIELD(NGX,NGY))
         CALL GRIB2WGRID (IU06, NPROMA_WAM,                             &
     &                    KGRIB_HANDLE, KGRIB_BUFR, ISIZE,              &
     &                    LLUNSTR,                                      &
     &                    NGY, IRGG, NLONRGG_LOC,                       &
     &                    NXS, NXE, NYS, NYE,                           &
     &                    XLON, YLAT,                                   &
     &                    ZMISS, PPREC, PPEPS,                          &
     &                    CDATE, IFORP, JPARAM, KZLEV, KKK, MMM, FIELD)

         DEALLOCATE(KGRIB_BUFR)


!        ENCODE AND OUTPUT INTEGRATED PARAMETER
!        --------------------------------------

         KKK=0
         MMM=0

         CALL WGRIBENOUT(IU06, ITEST, NGX, NGY, FIELD,                  &
     &                 ITABLE, IPARAM, KZLEV, KKK , MMM,                &
     &                 CDATE, IFORP, MARSTYPE, LFDB, IUOUT)

       ENDIF

8888  CONTINUE

      CALL IGRIB_RELEASE(KGRIB_HANDLE)
      CALL IGRIB_CLOSE_FILE(KFILE_HANDLE)
      IF (LLINTERPOL) CALL IGRIB_CLOSE_FILE(IUOUT)
 
      CALL MPL_END()
      WRITE (IU06,*) ' PROGRAM INTWAMINPUT: ALL DONE'

END PROGRAM intwaminput
