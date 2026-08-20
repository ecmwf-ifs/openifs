! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE READSTA (IUIN, CBDT, CEDT, ANALPD, FOREPD, IS,         &
     &                    CRSDT, CLSDT, CABDT, CAEDT, IASS, NF, ISTAT,  &
     &                    CCURA, LRSTPARAL, NPROC_RST )
                                                                        
! -------------------------------------------------------------------   

!****  *READSTA*  READS THE WAMINFO FILE.

!     PURPOSE.
!     --------

!         *READSTA* - READS THE WAMINFO FILE AND DETERMINES 
!                     BEGIN DATE TIME, END DATE TIME, ANALYSIS PERIOD,
!                     FORECAST PERIOD, SHIFT,
!                     PLOTTING DATES, IASS, NF, ISTATUS
!                     FROM THE PREVIOUS RUN AND WRITES THIS
!                     INFO ON ARRAY CARD
!                     AND THE LAST DATE FOR OUTPUT OF RESTART FILES
!                     IF SPECIFIED OTHERWISE IT WILL BE SET TO END DATE

!     PROGRAMMER.
!     -----------

!          P. JANSSEN
!          J. BIDLOT NOVEMBER 96  : add CRSDT and CLSDT definition
!          B. HANSEN APRIL    97  : move COMMON CARD INTO INCLUDE FILE.

!**   INTERFACE.
!     ----------

!          IUIN     INTEGER    DISC UNIT FOR WAMINFO.
!          CBDT     CHAR*14    BEGIN DATE-TIME-GROUP.
!          CEDT     CHAR*14    END DATE-TIME-GROUP.
!          ANALPD   INTEGER    ANALYSIS PERIOD IN SECONDS.
!          FOREPD   INTEGER    FORECAST PERIOD IN SECONDS.
!          IS       INTEGER    WIND INPUT TIME STEP IN SECONDS.
!          CRSDT    CHAR*14    DATE FOR RESTART FILES OUTPUT
!          CLSDT    CHAR*14    LAST DATE FOR SPECTRUM FILE OUTPUT
!          CABDT    CHAR*14    BEGIN ANALYSIS DATE.
!          CAEDT    CHAR*14    END ANALYSIS DATE.
!          IASS     INTEGER    ASSIMILATION CONDITION.                  
!          NF       INTEGER    NEW FORECAST CONDITION.                  
!          ISTAT    INTEGER    STATUS FLAGS FOR JOBS.                   
!         *CCURA*   CHAR*14    BEGIN DATE FOR USE OF CURRENTS. 
!          LRSTPARAL LOGICAL   TRUE IF RESTART FILES WRITTEN IN PARALLEL
!          NPROC_RST INTEGER   NUMBER OF PARALLEL MPI TASKS THAT HAVE BEEN USED
                                                                        


!     METHOD.
!     -------                                                           

!          READ THE WAMINFO FILE

!     EXTERNALS.
!     ----------

!          ABORT1    TERMINATES PROCESSING

!     REFERENCES.
!     -----------

!          NONE

!---------------------------------------------------------------------- 

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWCARD  , ONLY : JPCL     ,CARD
      USE YOWTEST  , ONLY : IU06 

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!---------------------------------------------------------------------- 
      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IUIN
      INTEGER(KIND=JWIM), INTENT(OUT) :: ANALPD, FOREPD, IS, IASS, NF 
      INTEGER(KIND=JWIM), INTENT(OUT) :: NPROC_RST
      INTEGER(KIND=JWIM), INTENT(OUT) :: ISTAT(3)
      CHARACTER(LEN=14), INTENT(OUT) :: CBDT, CEDT
      CHARACTER(LEN=14), INTENT(OUT) :: CRSDT, CLSDT, CABDT, CAEDT 
      CHARACTER(LEN=14), INTENT(OUT) :: CCURA 
      LOGICAL, INTENT(OUT) :: LRSTPARAL

      INTEGER(KIND=JWIM) :: I, L
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!---------------------------------------------------------------------- 

      IF (LHOOK) CALL DR_HOOK('READSTA',0,ZHOOK_HANDLE)


!*    1.0  READ WAMINFO FILE FROM IUIN.
!          ---------------------------- 

      CARD=' '
      CRSDT = ' '
      REWIND(IUIN)
      DO I=1,JPCL
        READ(IUIN,'(A72)', END=3000, ERR=3000) CARD(I)
      ENDDO
      REWIND(IUIN)

!*    2.0  CONVERT CHARACTER TO INTEGER.
!          -----------------------------

      READ (CARD(1), '(15X,A14,4X,A14)') CBDT, CEDT
      READ (CARD(4), '(18X,I7)') ANALPD
      READ (CARD(5), '(18X,I10)') FOREPD
      READ (CARD(6), '(28X,I7)') IS
      READ (CARD(7), '(14X,A14,4X,A14)') CABDT, CAEDT
      IASS=0
      IF (CARD(8)(14:14).EQ.'Y') IASS=1
      NF=0
      IF (CARD(9)(14:14).EQ.'Y') NF=1
      ISTAT(1)=0
      IF (CARD(10)(16:16).EQ.'F') ISTAT(1)=1
      ISTAT(2)=0
      IF (CARD(11)(16:16).EQ.'F') ISTAT(2)=1
      ISTAT(3)=0
      IF (CARD(12)(16:16).EQ.'F') ISTAT(3)=1
                  
      IF(CARD(13).NE.' ') THEN
        READ (CARD(13), '(40X,A14)') CRSDT
      ELSE
        CRSDT = '00000000000000'
      ENDIF

      IF(CARD(14).NE.' ') THEN
        READ (CARD(14), '(36X,A14)') CLSDT
      ELSE
        CLSDT = CEDT
      ENDIF

      IF(CARD(15).NE.' ') THEN
        READ (CARD(15), '(39X,A14)') CCURA
      ELSE
        CCURA = CBDT
      ENDIF

      LRSTPARAL=.FALSE.
      IF (CARD(16)(20:20).EQ.'Y') LRSTPARAL=.TRUE.

      READ (CARD(17), '(27X,I10)') NPROC_RST 

      IF (LHOOK) CALL DR_HOOK('READSTA',1,ZHOOK_HANDLE)
      RETURN                                                            

!*    3.0  ERROR EXIT
!          ----------

 3000 CONTINUE
      WRITE(IU06,*) ' *************************************************'
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *        FATAL ERROR IN SUB. READSTA            *'
      WRITE(IU06,*) ' *        ===========================            *'
      WRITE(IU06,*) ' * READ ERROR ON UNIT  IUIN = ', IUIN              
      WRITE(IU06,*) ' * WAMINFO FILE.                                 *'
      WRITE(IU06,*) ' * RECORD NUMBER IS I = ', I                       
      WRITE(IU06,*) ' * WAMINFO READ BEFORE ERROR IS:                 *'
      WRITE(IU06,*) ' *                                               *'
      DO L=1,I
        WRITE(IU06,'(5X,A72)') CARD(L)
      ENDDO
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *  PROGRAMM ABORTS  PROGRAMM ABORTS             *'
      WRITE(IU06,*) ' *                                               *'
      WRITE(IU06,*) ' *************************************************'
      CALL ABORT1

      IF (LHOOK) CALL DR_HOOK('READSTA',1,ZHOOK_HANDLE)

      END SUBROUTINE READSTA
