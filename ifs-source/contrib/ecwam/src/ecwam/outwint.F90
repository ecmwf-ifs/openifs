! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTWINT(BOUT)

! ----------------------------------------------------------------------

!**** *OUTWINT* -

!*    PURPOSE.
!     --------


!**   INTERFACE.
!     ----------

!      *BOUT*    - OUTPUT PARAMETERS BUFFER

!     METHOD.
!     -------

! ----------------------------------------------------------------------
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWDRVTYPE  , ONLY : WVGRIDGLO

      USE YOWCOUP  , ONLY : LWCOU, LIFS_IO_SERV_ENABLED, OUTINT_IO_SERV
      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU
      USE YOWCOUT  , ONLY : JPPFLAG ,NIPRMOUT , NINFOBOUT,              &
     &                      INFOBOUT,LWAM_USE_IO_SERV
      USE YOWGRID  , ONLY : NPROMA_WAM, NCHNK
      USE YOWSTAT  , ONLY : CDATEA, CDATEF, CDTPRO, MARSTYPE
      USE YOWTEST  , ONLY : IU06
      USE YOMHOOK   ,ONLY : LHOOK, DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include "abort1.intfb.h"
#include "difdate.intfb.h"
#include "incdate.intfb.h"
#include "outint.intfb.h"

      REAL(KIND=JWRB), DIMENSION(NPROMA_WAM, NIPRMOUT, NCHNK), INTENT(IN) :: BOUT


      INTEGER(KIND=JWIM) :: IFCST, INHOUR, ISHIFT, IY, IM, ID, IH, IMN, ISS

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=14) :: CDATE, CDATED, CDATE1, CDATE2

! ----------------------------------------------------------------------
      IF (LHOOK) CALL DR_HOOK('OUTWINT',0,ZHOOK_HANDLE)

!     DETERMINE THE FORECAST STEP INFORMATION
      IF (MARSTYPE == '4v') THEN
!*      THIS IS A 4V CASE
        CDATE=CDATEA
        CDATED=CDATEA
        READ(CDATEA,'(I4,5I2)')IY,IM,ID,IH,IMN,ISS
        CALL DIFDATE(CDATEA, CDTPRO, IFCST)
        IFCST=IFCST/3600
      ELSEIF (CDTPRO <= CDATEF) THEN
!*      THIS IS AN ANALYSIS DATE
        CDATE=CDTPRO
        CDATED=CDTPRO
        READ(CDTPRO,'(I4,5I2)') IY,IM,ID,IH,IMN,ISS
        IFCST = 0
      ELSE
!*      THIS IS A  FORECAST DATE
        CDATE=CDATEF
        CDATED=CDTPRO
        READ(CDATEF,'(I4,5I2)')IY,IM,ID,IH,IMN,ISS
!       FIND THE STEP IN HOURS ! 
        CDATE1=CDATEF
        CDATE2=CDATE1(1:8)//'010000'
!       MORE THAN A YEAR IN SECONDS 366*24*3600=31622400
        CALL INCDATE(CDATE2, 31622400) 
        CDATE2=CDATE2(1:8)//'010000'
        INHOUR=0
!       FIND THE NUMBER OF FULL YEARS AND CONVERT THEM IN HOURS
        DO WHILE (CDATE2 < CDTPRO)
          CALL DIFDATE(CDATE1, CDATE2, ISHIFT)
          INHOUR=INHOUR+ISHIFT/3600
          CDATE1=CDATE2
          CDATE2=CDATE1(1:8)//'010000'
          CALL INCDATE(CDATE2,31622400) 
          CDATE2=CDATE2(1:8)//'010000'
        ENDDO
        CALL DIFDATE(CDATE1, CDTPRO, IFCST)
        IF(MOD(IFCST,3600) == 0 ) THEN
          IFCST=IFCST/3600
        ELSE
          WRITE(IU06,*) ' -----------------------------------------'
          WRITE(IU06,*) ' ERROR in routine OUTINT :'
          WRITE(IU06,*) ' forecast step must be  multiple of hours!'
          WRITE(IU06,*) ' IFCST =', IFCST
          WRITE(IU06,*) ' CDATEF=', CDATEF 
          WRITE(IU06,*) ' CDATE1=', CDATE1 
          WRITE(IU06,*) ' CDTPRO=', CDTPRO 
          WRITE(IU06,*) ' INHOUR=', INHOUR 
          WRITE(IU06,*) ' -----------------------------------------'
          CALL ABORT1
        ENDIF
        IFCST=INHOUR+IFCST
      ENDIF

      ! Use IFS IO server?
      IF (LWCOU .AND. LIFS_IO_SERV_ENABLED .AND. LWAM_USE_IO_SERV) THEN
          CALL OUTINT_IO_SERV(NIPRMOUT, BOUT, INFOBOUT, MARSTYPE, CDATE, IFCST)
      ELSE
          CALL OUTINT(CDATE, CDATED, IFCST, BOUT)
      ENDIF

      WRITE(IU06,*) ' '
      WRITE(IU06,*) '  INTEGRATED PARAMETER DISPOSED AT... CDTPRO = ', CDTPRO
      WRITE(IU06,*) '  MARSTYPE = ',MARSTYPE
      WRITE(IU06,*) '  CDATE = ',CDATE
      WRITE(IU06,*) '  IFCST = ',IFCST
      WRITE(IU06,*) ' '
      CALL FLUSH(IU06)

    IF (LHOOK) CALL DR_HOOK('OUTWINT',1,ZHOOK_HANDLE)

    END SUBROUTINE OUTWINT
