! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

SUBROUTINE FPINT12(KASLB1,KFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF,KDER, &
 & KDMP,LDML,LDNRST,LDMONO,KBOX,KL0,PWXX,PWXY,LDMASK,PBUF,PROW,PUNDEF)

!**** *FPINT12*  - FULL-POS horizontal interpolator.

!     PURPOSE.
!     --------

!        Performs interpolations for each point of a row and for a set 
!        of fields.
!        Numbering of the points (I is the interpolation point):
!                     13       5       6      14

!                      7       1       2       8
!                                 (I)
!                      9       3       4      10

!                     15      11      12      16
!       Points number 1 to 12 are used if bicubic interpolations.
!       All points are used if multi-linear + bicubic interpolations.

!**   INTERFACE.
!     ----------
!       *CALL* *FPINT12*

!        EXPLICIT ARGUMENTS
!        --------------------

!        INPUT:
!         KASLB1   : size of interpolation buffer (core + halo)
!         KFPROW   : number of raw adresses needed for horizontal interpolations
!         KFIELDS  : number of fields in row.
!         KGPST    : first output point in row.
!         KGPEND   : last output point in row.
!         KFPROMA  : length of the output row.
!         KFLDBUF  : number of fields in input buffer
!         KDER     : number of "boxes" for each field (used for
!                    "horizontal derivative" in stretched mode). 
!         KDMP     : subdomain number for each field in input buffer (if relevant)
!         LDML     : multilinear interpolations if T.
!         LDNRST   : .TRUE. to replace interpolations by the search of the nearest point
!         LDMONO   : .TRUE. if monotonic interpolations
!         KBOX     : index of the subdomain each point belongs.
!                    Used for for pp. of stretched filtred derivatives.
!         KL0      : indexes of the western points, from north to south :
!                    4 points (nr 13, 7, 9, 15)
!         PWXX     : weights for bi-linear interpolations.
!         PWXY     : additional weights array for multi-linear interpolations.
!         LDMASK   : logical mask : .TRUE. => no computation
!         PBUF     : buffer which contain the fields to interpolate.
!         PUNDEF   : value for undefined points

!        OUTPUT:
!         PROW     : interpolated fields.

!        IMPLICIT ARGUMENTS
!        --------------------
!          None

!     METHOD.
!     -------
!        Horizontal interpolations : with 12 points.
!        Note that the adress of the input field may depend of the 
!        output point (in the case  field is considered as an 
!        "horizontal derivative".
!        Anyway, the computation is done only if the logical mask is .TRUE.

!     EXTERNALS.
!     ----------
!      None

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        See documentation about FULL-POS.

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 94-04-08

!     MODIFICATIONS.
!     --------------
!      R. El Khatib 01-10-03 Fix monotonic interpolations for land-sea dependant fields
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib  20-May-2005 YFPSTRUCT%NSLWIDE moved to YOMWFPDS
!      K. Yessad (oct 2010): multi-linear interpolations.
!     G.Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!     P.Marguinaud (Oct 2014): Accept undefined values in fields to be
!     interpolated, allow a wider halo
!      R. El Khatib : 18-Mar-2016 Vectorize again
!      R. El Khatib 27-Jul-2016 Cleaning
!     -----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

!     -----------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROW
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDBUF
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPST
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KDER(KFIELDS)
INTEGER(KIND=JPIM),INTENT(IN)    :: KDMP(KFLDBUF)
LOGICAL           ,INTENT(IN)    :: LDML
LOGICAL           ,INTENT(IN)    :: LDNRST(KFIELDS)
LOGICAL           ,INTENT(IN)    :: LDMONO(KFIELDS)
INTEGER(KIND=JPIM),INTENT(IN)    :: KBOX(KGPEND)
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KFPROMA,KFPROW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWXX(KFPROMA,12)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWXY(KFPROMA,5:16)
LOGICAL           ,INTENT(IN)    :: LDMASK(KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUF(KASLB1*KFLDBUF)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROW(KFPROMA,KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)  ,OPTIONAL :: PUNDEF

!     -----------------------------------------------------------------------
!     IADD : field index in input buffer.
!     INRST: location in input buffer of the points with highest weight 

INTEGER(KIND=JPIM) :: IADD(KFPROMA), INRST(KFPROMA,2)
INTEGER(KIND=JPIM) :: IADDFLD, IMAX0(KFPROMA,2)

INTEGER(KIND=JPIM) :: JDER, JF, JI, JW, IFLDBUF

REAL(KIND=JPRB) :: ZE, ZMX(KFPROMA), ZMN(KFPROMA)
LOGICAL :: LLIZNA
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------------

#include "abor1.intfb.h"

!     -----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPINT12',0,ZHOOK_HANDLE)
!     -----------------------------------------------------------------------

IF (KFPROW < 4) CALL ABOR1('FPINT12 : KFPROW TOO SMALL')

LLIZNA = PRESENT (PUNDEF)

!*    2. FIND NEAREST POINT AND A LAND-SEA COMPATIBLE POINT TO START
!        -----------------------------------------------------------

ZE=EPSILON(1.0_JPRB)*10._JPRB

IF (ANY(LDNRST).OR.ANY(LDMONO)) THEN
  ! Search the nearest land-sea compatible point
  DO JI=KGPST,KGPEND
    IMAX0(JI,1:1)=MAXLOC(PWXX(JI,1:12))
  ENDDO
  ! Search a land-sea compatible point among the 4 nearest ones first,
  ! then among all
  IMAX0(:,2)=0
  DO JW=1,12
    DO JI=KGPST,KGPEND
      IF (PWXX(JI,JW) >= ZE .AND. IMAX0(JI,2) == 0) IMAX0(JI,2)=JW
    ENDDO
  ENDDO
  DO JW=1,2
    DO JI=KGPST,KGPEND
      SELECT CASE (IMAX0(JI,JW))
      CASE (1)
        INRST(JI,JW)=KL0(JI,2)+1
      CASE (2)
        INRST(JI,JW)=KL0(JI,2)+2
      CASE (3)
        INRST(JI,JW)=KL0(JI,3)+1 
      CASE (4)
        INRST(JI,JW)=KL0(JI,3)+2
      CASE (5)
        INRST(JI,JW)=KL0(JI,1)+1
      CASE (6)
        INRST(JI,JW)=KL0(JI,1)+2
      CASE (7)
        INRST(JI,JW)=KL0(JI,2) 
      CASE (8)
        INRST(JI,JW)=KL0(JI,2)+3
      CASE (9)
        INRST(JI,JW)=KL0(JI,3)
      CASE (10)
        INRST(JI,JW)=KL0(JI,3)+3
      CASE (11)
        INRST(JI,JW)=KL0(JI,4)+1
      CASE (12)
        INRST(JI,JW)=KL0(JI,4)+2 
      END SELECT
    ENDDO
  ENDDO
ENDIF

!     -----------------------------------------------------------------------

!*    3. LOOP ON ALL FIELDS
!        ------------------

IADDFLD=0
DO JF = 1, KFIELDS

  IF (.NOT.LDMASK(JF)) THEN

    !*  3.1 COMPUTE INPUT FIELD ADRESS

    IADD(KGPST:KGPEND)=KASLB1*IADDFLD
    IF (KDER(JF) > 1) THEN
      ! Modify gp adresses to point onto the right subdomain
      DO JDER=1,KDER(JF)
        IFLDBUF=IADDFLD+JDER
        DO JI=KGPST,KGPEND
          IF (KBOX(JI) == KDMP(IFLDBUF)) THEN
            IADD(JI)=KASLB1*(IADDFLD+JDER-1)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    !*  3.2 COMPUTE OUTPUT POINT

    IF (.NOT.LDNRST(JF)) THEN
      IF (LLIZNA) THEN
        DO JI=KGPST,KGPEND
          IF (ANY (PWXX(JI,1:12) > 0._JPRB)) THEN
            PROW(JI,JF)= &
             & ( PWXX(JI, 1)*PBUF(IADD(JI)+KL0(JI,2)+1)&
             & + PWXX(JI, 2)*PBUF(IADD(JI)+KL0(JI,2)+2)&
             & + PWXX(JI, 3)*PBUF(IADD(JI)+KL0(JI,3)+1)&
             & + PWXX(JI, 4)*PBUF(IADD(JI)+KL0(JI,3)+2)&
             & + PWXX(JI, 5)*PBUF(IADD(JI)+KL0(JI,1)+1)&
             & + PWXX(JI, 6)*PBUF(IADD(JI)+KL0(JI,1)+2)&
             & + PWXX(JI, 7)*PBUF(IADD(JI)+KL0(JI,2)  )&
             & + PWXX(JI, 8)*PBUF(IADD(JI)+KL0(JI,2)+3)&
             & + PWXX(JI, 9)*PBUF(IADD(JI)+KL0(JI,3)  )&
             & + PWXX(JI,10)*PBUF(IADD(JI)+KL0(JI,3)+3)&
             & + PWXX(JI,11)*PBUF(IADD(JI)+KL0(JI,4)+1)&
             & + PWXX(JI,12)*PBUF(IADD(JI)+KL0(JI,4)+2) )  
          ELSE
            PROW(JI,JF)=PUNDEF
          ENDIF
        ENDDO
      ELSE
        DO JI=KGPST,KGPEND
          PROW(JI,JF)= &
           & ( PWXX(JI, 1)*PBUF(IADD(JI)+KL0(JI,2)+1)&
           & + PWXX(JI, 2)*PBUF(IADD(JI)+KL0(JI,2)+2)&
           & + PWXX(JI, 3)*PBUF(IADD(JI)+KL0(JI,3)+1)&
           & + PWXX(JI, 4)*PBUF(IADD(JI)+KL0(JI,3)+2)&
           & + PWXX(JI, 5)*PBUF(IADD(JI)+KL0(JI,1)+1)&
           & + PWXX(JI, 6)*PBUF(IADD(JI)+KL0(JI,1)+2)&
           & + PWXX(JI, 7)*PBUF(IADD(JI)+KL0(JI,2)  )&
           & + PWXX(JI, 8)*PBUF(IADD(JI)+KL0(JI,2)+3)&
           & + PWXX(JI, 9)*PBUF(IADD(JI)+KL0(JI,3)  )&
           & + PWXX(JI,10)*PBUF(IADD(JI)+KL0(JI,3)+3)&
           & + PWXX(JI,11)*PBUF(IADD(JI)+KL0(JI,4)+1)&
           & + PWXX(JI,12)*PBUF(IADD(JI)+KL0(JI,4)+2) )  
        ENDDO
      ENDIF
      IF (LDML) THEN
        DO JI=KGPST,KGPEND
          PROW(JI,JF)=PROW(JI,JF) &
           & + PWXY(JI, 5)*PBUF(IADD(JI)+KL0(JI,1)+1)&
           & + PWXY(JI, 6)*PBUF(IADD(JI)+KL0(JI,1)+2)&
           & + PWXY(JI, 7)*PBUF(IADD(JI)+KL0(JI,2)  )&
           & + PWXY(JI, 8)*PBUF(IADD(JI)+KL0(JI,2)+3)&
           & + PWXY(JI, 9)*PBUF(IADD(JI)+KL0(JI,3)  )&
           & + PWXY(JI,10)*PBUF(IADD(JI)+KL0(JI,3)+3)&
           & + PWXY(JI,11)*PBUF(IADD(JI)+KL0(JI,4)+1)&
           & + PWXY(JI,12)*PBUF(IADD(JI)+KL0(JI,4)+2)&
           & + PWXY(JI,13)*PBUF(IADD(JI)+KL0(JI,1)  )&
           & + PWXY(JI,14)*PBUF(IADD(JI)+KL0(JI,1)+3)&
           & + PWXY(JI,15)*PBUF(IADD(JI)+KL0(JI,4)  )&
           & + PWXY(JI,16)*PBUF(IADD(JI)+KL0(JI,4)+3)
        ENDDO
      ENDIF
      IF (LDMONO(JF)) THEN
        DO JI=KGPST,KGPEND
          ZMX(JI)=PBUF(IADD(JI)+INRST(JI,2))
          ZMN(JI)=PBUF(IADD(JI)+INRST(JI,2))
          IF (PWXX(JI, 1) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,2)+1))
          IF (PWXX(JI, 2) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,2)+2))
          IF (PWXX(JI, 3) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,3)+1))
          IF (PWXX(JI, 4) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,3)+2))
          IF (PWXX(JI, 1) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,2)+1))
          IF (PWXX(JI, 2) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,2)+2))
          IF (PWXX(JI, 3) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,3)+1))
          IF (PWXX(JI, 4) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,3)+2))
          IF (IMAX0(JI,2) > 4) THEN
            IF (PWXX(JI, 5) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,1)+1))
            IF (PWXX(JI, 6) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,1)+2))
            IF (PWXX(JI, 7) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,2)  ))
            IF (PWXX(JI, 8) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,2)+3))
            IF (PWXX(JI, 9) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,3)  ))
            IF (PWXX(JI,10) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,3)+3))
            IF (PWXX(JI,11) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,4)+1))
            IF (PWXX(JI,12) >= ZE) ZMX(JI)=MAX(ZMX(JI),PBUF(IADD(JI)+KL0(JI,4)+2))
            IF (PWXX(JI, 5) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,1)+1))
            IF (PWXX(JI, 6) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,1)+2))
            IF (PWXX(JI, 7) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,2)  ))
            IF (PWXX(JI, 8) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,2)+3))
            IF (PWXX(JI, 9) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,3)  ))
            IF (PWXX(JI,10) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,3)+3))
            IF (PWXX(JI,11) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,4)+1))
            IF (PWXX(JI,12) >= ZE) ZMN(JI)=MIN(ZMN(JI),PBUF(IADD(JI)+KL0(JI,4)+2))
          ENDIF
          PROW(JI,JF) = MIN(ZMX(JI),MAX(ZMN(JI),PROW(JI,JF)))
        ENDDO
      ENDIF
    ELSE
      DO JI=KGPST,KGPEND
        PROW(JI,JF) = PBUF(IADD(JI)+INRST(JI,1))
      ENDDO
    ENDIF

  ENDIF

  IADDFLD=IADDFLD+KDER(JF)
  IF (IADDFLD > KFLDBUF) THEN
    CALL ABOR1('FPINT12 : INTERNAL ERROR IADDFLD')
  ENDIF

ENDDO

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FPINT12',1,ZHOOK_HANDLE)
END SUBROUTINE FPINT12
