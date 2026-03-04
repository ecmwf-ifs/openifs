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

SUBROUTINE FPINT4X(KASLB1,KFPROW,KFIELDS,KGPST,KGPEND,KFPROMA,KFLDBUF, &
 & KL0,PWXX,KPSL,LDMASK,PBUF,PROW,PUNDEF)

!**** *FPINT4X*  - Interpolate fields with missing values using 4 points

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO FRANCE* 
!      Original : 01-10-2014

!     Description.
!     ------------
!     The meaning of arguments are the same as in FPINT4, except that
!     PWXX has an extra dimension of 16 corresponding the number of 
!     possibilities we have here; as we use 4 points, it is possible 
!     to have either the following combinations of valid points :
!     1,2,3,4     ,2,3,4
!     1,2,3,      ,2,3, 
!     1,2, ,4     ,2, ,4
!     1, ,3,4     , ,3,4
!     1,2, ,      ,2, , 
!     1, ,3,      , ,3, 
!     1, , ,4     , , ,4
!     1, , ,      , , , 
!     KPSL gives the connection between validity mask and the index in PWXX
!     KFPROW   : number of raw adresses needed for horizontal interpolations

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KASLB1
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROW
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS
INTEGER(KIND=JPIM),INTENT(IN)    :: KFPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDBUF
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPST
INTEGER(KIND=JPIM),INTENT(IN)    :: KGPEND
INTEGER(KIND=JPIM),INTENT(IN)    :: KL0(KFPROMA,KFPROW)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWXX(KFPROMA,4,16)
INTEGER(KIND=JPIM),INTENT(IN)    :: KPSL (0:1,0:1,0:1,0:1)
LOGICAL           ,INTENT(IN)    :: LDMASK(KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBUF(KASLB1*KFLDBUF)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PROW(KFPROMA,KFIELDS)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PUNDEF

INTEGER(KIND=JPIM) :: IADD(KFPROMA)
INTEGER(KIND=JPIM) :: IADDFLD

INTEGER(KIND=JPIM) :: JF, JI
INTEGER(KIND=JPIM) :: I1, I2, I3, I4, IPSL
REAL(KIND=JPRB)    :: Z1, Z2, Z3, Z4
LOGICAL :: LLINTP
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('FPINT4X',0,ZHOOK_HANDLE)

IF (KFPROW < 3) CALL ABOR1('FPINT4X : KFPROW TOO SMALL')

IADDFLD=0
DO JF = 1, KFIELDS

  IF (.NOT.LDMASK(JF)) THEN

    IADD(KGPST:KGPEND)=KASLB1*IADDFLD

    DO JI=KGPST,KGPEND

! Check which points are valid

      Z1 = PBUF(IADD(JI)+KL0(JI,2)+1); I1 = 0; IF (Z1 /= PUNDEF) I1 = 1
      Z2 = PBUF(IADD(JI)+KL0(JI,2)+2); I2 = 0; IF (Z2 /= PUNDEF) I2 = 1
      Z3 = PBUF(IADD(JI)+KL0(JI,3)+1); I3 = 0; IF (Z3 /= PUNDEF) I3 = 1
      Z4 = PBUF(IADD(JI)+KL0(JI,3)+2); I4 = 0; IF (Z4 /= PUNDEF) I4 = 1

! Index in weight array

      IPSL = KPSL (I1, I2, I3, I4)

      LLINTP = ANY (PWXX(JI,1:4,IPSL) > 0._JPRB)

      IF (LLINTP) THEN
        PROW(JI,JF)=&
         & ( PWXX(JI, 1, IPSL)*Z1 &
         & + PWXX(JI, 2, IPSL)*Z2 &
         & + PWXX(JI, 3, IPSL)*Z3 &
         & + PWXX(JI, 4, IPSL)*Z4 )  
      ELSE
        PROW(JI,JF)=PUNDEF
      ENDIF

    ENDDO

  ENDIF

  IADDFLD=IADDFLD+1
  IF (IADDFLD > KFLDBUF) THEN
    CALL ABOR1('FPINT4X : INTERNAL ERROR IADDFLD')
  ENDIF

ENDDO

IF (LHOOK) CALL DR_HOOK('FPINT4X',1,ZHOOK_HANDLE)

END SUBROUTINE FPINT4X

