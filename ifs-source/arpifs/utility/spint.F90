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

SUBROUTINE SPINT(PSPEC,KDIM,KSMAXI,KSMAXO)

!**** *SPINT* - simple spectral interpolator

!     Purpose. "Interpolate" from low to high resolution by filling the
!     --------  rest of the spectrum with zeros

!**   Interface.
!     ----------
!        *CALL* *SPINT(PSPEC,KDIM,KSMAXI,KSMAXO)*

!        Explicit arguments : PSPEC - spectral array
!        -------------------- KDIM  - dimension of PSPEC
!                             KSMAXI - input resolution
!                             KSMAXO - output resolution

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------
!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        M.Hamrud  *ECMWF* (copied from expernal program spint)

!     Modifications.
!     --------------
!        Original : 98-08-01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)    :: KDIM 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PSPEC(KDIM) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAXI 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSMAXO 
REAL(KIND=JPRB) :: ZSPECI((KSMAXI+1)*(KSMAXI+2))
INTEGER(KIND=JPIM) :: IASM0I(0:KSMAXI)
INTEGER(KIND=JPIM) :: IASM0O(0:KSMAXO)

INTEGER(KIND=JPIM) :: IPOS, ISEI, ISEO, ISMAX, ISPEC2I, ISPEC2O,&
 & ISPECI, ISPECO, JI, JM, JN  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*    1.    FILL WITH ZEROS
!           ---------------

IF (LHOOK) CALL DR_HOOK('SPINT',0,ZHOOK_HANDLE)
ISPECI=(KSMAXI+1)*(KSMAXI+2)/2
ISPEC2I=ISPECI*2
ISPECO=(KSMAXO+1)*(KSMAXO+2)/2
ISPEC2O=ISPECO*2
IPOS = 1
DO JM=0,KSMAXI
  IASM0I(JM) = IPOS
  IPOS = IPOS+(KSMAXI-JM+1)*2
ENDDO
IPOS = 1
DO JM=0,KSMAXO
  IASM0O(JM) = IPOS
  IPOS = IPOS+(KSMAXO-JM+1)*2
ENDDO
ZSPECI(:)=PSPEC(1:ISPEC2I)
PSPEC(:)=0.0_JPRB
ISMAX=MIN(KSMAXI,KSMAXO)
DO JM=0,ISMAX
  DO JI=1,2
    DO JN=JM,ISMAX
      ISEI=IASM0I(JM)+(JN-JM)*2+JI-1
      ISEO=IASM0O(JM)+(JN-JM)*2+JI-1
      PSPEC(ISEO)=ZSPECI(ISEI)
    ENDDO
  ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('SPINT',1,ZHOOK_HANDLE)
END SUBROUTINE SPINT
