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

SUBROUTINE PREPACKA1_MT (YDFA, KUNTIN, YDFLDSC, PFIELD, YDCPDSC, PVALCO)

!**** *PREPACKA1_MT*  - preliminar packing of gridpoint or spectral data
!                    before writing out to file; a single field is processed

!     Author. 
!     ------- 
!      Philippe Marguinaud *METEO-FRANCE*
!      Original : 12-04-2012 (from prepacka.F90)

!     Modifications :
!     P. Marguinaud : 10-10-2013 : Use IOFLDDESC_MOD & IOCPTDESC_MOD
!                                  compress fields using method 4
!     P. Marguinaud : 10-10-2014 : Use FACONO
!     P. Marguinaud : 04-10-2016 : Port to single precision



USE PARKIND1, ONLY : JPIM, JPRB
USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE FA_MOD,   ONLY : FA_COM, JPPRCM
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE IOCPTDESC_MOD, ONLY : IOCPTDESC

IMPLICIT NONE

TYPE (FA_COM),       INTENT (INOUT) :: YDFA
INTEGER (KIND=JPIM), INTENT (IN)    :: KUNTIN
TYPE (IOFLDDESC),    INTENT (IN)    :: YDFLDSC
REAL (KIND=JPRB),    INTENT (IN)    :: PFIELD (:)
TYPE (IOCPTDESC),    INTENT (INOUT) :: YDCPDSC
REAL (KIND=JPRB),    INTENT (OUT)   :: PVALCO (:)

#include "facono_mt.h"

! FAVEUR & FAGOTE arguments
INTEGER (KIND=JPIM) :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
INTEGER (KIND=JPIM) :: INGRIG
INTEGER (KIND=JPIM) :: ILONGD
INTEGER (KIND=JPIM) :: ICOD
INTEGER (KIND=JPIM) :: IBIT, IGRIB
INTEGER (KIND=JPIM) :: JBITS
INTEGER (KIND=JPIM) :: IREP

INTEGER(KIND=JPIM), PARAMETER :: JPMAXBIT = 31

REAL (KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK ('PREPACKA1_MT',0,ZHOOK_HANDLE)

! Prepare for packing 

IGRIB=YDFLDSC%NGRIBL

CALL FAVEUR_MT(YDFA,IREP,KUNTIN,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)

IF (IGRIB == 4) THEN

  CALL FAGOTE_MT(YDFA,IREP,KUNTIN,IGRIB,YDFLDSC%ICPLS,&
               & YDFLDSC%ICPLB,0_JPIM,0_JPIM,0_JPIM)

ELSEIF (IGRIB <= 3) THEN

  JBITS=YDFLDSC%JBITS

! data should be reordered only if packing type is not -1 and not 3
  IF (IGRIB==-1 .OR. IGRIB==3) THEN
    INGRIG=-1
  ELSE    
    INGRIG=0
  ENDIF   

  IF ((JBITS > JPMAXBIT) .AND. (IGRIB <= 3)) THEN
    ICOD=INGRIG
  ELSE    
    ICOD=IGRIB
  ENDIF   

!     Prevent abort in FAGOTE :
  IBIT=MIN(JPMAXBIT,JBITS)
  
  CALL FAGOTE_MT(YDFA,IREP,KUNTIN,ICOD,IBIT,IBIT,ISTRON,IPUILA,IDMOPL)

ELSE
  CALL FAGOTE_MT(YDFA,IREP,KUNTIN,YDFLDSC%NGRIBL,YDFLDSC%JBITS,YDFLDSC%JBITS,ISTRON,IPUILA,IDMOPL)
ENDIF

!*        2.3 Pack

ILONGD = SIZE (PVALCO) / JPPRCM
CALL FACONO_MT (YDFA, IREP, KUNTIN, YDFLDSC%CPREF, YDFLDSC%ILEVG, YDFLDSC%CSUFF, PFIELD, &
             &  YDFLDSC%LSPEC, YDCPDSC%CNOMA, YDCPDSC%ILNOMA, PVALCO, ILONGD,            &
             &  LDUNDF=YDFLDSC%LUNDF, PUNDF=YDFLDSC%XUNDF)
YDCPDSC%ILONGD = ILONGD

! Reset GRIB parameters 

CALL FAGOTE_MT(YDFA,IREP,KUNTIN,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)

IF (LHOOK) CALL DR_HOOK ('PREPACKA1_MT',1,ZHOOK_HANDLE)

END SUBROUTINE PREPACKA1_MT

