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

SUBROUTINE PREPACKA(YDFACTX, KDOM,YDFLDSC,KFLDS,KMASK,YDCPDSC, &
                  & KSTART,KLEN,YDFPOPH,PFIELD,PVALCO)

!**** *PREPACKA*  - preliminar packing of gridpoint or spectral data 
!                   before writing out to file

!     Purpose.
!     --------
!        To pack a chunk of global fields (either spectral or gridpoint) 
!        prior to writing out to a file Arpege/Aladin
!        WARNING : spectral data is expected to be ordered like in FA file

!**   Interface.
!     ----------
!        *CALL* *PREPACKA*

!        Explicit arguments :
!        --------------------

!            KDOM    : number of subdomains
!            PFIELD   : global fields packed (input/output)
!            KFLDS   : number of fields to pack
!            KMASK   : binary mask ; =0 to skip the field (no write out) ;
!                      else =1
!            KNOMA   : length of field name to write afterward
!            KSTART  : start adress of fields in PFIELD
!            KLEN    : actual size of fields in PFIELD
!            PVALCO  : If present, then this array will hold compressed data,
!                      and PFIELD will not be overwritten

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!           FA routines

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-12-18 from wrpxfa

!     Modifications.
!     --------------
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      D. Paradis & R. El Khatib : 04-07-22 GRIBEX
!      R. El Ouaraini : 03-Oct-06 GRIBEX: INGRIG is computed from KGRIB(NFPGRIB) instead of INGRIB(NVGRIB)
!      P.Marguinaud   : 26-Apr-2012 : Adapt for IO server (compression now performed
!                                     in PREPACKA1)
!      P.Marguinaud   : 06-Jun-2012 Set date in FA units used for compression
!      R. El Khatib : 20-Aug-2012 PSPECG + KRESOL and replace model dimensionning
!      P.Marguinaud   : 11-Sep-2012 Add extra argument PVALCO
!      R. El Khatib : 08-Feb-2013 Bugfix
!      P.Marguinaud   : 10-Oct-2013 Re-write using field descriptors and FACTX
!      R. El Khatib : 09-Aug-2013 Bugfix (no call to ellips in not LDCOSP)
!      P.Marguinaud   : 10-Oct-2014 Cleaning
!      P.Marguinaud   : 04-Oct-2016 Port to single precision
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE TYPE_FAOPH, ONLY : TFAOPH
USE IOFLDDESC_MOD, ONLY : IOFLDDESC
USE IOCPTDESC_MOD, ONLY : IOCPTDESC
USE FACTX_MOD, ONLY : FACTX, OPENCPFACTX
USE FA_MOD, ONLY : FA_COM, FA_COM_DEFAULT
USE OML_MOD, ONLY : OML_MY_THREAD

IMPLICIT NONE

TYPE(FACTX),       INTENT(INOUT),           TARGET :: YDFACTX (:)
INTEGER(KIND=JPIM),INTENT(IN)                      :: KDOM 
INTEGER(KIND=JPIM),INTENT(IN)                      :: KFLDS 
TYPE (IOFLDDESC)  ,INTENT(IN)                      :: YDFLDSC(KFLDS)
INTEGER(KIND=JPIM),INTENT(IN)                      :: KMASK(KFLDS,KDOM) 
TYPE (IOCPTDESC)  ,INTENT(INOUT)                   :: YDCPDSC(KFLDS,KDOM)
INTEGER(KIND=JPIM),INTENT(IN)                      :: KSTART(KDOM) 
INTEGER(KIND=JPIM),INTENT(IN)                      :: KLEN(KDOM) 
TYPE (TFAOPH)     ,INTENT(IN)                      :: YDFPOPH(KDOM)
REAL(KIND=JPRB)   ,INTENT(IN),              TARGET :: PFIELD(:,:) 
REAL (KIND=JPRB),  INTENT(INOUT),           TARGET :: PVALCO (:,:)

#include "abor1.intfb.h"
#include "prepacka1_mt.intfb.h"

LOGICAL :: LLOPENMP

INTEGER(KIND=JPIM) :: JFLD, JDOM, JDF
INTEGER(KIND=JPIM) :: IGRIB, ITID
INTEGER(KIND=JPIM) :: IMINCO, IMAXCO, IMINFL, IMAXFL
TYPE (FACTX),  POINTER :: YLFACTX
TYPE (FA_COM), POINTER :: YLFA

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PREPACKA',0,ZHOOK_HANDLE)
ASSOCIATE(NHEADMAX=>YDFPOPH%NHEADMAX,CFPCA=>YDFPOPH%CFPCA)
!      -----------------------------------------------------------

IF (ANY(YDFLDSC%LSPEC) .AND. (KDOM > 1)) THEN
  CALL ABOR1('PREPACKA : NOT SUITABLE FOR MULTI-SPECTRAL PACKING')
ENDIF

!*    1.    SETUP PACKING ENVIRONMENT WITHOUT OPENING FILE
!           ----------------------------------------------

DO JDOM = 1, KDOM
  IF (ANY (KMASK (:,JDOM) /= 0)) THEN
    CALL OPENCPFACTX (YDFACTX (JDOM),CFPCA(JDOM))
  ENDIF
ENDDO

LLOPENMP=ALL (YDFACTX%YIOOPTS%LFACPT_OPENMP)

!$OMP  PARALLEL PRIVATE (JDF, JDOM, JFLD, ITID, YLFA, IGRIB,        &
!$OMP&                   IMINCO, IMAXCO, IMINFL, IMAXFL, YLFACTX),  &
!$OMP& IF (LLOPENMP)

ITID = OML_MY_THREAD ()

!$OMP DO SCHEDULE(DYNAMIC, 1)
DO JDF = 1, KDOM * KFLDS

  JDOM = 1 + (JDF-1) / KFLDS
  JFLD = 1 + MODULO (JDF-1, KFLDS)

  IF (KMASK (JFLD,JDOM) /= 1) CYCLE

  YLFACTX => YDFACTX (JDOM)

  IF (LLOPENMP) THEN
    YLFA => YLFACTX%YFABYTID (ITID)
  ELSE
    YLFA => FA_COM_DEFAULT
  ENDIF

  IGRIB=YDFLDSC (JFLD)%NGRIBL

  IMINCO = KSTART(JDOM)+SUM (NHEADMAX(1:JDOM-1))
  IMAXCO = KSTART(JDOM)+SUM (NHEADMAX(1:JDOM))+KLEN(JDOM)-1
  IMINFL = KSTART(JDOM)
  IMAXFL = KSTART(JDOM)+KLEN(JDOM)-1

  CALL PREPACKA1_MT (YLFA, YLFACTX%IUNITCP, YDFLDSC (JFLD), PFIELD (IMINFL:IMAXFL,JFLD),   &
                   & YDCPDSC (JFLD,JDOM), PVALCO (IMINCO:IMAXCO,JFLD))

ENDDO
!$OMP END DO
  
!$OMP END PARALLEL

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('PREPACKA',1,ZHOOK_HANDLE)

END SUBROUTINE PREPACKA
