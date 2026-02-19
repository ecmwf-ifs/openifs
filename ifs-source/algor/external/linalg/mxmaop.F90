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

SUBROUTINE MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!**** *MXMAOP - Optimize call to SGEMMX

!     Purpose.
!     --------
!        Make sure SGEMMX is called in a way to insure maximum optimization.
!        Does matrix product PC=PA*PB

!**   Interface.
!     ----------
!        CALL MXMAOP(PA,KA,KAD,PB,KB,KBD,PC,KC,KCA,KAR,KAC,KBC)

!        Explicit arguments : See SGEMMX documentaion.
!        --------------------
!         PA     - input matrix PA                                     (in)
!         KA     - memory jump between two lines in PA (generally 1)   (in)
!         KAD    - memory jump between two columns in PA               (in)
!         PB     - input matrix PB                                     (in)
!         KB     - memory jump between two lines in PB (generally 1)   (in)
!         KBD    - memory jump between two columns in PB               (in)
!         PC     - output matrix PC                                    (out)
!         KC     - memory jump between two lines in PC (generally 1)   (in)
!         KCA    - memory jump between two columns in PC               (in)
!         KAR    - number of useful lines of PA                        (in)
!         KAC    - number of useful columns of PA (and lines of PB)    (in)
!         KBC    - number of useful columns of PB                      (in)

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.   SGEMMX in Cray library.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 88-01-28
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        K. Yessad (Apr 2009): add comments.
!        J. Hague  (Oct 2012): Parallelise call to SGEMMX
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
USE OML_MOD   ,ONLY : OML_MAX_THREADS, OML_IN_PARALLEL

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPRB)   ,INTENT(IN)    :: PA(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAD 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PB(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBD 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PC(*) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KCA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAR 
INTEGER(KIND=JPIM),INTENT(IN)    :: KAC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KBC 

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JJJ,JLEN,JT,JCHUNK

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MXMAOP',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       1.       PERFORM LEGENDRE TRANFORM.
!                 --------------------------

IF(OML_IN_PARALLEL()) THEN

  IF (KAR >= KBC) THEN
    CALL SGEMMX(KAR,KBC,KAC,1.0_JPRB,PA,KA,KAD,PB,KB,KBD,0.0_JPRB,PC,KC,KCA)
  ELSE
    CALL SGEMMX(KBC,KAR,KAC,1.0_JPRB,PB,KBD,KB,PA,KAD,KA,0.0_JPRB,PC,KCA,KC)
  ENDIF

ELSE

  JT=OML_MAX_THREADS()

  IF (KAR >= KBC) THEN
    JCHUNK=(KAR-1)/JT+1
!$OMP PARALLEL DO PRIVATE(JJJ,JLEN)
    DO JJJ=1,KAR,JCHUNK
      JLEN=MIN(JCHUNK,KAR-JJJ+1)
       IF(JLEN>0)  CALL SGEMMX(JLEN,KBC,KAC,1.0_JPRB,PA((JJJ-1)*KA+1),KA,KAD,PB,KB,KBD,0.0_JPRB,PC((JJJ-1)*KC+1),KC,KCA)
    ENDDO
!$OMP END PARALLEL DO 

  ELSE

    JCHUNK=(KBC-1)/JT+1
!$OMP PARALLEL DO PRIVATE(JJJ,JLEN)
    DO JJJ=1,KBC,JCHUNK
      JLEN=MIN(JCHUNK,KBC-JJJ+1)
      IF(JLEN>0) CALL SGEMMX(JLEN,KAR,KAC,1.0_JPRB,PB((JJJ-1)*KBD+1),KBD,KB,PA,KAD,KA,0.0_JPRB,PC((JJJ-1)*KCA+1),KCA,KC)
    ENDDO
!$OMP END PARALLEL DO 
  ENDIF

ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('MXMAOP',1,ZHOOK_HANDLE)
END SUBROUTINE MXMAOP

