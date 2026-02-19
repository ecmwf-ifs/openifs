! (C) Copyright 1996- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE SUWAVEDI_MOD
CONTAINS
SUBROUTINE SUWAVEDI(KSMAX,KTMAX,KPRTRW,KMYSETW,KASM0,KSPOLEGL,KPROCM,&
                    &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS,&
                    &KPTRMS,KALLMS,KDIM0G)

!**** *SUWAVEDI * - Routine to initialize spectral wave distribution

!     Purpose.
!     --------
!           Initialize arrays controlling spectral wave distribution

!**   Interface.
!     ----------
!        *CALL* *SUWAVEDI *

!        Explicit arguments : 
!        --------------------
!           KSMAX    - Spectral truncation limit (input)
!           KTMAX    - Overtruncation for KSMAX (input)
!           KPRTRW   - Number of processors in A-direction (input)
!           KMYSETW  - A-set for present processor (input)
!           KASM0    - Offsets for spectral waves (output)
!           KSPOLEGL - Local version of NSPOLEG (output)
!           KPROCM   - Where a certain spectral wave belongs  (output)
!           KUMPP    - Number of spectral waves on this PE (output)
!           KSPEC    - Local version on NSPEC (output)
!           KSPEC2   - Local version on NSPEC2 (output)
!           KSPEC2MX - Maximum KSPEC2 across PEs (output)
!           KPOSSP   - Global spectral fields partitioning (output)
!           KMYMS    - This PEs spectral zonal wavenumbers (output)
!           KPTRMS   - Pointer to the first wave number of a given a-set (output)
!           KALLMS   - Wave numbers for all wave-set concatenated together
!                      to give all wave numbers in wave-set order (output)

!        Implicit arguments : NONE
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        MPP Group *ECMWF*

!     Modifications.
!     --------------
!        Original : 96-01-10
!        L.Isaksen: 96-02-02 - Calculation of KSPEC2MX added
!        K.YESSAD : 97-02-18 - Add KTMAX, bug correction for KSPOLEGL.
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE


!     DUMMY 
INTEGER(KIND=JPIM),INTENT(IN)  :: KSMAX
INTEGER(KIND=JPIM),INTENT(IN)  :: KTMAX
INTEGER(KIND=JPIM),INTENT(IN)  :: KPRTRW
INTEGER(KIND=JPIM),INTENT(IN)  :: KMYSETW
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC2
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPEC2MX
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KSPOLEGL

INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KASM0(0:KSMAX)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPROCM(0:KSMAX)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KUMPP(KPRTRW)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KMYMS(KSMAX+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPOSSP(KPRTRW+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPTRMS(KPRTRW)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KALLMS(KSMAX+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KDIM0G(0:KSMAX)

!     LOCAL 
INTEGER(KIND=JPIM) :: IK, IL, IND, IPOS, ISPEC2P, JA, JM,JMLOC,IM
INTEGER(KIND=JPIM) :: ISPOLEGL,ISPEC2MX,IASM0(0:KSMAX),IPROCM(0:KSMAX)
INTEGER(KIND=JPIM) :: IUMPP(KPRTRW),IMYMS(KSMAX+1),IPOSSP(KPRTRW+1)
INTEGER(KIND=JPIM) :: IPTRMS(KPRTRW),IALLMS(KSMAX+1),IDIM0G(0:KSMAX)
INTEGER(KIND=JPIM) :: ISPEC(KPRTRW),IC(KPRTRW)


!      -----------------------------------------------------------------

!*       1.    Initialize partitioning of wave numbers to PEs
!              ----------------------------------------------

ISPEC(:) = 0

IUMPP(:) = 0
IASM0(:) = -99
ISPOLEGL = 0

IL  = 1
IND = 1
IK  = 0
IPOS = 1
DO JM=0,KSMAX
  IK = IK + IND
  IF (IK > KPRTRW) THEN
    IK = KPRTRW
    IND = -1
  ELSEIF (IK < 1) THEN
    IK = 1
    IND = 1
  ENDIF
  IPROCM(JM) = IK
  ISPEC(IK) = ISPEC(IK)+KSMAX-JM+1
  IUMPP(IK) = IUMPP(IK)+1
  IF (IK == KMYSETW) THEN
    ISPOLEGL = ISPOLEGL +KTMAX+1-JM+1
    IMYMS(IL) = JM
    IASM0(JM) = IPOS
    IPOS = IPOS+(KSMAX-JM+1)*2
    IL = IL+1
  ENDIF
ENDDO

IPOSSP(1) = 1
ISPEC2P = 2*ISPEC(1)
ISPEC2MX = ISPEC2P
IPTRMS(1) = 1
DO JA=2,KPRTRW
  IPOSSP(JA) = IPOSSP(JA-1)+ISPEC2P
  ISPEC2P = 2*ISPEC(JA)
  ISPEC2MX = MAX(ISPEC2MX,ISPEC2P)
! pointer to the first wave number of a given wave-set in NALLMS array
  IPTRMS(JA) = IPTRMS(JA-1)+IUMPP(JA-1)
ENDDO
IPOSSP(KPRTRW+1) = IPOSSP(KPRTRW)+ISPEC2P

!  IALLMS :  wave numbers for all wave-set concatenated together to give all
!            wave numbers in wave-set order.
IC(:) = 0
DO JM=0,KSMAX
  IALLMS(IC(IPROCM(JM))+IPTRMS(IPROCM(JM))) = JM
  IC(IPROCM(JM)) = IC(IPROCM(JM))+1
ENDDO

IPOS = 1
DO JA=1,KPRTRW
  DO JMLOC=1,IUMPP(JA)
    IM = IALLMS(IPTRMS(JA)+JMLOC-1)
    IDIM0G(IM) = IPOS
    IPOS = IPOS+(KSMAX+1-IM)*2
  ENDDO
ENDDO

IF(PRESENT(KSPEC))    KSPEC  = ISPEC(KMYSETW)
IF(PRESENT(KSPEC2))   KSPEC2 = 2*ISPEC(KMYSETW)
IF(PRESENT(KSPEC2MX)) KSPEC2MX = ISPEC2MX
IF(PRESENT(KSPOLEGL)) KSPOLEGL = ISPOLEGL

IF(PRESENT(KASM0))  KASM0(:)  = IASM0(:)
IF(PRESENT(KPROCM)) KPROCM(:) = IPROCM(:)
IF(PRESENT(KUMPP))  KUMPP(:)  = IUMPP(:)
IF(PRESENT(KMYMS))  KMYMS(:)  = IMYMS(:)
IF(PRESENT(KPOSSP)) KPOSSP(:) = IPOSSP(:)
IF(PRESENT(KPTRMS)) KPTRMS(:) = IPTRMS(:)
IF(PRESENT(KALLMS)) KALLMS(:) = IALLMS(:)
IF(PRESENT(KDIM0G)) KDIM0G(:) = IDIM0G(:)

END SUBROUTINE SUWAVEDI
END MODULE SUWAVEDI_MOD


