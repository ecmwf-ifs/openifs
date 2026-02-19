! (C) Copyright 2000- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE INI_SPEC_DIST(KSMAX,KTMAX,KPRTRW,KMYSETW,KASM0,KSPOLEGL,KPROCM,&
                    &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS,KPTRMS,KALLMS)


!**** *INI_SPEC_DIST* - Initialize spectral wave distribution

!     Purpose.
!     --------
!     Initialize arrays controlling spectral wave distribution

!**   Interface.
!     ----------
!     CALL INI_SPEC_DIST(...)

!     Explicit arguments :
!     --------------------
!           KSMAX    - spectral truncation required
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

!     Externals.  SUWAVEDI
!     ----------
!        Original : 00-03-03

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
!ifndef INTERFACE
USE SUWAVEDI_MOD    ,ONLY : SUWAVEDI
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!endif INTERFACE

IMPLICIT NONE

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
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPOSSP(KPRTRW+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KMYMS(KSMAX+1)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KPTRMS(KPRTRW)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KALLMS(KSMAX+1)

!ifndef INTERFACE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('INI_SPEC_DIST',0,ZHOOK_HANDLE)

CALL  SUWAVEDI(KSMAX,KTMAX,KPRTRW,KMYSETW,KASM0,KSPOLEGL,KPROCM,&
               &KUMPP,KSPEC,KSPEC2,KSPEC2MX,KPOSSP,KMYMS,KPTRMS,KALLMS)

IF (LHOOK) CALL DR_HOOK('INI_SPEC_DIST',1,ZHOOK_HANDLE)

!endif INTERFACE

END SUBROUTINE INI_SPEC_DIST
