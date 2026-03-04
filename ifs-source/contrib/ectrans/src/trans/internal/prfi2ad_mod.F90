! (C) Copyright 1987- ECMWF.
! (C) Copyright 2013- Meteo-France.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

MODULE PRFI2AD_MOD
CONTAINS
SUBROUTINE PRFI2AD(KM,KMLOC,KF_FS,PAIA,PSIA)

!**** *PRFI2AD* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!        *CALL* *PRFI2AD(..)

!        Explicit arguments :
!        --------------------  KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAIA - antisymmetric part of Fourier
!                              components for KM  (output)
!                              PSIA - symmetric part of Fourier
!                              components for KM  (output)

!        Implicit arguments :  The Grid point arrays of the model.
!        --------------------

!     Method.
!     -------

!     Externals.   PRFI2ADB   - basic copying routine
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 87-11-25
!        Modified : 91-07-01 Philippe Courtier/Mats Hamrud - Rewrite
!                            for uv formulation
!        Modified : 93-03-19 D. Giard - CDCONF='T'
!        Modified : 93-??-?? ???????? - CDCONF='C'--> bug if CDCONF='T'
!        Modified : 93-05-13 D. Giard - correction of the previous bug
!        Modified : 93-11-18 M. Hamrud - use only one Fourier buffer
!        Modified : 94-04-06 R. El Khatib - Full-POS configuration 'P'
!        M.Hamrud : 94-11-01 New conf 'G' - vor,div->vor,div
!                            instead of u,v->vor,div
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE PRFI2BAD_MOD    ,ONLY : PRFI2BAD
!

IMPLICIT NONE

INTEGER(KIND=JPIM) , INTENT(IN) :: KM
INTEGER(KIND=JPIM) , INTENT(IN) :: KMLOC
INTEGER(KIND=JPIM) , INTENT(IN) :: KF_FS

REAL(KIND=JPRB) , INTENT(IN) :: PSIA(:,:),   PAIA(:,:)

!     ------------------------------------------------------------------

!*       2.    EXTRACT SYM./ANTISYM. FIELDS FROM TIME T+1.
!              -------------------------------------------

CALL PRFI2BAD(KF_FS,KM,KMLOC,PAIA,PSIA)

!     ------------------------------------------------------------------

END SUBROUTINE PRFI2AD
END MODULE PRFI2AD_MOD
