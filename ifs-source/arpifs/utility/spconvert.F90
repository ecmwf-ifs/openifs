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

SUBROUTINE SPCONVERT(YDGEOMETRY,LDMODEL_TO_FILE, YDSPEC)

!**** *SPCONVERT*  - Spectral array conversion

!     Purpose.
!     --------
!        To convert model spectral fields into the fields actually written on file, 
!        NH variables are scaled.

!**   Interface.
!     ----------
!        *CALL* *SPCONVERT(LDMODEL_TO_FILE)

!        Explicit arguments :
!        --------------------
!           LDMODEL_TO_FILE : .TRUE. to convert from model to file
!                             .FALSE. to convert from file to model

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.    
!     ----------    

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        R. El Khatib  *METEO-FRANCE*

!     Modifications.
!     --------------
!        Original : 98-03-30
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE SPECTRAL_FIELDS_MOD, ONLY : SPECTRAL_FIELD, ASSIGNMENT(=)

IMPLICIT NONE

TYPE(GEOMETRY)      , INTENT(IN)    :: YDGEOMETRY
LOGICAL             , INTENT(IN)    :: LDMODEL_TO_FILE 
TYPE(SPECTRAL_FIELD), INTENT(INOUT) :: YDSPEC

INTEGER(KIND=JPIM) :: JJ, JLEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!      -----------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SPCONVERT',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, &
  & YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDLAP=>YDGEOMETRY%YRLAP)
ASSOCIATE(NSPEC2=>YDDIM%NSPEC2, &
 & NFLEVL=>YDDIMV%NFLEVL, &
 & NVALUE=>YDLAP%NVALUE, RLAPDI=>YDLAP%RLAPDI, RLAPIN=>YDLAP%RLAPIN)
IF (LDMODEL_TO_FILE) THEN

!*       1.    FROM MODEL FIELDS TO FILE
!              -------------------------

  DO JLEV=1,NFLEVL
!OCL NOVREC
    DO JJ=1,NSPEC2
!     VOR,DIV overwritten by PSI,KHY
      YDSPEC%VOR(JLEV,JJ)=YDSPEC%VOR(JLEV,JJ)*RLAPIN(NVALUE(JJ))
      YDSPEC%DIV(JLEV,JJ)=YDSPEC%DIV(JLEV,JJ)*RLAPIN(NVALUE(JJ))
    ENDDO
  ENDDO

ELSE

!*       2.    FROM FILE FIELDS TO MODEL
!              -------------------------

  DO JLEV=1,NFLEVL
!OCL NOVREC
    DO JJ=1,NSPEC2
      YDSPEC%VOR(JLEV,JJ)=YDSPEC%VOR(JLEV,JJ)*RLAPDI(NVALUE(JJ))
      YDSPEC%DIV(JLEV,JJ)=YDSPEC%DIV(JLEV,JJ)*RLAPDI(NVALUE(JJ))
    ENDDO
  ENDDO

ENDIF

!     ------------------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPCONVERT',1,ZHOOK_HANDLE)
END SUBROUTINE SPCONVERT
