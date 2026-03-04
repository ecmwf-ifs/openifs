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

SUBROUTINE WRMLPPLG(YDGEOMETRY,YDSURF,YDEPHY,YDEPHLI,YDRIP,YDDYN,YDPHY)

!**** *WRMLPPLG*  - writes out some model level fields in GRIB

!     Purpose.
!     --------
!     Write out model level fields in lagged mode, after going
!     into grid-point space and calling the physics.

!**   Interface.
!     ----------
!        *CALL* *WRMLPPLG

!        Implicit arguments :      The state variables of the model
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        MARS users guide

!     Author.
!     -------
!      Mats Hamrud *ECMWF*
!      Original : 93-10-25 (From WRMLPPL)

!     Modifications.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!     ------------------------------------------------------------------

USE YOMDYN             , ONLY : TDYN
USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF, GPOPER, TYPE_SFL_COMM
USE PARKIND1           , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT
USE YOMCT0             , ONLY : NCONF, LOBSC1   
USE YOMPHY             , ONLY : TPHY
USE YOEPHY             , ONLY : TEPHY
USE YOEPHLI            , ONLY : TEPHLI
USE YOMVAR             , ONLY : LTWBGV, LTWCGL  
USE YOM_GRIB_CODES     , ONLY :&
 & NGRBBV   ,NGRBSWVL1 ,NGRBSD   ,NGRBSDSL ,NGRBLSP  ,&
 & NGRBCP   ,NGRBSF    ,NGRBFZRA ,NGRBSDOR ,NGRBSWVL2,NGRBSR   ,&
 & NGRBE    ,NGRBSWVL3 ,NGRBSRC  ,NGRBRO   ,NGRBSWVL4,&
 & NGRBCSF  ,NGRBLSF   ,NGRBFSR  
USE YOMCST             , ONLY : RG
USE YOMPPC             , ONLY : M2DGGPL, NO2DGG, NO2DGGL, LRSUP  
USE YOMCT3   , ONLY : NSTEP
USE YOMOPH0  , ONLY : LINC

IMPLICIT NONE

TYPE(GEOMETRY),INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)   ,INTENT(INOUT) :: YDSURF
TYPE(TDYN)    ,INTENT(INOUT) :: YDDYN
TYPE(TEPHY)   ,INTENT(INOUT) :: YDEPHY
TYPE(TEPHLI)  ,INTENT(IN)    :: YDEPHLI
TYPE(TPHY)    ,INTENT(INOUT) :: YDPHY
TYPE(TRIP)    ,INTENT(INOUT) :: YDRIP

CHARACTER (LEN = 1) ::  CLMODE
CHARACTER (LEN = 30) :: CLFNGG
CHARACTER (LEN = 210) :: CLFNSH
REAL(KIND=JPRB), ALLOCATABLE :: ZREAL(:,:)
REAL(KIND=JPRB) :: ZFIELD(YDGEOMETRY%YRDIM%NPROMA,1)
INTEGER(KIND=JPIM) :: IGRIBIOGG(2,NO2DGG+1)

INTEGER(KIND=JPIM) ::  IEND,  IGGMAX, IGRIBCD,&
 & IST, J, JCODE, JF, JSTGLO ,IBLK

TYPE(TYPE_SFL_COMM) :: YLCOM
REAL(KIND=JPRB) :: ZD1, ZD2, ZD3, ZD4, ZSCALE
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fcttim.func.h"

#include "suecfname.intfb.h"
#include "wroutgpgb.intfb.h"

!     ------------------------------------------------------------------

!*       1.    PREPARATIONS.
!              --------------

IF (LHOOK) CALL DR_HOOK('WRMLPPLG',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM, YDGEM=>YDGEOMETRY%YRGEM)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, NRESOL=>YDDIM%NRESOL, &
 & LEPHYS=>YDEPHY%LEPHYS, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & LMPHYS=>YDPHY%LMPHYS, &
 & LTLEVOL=>YDEPHLI%LTLEVOL, &
 & TSTEP=>YDRIP%TSTEP)

CLMODE='a'

!      -----------------------------------------------------------

!*       4.0   WRITE OUT SURFACE GRID POINT FIELDS
!              -----------------------------------

IF ( .NOT.(LTWBGV.OR.LTWCGL) .AND.&
   & (LRSUP.AND.(LMPHYS.OR.LEPHYS).OR.LOBSC1.OR.&
   & NCONF == 131) ) THEN  

  ALLOCATE(ZREAL(NGPTOT,NO2DGGL))

!*       4.1   COPY FROM NPROMA FORMAT TO FIELD FORMAT


  IGGMAX=0
  DO JCODE=1,NO2DGGL
    IGRIBCD = M2DGGPL(JCODE)
    YLCOM%IGRBCODE = IGRIBCD
    CALL GPOPER(YDGEOMETRY%YRDIM,YDDYN,'GETGRIBPOS',YDSURF,KBL=1,YDCOM=YLCOM)
    IF(YLCOM%L_OK) THEN
      IGGMAX = IGGMAX+1
      IGRIBIOGG(1,IGGMAX) = IGRIBCD
      DO JSTGLO=1,NGPTOT,NPROMA
        IST=1
        IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
        IBLK=(JSTGLO-1)/NPROMA+1
        CALL GPOPER(YDGEOMETRY%YRDIM,YDDYN,'GETFIELD',YDSURF,KBL=IBLK,YDCOM=YLCOM,PFIELD=ZFIELD)
        DO J=JSTGLO,JSTGLO+IEND-IST
          ZREAL(J,IGGMAX)=ZFIELD(J-JSTGLO+1,1)
        ENDDO
      ENDDO
    ELSE
      IF(IGRIBCD /= NGRBBV) THEN
        WRITE(NULOUT,'(A,I4,A)') ' WRMLPPLG - GRIB CODE ',IGRIBCD,' NOT FOUND '
      ENDIF
    ENDIF
  ENDDO

  IGRIBIOGG(2,:) = 0


  DO JF=1,IGGMAX
    IGRIBCD = IGRIBIOGG(1,JF)
!*       4.2   SCALINGS FOR ECMWF MARS ARCHIVING

    ZD1 = 0.07_JPRB
    ZD2 = 0.21_JPRB
    ZD3 = 0.72_JPRB
    ZD4 = 1.89_JPRB
    IF    (IGRIBCD == NGRBSWVL1 ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBSD  ) THEN
      ZSCALE = 1._JPRB ! in kg m-2
    ELSEIF(IGRIBCD == NGRBSDSL  ) THEN
      ZSCALE = 1.0E-3_JPRB ! in m water eq.
    ELSEIF(IGRIBCD == NGRBLSP ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBCP  ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBSF  ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBFZRA) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBSWVL2) THEN
      ZSCALE = 1.0E-3_JPRB*ZD1/ZD2
    ELSEIF(IGRIBCD == NGRBSR  ) THEN
      ZSCALE = 1.0_JPRB/RG
    ELSEIF(IGRIBCD == NGRBE   ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBSWVL3) THEN
      ZSCALE = 1.0E-3_JPRB*ZD1/ZD3
    ELSEIF(IGRIBCD == NGRBSRC ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBRO  ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBSWVL4) THEN
      ZSCALE = 1.0E-3_JPRB*ZD1/ZD4
    ELSEIF(IGRIBCD == NGRBCSF ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBLSF ) THEN
      ZSCALE = 1.0E-3_JPRB
    ELSEIF(IGRIBCD == NGRBFSR ) THEN
      ZSCALE = 1.0_JPRB/RG
    ELSEIF(IGRIBCD == NGRBSDOR) THEN
      ZSCALE = 1.0_JPRB/RG
    ELSE
      ZSCALE = 1.0_JPRB
    ENDIF
    IF(ZSCALE /= 1.0_JPRB) THEN
      DO J=1,NGPTOT
        ZREAL(J,JF) = ZREAL(J,JF)*ZSCALE
      ENDDO
    ENDIF
    
  ENDDO

  CALL SUECFNAME(LINC,'s',NSTEP,YDRIP%TSTEP,YDRIP%NSTOP,CLFNGG,CLFNSH)
  CALL WROUTGPGB(ZREAL,NGPTOT,IGGMAX,IGRIBIOGG,'s',CLMODE,NRESOL,CLFNGG,YDRIP%TSTEP)
  DEALLOCATE(ZREAL)

ENDIF

!      -----------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('WRMLPPLG',1,ZHOOK_HANDLE)
END SUBROUTINE WRMLPPLG
