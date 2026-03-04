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

SUBROUTINE RESET_ACCFIE_VAREPS(YDGEOMETRY,YDSURF,YDMCC,YDEPHY,YDRIP,YDDYN,KSTEP)

USE YOEPHY             , ONLY : TEPHY
USE YOMDYN             , ONLY : TDYN
USE YOMMCC             , ONLY : TMCC
USE YOMRIP             , ONLY : TRIP
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE SURFACE_FIELDS_MIX , ONLY : TSURF
USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT
USE YOMCST   , ONLY : RDAY     ,REA      ,REPSM    ,&
 & RV       ,RCPV     ,RETV     ,&
 & RCW      ,RCS      ,RLVTT    ,RLSTT    ,RTT      ,&
 & RALPW    ,RBETW    ,RGAMW    ,RALPS    ,RBETS    ,&
 & RGAMS    ,RALPD    ,RBETD    ,RGAMD   ,RG
USE YOMVAREPS, ONLY : NST_FCLENGTH_INI, NAFVAREPS, NAFVAREPSGC
USE YOM_GRIB_CODES  , ONLY : &
 & NGRBSWVL1 ,NGRBSD   ,NGRBLSP  ,&
 & NGRBCP   ,NGRBSF    ,NGRBFZRA ,NGRBSDOR ,NGRBSWVL2 ,NGRBSR   ,&
 & NGRBE    ,NGRBSWVL3 ,NGRBSRC  ,NGRBRO   ,NGRBSWVL4 ,&
 & NGRBCSF  ,NGRBLSF  ,NGRBFSR   ,NGRBTP   ,&
 & NGRBSRO  ,NGRBSSRO ,NGRBPEV
USE YOMGRIB  , ONLY : NLEG
USE IOSTREAM_MIX , ONLY : SETUP_IOSTREAM, SETUP_IOREQUEST, Y_IOSTREAM_FDB,&
 & IO_PUT, CLOSE_IOSTREAM, TYPE_IOSTREAM , TYPE_IOREQUEST
USE YOMMP0   , ONLY : NOUTTYPE
USE YOMOPH0  , ONLY : CFNVAREPS

!**** *RESET_ACCFIE_VAREPS*  - UPDATE ACCUMULATED VAREPS FIELDS

!     Purpose.
!     --------
!     Update accumulated fields for VAREPS lower resolution leg

!**   Interface.
!     ----------
!        *CALL* *RESET_ACCFIE_VAREPS(KSTEP)

!        Explicit arguments :
!        --------------------
!        KSTEP   : TIME STEP INDEX

!        Implicit arguments :
!        --------------------
!        

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
!        Roberto Buizza *ECMWF*

!     Modifications.
!     --------------
!        Original : 24 June 2005
!        M.Hamrud      01-Jul-2006 Revised surface fields
!        M.Leutbecher  23-Mar-2011 Write interpolated accumulated fields in 
!                                  overlap stream, Table=230
!        R. Forbes    10-Jan-2015 Added FZRA 
!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(IN)    :: YDGEOMETRY
TYPE(TSURF)   , INTENT(INOUT) :: YDSURF
TYPE(TDYN)     ,INTENT(INOUT) :: YDDYN
TYPE(TEPHY)    ,INTENT(INOUT) :: YDEPHY
TYPE(TMCC)     ,INTENT(INOUT) :: YDMCC
TYPE(TRIP)     ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN) :: KSTEP 
INTEGER(KIND=JPIM), EXTERNAL  :: ISRCHEQ
INTEGER(KIND=JPIM) :: IFILE, JL, IFIELDS, JSTGLO, IEND, IBLK, IFLDP, JV, IGRBCODE
INTEGER(KIND=JPIM), ALLOCATABLE, DIMENSION(:) :: IGRIB2DX
LOGICAL         :: LLDSPOR
REAL(KIND=JPRB) :: ZR3(YDGEOMETRY%YRGEM%NGPTOT,NAFVAREPS,1)
REAL(KIND=JPRB) :: ZD1, ZD2, ZD3, ZD4, ZSCALE, ZSCALER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

TYPE(TYPE_IOSTREAM)  :: YL_IOSTREAM
TYPE(TYPE_IOREQUEST) :: YL_IOREQUEST

!     ------------------------------------------------------------------

!#include "gpnorm2.intfb.h"
#include "sugridg.intfb.h"

#include "fctast.func.h"
#include "fcttim.func.h"
#include "fcttrm.func.h"

IF (LHOOK) CALL DR_HOOK('RESET_ACCFIE_VAREPS',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT, &
 & SD_VD=>YDSURF%SD_VD, YSD_VD=>YDSURF%YSD_VD, YSD_VDD=>YDSURF%YSD_VDD)
!        Update accumulated fields for VAREPS

IF( (NLEG >= 2) .AND. (KSTEP == NST_FCLENGTH_INI) ) THEN
  WRITE(NULOUT,'(A,1X,I10)') ' VAREPS accumulated fields - UPDATE KSTEP=',KSTEP
  WRITE(NULOUT,'(A,1X,I10)') ' >> VAREPS accumulated fields - KSTEP=',KSTEP
  WRITE(NULOUT,'(A,1X,I10)') ' >> ......................... - NST_FCLENGTH_INI=',NST_FCLENGTH_INI
  WRITE(NULOUT,'(A,1X,I10)') ' >> ......................... - NLEG=',NLEG
  WRITE(NULOUT,'(A,1X,I10)') ' >> ......................... - NAFVAREPS=',NAFVAREPS
  WRITE(NULOUT,'(A,1X,I10)') ' >> ......................... - NGPTOT=',NGPTOT
  !debug WRITE(NULOUT,'(A)') ' VAREPS accumulated fields - norms BEFORE update'
  !debug CALL FLUSH(NULOUT)
  !debug  CALL GPNORM2(NVDIAG,MVDIAG,GPPBUF)

  IFILE=9
  IFIELDS=NAFVAREPS

  LLDSPOR=.FALSE.
  CALL SUGRIDG(YDGEOMETRY,YDSURF,YDEPHY,YDMCC,YDDYN,IFILE,LLDSPOR,PFIELD=ZR3(:,:,1))
  !   write interpolated accumulated fields from previous leg
  ALLOCATE (IGRIB2DX(IFIELDS))
  DO JV=1,NAFVAREPS
    IGRIB2DX(JV)=230*1000 + MOD(NAFVAREPSGC(JV),1000)
  ENDDO
  CALL SETUP_IOREQUEST(YL_IOREQUEST,'GRIDPOINT_FIELDS',LDGRIB=.TRUE.,&
   & KGRIB2D=IGRIB2DX,KPROMA=NGPTOT,PTSTEP=YDRIP%TSTEP)
  DEALLOCATE (IGRIB2DX)
  IF(NOUTTYPE == 2) THEN
    CALL IO_PUT(Y_IOSTREAM_FDB,YL_IOREQUEST,PR3=ZR3)
    WRITE(NULOUT,'(I3,A)') IFIELDS,&
         & ' lower resolution accumulated fields written to FDB ',KSTEP
  ELSE
    CALL SETUP_IOSTREAM(YL_IOSTREAM,'CIO',TRIM(CFNVAREPS)//'_lowres',CDMODE='w',KIOMASTER=1)
    CALL IO_PUT(YL_IOSTREAM,YL_IOREQUEST,PR3=ZR3)
    CALL CLOSE_IOSTREAM(YL_IOSTREAM)
    WRITE(NULOUT,'(I3,A)') IFIELDS,&
         & ' lower resolution accumulated fields written to  '//TRIM(CFNVAREPS)//'_lowres'
  ENDIF

  ZD1 = 0.07_JPRB
  ZD2 = 0.21_JPRB
  ZD3 = 0.72_JPRB
  ZD4 = 1.89_JPRB

  
  DO JV=1,YSD_VDD%NUMFLDS
    IGRBCODE = YSD_VD%YVD(JV)%IGRBCODE
    IF(IGRBCODE /= -999) THEN
      IFLDP=ISRCHEQ(IFIELDS,NAFVAREPSGC,1,IGRBCODE)
      !debug WRITE(NULOUT,'(">> RESET_ACCFIE_VAREPS: JV,IFLDP=",2(I10,1X))') JV,IFLDP
      IF(IFLDP <= IFIELDS) THEN

        ! Define re-scaling factor for fields read from MARS (as in wrmlpplg.F90)

        IF (IGRBCODE == NGRBSWVL1 ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSD  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBLSP ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBCP  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBTP  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSRO ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSSRO) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSF  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBFZRA  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSWVL2) THEN
          ZSCALE = 1.0E-3_JPRB*ZD1/ZD2
        ELSEIF(IGRBCODE == NGRBSR  ) THEN
          ZSCALE = 1.0_JPRB/RG
        ELSEIF(IGRBCODE == NGRBE   ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBPEV   ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSWVL3) THEN
          ZSCALE = 1.0E-3_JPRB*ZD1/ZD3
        ELSEIF(IGRBCODE == NGRBSRC ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBRO  ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBSWVL4) THEN
          ZSCALE = 1.0E-3_JPRB*ZD1/ZD4
        ELSEIF(IGRBCODE == NGRBCSF ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBLSF ) THEN
          ZSCALE = 1.0E-3_JPRB
        ELSEIF(IGRBCODE == NGRBFSR ) THEN
          ZSCALE = 1.0_JPRB/RG
        ELSEIF(IGRBCODE == NGRBSDOR) THEN
          ZSCALE = 1.0_JPRB/RG
        ELSE
          ZSCALE = 1.0_JPRB
        ENDIF
        ZSCALER = 1.0_JPRB / ZSCALE
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JSTGLO,IEND,IBLK,JL)
        DO JSTGLO=1,NGPTOT,NPROMA

          IEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
          IBLK=(JSTGLO-1)/NPROMA+1
          DO JL=1,IEND
            SD_VD(JL,JV,IBLK) = ZSCALER*ZR3(JSTGLO+JL-1,IFLDP,1)
          ENDDO
        ENDDO
!$OMP END PARALLEL DO
      ENDIF
    ENDIF

  ENDDO
  !debug WRITE(NULOUT,'(A)') ' VAREPS accumulated fields - norms AFTER update'
  !debug CALL FLUSH(NULOUT)
  !deubg CALL GPNORM2(NVDIAG,MVDIAG,GPPBUF)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RESET_ACCFIE_VAREPS',1,ZHOOK_HANDLE)

END SUBROUTINE RESET_ACCFIE_VAREPS
