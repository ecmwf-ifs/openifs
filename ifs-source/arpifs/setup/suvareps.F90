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

SUBROUTINE SUVAREPS(YDRIP)

!**** *SUVAREPS*   - Routine to initialize variables that control 
!                    VARiable Resolution EPS (VAREPS)

!     Purpose.
!     --------
!           Initialize VAREPS control common
!**   Interface.
!     ----------
!        *CALL* *SUVAREPS

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        COMMON YOMVAREPS

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Authors.
!     -------
!        Roberto Buizza  *ECMWF*
!        Original : 2004-01-29

!     Modifications.
!     --------------
!        MODIFIED : 2008-03-11 F. Vitart (add STRD, SSR and STR to the list of possible accumulated parameters)
!        K. Yessad (July 2014): Move some variables.
!        2015-01-10 R. Forbes Added freezing rain FZRA
!     -----------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT   ,NULNAM
USE YOMVAREPS, ONLY : NFCHO_TRUNC_INI, NFCLENGTH_INI, NST_FCLENGTH_INI,&
                    & NFCHO_TRUNC, NST_TRUNC, NST_TRUNC_INI,&
                    & NAFVAREPSMAX,NAFVAREPS,NAFVAREPSGC,&
                    & LVAREPS
USE YOM_GRIB_CODES , ONLY : NGRBES   ,NGRBSMLT ,NGRBLSPF ,NGRBLSP  ,NGRBCP   ,&
                    & NGRBSF   ,NGRBFZRA ,NGRBBLD  ,NGRBSSHF ,NGRBSLHF ,&
                    & NGRBSSRD ,NGRBSTRD  ,NGRBSSR  ,NGRBSTR  , NGRBTSR ,&
                    & NGRBTTR  ,NGRBEWSS ,NGRBNSSS ,&
                    & NGRBE    ,NGRBSUND ,NGRBLGWS ,NGRBMGWS ,NGRBGWD  ,&
                    & NGRBPEV  ,NGRBRO   ,NGRBTSRC ,NGRBTTRC ,NGRBSSRC ,&
                    & NGRBSTRC ,NGRBTP   ,NGRBCSF  ,NGRBLSF  ,&
                    & NGRBSRO  ,NGRBSSRO ,NGRBFDIR ,NGRBCDIR , NGRBDSRP,&
                    & NGRBPARCS,NGRBUVB  ,NGRBPAR  ,NGRBTISR ,&
                    & NGRBVIMD ,NGRBSSRDC,NGRBSTRDC,&
                    & NGRBNEE  ,NGRBGPP  ,NGRBREC
USE YOMRIP   , ONLY : TRIP
USE YOMFPC   , ONLY : TNAMFPL

!     -----------------------------------------------------------------

IMPLICIT NONE

TYPE(TRIP),INTENT(INOUT):: YDRIP
INTEGER(KIND=JPIM) :: ISECS
INTEGER(KIND=JPIM) :: JAF,JSFC,IAFAF,IAFGC(NAFVAREPSMAX)
TYPE(TNAMFPL) :: YLNAMFPC
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     -----------------------------------------------------------------

#include "sufpc.intfb.h"
#include "abor1.intfb.h"
#include "posnam.intfb.h"
#include "namvareps.nam.h"

!     -----------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUVAREPS',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP,MFPPHY=>YLNAMFPC%MFPPHY, NFPPHY=>YLNAMFPC%NFPPHY)
!     -----------------------------------------------------------------

!*    1. Set default values.
!     ----------------------

!     1.1 Set implicit default values

LVAREPS=.FALSE.
NFCHO_TRUNC_INI=0
NFCLENGTH_INI=0
NFCHO_TRUNC=0
NST_TRUNC_INI=0
NST_TRUNC=0
NAFVAREPS=0
DO JAF=1,NAFVAREPSMAX
  NAFVAREPSGC(JAF)=0
ENDDO

!     -----------------------------------------------------------------

!*    2. Read Namelist
!     ----------------

CALL POSNAM(NULNAM,'NAMVAREPS')
READ(NULNAM,NAMVAREPS)

CALL SUFPC(YDNAMFPL=YLNAMFPC,LDPRINT=.FALSE.)

ISECS=3600.0
IF( (NFCLENGTH_INI-NFCHO_TRUNC_INI) >= 0 ) THEN
  NST_FCLENGTH_INI=NINT(((NFCLENGTH_INI-NFCHO_TRUNC_INI)*ISECS)/TSTEP)
ELSE
  NST_FCLENGTH_INI=0
ENDIF

NST_TRUNC_INI=NINT((NFCHO_TRUNC_INI*ISECS)/TSTEP)
NST_TRUNC=NINT((NFCHO_TRUNC*ISECS)/TSTEP)

IF(LVAREPS) THEN

! List of accumulated fields (NB: IAFGC is dimensioned by NAFVAREPSMAX, set in YOMVAREPS)

  IAFGC(1)=NGRBES
  IAFGC(2)=NGRBSMLT
  IAFGC(3)=NGRBLSPF
  IAFGC(4)=NGRBLSP
  IAFGC(5)=NGRBCP
  IAFGC(6)=NGRBSF
  IAFGC(7)=NGRBBLD
  IAFGC(8)=NGRBSSHF
  IAFGC(9)=NGRBSLHF
  IAFGC(10)=NGRBSSRD
  IAFGC(11)=NGRBSTRD
  IAFGC(12)=NGRBSSR
  IAFGC(13)=NGRBSTR
  IAFGC(14)=NGRBTSR
  IAFGC(15)=NGRBTTR
  IAFGC(16)=NGRBEWSS
  IAFGC(17)=NGRBNSSS
  IAFGC(18)=NGRBE
  IAFGC(19)=NGRBSUND
  IAFGC(20)=NGRBLGWS
  IAFGC(21)=NGRBMGWS
  IAFGC(22)=NGRBGWD
  IAFGC(23)=NGRBPEV
  IAFGC(24)=NGRBRO
  IAFGC(25)=NGRBTSRC
  IAFGC(26)=NGRBTTRC
  IAFGC(27)=NGRBSSRC
  IAFGC(28)=NGRBSTRC
  IAFGC(29)=NGRBTP
  IAFGC(30)=NGRBCSF
  IAFGC(31)=NGRBLSF
  IAFGC(32)=NGRBSRO
  IAFGC(33)=NGRBSSRO
  IAFGC(34)=NGRBFDIR
  IAFGC(35)=NGRBCDIR
  IAFGC(36)=NGRBPARCS
  IAFGC(37)=NGRBVIMD
  IAFGC(38)=NGRBUVB
  IAFGC(39)=NGRBPAR
  IAFGC(40)=NGRBTISR
  IAFGC(41)=NGRBSSRDC
  IAFGC(42)=NGRBSTRDC
  IAFGC(43)=NGRBNEE
  IAFGC(44)=NGRBGPP
  IAFGC(45)=NGRBREC
  IAFGC(46)=NGRBFZRA
  IAFGC(47)=NGRBDSRP


! Check whether any of the post-processed surface fields is an accumulated one,
! and add it to the list NAFVAREPSGC

  IAFAF=0
  DO JSFC=1,NFPPHY
    DO JAF=1,NAFVAREPSMAX
      IF(MFPPHY(JSFC) == IAFGC(JAF)) THEN
        IAFAF=IAFAF+1
        IF (IAFAF > NAFVAREPSMAX) THEN
          WRITE(NULOUT,'('' >> Sub SUVAREPS ERROR - IAFAF gt NAFVAREPSMAX'')')
          WRITE(NULOUT,'(''    IAFAF, NAFVAREPSMAX ='',2(I10,1X))') IAFAF, NAFVAREPSMAX
          CALL ABOR1('SUVAREPS: ABOR1 CALLED')
        ENDIF
        NAFVAREPSGC(IAFAF)=MFPPHY(JSFC)
        WRITE(NULOUT,'('' >> Sub SUVAREPS found grid file'')')
        WRITE(NULOUT,'(''    IAFAF, JSFC, MFPPHY(JSFC)='',3(I10,1X))') IAFAF,JSFC,MFPPHY(JSFC)
      ENDIF
    ENDDO
  ENDDO
  NAFVAREPS=IAFAF
ENDIF

!     -----------------------------------------------------------------

!*    3. Print final values.
!     ----------------------

IF (LVAREPS) THEN 
 WRITE(NULOUT,'('' >> Sub SUVAREPS'')')
 WRITE(NULOUT,'('' VARiable Resolution EPS     = '',L1)') LVAREPS
 WRITE(NULOUT,'('' VAREPS variables in common YOMVAREPS'')')
 WRITE(NULOUT,'('' Fc-step used to define ICS NFCHO_TRUNC_INI= '',I10)') NFCHO_TRUNC_INI
 WRITE(NULOUT,'('' Fc length previous leg       NFCLENGTH_INI= '',I10)') NFCLENGTH_INI
 WRITE(NULOUT,'('' .. NST Fc len prev leg    NST_FCLENGTH_INI= '',I10)') NST_FCLENGTH_INI
 WRITE(NULOUT,'('' Truncation forecast time       NFCHO_TRUNC= '',I10)') NFCHO_TRUNC
 WRITE(NULOUT,'('' >> SUVAREPS: Accumulated fields to be re-set'')')
 WRITE(NULOUT,'(''    Tot number of Acc Fields  NAFVAREPS= '',I10)') NAFVAREPS
 DO JAF=1,NAFVAREPS
   WRITE(NULOUT,'(''    List of accum fields '',2(1X,I10))') JAF,NAFVAREPSGC(JAF)
 ENDDO
ENDIF

!     -----------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVAREPS',1,ZHOOK_HANDLE)
END SUBROUTINE SUVAREPS
