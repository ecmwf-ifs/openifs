SUBROUTINE SU0PHY1S(KULOUT)

USE PARKIND1  ,ONLY : JPIM     ,JPRB,   JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE YOMLUN1S , ONLY : NULNAM
USE YOEPHY   , ONLY : LERADS   ,LESICE   ,LESURF   ,LEVDIF ,RTHRFRTI,&
     &                LEVGEN   ,LESSRO   ,LESN09,&
     &                NALBEDOSCHEME,NEMISSSCHEME,LEOCWA   ,LEOCCO  ,LEOCSA ,LEOCLA, &
     &                LEFLAKE  ,LEOCML   ,LELAIV  ,LECTESSEL, LEAGS, RLAIINT, &
     &                LWCOU    ,LWCOU2W  ,LWCOUHMF, &
     &                LEFARQUHAR, LEC4MAP, LEAIRCO2COUP, LEOPTSURF, &
     &                LECLIM10D,LESNML   ,LEURBAN ,LEINTWIND, NSNMLWS, LECMF1WAY,LECMF2LAKEC
USE YOMLOG1S , ONLY : LWRLKE

#ifdef DOC
! (C) Copyright 1991- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SU0PHY*   - Initialize common YOxPHY controlling physics

!     Purpose.
!     --------
!           Initialize YOxPHY, the common that includes the
!           basic switches for the physics of the model.

!**   Interface.
!     ----------
!        *CALL* *SU0PHY(KULOUT) from SU0YOM1S

!        Explicit arguments :
!        --------------------
!        KULOUT : Logical unit for the output

!        Implicit arguments :
!        --------------------
!        COMMON YOMPHY, YOEPHY

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!        or

!        Documentation ARPEGE (depending on which physics will be used)

!     Author.
!     -------
!        J.-J. Morcrette                    *ECMWF*
!        J.-F. Geleyn for the ARPEGE rewriting.

!     Modifications.
!     --------------
!        Original : 91-11-12
!        Modified 92-02-22 by M. Deque (tests of consistency with *YOMDPHY*)
!        Modified by R. EL Khatib : 93-04-02 Set-up defaults controled by LECMWF
!        Modified 94-02-28 by M.  Deque  : Shallow convection clouds
!        Modified 93-10-28 by Ph. Dandin : FMR scheme with MF physics
!        Modified 93-08-24 by D. Giard (test of consistency with digital filter)
!        Modified by M. Hamrud    : 93-06-05 Make use of LECMWF for ECMWF
!        Modified 95-06-24 by Jean-Francois Mahfouf for 1D surface scheme only
!        Modified 95-03-30 by D. Giard (test of consistency NDPSFI,LTWOTL)
!        Modified 95-11-27 by M. Deque (2nd call to APLPAR)
!        Modified by F. Rabier    : 96-09-25 Full physics set-up for 801 job

!        P. Viterbo   ECMWF   03-12-2004  Include user-defined RTHRFRTI
!        Y. Takaya    ECMWF   07-10-2008  Implement ocean mixed layer model
!        E. Dutra             16-11-2009  snow 2009 cleaning
!        S. Boussetta/G.Balsamo May 2010  Include CTESSEL switch LECTESSEL
!        G.Balsamo/S. Boussetta June 2011 Include switch LEAGS (for modularity CO2&Evap)
!        R. Hogan             14-01-2019  Changed LE4ALB to NALBEDOSCHEME
!        A. Agusti-Panareda 18-11-2020 Include LEAIRCO2COUP to use variable air CO2 in photosynthesis
!        A. Agusti-Panareda 06-07-2021 Include LEFARQUHAR switch for Farquhar photosynthesis model
!        A. Agusti-Panareda (Jul 2021):Add LEC4MAP flag for C4 photosynthesis
!     ------------------------------------------------------------------
#endif

IMPLICIT NONE
INTEGER(KIND=JPIM) :: KULOUT, JSFC, IPATH

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "namphy1s.h"

IF (LHOOK) CALL DR_HOOK('SU0PHY1S',0,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

!*       1.    Set default values.
!              -------------------

LERADS=.TRUE.
LESICE=.TRUE.
LESURF=.TRUE.
LEVDIF=.TRUE.
LEVGEN=.TRUE.
LESSRO=.TRUE.
LESN09=.TRUE.
NALBEDOSCHEME=0
NEMISSSCHEME=1
LEOCWA=.FALSE. 
LEOCCO=.FALSE. 
LEOCSA=.FALSE.
LEOCLA=.FALSE.
LEFLAKE=.FALSE. 
LWCOU=.FALSE.
LWCOU2W=.FALSE.
LWCOUHMF=.FALSE.
LEOCML=.FALSE. 
LELAIV=.FALSE.
LECTESSEL=.FALSE.
LEAGS=.FALSE.
LEFARQUHAR=.FALSE.
LEOPTSURF=.FALSE.
LEC4MAP=.FALSE.
LEAIRCO2COUP=.FALSE.
LECLIM10D=.FALSE.
LESNML=.FALSE.
LEURBAN=.FALSE.
RLAIINT=0.0_JPRB
RTHRFRTI=0.0_JPRB
LEINTWIND=.FALSE.
NSNMLWS=1_JPIM
LECMF1WAY=.FALSE.
!     ------------------------------------------------------------------


!*       2.    Print final values.
!              -------------------

REWIND(NULNAM)
READ(NULNAM,NAMPHY)

WRITE(UNIT=KULOUT,FMT='('' COMMON YOEPHY '')')
WRITE(UNIT=KULOUT,FMT='('' LERADS = '',L5 &
     &     ,'' LESURF = '',L5,'' LEVDIF = '',L5 &
     &     ,'' LESICEF = '',L5 &
     &     )') &
     & LERADS,LESURF,LEVDIF,LESICE
WRITE(UNIT=KULOUT,FMT='('' LEVGEN = '',L5 &
      &     ,'' LESSRO  = '',L5 )') LEVGEN,LESSRO 
WRITE(UNIT=KULOUT,FMT='('' LEFLAKE = '',L5)') LEFLAKE  
WRITE(UNIT=KULOUT,FMT='('' LEOCML = '',L5)') LEOCML 
WRITE(UNIT=KULOUT,FMT='('' LESN09 = '',L5)') LESN09
WRITE(UNIT=KULOUT,FMT='('' NALBEDOSCHEME = '',I0)') NALBEDOSCHEME
WRITE(UNIT=KULOUT,FMT='('' NEMISSSCHEME = '',I0)') NEMISSSCHEME
WRITE(UNIT=KULOUT,FMT='('' LESNML = '',L5)') LESNML
WRITE(UNIT=KULOUT,FMT='('' LEURBAN = '',L5)') LEURBAN
WRITE(UNIT=KULOUT,FMT='('' NSNMLWS= '',I3)') NSNMLWS
WRITE(UNIT=KULOUT,FMT='('' LELAIV = '',L5)') LELAIV
WRITE(UNIT=KULOUT,FMT='('' LECTESSEL = '',L5)') LECTESSEL
WRITE(UNIT=KULOUT,FMT='('' LEAGS = '',L5)') LEAGS
WRITE(UNIT=KULOUT,FMT='('' LEFARQUHAR = '',L5)') LEFARQUHAR
WRITE(UNIT=KULOUT,FMT='('' LEOPTSURF = '',L5)') LEOPTSURF
WRITE(UNIT=KULOUT,FMT='('' LEC4MAP = '',L5)') LEC4MAP
WRITE(UNIT=KULOUT,FMT='('' LEAIRCO2COUP = '',L5)') LEAIRCO2COUP
WRITE(UNIT=KULOUT,FMT='('' LECLIM10D = '',L5)') LECLIM10D
WRITE(UNIT=KULOUT,FMT='('' RLAIINT = '',f4.2)') RLAIINT
WRITE(UNIT=KULOUT,FMT='('' LECMF1WAY = '',L5)') LECMF1WAY

!     ------------------------------------------------------------------
IF ( .NOT. LEFLAKE ) THEN
  IF ( LWRLKE ) THEN
    WRITE(UNIT=KULOUT,FMT='('' LEFLAKE AND LWRLKE INCOMPAT! - CHANGING LWRLKE TO FALSE '')')
    LWRLKE=.FALSE.
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SU0PHY1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SU0PHY1S
