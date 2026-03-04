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

SUBROUTINE SUCT1

!**** *SUCT1*  - Initialize module yomct1

!     Purpose.
!     --------
!           Initialize module yomct1
!
!**   Interface.
!     ----------
!        *CALL* *SUCT1

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :   Common YOMCT1
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
!      R. El Khatib *Meteo-France*
!      Original : 02-Apr-2015 from SU1YOM.

!     Modifications.
!     --------------
!        R. El Khatib 07-Mar-2016 Pruning of ISP
!      R. El Khatib  15-Sep-2016 Setup of LRFILAF
!      R. El Khatib  04-Jun-2018 refactor suct1 against monio
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMARG   , ONLY : NGRIBFILE
USE YOMCT0   , ONLY : LARPEGEF
USE YOMCT1   , ONLY : N1POS, N1HIS, N1GDI, N1SDI, N1DHP, N1RES, N1XFU, &
 & N1DHFG, N1DHFZ, N1DHFD, N1CFU, LRFILAF, N1SFXHIS, N1MASSCON, LWRSPEC

!     ------------------------------------------------------------------

IMPLICIT NONE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namct1.nam.h"

#include "posnam.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUCT1',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       2.    Initialize YOMCT1.
!              ------------------

!        2.1 Set implicit default values

N1POS=1
N1HIS=1
N1SDI=1
N1RES=1
N1SFXHIS=0
LWRSPEC=.TRUE.
LRFILAF=(NGRIBFILE==0)
N1GDI=1
N1DHP=1
N1DHFG=1
N1DHFZ=1
N1DHFD=1
N1CFU=1
N1XFU=1
N1MASSCON=1

!*       2.3   READ NAMELIST.

CALL POSNAM(NULNAM,'NAMCT1')
READ(NULNAM,NAMCT1)

LRFILAF=LRFILAF.AND.LARPEGEF

!*       2.5   PRINT VALUES

WRITE(UNIT=NULOUT,FMT='('' N1POS = '',I1,'' N1HIS = '',I1,&
 & '' N1GDI = '',I1,'' N1SDI = '',I1,'' N1DHP = '',I1,&
 & '' N1DHFG= '',I1,'' N1DHFZ= '',I1,'' N1DHFD= '',I1,&
 & '' N1CFU = '',I1,'' N1XFU = '',I1,&
 & '' N1RES = '',I1,'' N1SFXHIS = '',I1,'' N1MASSCON = '',I1)')&
 & N1POS,N1HIS,N1GDI,N1SDI,N1DHP,N1DHFG,N1DHFZ,N1DHFD,N1CFU,N1XFU,&
 & N1RES,N1SFXHIS,N1MASSCON 
WRITE(UNIT=NULOUT,FMT='('' LRFILAF = '',L2,'' LWRSPEC = '',L2)')&
 & LRFILAF, LWRSPEC

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUCT1',1,ZHOOK_HANDLE)
END SUBROUTINE SUCT1
