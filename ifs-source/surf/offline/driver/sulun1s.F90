SUBROUTINE SULUN1S

USE PARKIND1  ,ONLY : JPIM     ,JPRB, JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK
USE MPL_MODULE, ONLY: MPL_BARRIER, MPL_MYRANK
USE YOMLUN1S , ONLY : NULOUT   ,NULNAM   ,NULFOR   ,NPOSGG,NPOSRES,NPOSCLM ,&
            &NPOSDFO  ,NPOSDBD  ,NPOSDTI1 ,NPOSDTI2 ,NPOSDTI3 ,&
            &NPOSDTI4 ,NPOSDTI5 ,NPOSDTI6 ,NPOSDTI7 ,NPOSDTI8 ,&
            &NPOSDST  ,NPOSDT0  ,NPOSDT1  ,&
            &NPOSDT2  ,NPOSDT3  ,NPOSDT4  ,NPOSDSW  ,NPOSDW0  ,&
            &NPOSDW1  ,NPOSDW2  ,NPOSDW3  ,NPOSDW4  ,NULGP0   ,&
            &NULGPD0  ,NPOSDMO  ,NPOSRC,&
            &NPOSGGL  ,NPOSDTI9 ,& 
            &NPOSOCP  ,NPOSOCD  ,&
            &NPOSCO2  ,NPOSVEG, &
            &RMISS, ZMISS, DMISS, IMISS, CTUNITS,NCTYPE
USE YOMLOG1S , ONLY : LACCUMW  , CFOUT,NCDFTYPE
USE YOEPHY   , ONLY : LEOCML
USE YOMCST   , ONLY : RDAY
USE YOMRIP   , ONLY : NINDAT   ,NSSSSS
USE NETCDF
#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SULUN1S * - Routine to initialize the common YOMLUN1S

!     Purpose.
!     --------
!           Initialize and print the common YOMLUN1S

!**   Interface.
!     ----------
!        *CALL* *SULUN1S

!        Explicit arguments :
!        --------------------
!        none

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Pedro Viterbo and Jean-Francois Mahfouf  *ECMWF*

!     Modifications.
!     --------------
!        Original     : 95-03-08
!        J.F. Mahfouf : 95-03-14 (Logical unit of forcing data)
!        BART vd HURK (KNMI): preparation of NetCDF output (jul 2000)
!        Y. Takaya (ECMWF) : add ocean mixed layer model part
!        E Dutra : June 2014: update netcdf4 interface 
#endif
IMPLICIT NONE

INTEGER(KIND=JPIM) :: MYPROC

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER*2 CLI


#include "fcttim.h"
#include "sucdfres.intfb.h"
#include "supcdf.intfb.h"
#include "sudcdf.intfb.h"

IF (LHOOK) CALL DR_HOOK('SULUN1S',0,ZHOOK_HANDLE)


!        1.    Initialize YOMLUN1S.
!              --------------------

!NULOUT     =  6
NULFOR     =  7
MYPROC     = MPL_MYRANK()
WRITE(NULOUT,*) CFOUT
IF (CFOUT == 'netcdf')THEN
!*     open NetCDF files

  RMISS = 1.E20
  ZMISS = 1.E20
  DMISS = 1.E20
  IMISS = 2147483647  ! maximum allowed integer4
  IF ( NCDFTYPE == 4 ) THEN
    NCTYPE = NF90_NETCDF4
  ELSEIF ( NCDFTYPE == 44 ) THEN
    NCTYPE = NF90_CLASSIC_MODEL + NF90_NETCDF4
  ELSEIF (NCDFTYPE == 3 ) THEN
    NCTYPE = NF90_64BIT_OFFSET  
  ELSEIF (NCDFTYPE == 33 ) THEN
    NCTYPE = NF90_CLOBBER
  ENDIF
  IF( MYPROC == 1 ) THEN
    WRITE(CTUNITS,'(''seconds since '',i4,''-'',i2.2,''-'',i2.2,'' '',i2.2,'':'',i2.2,'':'',i2.2)') &
      NCCAA(NINDAT),NMM(NINDAT),NDD(NINDAT),MOD(NSSSSS,86400)/3600,MOD(NSSSSS,3600)/60,MOD(NSSSSS,60)/60
  ENDIF

  NPOSGG=-1
  NPOSRES=-1
  NPOSCLM=-1
  CALL SUCDFRES
  CALL SUPCDF
  CALL SUDCDF

ELSE

!     ------------------------------------------------------------------

  NPOSGG =20
  NPOSDFO=21
  NPOSDBD=22
  NPOSDT0=23
  NPOSDT1=24
  NPOSDT2=25
  NPOSDT3=26
  NPOSDT4=27
  NPOSDSW=28
  NPOSDW0=29
  NPOSDW1=30
  NPOSDW2=31
  NPOSDW3=32
  NPOSDW4=33
  NPOSDST=34
  NPOSDTI1=35
  NPOSDTI2=36
  NPOSDTI3=37
  NPOSDTI4=38
  NPOSDTI5=39
  NPOSDTI6=40
  NPOSDTI7=41
  NPOSDTI8=42
  NPOSDMO=43
  NPOSRC=44
  NPOSOCP=45
  NPOSOCD=46 
  NULGP0=50
  NULGPD0=51
  
  NPOSGGL=52 ! ENDUTRA 
  NPOSDTI9=53 ! ENDUTRA 

  NPOSCO2=54 !CTESSEL
  NPOSVEG=55  !CTESSEL
  
  IF (LACCUMW) THEN
    CLI=''
  ELSE
    CLI='_i'
  ENDIF
  OPEN(NPOSGG,FILE='o_gg')
  OPEN(NPOSDFO,FILE='o_fo')
  OPEN(NPOSDBD,FILE='o_bd'//CLI)
  OPEN(NPOSDT0,FILE='o_t0'//CLI)
  OPEN(NPOSDT1,FILE='o_t1'//CLI)
  OPEN(NPOSDT2,FILE='o_t2'//CLI)
  OPEN(NPOSDT3,FILE='o_t3'//CLI)
  OPEN(NPOSDT4,FILE='o_t4'//CLI)
  OPEN(NPOSDSW,FILE='o_sw'//CLI)
  OPEN(NPOSDW0,FILE='o_w0'//CLI)
  OPEN(NPOSDW1,FILE='o_w1'//CLI)
  OPEN(NPOSDW2,FILE='o_w2'//CLI)
  OPEN(NPOSDW3,FILE='o_w3'//CLI)
  OPEN(NPOSDW4,FILE='o_w4'//CLI)
  OPEN(NPOSDST,FILE='o_st'//CLI)
  OPEN(NPOSDMO,FILE='o_mo'//CLI)
  OPEN(NPOSDTI1,FILE='o_ti1'//CLI)
  OPEN(NPOSDTI2,FILE='o_ti2'//CLI)
  OPEN(NPOSDTI3,FILE='o_ti3'//CLI)
  OPEN(NPOSDTI4,FILE='o_ti4'//CLI)
  OPEN(NPOSDTI5,FILE='o_ti5'//CLI)
  OPEN(NPOSDTI6,FILE='o_ti6'//CLI)
  OPEN(NPOSDTI7,FILE='o_ti7'//CLI)
  OPEN(NPOSDTI8,FILE='o_ti8'//CLI)
  OPEN(NPOSRC,FILE='o_rc'//CLI)
  
  OPEN(NPOSGGL,FILE='o_gg_lake') ! ENDUTRA
  OPEN(NPOSDTI9,FILE='o_ti9'//CLI)  ! ENDUTRA 

  OPEN(NPOSOCP,FILE='o_ocp'//CLI)  !KPP
  OPEN(NPOSOCD,FILE='o_ocd'//CLI)  !KPP

  OPEN(NPOSCO2,FILE='o_co2'//CLI) !CTESSEL
  OPEN(NPOSVEG,FILE='o_veg'//CLI) !CTESSEL


!     ------------------------------------------------------------------

!        2.    Print YOMLUN1S.
!              ---------------

  WRITE(UNIT=NULOUT,FMT='(&
   &'' NULOUT='',I3,'' NULNAM='',I3,'' NULFOR='',I3,  &
   &'' NPOSGG='',I3,&
   &'' NPOSDFO='',I3,'' NPOSDBD='',I3,'' NPOSDT0='',I3,&
   &'' NPOSDT1='',I3,'' NPOSDT2='',I3,'' NPOSDT3='',I3,&
   &'' NPOSDT4='',I3,'' NPOSDSW='',I3,'' NPOSDW0='',I3,&
   &'' NPOSDW1='',I3,'' NPOSDW2='',I3,'' NPOSDW3='',I3,&
   &'' NPOSDW4='',I3,'' NPOSDST='',I3,'' NPOSDTI1='',I3,&
   &'' NPOSDTI2='',I3,'' NPOSDTI3='',I3,'' NPOSDTI4='',I3,&
   &'' NPOSDTI5='',I3,'' NPOSDTI6='',I3,'' NPOSDTI7='',I3,&
   &'' NPOSDTI8='',I3,'' NULGP0='',I3 ,'' NULGPD0='',I3,&
   &'' NPOSDMO='',I3, '' NPOSGGL='',I3, '' NPOSOCP='',I3,&
   &'' NPOSOCD='',I3,'' NPOSCO2='',I3'' NPOSVEG='',I3)')&
   &NULOUT,NULNAM,NULFOR,NPOSGG &
   &,NPOSDFO,NPOSDBD,NPOSDT0,NPOSDT1,NPOSDT2,NPOSDT3 &
   &,NPOSDT4,NPOSDSW,NPOSDW0,NPOSDW1,NPOSDW2,NPOSDW3 &
   &,NPOSDW4,NPOSDST,NPOSDTI1,NPOSDTI2,NPOSDTI3,NPOSDTI4 &
   &,NPOSDTI5,NPOSDTI6,NPOSDTI7,NPOSDTI8,NULGP0,NULGPD0 &
   &,NPOSDMO,NPOSGGL,NPOSOCP,NPOSOCD,NPOSCO2,NPOSVEG

!     ------------------------------------------------------------------

ENDIF

!     ------------------------------------------------------------------
CALL MPL_BARRIER()

IF (LHOOK) CALL DR_HOOK('SULUN1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE SULUN1S
