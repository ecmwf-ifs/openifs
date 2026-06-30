SUBROUTINE CNTEND
USE PARKIND1  ,ONLY : JPIM     ,JPRB

#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *CNTEND*  - Closes netCDF datafiles.

!     Purpose.
!     --------
!          Data file closure after finishing the run.

!***  Interface.
!     ----------
!        *CALL* *CNTEND

!        Explicit arguments :
!        --------------------
!        None

!        Implicit arguments :
!        --------------------
!        None

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd Hurk *KNMI*

!     Modifications.
!     --------------
!        Original : 14/7/2000
!
!        E. DUTRA - New outuput file for lakes  04/07/2008
!        E. Dutra : Coupling with cama-flood May 2019
!     ------------------------------------------------------------------
#endif
!     ------------------------------------------------------------------

USE YOMLUN1S , ONLY :  NPOSGG  ,NPOSEFL   ,NPOSWAT    ,NPOSRES , NPOSRESO, &  !KPP
                     & NPOSSUS ,NPOSSUB   ,NPOSEVA    ,NPOSCLD, NPOSCLM, &
                     & NPOSLKE ,NPOSOCD   ,NPOSOCP, &  !KPP
                     & NPOSCO2  ,NPOSVEG  ,NPOSEXT  ,NPOSBIO, &
                     & NPOSTIL  ,NPOSVTY  ,NPOSD2M,NPOSGGD

USE YOMLOG1S , ONLY : CFOUT &
           &,LWREFL   ,LWRWAT   ,LWRSUB   ,LWRSUS   ,LWREVA &
           &,LWRCLD   ,LWRGG    ,LWRCLM   ,LWRLKE   &
           &,LWROCD   ,LWROCP   ,LWROCR   & !KPP
           &,LWRCO2 ,LWRVEG ,LWREXT ,LWRBIO ,LWRTIL ,LWRVTY,LWRD2M,LWRGGD
USE YOEPHY,    ONLY : LECMF1WAY
USE CMF_DRV_CONTROL_MOD,     ONLY: CMF_DRV_END  

USE MPL_MODULE
IMPLICIT NONE
#include "netcdf.inc"

INTEGER(KIND=JPIM), PARAMETER :: JPNCDF=21   !KPP JPNCDF=9+3ocean+1lake+6ctessel
INTEGER(KIND=JPIM) :: IPOS(JPNCDF),J
INTEGER NPOS,IERR
LOGICAL LPOS(JPNCDF)
INTEGER(KIND=JPIM) :: MYPROC, NPROC

MYPROC = MPL_MYRANK()
NPROC  = MPL_NPROC()

IF(CFOUT == 'netcdf')THEN
  IPOS(1)=NPOSGG
  IPOS(2)=NPOSEFL
  IPOS(3)=NPOSWAT
  IPOS(4)=NPOSRES
  IPOS(5)=NPOSSUS
  IPOS(6)=NPOSSUB
  IPOS(7)=NPOSEVA
  IPOS(8)=NPOSCLD
  IPOS(9)=NPOSCLM
  IPOS(10)=NPOSLKE
  IPOS(11)=NPOSOCD !KPP
  IPOS(12)=NPOSOCP !KPP
  IPOS(13)=NPOSRESO!KPP
  IPOS(14)=NPOSCO2 !CTESSEL
  IPOS(15)=NPOSBIO !CTESSEL
  IPOS(16)=NPOSVEG !CTESSEL
  IPOS(17)=NPOSEXT !CTESSEL
  IPOS(18)=NPOSTIL !CTESSEL
  IPOS(19)=NPOSVTY !CTESSEL
  IPOS(20)=NPOSD2M 
  IPOS(21)=NPOSGGD 


  LPOS(1)=LWRGG
  LPOS(2)=LWREFL
  LPOS(3)=LWRWAT
  LPOS(4)=.TRUE.
  LPOS(5)=LWRSUS
  LPOS(6)=LWRSUB
  LPOS(7)=LWREVA
  LPOS(8)=LWRCLD
  LPOS(9)=LWRCLM
  LPOS(10)=LWRLKE
  LPOS(11)=LWROCD !KPP
  LPOS(12)=LWROCP !KPP
  LPOS(13)=LWROCR !KPP
  LPOS(14)=LWRCO2 !CTESSEL
  LPOS(15)=LWRBIO !CTESSEL
  LPOS(16)=LWRVEG !CTESSEL
  LPOS(17)=LWREXT !CTESSEL
  LPOS(18)=LWRTIL !CTESSEL
  LPOS(19)=LWRVTY !CTESSEL
  LPOS(20)=LWRD2M
  LPOS(21)=LWRGGD

  IF( MYPROC == 1 ) THEN
    DO J=1,JPNCDF
      NPOS=IPOS(J)
      IF(LPOS(J))CALL NCCLOS(NPOS,IERR)
    ENDDO
  ENDIF
ENDIF

IF (LECMF1WAY) THEN
  !* Close Cama-flood output files 
  IF( MYPROC == 1 ) THEN
    CALL CMF_DRV_END()
  ENDIF 
ENDIF 

RETURN
END SUBROUTINE CNTEND
