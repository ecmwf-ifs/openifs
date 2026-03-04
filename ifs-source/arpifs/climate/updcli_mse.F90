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

SUBROUTINE UPDCLI_MSE(YDGEOMETRY,YDRIP,KGP,PTEFRCL)

!**** *UPDCLI_MSE*

!     PURPOSE.
!     --------

!     Sends climatolgical data to SURFEX

!**   INTERFACE.
!     ----------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     FAITOU
!     FACILE
!     ARO_PUT_SST

!     AUTHORS.
!     --------


!     MODIFICATIONS.
!     --------------
!      R. El Khatib 17-Aug-2016 remove suarpio
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1     , ONLY : JPIM,   JPRB, JPRD
USE YOMHOOK      , ONLY : LHOOK,  DR_HOOK, JPHOOK
USE YOMLUN       , ONLY : NULOUT
USE YOMRIP       , ONLY : TRIP
USE YOMMP0       , ONLY : MYPROC
USE DISGRID_MOD  , ONLY : DISGRID_SEND, DISGRID_RECV
!     ------------------------------------------------------------------

IMPLICIT NONE

!     ------------------------------------------------------------------
TYPE(GEOMETRY)    ,INTENT(IN)    :: YDGEOMETRY
TYPE(TRIP)        ,INTENT(INOUT) :: YDRIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KGP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTEFRCL
INTEGER(KIND=JPIM) :: JPCCL
PARAMETER(JPCCL=20)
REAL(KIND=JPRB) ::    ZREALG(YDGEOMETRY%YRDIM%NDLON*(YDGEOMETRY%YRDIM%NDGLG+2))
REAL(KIND=JPRB) ::    ZFIELDBUF(YDGEOMETRY%YRGEM%NGPTOTG)
REAL(KIND=JPRB) ::    ZFIELD(YDGEOMETRY%YRGEM%NGPTOT,KGP)
CHARACTER :: CLNOMC*16, CLSTTO*7, CLNOMF*20, CLCST(JPCCL)*20, CLSTTF*6 
INTEGER(KIND=JPIM) :: IREP, INBARI, INDEX, IFIELD, INIMES, INBARP, INULCL1, IZFSST,&
 & IZCT1, IBC1, IINF, INULCL, JNULCL, IREG
LOGICAL :: LLERFA, LLIMST, LLNOMM, LLOP
REAL(KIND=JPRB)   :: ZEPS 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE 

#include "fcttim.func.h"
#include "aro_put_sst.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('UPDCLI_MSE',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NDGLG=>YDDIM%NDGLG, NDLON=>YDDIM%NDLON, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTG=>YDGEM%NGPTOTG, NLOENG=>YDGEM%NLOENG, &
 & RSTATI=>YDRIP%RSTATI)
!     ------------------------------------------------------------------

!*
!     1. SETTING CONSTANT VALUES.
!     ---------------------------

!*    1.1 Arpege files

CLCST(1)   =  'SURFTEMPERATURE '

LLNOMM=.TRUE.
CLSTTO='OLD'
LLERFA=.TRUE.
LLIMST=.FALSE.
INIMES=1
INBARP=14
CLNOMC='Cadre.Clim'
INULCL1=0
LLOP=.FALSE.
IZFSST=12
IZCT1=1

!*    2. READING MONTHLY VALUES.
!     --------------------------

IF (MYPROC == 1) THEN

  !Update SST with LB file instead of clim file.
  IBC1=RSTATI/PTEFRCL
  
  !*    2.1 Opening Arpege files
 
  WRITE(UNIT=CLNOMF,FMT='(''ELSCFHARMALBC'',I3.3)') IBC1
  DO JNULCL=1,99
    INQUIRE(UNIT=JNULCL,OPENED=LLOP)
    IF (.NOT.LLOP) EXIT
  ENDDO
  INULCL1=JNULCL

  CALL FAITOU (IREP,INULCL1,LLNOMM,CLNOMF,CLSTTO,&
   & LLERFA,LLIMST,INIMES,INBARP,INBARI,CLNOMC)  
  ZEPS=1.E-10_JPRB
  IINF=0

  IF ( MYPROC == 1 ) WRITE(NULOUT,*)' SUCESSFULLY OPENED ',TRIM(CLNOMF)

ENDIF

!*    2.4 Indices for horizontal loops

INDEX=NLOENG(1)+1

! Reading and Sending fields

IFIELD=0

!*       2.9 Reading SST

INULCL=INULCL1
IZCT1=IFIELD+IZCT1
IFIELD=IFIELD+1

IF (MYPROC == 1) THEN

  CALL FACILE (IREP,INULCL,CLCST(1)(1:4),1,CLCST(1)(5:),&
   & ZREALG,.FALSE.)  
  
  ZFIELDBUF=ZREALG(INDEX:NGPTOTG+INDEX-1)
  WRITE(NULOUT,'(2X,A,'' READ FROM ARPEGE FILE'',&
   & '' FIELD:'',I3)') CLCST(1),IFIELD  
  CALL DISGRID_SEND(YDGEOMETRY,1,ZFIELDBUF,IFIELD,ZFIELD)
ELSE
! Read distributed data
  CALL DISGRID_RECV(YDGEOMETRY,1,1,ZFIELD,IFIELD)  
ENDIF !MYPROC

CALL FLUSH(NULOUT)

CALL ARO_PUT_SST( NGPTOT, ZFIELD(:,IZCT1) )

!* Close file
IF (MYPROC == 1) THEN
  INQUIRE(UNIT=INULCL1,OPENED=LLOP)
    IF (LLOP) CALL FAIRME (IREG,INULCL1,CLSTTF)
ENDIF

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('UPDCLI_MSE',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCLI_MSE

