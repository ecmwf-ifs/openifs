SUBROUTINE SURFSEBS   (YDSURF,KIDIA,KFDIA,KLON,KLEVSN,KTILES,KTVL,KTVH,&
 & PTMST,PSSKM1M,PTSKM1M,PQSKM1M,PDQSDT,PRHOCHU,PRHOCQU,&
 & PALPHAL,PALPHAS,PSSRFL,PFRTI,PTSRF,&
 & PSNS,PRSN,PHLICE, & 
 & PSLRFL,PTSKRAD,PEMIS,PASL,PBSL,PAQL,PBQL,&
 & PTHKICE,PSNTICE,&
 !out
 & PJS,PJQ,PSSK,PTSK,PSSH,PSLH,PSTR,PG0,&
 & PSL,PQL,&
 & LNEMOLIMTHK)  

USE PARKIND1, ONLY : JPIM, JPRB
USE ISO_C_BINDING

!ifndef INTERFACE

USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOS_SURF, ONLY : TSURF, GET_SURF

USE ABORT_SURF_MOD
USE SURFSEBS_CTL_MOD

!endif INTERFACE

! (C) Copyright 2003- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!------------------------------------------------------------------------

!  PURPOSE:
!    Routine SURFSEB computes surface energy balance and skin temperature 
!    for each tile. 

!  SURFSEB is called by VDFDIFH

!  METHOD:
!    A linear relation between lowest model level dry static 
!    energy and moisture and their fluxes is specified as input. 
!    The surface energy balance equation is used to eliminate 
!    the skin temperature as in the derivation of the 
!    Penmann-Monteith equation. 

!    The routine can also be used in stand alone simulations by
!    putting PASL and PAQL to zero and by specifying for PBSL and PBQL 
!    the forcing with dry static energy and specific humidity. 

!  AUTHOR:
!    A. Beljaars       ECMWF April 2003   

!  REVISION HISTORY:
!    J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!    E. Dutra/G. Balsamo  01-05-08 lake tile
!    Linus Magnusson   10-09-28 Sea-ice 
!    G. Balsamo/A. Beljaars 09-08-2013 snow scheme stability fix
!      I. Sandu    24-02-2014  Lambda skin values by vegetation type instead of tile

!  INTERFACE: 

!    Integers (In):
!      KIDIA   :    Begin point in arrays
!      KFDIA   :    End point in arrays
!      KLON    :    Length of arrays
!      KTILES  :    Number of tiles
!      KTVL    :    Low vegetation type 
!      KTVH    :    High vegetation type

!    Reals with tile index (In): 
!      PTMST   :    Time Step
!      PSSKM1M :    Dry static energy of skin at T-1           (J/kg)
!      PTSKM1M :    Skin temperature at T-1                    (K)
!      PQSKM1M :    Saturation specific humidity at PTSKM1M    (kg/kg)
!      PDQSDT  :    dqsat/dT at PTSKM1M                        (kg/kg K)
!      PRHOCHU :    Rho*Ch*|U|                                 (kg/m2s)
!      PRHOCQU :    Rho*Cq*|U|                                 (kg/m2s)
!      PALPHAL :    multiplier of ql in moisture flux eq.      (-)
!      PALPHAS :    multiplier of qs in moisture flux eq.      (-)
!      PSSRFL  :    Net short wave radiation at the surface    (W/m2)
!      PFRTI   :    Fraction of surface area of each tile      (-)
!      PTSRF   :    Surface temp. below skin (e.g. Soil or SST)(K) 
!      PSNS    :    Snow mass per unit area                    (kg/m2)
!      PSNS    :    Snow density                               (kg/m3)
!      PHLICE  :    Lake ice thickness                         (m) 

!    Reals independent of tiles (In):
!      PSLRFL  :    Net long wave radiation at the surface     (W/m2) 
!      PTSKRAD :    Mean skin temp. at radiation time level    (K)
!      PEMIS   :    Surface emissivity                         (-)
!      PASL    :    Asl in Sl=Asl*Js+Bsl                       (m2s/kg)
!      PBSL    :    Bsl in Sl=Asl*Js+Bsl                       (J/kg)
!      PAQL    :    Aql in Ql=Aql*Jq+Bql                       (m2s/kg)
!      PBQL    :    Bql in Ql=Aql*Jq+Bql                       (kg/kg)
!      PTHKICE :    Sea ice thickness                          (m)
!      PSNTICE :    Thickness of snow layer on sea ice         (m)

!    Reals with tile index (Out):
!      PJS     :    Flux of dry static energy                  (W/m2)
!      PQS     :    Moisture flux                              (kg/m2s)
!      PSSK    :    New dry static energy of skin              (J/kg)
!      PTSK    :    New skin temperature                       (K)
!      PSSH    :    Surface sensible heat flux                 (W/m2)
!      PSLH    :    Surface latent heat flux                   (W/m2)
!      PSTR    :    Surface net thermal radiation              (W/m2)
!      PG0     :    Surface ground heat flux (solar radiation  (W/m2)
!                   leakage is not included in this term)

!    Reals independent of tiles (Out):
!      PSL     :    New lowest model level dry static energy   (J/kg)
!      PQL     :    New lowest model level specific humidity   (kg/kg)

!  DOCUMENTATION:
!    See Physics Volume of IFS documentation
!    This routine uses the method suggested by Polcher and Best
!    (the basic idea is to start with a linear relation between 
!     the lowest model level varibles and their fluxes, which is 
!     obtained after the downward elimination of the vertical 
!     diffusion tridiagonal matrix). 

!------------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR)       ,INTENT(IN)    :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVSN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVL(KLON)
INTEGER(KIND=JPIM),INTENT(IN)    :: KTVH(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTMST

REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSKM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQSKM1M(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDQSDT(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHOCHU(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRHOCQU(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHAL(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALPHAS(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSSRFL(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSRF(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSLRFL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTSKRAD(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBSL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAQL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PBQL(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNS(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRSN(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTHKICE(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSNTICE(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PJS(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PJQ(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSK(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTSK(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSSH(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSLH(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSTR(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PG0(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSL(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PQL(:) 
LOGICAL           ,INTENT(IN)    :: LNEMOLIMTHK

!ifndef INTERFACE

TYPE(TSURF), POINTER :: YSURF
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SURFSEB',0,ZHOOK_HANDLE)

YSURF => GET_SURF(YDSURF)

IF(UBOUND(PSLRFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PSLRFL TOO SHORT!')
ENDIF

IF(UBOUND(PTSKRAD,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PTSKRAD TOO SHORT!')
ENDIF

IF(UBOUND(PEMIS,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PEMIS TOO SHORT!')
ENDIF

IF(UBOUND(PASL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PASL TOO SHORT!')
ENDIF

IF(UBOUND(PBSL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PBSL TOO SHORT!')
ENDIF

IF(UBOUND(PAQL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PAQL TOO SHORT!')
ENDIF

IF(UBOUND(PBQL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB: PBQL TOO SHORT!')
ENDIF

IF(UBOUND(PSSKM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PSL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSL TOO SHORT!')
ENDIF

IF(UBOUND(PQL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PQL TOO SHORT!')
ENDIF

IF(UBOUND(PSSKM1M,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PTSKM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PTSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PTSKM1M,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PTSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PQSKM1M,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PQSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PQSKM1M,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PQSKM1M TOO SHORT!')
ENDIF

IF(UBOUND(PDQSDT,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PDQSDT TOO SHORT!')
ENDIF

IF(UBOUND(PDQSDT,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PDQSDT TOO SHORT!')
ENDIF

IF(UBOUND(PRHOCHU,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PRHOCHU TOO SHORT!')
ENDIF

IF(UBOUND(PRHOCHU,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PRHOCHU TOO SHORT!')
ENDIF

IF(UBOUND(PRHOCQU,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PRHOCQU TOO SHORT!')
ENDIF

IF(UBOUND(PRHOCQU,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PRHOCQU TOO SHORT!')
ENDIF

IF(UBOUND(PALPHAL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PALPHAL TOO SHORT!')
ENDIF

IF(UBOUND(PALPHAL,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PALPHAL TOO SHORT!')
ENDIF

IF(UBOUND(PALPHAS,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PALPHAS TOO SHORT!')
ENDIF

IF(UBOUND(PALPHAS,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PALPHAS TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFL,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSSRFL TOO SHORT!')
ENDIF

IF(UBOUND(PSSRFL,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSSRFL TOO SHORT!')
ENDIF

IF(UBOUND(PFRTI,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PFRTI,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PFRTI TOO SHORT!')
ENDIF

IF(UBOUND(PTSRF,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PTSRF TOO SHORT!')
ENDIF

IF(UBOUND(PTSRF,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PTSRF TOO SHORT!')
ENDIF

IF(UBOUND(PTHKICE,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PTHKICE TOO SHORT!')
ENDIF

IF(UBOUND(PSNTICE,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSNTICE TOO SHORT!')
ENDIF

IF(UBOUND(PJS,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PJS TOO SHORT!')
ENDIF

IF(UBOUND(PJS,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PJS TOO SHORT!')
ENDIF

IF(UBOUND(PJQ,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PJQ TOO SHORT!')
ENDIF

IF(UBOUND(PJQ,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PJQ TOO SHORT!')
ENDIF

IF(UBOUND(PSSK,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSSK TOO SHORT!')
ENDIF

IF(UBOUND(PSSK,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSSK TOO SHORT!')
ENDIF

IF(UBOUND(PTSK,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PTSK TOO SHORT!')
ENDIF

IF(UBOUND(PTSK,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PTSK TOO SHORT!')
ENDIF

IF(UBOUND(PSSH,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSSH TOO SHORT!')
ENDIF

IF(UBOUND(PSSH,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSSH TOO SHORT!')
ENDIF

IF(UBOUND(PSLH,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSLH TOO SHORT!')
ENDIF

IF(UBOUND(PSLH,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSLH TOO SHORT!')
ENDIF

IF(UBOUND(PSTR,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSTR TOO SHORT!')
ENDIF

IF(UBOUND(PSTR,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PSTR TOO SHORT!')
ENDIF

IF(UBOUND(PG0,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PG0 TOO SHORT!')
ENDIF

IF(UBOUND(PG0,2) < KTILES) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PG0 TOO SHORT!')
ENDIF

IF(UBOUND(PHLICE,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PHLICE TOO SHORT!')
ENDIF

IF(UBOUND(PSNS,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSNS TOO SHORT!')
ENDIF

IF(UBOUND(PSNS,2) < KLEVSN) THEN
  CALL ABORT_SURF('SURFSEB:: FIRST DIMENSION OF PSNS TOO SHORT!')
ENDIF

IF(UBOUND(PRSN,1) < KLON) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PRSN TOO SHORT!')
ENDIF
IF(UBOUND(PRSN,2) < KLEVSN) THEN
  CALL ABORT_SURF('SURFSEB:: SECOND DIMENSION OF PRSN TOO SHORT!')
ENDIF

CALL SURFSEBS_CTL      (KIDIA,KFDIA,KLON,KTILES,KTVL,KTVH,&
 & PTMST,PSSKM1M,PTSKM1M,PQSKM1M,PDQSDT,PRHOCHU,PRHOCQU,&
 & PALPHAL,PALPHAS,PSSRFL,PFRTI,PTSRF,&
 & PSNS,PRSN,PHLICE, & 
 & PSLRFL,PTSKRAD,PEMIS,PASL,PBSL,PAQL,PBQL,&
 & PTHKICE,PSNTICE,&
 & YSURF%YCST,YSURF%YEXC,YSURF%YVEG,YSURF%YFLAKE,YSURF%YSOIL,&
 & PJS,PJQ,PSSK,PTSK,PSSH,PSLH,PSTR,PG0,&
 & PSL,PQL,&
 & LNEMOLIMTHK)  
IF (LHOOK) CALL DR_HOOK('SURFSEB',1,ZHOOK_HANDLE)

!endif INTERFACE

!     ------------------------------------------------------------------

END SUBROUTINE SURFSEBS
