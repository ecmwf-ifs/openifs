INTERFACE
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
USE, INTRINSIC :: ISO_C_BINDING

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
!    Linus Magnusson   10-09-28 Sea-ice 

!  INTERFACE: 

!    Integers (In):
!      KIDIA   :    Begin point in arrays
!      KFDIA   :    End point in arrays
!      KLON    :    Length of arrays
!      KTILES  :    Number of tiles
!      KTVL    :    Dominant low vegetation type
!      KTVH    :    Dominant high vegetation type

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
!      PTSRF   :    Surface temp. below skin (e.g. Soil or SST)  (K) 
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


!     ------------------------------------------------------------------

END SUBROUTINE SURFSEBS
END INTERFACE
