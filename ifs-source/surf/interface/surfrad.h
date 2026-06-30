INTERFACE
SUBROUTINE SURFRAD    (YDSURF,KDD,KMM,KMON,KSECO,&
 & KIDIA,KFDIA,KLON,KTILES,KSW,KLW,&
 & LDNH,&
 & PALBF,PALBICEF,PTVH,&
 & PALCOEFF, PCURC,PCVH,&
 & PASN,PMU0,PTS,PWND,&
 & PWS1,KSOTY,PFRTI,PHLICE,PTLICE,&  
 & PALBD,PALBP,PALB,&
 & PSPECTRALEMISS,PEMIT,&
 & PALBTI,PCCNL,PCCNO,&
 & LNEMOLIMALB)

USE PARKIND1, ONLY : JPIM, JPRB
USE, INTRINSIC :: ISO_C_BINDING

! (C) Copyright 1991- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SURFRAD - COMPUTES RADIATIVE PROPERTIES OF SURFACE

!     PURPOSE.
!     --------

!**   INTERFACE.
!     ----------
!        CALL *SURFRAD* FROM *CALLPAR*

!        EXPLICIT ARGUMENTS :
!        --------------------
!     ==== INPUTS ===
! KDD,KMM : final date
! KMON   : length of the 12 months
! KSECO  : Number of seconds since 00 UTC on the initial day of the forecast
! LDNH   : LOGICAL  : .TRUE. FOR Northern Hemisphere
! PALBF  : REAL     : FIXED BACKGROUND SURFACE SHORTWAVE ALBEDO
! PALBICEF REAL     : FIXED SEA-ICE ALBEDO (FROM COUPLER)
! PALCOEFF : REAL   : MODIS albedo coefficients
! PCURC  : REAL     : URBAN COVER (PASSIVE)
! PCVH   : REAL     : HIGH VEGETATION COVER
! PASN   : REAL     : ALBEDO OF EXPOSED SNOW (TYPE 5)
! PMU0   : REAL     : COSINE OF SOLAR ZENITH ANGLE
! PTS    : REAL     : SURFACE TEMPERATURE
! PWND   : REAL     : WIND INTENSITY AT LOWEST LEVEL
! PWS1   : REAL     : TOP LAYER SOIL MOISTURE CONTENT
! KSOTY  : INTEGER  : SOIL TYPE                           (1-7)
! PFRTI  : REAL     : TILE FRACTIONS                      (0-1)
!            1 : WATER                  5 : SNOW ON LOW-VEG+BARE-SOIL
!            2 : ICE                    6 : DRY SNOW-FREE HIGH-VEG
!            3 : WET SKIN               7 : SNOW UNDER HIGH-VEG
!            4 : DRY SNOW-FREE LOW-VEG  8 : BARE SOIL
! PHLICE : REAL     : LAKE ICE THICKNESS (m)
! PTLICE : REAL     : LAKE ICE TEMPERATURE (K)

!     ==== OUTPUTS ===
! PALBD  : REAL     : SURFACE ALBEDO FOR DIFFUSE RADIATION
! PALBP  : REAL     : SURFACE ALBEDO FOR PARALLEL RADIATION
! PALB   : REAL     : AVERAGE SW ALBEDO (DIAGNOSTIC ONLY)
! PSPECTRALEMISS : REAL : SURFACE LONGWAVE SPECTRAL EMISSIVITY
! PEMIT  : REAL     : SURFACE LONGWAVE EMISSIVITY
! PALBTI : REAL     : BROADBAND ALBEDO FOR TILE FRACTIONS
! PCCNL  : REAL     : CCN CONCENTRATION OVER LAND
! PCCNO  : REAL     : CCN CONCENTRATION OVER OCEAN

!     ==== OUTPUTS ===

!        IMPLICIT ARGUMENTS :   NONE
!        --------------------

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCE.
!     ----------

!        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE "I.F.S"

!     AUTHOR.
!     -------
!     J.-J. MORCRETTE  E.C.M.W.F.    91/03/15

!     MODIFICATIONS.
!     --------------
!     J.-J. MORCRETTE  ECMWF 94/11/15  DIRECT/DIFFUSE ALBEDOS
!     J.-J. MORCRETTE 96/06/07  moisture dep. emissiv. / spectral alb.  
!     PJANSSEN/JJMORCRETTE ECMWF     96/11/07  WIND DEPENDENT SEA ALBEDO
!     PViterbo         ECMWF 99/03/03  Albedo for tile fractions
!     J.-J. Morcrette ECMWF 00/10/24 Spectral albedo for all surfaces
!     JJMorcrette     01-10-08  CCNs concentration over ocean
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!     M.Hamrud      01-Oct-2003 CY28 Cleaning
!     G.Balsamo     03-Jul-2006 Add soil type
!     JJMorcrette ECMWF 2006-05-10 MODIS albedo
!     E.Dutra/G.Balsamo 2008-05-01 Add lake tile 
!     Linus Magnusson   10-09-28 Sea-ice 
!-----------------------------------------------------------------------

IMPLICIT NONE

! Declaration of arguments

TYPE(C_PTR)       ,INTENT(IN)    :: YDSURF
INTEGER(KIND=JPIM),INTENT(IN)    :: KDD 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KMON(:) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSECO
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)    :: KTILES 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSW 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLW
INTEGER(KIND=JPIM),INTENT(IN)    :: KSOTY(:) 
LOGICAL           ,INTENT(IN)    :: LDNH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBF(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALBICEF(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTVH(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PALCOEFF(:,:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCURC(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCVH(:)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PASN(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PMU0(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWND(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PWS1(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFRTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PHLICE(:) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTLICE(:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALBD(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALBP(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALB(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PSPECTRALEMISS(:,:)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PALBTI(:,:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCNL(:) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PCCNO(:) 
LOGICAL           ,INTENT(IN)    :: LNEMOLIMALB


!     ------------------------------------------------------------------

END SUBROUTINE SURFRAD
END INTERFACE
