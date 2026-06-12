MODULE SRFROOTFR_MOD
CONTAINS
SUBROUTINE SRFROOTFR(KLEVS,KVTYPES,PSDEPTH,PROOTF)
 
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

! (C) Copyright 2000- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *SRFROOTFR* - Initializing root fraction

!     Purpose.
!     --------
!          Computes the fraction of roots for each soil layer, for each
!     vegetation type in the BATS classification.

!     Interface.
!     ----------
!          *SRFROOTFR* is called from *SUVEG* and *SUSURF*.

!     PARAMETER   DESCRIPTION                                    UNITS
!     ---------   -----------                                    -----
!     INPUT PARAMETERS (INTEGER):
!    *KLEVS*      Number of soil layers
!    *KVTYPES*    Number of vegetation (surface cover) types

!     INPUT PARAMETERS (REAL):
!    *PSDEPTH*    Soil layer depth (dimension KLEVS)

!     OUTPUT PARAMETERS (REAL):
!    *PROOTF*     Root fraction (dimension KLEVS,KVTYPES)

!     Method.
!     -------
!     Straightforward substitution of the ECMWF accumulated layer depths
!     into the integrated root fraction formula used by Canadell et al.
!     1996, Jackson et al. 1996, and adapted to the BATS classification
!     by Zeng et al. 1998.

!     References.
!     ----------

!     Authors:
!     --------
! Canadell, 1996:
! Jackson, 1996:
! Zeng, 1998:

!     Modifications:
!     --------------
!     C. Fischer 00-12-20 Meteo-France recode initialization for prootf
!                         to avoid memory overflow on SUN workstations
!     P.Viterbo       E.C.M.W.F.     12/01/1999.
!     J.F. Estrade *ECMWF* 03-10-01 move in surf vob
!        M.Hamrud      01-Oct-2003 CY28 Cleaning

IMPLICIT NONE

! Declaration of arguments

INTEGER(KIND=JPIM), INTENT(IN)   :: KLEVS
INTEGER(KIND=JPIM), INTENT(IN)   :: KVTYPES

REAL(KIND=JPRB),    INTENT(IN)   :: PSDEPTH(:)

REAL(KIND=JPRB),    INTENT(OUT)  :: PROOTF(:,:)

!      LOCAL VARIABLES

INTEGER(KIND=JPIM), PARAMETER :: JPTYP=20
REAL(KIND=JPRB) :: ZAA(JPTYP),ZBB(JPTYP),ZISD(KLEVS)

INTEGER(KIND=JPIM) :: IFLAG, J, JT, JZ

REAL(KIND=JPRB) :: ZSD
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SRFROOTFR_MOD:SRFROOTFR',0,ZHOOK_HANDLE)
IF (KVTYPES /= JPTYP) THEN
  WRITE(*,*) ' UNEXPECTED NUMBER OF VEGETATION TYPES'
  WRITE(*,*) ' NUMBER OF TYPES ASSUMED BY THIS ROUTINE= ',JPTYP
  WRITE(*,*) ' NUMBER OF TYPES GIVEN AS INPUT= ',JPTYP
ENDIF

! beta values from the 3 papers above (note those papers have depth in
!                                      cm)
! The following order is assumed
!   KTYPE        Description              Low veg/High veg/bare ground/not used
!     1       Crops, Mixed Farming                L
!     2       Short Grass                         L
!     3       Evergreen Needleleaf Trees          H
!     4       Deciduous Needleleaf Trees          H
!     5       Deciduous Broadleaf Trees           H
!     6       Evergreen Broadleaf Trees           H
!     7       Tall Grass                          L
!     8       Desert                              b
!     9       Tundra                              L
!    10       Irrigated Crops                     L
!    11       Semidesert                          L
!    12+      Ice Caps and Glaciers               b
!    13*      Bogs and Marshes                    L
!    14+      Inland Water                        n
!    15+      Ocean                               n
!    16       Evergreen Shrubs                    L
!    17       Deciduous Shrubs                    L
!    18       Mixed Forest                        H
!    19       Interrupted Forest                  H
!    20*      Water and Land Mixtures             L

! Notes: * values not given in Zeng et al. paper
!          13 (Bogs and marshes) exists over Northern latitudes, and
!            it could be initialized as for tundra (root depth 0.5 m).
!            However, it also exists in tropical areas (i.e., the 
!            everglades, west Africa (such as Guine-Bissau)), where
!            a root depth could soon turn the point into a desert-type
!            climate. Values for Type 6 (Evergreen Broadleaf Trees)
!            have been assigned, with a root depth of 3 m. 
!          20 (Water and Land Mixtures) has been assigned as for water
!            (but this is largely academic, because at the time of
!            writing this routine (early 1999) this type was absent from
!            the GLCC data set.
!        + values arbitrarily assigned to 1, to give 0 layer depth

! New table according to Xubing Zeng (adapted to multilayer configuration)
ZAA(1:jptyp)=(/5.558_JPRB,10.739_JPRB, 6.706_JPRB, 7.066_JPRB, 5.990_JPRB,&
 & 7.344_JPRB, 8.235_JPRB, 4.372_JPRB, 8.992_JPRB, 5.558_JPRB,&
 & 4.372_JPRB, 0.0_JPRB    , 7.344_JPRB, 0.0_JPRB    , 0.0_JPRB    ,&
 & 6.326_JPRB, 6.326_JPRB, 4.453_JPRB, 4.453_JPRB, 0.0_JPRB    /)  
ZBB(1:jptyp)=(/2.614_JPRB, 2.608_JPRB, 2.175_JPRB, 1.953_JPRB, 1.955_JPRB,&
 & 1.303_JPRB, 1.627_JPRB, 0.978_JPRB, 8.992_JPRB, 2.614_JPRB,&
 & 0.978_JPRB, 0.0_JPRB    , 1.303_JPRB,  0.0_JPRB   , 0.0_JPRB    ,&
 & 1.567_JPRB, 1.567_JPRB, 1.631_JPRB, 1.631_JPRB, 0.0_JPRB    /)  

ZISD(1)=PSDEPTH(1)
DO J=2,KLEVS
  ZISD(J)=ZISD(J-1)+PSDEPTH(J)
ENDDO
ZSD=ZISD(KLEVS)

DO JT=1,KVTYPES
  IF (ZAA(JT) <= 0.01_JPRB) THEN
!   Initialize "odd" types, for which root depth is not relevant
    PROOTF(1,JT)=0.0_JPRB
    DO J=2,KLEVS
      IF (KLEVS >= J) PROOTF(J,JT)=0.0_JPRB
    ENDDO
    IFLAG=2
  ELSE
! Layer 1
    IFLAG=0
    PROOTF(1,JT)=1.0_JPRB-0.5_JPRB*&
     & (EXP(-ZAA(JT)*ZISD(1))+EXP(-ZBB(JT)*ZISD(1)))   
!     Layers 2 to KLEVS-1
    DO JZ=2,KLEVS-1
      PROOTF(JZ,JT)=0.5_JPRB*(&
       & EXP(-ZAA(JT)*ZISD(JZ-1))+&
       & EXP(-ZBB(JT)*ZISD(JZ-1))-&
       & EXP(-ZAA(JT)*ZISD(JZ  ))-&
       & EXP(-ZBB(JT)*ZISD(JZ  ))&
       & )  
    ENDDO
!   Layer KLEVS
    PROOTF(KLEVS,JT)=1.-SUM(PROOTF(1:KLEVS-1,JT))
  ENDIF
!  write(*,'(A32,2X,I2,4(1X,F5.3))') 'ROOTFRACTION: VEGT/L1/L2/L3/L4: ',&
!    &jt,(prootf(jz,jt),jz=1,4)
ENDDO

! Redefine root depth of "funny" classes as depth of top soil layer

PROOTF(1,8) =1.0_JPRB
PROOTF(1,12)=1.0_JPRB
PROOTF(1,14)=1.0_JPRB
PROOTF(1,15)=1.0_JPRB
PROOTF(1,20)=1.0_JPRB
DO J=2,KLEVS
  IF (KLEVS >= J) THEN
    PROOTF(J,8) =0.0_JPRB
    PROOTF(J,12)=0.0_JPRB
    PROOTF(J,14)=0.0_JPRB
    PROOTF(J,15)=0.0_JPRB
    PROOTF(J,20)=0.0_JPRB
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('SRFROOTFR_MOD:SRFROOTFR',1,ZHOOK_HANDLE)

END SUBROUTINE SRFROOTFR
END MODULE SRFROOTFR_MOD
