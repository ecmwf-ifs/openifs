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

       SUBROUTINE SUVOLC(YDGEOMETRY)
!**** *SUVOLC * -  Setup Volcano emissions based on namelist input

!     Purpose.
!     --------
!           Initialization of YOMVOLCANO, and some prints

!**   Interface.
!     ----------
!        *CALL* *SUVOLC*

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Johannes Flemming   *ECMWF*

!     Modifications.
!     --------------
!      T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!----------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1  , ONLY : JPRB, JPIM
USE YOMHOOK   , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMVOLCANO, ONLY : SEMIVOC, SEMIVOCFLX, SEMIVOCENS, SEMIVOCFLXENS, LVOCMP ,&
 & IVOCGP, IVOCBLK, SVOCLAT, SVOCLON, SLVOC1, SLVOC2, ILVOCJB, NVOCDATES,&
 & IVOCSTART, SLVOCES1, SLVOCES2, LVOCENS  
USE MPL_MODULE, ONLY : MPL_ALLREDUCE 
USE YOMCST    , ONLY : RPI, RA
USE YOMLUN    , ONLY : NULOUT, NULNAM

IMPLICIT NONE

!  INPUT ARGUMENT
TYPE(GEOMETRY), INTENT(IN) :: YDGEOMETRY
!  LOCAL ARRYS
REAL(KIND=JPRB)                   :: ZLON,ZLAT, ZAREA    
REAL(KIND=JPRB)                   :: ZDIST, ZDISTBUF, ZDISTMIN
INTEGER(KIND=JPIM)                :: IK, IKMIN, IGLOMIN
INTEGER(KIND=JPIM)::  I, IT, JSTGLO, ICEND, ISTC, J  

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! 
#include "namvolcano.nam.h"
#include "abor1.intfb.h"
#include "posnam.intfb.h"


IF (LHOOK) CALL DR_HOOK('SUVOLC',0,ZHOOK_HANDLE) 
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV, &
  & YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP, YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB)
ASSOCIATE(NPROMA=>YDDIM%NPROMA, &
 & NGPTOT=>YDGEM%NGPTOT)


! default settings
SLVOC1(:)=-999.9_JPRB
SLVOC2(:)=-999.9_JPRB
LVOCENS=.FALSE.
SEMIVOC(:)=0.0_JPRB
SEMIVOCENS=1000.0_JPRB
IVOCGP=-999
IVOCBLK=-999
LVOCMP = .FALSE. 
SEMIVOCFLX(:)=-999.9
SEMIVOCFLXENS=-999.9
NVOCDATES=0
IVOCSTART(:)=-99999
ILVOCJB(:)=-9

! read namelist
CALL POSNAM(NULNAM,'NAMVOLCANO')
READ(NULNAM,NAMVOLCANO)

! if no namelistinput return  
IF (  NVOCDATES == 0 ) THEN
  IF (LHOOK) CALL DR_HOOK('SUVOLC',1,ZHOOK_HANDLE)
  RETURN
ENDIF
   
! find nearest gridpoint on mpi task       
ZDISTMIN=10.E10
IKMIN=-999 
IGLOMIN=-999 
DO JSTGLO=1,NGPTOT,NPROMA
  ICEND=MIN(NPROMA,NGPTOT-JSTGLO+1)
  ISTC=1
  DO J = ISTC,ICEND
    IK=J+JSTGLO-1  
    ZLON = (360.0_JPRB/(2.0_JPRB*RPI))*YDGSGEOM_NB%GELAM(IK) 
    IF (ZLON > 180.0_JPRB ) ZLON = ZLON - 360.0_JPRB 
    ZLAT = (360.0_JPRB/(2.0_JPRB*RPI))*YDGSGEOM_NB%GELAT(IK)
    ZDIST= SQRT ((SVOCLAT- ZLAT)**2.0 + (COS(YDGSGEOM_NB%GELAT(IK))*(SVOCLON- ZLON))**2.0)
    IF (ZDIST < ZDISTMIN ) THEN 
      ZDISTMIN=ZDIST
      IKMIN=J
      IGLOMIN=JSTGLO
    ENDIF  
  ENDDO
ENDDO  

! find smallest distance of all processes 
ZDISTBUF=ZDISTMIN
CALL MPL_ALLREDUCE(ZDISTBUF ,'min',LDREPROD=.FALSE.,CDSTRING='SUVOLC1')

! the one with the smallest distance
IF ( ZDISTBUF == ZDISTMIN ) THEN
  LVOCMP = .TRUE.  
  IVOCGP=IKMIN  
  IVOCBLK=IGLOMIN
! double check of GAW is already weighted with reduced gaussian grid 
  ZAREA = 4.0_JPRB * RPI * RA * RA * YDGSGEOM_NB%GAW(IVOCGP+IVOCBLK)  
! to be improved ! 
!  ZAREA = 4.0_JPRB * RPI * RA * RA / FLOAT(NGPTOTG)   
! flux in kg/m2s
  WRITE(NULOUT,*) ' VOLCANO EMISSIONS AT COOR ', SVOCLAT, SVOCLON
  WRITE(NULOUT,*) ' VOLCANO EMISSIONS AT I, IBLK ', IVOCGP, IVOCBLK 
  SEMIVOCFLXENS = SEMIVOCENS/ZAREA
  WRITE(NULOUT,*) ' VOLC SOURCE FLUX for ENSEMBLE', SEMIVOCENS, SEMIVOCFLXENS 
  DO IT = 1, NVOCDATES 
    SEMIVOCFLX(IT) = SEMIVOC(IT)/ZAREA
  
    WRITE(NULOUT,*) ' VOLC DATE SOURCE FLUX ', IVOCSTART(IT),SEMIVOC(IT),&
&            SEMIVOCFLX(IT)
  ENDDO
ENDIF

!  check if more than one process has minimum distance
ZDISTBUF = 0.0_JPRB
IF ( LVOCMP)  ZDISTBUF = 1.0 

CALL MPL_ALLREDUCE(ZDISTBUF ,'sum',LDREPROD=.FALSE.,CDSTRING='SUVOLC2')

IF ( ZDISTBUF > 1.0 ) THEN
  WRITE(NULOUT,'(a)') ' Volcano position not clear: change position by tiny increment or'
  WRITE(NULOUT,'(a)') ' Volcano position not clear: revist this code '
  CALL ABOR1(' VOLCANO POSITION ERROR ' )
ENDIF  


! set fixed levels for ensemble in km 
IF ( LVOCENS ) THEN
 SLVOCES2(1)=0.0_JPRB 
 SLVOCES2(2)=1.0_JPRB  
 SLVOCES2(3)=2.0_JPRB  
 SLVOCES2(4)=4.0_JPRB  
 SLVOCES2(5)=6.0_JPRB  
 SLVOCES2(6)=8.0_JPRB  
 SLVOCES2(7)=10.0_JPRB  
 SLVOCES2(8)=12.0_JPRB  
 SLVOCES2(9)=14.0_JPRB  
 SLVOCES2(10)=16.0_JPRB  

 SLVOCES1(1)=1.0_JPRB
 SLVOCES1(2)=2.0_JPRB
 SLVOCES1(3)=4.0_JPRB
 SLVOCES1(4)=6.0_JPRB
 SLVOCES1(5)=8.0_JPRB
 SLVOCES1(6)=10.0_JPRB
 SLVOCES1(7)=12.0_JPRB
 SLVOCES1(8)=14.0_JPRB
 SLVOCES1(9)=16.0_JPRB
 SLVOCES1(10)=18.0_JPRB
ENDIF


END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUVOLC',1,ZHOOK_HANDLE) 

END SUBROUTINE SUVOLC  
