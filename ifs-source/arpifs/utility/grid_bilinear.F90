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

SUBROUTINE GRID_BILINEAR(KLAT_IN,KLON_IN,PLAT_IN,KIN,&
                       & KLAT_OUT,KLON_OUT,PLAT_OUT,KOUT,&
                       & PFLD_IN,PFLD_OUT)  

!    Purpose.
!    --------
!      Bilinear interpolation of grid field from model resolution
!      to grid defined by input parameter. Assumes the whole
!      field is in memory.

!    Explicit arguments:
!    -------------------

!    Input:
!      KLAT_IN  - number of latitude rows of the input grid
!      KLON_IN  - number of longitudes for each row of the input grid
!      PLAT_IN  - latitude (radians) of each row of the input grid
!      KOUT     - size of input array
!      KLAT_OUT - number of latitude rows of the output grid
!      KLON_OUT - number of longitudes for each row of the output grid
!      PLAT_OUT - latitude (radians) of each row of the output grid
!      KOUT     - size of output array
!      PFLD_IN  - array of grid values at input resolution

!    Output:
!      PFLD_OUT - array of interpolated values

!    Author.
!    -------
!      Y.Tremolet

!    Modifications.
!    --------------
!      Original: 18/01/01
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        D.Salmond     01-Jan-2006 OpenMP

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST, ONLY : RPI

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT_IN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON_IN(KLAT_IN) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT_IN(KLAT_IN) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIN 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLAT_OUT 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON_OUT(KLAT_OUT) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PLAT_OUT(KLAT_OUT) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KOUT 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PFLD_IN(KIN) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLD_OUT(KOUT)

INTEGER(KIND=JPIM) :: JLAT, JLON, ILATS, ILATN, ISW, INW, ISE, INE, II, IJLAT
REAL(KIND=JPRB) :: ZLAT, ZLON, ZDLON, ZDLONN, ZDLONS, ZWGTLAT, ZWGTNW, ZWGTSW, &
 & ZLONNW, ZLONSW  
INTEGER(KIND=JPIM) :: IND_IN(KLAT_IN)
INTEGER(KIND=JPIM) :: IND_OUT(KLAT_OUT)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('GRID_BILINEAR',0,ZHOOK_HANDLE)

CALL GSTATS(1945,0)
II=0
DO IJLAT=1,KLAT_IN
  IND_IN(IJLAT)=II
  II=II+KLON_IN(IJLAT)
ENDDO
II=0
DO IJLAT=1,KLAT_OUT
  IND_OUT(IJLAT)=II
  II=II+KLON_OUT(IJLAT)
ENDDO
CALL GSTATS(1945,1)

CALL GSTATS(1472,0)
!$OMP PARALLEL&
!$OMP& PRIVATE(jlat,zlat,ilats,ilatn,zwgtlat,&
!$OMP& zdlon,zdlons,zdlonn,jlon,zlon,&
!$OMP& isw,inw,ise,ine,zlonsw,zlonnw,zwgtsw,zwgtnw)

!$OMP DO SCHEDULE(DYNAMIC,1)
DO JLAT=1,KLAT_OUT
  ZLAT = PLAT_OUT(JLAT)
  IF (ZLAT <= PLAT_IN(KLAT_IN)) THEN
    ILATS = KLAT_IN
    ILATN = KLAT_IN
  ELSE
    ILATS = 1
    DO WHILE (ZLAT < PLAT_IN(ILATS))
      ILATS = ILATS+1
    ENDDO
    ILATN = MAX(1,ILATS-1)
  ENDIF

  IF (ILATN == ILATS) THEN
    ZWGTLAT = 0.
  ELSE
    ZWGTLAT = (ZLAT-PLAT_IN(ILATN))/(PLAT_IN(ILATS)-PLAT_IN(ILATN))
  ENDIF

  ZDLON  = 2.0*RPI/KLON_OUT(JLAT)
  ZDLONS = 2.0*RPI/KLON_IN(ILATS)
  ZDLONN = 2.0*RPI/KLON_IN(ILATN)

  DO JLON=1,KLON_OUT(JLAT)
    ZLON=(JLON-1)*ZDLON

    ISW = 1+MOD(INT(ZLON/ZDLONS), KLON_IN(ILATS))
    INW = 1+MOD(INT(ZLON/ZDLONN), KLON_IN(ILATN))
    ISE = 1+MOD(ISW, KLON_IN(ILATS) )
    INE = 1+MOD(INW, KLON_IN(ILATN) )

    ZLONSW = ZDLONS*(ISW-1)
    ZLONNW = ZDLONN*(INW-1)
    ZWGTSW = (ZLON-ZLONSW)/ZDLONS
    ZWGTNW = (ZLON-ZLONNW)/ZDLONN

    PFLD_OUT(JLON+IND_OUT(JLAT)) = &
     & ZWGTLAT *(      ZWGTSW *PFLD_IN(IND_IN(ILATS)+ISE)   &
     & +(1.0-ZWGTSW)*PFLD_IN(IND_IN(ILATS)+ISW) ) &
     & + (1.0-ZWGTLAT)*(      ZWGTNW *PFLD_IN(IND_IN(ILATN)+INE)   &
     & +(1.0-ZWGTNW)*PFLD_IN(IND_IN(ILATN)+INW) )  
  ENDDO
ENDDO
!$OMP END DO
!$OMP END PARALLEL
CALL GSTATS(1472,1)


IF (LHOOK) CALL DR_HOOK('GRID_BILINEAR',1,ZHOOK_HANDLE)
END SUBROUTINE GRID_BILINEAR
