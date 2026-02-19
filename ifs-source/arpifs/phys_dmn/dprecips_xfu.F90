! (C) Copyright 1989- Meteo-France.

SUBROUTINE DPRECIPS_XFU(KIDIA,KFDIA,KLON,KDTPREC,KSTATS,PDPRECIPS,PXPTYPE)

!**** *DPRECIPS_XFU*   -  Compute precipitation type diagnostic

!     Purpose.
!     --------
!           Compute precipitation type diagnostics for fullpos 

!**   Interface.
!     ----------
!        *CALL* *DPRECIPS_XFU(...)

!        Explicit arguments :
!        --------------------
!----
! 0D :
!----
! KIDIA, KFDIA : START/END OF HORIZONTAL LOOP  (IST,IEND IN *CPG*).
! KLON : HORIZONTAL DIMENSION                  (ILONMNH IN *APL_AROME*)
! KSTATS : 0 if PTYPE FREQUENT, 2 if PTYPE SEVERE
! ------
! INOUT :
! ------
! PDPRECIPS   : precipitation type diagnostic :
!    0: no precipitation
!    1: rain   / pluie
!    3: freezing rain / pluie verglacante
!    5: dry snow / neige seche 
!    6: wet snow / neige humide
!    7: rain now mixture / pluie et neige melees
!    8: ice pellets/ granules de glace
!    9: graupel   / gresil ou petite grele
!   10: hail      / grele
!   11: drizzle/ bruine
!   12: freezing drizzle / bruine verglacante
!  193: moist snow / neige mouillee
!  201: Pluie intermittente
!  205: Neige sèche intermittente
!  206: Neige humide intermittente
!  207: Pluie et neige mêlées intermittentes
!  213: Neige mouillée intermittente

!                                 
!        Implicit arguments :
!        --------------------
!        COMMON YOMDPRECIPS

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation ARPEGE/AROME

!     Authors.
!     -------
!      I.Etchevers Y. Seity.
!      Original : 2018-09-14

!     Modifications.
!     --------------
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK    ,DR_HOOK,  JPHOOK
USE YOMLUN   , ONLY : NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON
INTEGER(KIND=JPIM),INTENT(IN)    :: KDTPREC
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTATS
REAL(KIND=JPRB)   ,INTENT(IN)    :: PDPRECIPS(KLON,KDTPREC)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PXPTYPE(KLON)

INTEGER (KIND=JPIM):: JLON, JPREC, JTYPE, INB
INTEGER (KIND=JPIM):: IDPRECIPS(KLON,11),INDTOT(KLON)
INTEGER (KIND=JPIM):: ITYPE(KLON)
INTEGER(KIND=JPIM ):: INDPTYPE(11)
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
!     ------------------------------------------------------------------


!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DPRECIPS_XFU',0,ZHOOK_HANDLE)

! Initialisations
IDPRECIPS(:,:)=0_JPIM
INDTOT(:)=0_JPIM
ITYPE(:)=0_JPIM

!ordering important : from less dangerous to most dangerous
INDPTYPE(1)=11
INDPTYPE(2)=1
INDPTYPE(3)=7
INDPTYPE(4)=8
INDPTYPE(5)=9
INDPTYPE(6)=5
INDPTYPE(7)=6
INDPTYPE(8)=193
INDPTYPE(9)=10
INDPTYPE(10)=12
INDPTYPE(11)=3
!

! Number of ptype occurences
DO JLON=KIDIA,KFDIA

  INB=0
  DO JTYPE=1,11
    DO JPREC=1,KDTPREC
      IF (PDPRECIPS(JLON,JPREC)==FLOAT(INDPTYPE(JTYPE))) THEN
        IDPRECIPS(JLON,JTYPE)=IDPRECIPS(JLON,JTYPE)+1_JPIM
      ENDIF
    ENDDO

   ! Most Frequent ptype
    IF (KSTATS==0) THEN
       IF (IDPRECIPS(JLON,JTYPE)>=INB) THEN
         ITYPE(JLON)=INDPTYPE(JTYPE)
         INB=IDPRECIPS(JLON,JTYPE)
       ENDIF
   ! Most severe ptype
    ELSEIF (KSTATS==2) THEN
       IF (IDPRECIPS(JLON,JTYPE)>0_JPIM) THEN
         ITYPE(JLON)=INDPTYPE(JTYPE)
         INB=IDPRECIPS(JLON,JTYPE)
       ENDIF
    ENDIF
    INDTOT(JLON)=INDTOT(JLON)+IDPRECIPS(JLON,JTYPE)
  ENDDO

ENDDO  

DO JLON=KIDIA,KFDIA

!set to 0 when only 1/12 NDTPERIOD is concerned by precipitations (except for hail)
  IF (INDTOT(JLON)< INT(FLOAT(KDTPREC)/12._JPRB).AND.ITYPE(JLON)/=10_JPIM) THEN
     ITYPE(JLON)=0_JPIM
  ENDIF   

! Intermittent character
  IF (INDTOT(JLON)< INT(FLOAT(KDTPREC)*5._JPRB/6._JPRB)) THEN
     IF (ITYPE(JLON)==1_JPIM) ITYPE(JLON)=201_JPRB
     IF (ITYPE(JLON)==7_JPIM) ITYPE(JLON)=207_JPRB
     IF (ITYPE(JLON)==6_JPIM) ITYPE(JLON)=206_JPRB
     IF (ITYPE(JLON)==5_JPIM) ITYPE(JLON)=205_JPRB
     IF (ITYPE(JLON)==193_JPIM) ITYPE(JLON)=213_JPRB
  ENDIF   

  PXPTYPE(JLON)=FLOAT(ITYPE(JLON)) 

ENDDO

IF (LHOOK) CALL DR_HOOK('DPRECIPS_XFU',1,ZHOOK_HANDLE)
END SUBROUTINE DPRECIPS_XFU
