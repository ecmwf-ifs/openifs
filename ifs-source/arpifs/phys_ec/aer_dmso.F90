! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE AER_DMSO &
  &( YDEAERSRC,KIDIA, KFDIA, KLON , &
  &  PCI  , PDMSI, PLDAY, PLSM, PSKT, PWND, &
  &  PDMSO, PLISS, PTDMS, PODMS &
  & ) 

!*** * AER_DMSO* - SOURCE TERM FOR OCEANIC DMS

!**   INTERFACE.
!     ----------
!          *AER_DMSO* IS CALLED FROM *AER_SRC*.

!* INPUTS:
!  -------
! PCI    (KON)  : fraction sea-ice                    [0-1]  
! PDMSI  (KLON) : DMS concentration                   (nmol l-1)
! PLDAY  (KLON) : "length of the day"                 [0-1]
! PLSM   (KLON) : land-sea mask                       [0-1]
! PSKT   (KLON) : ocean skin temperature              (K)
! PWND   (KLON) : wind strength at 10 m               (m s-1)

!* OUTPUTS:
!  --------
! PDMSO(KLON)   : Flux of D.M.S.               (molec cm-2 s-1)

!     AUTHOR.
!     -------
!        JEAN-JACQUES MORCRETTE  *ECMWF*
!        FROM E.COSME's liss (2002-11-01) 

!     SOURCE.
!     -------
!     LISS & MERLIVAT (The Role of Air-Sea Exchange in Geochemical Cycling
!     P. Buat-Menart, Ed., 113-127     



!     MODIFICATIONS.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE YOMCST    ,ONLY : RTT 
USE YOEAERSRC ,ONLY : TEAERSRC

!-----------------------------------------------------------------------

IMPLICIT NONE

TYPE(TEAERSRC)    ,INTENT(INOUT):: YDEAERSRC
INTEGER(KIND=JPIM),INTENT(IN)  :: KLON 
INTEGER(KIND=JPIM),INTENT(IN)  :: KIDIA 
INTEGER(KIND=JPIM),INTENT(IN)  :: KFDIA 

REAL(KIND=JPRB)   ,INTENT(IN)  :: PCI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PDMSI(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLDAY(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PLSM(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PSKT(KLON)
REAL(KIND=JPRB)   ,INTENT(IN)  :: PWND(KLON)

REAL(KIND=JPRB)   ,INTENT(OUT) :: PDMSO(KLON)    
REAL(KIND=JPRB)   ,INTENT(OUT) :: PLISS(KLON)
REAL(KIND=JPRB)   ,INTENT(OUT) :: PTDMS(KLON) 
REAL(KIND=JPRB)   ,INTENT(OUT) :: PODMS(KLON)

!-----------------------------------------------------------------------

INTEGER(KIND=JPIM) :: JL
REAL(KIND=JPRB) :: ZDMSSRC(KLON)
REAL(KIND=JPRB) :: ZFCI, ZLISS, ZSCHMIDT, ZSCHMIDT2, ZSO2MSS, Z_S_SO2
REAL(KIND=JPRB) :: ZTSCHM, Z1SCHMDT, Z1SCHMD2, ZDMSI, ZUNIT
  
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('AER_DMSO',0,ZHOOK_HANDLE)
ASSOCIATE(NPIST=>YDEAERSRC%NPIST, RDMSMIN=>YDEAERSRC%RDMSMIN)
!-----------------------------------------------------------------------

ZSO2MSS = 64.058E-03_JPRB
Z_S_SO2 = 0.5_JPRB
ZUNIT = 1.E-10_JPRB/36._JPRB

DO JL=KIDIA,KFDIA
  PDMSO(JL)=0._JPRB
  PODMS(JL)=0._JPRB
  PTDMS(JL)=0._JPRB
  PLISS(JL)=0._JPRB

  ZLISS = 0._JPRB

  IF (PLSM(JL) <= 0.01_JPRB) THEN
!-- Schmidt's mumber (NB: two variants)
    ZTSCHM = PSKT(JL)-RTT
    ZTSCHM = MAX(ZTSCHM, 35._JPRB)

    ZSCHMIDT=2674._JPRB -ZTSCHM* (147.12_JPRB -ZTSCHM* (3.726_JPRB -ZTSCHM* 0.038_JPRB ))
    Z1SCHMDT=600._JPRB/ZSCHMIDT
    
    ZSCHMIDT2=3652.047271_JPRB&
      & +ZTSCHM*( -246.99_JPRB +ZTSCHM*( 8.536397_JPRB -ZTSCHM*0.124397_JPRB) )
    Z1SCHMD2=600._JPRB/ZSCHMIDT2

    IF (NPIST == 1) THEN                           ! Liss & Merlivat, 1986
!-- Computing transfer speed (a la Curran & Jones, 2000, JGR 105, D16, 20451)
      IF (PWND(JL) <= 3.6_JPRB) THEN
        ZLISS = 0.17_JPRB*PWND(JL)
        ZLISS = ZLISS*Z1SCHMDT**0.66667_JPRB
      ELSEIF (PWND(JL) <= 13.0_JPRB) THEN
        ZLISS = 2.85_JPRB *PWND(JL) - 9.65_JPRB
        ZLISS = ZLISS*SQRT(MAX(0._JPRB,Z1SCHMDT))
      ELSE
        ZLISS = 5.9_JPRB*PWND(JL) - 49.3_JPRB
        ZLISS = ZLISS*SQRT(MAX(0._JPRB,Z1SCHMDT))
      ENDIF

    ELSEIF (NPIST == 2) THEN                        ! Wanninkhof, 1992
       ZLISS = 0.31_JPRB*PWND(JL)*PWND(JL)
       ZLISS = ZLISS*SQRT(Z1SCHMD2)

    ELSEIF (NPIST == 3) THEN                        ! Nightingale, 2000
      ZLISS = PWND(JL) * (0.333_JPRB + PWND(JL)*0.222_JPRB)
      ZLISS = ZLISS*SQRT(Z1SCHMD2)
    ENDIF

!-- scaling down if sea-ice cover from 60 to 100%
    ZFCI = 1.0_JPRB
    IF (PCI(JL) > 0.60_JPRB) THEN
      ZFCI = 1.0_JPRB - (PCI(JL)-0.60_JPRB)/0.40_JPRB
    ENDIF
    ZLISS = ZLISS*ZFCI
    PLISS(JL)=ZLISS

!-- the source (in principle related to chlorophyll content) is 
!   simply linked to surface temperature, and the flux emitted is 
!   linked to the transfer speed
    PTDMS(JL)=0._JPRB
    IF (PSKT(JL) <= 298._JPRB) THEN
      PTDMS(JL)=1._JPRB
    ENDIF 
    ZDMSSRC(JL)= RDMSMIN*PTDMS(JL)*ZLISS*PLDAY(JL)
    PDMSO(JL) = ZDMSSRC(JL) *ZSO2MSS *Z_S_SO2

    ZDMSI = MAX(0.2_JPRB, PDMSI(JL))
    PODMS(JL) = ZLISS *ZDMSI *ZSO2MSS *Z_S_SO2 *ZUNIT
  ENDIF
ENDDO
!IF (NSTEP == YRRIP%NSTART) THEN
!  WRITE(NULOUT,FMT='("dms-related DMS ",I4,F8.3,9E12.5)') KIDIA,PLSM(KIDIA),PDMSO(KIDIA),ZDMSSRC(KIDIA),&
!      & PSKT(KIDIA),PCI(KIDIA),PTDMS(KIDIA),PWND(KIDIA),PLISS(KIDIA),PLDAY(KIDIA),PODMS(KIDIA)
!ENDIF



! - DMS flux
!    ZCDMS(JL)=MAX(0.2_JPRB,PCDMS(JL))         ! correction on input source file   
!    PFLXDMS(JL)=ZLISS*ZCDMS(JL)               ! in nmol*cm dm-3 h-1
!    PFLXDMS(JL)=PFLXDMS(JL)/3600._JPRB        ! in nmol*cm dm-3 s-1
!    PFLXDMS(JL)=PFLXDMS(JL)*RNAVO*1.E-09_JPRB ! in molec*cm dm-3 s-1
!    PFLXDMS(JL)=PFLXDMS(JL)*0.001_JPRB        ! in molec cm-2 s-1
!      
!  ENDIF
!
!ENDDO

!-----------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('AER_DMSO',1,ZHOOK_HANDLE)
END SUBROUTINE AER_DMSO
