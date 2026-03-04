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

#ifdef RS6K
@PROCESS NOCHECK
#endif
SUBROUTINE SLEXTPOL(YDDIM,YDSL, &
 & KFLDSLB,KFIXFLD,KTYP,PARSLB,PB1A,KMASK2)

!****-------------------------------------------------------------------
!**** *SLEXTPOL* - Communication for horizontal interpolations: 
!                  Special Pole treatment.

!                  This routine is called only if distributed memory.
!                  Computations are DM-local.
!                  This routine assumes that there are exactly three 
!                  extra-longitudes: jlon=0, jlon=ndlon+1, jlon=ndlon+2.
!****-------------------------------------------------------------------
!     Purpose.
!     --------
!**   Interface.
!     ----------
!        *CALL* *SLEXTPOL(...)*

!        Explicit arguments :
!        --------------------
!        INPUT:
!        - YDSL       : SL_STRUCT definition
!        - KFLDSLB    : number of fields.
!        - KFIXFLD    : description of 'fixed' fields (to be processed)
!        - KTYP       : type of communication.
!                       KTYP=1: equivalent of the CY35 version of SLEXTPOL
!                       KTYP=2: equivalent of the CY35 version of SLEXTPOL2
!                       KTYP=3: equivalent of the CY35 version of SLEXTPOL1
!                               if KMASK2 is provided.
!                               equivalent of the CY35 version of SLEXTPOL1A
!        - PARSLB     : parity for extra-polar calculations.

!        INPUT/OUTPUT:
!        - PB1A       : fields to be computed at the extra-polar latitudes.

!        OPTIONAL OUTPUT:
!        At least one of these masks should be present if KTYP=3.
!        - KMASK2     : Mask2 for on demand-comms.

!        Implicit arguments :
!        --------------------

!     Externals.
!     ----------

!     Method.
!     -------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS
!        Documentation about distributed memory.

!     Author.
!     -------
!      Lars ISAKSEN (ECMWF), OCTOBER 1995.

!     Modifications.
!     --------------
!      02-10-01  G.Mozdzynski     support for radiation on-demand comms
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      K.Yessad  (Oct 2008) Unified version
!      G. Mozdzynski (Jan 2011): OOPS cleaning, use of derived type SL_STRUCT
!      G. Mozdzynski (Aug 2011): support higher order interpolation
!      G. Mozdzynski (May 2012): further cleaning
!     ------------------------------------------------------------------

USE YOMDIM   , ONLY : TDIM
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN_IFSAUX, ONLY : NULOUT
USE EINT_MOD , ONLY : SL_STRUCT

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIM)        ,INTENT(IN)    :: YDDIM
TYPE(SL_STRUCT)   ,INTENT(IN)    :: YDSL
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSLB 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIXFLD(2)
INTEGER(KIND=JPIM),INTENT(IN)    :: KTYP
REAL(KIND=JPRB)   ,INTENT(IN)    :: PARSLB(KFLDSLB) 
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PB1A(YDSL%NASLB1,KFLDSLB,1) 
INTEGER(KIND=JPIM),OPTIONAL,INTENT(OUT) :: KMASK2(YDSL%NASLB1+YDDIM%NSTENCILWIDE*2)

!     ------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IEND, IGLGLO, IGLLOC, IGLSGLO, IGLSLOC, ILFB,&
 & ILPT, ILSB, IPT, IPTF, IST, JFLD, JGLLOC, JLP  

LOGICAL :: LLNOR, LLSUD
INTEGER(KIND=JPIM) :: IJPT
INTEGER(KIND=JPIM) :: IPTA(YDSL%NASLB1),IJPTA(YDSL%NASLB1)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SLEXTPOL',0,ZHOOK_HANDLE)
ASSOCIATE(NSTENCILWIDE=>YDDIM%NSTENCILWIDE)
!     ------------------------------------------------------------------

! KTYP=1: equivalent of the old SLEXTPOL
!  (assumes KFIXFLD(1)=1 ,KFIXFLD(2)=KFLDSLB, no mask required).
! KTYP=2: equivalent of the old SLEXTPOL2.
! KTYP=3: equivalent of the old SLEXTPOL1 and SLEXTPOL1A
!  (on-demand communications); requires the calculation of KMASK2

IF (PRESENT(KMASK2)) KMASK2(1:YDSL%NASLB1)=0

!     ** Search through all latitudes within halo region.

ILSB=MAX(YDSL%NDGSAG,YDSL%MYFRSTACTLAT-YDSL%NSLWIDEN)-YDSL%NFRSTLOFF
ILFB=MIN(YDSL%NDGENG,YDSL%MYLSTACTLAT +YDSL%NSLWIDES)-YDSL%NFRSTLOFF

!     ** Fix extra data required that has not been communicated:
!        1) pole mimics from core data.
!        2) end of latitude star copies.

IF (KTYP == 1) CALL GSTATS(1133,0)
IF (KTYP == 2) CALL GSTATS(1131,0)
IF (KTYP == 3) CALL GSTATS(1834,0)

DO JGLLOC=ILSB,ILFB

  IGLGLO=JGLLOC+YDSL%NFRSTLOFF
  LLNOR=(IGLGLO >= YDSL%NDGSAG).AND.(IGLGLO <= 0)
  LLSUD=(IGLGLO >= YDSL%NDGLG+1).AND.(IGLGLO <= YDSL%NDGENG)

  IF (LLNOR.OR.LLSUD) THEN

    IF (LLNOR) IGLSGLO=1-IGLGLO
    IF (LLSUD) IGLSGLO=2*YDSL%NDGLG+1-IGLGLO

    ! * Local view also required.

    IGLLOC=IGLGLO-YDSL%NFRSTLOFF
    IGLSLOC=IGLSGLO-YDSL%NFRSTLOFF

    ! * Determine which mimic points come from this processors core region.

    DO JLP=YDSL%NSLSTA(IGLLOC),YDSL%NSLONL(IGLLOC)+YDSL%NSLSTA(IGLLOC)-1

      ILPT=MOD(JLP-1+YDSL%NLOENG(IGLGLO),YDSL%NLOENG(IGLGLO))+1
      IPT =YDSL%NSLOFF(IGLLOC)+&
       & MOD(ILPT-YDSL%NSLSTA(IGLLOC)+YDSL%NLOENG(IGLGLO),YDSL%NLOENG(IGLGLO))  
      IPTF=YDSL%NSLOFF(IGLSLOC)+&
       & MOD(ILPT-YDSL%NSLSTA(IGLSLOC)+YDSL%NLOENG(IGLSGLO),YDSL%NLOENG(IGLSGLO))  
      IF(IPT >= YDSL%NASLB1 .OR. IPT < 1 .OR. IPTF >= YDSL%NASLB1 .OR. IPTF < 1 ) THEN
        WRITE(NULOUT,'(A,14I8)') ' ***** SLEXTPOL ERROR ***** ',&
         & IPT,IPTF,ILPT,JLP,IGLLOC,IGLGLO,IGLSLOC,&
         & YDSL%NASLB1,YDSL%NSLOFF(IGLLOC),YDSL%NSLOFF(IGLSLOC),&
         & YDSL%NSLSTA(IGLLOC),YDSL%NSLSTA(IGLSLOC),&
         & YDSL%NLOENG(IGLGLO),YDSL%NLOENG(IGLSGLO)  
        CALL ABOR1(' SLEXTPOL ERROR ')
      ENDIF
     
      ! update only if the polar mimic points are required on this processor

      IF (KTYP == 1) THEN
        IF ( (YDSL%NSLOFF(IGLSLOC) <= IPTF).AND.&
         & (IPTF <= YDSL%NSLOFF(IGLSLOC)+YDSL%NSLONL(IGLSLOC)-1) ) THEN  
          DO JFLD=1,KFLDSLB
            PB1A(IPT+1,JFLD,1)=PB1A(IPTF+1,JFLD,1)
          ENDDO
        ENDIF
        ! The mirrored latitudes close to the poles require correction for
        !  symmetric and antisymmetric properties (PARSLB).
        DO JFLD=1,KFLDSLB
          PB1A(IPT+1,JFLD,1)=PARSLB(JFLD)*PB1A(IPT+1,JFLD,1)
        ENDDO
      ELSEIF (KTYP == 2 .OR. KTYP == 3) THEN
        IJPT=IPT
        IF ( (YDSL%NSLOFF(IGLSLOC) <= IPTF).AND.&
         & (IPTF <= YDSL%NSLOFF(IGLSLOC)+YDSL%NSLONL(IGLSLOC)-1) ) IJPT=IPTF
        IPTA(JLP-YDSL%NSLSTA(IGLLOC)+1)=IPT
        IJPTA(JLP-YDSL%NSLSTA(IGLLOC)+1)=IJPT
      ENDIF

    ENDDO

    ! The mirrored latitudes close to the poles require correction for 
    !  symmetric and antisymmetric properties (PARSLB).

    IF (KTYP == 2) THEN
!$OMP PARALLEL PRIVATE(JFLD,JLP)
!$OMP DO SCHEDULE(STATIC)
      DO JFLD=1,KFIXFLD(1)-1
        DO JLP=1,YDSL%NSLONL(IGLLOC)
          PB1A(IPTA(JLP)+1,JFLD,1)=PARSLB(JFLD)*PB1A(IJPTA(JLP)+1,JFLD,1)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
      DO JFLD=KFIXFLD(2)+1,KFLDSLB
        DO JLP=1,YDSL%NSLONL(IGLLOC)
          PB1A(IPTA(JLP)+1,JFLD,1)=PARSLB(JFLD)*PB1A(IJPTA(JLP)+1,JFLD,1)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP END PARALLEL
    ELSEIF (KTYP == 3) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JLP)
      DO JFLD=KFIXFLD(1),KFIXFLD(2)
        DO JLP=1,YDSL%NSLONL(IGLLOC)
          PB1A(IPTA(JLP)+1,JFLD,1)=PARSLB(JFLD)*PB1A(IJPTA(JLP)+1,JFLD,1)
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      IF (PRESENT(KMASK2)) THEN
        DO JLP=1,YDSL%NSLONL(IGLLOC)
          KMASK2(IJPTA(JLP)+1)=1
        ENDDO
      ENDIF
    ENDIF

  ENDIF

ENDDO
!     ** If complete latitude present we correct for additional points
!        required at start/finish for interpolation star.

IF (KTYP == 1) THEN

!$OMP PARALLEL PRIVATE(JFLD,JGLLOC,IST,IEND)
!$OMP DO SCHEDULE(STATIC)
  DO JGLLOC=ILSB,ILFB
    IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
      IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
      IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
      DO JFLD=1,KFLDSLB
        IF(     NSTENCILWIDE==2 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+2,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+3,JFLD,1)
        ELSEIF( NSTENCILWIDE==3 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND-1,JFLD,1)
          PB1A(IST +2,JFLD,1)=PB1A(IEND  ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+3 ,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+4 ,JFLD,1)
          PB1A(IEND+3,JFLD,1)=PB1A(IST+5 ,JFLD,1)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ELSEIF (KTYP == 2) THEN

!$OMP PARALLEL PRIVATE(JFLD,JGLLOC,IST,IEND)
!$OMP DO SCHEDULE(STATIC)
  DO JFLD=1,KFIXFLD(1)-1
    DO JGLLOC=ILSB,ILFB
      IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
        IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
        IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
        IF(     NSTENCILWIDE==2 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+2,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+3,JFLD,1)
        ELSEIF( NSTENCILWIDE==3 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND-1,JFLD,1)
          PB1A(IST +2,JFLD,1)=PB1A(IEND  ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+3 ,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+4 ,JFLD,1)
          PB1A(IEND+3,JFLD,1)=PB1A(IST+5 ,JFLD,1)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO JFLD=KFIXFLD(2)+1,KFLDSLB
    DO JGLLOC=ILSB,ILFB
      IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
        IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
        IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
        IF(     NSTENCILWIDE==2 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+2,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+3,JFLD,1)
        ELSEIF( NSTENCILWIDE==3 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND-1,JFLD,1)
          PB1A(IST +2,JFLD,1)=PB1A(IEND  ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+3 ,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+4 ,JFLD,1)
          PB1A(IEND+3,JFLD,1)=PB1A(IST+5 ,JFLD,1)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

ELSEIF (KTYP == 3) THEN

!$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(JFLD,JGLLOC,IST,IEND)
  DO JFLD=KFIXFLD(1),KFIXFLD(2)
    DO JGLLOC=ILSB,ILFB
      IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
        IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
        IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
        IF(     NSTENCILWIDE==2 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+2,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+3,JFLD,1)
        ELSEIF( NSTENCILWIDE==3 )THEN
          PB1A(IST +1,JFLD,1)=PB1A(IEND-1,JFLD,1)
          PB1A(IST +2,JFLD,1)=PB1A(IEND  ,JFLD,1)
          PB1A(IEND+1,JFLD,1)=PB1A(IST+3 ,JFLD,1)
          PB1A(IEND+2,JFLD,1)=PB1A(IST+4 ,JFLD,1)
          PB1A(IEND+3,JFLD,1)=PB1A(IST+5 ,JFLD,1)
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

  DO JGLLOC=ILSB,ILFB
    IF (YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF) == YDSL%NSLONL(JGLLOC)) THEN
      IST =YDSL%NSLOFF(JGLLOC)-(NSTENCILWIDE-1)
      IEND=YDSL%NSLOFF(JGLLOC)+YDSL%NLOENG(JGLLOC+YDSL%NFRSTLOFF)
      IF (PRESENT(KMASK2)) THEN
        IF(     NSTENCILWIDE==2 )THEN
          KMASK2(IEND )=1
          KMASK2(IST+2)=1
          KMASK2(IST+3)=1
        ELSEIF( NSTENCILWIDE==3 )THEN
          KMASK2(IEND-1)=1
          KMASK2(IEND  )=1
          KMASK2(IST +3)=1
          KMASK2(IST +4)=1
          KMASK2(IST +5)=1
        ENDIF
      ENDIF
    ENDIF
  ENDDO

ENDIF


IF (KTYP == 1) CALL GSTATS(1133,1)
IF (KTYP == 2) CALL GSTATS(1131,1)
IF (KTYP == 3) CALL GSTATS(1834,1)

!     ------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SLEXTPOL',1,ZHOOK_HANDLE)
END SUBROUTINE SLEXTPOL
