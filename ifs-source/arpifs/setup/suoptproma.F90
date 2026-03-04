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

SUBROUTINE SUOPTPROMA(YDGEOMETRY)

!**** *SUOPTPROMA * - Routine to initialize parallel environment

!     Purpose.
!     --------
!           Initialize NPROM* variables
!           NPROMA optimization

!**   Interface.
!     ----------
!        *CALL* *SUOPTPROMA *

!        Explicit arguments :
!        --------------------

!        Implicit arguments :
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.   NONE.
!     ----------

!     Reference.
!     ----------
!        METEO-FRANCE Research Department documentation of the ALADIN

!     Author.
!     -------
!      Bogdan Bochenek
!      Original : 2012-04-12

!     Modifications.
!     --------------
!      K. Yessad (July 2014): Move some variables.
!      K. Yessad (Feb 2018): remove deep-layer formulations.
!      M.Hamrud  (Jan 2022) remove NPROMM9 etc as they intruduce undesireble
!                           relience on LTWOTKL in geometry object
!     End Modifications
!     ------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE OML_MOD  , ONLY : OML_MAX_THREADS
USE YOMCT0   , ONLY : LELAM
USE YOMLUN   , ONLY : NULOUT
!!USE YOMPHY   , ONLY : YRPHY
!!USE YOMSIMPHL, ONLY : YRSIMPHL
USE YOMMP0   , ONLY : LOPT_SCALAR

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY), INTENT(INOUT) :: YDGEOMETRY
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: I, IPLAT1, IPLAT2, JJ,&
 & IGPBLKSMX, IMX, IPR, IPROMABEG, IPROMAEND, IPROMAREF, IPROMA

INTEGER(KIND=JPIM) :: ITHRMAX

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUOPTPROMA',0,ZHOOK_HANDLE)
ASSOCIATE(YDDIM=>YDGEOMETRY%YRDIM,YDDIMV=>YDGEOMETRY%YRDIMV,YDGEM=>YDGEOMETRY%YRGEM, YDMP=>YDGEOMETRY%YRMP)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, &
 & NPROMNH=>YDDIM%NPROMNH,  &
 & NPROMA=>YDDIM%NPROMA, &
 & LOPTPROMA=>YDDIM%LOPTPROMA, &
 & NPROMM=>YDDIM%NPROMM, &
 & NGPTOT=>YDGEM%NGPTOT, NGPTOTMX=>YDGEM%NGPTOTMX)

!     ------------------------------------------------------------------



!     ------------------------------------------------------------------

IF (LELAM) THEN

!  Adjust NPROMA to optimise for memory usage and vector length
!  Only do this if LOPTPROMA is true

  IF (LOPTPROMA) THEN
    IF (NGPTOTMX < NPROMA) THEN
      NPROMA=NGPTOTMX
    ELSE
      ! Platform dependent parameters (defining lower: 100-IPLAT1 % and
      !  upper: 100+IPLAT2 % bounds for NPROMA optimization)
      IF (LOPT_SCALAR) THEN
        IPLAT1=0
        IPLAT2=5
      ELSE
        IPLAT1=10
        IPLAT2=15
      ENDIF

#ifdef NECSX
      ! vector platforms require some minimal NPROMA to effectively feed
      !  registers (for SX9 this is 1800).
      IF (.NOT. LOPT_SCALAR) THEN
        IF (NPROMA < 1800) THEN
           NPROMA=1800
           IPLAT1=0
        ELSE
          IPLAT1=MIN((NPROMA-1800)/NPROMA,IPLAT1)
        ENDIF
      ENDIF
#endif

      ! now do some optimization within the given limits
      ITHRMAX=OML_MAX_THREADS()
      IF( ITHRMAX == 1 )THEN
        ! Allow NPROMA to drift upwards if this results
        ! in a reduction of the number of grid point blocks
        IPROMABEG=NPROMA+1
        IPROMAEND=IPROMABEG+NPROMA*IPLAT2/100
        IGPBLKSMX=(NGPTOTMX-1)/NPROMA+1
        NPROMA=(NGPTOTMX-1)/IGPBLKSMX+1
        DO JJ=IPROMABEG,IPROMAEND
          IMX=(NGPTOTMX-1)/JJ+1
          IPR=(NGPTOTMX-1)/IMX+1
          IF (IMX < IGPBLKSMX) THEN
            IGPBLKSMX=IMX
            NPROMA=IPR
          ENDIF
        ENDDO
      ELSE
        !   Optimise NPROMA to better balance the OMP threads
        IPROMABEG=MAX(1,NPROMA*(100-IPLAT1)/100)
        IPROMAEND=NPROMA*(100+IPLAT2)/100
        IPROMAREF=NPROMA
        NPROMA=(NGPTOT-1)/ITHRMAX+1
        I=1
        DO WHILE( NPROMA > IPROMAEND )
          I=I+1
          NPROMA=(NGPTOT-1)/(ITHRMAX*I)+1
        ENDDO
        !In case, there's no space for optimization within the given
        !  interval, keep the original value.
        IF (NPROMA < IPROMABEG) NPROMA=IPROMAREF
      ENDIF
      IPROMABEG=NPROMA

      ! Finally some platform dependent optimizations:
      IF (LOPT_SCALAR) THEN
        !  Force NPROMA to be odd
        IF( MOD(NPROMA,2) == 0 ) NPROMA=NPROMA+1
      ELSE
        ! Vector platforms require more constraints:
        !  Make sure NPROMA is odd, NPROMA+1 is not multiple of 4 and
        !  NPROMA= N * NBANKS + odd_offset
        DO WHILE ( (NPROMA >= IPROMABEG) .AND. (NPROMA <= IPROMAEND) .AND. &
         & ((MOD(NPROMA,2) == 0).OR.(MOD(NPROMA+1,4) == 0).OR. &
         &  (MOD(MOD(NPROMA,256),2) == 0) ) )
           NPROMA=NPROMA+1
        ENDDO
      ENDIF

      ! At the end of optimization re-check the basic things...
      IF (NGPTOTMX < NPROMA) NPROMA=NGPTOTMX
    ENDIF
    WRITE(NULOUT,'("SUOPTPROMA: NPROMA OPTIMISED =",I6)') NPROMA
  ENDIF

ELSE

  IF (LOPTPROMA) THEN
  !   Allow nproma to drift upwards by up to 5 percent if this results
  !   in a reduction of the number of grid point blocks
    ITHRMAX=OML_MAX_THREADS()
    IF( ITHRMAX > 1 )THEN
      !  Optimise NPROMA for maximum number of threads
      IPROMA=NPROMA*1.05
      NPROMA=(NGPTOTMX-1)/ITHRMAX+1
      I=1
      DO WHILE( NPROMA > IPROMA )
        I=I+1
        NPROMA=(NGPTOTMX-1)/(ITHRMAX*I)+1
      ENDDO
      ! Check if we have enough blocks to switch to a more optimal nproma
      ! i.e. more than 10 blocks per thread (for scalar only)
      IGPBLKSMX=(NGPTOTMX-1)/NPROMA+1
      IF( LOPT_SCALAR )THEN
        IF( IGPBLKSMX/ITHRMAX > 10 )THEN
          IF( NFLEVG <= 62 )THEN
            NPROMA=56
          ELSEIF( NFLEVG <= 91 )THEN
            NPROMA=32
          ELSE
            NPROMA=24
          ENDIF
        ENDIF
      ENDIF
    ELSE
      IPROMABEG=NPROMA+1
      IPROMAEND=IPROMABEG+NPROMA*5/100
      IGPBLKSMX=(NGPTOTMX-1)/NPROMA+1
      NPROMA=(NGPTOTMX-1)/IGPBLKSMX+1
      DO JJ=IPROMABEG,IPROMAEND
        IMX=(NGPTOTMX-1)/JJ+1
        IPR=(NGPTOTMX-1)/IMX+1
        IF (IMX < IGPBLKSMX) THEN
          IGPBLKSMX=IMX
          NPROMA=IPR
        ENDIF
      ENDDO
      WRITE(NULOUT,'("SUOPTPROMA: NPROMA OPTIMISED =",I6)') NPROMA
    ENDIF
  ENDIF
  
ENDIF


!! FIXME: LMPHYS and LSIMPH belong to MODEL and cannot be used during GEOMETRY setup.
NPROMM=NPROMA
IF(YDGEOMETRY%LNONHYD_GEOM) THEN
  NPROMNH=NPROMA
ELSE
  NPROMNH=1
ENDIF

WRITE(NULOUT,'(''SUOPTPROMA:  COLLOCATION GRID DIMENSIONING '')')
WRITE(UNIT=NULOUT,FMT='('' NPROMM='',I6,'' NPROMNH='',I6)') NPROMM,NPROMNH

!     -----------------------------------------------------------

END ASSOCIATE
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUOPTPROMA',1,ZHOOK_HANDLE)
END SUBROUTINE SUOPTPROMA
