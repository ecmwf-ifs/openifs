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

SUBROUTINE PPFLEV(KPROMA,KSTART,KPROF,KFLEV,KLEVP,KPPM,&
 & PRPRES,PRPRESH,PRPRESF,&
 & KLEVB,LDBELO,LDBELS,LDBLOW,LDBLES)  

!**** *PPFLEV* - FIND PRESSURE LEVEL UNDER SPECIFIED PRESSURES

!     PURPOSE.
!     --------
!       FINDS THE POINTER TO THE INPUT GRID LEVEL BELOW (HIGHER NUMBER)
!       THE SPECIFIED PRESSURES WE WANT TO INTERPOLATE TO.
!       FLAGS POINTS OUT OF RANGE OF INPUT GRID LEVELS.

!**   INTERFACE.
!     ----------
!        *CALL* *PPFLEV(...)*

!        EXPLICIT ARGUMENTS
!        --------------------

!        KPROMA                    - HORIZONTAL DIMENSION.             (INPUT)
!        KSTART                    - START OF WORK.                    (INPUT)
!        KPROF                     - DEPTH OF WORK.                    (INPUT)
!        KFLEV                     - NUMBER OF INPUT PRESSURE LEVELS   (INPUT)
!        KLEVP                     - NUMBER OF OUTPUT PRESSURE LEVELS  (INPUT)
!        KPPM             - Number of interpolation methods in post-processing

!        PRPRES(KPROMA,KLEVP)       - POST-PROCESSING LEVEL PRESSURES.  (INPUT)
!        PRPRESH(KPROMA,0:KFLEV)    - INPUT HALF LEVEL PRESSURES        (INPUT)
!        PRPRESF(KPROMA,KFLEV)      - INPUT FULL LEVEL PRESSURES        (INPUT)

!        KLEVB(KPROMA,KLEVP,KPPM)  - INPUT LEVEL BELOW PRPRES           (OUTPUT)
!        LDBELO(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     LOWEST (FULL) MODEL LEVEL        (OUTPUT)
!        LDBELS(KPROMA,KLEVP)      - .TRUE. IF PRESSURE IS UNDER
!                                     MODEL SURFACE                    (OUTPUT)
!        LDBLOW(KLEVP)             - .TRUE. IF LDBELO(J) IS CONTAINING
!                                    AT LEAST ONE .TRUE.               (OUTPUT)
!        LDBLES(KLEVP)             - .TRUE. IF LDBELS(J) IS CONTAINING
!                                    AT LEAST ONE .TRUE.               (OUTPUT)

!        IMPLICIT ARGUMENTS :  NONE.
!        --------------------

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.  NONE
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      MATS HAMRUD AND PHILIPRE COURTIER  *ECMWF*
!      ORIGINAL : 89-01-25

!     MODIFICATIONS.
!     --------------
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R. El Khatib     25-Apr-2007 Optimization (reduce memory conflicts)
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KPROMA
INTEGER(KIND=JPIM),INTENT(IN)    :: KSTART 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPROF 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLEV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEVP 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPPM
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRES(KPROMA,KLEVP) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESH(KPROMA,0:KFLEV) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRPRESF(KPROMA,KFLEV) 
INTEGER(KIND=JPIM),INTENT(INOUT) :: KLEVB(KPROMA,KLEVP,KPPM) 
LOGICAL           ,INTENT(OUT)   :: LDBELO(KPROMA,KLEVP) 
LOGICAL           ,INTENT(OUT)   :: LDBELS(KPROMA,KLEVP) 
LOGICAL           ,INTENT(OUT)   :: LDBLOW(KLEVP) 
LOGICAL           ,INTENT(OUT)   :: LDBLES(KLEVP) 
INTEGER(KIND=JPIM) :: ILB,JLEV,JLEVP,JL,IDONE
REAL(KIND=JPRB) :: ZRPRESHMAX(KFLEV),ZRPRESMIN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

!*       1.    FIND INPUT GRID LEVELS BELOW OUTPUT GRID LEVELS.
!              ------------------------------------------------

IF (LHOOK) CALL DR_HOOK('PPFLEV',0,ZHOOK_HANDLE)

DO JLEV=1,KFLEV
  ZRPRESHMAX(JLEV)=MAXVAL(PRPRESH(KSTART:KPROF,JLEV))
ENDDO

!*       1.1   FIND HALF AND FULL LEVEL BELOW SPECIFIED PRESSURE.

DO JLEVP=1,KLEVP
  IDONE=0
  ZRPRESMIN=MINVAL(PRPRES(KSTART:KPROF,JLEVP))
  DO JL=KSTART,KPROF
    KLEVB(JL,JLEVP,1)=KFLEV+1  
  ENDDO
  DO JLEV=1,KFLEV
    IF(ZRPRESHMAX(JLEV) > ZRPRESMIN) THEN
      DO JL=KSTART,KPROF
        IF( (PRPRESH(JL,JLEV) > PRPRES(JL,JLEVP))&
           & .AND. KLEVB(JL,JLEVP,1) == KFLEV+1) THEN  
          KLEVB(JL,JLEVP,1)=JLEV
          IDONE=IDONE+1
        ENDIF
      ENDDO
    ENDIF
    IF(IDONE == KPROF-KSTART+1) EXIT
  ENDDO
ENDDO

DO JLEVP=1,KLEVP
  DO JL=KSTART,KPROF
    KLEVB(JL,JLEVP,1)=MIN(KLEVB(JL,JLEVP,1),KFLEV)
    ILB=KLEVB(JL,JLEVP,1)
    IF ( PRPRESF(JL,ILB)-PRPRES(JL,JLEVP)  >=  0.0_JPRB ) THEN
      KLEVB(JL,JLEVP,2)=ILB
    ELSE
      KLEVB(JL,JLEVP,2)=ILB+1
    ENDIF
    KLEVB(JL,JLEVP,2)=MIN(KLEVB(JL,JLEVP,2),KFLEV)
    KLEVB(JL,JLEVP,3)=KLEVB(JL,JLEVP,1)
    KLEVB(JL,JLEVP,4)=KLEVB(JL,JLEVP,2)
  ENDDO
ENDDO

!*       1.2   FLAG VALUES OUTSIDE RANGE OF MODEL LEVELS.

DO JLEVP=1,KLEVP
  DO JL=KSTART,KPROF
    LDBELO(JL,JLEVP)=PRPRES(JL,JLEVP) > PRPRESF(JL,KFLEV)
    LDBELS(JL,JLEVP)=PRPRES(JL,JLEVP) > PRPRESH(JL,KFLEV)
  ENDDO
  LDBLOW(JLEVP) = ANY(LDBELO(KSTART:KPROF,JLEVP))
  LDBLES(JLEVP) = ANY(LDBELS(KSTART:KPROF,JLEVP))
ENDDO

IF (LHOOK) CALL DR_HOOK('PPFLEV',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE PPFLEV
