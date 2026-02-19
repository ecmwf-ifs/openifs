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

SUBROUTINE UPDCELAUT(YDML_PHY_STOCH,KCELL,PWGHT,KIIP,KJJP)

!**** *CA(KCELL,PWGHT,KIIP,KJJP)*  - Cellular Automaton
 
!     Purpose.
!     --------
!           Cellular Automaton update for stochastic physics

!**   Interface.
!     ----------
!        *CALL* *CA(KCELL,PWGHT,KIIP,KJJP)

!        Implicit arguments :
!        --------------------
!        COMMON STOPH_MIX

!     Author. 
!     ------- 
!        Glenn Shutts 

!     Modifications.
!     --------------
!        Original : 2004
!        Judith Berner:  Cleanup


USE MODEL_PHYSICS_STOCHAST_MOD , ONLY : MODEL_PHYSICS_STOCHAST_TYPE
USE PARKIND1  , ONLY : JPIM, JPRB
USE STOPH_MIX , ONLY : FERTILE_NEIGHBOUR_COUNT,NEIGHBOUR_COUNT,UPDATE_CELLS,WRAP_CELLS,WEIGHTING_FIELD
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK 
IMPLICIT NONE

TYPE(MODEL_PHYSICS_STOCHAST_TYPE),INTENT(INOUT):: YDML_PHY_STOCH
INTEGER(KIND=JPIM),INTENT(IN)    :: KIIP
INTEGER(KIND=JPIM),INTENT(IN)    :: KJJP
INTEGER(KIND=JPIM),INTENT(INOUT) :: KCELL(KIIP*4,KJJP*4)
REAL(KIND=JPRB)   ,INTENT(INOUT) :: PWGHT(KIIP,KJJP)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM)               :: IC,JC
INTEGER(KIND=JPIM), PARAMETER    :: ILIVES=10

!     based on the `generations' rule class of CELLular automata

INTEGER(KIND=JPIM), PARAMETER    :: INEIGHBOURHOOD= 8
INTEGER(KIND=JPIM)               :: ICELLS(0:KIIP*4+1,0:KJJP*4+1)

!     ICELLS() holds the state of the cellular automaton
!     KCELL() is used to compute the new state and this array is
!     input to the image-making subroutine MAKE_IFF_IMAGE()

INTEGER(KIND=JPIM)               :: ISURVIVAL_RULES(INEIGHBOURHOOD)
INTEGER(KIND=JPIM)               :: IBIRTH_RULES(INEIGHBOURHOOD)
INTEGER(KIND=JPIM)               :: I,J,ICOUNT,IN,IFERT
LOGICAL                          :: LLKEEPALIVE,LLCREATELIFE
DATA ISURVIVAL_RULES/3,4,5,0,0,0,0,0/
DATA IBIRTH_RULES   /2,3,0,0,0,0,0,0/

IF (LHOOK) CALL DR_HOOK('UPDCELAUT',0,ZHOOK_HANDLE)

IC=KIIP*4
JC=KJJP*4
CALL UPDATE_CELLS(KCELL,ICELLS,KIIP,KJJP)      
CALL WRAP_CELLS(ICELLS,KIIP,KJJP)

DO I= 1,IC
  DO J= 1,JC
    LLKEEPALIVE=  .FALSE.
    LLCREATELIFE= .FALSE.
    IN= NEIGHBOUR_COUNT(I,J,ICELLS,KIIP,KJJP)
    IFERT= FERTILE_NEIGHBOUR_COUNT(I,J,ICELLS,KIIP,KJJP)

!       check  survival  condition

    IF (ICELLS(I,J)/=0) THEN    ! check for living cell
      ICOUNT= 1
      DO WHILE(ISURVIVAL_RULES(ICOUNT)/=0)
        IF (IFERT==ISURVIVAL_RULES(ICOUNT)) LLKEEPALIVE= .TRUE.


!  N.B. the survival condition was earlier defined in terms of the 
!       age-independent neighbour count N.
!      It is now set to the `fertile' neighbour count IFERT i.e. w.r.t 
!      ICELLs enjoying their first year of life

        ICOUNT=  ICOUNT +1 
      ENDDO
    ELSE                        ! dead KCELL

!      check for birth condition  (only count `fertile' neighbours)  


      ICOUNT= 1
      DO WHILE(IBIRTH_RULES(ICOUNT)/=0)
        IF (IFERT==IBIRTH_RULES(ICOUNT)) LLCREATELIFE= .TRUE.
        ICOUNT=  ICOUNT +1
      ENDDO
    ENDIF   
    IF (.NOT.LLKEEPALIVE.AND.ICELLS(I,J)/=0) THEN
      KCELL(I,J)=  ICELLS(I,J) - 1  !  cells failing the survival
                                        !  test lose a life
      ELSEIF (LLCREATELIFE) THEN
        KCELL(I,J)=  ILIVES
      ELSE
        KCELL(I,J)=  ICELLS(I,J)
      ENDIF
  ENDDO
ENDDO

CALL UPDATE_CELLS(KCELL,ICELLS,KIIP,KJJP)
CALL WEIGHTING_FIELD(KCELL,PWGHT,KIIP,KJJP) !  average blocks of cells to a coarser grid

IF (LHOOK) CALL DR_HOOK('UPDCELAUT',1,ZHOOK_HANDLE)
END SUBROUTINE UPDCELAUT
