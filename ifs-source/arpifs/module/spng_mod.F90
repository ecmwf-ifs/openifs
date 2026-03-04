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

!OCL  NOEVAL
MODULE SPNG_MOD

! Purpose :
! -------
!    Module for sponge applied at the top of the models.
!    This new sponge scheme can be used in 3D models.
!    - sponge is done in spectral space for GMV and spectral GFL.
!    - sponge is done in grid-point space for grid-point GFL.
!    Replaces the old sponge under NSPONGE=2

! Interface :
! ---------
!    Empty.

! External :  
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    K. YESSAD (CNRM/GMAP), after old sponge
!    Original : October 2011

! Modifications :
! -------------
!    T. Wilhelmsson (Sept 2013) Geometry and setup refactoring.
!    K. Yessad (July 2014): Move some variables.
!-----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE YOMLUN   , ONLY : NULNAM, NULOUT

IMPLICIT NONE
SAVE

!=============================================================================

!      1.    TYPE DEFINITIONS
!            ----------------

! none for the time being.

!=============================================================================

!      2.    DECLARATIONS
!            ------------

! Namelist variables:
! * LNSPONGE : to switch on new sponge:
! * RSPONBT  : altitude of the bottom of the sponge layer.
! * RSPONTP  : altitude of the top of the sponge layer.
! * RSPONTAU : time scale of absorption.
! * RSPONGN  : minimum value for RSPONGF.
! Other variables:
! * RSPONGF  : sponge function at full levels.
! Remarks:
! * for RSPONBT,RSPONTP,RSPONTAU,RSPONGN,RSPONGF: index 1 for GMV, 2 for sp. GFL, 3 for gp. GFL.

TYPE :: TSPNG
LOGICAL :: LNSPONGE
REAL(KIND=JPRB), ALLOCATABLE :: RSPONGF(:,:)

CONTAINS
  
  PROCEDURE, PASS :: PRINT => PRINT_CONFIGURATION 
    
END TYPE TSPNG

!!TYPE(TSPNG), POINTER :: YRSPNG => NULL()

#include "posnam.intfb.h"

!=============================================================================

CONTAINS

!      3.    ENCAPSULATED ROUTINES (SETUP, ALLOC, DEALLOC)
!            ---------------------------------------------

SUBROUTINE SUSPNG(YDSPNG,YDRIP,YDDYNA,NFLEVG,PSTZ)

!--------------------------------------------------------------------------
! Sets-up new sponge variables
!--------------------------------------------------------------------------

USE YOMRIP , ONLY : TRIP
USE YOMDYNA, ONLY : TDYNA
TYPE(TSPNG)       , TARGET, INTENT(INOUT) :: YDSPNG
TYPE(TRIP)        , INTENT(INOUT) :: YDRIP
TYPE(TDYNA)       , INTENT(INOUT) :: YDDYNA
INTEGER(KIND=JPIM), INTENT(IN)    :: NFLEVG
REAL(KIND=JPRB)   , INTENT(IN)    :: PSTZ(:)
CHARACTER(LEN = 80) :: CLFMT
INTEGER(KIND=JPIM) :: JLEV,JJ
REAL(KIND=JPRB) :: ZARG, ZTSTEP, ZTSTEP_REF, ZEXPO, ZRSPONGN
REAL(KIND=JPRB) :: ZALTF(NFLEVG)

REAL(KIND=JPRB) :: RSPONBT(3)
REAL(KIND=JPRB) :: RSPONTP(3)
REAL(KIND=JPRB) :: RSPONTAU(3)
REAL(KIND=JPRB) :: RSPONGN(3)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

LOGICAL, POINTER :: LNSPONGE
NAMELIST /NAMSPNG/ LNSPONGE,RSPONBT,RSPONTP,RSPONTAU,RSPONGN

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPNG_MOD:SUSPNG',0,ZHOOK_HANDLE)
ASSOCIATE(TSTEP=>YDRIP%TSTEP)
! Associate pointers for variables in namelist
LNSPONGE => YDSPNG%LNSPONGE
!--------------------------------------------------------------------------

!      3.1   Set default values:

IF ( YDDYNA%LTWOTL ) THEN
  ZTSTEP = TSTEP
ELSE
  ZTSTEP = 2.0_JPRB*TSTEP
ENDIF

! ZTSTEP_REF and ZEXPO can go later in namelist;
!  they allow to control sponge for "evanescent timesteps".
! ZTSTEP << ZTSTEP_REF => no sponge (ZRSPONGN=1).
! ZTSTEP >> ZTSTEP_REF => ZRSPONGN=RSPONGN
ZTSTEP_REF=1.0_JPRB
ZEXPO=2.0_JPRB

LNSPONGE=.FALSE.

RSPONBT(1) = 10000.0_JPRB
RSPONBT(2) = 10000.0_JPRB
RSPONBT(3) = 10000.0_JPRB
RSPONTP(1) = 20000.0_JPRB
RSPONTP(2) = 20000.0_JPRB
RSPONTP(3) = 20000.0_JPRB
RSPONTAU(1)= 5.0_JPRB*ZTSTEP
RSPONTAU(2)= 5.0_JPRB*ZTSTEP
RSPONTAU(3)= 5.0_JPRB*ZTSTEP
RSPONGN(1) = 0.0_JPRB
RSPONGN(2) = 0.0_JPRB
RSPONGN(3) = 0.0_JPRB

!      3.2   Read namelist:

CALL POSNAM(NULNAM,'NAMSPNG')
READ(NULNAM,NAMSPNG)

!      3.3   Write namelist variables:

WRITE(NULOUT, '(A)') '----- SPNG_MOD PRINTINGS:'
WRITE(NULOUT, '(A, L2)'   ) 'LNSPONGE  = ', LNSPONGE

IF (LNSPONGE) THEN
  WRITE(NULOUT, '(A)') 'NEW SPONGE IS ON'
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTAU for GMV     = ', RSPONTAU(1)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTAU for sp. GFL = ', RSPONTAU(2)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTAU for gp. GFL = ', RSPONTAU(3)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTP  for GMV     = ', RSPONTP(1)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTP  for sp. GFL = ', RSPONTP(2)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONTP  for gp. GFL = ', RSPONTP(3)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONBT  for GMV     = ', RSPONBT(1)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONBT  for sp. GFL = ', RSPONBT(2)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONBT  for gp. GFL = ', RSPONBT(3)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONGN  for GMV     = ', RSPONGN(1)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONGN  for sp. GFL = ', RSPONGN(2)
  WRITE(NULOUT, '(A, F10.2)') 'RSPONGN  for gp. GFL = ', RSPONGN(3)
  WRITE(NULOUT, '(A)') ' '
ENDIF

!      3.4   Allocate and compute RSPONGF:

IF (LNSPONGE) THEN

  ALLOCATE(YDSPNG%RSPONGF(NFLEVG,3))
  CLFMT = '(1X,''ARRAY '',A10,'' ALLOCATED '',8I8)'
  WRITE(NULOUT, CLFMT) 'RSPONGF', SIZE(YDSPNG%RSPONGF), SHAPE(YDSPNG%RSPONGF)

  ZALTF(1:NFLEVG)=PSTZ(1:NFLEVG)
  DO JJ=1,3
    ZRSPONGN=1.0_JPRB-TANH((ZTSTEP/ZTSTEP_REF)**ZEXPO)*(1.0_JPRB-RSPONGN(JJ))
    IF (ABS(RSPONTP(JJ)-RSPONBT(JJ)) > 1.E-10) THEN
      ! RSPONTP>RSPONBT (ZALTF may be between RSPONBT and RSPONTP).
      DO JLEV=1,NFLEVG
        IF ( ZALTF(JLEV) >= RSPONTP(JJ) ) THEN
          YDSPNG%RSPONGF(JLEV,JJ)=ZRSPONGN
        ELSEIF ( ZALTF(JLEV) <= RSPONBT(JJ) ) THEN
          YDSPNG%RSPONGF(JLEV,JJ)=1.0_JPRB
        ELSE
          ZARG=RPI*(RSPONTP(JJ)-ZALTF(JLEV))/(2.0_JPRB*(RSPONTP(JJ)-RSPONBT(JJ)))
          YDSPNG%RSPONGF(JLEV,JJ)=ZRSPONGN&
           & +(1.0_JPRB-ZRSPONGN)*RSPONTAU(JJ)*TAN(ZARG)/(ZTSTEP+RSPONTAU(JJ)*TAN(ZARG))
        ENDIF
      ENDDO
    ELSE
      ! RSPONTP=RSPONBT.
      DO JLEV=1,NFLEVG
        IF ( ZALTF(JLEV) >= RSPONTP(JJ) ) THEN
          YDSPNG%RSPONGF(JLEV,JJ)=ZRSPONGN
        ELSE
          YDSPNG%RSPONGF(JLEV,JJ)=1.0_JPRB
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  WRITE(NULOUT, '(A)') ' '
  WRITE(NULOUT, '(A)') 'ABSORBANT COEFFICIENT ON FULL LEVELS:'
  WRITE(NULOUT, '(A4,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9)') ' LEV',' ALTI (m)',&
   & '1-RSPONGF',' RSPONGF ','1-RSPONGF',' RSPONGF ','1-RSPONGF',' RSPONGF '
  WRITE(NULOUT, '(A4,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9,1X,A9)') '    ','         ',&
   & ' for GMV ',' for GMV ','for spGFL','for spGFL','for gpGFL','for gpGFL'
  DO JLEV=1,NFLEVG
    WRITE(NULOUT, '(I4,1X,F9.2,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6,1X,F9.6)')&
     & JLEV,ZALTF(JLEV),&
     & 1.0_JPRB-YDSPNG%RSPONGF(JLEV,1),YDSPNG%RSPONGF(JLEV,1),&
     & 1.0_JPRB-YDSPNG%RSPONGF(JLEV,2),YDSPNG%RSPONGF(JLEV,2),&
     & 1.0_JPRB-YDSPNG%RSPONGF(JLEV,3),YDSPNG%RSPONGF(JLEV,3)
  ENDDO
  WRITE(NULOUT, *)

ENDIF

!--------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SPNG_MOD:SUSPNG',1,ZHOOK_HANDLE)
END SUBROUTINE SUSPNG

!--------------------------------------------------------------------------
SUBROUTINE PRINT_CONFIGURATION(SELF, KDEPTH, KOUTNO)
IMPLICIT NONE
CLASS(TSPNG), INTENT(IN) :: SELF
INTEGER     , INTENT(IN) :: KDEPTH
INTEGER     , INTENT(IN) :: KOUTNO

INTEGER :: IDEPTHLOC 

IDEPTHLOC = KDEPTH + 2
WRITE(KOUTNO,*) REPEAT(' ',KDEPTH   ) // 'model%yrml_dyn%yrspng : '
WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'LNSPONGE = ', SELF%LNSPONGE
IF (ALLOCATED(SELF%RSPONGF)) THEN
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSPONGF allocated of shape ',SHAPE(SELF%RSPONGF)
ELSE
  WRITE(KOUTNO,*) REPEAT(' ',IDEPTHLOC) // 'RSPONGF not allocated'
ENDIF

END SUBROUTINE PRINT_CONFIGURATION
!--------------------------------------------------------------------------

END MODULE SPNG_MOD
