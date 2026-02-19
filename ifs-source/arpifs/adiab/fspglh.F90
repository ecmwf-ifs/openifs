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

SUBROUTINE FSPGLH(KM,KSL,KDGL,KFIELDS,PR1MU2,PFIELD,&
 & KPTRU,KFLDUV,KFLDSC,&
 & KFLDPTRUV )  

!**** *FSPGLH* - Dynamics calculations in Fourier space

!     Purpose.
!     --------
!         Perform dynamics calculations in Fourier space, specifically in
!         Fourier space before transposition (where you have access to all
!         latitudes for a given zonal wave number M).

!**   Interface.
!     ----------
!        CALL FSPGLH(...)
!        Explicit arguments : 
!        -------------------- 
!        KM - zonal wavenumer
!        KMLOC - local zonal wavenumber                         !!  REMOVED, CY45
!!       KSL - first latitude for given KM
!        KDGL - number of Gaussian latitudes
!        KFIELDS - number of fields in PFIELD
!        PR1MU2 - gaussian grid geometry  1.-MU*MU  , cos(theta)**2
!        PFIELD - fourier data
!        KPTRVOR - pointer to first vort.field in PFIELD        !!  REMOVED, CY45
!        KPTRDIV - pointer to first div. field in PFIELD        !!  REMOVED, CY45
!        KPTRU - pointer to first u field in PFIELD
!        KPTRV - pointer to first v field in PFIELD             !!  REMOVED, CY45
!        KPTRSC  - pointer to first scalar field in PFIELD      !!  REMOVED, CY45
!        KPTRNSD  - pointer to first N.S der.  field in PFIELD  !!  REMOVED, CY45        
!        KFLDUV - number of u/v fields (also for vor. and div.)
!        KFLDSC - number of scalar fields
!        KFLDNSD - number of N.S der.  field                    !!  REMOVED, CY45
!        KFLDPTRUV - level list for U-V type fields
!        KFLDPTRSC - level list for scalar fields               !!  REMOVED, CY45
!
!     Method.
!     -------
!       This routine is passed as an external to the transform package
!       when calling INV_TRANS. Inside the transform package the routine is 
!       invoked from routine FSPGL_INT. Tha argument list can thus not be
!       changed without a corrosponding change to this routine.
  
!     Externals.   None.
!     ----------
!     Author.
!     -------
!        Mats Hamrud *ECMWF*
!        Original : 00-10-26

!     Modifications.
!     --------------
!        P. Bechtold 26-Oct-2008  non-oro GWD (LEGWWMS)
!        G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!        K. Yessad  (July 2014): Move some variables.
!        O. Marsden (Aug 2016) : Pass R1MU2 and MYLEVS explicitly to routine as arguments
!        O. Marsden (June 2017): Store MYLEVS, R1MU2, and TSTEP in FSPGL_MOD module before going into 
!                                spectral space, and use stored values when called from trans library
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE FSPGL_MOD, ONLY : NFSPGL_MYLEVS, RFSPGL_KRF, RFSPGL_TSTEP
!!USE YOMDYNA  , ONLY : LRFRIC



!     ------------------------------------------------------------------

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)    :: KM 
INTEGER(KIND=JPIM),INTENT(IN)    :: KSL 
INTEGER(KIND=JPIM),INTENT(IN)    :: KDGL 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PR1MU2(KDGL)
INTEGER(KIND=JPIM),INTENT(IN)    :: KFIELDS 
REAL(KIND=JPRB)   ,TARGET,INTENT(INOUT) :: PFIELD(2*KFIELDS,0:KDGL+1) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KPTRU 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDUV 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDSC 
INTEGER(KIND=JPIM),INTENT(IN)    :: KFLDPTRUV(KFLDUV) 

!     ------------------------------------------------------------------

REAL(KIND=JPRB) :: ZCOS2
INTEGER(KIND=JPIM) :: JGL,JLEV,IR,ILEV
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FSPGLH',0,ZHOOK_HANDLE)
!!ASSOCIATE(RKRF=>YDDYN%RKRF, &
!! & LEGWWMS=>YDEPHY%LEGWWMS, &
!! & TSTEP=>YDRIP%TSTEP)

  
IF (.NOT. ALLOCATED(RFSPGL_KRF) ) CALL ABOR1("Trying to run FSPGLH without first having saved the necessary values")

!     ------------------------------------------------------------------

!*       1.1.1  RAYLEIGH FRICTION ON U ( FOR M=0 ONLY)

!!IF (LEGWWMS) THEN
!!  IF(YDDYN%LRFRIC .AND. LSLAG) THEN
    IF(KM == 0) THEN
      DO JGL=KSL,KDGL-KSL+1
        ZCOS2 = PR1MU2(JGL)**6
        DO JLEV=1,KFLDUV
          IR = KPTRU+2*JLEV-2
          ILEV = NFSPGL_MYLEVS(KFLDPTRUV(JLEV))
  !       IF(PFIELD(IR,JGL) > 0.0_JPRB ) THEN
            PFIELD(IR,JGL) = PFIELD(IR,JGL)/(1.0_JPRB+RFSPGL_KRF(ILEV)*ZCOS2*RFSPGL_TSTEP)
  !       ENDIF
        ENDDO
      ENDDO
    ENDIF
!!  ENDIF
!!ENDIF

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('FSPGLH',1,ZHOOK_HANDLE)
END SUBROUTINE FSPGLH
