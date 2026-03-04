! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE N2O_O1D_INTERP(YGFL,PJ_VAL,PRSF1,PO1D_MMR)


!**   DESCRIPTION 
!     ----------
!
!   Part of TM5 routines for IFS chemistry: 
!--------------------------------------------------------------------------
!
!*** Surface boundary conditions for CH4 as a function of latitude
!*** Note: currently no seasonal cycle nor annual changes are implemented
!
!--------------------------------------------------------------------------
!
!
!**   INTERFACE.
!     ----------
!          *N2O_O1D_INTERP  IS CALLED FROM *CHEM_N2O*.

! INPUTS:
! -------
! PJ_VAL              : N2O and O3 Photolysis rate           
! PRSF1               : Midlevel pressure
!
!
! OUTPUTS:
! -------
! PTENCH4 (KLON)   : O1D_MMR based on full-BASCOE simulation
!
! LOCAL:
! -------
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-05-01



USE PARKIND1 , ONLY : JPIM, JPRB
USE YOM_YGFL , ONLY : TYPE_GFLD
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------

TYPE(TYPE_GFLD),INTENT(INOUT):: YGFL
REAL(KIND=JPRB),INTENT(IN)  :: PJ_VAL(2)       ! N2O and O3_O1D photolysis 
REAL(KIND=JPRB),INTENT(IN)  :: PRSF1           ! Pressure
REAL(KIND=JPRB),INTENT(OUT) :: PO1D_MMR       ! Conc

! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE
INTEGER(KIND=JPIM), PARAMETER :: J_NJN2O = 10
INTEGER(KIND=JPIM), PARAMETER :: J_NJO3  = 10
INTEGER(KIND=JPIM), PARAMETER :: J_NPRES = 13
REAL(KIND=JPRB),PARAMETER, DIMENSION(J_NPRES,J_NJN2O)  :: ZO1D_PARAM_2D=RESHAPE((/&
& 2.962E-17, 1.863E-17, 1.576E-17, 1.338E-17, 1.264E-17, 1.230E-17, 1.14804E-17, 9.819E-18, 5.953E-18, 3.421E-18 , 2.28E-19, &
& 1.62E-20, 4.75E-21 ,&
& 2.405E-15, 1.501E-15, 1.105E-15, 8.104E-16, 6.325E-16, 4.724E-16, 2.63643E-16, 1.559E-16, 7.529E-17, 2.304E-17, 0.0, 0.0, 0.0 ,&
& 4.121E-15, 3.004E-15, 2.205E-15, 1.574E-15, 1.204E-15, 8.253E-16, 5.19943E-16, 2.977E-16, 1.395E-16, 0.0 , 0.0, 0.0, 0.0 ,&
& 5.957E-15, 4.711E-15, 3.496E-15, 2.529E-15, 1.848E-15, 1.337E-15, 8.27981E-16, 4.916E-16, 0.00000, 0.000 , 0.0, 0.0, 0.0 ,&
& 7.893E-15, 6.596E-15, 5.072E-15, 3.721E-15, 2.679E-15, 1.989E-15, 1.26876E-15, 7.304E-16, 0.00000, 0.000 , 0.0, 0.0, 0.0 ,&
& 9.711E-15, 8.556E-15, 6.862E-15, 5.035E-15, 3.792E-15, 2.874E-15, 1.89136E-15, 0.00000  , 0.00000, 0.000 , 0.0, 0.0, 0.0 ,&
& 1.126E-14, 1.008E-14, 8.384E-15, 6.545E-15, 5.122E-15, 3.999E-15, 2.43533E-15, 0.00000  , 0.00000, 0.000 , 0.0, 0.0, 0.0 ,&
& 1.252E-14, 1.128E-14, 9.639E-15, 7.886E-15, 6.467E-15, 5.058E-15,  0.00000, 0.00000 , 0.00000 , 0.00000  , 0.0, 0.0, 0.0 ,&
& 1.319E-14, 1.212E-14, 1.054E-14, 8.956E-15, 7.424E-15,   0.00000,  0.00000, 0.00000 , 0.00000 , 0.00000  , 0.0, 0.0, 0.0 ,&
& 1.293E-14, 1.276E-14, 1.120E-14, 9.402E-15, 0.000    ,   0.00000,  0.00000, 0.00000 , 0.00000 , 0.00000  , 0.0, 0.0, 0.0 /),&
& (/J_NPRES,J_NJN2O/))

REAL(KIND=JPRB),PARAMETER, DIMENSION(J_NPRES,J_NJO3)  :: ZO3_O1D_PARAM_2D=RESHAPE((/&
&  2.767E-17, 1.786E-17, 1.667E-17, 1.719E-17, 2.054E-17, 2.662E-17, 3.721E-17, 3.898E-17, 2.828E-17, 4.311E-18, 2.28176E-19, &
&  1.62E-20,  4.75100E-21 ,&
&  2.338E-15, 1.581E-15, 1.334E-15, 1.202E-15, 1.212E-15, 9.647E-16, 7.995E-16, 5.434E-16,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  4.323E-15, 3.280E-15, 2.837E-15, 2.605E-15, 2.330E-15, 2.122E-15, 1.721E-15, 0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  5.815E-15, 4.786E-15, 4.265E-15, 3.920E-15, 3.516E-15, 3.202E-15, 2.361E-15, 0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  7.061E-15, 6.205E-15, 5.633E-15, 5.106E-15, 4.649E-15, 4.172E-15, 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  8.195E-15, 7.534E-15, 6.966E-15, 6.254E-15, 5.693E-15, 5.012E-15, 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  9.192E-15, 8.765E-15, 8.094E-15, 7.332E-15, 6.664E-15, 5.360E-15, 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  1.024E-14, 9.941E-15, 9.210E-15, 8.334E-15, 7.386E-15, 0.000    , 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  1.150E-14, 1.127E-14, 1.041E-14, 9.144E-15, 0.000    , 0.000    , 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000 ,&
&  1.241E-14, 1.253E-14, 1.121E-14, 0.000    , 0.000    , 0.000    , 0.00000,   0.00000  ,0.0, 0.00,0.00, 0.00000,  0.00000/),&
&  (/J_NPRES,J_NJO3/))
  
! use either J_N2O or J_O3_O1D as proxy for O1D
INTEGER(KIND=JPIM), PARAMETER :: IMODE = 2

REAL(KIND=JPRB),PARAMETER,DIMENSION(J_NJN2O) :: ZJN2O_PARAM=(/&
& 0.00000, 6.50000E-08, 1.30000E-07, 1.95000E-07, 2.60000E-07, 3.25000E-07, 3.90000E-07, 4.55000E-07, 5.20000E-07, 5.850E-07 /)
REAL(KIND=JPRB),PARAMETER,DIMENSION(J_NJO3) :: ZJO3_PARAM=(/ & 
& 0.0 ,  0.0008 ,  0.0016  , 0.0024  , 0.0032,  0.004  , 0.0048 ,  0.0056 ,  0.0064 ,  0.0072 /)
      
REAL(KIND=JPRB),PARAMETER,DIMENSION(J_NPRES) :: ZPRES_PARAM=(/&
& 10.0, 30.0,  50.0,  80.0,  120.0,  180.0,  300.0,  500.0,  1000.0,  2000.0, 5000.0, 1E4, 1.2E5 /)



! * counters
INTEGER(KIND=JPIM) :: JLEV,JLEVM,JN2O,JN2OM,JO3,JO3M
REAL(KIND=JPRB)    :: ZPFRAC,ZJFRAC

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('N2O_O1D_INTERP',0,ZHOOK_HANDLE )

ASSOCIATE(YCHEM=>YGFL%YCHEM)
!-----------------------------------------------------------------------
!  Prepare pressure interpolations: find the index jlev to use in O1D table
!-----------------------------------------------------------------------
      JLEV = 1
      DO WHILE( ZPRES_PARAM(JLEV) < PRSF1 )
         JLEV = JLEV + 1
      ENDDO
!-----------------------------------------------------------------------
!  Find interpolation weights
!-----------------------------------------------------------------------
      JLEVM=JLEV-1
      IF (JLEVM == 0 ) THEN 
        ZPFRAC = 1.0
        JLEVM = 1
      ELSE
        ZPFRAC=(PRSF1-ZPRES_PARAM(JLEVM)) / (ZPRES_PARAM(JLEV)-ZPRES_PARAM(JLEVM))
      ENDIF

    IF (IMODE == 1) THEN
!-----------------------------------------------------------------------
!  Prepare J N2O interpolations: find the index jN2O to use in O1D table
!-----------------------------------------------------------------------
      JN2O = 1
      DO WHILE( ZJN2O_PARAM(JN2O) < PJ_VAL(1) .AND. JN2O < J_NJN2O )
         JN2O = JN2O + 1
      ENDDO

!-----------------------------------------------------------------------
!  Find interpolation weights
!-----------------------------------------------------------------------

      JN2OM=JN2O-1
      IF (JN2OM == 0 ) THEN 
        ZJFRAC = 1.0
        JN2OM = 1
      ELSE
        ZJFRAC=(PJ_VAL(1)-ZJN2O_PARAM(JN2OM)) / (ZJN2O_PARAM(JN2O)-ZJN2O_PARAM(JN2OM))
      ENDIF


!-----------------------------------------------------------------------
!  Bilinear horizontal interp.
!-----------------------------------------------------------------------
      PO1D_MMR =      (1.-ZPFRAC) * (1.-ZJFRAC) * ZO1D_PARAM_2D(JLEVM,JN2OM)&
     &              +     ZPFRAC  * (1.-ZJFRAC) * ZO1D_PARAM_2D(JLEV ,JN2OM)&
     &              + (1.-ZPFRAC) *     ZJFRAC  * ZO1D_PARAM_2D(JLEVM,JN2O )&
     &              +     ZPFRAC  *     ZJFRAC  * ZO1D_PARAM_2D(JLEV ,JN2O )

    ELSEIF (IMODE == 2) THEN
    
!-----------------------------------------------------------------------
!  Prepare J O3 interpolations: find the index jN2O to use in O1D table
!-----------------------------------------------------------------------
      JO3 = 1
      DO WHILE( ZJO3_PARAM(JO3) < PJ_VAL(2) .AND. JO3 < J_NJO3 )
         JO3 = JO3 + 1
      ENDDO

!-----------------------------------------------------------------------
!  Find interpolation weights
!-----------------------------------------------------------------------

      JO3M=JO3-1
      IF (JO3M == 0 ) THEN 
        ZJFRAC = 1.0
        JO3M = 1
      ELSE
        ZJFRAC=(PJ_VAL(2)-ZJO3_PARAM(JO3M)) / (ZJO3_PARAM(JO3)-ZJO3_PARAM(JO3M))
      ENDIF


!-----------------------------------------------------------------------
!  Bilinear horizontal interp.
!-----------------------------------------------------------------------
      PO1D_MMR =      (1.-ZPFRAC) * (1.-ZJFRAC) * ZO3_O1D_PARAM_2D(JLEVM,JO3M)&
     &              +     ZPFRAC  * (1.-ZJFRAC) * ZO3_O1D_PARAM_2D(JLEV ,JO3M)&
     &              + (1.-ZPFRAC) *     ZJFRAC  * ZO3_O1D_PARAM_2D(JLEVM,JO3 )&
     &              +     ZPFRAC  *     ZJFRAC  * ZO3_O1D_PARAM_2D(JLEV ,JO3 )

    ENDIF




END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('N2O_O1D_INTERP',1,ZHOOK_HANDLE )
END SUBROUTINE N2O_O1D_INTERP

