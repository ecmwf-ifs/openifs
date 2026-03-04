! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_PHO_CHAPMAN( PZKM, PSZA, PCHAPMAN, PH )
!-----------------------------------------------------------------------
!**   DESCRIPTION
!     ----------
!
!
!  Calculate the Chapman function.  Chapman function is used when the solar
!   zenith angle exceeds 65 degrees, but is not greater than 95 degrees.
!   The value is calculated following Smith and Smith (JGR,77,3592,1972).
!
!  In the heterosphere (i.e. above ~80km), molecular diffusion becomes
!   non-negligible -> different molar weights imply different scale heights
!   -> results can differ per species. To take this into account, give the
!   optional argument "PH"
!
!  Important modification: to correct the discontinuity when going from
!   secant (sza<=65) to chapman, we calc the secant for sza<75, use
!   it as fct value if sza<=65, and do a linear interp from secant to
!   chapman when going fromsza=65 to sza=75. This reduces, but does not
!   fix, the small discontinuity intrinsic in the Smith/Smith-1972 param
!   at y=8 (i.e. sza=71.58 for zkm=85 & h=5, see D:\...\chapman.xls)
!                                 simonc@oma.be, BASCOE v3s04, Nov 2003
!
! *** THIS WAS ORIGINALLY A FUNCTION, BUT gen_interfaces IS UNABLE
! ***   TO GENERATE INTERFACES FOR FUNCTIONS (ONLY FOR SUBROUTINES)
!           FUNCTION BASCOE_PHO_CHAPMAN( PZKM, PSZA )
!
!-----------------------------------------------------------------------
    USE PARKIND1 , ONLY : JPIM     ,JPRB
    USE YOMCST,    ONLY : RPI
    USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
    USE YOMLUN   , ONLY : NULOUT

    IMPLICIT NONE


! * INPUTS:
    REAL(KIND=JPRB), INTENT(IN)   :: PZKM     ! altitude (km)
    REAL(KIND=JPRB), INTENT(IN)   :: PSZA     ! solar zenith angle (degrees)
    REAL(KIND=JPRB), INTENT(OUT)  :: PCHAPMAN ! Chapman function
    REAL(KIND=JPRB), INTENT(IN), OPTIONAL   :: PH     ! scaled height

! * LOCAL
    CHARACTER(LEN=*), PARAMETER   :: CL_MY_NAME   = 'BASCOE_PHO_CHAPMAN'
    REAL(KIND=JPRB), PARAMETER :: ZD2R =  3.14159265/180._JPRB       ! degree to radian
    REAL(KIND=JPRB), PARAMETER :: ZR0 = 6371.0d0              ! effective earth radius (km)
    REAL(KIND=JPRB), PARAMETER :: ZHMEAN = 7._JPRB            ! mean scale height (km)

    INTEGER(KIND=JPIM) :: ierr
    REAL(KIND=JPRB)    :: Zxp, Zy, Zh, Zcoskhi, Zseckhi, Zsinkhi, Zerfcexp, Zchap_fct, Zresult
    REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


    IF (LHOOK) CALL DR_HOOK('BASCOE_PHO_CHAPMAN',0,ZHOOK_HANDLE )

    ierr = 0_JPIM
    Zcoskhi = COS( PSZA * ZD2R )
    Zseckhi = 1._JPRB / Zcoskhi
    IF( PSZA <= 65._JPRB ) THEN
       PCHAPMAN = Zseckhi
      IF (LHOOK) CALL DR_HOOK('BASCOE_PHO_CHAPMAN',1,ZHOOK_HANDLE )
      RETURN
    ENDIF

    Zh = ZHMEAN
    IF( PRESENT( PH ) ) ZH = PH

    Zxp = ( ZR0 + PZKM ) / ZH
    IF( Zxp < 50._JPRB .or. Zxp > 1.E9_JPRB ) ierr = 1_JPIM

    Zy = SQRT( 0.5_JPRB*Zxp ) * ABS( Zcoskhi )

    IF( Zy < 8._JPRB) THEN
      Zerfcexp = ( 1.0606963_JPRB + 0.55643831_JPRB*ZY ) &
   &          /  ( 1.0619896_JPRB + Zy*(1.7245609_JPRB + ZY) )
    ELSEIF( ZY <= 100._JPRB ) THEN
      Zerfcexp = 0.56498823_JPRB / (0.06651874_JPRB + Zy)
    ELSE
      Zerfcexp = 0.56498823_JPRB / (0.06651874_JPRB + Zy)
      ierr = 2_JPIM
    ENDIF

    IF( Zcoskhi >= 0._JPRB ) THEN
      Zchap_fct = SQRT(RPI*Zxp*0.5_JPRB ) * Zerfcexp
      IF( PSZA < 75._JPRB ) THEN
        ! Linear interp from seckhi at sza=65 to chap_fct at sza=75 ; 0.1 = 1 / (75-65)
         Zresult = (75._JPRB-PSZA)*Zseckhi*0.1_JPRB + (PSZA-65._JPRB)*Zchap_fct*0.1_JPRB
       ELSE
         Zresult = Zchap_fct
      ENDIF
    ELSE
      Zsinkhi = SQRT( 1._JPRB - Zcoskhi*Zcoskhi )
      Zchap_fct = SQRT( 2._JPRB*RPI*Zxp ) &
   &    * ( SQRT(Zsinkhi) * EXP( Zxp*(1._JPRB - Zsinkhi)) - 0.5_JPRB*Zerfcexp )
      Zresult = Zchap_fct
    ENDIF
    IF( Zchap_fct < 1. ) ierr = 3

    PCHAPMAN = Zresult

    IF( ierr > 0 ) THEN
       WRITE(NULOUT,*) CL_MY_NAME//' error, ierr= ',ierr
       WRITE(NULOUT,'(3(a,es12.3))') '   zkm= ',PZKM,' ; h= ',ZH,' ; xp= ',Zxp
       WRITE(NULOUT,'(3(a,es12.3))') '   sza= ',PSZA,' ; y= ',zy
       WRITE(NULOUT,'(3(a,es12.3))') '  erfcexp= ',Zerfcexp,' ; seckhi= ' &
   &           ,zseckhi,' ; chap_fct= ',zchap_fct
       WRITE(NULOUT,'(3(a,es12.3))')'   '//CL_MY_NAME//'= ',Zresult
       CALL ABOR1(CL_MY_NAME//': fatal error')

    ENDIF

    IF (LHOOK) CALL DR_HOOK('BASCOE_PHO_CHAPMAN',1,ZHOOK_HANDLE )

END SUBROUTINE BASCOE_PHO_CHAPMAN

