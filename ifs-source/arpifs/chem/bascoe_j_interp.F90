! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_J_INTERP(PSZA,PZKM,PO3_OVER, PAJVAL )


!**   DESCRIPTION 
!     ----------
!
!   Part of BASCOE / TM5 routines for IFS chemistry: 
!   The evaluation of photolysis rates
!
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_J_INTERP* IS CALLED FROM *CHEM_tm5 / CHEM_BASCOE*.

! INPUTS:
! -------
!
! PSZA                ! Solar zenith angle (degrees)
! PZKM                ! Geometric altitude (km)
! PO3VMR, PO3VMR_ABOVE! ozone volume mixing ratio
! PRES, PRES_ABOVE    ! pressure (Pa)
! PO3_OVER            ! O3 column above model point (DU)
!
! OUTPUTS:
! -------
! PAJVAL  (NPHOTO)     : photolysis rates
!
!     AUTHOR.
!     -------
!        Coded in C-IFS by VINCENT HUIJNEN    *KNMI*
!        Original code from BASCOE_CTM v4s09, simonc@oma.be, June 2008
!
!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2014-02-01
!
! LOCAL:
! -------
USE PARKIND1 , ONLY : JPIM     ,JPRB
USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RPI
USE YOMLUN   , ONLY : NULOUT
USE BASCOE_J_MODULE , ONLY :  NDISS
USE BASCOE_J_TABLES_MODULE , ONLY :  NZEN =>NSZA, NJO3, NJLEV, &
                            & ZJLEV_TABLES, SZA_TABLES, O3COL_TABLES, ALJ_TABLES

IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------
REAL(KIND=JPRB),INTENT(IN) :: PSZA                ! Solar zenith angle (degrees)
REAL(KIND=JPRB),INTENT(IN) :: PZKM                ! Geometric altitude (km)
REAL(KIND=JPRB),INTENT(IN) :: PO3_OVER            ! O3 column above model point (DU)
REAL(KIND=JPRB),DIMENSION(NDISS), INTENT(out) :: PAJVAL  ! J-values


! * LOCAL 
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

LOGICAL, PARAMETER :: lldebug = .false.
REAL(KIND=JPRB), PARAMETER :: Zvery_large = 9.E29
!REAL(KIND=JPRB), PARAMETER :: cst = (1.E-4/RG)*( RNAVO / (1.E-3* RMD) ) 1.E3 / 2.687E19  ! molec/cm2 -> DU
REAL(KIND=JPRB), PARAMETER :: ZD2R =  3.14159265/180._JPRB       ! degree to radian


INTEGER(KIND=JPIM) :: izen, JLEV, jn, idx1(1), io3lo, io3hi, idx(4)
INTEGER(KIND=JPIM) :: iok, ierr, itmp, JO3
REAL(KIND=JPRB) :: Zdlev, Zscsza, Zrd(2), Zalj_o3lo, Zalj_o3hi, Zdo3, Zc0, Zc1
REAL(KIND=JPRB), DIMENSION(NJO3) :: Zb1, Zb0, Zo3col_modlev, Zxtmp
REAL(KIND=JPRB), DIMENSION(nzen) :: Zsct_modlev


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "bascoe_interp2d_loc.intfb.h"
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('BASCOE_J_INTERP',0,ZHOOK_HANDLE )
PAJVAL(1:NDISS)=0._JPRB
!-----------------------------------------------------------------------
!  Calculate the secant/chapman of the sza_tables, at model altitudes
!-----------------------------------------------------------------------
      do izen = 1, nzen
         IF( sza_tables(izen) < 65._JPRB ) THEN
            Zsct_modlev(izen) = 1._JPRB/ COS( sza_tables(izen)*ZD2R )
         ELSE
            Zsct_modlev(izen) = PHO_CHAPMAN( PZKM, sza_tables(izen) )
         ENDIF
      ENDDO

!-----------------------------------------------------------------------
!  Secant/Chapman of PSZA at model gridpoint
!-----------------------------------------------------------------------
      IF( PSZA < 65._JPRB ) THEN
        Zscsza = 1.d0 / COS( PSZA*ZD2R )
      ELSE
        Zscsza = PHO_CHAPMAN( PZKM, PSZA )
      endif

!-----------------------------------------------------------------------
!  Find 2D interpolation weights: X-axis = scsza ; Y-axis = altitude
!-----------------------------------------------------------------------
      CALL BASCOE_INTERP2D_LOC( nzen, NJLEV, Zsct_modlev, ZJLEV_TABLES, .false.,  &
    &                   Zscsza, PZKM, idx, Zrd, iok )
      IF( iok /= 0 ) THEN
         WRITE(NULOUT,*)nzen, NJLEV, Zscsza, PZKM, iok
         WRITE(NULOUT,*)Zsct_modlev
         WRITE(NULOUT,*)'J_INTERP: INTERP2D_LOC failed'
         CALL ABOR1(" BASCOE_J_INTERP: INTERP2D_LOC failed ")
      ENDIF

!-----------------------------------------------------------------------
!  Prepare altitude interpolations: find the index JLEV to use in J tables
!-----------------------------------------------------------------------
      JLEV = 2_JPIM
      do while( ZJLEV_TABLES(JLEV) < PZKM )
         JLEV = JLEV + 1_JPIM
      ENDDO

!-----------------------------------------------------------------------
!  Interpolate o3col_tables as a function of model altitude
!-----------------------------------------------------------------------
      IF( ( ZJLEV_TABLES(JLEV) - PZKM ) >= 1.E-3_JPRB ) THEN
         Zdlev = ZJLEV_TABLES(JLEV) - ZJLEV_TABLES(JLEV-1)
         Zb1(1:NJO3) = (o3col_tables(JLEV,1:NJO3)-o3col_tables(JLEV-1,1:NJO3)) / Zdlev
         Zb0(1:NJO3) = o3col_tables(JLEV,1:NJO3) - zb1(1:NJO3)*ZJLEV_TABLES(JLEV)
         Zo3col_modlev(1:NJO3) = Zb1(1:NJO3)*PZKM + Zb0(1:NJO3)
       ELSE
         Zo3col_modlev(1:NJO3) = o3col_tables(JLEV,1:NJO3)
      ENDIF

!-----------------------------------------------------------------------
! Calculate overhead ozone columns *at* levels:  cst * SUM( vmr * delta_p )
! where delta_p is between levels and vmr are at mid-levels
!
! VH - This part is moved to chem_bascoe / chem_tm5, as it is prone to
!    - coding errors  (specifically use of INTENT(INOUT) )
! 
!-----------------------------------------------------------------------
!      IF( Po3_over == 0.0_JPRB ) THEN
!         Po3_over = cst * PO3VMR * ( PRES - 0._JPRB )
!       ELSE
!         Po3_over = Po3_over &
!     &           + cst * ( PRES - PRES_ABOVE ) * ( PO3VMR + PO3VMR_ABOVE )*0.5d0
!      ENDIF

!-----------------------------------------------------------------------
!  Find [io3_lo,io3_hi] indices in tables around model overhead O3 columns
!-----------------------------------------------------------------------
      IF( PO3_OVER <= Zo3col_modlev(1) ) THEN
         io3lo = 1_JPIM
         io3hi = 1_JPIM
       ELSEIF( PO3_OVER >= Zo3col_modlev(NJO3) ) THEN
         io3lo = NJO3
         io3hi = NJO3
       ELSE
         Zxtmp(1:NJO3) = Zo3col_modlev(1:NJO3)
         ! where( Zxtmp < PO3_OVER ) Zxtmp = Zvery_large
         DO JO3=1,NJO3
            IF (ZXTMP(JO3) < PO3_OVER) THEN
              ZXTMP(JO3)=ZVERY_LARGE
            ENDIF
         ENDDO
         idx1 = MINLOC( Zxtmp )
         io3hi = idx1(1)
         io3lo = idx1(1) - 1_JPIM
      ENDIF
      IF( io3lo < 1 .OR. io3hi > NJO3 ) THEN
        WRITE(NULOUT,*)'INTERP: io3 failed'
        CALL ABOR1(" BASCOE_J_INTERP: INTERP: io3 failed ")
      ENDIF

!-----------------------------------------------------------------------
!  2D Interp for tables with O3col corresponding to io3_lo and io3_hi
!-----------------------------------------------------------------------
      do jn = 1, ndiss
         Zalj_o3lo = INTERP2D_FAST( nzen, NJLEV, alj_tables(1:NZEN,1:njlev,io3lo,jn), idx, Zrd )
         IF( io3hi == io3lo ) THEN
            PAJVAL(jn) = EXP( Zalj_o3lo )
          ELSE
            Zalj_o3hi = INTERP2D_FAST( nzen, NJLEV, alj_tables(1:NZEN,1:njlev,io3hi,jn), idx, Zrd )

!-----------------------------------------------------------------------
!  Interpolate as a function of the overhead O3 columns
!-----------------------------------------------------------------------
            Zdo3 = Zo3col_modlev(io3hi) - Zo3col_modlev(io3lo)
            Zc1 = ( Zalj_o3hi - Zalj_o3lo ) / Zdo3
            Zc0 = Zalj_o3hi - Zc1*Zo3col_modlev(io3hi)
            PAJVAL(jn) = EXP( Zc1* PO3_OVER + Zc0 )
         ENDIF

      ENDDO

!-----------------------------------------------------------------------
!  Diagnostics and error reporting. Tests on alj0,alj1 to catch errors
!    when alj_global has still its default value as set in J_TABLES_READ
!-----------------------------------------------------------------------
      IF( lldebug ) THEN
         ierr = 0
         IF( PZKM<ZJLEV_TABLES(1) .OR. PZKM>ZJLEV_TABLES(NJLEV) ) ierr=1
         IF( JLEV > NJLEV ) ierr = 2
         IF( Zsct_modlev(1) /= 1. ) ierr = 3
         do izen = 2, nzen
           IF( Zsct_modlev(izen) <= Zsct_modlev(izen-1) ) ierr = 4
         ENDDO
         IF( izen > nzen ) ierr = 1
         IF( Zsct_modlev(1) /= 1. ) ierr = ierr + 1
         do itmp = 2, nzen
            IF( Zsct_modlev(itmp)<=Zsct_modlev(itmp-1) ) ierr = ierr+10
         ENDDO
         IF( PAJVAL(jn) < 0. ) ierr = ierr + 10000
         IF( PAJVAL(jn) >= 1. ) ierr = ierr + 10000
      ENDIF
      IF(  lldebug .and. ierr>0 ) THEN
         WRITE(NULOUT,'(a,3i3,a,i3)')' J_INTERP for process ',jn
         IF(ierr>0)WRITE(NULOUT,*)'   ERROR !!!! ierr= ',ierr
         WRITE(NULOUT,'(2(a,f9.3),a,i3)') '   PSZA= ',PSZA &
     &     ,' ; PZKM= ',PZKM,' ; izen= ',izen
         WRITE(NULOUT,'(a,f6.2,a,i3,a,f9.3)') ' J_INTERP,   PZKM= ',PZKM &
     &    ,' ; JLEV= ',JLEV,' ; ZJLEV_TABLES(JLEV)= ',ZJLEV_TABLES(JLEV)
         WRITE(NULOUT,'(2(a,f9.3))') &
     &     ' ; ZJLEV_TABLES(1)= ',ZJLEV_TABLES(1) &
     &    ,' ; ZJLEV_TABLES(NJLEV)= ',ZJLEV_TABLES(NJLEV)
         do izen = 1, nzen
            WRITE(NULOUT,'(a,i3,a,f9.3,a,es12.5)')'   izen= ',izen &
     &       ,' ; sza_tables= ',sza_tables(izen) &
     &       ,' ; sct_modlev= ',Zsct_modlev(izen)
         ENDDO
         WRITE(NULOUT,'(3(a,es12.3))') '   sct_modlev(izen-1)= ' &
     &     ,Zsct_modlev(izen-1),' ; scsza= ',Zscsza &
     &     ,' ; sct_modlev(izen)= ',Zsct_modlev(izen)
         WRITE(NULOUT,'(3(a,es12.3))') '   ajval= ',PAJVAL(jn)
         IF( lldebug .and. ierr > 0 ) THEN
            do itmp = 1, nzen
               WRITE(NULOUT,'(a,i3,a,es12.5)')'   izen= ',itmp &
     &             ,' ; sc_tab= ',Zsct_modlev(itmp)
            ENDDO
            WRITE(NULOUT,*) 'J_INTERP fatal error'
            CALL ABOR1(" BASCOE_J_INTERP: fatal error ")
         ENDIF
      ENDIF


IF (LHOOK) CALL DR_HOOK('BASCOE_J_INTERP',1,ZHOOK_HANDLE )

CONTAINS

  FUNCTION PHO_CHAPMAN( PZKM, PSZA )

!-----------------------------------------------------------------------
!  Calculate the Chapman function.  Chapman function is used when the solar
!  zenith angle exceeds 65 degrees, but is not greater than 95 degrees.
!  The fct value is calculated following Smith and Smith (JGR,77,3592,1972).
!
!  Important modification: to correct the discontinuity when going from
!   secant (sza<=65) to chapman, we calc the secant for sza<75, use
!   it as fct value if sza<=65, and do a linear interp from secant to
!   chapman when going fromsza=65 to sza=75. This reduces, but does not
!   fix, the small discontinuity intrinsic in the Smith/Smith-1972 param
!   at y=8 (i.e. sza=71.58 for zkm=85 & h=5, see D:\...\chapman.xls)
!                                 simonc@oma.be, BASCOE v3s04, Nov 2003
!-----------------------------------------------------------------------
      USE PARKIND1 , ONLY : JPIM     ,JPRB
      USE YOMCST,    ONLY : RPI
      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOMLUN   , ONLY : NULOUT

      IMPLICIT NONE

      REAL(KIND=JPRB)            :: PHO_CHAPMAN
      REAL(KIND=JPRB), PARAMETER :: ZD2R =  3.14159265/180._JPRB       ! degree to radian
      REAL(KIND=JPRB), PARAMETER :: ZR0 = 6356.7                ! effective earth radius (km)

! * INPUTS:
      REAL(KIND=JPRB), INTENT(IN)   :: PZKM     ! altitude (km)
      REAL(KIND=JPRB), INTENT(IN)   :: PSZA     ! solar zenith angle (degrees)

! * LOCAL 
      REAL(KIND=JPRB), parameter :: ZH = 7._JPRB            ! mean scale height (km)
      INTEGER(KIND=JPIM) :: ierr
      REAL(KIND=JPRB)    :: Zxp, Zy, Zcoskhi, Zseckhi, Zsinkhi, Zerfcexp, Zchap_fct, Zresult
#ifdef WITH_COMPO_DR_HOOK
      REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE


      IF (LHOOK) CALL DR_HOOK('BASCOE_J_INTERP:PHO_CHAPMAN',0,ZHOOK_HANDLE )
#endif
      ierr = 0_JPIM 
      Zcoskhi = COS( PSZA * ZD2R )
      Zseckhi = 1._JPRB / Zcoskhi
      IF( PSZA <= 65._JPRB ) THEN
         PHO_CHAPMAN = Zseckhi
#ifdef WITH_COMPO_DR_HOOK
        IF (LHOOK) CALL DR_HOOK('BASCOE_J_INTERP:PHO_CHAPMAN',1,ZHOOK_HANDLE )
#endif
        RETURN
      ENDIF

      Zxp = ( ZR0 + PZKM ) / ZH
      IF( Zxp < 50._JPRB .OR. Zxp > 1.E9_JPRB ) ierr = 1_JPIM

      ZY = SQRT( 0.5_JPRB*Zxp ) * ABS( Zcoskhi )

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

      PHO_CHAPMAN = Zresult

      IF( ierr > 0 ) THEN
         WRITE(NULOUT,*) 'PHO_CHAPMAN error, ierr= ',ierr
         WRITE(NULOUT,'(3(a,es12.3))') '   zkm= ',PZKM,' ; h= ',ZH,' ; xp= ',Zxp
         WRITE(NULOUT,'(3(a,es12.3))') '   sza= ',PSZA,' ; y= ',zy
         WRITE(NULOUT,'(3(a,es12.3))') '  erfcexp= ',Zerfcexp,' ; seckhi= ' &
     &           ,zseckhi,' ; chap_fct= ',zchap_fct
         WRITE(NULOUT,'(3(a,es12.3))')'   PHO_CHAPMAN= ',Zresult
         CALL ABOR1(" PHO_CHAPMAN: fatal error ")
 
      ENDIF

#ifdef WITH_COMPO_DR_HOOK
      IF (LHOOK) CALL DR_HOOK('BASCOE_J_INTERP:PHO_CHAPMAN',1,ZHOOK_HANDLE )
#endif

  END FUNCTION PHO_CHAPMAN

  FUNCTION INTERP2D_FAST( KX, KY, P_IN, KDX, PRD )

      USE PARKIND1 , ONLY : JPIM,JPRB
      IMPLICIT NONE

!-----------------------------------------------------------------------
!   ... Dummy args.
!  If lxlon==.true., x is longitudes from 0 to (360-dx) degrees by dx
!-----------------------------------------------------------------------
      REAL(KIND=JPRB)                 :: INTERP2D_FAST
      INTEGER(KIND=JPIM), INTENT(IN)  :: KX, KY       ! dims of 2D grid
      REAL(KIND=JPRB)   , INTENT(IN)  :: P_IN(KX,KY)  ! input matrix defined on 2D grid
      INTEGER(KIND=JPIM), INTENT(IN)  :: KDX(4)       ! vector of 4 indexes of neighbour gridpoints
      REAL(KIND=JPRB)   , INTENT(IN)  :: PRD(2)       ! ratios of x distance and y distance

      ! No DR_HOOK for performance

!-----------------------------------------------------------------------
!  Bilinear horizontal interp.
!-----------------------------------------------------------------------
      INTERP2D_FAST = (1.-PRD(1)) * (1.-PRD(2)) * P_IN(KDX(1),KDX(2)) &
     &              +     PRD(1)  * (1.-PRD(2)) * P_IN(KDX(3),KDX(2)) &
     &              + (1.-PRD(1)) *      PRD(2) * P_IN(KDX(1),KDX(4)) &
     &              +     PRD(1)  *      PRD(2) * P_IN(KDX(3),KDX(4))

  END FUNCTION INTERP2D_FAST

END SUBROUTINE BASCOE_J_INTERP

