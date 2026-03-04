SUBROUTINE RRTM_RRTM_140GP_MCICA&
 &(YDDIMV, YDEAERATM,YDERAD,YGFL,KIDIA , KFDIA , KLON , KLEV, KCOLS, KCLDCOL,&
 &PRAER, PAER  , PAPH  , PAP  , PAERTAUL, PAERASYL, PAEROMGL,&
 &PTS   , PTH   , PT   ,&
 &PEMIS , PEMIW ,&
 &PQ    , PCO2  , PCH4 , PN2O, PNO2 , PC11, PC12, PC22, PCL4, POZN   ,&
 &PCLDF , PTAUCLD, PCLFR,&
 &PEMIT , PFLUX , PFLUC,&
 &  PLWDERIVATIVE)  

! *** This program is the driver for the McICA version of RRTM_LW, 
!     the AER rapid model.  

!     For each atmosphere the user wishes to analyze, this routine
!     a) calls ECRTATM to read in the atmospheric profile 
!     b) calls SETCOEF to calculate various quantities needed for 
!        the radiative transfer algorithm
!     c) calls RTRN to do the radiative transfer calculation for
!        clear or cloudy sky
!     d) writes out the upward, downward, and net flux for each
!        level and the heating rate for each layer

!     JJMorcrette 20050110 McICA version revisited (changes in RRTM_ECRT, RRTM_RTRN)
!        NEC           25-Oct-2007 Optimisations
!     JJMorcrette 20080424 3D fields of CO2, CH4, N2O, NO2, CFC11, 12, 22 and CCL4
!     JJMorcrette 20110613 flexible number of g-points
!     P Bechtold 14/05/2012 replace ZHEATF by core constants RG*RDAY/RCPD
!                           and put arrays to scalars
!     R Hogan    20/05/2014 pass partial derivatives back to calling function
!-----------------------------------------------------------------------

USE YOEAERATM , ONLY : TEAERATM
USE YOERAD    , ONLY : TERAD
USE YOM_YGFL  , ONLY : TYPE_GFLD
USE YOMDIMV   , ONLY : TDIMV
USE PARKIND1  , ONLY : JPIM, JPRB
USE YOMHOOK   , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM   , ONLY : JPBAND, JPXSEC, JPINPX 
USE YOERRTM   , ONLY : JPGPT  
USE YOMCST    , ONLY : RG, RDAYI, RCPD

IMPLICIT NONE

!------------------------------Arguments--------------------------------

! Input arguments

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TEAERATM)    ,INTENT(INOUT) :: YDEAERATM
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
TYPE(TYPE_GFLD)   ,INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA ! First atmosphere index
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA ! Last atmosphere index
INTEGER(KIND=JPIM),INTENT(IN)    :: KCOLS ! Number of columns on which to perform RT
                                          ! should be the same as number of g-points, JPGPT
INTEGER(KIND=JPIM),INTENT(IN)    :: KCLDCOL(KLON) ! cloudy column index: 1=cloud, 0: clear    

REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAER(6,6) ! LW absorption coefficients
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAERTAUL(KLON,KLEV,16), PAERASYL(KLON,KLEV,16), PAEROMGL(KLON,KLEV,16)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (JK)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (JK)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (JK)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIS(KLON) ! Non-window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PEMIW(KLON) ! Window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: PQ(KLON,KLEV) ! H2O specific humidity (mmr)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCO2(KLON,KLEV) ! CO2 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCH4(KLON,KLEV) ! CH4 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PN2O(KLON,KLEV) ! N2O mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PNO2(KLON,KLEV) ! NO2 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC11(KLON,KLEV) ! CFC11 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC12(KLON,KLEV) ! CFC12 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PC22(KLON,KLEV) ! CFC22 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCL4(KLON,KLEV) ! CCL4  mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: POZN(KLON,KLEV) ! O3 mass mixing ratio
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLFR(KLON,KLEV)

REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) ! Surface LW emissivity
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KLON,2,KLEV+1) ! LW total sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) ! LW clear sky flux (1=up, 2=down)

! Partial derivative of total upward flux at each level with respect
! to upward flux at surface, used to correct heating rates at
! gridpoints/timesteps between calls to the full radiation scheme:
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PLWDERIVATIVE(KLON,KLEV+1)

!-- McICA ----------
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KCOLS,KLEV)    ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLD(KLON,KLEV,KCOLS)  ! Cloud optical depth

REAL(KIND=JPRB) :: ZCLDFRAC(KIDIA:KFDIA,KCOLS,YDDIMV%NFLEVG)   ! Cloud fraction
REAL(KIND=JPRB) :: ZTAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,KCOLS)    ! Spectral optical thickness
!-- McICA ----------

REAL(KIND=JPRB) :: ZATR1(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZOD(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZTF1(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZCOLDRY(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZWBRODL(KIDIA:KFDIA,YDDIMV%NFLEVG) !BROADENING GASES,column density (mol/cm2)
REAL(KIND=JPRB) :: ZCOLBRD(KIDIA:KFDIA,YDDIMV%NFLEVG) !BROADENING GASES, column amount
REAL(KIND=JPRB) :: ZWKL(KIDIA:KFDIA,JPINPX,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: ZWX(KIDIA:KFDIA,JPXSEC,YDDIMV%NFLEVG)         ! Amount of trace gases

REAL(KIND=JPRB) :: ZTOTDFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTOTDFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTOTUFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTOTUFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: JL, JK
INTEGER(KIND=JPIM) :: ISTART
INTEGER(KIND=JPIM) :: IEND

REAL(KIND=JPRB) :: ZFLUXFAC, ZHEATFAC, ZPI
REAL(KIND=JPRB) :: ZEPSEC

!- from AER
REAL(KIND=JPRB) :: ZTAUAERL(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND)

!- from INTFAC      
REAL(KIND=JPRB) :: ZFAC00(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFAC01(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFAC10(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFAC11(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from FOR
REAL(KIND=JPRB) :: ZFORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZFORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 

!- from MINOR
INTEGER(KIND=JPIM)  :: INDMINOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: ZSCALEMINOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: ZSCALEMINORN2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: ZMINORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 

REAL(KIND=JPRB)     ::&
                    &  ZRAT_H2OCO2(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2OCO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                    &  ZRAT_H2OO3(KIDIA:KFDIA,YDDIMV%NFLEVG) ,ZRAT_H2OO3_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                    &  ZRAT_H2ON2O(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2ON2O_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                    &  ZRAT_H2OCH4(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2OCH4_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                    &  ZRAT_N2OCO2(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_N2OCO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                    &  ZRAT_O3CO2(KIDIA:KFDIA,YDDIMV%NFLEVG) ,ZRAT_O3CO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from INTIND
INTEGER(KIND=JPIM) :: JP(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: JT(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: JT1(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from PRECISE             
REAL(KIND=JPRB) :: ZONEMINUS

!- from PROFDATA             
REAL(KIND=JPRB) :: ZCOLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOLCO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOLO3(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOLN2O(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOLCH4(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCOLO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZCO2MULT(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: ILAYTROP(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: ILAYSWTCH(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: ILAYLOW(KIDIA:KFDIA)

!- from PROFILE             
REAL(KIND=JPRB) :: ZPAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZPZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZTBOUND(KIDIA:KFDIA)

!- from SELF             
REAL(KIND=JPRB) :: ZSELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: ZSELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from SP             
REAL(KIND=JPRB) :: ZPFRAC(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

!- from SURFACE             
REAL(KIND=JPRB) :: ZSEMISS(KIDIA:KFDIA,JPBAND)
REAL(KIND=JPRB) :: ZSEMISLW(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: IREFLECT(KIDIA:KFDIA)

! Local variable required in case KFDIA /= KLON
REAL(KIND=JPRB) :: ZLWDERIVATIVE(KIDIA:KFDIA,KLEV+1)

LOGICAL            :: LLPRINT

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "rrtm_ecrt_140gp_mcica.intfb.h"
#include "rrtm_gasabs1a_140gp.intfb.h"
#include "rrtm_rtrn1a_140gp_mcica.intfb.h"
#include "rrtm_setcoef_140gp.intfb.h"

!     HEATFAC is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)

IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP_MCICA',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, RAER=>PRAER)
ZEPSEC = 1.E-06_JPRB
ZONEMINUS = 1.0_JPRB - ZEPSEC
ZPI = 2.0_JPRB*ASIN(1.0_JPRB)
ZFLUXFAC = ZPI * 2.E+4
!ZHEATFAC = 8.4391_JPRB
ZHEATFAC = RG*RDAYI/RCPD*1.E-2_JPRB

! *** mji ***

! For use with ECRT, this loop is over atmospheres (or longitudes)
LLPRINT=.TRUE.

!  do JK=1,KLEV
!    print 9901,JK,PT(JL,JK),PQ(JL,JK),POZN(JL,JK),PCLDF(JL,JK,1),PTAUCLD(JL,JK,1)
!  enddo

! *** mji ***
!- Prepare atmospheric profile from ECRT for use in RRTM, and define
!  other RRTM input parameters.  Arrays are passed back through the
!  existing RRTM commons and arrays.

  CALL RRTM_ECRT_140GP_MCICA&
   &(YDDIMV, YDEAERATM,YDERAD,YGFL, KIDIA, KFDIA, KLON , KLEV, KCOLS ,&
   &RAER, PAER , PAPH , PAP , PAERTAUL, PAERASYL, PAEROMGL,&
   &PTS  , PTH  , PT  ,&
   &PEMIS, PEMIW,&
   &PQ   , PCO2 , PCH4, PN2O, PNO2, PC11, PC12, PC22, PCL4, POZN , PCLDF, PTAUCLD,&
   &ZCLDFRAC, ZTAUCLD, ZCOLDRY, ZWBRODL,ZWKL, ZWX,&
   &  ZTAUAERL, ZPAVEL , ZTAVEL , ZPZ , ZTZ, ZTBOUND, ZSEMISS, IREFLECT)  

  ISTART = 1
  IEND   = 16

!  Calculate information needed by the radiative transfer routine
!  that is specific to this atmosphere, especially some of the 
!  coefficients and indices needed to compute the optical depths
!  by interpolating data from stored reference atmospheres. 

  CALL RRTM_SETCOEF_140GP&
   &(YDDIMV, KIDIA  , KFDIA    , KLEV   , ZCOLDRY  , ZWBRODL , ZWKL ,&
   &ZFAC00 , ZFAC01   , ZFAC10 , ZFAC11 , ZFORFAC,ZFORFRAC,INDFOR, JP, JT, JT1 ,&
   &ZCOLH2O, ZCOLCO2  , ZCOLO3 , ZCOLN2O, ZCOLCH4, ZCOLO2,ZCO2MULT , ZCOLBRD,&
   &ILAYTROP,ILAYSWTCH, ILAYLOW, ZPAVEL , ZTAVEL , ZSELFFAC, ZSELFFRAC, INDSELF,&
   &INDMINOR,ZSCALEMINOR,ZSCALEMINORN2,ZMINORFRAC,&
   &ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1,&
   &ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1,&
   &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)   

  CALL RRTM_GASABS1A_140GP&
   &(YDDIMV, KIDIA   , KFDIA  , KLEV, ZATR1, ZOD, ZTF1, ZPAVEL, ZCOLDRY, ZCOLBRD, ZWX ,&
   &ZTAUAERL, ZFAC00 , ZFAC01, ZFAC10 , ZFAC11 , ZFORFAC,ZFORFRAC,INDFOR, JP, JT, JT1, ZONEMINUS ,&
   &ZCOLH2O , ZCOLCO2, ZCOLO3, ZCOLN2O, ZCOLCH4, ZCOLO2,ZCO2MULT ,&
   &ILAYTROP, ILAYSWTCH,ILAYLOW, ZSELFFAC, ZSELFFRAC, INDSELF, ZPFRAC,&
   &INDMINOR,ZSCALEMINOR,ZSCALEMINORN2,ZMINORFRAC,&
   &ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1,&
   &ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1,&
   &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)      

!- Call the radiative transfer routine.

!  Clear and cloudy parts of column are treated together in RTRN.

!  print 9901,JL,ZTBOUND

  CALL RRTM_RTRN1A_140GP_MCICA&
   &(YDDIMV, KIDIA, KFDIA, KLEV, ISTART, IEND, KCOLS ,&
   &ZCLDFRAC, ZTAUCLD, ZATR1 ,&
   &ZOD , ZTF1 ,&
   &ZTOTDFLUC, ZTOTDFLUX, ZTOTUFLUC, ZTOTUFLUX ,&
   &ZTAVEL, ZTZ, ZTBOUND, ZPFRAC, ZSEMISS, ZSEMISLW ,&
   &  ZLWDERIVATIVE )  

! ***   Pass clear sky and total sky up and down flux profiles to ECRT
!       output arrays (zflux, zfluc). Array indexing from bottom to top 
!       is preserved for ECRT.
!       Invert down flux arrays for consistency with ECRT sign conventions.

DO JL = KIDIA,KFDIA

  PEMIT(JL) = ZSEMISLW(JL)
  DO JK = 0, KLEV
    PFLUC(JL,1,JK+1) =  ZTOTUFLUC(JL,JK)*ZFLUXFAC
    PFLUC(JL,2,JK+1) = -ZTOTDFLUC(JL,JK)*ZFLUXFAC
    PFLUX(JL,1,JK+1) =  ZTOTUFLUX(JL,JK)*ZFLUXFAC
    PFLUX(JL,2,JK+1) = -ZTOTDFLUX(JL,JK)*ZFLUXFAC
  ENDDO

  ! Copy to output array, noting that they may be dimensioned
  ! differently
  PLWDERIVATIVE(JL,:) = ZLWDERIVATIVE(JL,:)

ENDDO

9901 FORMAT(1X,'rrtm:',I4,12E12.5)

!------------------------------------------------------------------------

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP_MCICA',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RRTM_140GP_MCICA
