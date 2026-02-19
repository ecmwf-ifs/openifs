!***************************************************************************
!                                                                          *
!                RRTM :  RAPID RADIATIVE TRANSFER MODEL                    *
!                                                                          *
!             ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                 *
!                        840 MEMORIAL DRIVE                                *
!                        CAMBRIDGE, MA 02139                               *
!                                                                          *
!                           ELI J. MLAWER                                  *
!                         STEVEN J. TAUBMAN~                               *
!                         SHEPARD A. CLOUGH                                *
!                                                                          *
!                        ~currently at GFDL                                *
!                                                                          *
!                       email:  mlawer@aer.com                             *
!                                                                          *
!        The authors wish to acknowledge the contributions of the          *
!        following people:  Patrick D. Brown, Michael J. Iacono,           *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                    *
!                                                                          *
!***************************************************************************
!     Reformatted for F90 by JJMorcrette, ECMWF, 980714                    * 
!                                                                          *
!***************************************************************************
! *** mji ***
! *** This version of RRTM has been altered to interface with either
!     the ECMWF numerical weather prediction model or the ECMWF column 
!     radiation model (ECRT) package. 

!     Revised, April, 1997;  Michael J. Iacono, AER, Inc.
!          - initial implementation of RRTM in ECRT code
!     Revised, June, 1999;  Michael J. Iacono and Eli J. Mlawer, AER, Inc.
!          - to implement generalized maximum/random cloud overlap
!        NEC           25-Oct-2007 Optimisations
!        D. Salmond    11-Dec-2007 Optimizations
!        NEC/FC        05-Oct-2009 Optimisation

SUBROUTINE RRTM_RRTM_140GP&
 & (YDDIMV, YDERAD,KIDIA , KFDIA , KLON , KLEV,&
 & PRAER, PAER  , PAPH  , PAP,&
 & PTS   , PTH   , PT,&
 & P_ZEMIS , P_ZEMIW,&
 & PQ    , PCO2  , PCH4, PN2O, PNO2, PC11, PC12, PC22, PCL4, POZN,&
 & PCLDF , PTAUCLD,&
 & PEMIT , PFLUX , PFLUC, PTCLEAR&
 & )  

! *** This program is the driver for RRTM, the AER rapid model.  
!     For each atmosphere the user wishes to analyze, this routine
!     a) calls ECRTATM to read in the atmospheric profile 
!     b) calls SETCOEF to calculate various quantities needed for 
!        the radiative transfer algorithm
!     c) calls RTRN to do the radiative transfer calculation for
!        clear or cloudy sky
!     d) writes out the upward, downward, and net flux for each
!        level and the heating rate for each layer

!     JJMorcrette 20080424 3D fields of CO2, CH4, N2O, NO2, CFC11, 12, 22, and CCL4
!     JJMorcrette 20110613 flexible number of g-points
!     P Bechtold 14/05/2012 replace ZHEATF by core constants RG*RDAY/RCPD

USE YOERAD   , ONLY : TERAD
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARRRTM  , ONLY : JPBAND, JPXSEC, JPINPX
USE YOERRTM  , ONLY : JPGPT   
USE YOMCST   , ONLY : RG, RDAYI, RCPD

IMPLICIT NONE

!------------------------------Arguments--------------------------------

! Input arguments

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
TYPE(TERAD)       ,INTENT(INOUT) :: YDERAD
INTEGER(KIND=JPIM),INTENT(IN)    :: KLON! Number of atmospheres (longitudes) 
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV! Number of atmospheric layers 
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA ! First atmosphere index
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA ! Last atmosphere index
REAL(KIND=JPRB)   ,INTENT(IN)    :: PRAER(6,6) ! LW absorption coefficients
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAER(KLON,6,KLEV) ! Aerosol optical thickness
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAPH(KLON,KLEV+1) ! Interface pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAP(KLON,KLEV) ! Layer pressures (Pa)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTS(KLON) ! Surface temperature (JLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTH(KLON,KLEV+1) ! Interface temperatures (JLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: PT(KLON,KLEV) ! Layer temperature (JLEV)
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIS(KLON) ! Non-window surface emissivity
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_ZEMIW(KLON) ! Window surface emissivity
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
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCLDF(KLON,KLEV) ! Cloud fraction
REAL(KIND=JPRB)   ,INTENT(IN)    :: PTAUCLD(KLON,KLEV,JPBAND) ! Cloud optical depth
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PEMIT(KLON) ! Surface LW emissivity
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUX(KLON,2,KLEV+1) ! LW total sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFLUC(KLON,2,KLEV+1) ! LW clear sky flux (1=up, 2=down)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PTCLEAR(KLON) ! clear-sky fraction of column
INTEGER(KIND=JPIM) :: ICLDLYR(KIDIA:KFDIA,YDDIMV%NFLEVG)        ! Cloud indicator
REAL(KIND=JPRB) :: Z_CLDFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)           ! Cloud fraction
REAL(KIND=JPRB) :: Z_TAUCLD(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND)     ! Spectral optical thickness

REAL(KIND=JPRB) :: Z_ATR1  (KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

!REAL(KIND=JPRB) :: Z_OD    (JPGPT,NFLEVG,KIDIA:KFDIA)
REAL(KIND=JPRB) :: Z_OD    (KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_TF1   (KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_COLDRY(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_WBRODL(KIDIA:KFDIA,YDDIMV%NFLEVG) !BROADENING GASES
REAL(KIND=JPRB) :: Z_COLBRD(KIDIA:KFDIA,YDDIMV%NFLEVG) !BROADENING GASES,column amount 
REAL(KIND=JPRB) :: Z_WKL(KIDIA:KFDIA,JPINPX,YDDIMV%NFLEVG)

REAL(KIND=JPRB) :: Z_WX(KIDIA:KFDIA,JPXSEC,YDDIMV%NFLEVG)         ! Amount of trace gases

REAL(KIND=JPRB) :: Z_TOTDFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TOTDFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TOTUFLUC(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TOTUFLUX(KIDIA:KFDIA,0:YDDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: ICLD(KIDIA:KFDIA), JLON, JLEV, JI
INTEGER(KIND=JPIM) :: ISTART
INTEGER(KIND=JPIM) :: IEND

REAL(KIND=JPRB) :: Z_FLUXFAC, Z_HEATFAC, Z_PI, ZEPSEC, ZTCLEAR(KIDIA:KFDIA)

!- from AER
REAL(KIND=JPRB) :: Z_TAUAERL(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND)

!- from INTFAC      
REAL(KIND=JPRB) :: Z_FAC00(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_FAC01(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_FAC10(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_FAC11(KIDIA:KFDIA,YDDIMV%NFLEVG)
!- from FOR
REAL(KIND=JPRB) :: Z_FORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_FORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 

!- from MINOR
INTEGER(KIND=JPIM)  :: INDMINOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: Z_SCALEMINOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: Z_SCALEMINORN2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)     :: Z_MINORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 

REAL(KIND=JPRB)     ::&
                     & ZRAT_H2OCO2(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2OCO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                     & ZRAT_H2OO3(KIDIA:KFDIA,YDDIMV%NFLEVG) ,ZRAT_H2OO3_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                     & ZRAT_H2ON2O(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2ON2O_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                     & ZRAT_H2OCH4(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_H2OCH4_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                     & ZRAT_N2OCO2(KIDIA:KFDIA,YDDIMV%NFLEVG),ZRAT_N2OCO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG),&
                     & ZRAT_O3CO2(KIDIA:KFDIA,YDDIMV%NFLEVG) ,ZRAT_O3CO2_1(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from INTIND
INTEGER(KIND=JPIM) :: JP(KIDIA,KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: JT(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: JT1(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from PRECISE             
REAL(KIND=JPRB) :: Z_ONEMINUS

!- from PROFDATA             
REAL(KIND=JPRB) :: Z_COLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_COLCO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_COLO3 (KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_COLN2O(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_COLCH4(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_COLO2(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_CO2MULT(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: I_LAYTROP(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: I_LAYSWTCH(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: I_LAYLOW(KIDIA:KFDIA)

!- from PROFILE             
REAL(KIND=JPRB) :: Z_PAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_PZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TZ(KIDIA:KFDIA,0:YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_TBOUND(KIDIA:KFDIA)
INTEGER(KIND=JPIM) :: I_NLAYERS(KIDIA:KFDIA)

!- from SELF             
REAL(KIND=JPRB) :: Z_SELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
REAL(KIND=JPRB) :: Z_SELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG)

!- from SP             
REAL(KIND=JPRB) :: Z_PFRAC(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG)

!- from SURFACE             
REAL(KIND=JPRB) :: Z_SEMISS(KIDIA:KFDIA,JPBAND)
REAL(KIND=JPRB) :: Z_SEMISLW(KIDIA:KFDIA)
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "rrtm_ecrt_140gp.intfb.h"
#include "rrtm_gasabs1a_140gp.intfb.h"
#include "rrtm_rtrn1a_140gp.intfb.h"
#include "rrtm_setcoef_140gp.intfb.h"

!     HEATFAC is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)

IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG, RAER=>PRAER)
ZEPSEC = 1.E-06_JPRB
Z_ONEMINUS = 1.0_JPRB - ZEPSEC
Z_PI = 2.0_JPRB*ASIN(1.0_JPRB)
Z_FLUXFAC = Z_PI * 2.D4
!Z_HEATFAC = 8.4391_JPRB
Z_HEATFAC = RG*RDAYI/RCPD*1.E-2_JPRB

! *** mji ***
! For use with ECRT, this loop is over atmospheres (or longitudes)
DO JLON = KIDIA,KFDIA

! *** mji ***
!- Prepare atmospheric profile from ECRT for use in RRTM, and define
!  other RRTM input parameters.  Arrays are passed back through the
!  existing RRTM commons and arrays.
  ZTCLEAR(JLON)=1.0_JPRB
ENDDO

  CALL RRTM_ECRT_140GP&
   & (YDDIMV,YDERAD,KLON,KIDIA,KFDIA, KLEV, ICLD,&
   & RAER, PAER , PAPH , PAP,&
   & PTS  , PTH  , PT,&
   & P_ZEMIS, P_ZEMIW,&
   & PQ   , PCO2, PCH4, PN2O, PNO2, PC11, PC12, PC22, PCL4, POZN, PCLDF, PTAUCLD, ZTCLEAR,&
   & Z_CLDFRAC,Z_TAUCLD,Z_COLDRY, Z_WBRODL,Z_WKL,Z_WX,&
   & Z_TAUAERL,Z_PAVEL,Z_TAVEL,Z_PZ,Z_TZ,Z_TBOUND,I_NLAYERS,Z_SEMISS)  

DO JLON = KIDIA,KFDIA
  PTCLEAR(JLON)=ZTCLEAR(JLON)
ENDDO

  ISTART = 1
  IEND   = 16

!  Calculate information needed by the radiative transfer routine
!  that is specific to this atmosphere, especially some of the 
!  coefficients and indices needed to compute the optical depths
!  by interpolating data from stored reference atmospheres. 

  CALL RRTM_SETCOEF_140GP(YDDIMV,KIDIA,KFDIA,KLEV,Z_COLDRY,Z_WBRODL, Z_WKL,&
   & Z_FAC00,Z_FAC01,Z_FAC10,Z_FAC11,Z_FORFAC,Z_FORFRAC,INDFOR,JP,JT,JT1,&
   & Z_COLH2O,Z_COLCO2,Z_COLO3,Z_COLN2O,Z_COLCH4,Z_COLO2,Z_CO2MULT, Z_COLBRD,&
   & I_LAYTROP,I_LAYSWTCH,I_LAYLOW,Z_PAVEL,Z_TAVEL,Z_SELFFAC,Z_SELFFRAC,INDSELF,&
   & INDMINOR,Z_SCALEMINOR,Z_SCALEMINORN2,Z_MINORFRAC,&
   & ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1,&
   & ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1,&
   &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)     

  CALL RRTM_GASABS1A_140GP(YDDIMV,KIDIA,KFDIA,KLEV,Z_ATR1,Z_OD,Z_TF1,Z_PAVEL,Z_COLDRY,Z_COLBRD,Z_WX,&
   & Z_TAUAERL,Z_FAC00,Z_FAC01,Z_FAC10,Z_FAC11,Z_FORFAC,Z_FORFRAC,INDFOR,JP,JT,JT1,Z_ONEMINUS,&
   & Z_COLH2O,Z_COLCO2,Z_COLO3,Z_COLN2O,Z_COLCH4,Z_COLO2,Z_CO2MULT,&
   & I_LAYTROP,I_LAYSWTCH,I_LAYLOW,Z_SELFFAC,Z_SELFFRAC,INDSELF,Z_PFRAC,&
   & INDMINOR,Z_SCALEMINOR,Z_SCALEMINORN2,Z_MINORFRAC,&
   & ZRAT_H2OCO2, ZRAT_H2OCO2_1, ZRAT_H2OO3, ZRAT_H2OO3_1,&
   & ZRAT_H2ON2O, ZRAT_H2ON2O_1, ZRAT_H2OCH4, ZRAT_H2OCH4_1,&
   &  ZRAT_N2OCO2, ZRAT_N2OCO2_1, ZRAT_O3CO2, ZRAT_O3CO2_1)      
 

!- Call the radiative transfer routine.

! *** mji ***
!  Check for cloud in column.  Use ECRT threshold set as flag icld in
!  routine ECRTATM.  If icld=1 then column is cloudy, otherwise it is
!  clear.  Also, set up flag array, icldlyr, for use in radiative
!  transfer.  Set icldlyr to one for each layer with non-zero cloud
!  fraction.

DO JLON = KIDIA,KFDIA
  DO JLEV = 1, KLEV
    IF (ICLD(JLON) == 1.AND.Z_CLDFRAC(JLON,JLEV) > ZEPSEC) THEN
      ICLDLYR(JLON,JLEV) = 1
    ELSE
      ICLDLYR(JLON,JLEV) = 0
    ENDIF
  ENDDO
ENDDO

!  Clear and cloudy parts of column are treated together in RTRN.
!  Clear radiative transfer is done for clear layers and cloudy radiative
!  transfer is done for cloudy layers as identified by icldlyr.

  CALL RRTM_RTRN1A_140GP(YDDIMV,KIDIA,KFDIA,KLEV,ISTART,IEND,ICLDLYR,Z_CLDFRAC,Z_TAUCLD,Z_ATR1,&
   & Z_OD,Z_TF1,Z_TOTDFLUC,Z_TOTDFLUX,Z_TOTUFLUC,Z_TOTUFLUX,&
   & Z_TAVEL,Z_TZ,Z_TBOUND,Z_PFRAC,Z_SEMISS,Z_SEMISLW)

! ***   Pass clear sky and total sky up and down flux profiles to ECRT
!       output arrays (zflux, zfluc). Array indexing from bottom to top 
!       is preserved for ECRT.
!       Invert down flux arrays for consistency with ECRT sign conventions.

DO JLON = KIDIA,KFDIA
  PEMIT(JLON) = Z_SEMISLW(JLON)
ENDDO

DO JLON = KIDIA,KFDIA
  DO JI = 0, KLEV
    PFLUC(JLON,1,JI+1) =  Z_TOTUFLUC(JLON,JI)*Z_FLUXFAC
    PFLUC(JLON,2,JI+1) = -Z_TOTDFLUC(JLON,JI)*Z_FLUXFAC
    PFLUX(JLON,1,JI+1) =  Z_TOTUFLUX(JLON,JI)*Z_FLUXFAC
    PFLUX(JLON,2,JI+1) = -Z_TOTDFLUX(JLON,JI)*Z_FLUXFAC
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_RRTM_140GP',1,ZHOOK_HANDLE)
END SUBROUTINE RRTM_RRTM_140GP
