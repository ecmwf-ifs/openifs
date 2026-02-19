SUBROUTINE RRTM_TAUMOL1(YDDIMV,KIDIA,KFDIA,KLEV,P_TAU,PAVEL,&
 & P_TAUAERL,P_FAC00,P_FAC01,P_FAC10,P_FAC11,P_FORFAC,P_FORFRAC,K_INDFOR,K_JP,K_JT,K_JT1,&
 & P_COLH2O,K_LAYTROP,P_SELFFAC,P_SELFFRAC,K_INDSELF,PFRAC,P_MINORFRAC,K_INDMINOR,PSCALEMINORN2,PCOLBRD)  

!******************************************************************************
!                                                                             *
!                  Optical depths developed for the                           *
!                                                                             *
!                RAPID RADIATIVE TRANSFER MODEL (RRTM)                        *
!                                                                             *
!            ATMOSPHERIC AND ENVIRONMENTAL RESEARCH, INC.                     *
!                        840 MEMORIAL DRIVE                                   *
!                        CAMBRIDGE, MA 02139                                  *
!                                                                             *
!                           ELI J. MLAWER                                     *
!                         STEVEN J. TAUBMAN                                   *
!                         SHEPARD A. CLOUGH                                   *
!                                                                             *
!                       email:  mlawer@aer.com                                *
!                                                                             *
!        The authors wish to acknowledge the contributions of the             *
!        following people:  Patrick D. Brown, Michael J. Iacono,              *
!        Ronald E. Farren, Luke Chen, Robert Bergstrom.                       *
!                                                                             *
!******************************************************************************
!     TAUMOL                                                                  *
!                                                                             *
!     This file contains the subroutines TAUGBn (where n goes from            *
!     1 to 16).  TAUGBn calculates the optical depths and Planck fractions    *
!     per g-value and layer for band n.                                       *
!                                                                             *
!  Output:  optical depths (unitless)                                         *
!           fractions needed to compute Planck functions at every layer       *
!               and g-value                                                   *
!                                                                             *
!     COMMON /TAUGCOM/  TAUG(MXLAY,MG)                                        *
!     COMMON /PLANKG/   FRACS(MXLAY,MG)                                       *
!                                                                             *
!  Input                                                                      *
!                                                                             *
!     COMMON /FEATURES/ NG(NBANDS),NSPA(NBANDS),NSPB(NBANDS)                  *
!     COMMON /PRECISE/  ONEMINUS                                              *
!     COMMON /PROFILE/  NLAYERS,PAVEL(MXLAY),TAVEL(MXLAY),                    *
!    &                  PZ(0:MXLAY),TZ(0:MXLAY),TBOUND                        *
!     COMMON /PROFDATA/ LAYTROP,LAYSWTCH,LAYLOW,                              *
!    &                  COLH2O(MXLAY),COLCO2(MXLAY),                          *
!    &                  COLO3(MXLAY),COLN2O(MXLAY),COLCH4(MXLAY),             *
!    &                  COLO2(MXLAY),CO2MULT(MXLAY)                           *
!     COMMON /INTFAC/   FAC00(MXLAY),FAC01(MXLAY),                            *
!    &                  FAC10(MXLAY),FAC11(MXLAY)                             *
!     COMMON /INTIND/   JP(MXLAY),JT(KIDIA:KFDIA,MXLAY),JT1(KIDIA:KFDIA,MXLAY)                        *
!     COMMON /SELF/     SELFFAC(MXLAY), SELFFRAC(MXLAY), INDSELF(KIDIA:KFDIA,MXLAY)       *
!                                                                             *
!     Description:                                                            *
!     NG(IBAND) - number of g-values in band IBAND                            *
!     NSPA(IBAND) - for the lower atmosphere, the number of reference         *
!                   atmospheres that are stored for band IBAND per            *
!                   pressure level and temperature.  Each of these            *
!                   atmospheres has different relative amounts of the         *
!                   key species for the band (i.e. different binary           *
!                   species parameters).                                      *
!     NSPB(IBAND) - same for upper atmosphere                                 *
!     ONEMINUS - since problems are caused in some cases by interpolation     *
!                parameters equal to or greater than 1, for these cases       *
!                these parameters are set to this value, slightly < 1.        *
!     PAVEL - layer pressures (mb)                                            *
!     TAVEL - layer temperatures (degrees K)                                  *
!     PZ - level pressures (mb)                                               *
!     TZ - level temperatures (degrees K)                                     *
!     LAYTROP - layer at which switch is made from one combination of         *
!               key species to another                                        *
!     COLH2O, COLCO2, COLO3, COLN2O, COLCH4 - column amounts of water         *
!               vapor,carbon dioxide, ozone, nitrous ozide, methane,          *
!               respectively (molecules/cm**2)                                *
!     CO2MULT - for bands in which carbon dioxide is implemented as a         *
!               trace species, this is the factor used to multiply the        *
!               band's average CO2 absorption coefficient to get the added    *
!               contribution to the optical depth relative to 355 ppm.        *
!     FACij(JLAY) - for layer JLAY, these are factors that are needed to        *
!                  compute the interpolation factors that multiply the        *
!                  appropriate reference k-values.  A value of 0 (1) for      *
!                  i,j indicates that the corresponding factor multiplies     *
!                  reference k-value for the lower (higher) of the two        *
!                  appropriate temperatures, and altitudes, respectively.     *
!     JP - the index of the lower (in altitude) of the two appropriate        *
!          reference pressure levels needed for interpolation                 *
!     JT, JT1 - the indices of the lower of the two appropriate reference     *
!               temperatures needed for interpolation (for pressure           *
!               levels JP and JP+1, respectively)                             *
!     SELFFAC - scale factor needed to water vapor self-continuum, equals     *
!               (water vapor density)/(atmospheric density at 296K and        *
!               1013 mb)                                                      *
!     SELFFRAC - factor needed for temperature interpolation of reference     *
!                water vapor self-continuum data                              *
!     INDSELF - index of the lower of the two appropriate reference           *
!               temperatures needed for the self-continuum interpolation      *
!                                                                             *
!  Data input                                                                 *
!     COMMON /Kn/ KA(NSPA(n),5,13,MG), KB(NSPB(n),5,13:59,MG), SELFREF(10,MG) *
!        (note:  n is the band number)                                        *
!                                                                             *
!     Description:                                                            *
!     KA - k-values for low reference atmospheres (no water vapor             *
!          self-continuum) (units: cm**2/molecule)                            *
!     KB - k-values for high reference atmospheres (all sources)              *
!          (units: cm**2/molecule)                                            *
!     SELFREF - k-values for water vapor self-continuum for reference         *
!               atmospheres (used below LAYTROP)                              *
!               (units: cm**2/molecule)                                       *
!                                                                             *
!     DIMENSION ABSA(65*NSPA(n),MG), ABSB(235*NSPB(n),MG)                     *
!     EQUIVALENCE (KA,ABSA),(KB,ABSB)                                         *
!                                                                             *
!******************************************************************************

!     BAND 1:  10-250 cm-1 (low - H2O; high - H2O)
 
!     AUTHOR.
!     -------
!      JJMorcrette, ECMWF, from
!      Eli J. Mlawer, Atmospheric & Environmental Research.
!      (Revised by Michael J. Iacono, Atmospheric & Environmental Research.)

!     MODIFICATIONS.
!     --------------
!      D Salmond   2000-05-15 speed-up
!      JJMorcrette 2000-05-17 speed-up
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      NEC           25-Oct-2007 Optimisations
!      JJMorcrette 20110613 flexible number of g-points
!      ABozzo 200130517 updated to rrtmg_lw_v4.85:
!*********
!     band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                          (high key - h2o; high minor - n2)
!
!     note: previous versions of rrtm band 1: 
!           10-250 cm-1 (low - h2o; high - h2o)
! ---------------------------------------------------------------------------

USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK

USE PARRRTM  , ONLY : JPBAND
USE YOERRTM  , ONLY : JPGPT  ,NG1
USE YOERRTWN , ONLY : NSPA   ,NSPB
USE YOERRTA1 , ONLY : ABSA   ,ABSB   ,FRACREFA, FRACREFB,&
 & FORREF   ,SELFREF,  KA_MN2, KB_MN2   

IMPLICIT NONE

TYPE(TDIMV)       ,INTENT(IN)    :: YDDIMV
INTEGER(KIND=JPIM),INTENT(IN)    :: KIDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KFDIA
INTEGER(KIND=JPIM),INTENT(IN)    :: KLEV 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PAVEL(KIDIA:KFDIA,YDDIMV%NFLEVG) ! Layer pressures (hPa)
REAL(KIND=JPRB)   ,INTENT(OUT)   :: P_TAU(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_TAUAERL(KIDIA:KFDIA,YDDIMV%NFLEVG,JPBAND) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC00(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC01(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC10(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FAC11(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_FORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JP(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_JT1(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_COLH2O(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_LAYTROP(KIDIA:KFDIA) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_SELFFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: P_MINORFRAC(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDSELF(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(OUT)   :: PFRAC(KIDIA:KFDIA,JPGPT,YDDIMV%NFLEVG) 

INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDFOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
INTEGER(KIND=JPIM),INTENT(IN)    :: K_INDMINOR(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PSCALEMINORN2(KIDIA:KFDIA,YDDIMV%NFLEVG) 
REAL(KIND=JPRB)   ,INTENT(IN)    :: PCOLBRD(KIDIA:KFDIA,YDDIMV%NFLEVG)         
! ---------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: IND0(YDDIMV%NFLEVG),IND1(YDDIMV%NFLEVG),INDS(YDDIMV%NFLEVG)
INTEGER(KIND=JPIM) :: INDF(YDDIMV%NFLEVG),INDM(YDDIMV%NFLEVG)

INTEGER(KIND=JPIM) :: IG, JLAY
INTEGER(KIND=JPIM) :: JLON
REAL(KIND=JPRB) :: ZTAUFOR,ZTAUSELF,ZTAUN2,ZCORRADJ,ZPP,ZSCALEN2
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

! Minor gas mapping levels:
!     lower - n2, p = 142.5490 mbar, t = 215.70 k
!     upper - n2, p = 142.5490 mbar, t = 215.70 k

!     Compute the optical depth by interpolating in ln(pressure) and 
!     temperature.  Below LAYTROP, the water vapor self-continuum and
!     foreign continuum is interpolated (in temperature) separately.

IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL1',0,ZHOOK_HANDLE)
ASSOCIATE(NFLEVG=>YDDIMV%NFLEVG)

DO JLAY = 1, KLEV
  DO JLON = KIDIA, KFDIA
  IF (JLAY <= K_LAYTROP(JLON)) THEN
    IND0(JLAY) = ((K_JP(JLON,JLAY)-1)*5+(K_JT(JLON,JLAY)-1))*NSPA(1) + 1
    IND1(JLAY) = (K_JP(JLON,JLAY)*5+(K_JT1(JLON,JLAY)-1))*NSPA(1) + 1
    INDS(JLAY) = K_INDSELF(JLON,JLAY)
    INDF(JLAY) = K_INDFOR(JLON,JLAY)
    INDM(JLAY) = K_INDMINOR(JLON,JLAY)
    ZPP = PAVEL(JLON,JLAY) !hPa(mb)
    ZCORRADJ =  1.
    IF (ZPP < 250._JPRB) THEN
        ZCORRADJ = 1._JPRB - 0.15_JPRB * (250._JPRB-ZPP) / 154.4_JPRB
    ENDIF

    ZSCALEN2 = PCOLBRD(JLON,JLAY) * PSCALEMINORN2(JLON,JLAY)

!CDIR UNROLL=NG1
    DO IG = 1, NG1
!-- DS_000515  
           ZTAUSELF = P_SELFFAC(JLON,JLAY) * (SELFREF(INDS(JLAY),IG) + &
     &           P_SELFFRAC(JLON,JLAY) * &
     &           (SELFREF(INDS(JLAY)+1,IG) - SELFREF(INDS(JLAY),IG))) 

            ZTAUFOR =  P_FORFAC(JLON,JLAY) * (FORREF(INDF(JLAY),IG) + &
     &           P_FORFRAC(JLON,JLAY) * (FORREF(INDF(JLAY)+1,IG) - &
     &           FORREF(INDF(JLAY),IG)))  

            ZTAUN2 = ZSCALEN2*(KA_MN2(INDM(JLAY),IG) +   &
     &           P_MINORFRAC(JLON,JLAY) * &
     &           (KA_MN2(INDM(JLAY)+1,IG) - KA_MN2(INDM(JLAY),IG))) 

            P_TAU(JLON,IG,JLAY) = ZCORRADJ * (P_COLH2O(JLON,JLAY) *  &
     &          (P_FAC00(JLON,JLAY) * ABSA(IND0(JLAY),IG) + &
     &           P_FAC10(JLON,JLAY) * ABSA(IND0(JLAY)+1,IG) + &
     &           P_FAC01(JLON,JLAY) * ABSA(IND1(JLAY),IG) +  &
     &           P_FAC11(JLON,JLAY) * ABSA(IND1(JLAY)+1,IG))  &
     &           + ZTAUSELF + ZTAUFOR &
     &           + ZTAUN2) + P_TAUAERL(JLON,JLAY,1)

            PFRAC(JLON,IG,JLAY) = FRACREFA(IG)

    ENDDO
  ENDIF

  IF (JLAY > K_LAYTROP(JLON)) THEN
    IND0(JLAY) = ((K_JP(JLON,JLAY)-13)*5+(K_JT(JLON,JLAY)-1))*NSPB(1) + 1
    IND1(JLAY) = ((K_JP(JLON,JLAY)-12)*5+(K_JT1(JLON,JLAY)-1))*NSPB(1) + 1
    INDF(JLAY) = K_INDFOR(JLON,JLAY)
    INDM(JLAY) = K_INDMINOR(JLON,JLAY)
    ZPP = PAVEL(JLON,JLAY)  !hPa(mb)
    ZCORRADJ =  1._JPRB - 0.15_JPRB * (ZPP / 95.6_JPRB)

    ZSCALEN2 = PCOLBRD(JLON,JLAY) * PSCALEMINORN2(JLON,JLAY)

!-- JJM000517
!CDIR UNROLL=NG1
    DO IG = 1, NG1
!-- JJM000517
            ZTAUFOR = P_FORFAC(JLON,JLAY) * (FORREF(INDF(JLAY),IG) + &
     &           P_FORFRAC(JLON,JLAY) * &
     &           (FORREF(INDF(JLAY)+1,IG) - FORREF(INDF(JLAY),IG))) 

            ZTAUN2 = ZSCALEN2*(KB_MN2(INDM(JLAY),IG) +  &
     &           P_MINORFRAC(JLON,JLAY) * &
     &           (KB_MN2(INDM(JLAY)+1,IG) - KB_MN2(INDM(JLAY),IG)))

            P_TAU(JLON,IG,JLAY) = ZCORRADJ * (P_COLH2O(JLON,JLAY) * &
     &          (P_FAC00(JLON,JLAY) * ABSB(IND0(JLAY),IG) + &
     &           P_FAC10(JLON,JLAY) * ABSB(IND0(JLAY)+1,IG) + &
     &           P_FAC01(JLON,JLAY) * ABSB(IND1(JLAY),IG) + &
     &           P_FAC11(JLON,JLAY) * ABSB(IND1(JLAY)+1,IG)) &  
     &           + ZTAUFOR &
     &           + ZTAUN2)+ P_TAUAERL(JLON,JLAY,1)
      PFRAC(JLON,IG,JLAY) = FRACREFB(IG)


    ENDDO
  ENDIF
  ENDDO
ENDDO

END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('RRTM_TAUMOL1',1,ZHOOK_HANDLE)

END SUBROUTINE RRTM_TAUMOL1
