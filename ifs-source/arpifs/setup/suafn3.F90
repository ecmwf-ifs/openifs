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

SUBROUTINE SUAFN3(YDAFN)

!**** *SUAFN3*  - PRINT OUT ARPEGE FIELD DESCRIPTORS

!     PURPOSE.
!     --------
!        TO PRINT OUT ARPEGE FIELD DESCRIPTORS

!**   INTERFACE.
!     ----------
!       *CALL* *SUAFN3*

!        EXPLICIT ARGUMENTS
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS
!        --------------------
!        See #include below.

!     METHOD.
!     -------
!        See documentation about FULL-POS.

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     AUTHOR.
!     -------
!      RYAD EL KHATIB *METEO-FRANCE*
!      ORIGINAL : 98-04-08 from SUAFN

!     MODIFICATIONS.
!     --------------
!      R. El Khatib : 01-05-04 Duration of total precipitations,HUn,HUx
!      R. El Khatib : 01-08-07 Pruning options
!      R. El Khatib : 02-21-20 Fullpos B-level distribution + remove IO scheme
!      F. Vana      : 02-07-03 Virtual potential temperature
!      R. El Khatib : 02-09-04 Iso -10 Celsius
!      R. El Khatib : 03-04-17 Fullpos improvemnts
!      M.Hamrud : 03-10-01 CY28 Cleaning
!      F. Taillefer : 04-05-26 Add ozone and aerosols fields
!      Y. Bouteloup : 04-03-26 Introduction of MTS radiance
!      A. Tompkins  : 11-02-04 TCLW, TCLW
!      R. Engelen   : 05-03-15 Greenhouse gases (ECMWF)
!      A. Untch     : 05-03-15 Aerosols (ECMWF)
!      J. Flemming  : 05-04-11 Aerosols replaced with reactive gases
!      Y. Bouteloup : 05-01-20 Introduction of atmospheric rain
!      G. Hello     : 05-01-25 Adding snow, rain,graupel, tke
!      R. El Khatib 10-May-2005 IMASK instead of LLLSM
!      R. El Ouaraini: 05-07-28 Adding fields for the monitoring of update frequency of the coupling files for Aladin
!      A. Alias     : 05-12-21 Add sulfate and volcano aerosols fields
!      M. Bellus    : 28-Sep-2006 Add prognostic convection fields (ALARO-0)
!      JJMorcrette  : 10-05-2006 MODIS albedo
!      S. Serrar    : 05-06-23 CO2 surface fluxes added
!      S. Serrar    : 05-09-26 Add CO2 and SF6 total column
!      R. Engelen   : 02-Jun-06  CO2  replaced by generic GHG
!      JJMorcrette 20060721 PP of clear-sky PAR and TOA incident solar radiation
!      JJMorcrette 20060807 PP of vertically integerated mass divergence VIMD
!      JJMorcrette 20060925 DU, BC, OM, SO2, VOL, SOA climatol.surf.fields 
!      G. Balsamo  20070115 Soil type
!      Y. Seity    : 07-03-06 Add Graupel and Hail
!      P. Lopez      : 26-06-2009 Added 100m wind components
!      P. Lopez      : 22-10-2008 Added ducting diagnostic fields
!      Y. Bouteloup  : 01-05-2009 Add Total (radiative) cloud water and ice
!      S. Boussetta/G.Balsamo      May 2009   Added variable LAI field
!      G. Balsamo  20090807 Add  surface/subsurface runoff
!      K. Yessad (Aug 2009): name TFP_EXT for extra-GFL.
!      JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface 
!      H. Hersbach 04-Dec-2009 10-m neutral wind and friction velocity
!      JJMorcrette 20100212 PP of CBASE, 0DEGL and VISIH
!      R. Forbes   01-Mar-10 Added TCRW, TCSW diagnostics
!      T. Wilhelmsson: 25-Mar-2010 Add 6 hourly min/max fields
!      T. Wilhelmsson: 20-Dec-2010 Add 3 hourly min/max fields
!      G.Balsamo/S.Boussetta: 17-Apr-2011 Add land carbon dioxide
!      P.Bechtold    : 09-Aug-2011 Add CIN, convective Indices
!      M Ahlgrimm  31 Oct 2011 Clear-sky downward radiation at surface
!      R. El Khatib : 01-Mar-2012 LFPOS => NFPOS
!      A. Inness  23 Mar 2012 Add total column chemistry fields TCCHEM
!      JJMorcrette 20130213 PP optical depths GEMS/MACC aerosols
!      R. El Khatib : 01-Mar-2012 LFPOS => NFPOS
!      G. Balsamo   : 11-Jan-2014 Add lake fields
!      R. Forbes :   01-Mar-2014 Added precip rates/type, TCSLW,I10FG,PEV
!      JJMorcrette 20130730 15-variable aerosol model + oceanic DMS
!      R. El Khatib 17-Jul-2013 FABEC post-processing
!      Y. Bouteloup : 03-Oct-2013 Add LCONV, ICONV, RCONV and SCONV
!      A.Agusti-Panareda: 30-Oct-2013 Add GPP/REC flux adjustment coeffs 
!      P. Marguinaud: 01-Oct-2014 Setup of PREP fields to be interpolated, move
!      PRINT_GFP in a separate routine of re-use
!      P. Marguinaud: 10-Oct-2014 More fields
!      K. Yessad (July 2014): Move some variables.
!      R. Forbes :   10-Jan-2015 Add freezing rain FZRA
!      S.Remy : 28-Jan-2015 Add injection height for bb emissions  
!      P. Lopez 16-Nov-2015 Added lightning density field
!      R. El Khatib 08-Dec-2015 Interoperability GRIB2 vs FA
!      A. Bozzo Jan 2016 Add surface shortwave direct beam
!      S.Remy : 09-Mar-2016 Add aer SO2 dry dep velocity
!      S.Remy : 21-Apr-2017 Add altitude of volcanoes
!      S. Rémy, 25/4/2017, add fraction of calcite over dust 1 and 2
!      R. El Khatib 17-Aug-2016 Better interoperability GRIB2 vs FA
!      R. El Khatib 12-Mar-2018 total albedo
!      R. Brozkova Sep-2018 Added convective temperature, global normal
!                           irradiance and mean radiant temperature
!      R. Hogan     15-Jan-2019 6-component MODIS albedo
!      A.Agusti-Panareda : 23-06-2021 Add CO2 photosynthesis type (C3/C4)
!      R. Forbes    May-2022  Added precip-type most-frequent and most-severe
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARFPOS  , ONLY : JPOSSCVA, JPOSVX2, JPOSFSU, JPOSAERO, JPOSGHG, JPOSCHEM, &
 & JPOSCHEMFLX, JPOSEZDIAG, JPOSAERODIAG, JPOSAERO_WVL_DIAG, JPOSAERO_WVL_DIAG_TYPES, &
 & JPOSEMIS2D, JPOSEMIS2DAUX, JPOSEMIS3D
USE YOMLUN   , ONLY : NULOUT
USE YOMMP0   , ONLY : NPRINTLEV
USE YOMCT0   , ONLY : LARPEGEF
USE YOMAFN, ONLY : TAFN
USE FULLPOS_MIX,ONLY: FULLPOS_TYPE
USE YOE_AERODIAG, ONLY : CPAERODIAG_LABEL, CPAERO_WVL_DIAG_LABEL

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE (TAFN),  INTENT(IN) :: YDAFN

INTEGER(KIND=JPIM) :: J,JVAR,JDIAG, IERR, ISIZE, ILIMIT

CHARACTER(LEN=25) :: CLTEXT
CHARACTER(LEN=34) :: CLTEXU
CHARACTER(LEN=12) :: CLVARN
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "abor1.intfb.h"
#include "print_gfp.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUAFN3',0,ZHOOK_HANDLE)
ASSOCIATE(TFP=>YDAFN%TFP, TFP_DYNDS=>YDAFN%TFP_DYNDS, GFP=>YDAFN%GFP, GFP_PHYDS=>YDAFN%GFP_PHYDS, &
 & NSFXPRE_CNT=>YDAFN%NSFXPRE_CNT)
!     ------------------------------------------------------------------
!*       1. Print out descriptors
!           ---------------------

IF (NPRINTLEV >=1) THEN

  ILIMIT=3

  WRITE(UNIT=NULOUT,FMT='('' MODULE YOMAFN'')')

  WRITE(UNIT=NULOUT,FMT='(/,'' UPPER AIR FIELDS FROM DYNAMIC :'')')

  CALL PRINT_TFP

  CALL PRINT_TFP('Geopotential.............','TFP_Z   ',TFP%Z   )
  CALL PRINT_TFP('Temperature..............','TFP_T   ',TFP%T   )
  CALL PRINT_TFP('U-momentum of wind.......','TFP_U   ',TFP%U   )
  CALL PRINT_TFP('V-momentum of wind.......','TFP_V   ',TFP%V   )
  CALL PRINT_TFP('Specific humidity........','TFP_Q   ',TFP%Q   )
  CALL PRINT_TFP('Relative humidity........','TFP_HU  ',TFP%HU  )
  CALL PRINT_TFP('Ver velocity omega=DPi/Dt','TFP_VV  ',TFP%VV  )
  CALL PRINT_TFP('Etadot...................','TFP_ETAD',TFP%ETAD)
  CALL PRINT_TFP('Vorticity................','TFP_VOR ',TFP%VOR )
  CALL PRINT_TFP('Divergence...............','TFP_DIV ',TFP%DIV )
  CALL PRINT_TFP('Velocity potential.......','TFP_PSI ',TFP%PSI )
  CALL PRINT_TFP('Stream function..........','TFP_KHI ',TFP%KHI )
  CALL PRINT_TFP('Atmospheric liquid water.','TFP_L   ',TFP%L   )
  CALL PRINT_TFP('Atmospheric solid water..','TFP_I   ',TFP%I   )
  CALL PRINT_TFP('Atmos total liquid water.','TFP_LRAD',TFP%LRAD)
  CALL PRINT_TFP('Atmos total solid water..','TFP_IRAD',TFP%IRAD)
  CALL PRINT_TFP('Atmospheric snow.........','TFP_SN  ',TFP%SN  )
  CALL PRINT_TFP('Atmospheric rain.........','TFP_RR  ',TFP%RR  )
  CALL PRINT_TFP('Atmospheric graupel......','TFP_GR  ',TFP%GR  )
  CALL PRINT_TFP('Atmospheric hail.........','TFP_HL  ',TFP%HL  )
  CALL PRINT_TFP('Turbulent kinetic energy.','TFP_TKE ',TFP%TKE )
  CALL PRINT_TFP('Mass change rate in phys.','TFP_PHYCTY ',TFP%PHYCTY )
  CALL PRINT_TFP('First variable EFB scheme','TFP_EFB1',TFP%EFB1)
  CALL PRINT_TFP('Secon variable EFB scheme','TFP_EFB2',TFP%EFB2)
  CALL PRINT_TFP('Third variable EFB scheme','TFP_EFB3',TFP%EFB3)
  CALL PRINT_TFP('Potential temperature....','TFP_TH  ',TFP%TH  )
  CALL PRINT_TFP('Wet bulb pot. temperature','TFP_THPW',TFP%THPW)
  CALL PRINT_TFP('Wet bulb temperature.....','TFP_TPW' ,TFP%TPW )
  CALL PRINT_TFP('Cloud fraction...........','TFP_CLF ',TFP%CLF )
  CALL PRINT_TFP('Ozone....................','TFP_O3MX',TFP%O3MX)
  CALL PRINT_TFP('Convective precip. flux..','TFP_CPF ',TFP%CPF)
  CALL PRINT_TFP('Stratiform precip. flux..','TFP_SPF ',TFP%SPF)
  CALL PRINT_TFP('Wind velocity............','TFP_WND ',TFP%WND )
  CALL PRINT_TFP('Equiv. pot. temperature..','TFP_ETH ',TFP%ETH )
  CALL PRINT_TFP('Isobaric equivalent tempe','TFP_IET ',TFP%IET )
  CALL PRINT_TFP('Simulated reflect. in mmh','TFP_SRE ',TFP%SRE )
  CALL PRINT_TFP('Simulated reflect. in dBz','TFP_SREDB ',TFP%SREDB )
  CALL PRINT_TFP('Virtual Theta............','TFP_THV ',TFP%THV )
  CALL PRINT_TFP('Absolute Vorticity.......','TFP_ABS ',TFP%ABS )
  CALL PRINT_TFP('Stretching Deformation...','TFP_STD ',TFP%STD )
  CALL PRINT_TFP('Potential Vorticity......','TFP_PV  ',TFP%PV  )
  CALL PRINT_TFP('Shearing Deformation.....','TFP_SHD ',TFP%SHD )
  CALL PRINT_TFP('Pressure.................','TFP_P   ',TFP%P   )
  CALL PRINT_TFP('Pressure Departure.......','TFP_PD  ',TFP%PD  )
  CALL PRINT_TFP('Montgomery potential.....','TFP_MG  ',TFP%MG  )
  CALL PRINT_TFP('Pseudo Vertic. Divergence','TFP_VD  ',TFP%VD  )
  CALL PRINT_TFP('Vert velocity w=Dz/Dt....','TFP_VW  ',TFP%VW  )
  CALL PRINT_TFP('Pressure of iso-t........','TFP_PTB ',TFP%PTB )
  CALL PRINT_TFP('altitude of iso-t........','TFP_HTB ',TFP%HTB )
  CALL PRINT_TFP('Humid density............','TFP_RHO ',TFP%RHO )

! Prognostic convection
  CALL PRINT_TFP('DD Mesh Fraction.........','TFP_DAL ',TFP%DAL )
  CALL PRINT_TFP('DD Vertical Velocity.....','TFP_DOM ',TFP%DOM )
  CALL PRINT_TFP('UD Mesh Fraction.........','TFP_UAL ',TFP%UAL )
  CALL PRINT_TFP('UD Vertical Velocity.....','TFP_UOM ',TFP%UOM )
  CALL PRINT_TFP('UD pseudohist entrainment','TFP_UEN ',TFP%UEN )
  CALL PRINT_TFP('Cnv pseudohist cloud frac','TFP_UNEBH ',TFP%UNEBH )
  CALL PRINT_TFP('Total turbulent energy...','TFP_TTE ',TFP%TTE )
  CALL PRINT_TFP('Prognostic mixing length.','TFP_MXL ',TFP%MXL )
  CALL PRINT_TFP('Shear source term turbul.','TFP_SHTUR',TFP%SHTUR )
  CALL PRINT_TFP('Flux Q source term turbul','TFP_FQTUR',TFP%FQTUR)
  CALL PRINT_TFP('Flux S source term turbul','TFP_FSTUR',TFP%FSTUR)
  CALL PRINT_TFP('RK scheme temper tend ...','TFP_RKTH ',TFP%RKTH )
  CALL PRINT_TFP('RK scheme vapor  tend ...','TFP_RKTQV',TFP%RKTQV)
  CALL PRINT_TFP('RK scheme condens tend ..','TFP_RKTQC',TFP%RKTQC)
  CALL PRINT_TFP('Lifting Condensation Level','TFP_LCL ',TFP%LCL )
  CALL PRINT_TFP('Free Convection Level     ','TFP_FCL ',TFP%FCL )
  CALL PRINT_TFP('Equilibrium Level         ','TFP_EL  ',TFP%EL )
  CALL PRINT_TFP('Convective Liquid Water...','TFP_LCONV',TFP%LCONV )
  CALL PRINT_TFP('Convective Solid Water....','TFP_ICONV',TFP%ICONV )
  CALL PRINT_TFP('Convective Rain...........','TFP_RCONV',TFP%RCONV )
  CALL PRINT_TFP('Convective Snow...........','TFP_SCONV',TFP%SCONV )
  CALL PRINT_TFP('Eddy diffusivity rate.....','TFP_EDR  ',TFP%EDR)
  DO JVAR=1,JPOSEMIS3D
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'EMIS3D nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A8,I2.2,A2)') 'TFP_EMIS3D(',JVAR,') '
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%EMIS3D(JVAR))
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSEMIS3D  
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSGHG
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'GHG nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A8,I2.2,A2)') 'TFP_GHG(',JVAR,') '
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%GHG(JVAR))
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSGHG  
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSEZDIAG
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'EZDIAG nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') 'TFP_EZDIAG(',JVAR,')'
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%EZDIAG(JVAR))
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') 100
      EXIT
    ENDIF
  ENDDO
  DO JVAR=1,JPOSAERO
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'Aerosol nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') 'TFP_AERO(',JVAR,')'
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%AERO(JVAR))
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSAERO
      EXIT
    ENDIF
  ENDDO
  CALL PRINT_TFP('CH4 loss rate.............','TFP_LRCH4',TFP%LRCH4)
  DO JVAR=1,JPOSCHEM
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'CHEM var nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') 'TFP_CHEM(',JVAR,')'
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%CHEM(JVAR))
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSCHEM 
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSSCVA
    WRITE(CLTEXT,FMT='(A18,I2.2,A5)') 'Passive scalar nr ',JVAR,'.....'
    WRITE(CLVARN,FMT='(A8,I2.2,A2)') 'TFP_EXT(',JVAR,') '
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%EXT(JVAR))
    IF (JVAR == ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSSCVA  
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSVX2
    WRITE(CLTEXT,FMT='(A23,I2.2)') 'Free upper air field n ',JVAR
    WRITE(CLVARN,FMT='(A8,I2.2,A2)') 'TFP_FUA(',JVAR,') '
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%FUA(JVAR))
    IF (JVAR == ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSVX2  
      EXIT
    ENDIF    
  ENDDO

  WRITE(UNIT=NULOUT,FMT='(/,'' SURFACE FIELDS FROM DYNAMIC :'')')

  CALL PRINT_TFP('Surface pressure.........','TFP_SP  ',TFP%SP  )
  CALL PRINT_TFP('Surface NH pressure......','TFP_SPNH',TFP%SPNH)
  CALL PRINT_TFP('Mean sea level pressure..','TFP_MSL ',TFP%MSL )
  CALL PRINT_TFP('Filtred ln(Ps),NCUFNR=1..','TFP_CUF1',TFP%CUF1)
  CALL PRINT_TFP('Filtred ln(Ps),NCUFNR=2  ','TFP_CUF2',TFP%CUF2)
  CALL PRINT_TFP('Filtred ln(Ps),NCUFNR=3..','TFP_CUF3',TFP%CUF3)
  CALL PRINT_TFP('Filtred ln(Ps),NCUFNR=4..','TFP_CUF4',TFP%CUF4)
  CALL PRINT_TFP('Filtred ln(Ps),NCUFNR=5..','TFP_CUF5',TFP%CUF5)
  CALL PRINT_TFP('Surface geopotential.....','TFP_FIS ',TFP%FIS )
  CALL PRINT_TFP('Map factor...............','TFP_GM  ',TFP%GM  )
  CALL PRINT_TFP('Tropo. Folding Indicator.','TFP_FOL ',TFP%FOL )
  CALL PRINT_TFP('U-momentum of ICAO jet...','TFP_UJET',TFP%UJET)
  CALL PRINT_TFP('V-momentum of ICAO jet...','TFP_VJET',TFP%VJET)
  CALL PRINT_TFP('ICAO jet pressure........','TFP_PJET',TFP%PJET)
  CALL PRINT_TFP('ICAO Tropopause pressure.','TFP_PCAO',TFP%PCAO)
  CALL PRINT_TFP('ICAO Tropo. temperature..','TFP_TCAO',TFP%TCAO)
  CALL PRINT_TFP('altitude of iso-tprimw=0 ','TFP_HTPW',TFP%HTPW)
  CALL PRINT_TFP('altitude of iso-tprimw=1 ','TFP_HTPW1',TFP%HTPW1)
  CALL PRINT_TFP('altitude of iso-tpw=1P5..','TFP_HTPW2',TFP%HTPW2)
  CALL PRINT_TFP('Surface Vertical Velocity','TFP_WWS ',TFP%WWS )
  CALL PRINT_TFP('Log. of Surface pressure.','TFP_LNSP',TFP%LNSP)
  CALL PRINT_TFP('U cls....................','TFP_UCLS',TFP%UCLS)
  CALL PRINT_TFP('V cls....................','TFP_VCLS',TFP%VCLS)
  CALL PRINT_TFP('T cls....................','TFP_TCLS',TFP%TCLS)
  CALL PRINT_TFP('Q cls....................','TFP_QCLS',TFP%QCLS)
  CALL PRINT_TFP('HU cls...................','TFP_RCLS',TFP%RCLS)
  CALL PRINT_TFP('Module of wind cls.......','TFP_FCLS',TFP%FCLS)
  CALL PRINT_TFP('Maxi. temperature in cls.','TFP_TX  ',TFP%TX  )
  CALL PRINT_TFP('Mini. temperature in cls.','TFP_TN  ',TFP%TN  )
  CALL PRINT_TFP('Maxi. rel. moist. in cls.','TFP_HUX ',TFP%HUX )
  CALL PRINT_TFP('Mini. rel. moist. in cls.','TFP_HUN ',TFP%HUN )
  CALL PRINT_TFP('Maxi. simu. reflect(mm/h)','TFP_SREX',TFP%SREX)
  CALL PRINT_TFP('Maxi. simu. reflect (dBZ)','TFP_SREDBX',TFP%SREDBX)
  CALL PRINT_TFP('Pressure of echotop......','TFP_TOPR',TFP%TOPR)
  CALL PRINT_TFP('CAPE.....................','TFP_CAPE',TFP%CAPE)
  CALL PRINT_TFP('CIEN.....................','TFP_CIEN',TFP%CIEN)
  CALL PRINT_TFP('Temperature of convection','TFP_TCVS',TFP%TCVS)
  CALL PRINT_TFP('Storm relative helicity..','TFP_SRH',TFP%SRH)
  CALL PRINT_TFP('Storm motion, u component','TFP_STRMMU',TFP%STRMMU)
  CALL PRINT_TFP('Storm motion, v component','TFP_STRMMU',TFP%STRMMU)
  CALL PRINT_TFP('MOCON....................','TFP_MOCO',TFP%MOCO)
  CALL PRINT_TFP('Total water vapour.......','TFP_TWV ',TFP%TWV )
  CALL PRINT_TFP('U gusts..................','TFP_UGST',TFP%UGST)
  CALL PRINT_TFP('V gusts..................','TFP_VGST',TFP%VGST)
  CALL PRINT_TFP('Module of gusts..........','TFP_FGST',TFP%FGST)
  CALL PRINT_TFP('Height of the PBL........','TFP_HCLP',TFP%HCLP)
  CALL PRINT_TFP('Ventilation Index........','TFP_VEIN',TFP%VEIN)
  CALL PRINT_TFP('Analysed surface MOCON...','TFP_ASMC',TFP%ASMC)
  CALL PRINT_TFP('Surface MOCON for VARPACK','TFP_VSMC',TFP%VSMC)
  CALL PRINT_TFP('Forecasted surface MOCON.','TFP_SMC' ,TFP%SMC )
  CALL PRINT_TFP('QNH pressure.............','TFP_QNH' ,TFP%QNH )
  CALL PRINT_TFP('MSL NH pressure..........','TFP_MSLNH' ,TFP%MSLNH )
  CALL PRINT_TFP('MSAT7 MVIRI CH 1.........','TFP_MSAT7C1',TFP%MSAT7C1)
  CALL PRINT_TFP('MSAT7 MVIRI CH 2.........','TFP_MSAT7C2',TFP%MSAT7C2)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 1........','TFP_MSAT8C1',TFP%MSAT8C1)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 2........','TFP_MSAT8C2',TFP%MSAT8C2)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 3........','TFP_MSAT8C3',TFP%MSAT8C3)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 4........','TFP_MSAT8C4',TFP%MSAT8C4)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 5........','TFP_MSAT8C5',TFP%MSAT8C5)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 6........','TFP_MSAT8C6',TFP%MSAT8C6)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 7........','TFP_MSAT8C7',TFP%MSAT8C7)
  CALL PRINT_TFP('MSAT8 SEVIRI CH 8........','TFP_MSAT8C8',TFP%MSAT8C8)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 1........','TFP_MSAT9C1',TFP%MSAT9C1)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 2........','TFP_MSAT9C2',TFP%MSAT9C2)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 3........','TFP_MSAT9C3',TFP%MSAT9C3)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 4........','TFP_MSAT9C4',TFP%MSAT9C4)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 5........','TFP_MSAT9C5',TFP%MSAT9C5)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 6........','TFP_MSAT9C6',TFP%MSAT9C6)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 7........','TFP_MSAT9C7',TFP%MSAT9C7)
  CALL PRINT_TFP('MSAT9 SEVIRI CH 8........','TFP_MSAT9C8',TFP%MSAT9C8)
  CALL PRINT_TFP('GOES11 IMAGER CH 1.......','TFP_GOES11C1',TFP%GOES11C1)
  CALL PRINT_TFP('GOES11 IMAGER CH 2.......','TFP_GOES11C2',TFP%GOES11C2)
  CALL PRINT_TFP('GOES11 IMAGER CH 3.......','TFP_GOES11C3',TFP%GOES11C3)
  CALL PRINT_TFP('GOES11 IMAGER CH 4.......','TFP_GOES11C4',TFP%GOES11C4)
  CALL PRINT_TFP('GOES12 IMAGER CH 1.......','TFP_GOES12C1',TFP%GOES12C1)
  CALL PRINT_TFP('GOES12 IMAGER CH 2.......','TFP_GOES12C2',TFP%GOES12C2)
  CALL PRINT_TFP('GOES12 IMAGER CH 3.......','TFP_GOES12C3',TFP%GOES12C3)
  CALL PRINT_TFP('GOES12 IMAGER CH 4.......','TFP_GOES12C4',TFP%GOES12C4)
  CALL PRINT_TFP('MTSAT 1 IMAGER CH 1......','TFP_MTSAT1C1',TFP%MTSAT1C1)
  CALL PRINT_TFP('MTSAT 1 IMAGER CH 2......','TFP_MTSAT1C2',TFP%MTSAT1C2)
  CALL PRINT_TFP('MTSAT 1 IMAGER CH 3......','TFP_MTSAT1C3',TFP%MTSAT1C3)
  CALL PRINT_TFP('MTSAT 1 IMAGER CH 4......','TFP_MTSAT1C4',TFP%MTSAT1C4)
  CALL PRINT_TFP('Mask extra domain........','TFP_MSK' ,TFP%MSK ) 
  DO JVAR=1,JPOSFSU
    WRITE(CLTEXT,FMT='(A22,I2.2,A1)') 'Free surface field nr ',JVAR,'.'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') ' TFP_FSU(',JVAR,')'
    CALL PRINT_TFP(CLTEXT,CLVARN,TFP%FSU(JVAR))
    IF (JVAR == ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSFSU  
      EXIT
    ENDIF
  ENDDO

  WRITE(UNIT=NULOUT,FMT='(/,'' SURFACE FIELDS FROM PHYSIC :'')')

  CALL PRINT_GFP

  CALL PRINT_GFP('SD FILTERED OROG .................','GFP_SDFOR',GFP%SDFOR)
  CALL PRINT_GFP('GPP flux adjustment coefficients. ','GFP_CGPP ',GFP%CGPP )
  CALL PRINT_GFP('REC flux adjustment coefficients. ','GFP_CREC ',GFP%CREC )
  CALL PRINT_GFP('Land/sea mask.....................','GFP_LSM  ',GFP%LSM  )
  CALL PRINT_GFP('OUTPUT gp orography (*g)..........','GFP_GFIS ',GFP%GFIS )
  CALL PRINT_GFP('INTERPOLATED gp orography (*g)....','GFP_SFIS ',GFP%SFIS )
  CALL PRINT_GFP('Surface temperature...............','GFP_ST   ',GFP%ST   )
  CALL PRINT_GFP('Deep soil temperature.............','GFP_DST  ',GFP%DST  )
  CALL PRINT_GFP('INTERPOLATED Ts...................','GFP_RDST ',GFP%RDST )
  CALL PRINT_GFP('Surface soil wetness..............','GFP_SSW  ',GFP%SSW  )
  CALL PRINT_GFP('Deep soil wetness.................','GFP_DSW  ',GFP%DSW  )
  CALL PRINT_GFP('Relax. deep soil wetness..........','GFP_RDSW ',GFP%RDSW )
  CALL PRINT_GFP('Clim. rel. surf soil wetness......','GFP_CSSW ',GFP%CSSW )
  CALL PRINT_GFP('Clim. rel. deep soil wetness......','GFP_CDSW ',GFP%CDSW )
  CALL PRINT_GFP('Snow depth........................','GFP_SD   ',GFP%SD   )
  CALL PRINT_GFP('Snow depth tot....................','GFP_SDSL ',GFP%SDSL )
  CALL PRINT_GFP('Surface roughness (*g)............','GFP_SR   ',GFP%SR   )
  CALL PRINT_GFP('Roughness length (*g).............','GFP_BSR  ',GFP%BSR  )
  CALL PRINT_GFP('Albedo............................','GFP_AL   ',GFP%AL   )
  CALL PRINT_GFP('Emissivity........................','GFP_EMIS ',GFP%EMIS )
  CALL PRINT_GFP('Soil type.........................','GFP_SOTY ',GFP%SOTY )
  CALL PRINT_GFP('Std. dev. of orography (*g).......','GFP_SDOG ',GFP%SDOG )
  CALL PRINT_GFP('Lake cover........................','GFP_CLK  ',GFP%CLK  )
  CALL PRINT_GFP('Lake depth........................','GFP_DL   ',GFP%DL   )
  CALL PRINT_GFP('Lake mix layer temperature........','GFP_LMLT ',GFP%LMLT )
  CALL PRINT_GFP('Lake mix layer depth..............','GFP_LMLD ',GFP%LMLD )
  CALL PRINT_GFP('Lake bottom layer temperature.....','GFP_LBLT ',GFP%LBLT )
  CALL PRINT_GFP('Lake total layer temperature......','GFP_LTLT ',GFP%LTLT )
  CALL PRINT_GFP('Lake shape factor.................','GFP_LSHF ',GFP%LSHF )
  CALL PRINT_GFP('Lake ice temperature..............','GFP_LICT ',GFP%LICT )    
  CALL PRINT_GFP('Lake ice depth....................','GFP_LICD ',GFP%LICD )
  CALL PRINT_GFP('percentage of vegetation..........','GFP_VEG  ',GFP%VEG  )
  CALL PRINT_GFP('percentage of land................','GFP_LAN  ',GFP%LAN  )
  CALL PRINT_GFP('Anisotropy of topography..........','GFP_ACOT ',GFP%ACOT )
  CALL PRINT_GFP('Topography main direction.........','GFP_DPAT ',GFP%DPAT )
  CALL PRINT_GFP('Index of vegetation...............','GFP_IVEG ',GFP%IVEG )
  CALL PRINT_GFP('Stomatal minimum resistance.......','GFP_RSMIN',GFP%RSMIN)
  CALL PRINT_GFP('Percentage of clay in soil........','GFP_ARG  ',GFP%ARG  )
  CALL PRINT_GFP('Percentage of sand in soil........','GFP_SAB  ',GFP%SAB  )
  CALL PRINT_GFP('Soil depth........................','GFP_D2   ',GFP%D2   )
  CALL PRINT_GFP('Leaf area index...................','GFP_LAI  ',GFP%LAI  )
  CALL PRINT_GFP('Resist. to evapotranspir..........','GFP_HV   ',GFP%HV   )
  CALL PRINT_GFP('Heat roughness length (*g)........','GFP_Z0H  ',GFP%Z0H  )
  CALL PRINT_GFP('Interception content..............','GFP_IC   ',GFP%IC   )
  CALL PRINT_GFP('Surface snow albedo...............','GFP_ALSN ',GFP%ALSN )
  CALL PRINT_GFP('Surface snow density..............','GFP_SNDE ',GFP%SNDE )
  CALL PRINT_GFP('Surf. albedo for non snowed areas.','GFP_BAAL ',GFP%BAAL )
  CALL PRINT_GFP('Surf. albedo total................','GFP_ALBHIS',GFP%ALBHIS)
  CALL PRINT_GFP('Surf. albedo for bare ground......','GFP_ALS  ',GFP%ALS  )
  CALL PRINT_GFP('Surf. albedo for vegetation.......','GFP_ALV  ',GFP%ALV  )
  CALL PRINT_GFP('Frozen deep soil wetness..........','GFP_FDSW ',GFP%FDSW )
  CALL PRINT_GFP('Frozen surface soil wetness.......','GFP_FSSW ',GFP%FSSW )
  CALL PRINT_GFP('Surface relative moisture.........','GFP_PSRHU',GFP%PSRHU)
  CALL PRINT_GFP('Anisotropy vector/U-momentum......','GFP_PADOU',GFP%PADOU)
  CALL PRINT_GFP('Anisotropy vector/V-momentum......','GFP_PADOV',GFP%PADOV)
  CALL PRINT_GFP('Analysed RMS of PHI (CANARI)......','GFP_PCAAG',GFP%PCAAG)
  CALL PRINT_GFP('Forecasted RMS of PHI (CANARI)....','GFP_PCAPG',GFP%PCAPG)
  CALL PRINT_GFP('Interpolated dynamic surf. g*Z0...','GFP_IDZ0 ',GFP%IDZ0 )
  CALL PRINT_GFP('Interpolated thermal surf. g*Z0...','GFP_ITZ0 ',GFP%ITZ0 )
  CALL PRINT_GFP('Maximum proportion of vegetation..','GFP_PVGMX',GFP%PVGMX)
  CALL PRINT_GFP('Vegetation roughness length (*g)..','GFP_Z0V  ',GFP%Z0V  )
  CALL PRINT_GFP('Proportion of urbanisation........','GFP_PURB ',GFP%PURB )
  CALL PRINT_GFP('Maximum soil depth................','GFP_D2MX ',GFP%D2MX )
  CALL PRINT_GFP('Marine aerosols...................','GFP_ASEA ',GFP%ASEA )
  CALL PRINT_GFP('Continental aerosols..............','GFP_ALAN ',GFP%ALAN )
  CALL PRINT_GFP('Carbone aerosols..................','GFP_ASOO ',GFP%ASOO )
  CALL PRINT_GFP('Desert aerosols...................','GFP_ADES ',GFP%ADES )
  CALL PRINT_GFP('Sulfate aerosols..................','GFP_ASUL ',GFP%ASUL )
  CALL PRINT_GFP('Volcano aerosols..................','GFP_AVOL ',GFP%AVOL )
  CALL PRINT_GFP('First ozone profile (A)...........','GFP_O3A  ',GFP%O3A  )
  CALL PRINT_GFP('Second ozone profile (B)..........','GFP_O3B  ',GFP%O3B  )
  CALL PRINT_GFP('Third ozone profile (C)...........','GFP_O3C  ',GFP%O3C  )
  CALL PRINT_GFP('Vertically-integr. mass divergence','GFP_VIMD ',GFP%VIMD )

  DO JVAR=1,JPOSVX2
    WRITE(CLTEXT,FMT='(A14,I2.2,A9)') 'Extra field n ',JVAR,'.........'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') ' GFP_VX2(',JVAR,')'
    CALL PRINT_GFP(CLTEXT,CLVARN,GFP%VX2(JVAR))
    IF (JVAR == ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSVX2  
      EXIT
    ENDIF
  ENDDO
  DO JVAR=1,JPOSFSU
    WRITE(CLTEXT,FMT='(A14,I2.2,A9)') 'Free field nr ',JVAR,'.........'
    WRITE(CLVARN,FMT='(A9,I2.2,A1)') ' GFP_FSU(',JVAR,')'
    CALL PRINT_GFP(CLTEXT,CLVARN,GFP%FSU(JVAR))
    IF (JVAR == ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSFSU  
      EXIT
    ENDIF
  ENDDO

  CALL PRINT_GFP('WIND GUST 10 M  ..................','GFP_10FG ',GFP%EC10FG )
  CALL PRINT_GFP('WIND GUST 10 M SINCE LAST 3 HOURS.','GFP_10FG3',GFP%EC10FG3)
  CALL PRINT_GFP('WIND GUST 10 M SINCE LAST 6 HOURS.','GFP_10FG6',GFP%EC10FG6)
  CALL PRINT_GFP('WIND GUST 10 M (INSTANTANEOUS)....','GFP_I10FG',GFP%I10FG)
  CALL PRINT_GFP('10M.UWIND       ..................','GFP_10U  ',GFP%EC10U  )
  CALL PRINT_GFP('10M.VWIND       ..................','GFP_10V  ',GFP%EC10V  )
  CALL PRINT_GFP('10M.WIND SPEED  ..................','GFP_10SI ',GFP%EC10SI )
  CALL PRINT_GFP('100M.UWIND      ..................','GFP_100U ',GFP%EC100U )
  CALL PRINT_GFP('100M.VWIND      ..................','GFP_100V ',GFP%EC100V )
  CALL PRINT_GFP('100M.WIND SPEED ..................','GFP_100SI',GFP%EC100SI)
  CALL PRINT_GFP('200M.UWIND      ..................','GFP_200U ',GFP%EC200U )
  CALL PRINT_GFP('200M.VWIND      ..................','GFP_200V ',GFP%EC200V )
  CALL PRINT_GFP('200M.WIND SPEED ..................','GFP_200SI',GFP%EC200SI)
  CALL PRINT_GFP('FRICTION VELOCITY      ...........','GFP_ZUST ',GFP%ZUST )
  CALL PRINT_GFP('10M.NEUTRAL UWIND      ...........','GFP_10NU ',GFP%EC10NU )
  CALL PRINT_GFP('10M.NEUTRAL VWIND      ...........','GFP_10NV ',GFP%EC10NV )
  CALL PRINT_GFP('2M.DEWPOINT     ..................','GFP_2D   ',GFP%EC2D   )
  CALL PRINT_GFP('2M.SPECIFIC HUMIDITY..............','GFP_2SH  ',GFP%EC2SH  )
  CALL PRINT_GFP('2M.TEMPERATURE  ..................','GFP_2T   ',GFP%EC2T   )
  CALL PRINT_GFP('ALBEDO          ..................','GFP_ALB  ',GFP%ALB  )

  CALL PRINT_GFP('MODIS Albedo UVis parallel rad....','GFP_ALUVP',GFP%ALUVP)
  CALL PRINT_GFP('MODIS Albedo UVis diffuse rad.....','GFP_ALUVD',GFP%ALUVD)
  CALL PRINT_GFP('MODIS Albedo N-IR parallel rad....','GFP_ALNIP',GFP%ALNIP)
  CALL PRINT_GFP('MODIS Albedo N-IR diffuse rad.....','GFP_ALNID',GFP%ALNID)

  CALL PRINT_GFP('MODIS Albedo UVis parallel rad.iso','GFP_ALUVI',GFP%ALUVI)
  CALL PRINT_GFP('MODIS Albedo N-IR parallel rad.iso','GFP_ALNII',GFP%ALNII)
  CALL PRINT_GFP('MODIS Albedo UVis parallel rad.vol','GFP_ALUVV',GFP%ALUVV)
  CALL PRINT_GFP('MODIS Albedo N-IR parallel rad.vol','GFP_ALNIV',GFP%ALNIV)
  CALL PRINT_GFP('MODIS Albedo UVis parallel rad.geo','GFP_ALUVG',GFP%ALUVG)
  CALL PRINT_GFP('MODIS Albedo N-IR parallel rad.geo','GFP_ALNIG',GFP%ALNIG)
!--
  CALL PRINT_GFP('Dust emission potential...........','GFP_AERDEP',GFP%AERDEP)
  CALL PRINT_GFP('Lifting Threshold Speed...........','GFP_AERLTS',GFP%AERLTS)
  CALL PRINT_GFP('Soil Clay Content.................','GFP_AERSCC',GFP%AERSCC)
  CALL PRINT_GFP('Dust Source Function..............','GFP_DSF   ',GFP%DSF)
  CALL PRINT_GFP('Dust size distrib. modulation.....','GFP_DSZ   ',GFP%DSZ)

  CALL PRINT_GFP('SO2 dry deposition velocity.......','GFP_SO2DD',GFP%SO2DD )
  CALL PRINT_GFP('Oceanic DMS.......................','GFP_DMSO ',GFP%DMSO )
  CALL PRINT_GFP('Urban fraction....................','GFP_URBF ',GFP%URBF )
  CALL PRINT_GFP('Calcite Fraction 1 (small/medium).','GFP_FCA1 ',GFP%FCA1 )
  CALL PRINT_GFP('Calcite Fraction 2 (large)........','GFP_FCA2 ',GFP%FCA2 )

  CALL PRINT_GFP('Optical depth sea salt aerosols...','GFP_ODSS ',GFP%ODSS )
  CALL PRINT_GFP('Optical depth dust aerosols.......','GFP_ODDU ',GFP%ODDU )
  CALL PRINT_GFP('Optical depth organic M aerosols..','GFP_ODOM ',GFP%ODOM )
  CALL PRINT_GFP('Optical depth black C aerosols....','GFP_ODBC ',GFP%ODBC )
  CALL PRINT_GFP('Optical depth sulphate aerosols...','GFP_ODSU ',GFP%ODSU )
  CALL PRINT_GFP('Optical depth nitrate aerosols....','GFP_ODNI ',GFP%ODNI )
  CALL PRINT_GFP('Optical depth ammonium aerosols...','GFP_ODAM ',GFP%ODAM )
  CALL PRINT_GFP('Optical depth secondary org. aer..','GFP_ODSOA',GFP%ODSOA)
  CALL PRINT_GFP('Optical depth volc.fly ash aeros..','GFP_ODVFA',GFP%ODVFA)
  CALL PRINT_GFP('Optical depth volc.sulfate aeros..','GFP_ODVSU',GFP%ODVSU)
  CALL PRINT_GFP('Optical depth tot. aer. acc.......','GFP_ODTOACC',GFP%ODTOACC)
  CALL PRINT_GFP('Particulate matter le 1 um .......','GFP_AEPM1 ',GFP%AEPM1 )
  CALL PRINT_GFP('Particulate matter le 2.5um.......','GFP_AEPM25',GFP%AEPM25)
  CALL PRINT_GFP('Particulate matter le 10um .......','GFP_AEPM10',GFP%AEPM10)
  CALL PRINT_GFP('UV Biologically Effective Dose....','GFP_UVBED',GFP%UVBED)
  CALL PRINT_GFP('UV Biologically Effective Dose CS.','GFP_UVBEDCS',GFP%UVBEDCS)
!-- 
  CALL PRINT_GFP('ANGLE.SUB.SOROG ..................','GFP_ANOR ',GFP%ANOR )
  CALL PRINT_GFP('SNOW ALBEDO     ..................','GFP_ASN  ',GFP%ASN  )
  CALL PRINT_GFP('PBL.DISSIP.     ..................','GFP_BLD  ',GFP%BLD  )
  CALL PRINT_GFP('BOUND.LAY.HEIGHT..................','GFP_BLH  ',GFP%BLH  )
  CALL PRINT_GFP('BUDGET VALUES   ..................','GFP_BV   ',GFP%BV   )
  CALL PRINT_GFP('CLOUD BASE LEVEL..................','GFP_CBASE',GFP%CBASE)
  CALL PRINT_GFP('CONV.CC         ..................','GFP_CCC  ',GFP%CCC  )
  CALL PRINT_GFP('CHARNOCK.PARAM  ..................','GFP_CHAR ',GFP%CHAR )
  CALL PRINT_GFP('SEA ICE COVER   ..................','GFP_CI   ',GFP%CI   )
  CALL PRINT_GFP('CONV.PRECIP.    ..................','GFP_CP   ',GFP%CP   )
  CALL PRINT_GFP('CONV.SNOWFALL   ..................','GFP_CSF  ',GFP%CSF  )
  CALL PRINT_GFP('HIGH VEG. COVER ..................','GFP_CVH  ',GFP%CVH  )
  CALL PRINT_GFP('LOW VEG. COVER  ..................','GFP_CVL  ',GFP%CVL  )
  CALL PRINT_GFP('URBAN COVER     ..................','GFP_CUR  ',GFP%CUR  )
  CALL PRINT_GFP('CO2 PHOTOS. TYPE..................','GFP_COTYP',GFP%CO2TYP)
  CALL PRINT_GFP('WETLAND FRACTION..................','GFP_FWET ',GFP%FWET )
  CALL PRINT_GFP('LOW VEG. LAI    ..................','GFP_LAIL ',GFP%LAIL ) 
  CALL PRINT_GFP('HIGH VEG. LAI   ..................','GFP_LAIH ',GFP%LAIH ) 
  CALL PRINT_GFP('EVAPORATION     ..................','GFP_E    ',GFP%E    )
  CALL PRINT_GFP('EVA. SNOW       ..................','GFP_ES   ',GFP%ES   )
  CALL PRINT_GFP('U-STRESS        ..................','GFP_EWSS ',GFP%EWSS )
  CALL PRINT_GFP('EVA. SNOW       ..................','GFP_ES   ',GFP%ES   )
  CALL PRINT_GFP('U-STRESS        ..................','GFP_EWSS ',GFP%EWSS )
  CALL PRINT_GFP('FREEZINGRAIN    ..................','GFP_FZRA ',GFP%FZRA )
  CALL PRINT_GFP('GRAV.WAVE.DISSIP..................','GFP_GWD  ',GFP%GWD  )
  CALL PRINT_GFP('HIGH.CC         ..................','GFP_HCC  ',GFP%HCC  )
  CALL PRINT_GFP('INST.S.MOI.FLUX ..................','GFP_IE   ',GFP%IE   )
  CALL PRINT_GFP('INST.X-SURF.ST. ..................','GFP_IEWSS',GFP%IEWSS)
  CALL PRINT_GFP('INST.X-SURF.ST. ..................','GFP_IEWSS',GFP%IEWSS)
  CALL PRINT_GFP('INST.Y-SURF.ST. ..................','GFP_INSSS',GFP%INSSS)
  CALL PRINT_GFP('ANISO.SUB.SOROG ..................','GFP_ISOR ',GFP%ISOR )
  CALL PRINT_GFP('INST.S.HEAT.FLUX..................','GFP_ISSHF',GFP%ISSHF)
  CALL PRINT_GFP('ICE SURF.TEMP. 1..................','GFP_ISTL1',GFP%ISTL1)
  CALL PRINT_GFP('ICE SURF.TEMP. 2..................','GFP_ISTL2',GFP%ISTL2)
  CALL PRINT_GFP('ICE SURF.TEMP. 3..................','GFP_ISTL3',GFP%ISTL3)
  CALL PRINT_GFP('ICE SURF.TEMP. 4..................','GFP_ISTL4',GFP%ISTL4)
  CALL PRINT_GFP('LOW.CC          ..................','GFP_LCC  ',GFP%LCC  )
  CALL PRINT_GFP('LAT.GRAV.WAVE.ST..................','GFP_LGWS ',GFP%LGWS )
  CALL PRINT_GFP('TOTAL PRECIPITATION RATE .........','GFP_TPR  ',GFP%TPR  )
  CALL PRINT_GFP('LARGE SCALE RAINFALL RATE ........','GFP_LSRR ',GFP%LSRR )
  CALL PRINT_GFP('CONVECTIVE RAINFALL RATE .........','GFP_CRR  ',GFP%CRR  )
  CALL PRINT_GFP('LARGE SCALE SNOWFALL RATE ........','GFP_LSSFR',GFP%LSSFR)
  CALL PRINT_GFP('CONVECTIVE SNOWFALL RATE .........','GFP_CSFR ',GFP%CSFR )
  CALL PRINT_GFP('LARGE.SCALE.SNOW..................','GFP_LSF  ',GFP%LSF  )
  CALL PRINT_GFP('LARGESCAPRECIP. ..................','GFP_LSP  ',GFP%LSP  )
  CALL PRINT_GFP('LSC PREC. FRACT. ACCUMULATED......','GFP_LSPF ',GFP%LSPF )
  CALL PRINT_GFP('LSC PREC. FRACT. INSTANTANEOUS....','GFP_ILSPF',GFP%ILSPF)
  CALL PRINT_GFP('LOG.SURF.ROUGH  ..................','GFP_LSRH ',GFP%LSRH )
  CALL PRINT_GFP('LOG.S.ROUGH.HEAT..................','GFP_LZ0H ',GFP%LZ0H )
  CALL PRINT_GFP('MEDIUM.CC       ..................','GFP_MCC  ',GFP%MCC  )
  CALL PRINT_GFP('MER.GRAV.WAVE.ST..................','GFP_MGWS ',GFP%MGWS )
  CALL PRINT_GFP('MIN.TEMP.       ..................','GFP_MN2T ',GFP%MN2T )
  CALL PRINT_GFP('MIN.TEMP. SINCE LAST 3 HOURS......','GFP_MN2T3',GFP%MN2T3)
  CALL PRINT_GFP('MIN.TEMP. SINCE LAST 6 HOURS......','GFP_MN2T6',GFP%MN2T6)
  CALL PRINT_GFP('MIN.TOT.PRECIP. ..................','GFP_MNTPR ',GFP%MNTPR )
  CALL PRINT_GFP('MIN.TOT.PRECIP.SINCE LAST 6 HOURS.','GFP_MNTPR6',GFP%MNTPR6)
  CALL PRINT_GFP('MSL             ..................','GFP_MSLD ',GFP%MSLD )
  CALL PRINT_GFP('SURFACE PRESSURE..................','GFP_SP   ',GFP%SP   )
  CALL PRINT_GFP('MAX.TEMP.       ..................','GFP_MX2T ',GFP%MX2T )
  CALL PRINT_GFP('MAX.TEMP. SINCE LAST 3 HOURS......','GFP_MX2T3',GFP%MX2T3)
  CALL PRINT_GFP('MAX.TEMP. SINCE LAST 6 HOURS......','GFP_MX2T6',GFP%MX2T6)
  CALL PRINT_GFP('MAX.TOT.PRECIP. ..................','GFP_MXTPR ',GFP%MXTPR )
  CALL PRINT_GFP('MAX.TOT.PRECIP.SINCE LAST 3 HOURS.','GFP_MXTPR3',GFP%MXTPR3)
  CALL PRINT_GFP('MAX.TOT.PRECIP.SINCE LAST 6 HOURS.','GFP_MXTPR6',GFP%MXTPR6)
  CALL PRINT_GFP('V-STRESS        ..................','GFP_NSSS ',GFP%NSSS )
  CALL PRINT_GFP('POTENTIAL EVAPORATION  ...........','GFP_PEV  ',GFP%PEV  )
  CALL PRINT_GFP('PRECIPITATION TYPE  ..............','GFP_PTYPE',GFP%PTYPE)
  CALL PRINT_GFP('RUNOFF          ..................','GFP_RO   ',GFP%RO   )
  CALL PRINT_GFP('SURFACE RUNOFF  ..................','GFP_SRO  ',GFP%SRO  )
  CALL PRINT_GFP('SUB-SURFACE RUNOFF................','GFP_SSRO ',GFP%SSRO )
  CALL PRINT_GFP('SNOW DENSITY    ..................','GFP_RSN  ',GFP%RSN  )
  CALL PRINT_GFP('SNOWFALL        ..................','GFP_SF   ',GFP%SF   )
  CALL PRINT_GFP('SKINTEMPERATURE ..................','GFP_SKT  ',GFP%SKT  )
  CALL PRINT_GFP('S.L.HEAT.FLUX.  ..................','GFP_SLHF ',GFP%SLHF )
  CALL PRINT_GFP('SLOPE.SUB.SOROG ..................','GFP_SLOR ',GFP%SLOR )
  CALL PRINT_GFP('SLOPE.SUB.SOROG ..................','GFP_SLOR ',GFP%SLOR )
  CALL PRINT_GFP('SNOW MELT       ..................','GFP_SMLT ',GFP%SMLT )
  CALL PRINT_GFP('SKINRESERV.EAU  ..................','GFP_SRC  ',GFP%SRC  )
  CALL PRINT_GFP('S.S.HEAT.FLUX.  ..................','GFP_SSHF ',GFP%SSHF )
  CALL PRINT_GFP('SURF SOLAR CLEAR..................','GFP_SSRC ',GFP%SSRC )
  CALL PRINT_GFP('S.SOLAR.RAD.DOWN..................','GFP_SSRD ',GFP%SSRD )
  CALL PRINT_GFP('S.SOLAR.RAD.DOWN CLEAR............','GFP_SSRDC',GFP%SSRDC)
  CALL PRINT_GFP('SEA SURF. TEMP. ..................','GFP_SST  ',GFP%SST  )
  CALL PRINT_GFP('OCEAN U-CURRENT ..................','GFP_UCUR ',GFP%UCUR )
  CALL PRINT_GFP('OCEAN V-CURRENT ..................','GFP_VCUR ',GFP%VCUR )
  CALL PRINT_GFP('LEV1TEMPERATURE ..................','GFP_STL1 ',GFP%STL1 )
  CALL PRINT_GFP('LEV2TEMPERATURE ..................','GFP_STL2 ',GFP%STL2 )
  CALL PRINT_GFP('LEV3TEMPERATURE ..................','GFP_STL3 ',GFP%STL3 )
  CALL PRINT_GFP('LEV4TEMPERATURE ..................','GFP_STL4 ',GFP%STL4 )
  CALL PRINT_GFP('S.THERMAL.RAD.  ..................','GFP_STR  ',GFP%STR  )
  CALL PRINT_GFP('SURF THERM CLEAR..................','GFP_STRC ',GFP%STRC )
  CALL PRINT_GFP('S.THERM.RAD.DOWN..................','GFP_STRD ',GFP%STRD )
  CALL PRINT_GFP('S.THERM.RAD.DOWN CLEAR............','GFP_STRDC',GFP%STRDC)
  CALL PRINT_GFP('SUNSHI. DURATION..................','GFP_SUND ',GFP%SUND )
  CALL PRINT_GFP('LEV1RESERV.EAU  ..................','GFP_SWL1 ',GFP%SWL1 )
  CALL PRINT_GFP('LEV2RESERV.EAU  ..................','GFP_SWL2 ',GFP%SWL2 )
  CALL PRINT_GFP('LEV3RESERV.EAU  ..................','GFP_SWL3 ',GFP%SWL3 )
  CALL PRINT_GFP('LEV4RESERV.EAU  ..................','GFP_SWL4 ',GFP%SWL4 )
  CALL PRINT_GFP('TOTAL.CC        ..................','GFP_TCC  ',GFP%TCC  )
  CALL PRINT_GFP('TOT.COL.OZON    ..................','GFP_TCO3 ',GFP%TCO3 )
  CALL PRINT_GFP('NET ECOSYSTEM EXCHANGE OF CO2.....','GFP_NEE  ',GFP%NEE  )
  CALL PRINT_GFP('GROSS PRIMARY PRODUCTION OF CO2...','GFP_GPP  ',GFP%GPP  )
  CALL PRINT_GFP('ECOSYSTEM RESPIRATION OF CO2......','GFP_REC  ',GFP%REC  )
  CALL PRINT_GFP('INST.NET ECOSYSTEM EXCHANGE CO2...','GFP_INEE ',GFP%INEE )
  CALL PRINT_GFP('INST.GROSS PRIMARY PRODUCTION CO2.','GFP_IGPP ',GFP%IGPP )
  CALL PRINT_GFP('INST.ECOSYSTEM RESPIRATION CO2....','GFP_IREC ',GFP%IREC )
  CALL PRINT_GFP('INST.WETLAND CH4 FLUX.............','GFP_ICH4 ',GFP%ICH4 )
  DO JVAR=1,JPOSGHG
    WRITE(CLTEXU,FMT='(A24,I2.2,A8)') 'Total column GHG var nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_TCGHG',GFP%TCGHG(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSGHG  
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSCHEM
    WRITE(CLTEXU,FMT='(A24,I2.2,A8)') 'Total column CHEM var nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_TCCHEM',GFP%TCCHEM(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSCHEM
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSCHEMFLX
    WRITE(CLTEXU,FMT='(A24,I2.2,A8)') 'Total flux CHEM var nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_CHEMFLXO',GFP%CHEMFLXO(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSCHEMFLX
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSCHEMFLX
    WRITE(CLTEXU,FMT='(A24,I2.2,A8)') 'Wet dep flux CHEM var nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_CHEMWDFLX',GFP%CHEMWDFLX(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSCHEMFLX
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSCHEMFLX
    WRITE(CLTEXU,FMT='(A24,I2.2,A8)') 'Dry dep flux CHEM var nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_CHEMDDFLX',GFP%CHEMDDFLX(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSCHEMFLX
      EXIT
    ENDIF    
  ENDDO
  DO JDIAG=1,JPOSAERODIAG
    DO JVAR=1,JPOSAERO
      WRITE(CLTEXU,FMT='(A10,A3,A8,I2.2,A11)') 'AERO diag ', &
        & CPAERODIAG_LABEL(JDIAG), ' var nr ',JVAR,'........'
        CALL PRINT_GFP(CLTEXU,'GFP_AERODIAG',GFP%AERODIAG(JVAR,JDIAG) )
        IF (JVAR*JDIAG > ILIMIT*ILIMIT) THEN
          WRITE(UNIT=NULOUT, &
           & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSAERO*JPOSAERODIAG
          EXIT
        ENDIF    
    ENDDO
  ENDDO
  DO JDIAG=1,JPOSAERO_WVL_DIAG_TYPES
    DO JVAR=1,JPOSAERO_WVL_DIAG
      WRITE(CLTEXU,FMT='(A10,A3,A8,I2.2,A11)') 'AERO diag ', &
        & CPAERO_WVL_DIAG_LABEL(JDIAG), ' wvl nr ',JVAR,'........'
        CALL PRINT_GFP(CLTEXU,'GFP_AERO_WVL_DIAG',GFP%AERO_WVL_DIAG(JVAR,JDIAG) )
        IF (JVAR*JDIAG > ILIMIT*ILIMIT) THEN
          WRITE(UNIT=NULOUT, &
           & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSAERO_WVL_DIAG*JPOSAERO_WVL_DIAG_TYPES
          EXIT
        ENDIF    
    ENDDO
  ENDDO
  DO JVAR=1,JPOSEMIS2D
    WRITE(CLTEXU,FMT='(A24,I3.3,A8)') '2D emission flux nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_EMIS2D',GFP%EMIS2D(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSEMIS2D
      EXIT
    ENDIF    
  ENDDO
  DO JVAR=1,JPOSEMIS2DAUX
    WRITE(CLTEXU,FMT='(A24,I3.3,A8)') '2D emission aux fld nr ',JVAR,'........'
    CALL PRINT_GFP(CLTEXU,'GFP_EMIS2DAUX',GFP%EMIS2DAUX(JVAR) )
    IF (JVAR > ILIMIT) THEN
      WRITE(UNIT=NULOUT, &
       & FMT='(29X,''(truncated list - '',I3,'' variables)'')') JPOSEMIS2DAUX
      EXIT
    ENDIF    
  ENDDO
  CALL PRINT_GFP('Eastward water vapour flux........','GFP_VIWVE',GFP%VIWVE)
  CALL PRINT_GFP('Northward water vapour flux.......','GFP_VIWVN',GFP%VIWVN)
  CALL PRINT_GFP('TOT.COL.WATER   ..................','GFP_TCW  ',GFP%TCW  )
  CALL PRINT_GFP('TOT.COL.WATVAPOU..................','GFP_TCWV ',GFP%TCWV )
  CALL PRINT_GFP('TOT.COL.LIQ.WAT ..................','GFP_TCLW ',GFP%TCLW )
  CALL PRINT_GFP('TOT.COL.ICE.WAT ..................','GFP_TCIW ',GFP%TCIW )
  CALL PRINT_GFP('TOT.COL.RAIN.WAT .................','GFP_TCRW ',GFP%TCRW )
  CALL PRINT_GFP('TOT.COL.SNOW.WAT .................','GFP_TCSW ',GFP%TCSW )
  CALL PRINT_GFP('TOT.COL.SUPERCOOLED.LIQ.WAT ......','GFP_TCSLW',GFP%TCSLW )
  CALL PRINT_GFP('TOTAL.PRECIP.   ..................','GFP_TP   ',GFP%TP   )
  CALL PRINT_GFP('TEMP.SNOWLAYER  ..................','GFP_TSN  ',GFP%TSN  )
  CALL PRINT_GFP('T.SOLAR.RAD.    ..................','GFP_TSR  ',GFP%TSR  )
  CALL PRINT_GFP('TOP SOLAR  CLEAR..................','GFP_TSRC ',GFP%TSRC )
  CALL PRINT_GFP('T.THERMAL.RAD.  ..................','GFP_TTR  ',GFP%TTR  )
  CALL PRINT_GFP('TOP THERMALCLEAR..................','GFP_TTRC ',GFP%TTRC )
  CALL PRINT_GFP('HIG VEG. TYPE   ..................','GFP_TVH  ',GFP%TVH  )
  CALL PRINT_GFP('LOW VEG. TYPE   ..................','GFP_TVL  ',GFP%TVL  )
  CALL PRINT_GFP('HORIZ.VISIBILITY..................','GFP_VISIH',GFP%VISIH)
  CALL PRINT_GFP('ZERO DEG.LEVEL  ..................','GFP_0DEGL',GFP%EC0DEGL)
  CALL PRINT_GFP('-10 DEG.LEVEL  ...................','GFP_0DEGL',GFP%ECM10DEGL)
  CALL PRINT_GFP('SURF.ROUGHNESS  ..................','GFP_Z0F  ',GFP%Z0F  )
  CALL PRINT_GFP('CONV.AVAIL.POTENTIAL ENERGY.......','GFP_CAPE ',GFP%CAPE )
  CALL PRINT_GFP('MAX UNSTABLE.CAPE.................','GFP_CAPE ',GFP%MUCAPE )
  CALL PRINT_GFP('MIX.LAYER.50HPA  CAPE.............','GFP_CAPE ',GFP%MLCAPE50)
  CALL PRINT_GFP('MIX.LAYER.100HPA CAPE.............','GFP_CAPE ',GFP%MLCAPE100)
  CALL PRINT_GFP('DEPARTLEVEL.PRESSURE.MUCAPE.......','GFP_CAPE ',GFP%PDEPL)
  CALL PRINT_GFP('CAPE-SHEAR      ..................','GFP_CAPES',GFP%CAPES)
  CALL PRINT_GFP('CONV. INHIBITION .................','GFP_CIN  ',GFP%CIN  )
  CALL PRINT_GFP('CONV. INHIBITION 50HPA MIXLAYER...','GFP_CIN  ',GFP%MLCIN50  )
  CALL PRINT_GFP('CONV. INHIBITION 100HPA MIXLAYER..','GFP_CIN  ',GFP%MLCIN100  )
  CALL PRINT_GFP('CONV. K-INDEX    .................','GFP_KINDX',GFP%KINDEX )
  CALL PRINT_GFP('CONV. TT-INDEX   .................','GFP_TTIND',GFP%TTINDEX )
  CALL PRINT_GFP('CEILING CLOUD BASE................','GFP_CBASA',GFP%CBASEA)
  CALL PRINT_GFP('CLOUD TOP LEVEL CONVECTION........','GFP_CTOPC',GFP%CTOPC)
  CALL PRINT_GFP('MAX CAPE LAST 6H .................','GFP_CAPE6',GFP%MXCAP6)
  CALL PRINT_GFP('MAX CAPE-SHEAR LAST 6H ...........','GFP_CAPS6',GFP%MXCAPS6)
  CALL PRINT_GFP('PRESSURE THERMAL TROPOPAUSE ......','GFP_TROPT',GFP%TROPOTP)
  CALL PRINT_GFP('ZERO DEG.LEVEL WETBULB-T..........','GFP_ZTWB0',GFP%ZTWETB0)
  CALL PRINT_GFP('ONE DEG.LEVEL WETBULB-T. .........','GFP_ZTWB1',GFP%ZTWETB1)

  WRITE(UNIT=NULOUT,FMT='(/,'' CUMULATED FLUXES :'')')

  CALL PRINT_GFP('Large Scale Precipitation.........','GFP_CLSP ',GFP%CLSP )
  CALL PRINT_GFP('Convective precipitation..........','GFP_CCP  ',GFP%CCP  )
  CALL PRINT_GFP('Large Scale Snow fall.............','GFP_CLSS ',GFP%CLSS )
  CALL PRINT_GFP('Convective Snow Fall..............','GFP_CCSF ',GFP%CCSF )
  CALL PRINT_GFP('Large Scale Graupel fall .........','GFP_CLSG ',GFP%CLSG )
  CALL PRINT_GFP('Convective Graupel Fall...........','GFP_CCSG ',GFP%CCSG )
  CALL PRINT_GFP('Large Scale Hail fall.............','GFP_CLSH ',GFP%CLSH )
  CALL PRINT_GFP('Convective Hail Fall..............','GFP_CCSH ',GFP%CCSH )
  CALL PRINT_GFP('U-stress..........................','GFP_CUSS ',GFP%CUSS )
  CALL PRINT_GFP('V-stress..........................','GFP_CVSS ',GFP%CVSS )
  CALL PRINT_GFP('Surface Sensible Heat Flux........','GFP_CSSH ',GFP%CSSH )
  CALL PRINT_GFP('Surface Latent Heat Flux..........','GFP_CSLH ',GFP%CSLH )
  CALL PRINT_GFP('Tendency of Surface pressure......','GFP_CTSP ',GFP%CTSP )
  CALL PRINT_GFP('Total Cloud cover.................','GFP_CTCC ',GFP%CTCC )
  CALL PRINT_GFP('Boundary Layer Dissipation........','GFP_CBLD ',GFP%CBLD )
  CALL PRINT_GFP('Surface solar radiation...........','GFP_CSSR ',GFP%CSSR )
  CALL PRINT_GFP('Surface Thermal radiation.........','GFP_CSTR ',GFP%CSTR )
  CALL PRINT_GFP('Top Solar radiation...............','GFP_CTSR ',GFP%CTSR )
  CALL PRINT_GFP('Top Thermal radiation.............','GFP_CTTR ',GFP%CTTR )
  CALL PRINT_GFP('Convective Cloud Cover............','GFP_CCCC ',GFP%CCCC )
  CALL PRINT_GFP('High Cloud Cover..................','GFP_CHCC ',GFP%CHCC )
  CALL PRINT_GFP('Medium Cloud Cover................','GFP_CMCC ',GFP%CMCC )
  CALL PRINT_GFP('Low Cloud Cover...................','GFP_CLCC ',GFP%CLCC )
  CALL PRINT_GFP('U-Gravity-Wave Stress.............','GFP_CUGW ',GFP%CUGW )
  CALL PRINT_GFP('V-Gravity-Wave Stress.............','GFP_CVGW ',GFP%CVGW )
  CALL PRINT_GFP('U-Total Stress....................','GFP_CUTO ',GFP%CUTO )
  CALL PRINT_GFP('V-Total Stress....................','GFP_CVTO ',GFP%CVTO )
  CALL PRINT_GFP('Water Evaporation.................','GFP_CE   ',GFP%CE   )
  CALL PRINT_GFP('Snow Sublimation..................','GFP_CS   ',GFP%CS   )
  CALL PRINT_GFP('Water + Snow Evaporation..........','GFP_CT   ',GFP%CT   )
  CALL PRINT_GFP('Latent Heat Evaporation...........','GFP_CLHE ',GFP%CLHE )
  CALL PRINT_GFP('Latent Heat Sublimation...........','GFP_CLHS ',GFP%CLHS )
  CALL PRINT_GFP('Total Latent Heat (Water + Snow)..','GFP_CLHT ',GFP%CLHT )
  CALL PRINT_GFP('Soil Moisture.....................','GFP_CWS  ',GFP%CWS  )
  CALL PRINT_GFP('Snow mass.........................','GFP_CSNS ',GFP%CSNS )
  CALL PRINT_GFP('Total precipitable water..........','GFP_CQTO ',GFP%CQTO )
  CALL PRINT_GFP('Total Ozone.......................','GFP_CTO3 ',GFP%CTO3 )
  CALL PRINT_GFP('Top mesospheric enthalpy..........','GFP_CTME ',GFP%CTME )
  CALL PRINT_GFP('Solid specific moisture...........','GFP_CICE ',GFP%CICE )
  CALL PRINT_GFP('Liquid specific moisture..........','GFP_CLI  ',GFP%CLI  )
  CALL PRINT_GFP('Contrib. of Convection to U.......','GFP_CCVU ',GFP%CCVU )
  CALL PRINT_GFP('Contrib. of Convection to V.......','GFP_CCVV ',GFP%CCVV )
  CALL PRINT_GFP('Contrib. of Convection to Q.......','GFP_CCVQ ',GFP%CCVQ )
  CALL PRINT_GFP('Contrib. of Convection to Cp.T....','GFP_CCVS ',GFP%CCVS )
  CALL PRINT_GFP('Contrib. of Turbulence to Q.......','GFP_CTUQ ',GFP%CTUQ )
  CALL PRINT_GFP('Contrib. of Turbulence to Cp.T....','GFP_CTUS ',GFP%CTUS )
  CALL PRINT_GFP('Clear sky shortwave radiation.....','GFP_CSOC ',GFP%CSOC )
  CALL PRINT_GFP('Clear sky longwave radiation......','GFP_CTHC ',GFP%CTHC )
  CALL PRINT_GFP('Surface direct normal iradiance...','GFP_CDNI ',GFP%CDNI )
  CALL PRINT_GFP('Surface global normal iradiance...','GFP_CGNI ',GFP%CGNI )
  CALL PRINT_GFP('Surface parallel solar flux.......','GFP_CSOP ',GFP%CSOP )
  CALL PRINT_GFP('Top parallel solar flux...........','GFP_CTOP ',GFP%CTOP )
  CALL PRINT_GFP('Surface down solar flux...........','GFP_CSOD ',GFP%CSOD )
  CALL PRINT_GFP('Surface down thermic flux.........','GFP_CTHD ',GFP%CTHD )
  CALL PRINT_GFP('melt snow.........................','GFP_CFON ',GFP%CFON )
  CALL PRINT_GFP('flux de chaleur dans le sol.......','GFP_CCHS ',GFP%CCHS )
  CALL PRINT_GFP('flux d eau dans le sol............','GFP_CEAS ',GFP%CEAS )
  CALL PRINT_GFP('Surf. water content run-off.......','GFP_CSRU ',GFP%CSRU )
  CALL PRINT_GFP('Deep soil water content run-off...','GFP_CDRU ',GFP%CDRU )
  CALL PRINT_GFP('Interception water content run-off','GFP_CIRU ',GFP%CIRU )
  CALL PRINT_GFP('Evapotranspiration................','GFP_CETP ',GFP%CETP )
  CALL PRINT_GFP('Transpiration.....................','GFP_CTP  ',GFP%CTP  )
  CALL PRINT_GFP('Surface moon radiation............','GFP_CSMR ',GFP%CSMR )
  CALL PRINT_GFP('Top clear sky shortwave radiation.','GFP_CTSOC',GFP%CTSOC)
  CALL PRINT_GFP('Top clear sky longwave radiation..','GFP_CTTHC',GFP%CTTHC)
  CALL PRINT_GFP('Duration of total precipitations..','GFP_CDUTP',GFP%CDUTP)
  CALL PRINT_GFP('Surface PARadiation...............','GFP_SPAR ',GFP%SPAR )
  CALL PRINT_GFP('Surface UV-B radiation............','GFP_SUVB ',GFP%SUVB )
  CALL PRINT_GFP('Surface clear-sky PARadiation.....','GFP_SPARC',GFP%SPARC)
  CALL PRINT_GFP('TOA incident solar radiation......','GFP_STINC',GFP%STINC)
  CALL PRINT_GFP('Surf.total sky direct SW radiation','GFP_SFDIR',GFP%SFDIR)
  CALL PRINT_GFP('Surf.clear-sky direct SW radiation','GFP_SCDIR',GFP%SCDIR)
  CALL PRINT_GFP('Cumulative lightning density......','GFP_CFLASH',GFP%CFLASH)

  WRITE(UNIT=NULOUT,FMT='(/,'' INSTANTANEOUS DIAGNOSTICS :'')')

  CALL PRINT_GFP('Total Cloud cover.................','GFP_XTCC ',GFP%XTCC )
  CALL PRINT_GFP('U-component of wind at PBL height.','GFP_X10U ',GFP%X10U )
  CALL PRINT_GFP('V-component of wind at PBL height.','GFP_X10V ',GFP%X10V )
  CALL PRINT_GFP('U-neutral wind at PBL height......','GFP_X10NU',GFP%X10NU)
  CALL PRINT_GFP('V-neutral wind at PBL height......','GFP_X10NV',GFP%X10NV)
  CALL PRINT_GFP('Temperature at PBL height.........','GFP_X2T  ',GFP%X2T  )
  CALL PRINT_GFP('Wet bulb temperature at 2meters...','GFP_X2TPW',GFP%X2TPW)
  CALL PRINT_GFP('Specific Humidity at PBL height...','GFP_X2SH ',GFP%X2SH )
  CALL PRINT_GFP('Relative Humidity at PBL height...','GFP_X2RH ',GFP%X2RH )
  CALL PRINT_GFP('Convective Cloud Cover............','GFP_XCCC ',GFP%XCCC )
  CALL PRINT_GFP('High Cloud Cover..................','GFP_XHCC ',GFP%XHCC )
  CALL PRINT_GFP('Medium Cloud Cover................','GFP_XMCC ',GFP%XMCC )
  CALL PRINT_GFP('Low Cloud Cover...................','GFP_XLCC ',GFP%XLCC )
  CALL PRINT_GFP('Maximum temperature at PBL height.','GFP_XX2T ',GFP%XX2T )
  CALL PRINT_GFP('Minimum temperature at PBL height.','GFP_XN2T ',GFP%XN2T )
  CALL PRINT_GFP('Contribution of Convection to U...','GFP_XCVU ',GFP%XCVU )
  CALL PRINT_GFP('Contribution of Convection to V...','GFP_XCVV ',GFP%XCVV )
  CALL PRINT_GFP('Contribution of Convection to Q...','GFP_XCVQ ',GFP%XCVQ )
  CALL PRINT_GFP('Contribution of Convection to Cp.T','GFP_XCVS ',GFP%XCVS )
  CALL PRINT_GFP('Contribution of Turbulence to U...','GFP_XTUU ',GFP%XTUU )
  CALL PRINT_GFP('Contribution of Turbulence to V...','GFP_XTUV ',GFP%XTUV )
  CALL PRINT_GFP('Contribution of Turbulence to Q...','GFP_XTUQ ',GFP%XTUQ )
  CALL PRINT_GFP('Contribution of Turbulence to Cp.T','GFP_XTUS ',GFP%XTUS )
  CALL PRINT_GFP('Contribution of GWD to U..........','GFP_XGDU ',GFP%XGDU )
  CALL PRINT_GFP('Contribution of GWD to V..........','GFP_XGDV ',GFP%XGDV )
  CALL PRINT_GFP('Large Scale Precipitation.........','GFP_XLSP ',GFP%XLSP )
  CALL PRINT_GFP('Convective precipitation..........','GFP_XCP  ',GFP%CCP  )
  CALL PRINT_GFP('Large Scale Snow fall.............','GFP_XLSS ',GFP%XLSS )
  CALL PRINT_GFP('Convective Snow Fall..............','GFP_XCSF ',GFP%XCSF )
  CALL PRINT_GFP('Large Scale Graupel fall..........','GFP_XLSG ',GFP%XLSG )
  CALL PRINT_GFP('Convective Graupel Fall...........','GFP_XCSG ',GFP%XCSG )
  CALL PRINT_GFP('Large Scale Hail fall.............','GFP_XLSH ',GFP%XLSH )
  CALL PRINT_GFP('Convective Hail Fall..............','GFP_XCSH ',GFP%XCSH )
  CALL PRINT_GFP('Surface solar radiation...........','GFP_XSSR ',GFP%XSSR )
  CALL PRINT_GFP('Surface Thermal radiation.........','GFP_XSTR ',GFP%XSTR )
  CALL PRINT_GFP('Top Solar radiation...............','GFP_XTSR ',GFP%XTSR )
  CALL PRINT_GFP('Top Thermal radiation.............','GFP_XTTR ',GFP%XTTR )
  CALL PRINT_GFP('wind intensity at PBL height......','GFP_X10FF',GFP%X10FF)
  CALL PRINT_GFP('CAPE from the model...............','GFP_XCAPE',GFP%XCAPE)
  CALL PRINT_GFP('Pressure of top of convection.....','GFP_XCTOP',GFP%XCTOP)
  CALL PRINT_GFP('MOCON from the model..............','GFP_XMOCO',GFP%XMOCO)
  CALL PRINT_GFP('Height of the PBL.................','GFP_XCLPH',GFP%XCLPH)
  CALL PRINT_GFP('Ventilation Index.................','GFP_XVEIN',GFP%XVEIN)
  CALL PRINT_GFP('Gusts from of the model...........','GFP_XGUST',GFP%XGUST)
  CALL PRINT_GFP('U-momentum of gusts from the model','GFP_XUGST',GFP%XUGST)
  CALL PRINT_GFP('V-momentum of gusts from the model','GFP_XVGST',GFP%XVGST)
  CALL PRINT_GFP('Max. rel. moisture at PBL height..','GFP_XX2HU',GFP%XX2HU)
  CALL PRINT_GFP('Min. rel. moisture at PBL height..','GFP_XN2HU',GFP%XN2HU)
  CALL PRINT_GFP('U-momentum of gusts2 from model...','GFP_XUGST2',GFP%XUGST2)
  CALL PRINT_GFP('V-momentum of gusts2 from model...','GFP_XVGST2',GFP%XVGST2)
  CALL PRINT_GFP('Increment to mini temperature.....','GFP_INCTN',GFP%INCTN)
  CALL PRINT_GFP('Increment to maxi temperature.....','GFP_INCTX',GFP%INCTX)
  CALL PRINT_GFP('Increment to mini rel. moisture...','GFP_INCHN',GFP%INCHN)
  CALL PRINT_GFP('Increment to maxi rel. moisture...','GFP_INCHX',GFP%INCHX)
  CALL PRINT_GFP('Thetaprimwprim surface flux.......','GFP_XTHW',GFP%XTHW)
  CALL PRINT_GFP('Hail Diagnostic (AROME)...........','GFP_XXDIAGH',GFP%XXDIAGH)
  CALL PRINT_GFP('Dummy hail field (AROME)..........','GFP_ACCGREL',GFP%ACCGREL)
  CALL PRINT_GFP('Mean radiant temperature..........','GFP_XMRT',GFP%XMRT)
  CALL PRINT_GFP('Min. Refractivity Gradient in TPL.','GFP_DNDZN',GFP%DNDZN)
  CALL PRINT_GFP('Mean Refractivity Gradient in TPL.','GFP_DNDZA',GFP%DNDZA)
  CALL PRINT_GFP('Duct Base Height..................','GFP_DCTB ',GFP%DCTB )
  CALL PRINT_GFP('Trapping Layer Base Height........','GFP_TPLB ',GFP%TPLB )
  CALL PRINT_GFP('Trapping Layer Top Height.........','GFP_TPLT ',GFP%TPLT )
  CALL PRINT_GFP('Inst. total lightning density.....','GFP_LITOTI',GFP%LITOTI) 
  CALL PRINT_GFP('1h avg. total lightning density...','GFP_LITOTA1',GFP%LITOTA1)
  CALL PRINT_GFP('3h avg. total lightning density...','GFP_LITOTA3',GFP%LITOTA3)
  CALL PRINT_GFP('6h avg. total lightning density...','GFP_LITOTA6',GFP%LITOTA6)
  CALL PRINT_GFP('Inst. CTG lightning density.......','GFP_LICGI',GFP%LICGI)
  CALL PRINT_GFP('1h avg. CTG lightning density.....','GFP_LICGA1',GFP%LICGA1)
  CALL PRINT_GFP('3h avg. CTG lightning density.....','GFP_LICGA3',GFP%LICGA3)
  CALL PRINT_GFP('6h avg. CTG lightning density.....','GFP_LICGA6',GFP%LICGA6)
  CALL PRINT_GFP('1h precip type most-frequent(mode)','GFP_PTYPEMODE1',GFP%PTYPEMODE1)
  CALL PRINT_GFP('3h precip type most-frequent(mode)','GFP_PTYPEMODE3',GFP%PTYPEMODE3)
  CALL PRINT_GFP('6h precip type most-frequent(mode)','GFP_PTYPEMODE6',GFP%PTYPEMODE6)
  CALL PRINT_GFP('1h precip type most-severe........','GFP_PTYPESEVR1',GFP%PTYPESEVR1)
  CALL PRINT_GFP('3h precip type most-severe........','GFP_PTYPESEVR3',GFP%PTYPESEVR3)
  CALL PRINT_GFP('6h precip type most-severe........','GFP_PTYPESEVR6',GFP%PTYPESEVR6)
  CALL PRINT_GFP('Surf.tot sky direct beam SW rad...','GFP_SDSRP',GFP%SDSRP)
  CALL PRINT_GFP('Cloudy brightness temperature.....','GFP_CLBT ',GFP%CLBT )
  CALL PRINT_GFP('Clear-sky brightness temperature..','GFP_CSBT ',GFP%CSBT )
  CALL PRINT_GFP('Visibility due to water/ice cloud.','GFP_VISICLD',GFP%VISICLD)
  CALL PRINT_GFP('Visibility due to Precipitations..','GFP_VISIHYD',GFP%VISIHYD)
  CALL PRINT_GFP('Maximum of CLWC ..................','GFP_MXCLWC',GFP%MXCLWC)
  CALL PRINT_GFP('Visibility due to water/ice cloud.','GFP_VISICLD2',GFP%VISICLD2)
  CALL PRINT_GFP('Visibility due to Precipitations..','GFP_VISIHYD2',GFP%VISIHYD2)
  CALL PRINT_GFP('Maximum of CLWC ..................','GFP_MXCLWC2',GFP%MXCLWC2)
  CALL PRINT_GFP('Precipitation Type................','GFP_XPTYPE',GFP%XPTYPE)
  CALL PRINT_GFP('Severe Precipitation Type.........','GFP_XPTYPESEV',GFP%XPTYPESEV)
  CALL PRINT_GFP('Precipitation Type................','GFP_XPTYPE2',GFP%XPTYPE2)
  CALL PRINT_GFP('Severe Precipitation Type.........','GFP_XPTYPESEV2',GFP%XPTYPESEV2)

  WRITE(UNIT=NULOUT,FMT='(/,'' OCEAN MODEL DIAGNOSTICS :'')')

  CALL PRINT_GFP('Sea ice thickness.................','GFP_ICTH ',GFP%ICTH )
  CALL PRINT_GFP('Ocean sea surface height..........','GFP_SSH  ',GFP%SSH  )
  CALL PRINT_GFP('Ocean depth of 20degC isotherm....','GFP_EC20D',GFP%EC20D)
  CALL PRINT_GFP('Ocean mixed layer depth-..........','GFP_MLD  ',GFP%MLD  )
  CALL PRINT_GFP('Ocean sea surface salinity........','GFP_SSS  ',GFP%SSS  )
  CALL PRINT_GFP('Ocean temperature avg 300m........','GFP_TEM3 ',GFP%TEM3 )
  CALL PRINT_GFP('Ocean salnity avg 300m............','GFP_SAL3 ',GFP%SAL3 )

  DO J=1,NSFXPRE_CNT
    CALL PRINT_GFP ('SURFEX Field......................','GFP_SFXPRE        ',GFP%SFXPRE (J))
  ENDDO

ENDIF

!*       2. Control duplicated names/gribcodes
!           ----------------------------------

IERR=0
IF (LARPEGEF) THEN
! We should control that all the FA names are differents
  ISIZE=SIZE(TFP_DYNDS)
  DO J=1,ISIZE
    DO JVAR=J+1,ISIZE
      IF (TFP_DYNDS(J)%CLNAME /= ' ') THEN
        IF (TFP_DYNDS(J)%CLNAME==TFP_DYNDS(JVAR)%CLNAME) THEN
          IERR=IERR+1
          WRITE(NULOUT,*) 'SAME *FA* NAME DETECTED FOR TFP_*%ICOD = ',J, &
           & ' AND ',JVAR, ' : ',TFP_DYNDS(J)%CLNAME  
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  ISIZE=SIZE(GFP_PHYDS)

  DO J=1,ISIZE
    DO JVAR=J+1,ISIZE
      IF (GFP_PHYDS(J)%CLNAME /= ' ') THEN
        IF (GFP_PHYDS(J)%CLNAME==GFP_PHYDS(JVAR)%CLNAME) THEN
          IERR=IERR+1
          WRITE(NULOUT,*) 'SAME *FA* NAME DETECTED FOR GFP_*%ICOD = ',J, &
           & ' AND ',JVAR, ' : ',GFP_PHYDS(J)%CLNAME  
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ELSE
! We should control that all the grib codes (if valid) are different for upper
! air field on one side, and surface fields on the other side.
  ISIZE=SIZE(TFP_DYNDS)
  DO J=1,ISIZE
    DO JVAR=J+1,ISIZE
      IF (TFP_DYNDS(J)%IGRIB > 0) THEN
        IF (TFP_DYNDS(J)%IGRIB==TFP_DYNDS(JVAR)%IGRIB) THEN
          SELECT CASE (TFP_DYNDS(J)%LLSRF)
          CASE (.TRUE.)
            IF (TFP_DYNDS(JVAR)%LLSRF) THEN
              IERR=IERR+1
              WRITE(NULOUT,*) 'SAME GRIB CODE DETECTED FOR TFP_*%ICOD = ', &
               & J, ' AND ',JVAR, ' : ',TFP_DYNDS(J)%IGRIB  
            ENDIF 
          CASE (.FALSE.)
            IF (.NOT.TFP_DYNDS(JVAR)%LLSRF) THEN
              IERR=IERR+1
              WRITE(NULOUT,*) 'SAME GRIB CODE DETECTED FOR TFP_*%ICOD = ', &
               & J, ' AND ',JVAR, ' : ',TFP_DYNDS(J)%IGRIB  
            ENDIF
          END SELECT
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  ISIZE=SIZE(GFP_PHYDS)
  DO J=1,ISIZE
    DO JVAR=J+1,ISIZE
      IF (GFP_PHYDS(J)%IGRIB > 0) THEN
        IF (GFP_PHYDS(J)%IGRIB==GFP_PHYDS(JVAR)%IGRIB) THEN
          SELECT CASE (GFP_PHYDS(J)%LLSRF)
          CASE (.TRUE.)
            IF (GFP_PHYDS(JVAR)%LLSRF) THEN
              IERR=IERR+1
              WRITE(NULOUT,*) 'SAME GRIB CODE DETECTED FOR GFP_*%ICOD = ', &
               & J, ' AND ',JVAR, ' : ',GFP_PHYDS(J)%IGRIB  
            ENDIF
          CASE (.FALSE.)
            IF (.NOT.GFP_PHYDS(JVAR)%LLSRF) THEN
              IERR=IERR+1
              WRITE(NULOUT,*) 'SAME GRIB CODE DETECTED FOR GFP_*%ICOD = ', &
               & J, ' AND ',JVAR, ' : ',GFP_PHYDS(J)%IGRIB  
            ENDIF
          END SELECT 
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDIF
IF (IERR /= 0) THEN
  CALL ABOR1('SUAFN3 : ABOR1 CALLED')
ENDIF

!------------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('SUAFN3',1,ZHOOK_HANDLE)

!------------------------------------------------------------------------------

CONTAINS

SUBROUTINE PRINT_TFP(CDTEXT,CDNAME,YD)

! This tiny subroutine enables any sequence of the type.

CHARACTER(LEN=*),   OPTIONAL, INTENT(IN) :: CDTEXT
CHARACTER(LEN=*),   OPTIONAL, INTENT(IN) :: CDNAME
TYPE(FULLPOS_TYPE), OPTIONAL, INTENT(IN) :: YD

CHARACTER(LEN=12) :: CLNAME
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SUAFN3:PRINT_TFP',0,ZHOOK_HANDLE)

IF (LARPEGEF) THEN
  IF (PRESENT(CDTEXT).AND.PRESENT(CDNAME).AND.PRESENT(YD)) THEN
    CLNAME=CDNAME//' '
    WRITE(NULOUT,FMT='(1X,A25,'' : '',A12,1X,A16,3X,I2,5X,I2,4X,I4,4X,A8,  &
     & L1,4X,I2,2X,F4.1,2X,A10)')                                          &
     & CDTEXT,CLNAME, YD%CLNAME, YD%IBITS, YD%INTER, YD%IORDR, &
     & YD%CLPAIR, YD%LLGP, YD%ISF, YD%ZFK, YD%CLNIL  
  ELSE
    WRITE(UNIT=NULOUT,FMT='(3X,''Description'',12X,'' : '','' TYPE NAME  '', &
     & 1X,''    %CLNAME     '',1X,''%IBITS'',1X,''%INTER'',1X,''%IORDR'',1X, &
     & ''%CLPAIR'',1X,''%LLGP'',1X,''%ISF'',1X,''%ZFK'',1X,''%CLNIL'')')  
  ENDIF
ELSE
  IF (PRESENT(CDTEXT).AND.PRESENT(CDNAME).AND.PRESENT(YD)) THEN
    IF (YD%IGRIB > 0) THEN
      CLNAME=CDNAME//' '
      WRITE(NULOUT,FMT='(1X,A25,'' : '',A12,1X,I6,3X,I2,5X,I2,4X,I4,4X,A8,  &
       & L1,4X,I2,2X,F4.1,2X,A10)')                                         &
       & CDTEXT,CLNAME, YD%IGRIB, YD%IBITS, YD%INTER, YD%IORDR, &
       & YD%CLPAIR, YD%LLGP, YD%ISF, YD%ZFK, YD%CLNIL  
    ENDIF
  ELSE
    WRITE(UNIT=NULOUT,FMT='(3X,''Description'',12X,'' : '','' TYPE NAME  '', &
     & 1X,''%IGRIB'',1X,''%IBITS'',1X,''%INTER'',1X,''%IORDR'',1X,           &
     & ''%CLPAIR'',1X,''%LLGP'',1X,''%ISF'',1X,''%ZFK'',1X,''%CLNIL'')')  
  ENDIF
ENDIF

IF (LHOOK) CALL DR_HOOK('SUAFN3:PRINT_TFP',1,ZHOOK_HANDLE)

END SUBROUTINE PRINT_TFP

END SUBROUTINE SUAFN3


