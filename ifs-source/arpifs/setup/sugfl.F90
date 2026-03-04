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

SUBROUTINE SUGFL(YDDIMV,YDMODEL,KGFLCONF)

!**** *SUGFL*  - Initialize definition of unified_treatment grid_point fields

!     Purpose.
!     --------
!           Initialize definition of unified_treatment fields (GFL)
!           The GFL ordering in this routine is the ordering for "t" and
!           "t+dt" values of the GFL in grid-point arrays (GFL and GFLT1),
!           and also the ordering in the spectral array SPGFL.

!**   Interface.
!     ----------
!        *CALL* *SUGFL

!        Explicit arguments :
!        --------------------
!           KGFLCONF : default configuration of GFL attributes

!        Implicit arguments :
!        --------------------
!        MODULE YOM_YGFL

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------
!        None

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud  *ECMWF*

!     Modifications.
!     --------------
!        Original : 2002-03-11
!        M.Hamrud : 2002-08-01 Extensive mods
!        R. El Khatib : 2003-08-19 ARPEGE field names
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        M. Tudor : 2003-10-17 introduce YCPF and YSPF in GFL
!        A. Untch : 2005-03-11 introduce aerosols
!        J. Flemming : 2005-04-11 replace aerosols with reactive gases
!        Y. Seity  : 2004-11-16 AROME GFL fields
!        C. Fischer: 23-May-2005 No Ozone in the Meteo-France control variable
!        Y. Bouteloup: 13-Oct-2005 introduce YCVGQ (Moisture convergence for French physics)
!        J. Haseler : 2005-10-11 introduce LTRAJIO
!        JJMorcrette  20060512 GEMS variables in cy31
!        S. Serrar    20060907 tracers used for diagnostics only 
!        A. Alias : 2006-10-13 introduction of YSDSAT and YCVV
!        B. Sass :  HIRLAM pseudo-prognostic field YQVA
!        M. Bellus : 27-Sep-2006 introduce YUOM, YUAL, YDOM, YDAL, YUEN and
!                    YUNEBH (ALARO-0 prognostic convection) + spotted/corrected bug
!                    in index incrementation after YO3
!        J.Haseler : 27-Feb-2007 Generalize definition of cloud fields
!        JJMorcrette  20070312 Aerosol variables in cy32R1
!        S. Serrar 2007-03-22 introduce ERA40 GFL fields (YERA40)
!        S. Serrar 20070717 new GEMS GFL fields
!        Y. Seity  2007-12-18 add GFL Hail (YH)  
!        JJMorcrette  20090217 PP of aerosol and UV processor outputs
!        Y. Bouteloup 2009-10-15 add Total cloud water and ice (YIRAD and YLRAD)
!        K. YESSAD (Nov 2011) gather all GFL set-up
!        R. El Khatib 18-Jan-2018 optional argument KGFLCONF
!     ------------------------------------------------------------------

USE TYPE_MODEL , ONLY : MODEL
USE YOMDIMV  , ONLY : TDIMV
USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(TDIMV), INTENT(IN) :: YDDIMV
TYPE(MODEL) ,INTENT(INOUT):: YDMODEL
INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: KGFLCONF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------

#include "sugfl1.intfb.h"
#include "sugfl2.intfb.h"
#include "sugfl3.intfb.h"

!-------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGFL',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------

! code taken from SUDIM1 (reads NAMGFL)
CALL SUGFL1(YDMODEL,KGFLCONF)

! former content of SUGFL
CALL SUGFL2(YDMODEL%YRML_PHY_EC%YREPHY,YDMODEL%YRML_GCONF, YDMODEL%YRML_CHEM%YRCOMPO)

! former content of SUDYN_SETGFLATTR
CALL SUGFL3(YDDIMV,YDMODEL)

!-------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUGFL',1,ZHOOK_HANDLE)
END SUBROUTINE SUGFL
