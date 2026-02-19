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

SUBROUTINE SUOPH0(CDEXP)

!**** *SUOPH0* - Routine to prepare input/output handling

!     Purpose.   To set up common block YOMOPH which contains file-
!     --------   handling parameters.

!**   Interface.
!     ----------
!        *CALL* *SUOPH*
!           CDEXP : name of the experiment

!        Explicit arguments :
!        --------------------
!           CDEXP : name of experiment (expected 4 letters long) - input

!        Implicit arguments :
!        --------------------

!     Method.
!     -------

!     Externals :
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      Mats Hamrud and Philippe Courtier  *ECMWF*
!      Original : 88-03-01

!     Modifications.
!     --------------
!      M. Janousek: 01-11-26 : NCADFORM - format of FA file frame
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      R.Buizza      23-May-2005 Added CFNVAREPS (input file with VAREPS accum fields)
!      G. Desroziers and K. Yessad (sept 2005):
!       - split option LTRAJHR into LTRAJHR_ALTI and LTRAJHR_SURF.
!       - adapt option LTRAJHR to METEO-FRANCE configurations.
!      JD Gril : 03-feb-2006 : default NCADFORM=1
!      R. El Khatib : 23-Oct-2008 : CFPATH now a directory ending with a '/'
!        and used also for CFNHWF + merge with sueoph + cleanings + CEFNLSH in namelist
!      E.Holm        13-Nov-2008 Added CFNBGHRSH/GG, background
!                                at outer loop resolution.
!      K. Yessad (Jan 2010): remove useless variables.
!      G.Mozdzynski (Feb 2011): OOPS cleaning, use of derived type TCSGLEG
!      P.Marguinaud : 26-Apr-2012 : Handle new parameter NTIMEFMT
!      P. Bechtold 14/05/2012 replace 86400 by INT(RDAY)
!      R. El Khatib 27-Sep-2013 Boyd window in frame
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      M.Hamrud  (Jan 2015) ! Remove FDB stuff from here
!      R. El Khatib : 03-Dec-2014 skeleton of the configuration 903
!      R. El Khatib 08-Dec-2015 CDEXP and KCONF as argument ; remove LECMWF (useless)
!      R. El Khatib 22-mar-2016  CFNCLIMIN, CFNCLIMOUT, CCLIMINC
!      R. El Khatib 13-Sep-2016 CFNIUA name of Upper air gp GRIB data for conf 903
!     ------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMOPH0  , ONLY : CFNSH, CFNGG, CFNUA, CFNISH, CFNIGG, CFNIUA,CFNPL, CFNVAREPS, &
 & CFNRF, CFNAN, CFNGR, CFNINSH, CFNINGG, CFNHWF, CFNGSH, CFNGGG, CFNDDH, CFPEXTSFX, &
 & CFPATH, CFNBS, CFNBGV, CFNCGL, CFNFGI, CFNANI, CFNFGIG, CFNANIG, LINC, NTIMEFMT, &
 & CFANS, NCADFORM, CEFLS, CEFNLSH, LBCINC, CFNCLIMIN, CFNCLIMOUT, CCLIMINC, &
 & CFNTRAJHRGRID, CFNTRAJHRSPEC, CFNTRAJHRSURF, CFNTRAJHR, &
 & CFNTRAJBGGRID, CFNTRAJBGSPEC, CFNTRAJBG, CFNBGHRSH, CFNBGHRGG, &
 & CFNRFBDP, CFNGSHBDP, CFNPANSH, CFNPANGG
USE YOMLUN   , ONLY : NULNAM, NULOUT

!     ------------------------------------------------------------------

IMPLICIT NONE

CHARACTER(LEN=*),   INTENT(IN) :: CDEXP

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

#include "namoph.nam.h"

#include "posnam.intfb.h"

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SUOPH0',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

!*       0.    READ NAMELIST
!              -------------

!        0.1 Set implicit default values

CFPATH='./'
LINC=.FALSE.
LBCINC=.FALSE.
NCADFORM=1
NTIMEFMT=0
CFNHWF='hiwif'
CEFNLSH='ELSCF'//CDEXP(1:4)//'ALBC'
CFNCLIMIN='Const.Clim'
CFNCLIMOUT='const.clim.'
CCLIMINC=' '
CFPEXTSFX='.sfx'


CALL POSNAM(NULNAM,'NAMOPH')
READ(NULNAM,NAMOPH)

WRITE(NULOUT,'('' LINC = '',L2,'' CFNHWF = '',A,'' LBCINC = '',L2)') LINC, TRIM(CFNHWF), LBCINC
WRITE(NULOUT,'('' NCADFORM = '',I2, '' NTIMEFMT = '',I2)') NCADFORM, NTIMEFMT
WRITE(NULOUT,'('' CEFNLSH = '',A)') TRIM(CEFNLSH)
WRITE(NULOUT,'('' CFNCLIMIN = '',A,'' CFNCLIMOUT = '',A, '' CCLIMINC = '',A8, '' CFPEXTSFX = '',A8)') &
 & TRIM(CFNCLIMIN), TRIM(CFNCLIMOUT), TRIM(CCLIMINC), TRIM(CFPEXTSFX)

!     ------------------------------------------------------------------

!*       1.    SET FILE NAME FORMATS
!              ---------------------

!*       1.1  WORK FILES AND DIAGNOSTIC OUTPUT

CFNSH=TRIM(CFPATH)//'ICMSH'//CDEXP(1:4)//'0000'//'   '
CFNGG='ICMGG'//CDEXP(1:4)//'0000'//'   '
CFNUA='ICMUA'//CDEXP(1:4)//'0000'//'   '
CFNRF  ='ICMRF'//CDEXP(1:4)//'0000'//'   '
CFNRFBDP = 'ICMRF'//CDEXP(1:4)//'0000'//'   '

!*       1.2  DIAGNOSTIC OUTPUT FILE NAMES
!*            PX: PRESSURE LEVEL, CONTROL VARIABLE
!*            MX: MODEL LEVEL, CONTROL VARIABLE
!*            MG: MODEL LEVEL, GRADIENT
!*            MD: MODEL LEVEL, DEPARTURES (for 151)
!*            MB: MODEL LEVEL, RANDOM BACKGROUND VECTOR
!*            ME: MODEL LEVEL, (SCALED) EIGENVECTOR OF THE HESSIAN
CFNPL ='PX'//CDEXP(1:4)//'000'//'+'//'0000'
CFNAN ='MX'//CDEXP(1:4)//'000'//'+'//'0000'
CFNGR ='MG'//CDEXP(1:4)//'000'//'+'//'0000'
CFNINSH ='MISH'//CDEXP(1:4)//'000'//'+'//'0000'
CFNINGG ='MIGG'//CDEXP(1:4)//'000'//'+'//'0000'
CFNBGV='MB'//CDEXP(1:4)//'000'
CFNCGL='ME'//CDEXP(1:4)//'000'

!*       1.4 HIGH RESOLUTION TRAJECTORY FILE NAME

CFNTRAJHRGRID='TRAJHR00/trajgrid'
CFNTRAJHRSPEC='TRAJHR00/trajspec'
CFNTRAJHRSURF='TRAJHR00/trajsurf000'
CFNTRAJHR='TRAJ'//CDEXP(1:4)//'0000    '

!*       1.5 BACKGROUND FILE NAME

CFNTRAJBGGRID='TRAJBG00/bckggrid'
CFNTRAJBGSPEC='TRAJBG00/bckgspec'

! ky: I am not yet sure of the name of this file, to be validated
!     and maybe fixed later according to the validations we expect to
!     do with a 4DVAR-LTRAJHR OLIVE script.
CFNTRAJBG=CFNRF

!*       1.6 BACKGROUND AT OUTER LOOP RESOLUTION
!            The high resolution background is needed at all times
!            where the increments are added.

CFNBGHRSH='ICMSH'//CDEXP(1:4)//'BGHR'//'   '
CFNBGHRGG='ICMGG'//CDEXP(1:4)//'BGHR'//'   '



!     ------------------------------------------------------------------

!*       2.    PREPARE INITIAL DATA UNITS.
!              ----------------------------

CFNISH='ICMSH'//CDEXP(1:4)//'INIT'//'   '
CFNIGG='ICMGG'//CDEXP(1:4)//'INIT'//'   '
CFNIUA='ICMUA'//CDEXP(1:4)//'INIT'//'   '
CFNGSH='ICMSH'//CDEXP(1:4)//'IMIN'//'   '
CFNGGG='ICMGG'//CDEXP(1:4)//'IMIN'//'   '
CFNDDH='DHFTT'//CDEXP(1:4)//'IMIN'//'   '
CFNPANSH='ICMSH'//CDEXP(1:4)//'PANA'//'   '
CFNPANGG='ICMGG'//CDEXP(1:4)//'PANA'//'   '
CFNBS ='ICMSH'//CDEXP(1:4)//'BIAS'//'   '
CFANS ='ICMSH'//CDEXP(1:4)//'SURF'//'   '

CFNGSHBDP='ICMSH'//CDEXP(1:4)//'IBDP'//'   '

! Initial condition file with accum fields for VAREPS
CFNVAREPS ='ICVEP'//CDEXP(1:4)//'INIT'//'   '

CFNFGI ='ICMSH'//CDEXP(1:4)//'FGIN'//'   '
CFNANI ='ICMSH'//CDEXP(1:4)//'ANIN'//'   '
CFNFGIG='ICMGG'//CDEXP(1:4)//'FGIN'//'   '
CFNANIG='ICMGG'//CDEXP(1:4)//'ANIN'//'   '

! Lateral Boundary coupling data units

CEFLS='ELSCF'//CDEXP(1:4)//'COMP'//'    '

! -------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('SUOPH0',1,ZHOOK_HANDLE)
END SUBROUTINE SUOPH0
