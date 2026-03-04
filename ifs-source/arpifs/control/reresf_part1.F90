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

SUBROUTINE RERESF_PART1(YDGEOMETRY,YDFIELDS,YDMODEL,KRSTEP)

!**** *RERESF_PART1*  - read  restart files

!     Purpose.
!     --------
!     Read  a set of restart files if they exist:
!       (1) spectral fields
!       (2) upper-air grid-point fields (t-dt)
!       (3) surface fields
!       (4) upper-air grid-point fields (t)
!       (5) cumulated fluxes
!       (6) instantaneous fluxes
!       (7) restart control file

!**   Interface.
!     ----------
!        *CALL* *RERESF_PART1

!        Explicit arguments :      None.
!        --------------------

!        Implicit arguments :      None.
!        --------------------

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!      David Dent *ECMWF*
!      Original : 92-04-06

!     Modifications.
!     --------------
!      Modified : 01-06-13 R. El Khatib (ISP restart)
!      R. El Khatib : 01-08-07 Pruning options
!      Modified : 02-06-10 by D.Dent - Additional file for GP fields
!                          LFASTRES option and improved DDH handling
!      Modified : 08-2002 C. Smith : use "w" as prognostic variable in the
!       semi-lag advection of vertical divergence in the NH model.
!      M.Hamrud      01-Oct-2003 CY28 Cleaning
!      M.Hamrud      10-Jan-2004 CY28R1 Cleaning
!      Modified : 02-07-17 by P. Marquet : add VCLIA for aerosol files
!      Modified : 04-08-12 R. El Khatib Bugfixes + COPY
!      JJMorcrette 20060721 PP clear-sky PAR and TOA incident solar radiation
!      M.Hamrud      01-Jul-2006  Revised surface fields
!      JJMorcrette 20070321 Prognostic aerosols in rad and clouds
!      M.Hamrud      13-Oct-2008 Use IOSTREAM, one file per task, only FASTRES supported
!      JJMorcrette 20091201 Total and clear-sky direct SW radiation flux at surface 
!      M. Steinheimer Oct-2009  improved for stochastic physics options LSTOPH and LSTOPH_SPBS
!      M. Ahlgrimm 31 Oct 2011 Surface downward clear-sky LW and SW fluxes
!      K. Mogensen: 26-Jan-2012 Ocean fields when coupling to LIM2
!      T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
!      K. Yessad (July 2014): Move some variables.
!      F. Vana  05-Mar-2015  Support for single precision
!      F. Vana & M. Kharoutdinov 06-Feb-2015: Super-parametrization fields restart.
!      R. El Khatib 07-Mar-2016 Pruning of ISP
!      SJ Lock:    Jan-2016  Removed LSTOPH option
!      O. Marsden  Aug 2016  Removed use of SPA3, replaced by YDFIELDS%YRSPEC
!     ------------------------------------------------------------------

USE TYPE_MODEL         , ONLY : MODEL
USE GEOMETRY_MOD       , ONLY : GEOMETRY
USE FIELDS_MOD         , ONLY : FIELDS
USE PARKIND1           , ONLY : JPIM, JPRB
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMLUN             , ONLY : NULOUT, NULRCF, NULERR
USE YOMCT0             , ONLY : LNF
USE YOMCT2             , ONLY : NSTAR2
USE YOMIOS             , ONLY : CFRCF
USE YOMRES             , ONLY : CTIME, CSTEP
USE YOMGRIB            , ONLY : NSTEPLPP
USE YOMMP0             , ONLY : MYPROC, NPRTRW, NPRTRV, NPRGPNS, NPRGPEW, MYSETV

!     ------------------------------------------------------------------

IMPLICIT NONE

TYPE(GEOMETRY),                INTENT (INOUT)        :: YDGEOMETRY
TYPE(FIELDS),                  INTENT (INOUT)        :: YDFIELDS
TYPE(MODEL),  TARGET,          INTENT (INOUT)        :: YDMODEL
INTEGER(KIND=JPIM), OPTIONAL,  INTENT (OUT)          :: KRSTEP

LOGICAL :: LLEXIST
REAL(KIND=JPRB), POINTER  :: GMASS0, GMASSI 

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE


!     ------------------------------------------------------------------


#include "namrcf.nam.h"

!     ------------------------------------------------------------------

#include "abor1.intfb.h"

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RERESF_PART1',0,ZHOOK_HANDLE)

!*       1.    READ  RESTART CONTROL FILE
!              --------------------------

INQUIRE(FILE=CFRCF,EXIST=LLEXIST)

IF (PRESENT (KRSTEP)) KRSTEP = 0

IF (LLEXIST) THEN
  GMASS0 => YDMODEL%YRML_DYN%YRDYN%GMASS0
  GMASSI => YDMODEL%YRML_DYN%YRDYN%GMASSI

  LNF = .FALSE.
  OPEN  (NULRCF,FILE=CFRCF,DELIM="QUOTE")
  READ  (NULRCF,NAMRCF)
  CLOSE (NULRCF,STATUS='KEEP')

!        1.1   PERFORM CROSS-CHECKS ON NAMELIST VARIABLES

  IF(NPRGPNS /= IPRGPNSRES.OR.&
     & NPRGPEW /= IPRGPEWRES.OR.&
     & NPRTRW /= IPRTRWRES.OR.&
     & NPRTRV /= IPRTRVRES&
     & ) THEN  
!     must restart with same parallel configuration
    WRITE(NULOUT,'(A/A,I4.4/A,I4.4/A,I4.4/A,I4.4/A,L1/)')&
     & ' restart with '//&
     & ' different parallel configuration',&
     & 'NPRGPNS=',IPRGPNSRES,&
     & 'NPRGPEW=',IPRGPEWRES,&
     & 'NPRTRW=',IPRTRWRES,&
     & 'NPRTRV=',IPRTRVRES
    
    CALL ABOR1('RERESF_PART1: PROBLEM IN PARALLEL RESTART')
  ENDIF

!        1.2   MODIFY CONTROL VARIABLES FROM CONTENTS OF RCF

  READ(CSTEP,'(I8)') NSTAR2
  IF (PRESENT (KRSTEP)) KRSTEP = NSTAR2

ELSE

!*       6.    RESTART CONTROL FILE DOES NOT EXIST
!              -----------------------------------

  LNF = .TRUE.

ENDIF

IF (LHOOK) CALL DR_HOOK('RERESF_PART1',1,ZHOOK_HANDLE)
END SUBROUTINE RERESF_PART1

