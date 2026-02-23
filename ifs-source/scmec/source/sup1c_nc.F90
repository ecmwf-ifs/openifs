! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

      SUBROUTINE SUP1C_NC(YDDIMV,YDMODEL,YDSURF)

#ifdef DOC

!**** *SUP1C_NC * - Routine to initialize prognostic NetCDF output

!     Purpose.
!     --------
!           Initialize NetCDF outputfile progvar.nc

!**   Interface.
!     ----------
!        *CALL* *SUP1C_NC

!     Explicit arguments :
!     --------------------
!        none

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------

!     Author.
!     -------
!        Bart vd HURK (KNMI)

!     Modifications.
!     --------------
!        Original     2000-7-14
!        M. Ko"hler      9-2000  adoption to single column model
!        M. Ko"hler    6-6-2006  Single Column Model integration within IFS 
!        G. Carver     Aug/2012  Bugfixes for gfortran compile.
!     ------------------------------------------------------------------
#endif

USE YOMDIMV  , ONLY : TDIMV
USE TYPE_MODEL, ONLY : MODEL
USE PARKIND1 , ONLY : JPIM
USE YOMLOG1C , ONLY : CMODID, CSIMID, NPOSPRG
USE YOMRIP0  , ONLY : NINDAT, NSSSSS
USE SURFACE_FIELDS_MIX , ONLY : TSURF

IMPLICIT NONE

TYPE(TDIMV), INTENT(INOUT) :: YDDIMV
TYPE(MODEL), INTENT(INOUT) :: YDMODEL
TYPE(TSURF), INTENT(INOUT) :: YDSURF
INTEGER(KIND=JPIM) :: ICSS

INTEGER(KIND=JPIM) :: INCID, I,&
     &  NLEVDID, NLEVP1DID, NLEVSDID, NTIMDID, NVARID,&
     &  IDIMID2(2), IDIMID3(2), IDIMID4(2), IDIMID5(2),&
     &  ILEV(YDDIMV%NFLEVG), ILEVP1(YDDIMV%NFLEVG+1),&
     &  ILEVS(YDSURF%YSP_SBD%NLEVS), ISTATUS, IACCUR

CHARACTER (LEN=40) :: TITLE

#include "netcdf.inc"

!     ------------------------------------------------------------------
#include "handle_err_nc.intfb.h"
#include "varsetup1c_nc.intfb.h"
!     ------------------------------------------------------------------
ASSOCIATE(YDEPHY=>YDMODEL%YRML_PHY_EC%YREPHY)

ICSS=YDSURF%YSP_SBD%NLEVS

!        1.    open NetCDF files.
!              ------------------

ISTATUS = NF_CREATE ('progvar.nc', NF_CLOBBER, INCID)
CALL HANDLE_ERR_NC(ISTATUS)
WRITE(*,*) 'NETCDF-FILE progvar.nc OPENED ON UNIT ',INCID
NPOSPRG = INCID        !store unit number for later use

!        output accuracy
 IACCUR=NF_FLOAT
!iaccur=NF_DOUBLE


!        2.    meta data set up.
!              -----------------

!        title
TITLE = 'SCM: ' // TRIM(CMODID) // '  Sim: ' // TRIM(CSIMID)
ISTATUS = NF_PUT_ATT_TEXT (INCID, NF_GLOBAL, 'title', LEN(TRIM(TITLE)), TRIM(TITLE))
CALL HANDLE_ERR_NC(ISTATUS)

!        model identification
ISTATUS = NF_PUT_ATT_TEXT (INCID, NF_GLOBAL, 'modelID', LEN(TRIM(CMODID)), TRIM(CMODID))
CALL HANDLE_ERR_NC(ISTATUS)

!        simulation identification
ISTATUS = NF_PUT_ATT_TEXT (INCID, NF_GLOBAL, 'simulationID', LEN(TRIM(CSIMID)), TRIM(CSIMID))
CALL HANDLE_ERR_NC(ISTATUS)

!        dataID for MetView
ISTATUS = NF_PUT_ATT_TEXT (INCID, NF_GLOBAL, 'dataID', 10, 'SCM_OUTPUT' ) 
CALL HANDLE_ERR_NC(ISTATUS)

!        model startup
ISTATUS = NF_PUT_ATT_INT (INCID, NF_GLOBAL,  'start_day',  NF_INT, 1, NINDAT)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_PUT_ATT_INT (INCID, NF_GLOBAL,  'start_hour', NF_INT, 1, NSSSSS)
CALL HANDLE_ERR_NC(ISTATUS)

!        create dimensions
ISTATUS = NF_DEF_DIM (INCID, 'nlev',   YDDIMV%NFLEVG,       NLEVDID)   !atmosphere
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_DEF_DIM (INCID, 'nlevp1', YDDIMV%NFLEVG+1,     NLEVP1DID) !atmosphere + 1
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_DEF_DIM (INCID, 'nlevs',  ICSS,         NLEVSDID)  !land/sea-ice
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_DEF_DIM (INCID, 'time',   NF_UNLIMITED, NTIMDID)   !time
CALL HANDLE_ERR_NC(ISTATUS)

!        2.1   create dimensional variables.

!        atmospheric model levels
CALL VARSETUP1C_NC (INCID, NF_INT, 1, (/ NLEVDID /), &
   'nlev',    'Atmospheric Model Levels',  'count')
!        atmospheric model levels plus 1
CALL VARSETUP1C_NC (INCID, NF_INT, 1, (/ NLEVP1DID /), &
   'nlevp1',  'Atmospheric Model Levels',  'count')
!        soil/sea-ice model levels
CALL VARSETUP1C_NC (INCID, NF_INT, 1, (/ NLEVSDID /), &
   'nlevs',   'Soil/Sea-Ice Model Levels', 'count')
!        time
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'time',    'Time',                      'seconds')

IDIMID2(1) = NLEVDID
IDIMID2(2) = NTIMDID
IDIMID3(1) = NLEVP1DID
IDIMID3(2) = NTIMDID
IDIMID4(1) = NLEVSDID
IDIMID4(2) = NTIMDID

!        2.2   create atmospheric variables.

CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'pressure_f',        'Pressure - full levels',         'Pa')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID3, &
   'pressure_h',        'Pressure - half levels',         'Pa')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'height_f',          'Height - full levels',           'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID3, &
   'height_h',          'Height - half levels',           'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'relative_humidity', 'Relative Humidity',              '0-1')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   't',                 'Temperature',                    'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'pot_temperature',   'Potential Temperature',          'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'pot_temp_e',        'Equivalent Potential Temperature', 'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'u',                 'U Wind',                         'm/s')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'v',                 'V Wind',                         'm/s')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'dry_st_energy',     'Dry Static Energy',              'J/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'moist_st_energy',   'Moist Static Energy',            'J/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'q',                 'Water Vapor Mixing Ratio',       'kg/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'ql',                'Liquid Water Mixing Ratio',      'kg/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'qi',                'Ice Water Mixing Ratio',         'kg/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'cloud_fraction',    'Cloud Fraction',                 '0-1')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'qr',                'Rain Water Mixing Ratio',        'kg/kg')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'qsn',               'Snow Water Mixing Ratio',        'kg/kg')


!        2.3   create land/ocean/sea-ice variables.

CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   't_skin',            'Skin Temperature',               'Kelvin')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'q_skin',            'Skin Reservoir Water Content',   'scaled m')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID4, &
   't_soil',            'Soil Temperature ',              'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID4, &
   'q_soil',            'Soil Moisture',                  'm^3/m^3')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'snow',              'Snow Depth (liquid equivalent)', 'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   't_snow',            'Snow Temperature',               'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'albedo_snow',       'Snow Albedo',                    '0-1')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'density_snow',      'Snow Density',                   'kg/m^3')


!        2.4   create lake_variables

IF (YDEPHY%LEFLAKE) THEN

CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'dl',                'Lake Depth',                 'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'lmlt',              'Lake Mix-layer Temperature', 'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'lmld',              'Lake Mix-layer Depth',       'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'lblt',              'Lake Bottom Temperature',    'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'ltlt',              'Lake Total Layer Temperature', 'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'lshf',              'Lake Shape Factor', '-')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'lict',              'Lake Ice Temperature',       'K')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'licd',              'Lake Ice Depth',             'm')
CALL VARSETUP1C_NC (INCID, IACCUR, 1, (/ NTIMDID /), &
   'cl',                'Lake Cover',                 '0-1')

ENDIF

!        2.5   create q_sat

CALL VARSETUP1C_NC (INCID, IACCUR, 2, IDIMID2, &
   'q_sat',             'Saturation Specific Humidity',   'kg/kg')

!        3.    write model level data
!              ----------------------

!        return to data mode
ISTATUS = NF_ENDDEF(INCID)
CALL HANDLE_ERR_NC(ISTATUS)

!        atmospheric model levels
DO I=1,YDDIMV%NFLEVG
  ILEV(I)=I
ENDDO
ISTATUS = NF_INQ_VARID   (INCID, 'nlev', NVARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_PUT_VAR_INT (INCID, NVARID, ILEV)
CALL HANDLE_ERR_NC(ISTATUS)

!        atmospheric model levels + 1
DO I=1,YDDIMV%NFLEVG+1
  ILEVP1(I)=I
ENDDO
ISTATUS = NF_INQ_VARID   (INCID, 'nlevp1', NVARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_PUT_VAR_INT (INCID, NVARID,   ILEVP1)
CALL HANDLE_ERR_NC(ISTATUS)

!        soil/sea-ice model levels
DO I=1,ICSS
  ILEVS(I)=I
ENDDO
ISTATUS = NF_INQ_VARID   (INCID, 'nlevs', NVARID)
CALL HANDLE_ERR_NC(ISTATUS)
ISTATUS = NF_PUT_VAR_INT (INCID, NVARID, ILEVS)
CALL HANDLE_ERR_NC(ISTATUS)

END ASSOCIATE
END SUBROUTINE SUP1C_NC
