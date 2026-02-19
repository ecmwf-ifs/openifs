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

module obsop_sets
!-----------------------------------------------------------------------
!   Looks after observation "sets" which  are groups of similar observations;
!   there is one call to the top level observation operator "hop" per set.
!   The set is the basic unit of observations for parallelism and
!   load-balancing.

!   The basic properties of a set (e.g. obstype, timeslot) are intended to be
!   consistent across all members of the set. For satellites, although
!   not conventional data, there is also a single satgrp table "group",
!   codetype and sensor. But conventional can have multiple codetypes
!   in one set (hence "codetype_sat") and does not have a "group" or "sensor"

!   Author.
!   -------
!   Alan Geer 09-Oct-2015 - Call the old ECSET functions from this object;
!                           private ownership of set-related former global data
!                           variables and encapsulate a lot of copy-and-paste code

!   Modifications.
!   --------------

!   Niels Bormann    01-Feb-2016 Include zenith angle when available
!   A. Geer          04-Apr-2016 Use type-bound procedures for a tidier object
!   F. Suzat         08-Apr-2018 Report 43t2bf to make canari running
!   P. Lean          23-Mar-2017 Added set boundaries for OOPS Jo-table
!   N. Bormann       05-Aug-2020 Support for different types of 2dGOM treatments
!------------------------------------------------------------------------

use parkind1, only : jpim, jprb
use yomhook , only : lhook, dr_hook, jphook
use yomct0   , only : lcanari
use yomcoctp , only : nsatob, nsatem, nlimb, nlidar,ndwlpcd, nallsky, &
  & ngthrb,nssmi,ntcwc
use yomsats  , only : satgrp_table, satobgrp_table, chmethod, mxmethod
use yomlimb,   only : y_limbgrp_table
use gom2d_support, only : jp_treat_2dgom_unspec, jp_treat_2dgom_limb_plane
use rttov_const  , only : ninst, inst_name
use yomdimo  , only : nmxset
use dbase_mod, only : dbase

implicit none
private
save

integer(kind=jpim), parameter, public :: len_set_name=60

type, public :: type_set_info
  integer(kind=jpim) :: id            ! ID number of the set
  integer(kind=jpim) :: lnset         ! number of observations (header entries)
  integer(kind=jpim) :: mxbdy         ! maximum body length (body entries per header)
  integer(kind=jpim) :: sumbdy        ! total number of active body entries in set
  integer(kind=jpim) :: tslot         ! timeslot
  integer(kind=jpim) :: obstype       ! a grouping used for jo table and gom settings
  integer(kind=jpim) :: group         ! index into "satgrp" and "satobgrp" tables
  integer(kind=jpim) :: codetype_sat  ! limited to satellites, used for odb access
  integer(kind=jpim) :: sensor        ! true sensor ID, satellites only (else -1)
  integer(kind=jpim) :: sensor_safe   ! sensor ID range limited to 0:NINST-1 for use with e.g. VarQC arrays
  integer(kind=jpim) :: start_global_row ! index of first row in set in the global observation vector
  integer(kind=jpim) :: end_global_row   ! index of last row in set in the global observation vector
  integer(kind=jpim), pointer :: mapomm(:)=>null() ! points to the GOM indices needed to create a GOM_PLUS
  integer(kind=jpim) :: treat_2dgom_id ! ID indicating 2d-GOM treatment (see gom2d_support)
  real(kind=jprb), pointer :: zenith(:)=>null()    ! zenith angle in degrees for SATEMS
end type

type, public :: class_obsop_sets
  private

  integer(kind=jpim) :: ntotal_sets = 0_jpim

  type(type_set_info), pointer :: yset_info(:)=>null()
  integer(kind=jpim), pointer :: mmapomm(:,:)=>null()  ! Indexes into the GOM
  real(kind=jprb), pointer :: zenith(:,:)=>null() ! zenith angle in degrees for SATEMS
end type
  !FORECAST_ONLY or OIFS BUILD, all obsop_set routines removed. Retain type definitions

end module obsop_sets
