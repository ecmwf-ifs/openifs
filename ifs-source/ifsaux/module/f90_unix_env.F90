! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
#ifdef NAG
!-- The NAG-compiler has its own proper F90_UNIX_* modules
!   We need this for non-NAG compilers to satisfy some
!   make/build-environments only
#define F90_UNIX_ENV NAGDUMMY_F90_UNIX_ENV
#endif
MODULE F90_UNIX_ENV
LOGICAL L_NAGDUMMY
END MODULE F90_UNIX_ENV
#ifdef NAG
#undef F90_UNIX_ENV
#endif
