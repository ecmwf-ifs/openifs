! (C) Copyright 2005- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

REAL(KIND=JPRB) :: VALBF(12)
REAL(KIND=JPRB) :: VLAIL(12)
REAL(KIND=JPRB) :: VLAIH(12)


REAL(KIND=JPRB) ::  VR0L
REAL(KIND=JPRB) ::  VR0H
REAL(KIND=JPRB) :: VRSML
REAL(KIND=JPRB) :: VRSMH


REAL(KIND=JPRB) :: VZ0F
REAL(KIND=JPRB) :: VITM
REAL(KIND=JPRB) :: VGEO
REAL(KIND=JPRB) :: VZ0H
REAL(KIND=JPRB) :: VTVL
REAL(KIND=JPRB) :: VTVH
REAL(KIND=JPRB) :: VCVL
REAL(KIND=JPRB) :: VCVH
REAL(KIND=JPRB) :: VCUR
REAL(KIND=JPRB) :: VSOTY
REAL(KIND=JPRB) :: VSDOR
REAL(KIND=JPRB) :: VCO2TYP


REAL(KIND=JPRB) :: VSST
REAL(KIND=JPRB) :: VCI
REAL(KIND=JPRB) :: VEMISF
REAL(KIND=JPRB) :: VVEGF
REAL(KIND=JPRB) :: VLAIF
REAL(KIND=JPRB) :: VLDEPTH
REAL(KIND=JPRB) :: VCLAKE



!*    ------------------------------------------------------------------
NAMELIST/NAMGPD1S/VALBF,VLAIL,VLAIH,&
                 &VRSML, VRSMH,VVEGF,VLAIF,& 
                 &VZ0F ,VITM ,VZ0H ,VCVL ,VCVH ,VCUR  ,VTVL ,VCO2TYP, VTVH ,&
                 &VSST ,VCI  ,VEMISF,VSOTY,VSDOR,VLDEPTH,VCLAKE, VR0L, VR0H

!     ------------------------------------------------------------------
