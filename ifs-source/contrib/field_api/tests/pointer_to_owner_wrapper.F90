! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!CHECK THAT WE CAN CREATE A POINTER OF TYPE FIELD_2RB THAN CAN POINT
!TO A FIELD_2RB_OWNER OR A FIELD_2RB_WRAPPER BUT STILL CALL THE
!DESTRUCTOR OF THOSE EXTENDED CLASSES
MODULE POINT_TO
        IMPLICIT NONE
CONTAINS
        SUBROUTINE POINT_TO_WRAPPER(PTR)
                USE FIELD_MODULE
                USE FIELD_FACTORY_MODULE
                USE PARKIND1
                IMPLICIT NONE
                CLASS(FIELD_2RB), POINTER, INTENT(INOUT) :: PTR
                CLASS(FIELD_2RB), POINTER :: W => NULL()
                REAL(KIND=JPRB), POINTER :: D(:,:)
                ALLOCATE(D(10,10))
                CALL FIELD_NEW(W, DATA=D)
                PTR => W                                        
        END SUBROUTINE
        SUBROUTINE POINT_TO_OWNER(PTR)
                USE FIELD_MODULE
                USE FIELD_FACTORY_MODULE
                USE PARKIND1
                IMPLICIT NONE
                CLASS(FIELD_2RB), POINTER, INTENT(INOUT) :: PTR
                CLASS(FIELD_2RB), POINTER :: O => NULL()
                CALL FIELD_NEW(O, LBOUNDS=[1,1],UBOUNDS=[10,10])
                PTR => O                                        
        END SUBROUTINE
END MODULE POINT_TO

PROGRAM POINTER_TO_OWNER_WRAPPER
        USE FIELD_MODULE
        USE POINT_TO
        USE PARKIND1
        IMPLICIT NONE

        TYPE(FIELD_2RB_PTR) :: PTR
        CALL POINT_TO_WRAPPER(PTR%PTR)
        CALL PTR%PTR%FINAL()
        CALL POINT_TO_OWNER(PTR%PTR)
        CALL PTR%PTR%FINAL()
END PROGRAM POINTER_TO_OWNER_WRAPPER
