! (C) Copyright 2022- ECMWF.
! (C) Copyright 2022- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM GATHER_SCATTER
        !TEST THAT FIELD_GATHSCAT ONLY MODIFY VALUES THAT HAVE BEEN FILTERED
        USE FIELD_MODULE
        USE FIELD_FACTORY_MODULE
        USE FIELD_ACCESS_MODULE
        USE FIELD_GATHSCAT_MODULE
        USE PARKIND1
        IMPLICIT NONE

        TYPE(FIELD_GATHSCAT):: FGS
        CLASS(FIELD_2LM), POINTER :: FTRIG => NULL()
        INTEGER, PARAMETER :: NPROMA = 10, NGPBLKS = 4
        LOGICAL :: TRIG(NPROMA,NGPBLKS)
        INTEGER ::I,J

        INTEGER(KIND=JPIM), ALLOCATABLE :: D(:,:)
        CLASS(FIELD_2IM),POINTER :: FD => NULL()
        INTEGER(KIND=JPIM),POINTER :: FILTERED_D(:,:) => NULL()


        !CREATE A FILTER TO USE WITH THE GATHSCAT
        TRIG=.FALSE.
        DO I=1,NPROMA
        IF(MOD(I,2)==0)THEN !ONLY MODIFY BLOCKS
                TRIG(I,:)=.TRUE.
        ENDIF
        ENDDO
        CALL FIELD_NEW(FTRIG, DATA=TRIG)

        !CREATE THE FIELD TO BE FILTERED BY GATHSCAT
        ALLOCATE(D(NPROMA,NGPBLKS))
        D=1
        CALL FIELD_NEW(FD, DATA=D)

        !FILTER DATA, WE GET A POINTER TO A CONTIGUOUS ARRAY CONTAINING ONLY THE FILTERED DATA
        CALL FGS%INIT(FTRIG, NPROMA*NGPBLKS)
        FILTERED_D=>GATHER_HOST_DATA_RDWR(FGS, FD)
        FILTERED_D=2 !NOT ALL OF D WILL BE MODIFIED, ONLY THE FILTERED DATA

        !ACTUALLY UPDATE THE D ARRAY WITH THE MODIFIED DATA
        CALL FGS%SCATTER()

        DO I=1,NPROMA
        IF(MOD(I,2)==0)THEN
                IF (.NOT. ALL(D(I,:)==2))THEN
                        ERROR STOP
                ENDIF
        ELSE
                IF (.NOT. ALL(D(I,:)==1))THEN
                        ERROR STOP
                ENDIF
        ENDIF
        ENDDO

        CALL FIELD_DELETE(FTRIG)

END PROGRAM
