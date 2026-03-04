! (C) Copyright 2005- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE GETHEAPSTAT(KOUT, CDLABEL)

USE EC_PARKIND  ,ONLY : JPIM     ,JPRD     ,JPIB

USE MPL_MODULE

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KOUT
CHARACTER(LEN=*), INTENT(IN) :: CDLABEL
INTEGER(KIND=JPIM) :: I, IMYPROC, INPROC, IRET, IOFFSET, II
INTEGER(KIND=JPIM), PARAMETER :: JP_NPROFILE = 9 ! pls. consult ifsaux/utilities/getcurheap.c
INTEGER(KIND=JPIM), PARAMETER :: ISIZE = JP_NPROFILE+1
INTEGER(KIND=JPIB) ILIMIT(ISIZE)
INTEGER(KIND=JPIB) ICNT(ISIZE)
REAL(KIND=JPRD), ALLOCATABLE :: ZSEND(:), ZRECV(:)
INTEGER(KIND=JPIM), ALLOCATABLE :: ICOUNTS(:)
CHARACTER(LEN=1) CLENV
CHARACTER(LEN=80) CLTEXT(0:4)

CALL GET_ENVIRONMENT_VARIABLE("EC_PROFILE_HEAP", CLENV) ! turn ON by export EC_PROFILE_HEAP=1

IF (KOUT >= 0 .AND. CLENV == '1') THEN
  IMYPROC = MPL_MYRANK()
  INPROC  = MPL_NPROC()

  DO I=1,ISIZE
    ILIMIT(I) = I ! power of 10's ; pls. consult ifsaux/utilities/getcurheap.c
  ENDDO

  ALLOCATE(ZSEND(ISIZE))
  ALLOCATE(ZRECV(ISIZE * INPROC))
  ALLOCATE(ICOUNTS(INPROC))

  CLTEXT(0) = "free()/DEALLOCATE -hits per byte range"
  CLTEXT(1) = "malloc()/ALLOCATE -hits per byte range"
  CLTEXT(2) = "Outstanding malloc()/ALLOCATE -hits per byte range"
  CLTEXT(3) = "Outstanding amount of malloc()/ALLOCATE -bytes per byte range"
  CLTEXT(4) = "Average amount of outstanding malloc()/ALLOCATE -bytes per byte range"

  DO II=0,4
    ICNT(:) = 0
    CALL PROFILE_HEAP_GET(ICNT, ISIZE, II, IRET)

    ZSEND(:) = 0
    DO I=1,IRET
      ZSEND(I) = ICNT(I)
    ENDDO
    ZRECV(:) = -1

    ICOUNTS(:) = ISIZE
    CALL MPL_GATHERV(ZSEND(:), KROOT=1, KRECVCOUNTS=ICOUNTS(:), &
                    &PRECVBUF=ZRECV, CDSTRING='GETHEAPSTAT:')

    IF (IMYPROC == 1) THEN
!     Not more than 132 columns, please :-)
      WRITE(KOUT,9000) TRIM(CLTEXT(II)),TRIM(CDLABEL), "Task", &
                     & (ILIMIT(I),I=1,MIN(JP_NPROFILE,9)), "Larger"
9000  FORMAT(/,"Heap Utilization Profile (",A,"): ",A,&
            &/,126("="),&
            &//,(1X,A4,2X,9(:,2X,4X,"< 10^",I1),:,2X,A10))
      WRITE(KOUT,9001)
9001  FORMAT(1X,4("="),2X,10(2X,10("="))/)
      IOFFSET = 0
      DO I=1,INPROC
        ICNT(:) = ZRECV(IOFFSET+1:IOFFSET+ISIZE)
        WRITE(KOUT,'(i5,2x,(10(:,2x,i10)))') I,ICNT(:)
        IOFFSET = IOFFSET + ISIZE
      ENDDO
    ENDIF
  ENDDO

  IF (IMYPROC == 1) THEN
    WRITE(KOUT,'(/,a,/)') 'End of Heap Utilization Profile'
  ENDIF

  DEALLOCATE(ZSEND)
  DEALLOCATE(ZRECV)
  DEALLOCATE(ICOUNTS)
ENDIF
END SUBROUTINE GETHEAPSTAT
