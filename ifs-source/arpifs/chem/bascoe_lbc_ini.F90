! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE BASCOE_LBC_INI( KNBC, PBCVAL, CDBCNAME )
!**   DESCRIPTION
!     ----------
!
!   Read a set of tables containing lower boundary conditions
!       i.e. volume mixing ratio at lowest model level for stratospheric species
!   The data are typically a subset derived from CMIP6 historical GHG mole
!       fractions (Meinshausen et al. 2016), and encoded in an ASCII file with
!       the following format:
!        - Lines starting with '#' are comments and may be present in header and
!           between tables
!        - The file header contains 2 lines defining the coordinates:
!            - identify coordinates:    'COORDS species month latbounds'
!            - specify dimensions:      'DIMS' <nspecies> <nmonths> <nlatbounds>
!        - One table per species
!            - 1st header row: <species>
!            - 2nd header row: <nlatbounds latbounds>
!            - <nmonths> data rows: '<YYYYmm> <nlatbounds-1 values>
!
!   Part of BASCOE / TM5 routines for IFS chemistry:
!
!
!**   INTERFACE.
!     ----------
!          *BASCOE_LBC_INI* IS CALLED FROM *BASCOE_CHEM_INI*.
!     Input:
!         KNBC:                     Number of defined LBC species
!         PBCVAL(KNBC): (mol/mol)   Default global value
!         CDBCNAME(KNBC):            Name of species
!
!
!**   AUTHOR.
!     -------
!        Yves Christophe, 2018-08-30
!
!-----------------------------------------------------------------------
  USE PARKIND1           , ONLY : JPIM     ,JPRB
  USE YOMHOOK            , ONLY : LHOOK,   DR_HOOK, JPHOOK
  USE YOMLUN             , ONLY : NULOUT
  USE BASCOE_LBC_MODULE  , ONLY : NMONTH_LBC, NLATBOUND_LBC, &
        &                           MONTH_LBC, XLATBOUND_LBC, VALUES_LBC


  IMPLICIT NONE
!-----------------------------------------------------------------------
!   Arguments
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), INTENT(IN)                    :: KNBC
  REAL(KIND=JPRB), INTENT(IN), DIMENSION(KNBC)      :: PBCVAL
  CHARACTER (LEN = 12), INTENT(IN), DIMENSION(KNBC) :: CDBCNAME

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
  INTEGER(KIND=JPIM), PARAMETER :: IUTMP=77

  LOGICAL                       :: ll1              ! reading 1st table
  LOGICAL                       :: LL2

  INTEGER(KIND=JPIM)            :: IOS, ISPC, IMONTH
  INTEGER(KIND=JPIM)            :: JSPC, JBC, JMONTH, JLAT
  INTEGER(KIND=JPIM), DIMENSION(KNBC) :: JDBNAME_READ

  REAL(KIND=JPHOOK)             :: ZHOOK_HANDLE
  !VH REAL(KIND=JPRB), DIMENSION(NLATBOUND_LBC) :: ZLATBOUND_LBC
  REAL(KIND=JPRB), DIMENSION(:), ALLOCATABLE        :: ZLATBOUND_LBC

  character(len=128)            :: clflnm, clline

#include "abor1.intfb.h"
#include "bascoe_cskip.intfb.h"
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('BASCOE_LBC_INI',0,ZHOOK_HANDLE)

!-----------------------------------------------------------------------
!       ... Set filename and open input file
!-----------------------------------------------------------------------
      clflnm = 'BASCOE_LBC.dat'
      OPEN( unit = IUTMP, &
     &      file = clflnm, &
     &      status = 'OLD', form = 'formatted', iostat = IOS )
      IF( IOS /= 0 ) THEN
        CALL ABOR1('BASCOE_LBC_INI'//' Error opening '//TRIM(clflnm))
      ENDIF

      JDBNAME_READ(1:KNBC)=0_JPIM

!-----------------------------------------------------------------------
!       ... Read coordinates, allocate  and initialize arrays
!-----------------------------------------------------------------------
      CALL BASCOE_CSKIP( '#', IUTMP)
      READ( IUTMP, '(a)' , IOSTAT=IOS) clline
      IF( IOS /= 0 .OR. ADJUSTL(clline) /= 'COORDS species month latbounds') THEN
        CALL ABOR1('BASCOE_LBC_INI'//' Error reading COORDS header in '//TRIM(clflnm))
      ENDIF

      CALL BASCOE_CSKIP( '#', IUTMP)
      READ( IUTMP, '(a)' , IOSTAT=IOS) clline
      clline = ADJUSTL(clline)
      IF( IOS /= 0 .OR. (clline(1:4)) /= 'DIMS' ) THEN
        CALL ABOR1('BASCOE_LBC_INI'//' Error reading DIMS header in '//TRIM(clflnm))
      ENDIF
      READ( clline(5:), *, IOSTAT=IOS) ISPC, NMONTH_LBC, NLATBOUND_LBC
      WRITE (NULOUT,*) 'BASCOE_LBC_INI'//': Input dimensions: ',ISPC, NMONTH_LBC, NLATBOUND_LBC
      IF( IOS /= 0  ) THEN
        CALL ABOR1('BASCOE_LBC_INI'//' Error reading dimensions from DIMS header in '//TRIM(clflnm))
      ENDIF

      ALLOCATE( VALUES_LBC(KNBC, NMONTH_LBC, NLATBOUND_LBC-1),   &
              & MONTH_LBC(NMONTH_LBC),                           &
              & XLATBOUND_LBC(NLATBOUND_LBC), ZLATBOUND_LBC(NLATBOUND_LBC))
      DO JBC = 1,KNBC
        VALUES_LBC(JBC,1:NMONTH_LBC,1:NLATBOUND_LBC-1) = PBCVAL(JBC)
      ENDDO

!-----------------------------------------------------------------------
!       ... Read LBC values
!-----------------------------------------------------------------------
      ll1 = .TRUE.
      DO JSPC = 1,ISPC
        LL2= .FALSE.
        ! read table header
        !
        CALL BASCOE_CSKIP( '#', IUTMP)
        READ( IUTMP, * , IOSTAT=IOS) clline
        clline = ADJUSTL(clline)
        DO JBC = 1,KNBC
          IF( TRIM(clline) == TRIM(CDBCNAME(JBC)) )THEN

            ! species name matches one of the declared LBC species
            !
            JDBNAME_READ(JBC)=1_JPIM
            LL2=.TRUE.
            WRITE(NULOUT,*) 'BASCOE_LBC_INI'//': reading LBC for '//TRIM(clline)
            CALL BASCOE_CSKIP( '#', IUTMP)
            READ( IUTMP, * , IOSTAT=IOS) (ZLATBOUND_LBC(JLAT), JLAT=1,NLATBOUND_LBC)
            IF (ll1) THEN
              XLATBOUND_LBC = ZLATBOUND_LBC
            ELSEIF( ANY( XLATBOUND_LBC /= ZLATBOUND_LBC ) ) THEN
              CALL ABOR1('BASCOE_LBC_INI'//' Inconsistent latitude bounds in '//TRIM(clflnm))
            ENDIF

            ! read table content
            !
            DO JMONTH = 1,NMONTH_LBC
              READ( IUTMP, * , IOSTAT=IOS) IMONTH, (VALUES_LBC(JBC,JMONTH,JLAT), JLAT=1,NLATBOUND_LBC-1)
              IF( ll1 ) THEN
                MONTH_LBC(JMONTH) = IMONTH
              ELSEIF( MONTH_LBC(JMONTH) /= IMONTH ) THEN
                CALL ABOR1('BASCOE_LBC_INI'//' Inconsistent YearMonth in '//TRIM(clflnm))
              ENDIF
            ENDDO
            ll1 = .FALSE.

          ENDIF

        ENDDO
        IF (.NOT. LL2) THEN
            ! Skip and go to next entry
           CALL BASCOE_CSKIP( '#', IUTMP)
           READ( IUTMP, * , IOSTAT=IOS) CLLINE
           DO JMONTH = 1,NMONTH_LBC
              READ( IUTMP, * , IOSTAT=IOS) CLLINE
           ENDDO
        ENDIF
      ENDDO

    ! Check for missing LBC
    DO JBC = 1,KNBC
      IF (JDBNAME_READ(JBC) == 0_JPIM) THEN
          ! explicit special treatment for tracers that may currently be missing:
          IF ((TRIM(CDBCNAME(JBC)) == "ocs") .OR.(TRIM(CDBCNAME(JBC)) == "ch2br2")  &
        & .OR.(TRIM(CDBCNAME(JBC)) == "bro") .OR.(TRIM(CDBCNAME(JBC)) == "hbr"   )  &
        & .OR.(TRIM(CDBCNAME(JBC)) == "chbr3")) THEN
            VALUES_LBC(JBC,1:NMONTH_LBC,1:NLATBOUND_LBC-1)= PBCVAL(JBC)
          ELSE
            CALL ABOR1('BASCOE_LBC_INI'//' No LBC available for '//TRIM(CDBCNAME(JBC)))
          ENDIF      
      ENDIF
    ENDDO
    

    CLOSE(IUTMP)

    IF( ll1 ) then
        CALL ABOR1('BASCOE_LBC_INI'//' Error no matching species for LBC in '//TRIM(clflnm))
    ENDIF

    DEALLOCATE( ZLATBOUND_LBC )

IF (LHOOK) CALL DR_HOOK('BASCOE_LBC_INI',1,ZHOOK_HANDLE)

END SUBROUTINE BASCOE_LBC_INI
