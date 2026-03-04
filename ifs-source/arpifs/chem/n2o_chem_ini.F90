! (C) Copyright 2009- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction

SUBROUTINE N2O_CHEM_INI 

!**   DESCRIPTION 
!     ----------
!
!   BASCOE routine for IFS chemistry : Initialization of photolysis lookup-table
!
!
!
!**   INTERFACE.
!     ----------
!          *N2O_CHEM_INI* IS CALLED FROM *CHEM_INIT*.
!
!
!     AUTHOR.
!     -------
!        VINCENT HUIJNEN    *KNMI*
!        Simon Chabrillat   *BIRA*
!        This source code is based on BASCOE

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 2015-05-08



USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK


IMPLICIT NONE

!-----------------------------------------------------------------------
!*       0.1  ARGUMENTS
!             ---------


! * LOCAL 
REAL(KIND=JPHOOK)                 :: ZHOOK_HANDLE


!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI',0,ZHOOK_HANDLE )

  ! Prepare / read in photolysis table
 CALL N2O_J_TABLES_READ

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI',1,ZHOOK_HANDLE )

CONTAINS

!       Read file with infor for BASCOE photolysis calculations
!       -----------------------------------------------------------
SUBROUTINE N2O_J_TABLES_READ
USE BASCOE_J_MODULE    , ONLY :  NDISS, JNAMES
USE BASCOE_J_TABLES_MODULE    , ONLY :  NZEN => NSZA, NJO3,NJLEV, ZJLEV_TABLES, &
           & O3COL_TABLES, ALJ_TABLES
USE YOMLUN             , ONLY : NULOUT
USE PARKIND1           , ONLY : JPIM     ,JPRB
USE YOMHOOK            ,ONLY  : LHOOK,   DR_HOOK, JPHOOK

IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
INTEGER(KIND=JPIM), PARAMETER  :: IUTMP=77
LOGICAL, PARAMETER             :: LLDEBUG = .true.
INTEGER(KIND=JPIM)             :: ijlev, ijo3, idiss, izen, ios, index
REAL(KIND=JPRB)                :: ZNEWVAL, ZXYJO3(njo3),ZSZA_TABLES(nzen),  ZNEW_SZA(nzen), ZTMP
CHARACTER(len=128)             :: CL_FILENM


#include "abor1.intfb.h"

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI:N2O_J_TABLES_READ',0,ZHOOK_HANDLE )


!-----------------------------------------------------------------------
!  Read the O3 overhad columns corresponding to the ozone standard profiles
!-----------------------------------------------------------------------


      WRITE(*,*) 'J_TABLES_READ: reading the J tables'
      ! if (lj_tmls) THEN
      !   CL_FILENM = '../data_in/j_tables/t_mls/info_j.dat'
      ! else
      !   CL_FILENM = '../data_in/j_tables/t_msis/info_j.dat'
      ! ENDIF
      CL_FILENM = 'info_j.dat'
      OPEN( IUTMP, file=TRIM(ADJUSTL(CL_FILENM)), FORM='FORMATTED',ACTION='READ',status='OLD', iostat=ios )
      IF( ios /= 0 ) THEN
         write(NULOUT,*)'J_TABLES_READ: Error opening >'//TRIM(ADJUSTL(CL_FILENM))//'< , ios=',ios
         CALL ABOR1("J_TABLES_READ: Error opening a J-data input file")
      ENDIF
      CALL CSKIP_N2O( '!', IUTMP )
      READ(IUTMP,'(42x,i6)') index
      IF( index/=njo3 ) THEN
        WRITE(NULOUT,*)'J_TABLES_READ fatal error: njo3 & file mismatch'
        WRITE(NULOUT,*)'njo3, i',njo3, index
        CALL ABOR1("ERROR in J_TABLES_READ")
      ENDIF
      READ(IUTMP,'(42x,i6)') index
      IF(index/=njlev) THEN
        CALL ABOR1("J_TABLES_READ fatal error: njlev & file mismatch")
      ENDIF
      READ(IUTMP,'(42x,i6)') index
      IF( index/=nzen ) THEN
        CALL ABOR1("J_TABLES_READ fatal error: nzen & file mismatch")
      ENDIF
      CALL CSKIP_N2O( '!', IUTMP )
      DO ijlev = 1, njlev
         READ(IUTMP,*) ZJLEV_TABLES(ijlev) &
     &                    , ( O3COL_TABLES(ijlev,ijo3), ijo3 = 1, njo3 )
      ENDDO
      CLOSE( IUTMP )

!-----------------------------------------------------------------------
!  Read the J tables
!-----------------------------------------------------------------------
      ZSZA_TABLES(:) = -999._JPRB
      ZXYJO3(:) = -999._JPRB
      DO idiss = 1, ndiss
         !filenm = STRF_TOLOW( jnames(idiss) )   ! jname
         !filenm = '../data_in/j_tables/sb14a_svn6720/j_' // TRIM(ADJUSTL(filenm)) //'.dat'
         !CL_FILENM = STRF_TOLOW( jnames(idiss) )   ! jname - to lower-case . This doesn't work easily
         CL_FILENM = jnames(idiss)                  ! jname - w/o capitals !!!
         CL_FILENM = 'j_' // TRIM(ADJUSTL(CL_FILENM)) //'.dat'
!         write(NULOUT,*)'J_TABLES_READ fname:',CL_FILENM
         OPEN( IUTMP, file=TRIM(ADJUSTL(CL_FILENM)),  FORM='FORMATTED',ACTION='READ',status = 'OLD', iostat = ios )
         IF( ios /= 0 ) THEN
           write(NULOUT,*)'J_TABLES_READ: Error opening >'//TRIM(ADJUSTL(CL_FILENM))//'< , ios=',ios
           CALL ABOR1("J_TABLES_READ: Error opening a J-data input file")
         ENDIF
         CALL CSKIP_N2O( '!', IUTMP )
         DO ijo3 = 1, njo3
            CALL CSKIP_N2O( '!', IUTMP )
            READ(IUTMP,'(42x,E5.2)') ZNEWVAL
            IF( ZXYJO3(ijo3) == -999._JPRB ) THEN
               ZXYJO3(ijo3) = ZNEWVAL
             ELSEIF( ZNEWVAL /= ZXYJO3(ijo3) .or. &
     &                ZNEWVAL /= O3COL_TABLES(1,ijo3) ) THEN
               CALL ABOR1("J_TABLES_READ error: ZXYJO3 (surf O3col) mismatch")
            ENDIF
            CALL CSKIP_N2O( '!', IUTMP )
            READ(IUTMP,'(32x,99(10x,e6.1))') ( ZNEW_SZA(izen), izen=1,nzen )
            IF( ALL( ZSZA_TABLES(:) == -999._JPRB ) ) THEN
               ZSZA_TABLES(:) = ZNEW_SZA(:)
             ELSEIF( ANY( ZNEW_SZA /=  ZSZA_TABLES) ) THEN
               CALL ABOR1("J_TABLES_READ error: sza (zenith angles) mismatch")
            ENDIF
            CALL CSKIP_N2O( '!', IUTMP )
            DO ijlev = 1, njlev
               READ(IUTMP,'(2F16.1,22(E16.5))') ZNEWVAL, ZTMP, &
     &                   ( alj_tables(izen,ijlev,ijo3,idiss), izen=1, nzen )
               IF( ZNEWVAL /= ZJLEV_TABLES(ijlev) ) THEN
                 CALL ABOR1("J_TABLES_READ error: zjlev (vert level) mismatch")
               ENDIF
            ENDDO
         ENDDO
         CLOSE( IUTMP )
         WRITE(NULOUT,*)'>'//TRIM(ADJUSTL(CL_FILENM))//'< read'
      ENDDO
      WRITE(NULOUT,*) 'J_TABLES_READ: read successfully ',ndiss,' J tables'

!-----------------------------------------------------------------------
!  Store the (natural) LOG of the J tables - set alj_global to default
!-----------------------------------------------------------------------
      alj_tables = MAX( 1.d-30, alj_tables )
      alj_tables = LOG( alj_tables )

!-----------------------------------------------------------------------
! Diagnostics and error reporting
!-----------------------------------------------------------------------
      IF( LLDEBUG ) THEN
         DO izen = 1, nzen
            write(NULOUT,'(a,i3,a,f9.3)') 'J_TABLES_READ: izen= ',izen &
     &          ,' ; sza_tables= ',ZSZA_TABLES(izen)
         ENDDO
         write(NULOUT,*) '   alj_tables(jlev=31=30km,idiss=1=O2,' &
     &       ,'izen=5=60deg,ijo3=2=260DU)= ',alj_tables(5,31,2,1)
      ENDIF

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI:N2O_J_TABLES_READ',1,ZHOOK_HANDLE )

END SUBROUTINE N2O_J_TABLES_READ



SUBROUTINE CSKIP_N2O( CD , KUNIT )
!--------------------------------------------------------------------
!   ... Skip lines starting with char c in opened file iunit
!--------------------------------------------------------------------
USE YOMHOOK            , ONLY : LHOOK, DR_HOOK, JPHOOK
USE PARKIND1           , ONLY : JPIM,  JPRB

      implicit none
!-----------------------------------------------------------------------
!  Local variables
!-----------------------------------------------------------------------
REAL(KIND=JPHOOK)                :: ZHOOK_HANDLE
CHARACTER(len=*), intent(in)   :: CD
INTEGER(KIND=JPIM), intent(in) :: KUNIT
CHARACTER(len=5)               :: CL_C
logical                        :: LLOK
INTEGER(KIND=JPIM)             :: ILEN

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI:CSKIP_N2O',0,ZHOOK_HANDLE )

LLOK = .false.
ILEN = LEN_TRIM(CD)
DO WHILE (.not. LLOK)
   READ(KUNIT,'(a)')CL_C
!   print*, CL_C
   IF (CL_C(:ILEN)/=CD)  LLOK = .true.
ENDDO
BACKSPACE(KUNIT)

IF (LHOOK) CALL DR_HOOK('N2O_CHEM_INI:CSKIP_N2O',1,ZHOOK_HANDLE )
END SUBROUTINE CSKIP_N2O

END SUBROUTINE N2O_CHEM_INI 
