! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE OUTMDLDCP

! ----------------------------------------------------------------------

!**** *OUTMDLDCP* -
!     ----------

!  OUPUT THE MODEL DECOMPOSITON 
!    MPI RANK (IRANK) IN GRIB AS EXTRA PARAMETER 80
!    GRID POINT INDEX (IJ) AS EXTRA PARAMETER 81  

!  THE GRIB DATA WILL BE SAVED TO FILE (see OUTFILE)

!  It will temporarilly reset the output set up to only output what is necessary
!  All is reset at the end.

! ----------------------------------------------------------------------

USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

USE YOWCOUT  , ONLY : JPPFLAG, FFLAG, GFLAG, NFLAG, NIPRMOUT, LFDB
USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_WAM_I
USE YOWGRIBHD, ONLY : PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, NGRBRESS, LNEWLVTP
USE YOWGRID  , ONLY : IJS, IJL, NPROMA_WAM, NCHNK, KIJL4CHNK, IJFROMCHNK, IPRMFROMIJ
USE YOWINTP  , ONLY : GOUT
USE YOWMAP   , ONLY : IRGG, AMONOP, AMOSOP, XDELLA, NLONRGG, BLK2GLO, NIBLO, NGX, NGY, CLDOMAIN 
USE YOWMPP   , ONLY : IRANK, NPRECI 
USE YOWPARAM , ONLY : LL1D
USE YOWPCONS , ONLY : ZMISS
USE YOWSPEC  , ONLY : IJ2NEWIJ
USE YOWTEST  , ONLY : IU06, ITEST
USE YOWTEXT  , ONLY : ICPLEN, CPATH



USE YOMHOOK  , ONLY : LHOOK ,DR_HOOK, JPHOOK
USE YOWGRIB  , ONLY : JPKSIZE_T, &
                    & IGRIB_OPEN_FILE, &
                    & IGRIB_CLOSE_FILE, &
                    & IGRIB_WRITE_BYTES, &
                    & IGRIB_GET_MESSAGE, &
                    & IGRIB_GET_MESSAGE_SIZE, &
                    & IGRIB_CLONE, &
                    & IGRIB_RELEASE

! ----------------------------------------------------------------------

IMPLICIT NONE
#include "abort1.intfb.h"
#include "mpcrtbl.intfb.h"
#include "outgrid.intfb.h"
#include "preset_wgrib_template.intfb.h"
#include "wgribencode.intfb.h"


INTEGER(KIND=JWIM), PARAMETER :: NBITSPERVALUE = 24  !! need a higher precision to code the grid point index 

INTEGER(KIND=JWIM) :: IJ, IPRM, IJGLO, IK, IM, IR, IFLAG, IDT, KSTEP, IX, JSN, I,J,K ,IK1, IK2
INTEGER(KIND=JWIM) :: LFILE, IUOUT, ICOUNT
INTEGER(KIND=JWIM) :: IFCST, IT, IPARAM, ITABLE, IZLEV
INTEGER(KIND=JWIM) :: ICHNK, KIJS, KIJL, IJSB, IJLB
INTEGER(KIND=JWIM) :: IGRIB_HANDLE
INTEGER(KIND=JWIM) :: ISIZE
INTEGER(KIND=JPKSIZE_T) :: KBYTES
INTEGER(KIND=JWIM), ALLOCATABLE :: KGRIB_BUFR(:)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

CHARACTER(LEN= 2) :: MARSTYPE_DUM
CHARACTER(LEN=14) :: CDATE
CHARACTER(LEN=296) :: OUTFILE

LOGICAL :: LLCREATE, LFDBBAK
LOGICAL :: LGRHDIFS_DUM, LPADPOLES_DUM, LRSTST0_DUM
LOGICAL, DIMENSION(JPPFLAG) :: FFLAGBAK, GFLAGBAK, NFLAGBAK
! TEMPORARY ARRAYS FOR LOCAL MASK AND GLOBAL INDICES
INTEGER(KIND=JWIM), DIMENSION(NGX,NGY) :: IGLOBAL
INTEGER(KIND=JWIM), DIMENSION(IJS:IJL) :: NLOCMSK, NGLOIND
REAL(KIND=JWRB), ALLOCATABLE, DIMENSION(:, :, :) :: BOUT

! ----------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',0,ZHOOK_HANDLE)

! Output will be save in 
IF ( ICPLEN > 0 ) THEN
  OUTFILE = CPATH(1:ICPLEN)//'/wam_model_mpi_decomposition.grb'
ELSE
  OUTFILE = 'wam_model_mpi_decomposition.grb'
ENDIF
LFILE=LEN_TRIM(OUTFILE)

WRITE(IU06,*) ' OUTMDLDCP : The MPI decomposition is written to ', OUTFILE(1:LFILE)

! save output parameter selection (it will be overwritten and then reset)
FFLAGBAK(:) = FFLAG(:)
GFLAGBAK(:) = GFLAG(:)
NFLAGBAK(:) = NFLAG(:)
LFDBBAK = LFDB

! Create a new output parameter selection that will only output
! parameters relevant for the description of the model decomposition
FFLAG(:) = .FALSE.
GFLAG(:) = .FALSE.
NFLAG(:) = .FALSE.
LFDB = .FALSE. 


! Use the extra field codes (currently the last 5 fields of the list of potential output parameters:
!! see below that everyhting is hardcoded from JPPFLAG-8 to JPPFLAG
!! We have to use to parameters because of the limitation of grib)  :
! MPI RANK:
GFLAG(JPPFLAG-8) = .TRUE.
GFLAG(JPPFLAG-7) = .TRUE.

! LAND sea mask
GFLAG(JPPFLAG-6) = .TRUE.

! GLOBAL GRID POINT INDEX (land and sea) (as used inside the wave model):
GFLAG(JPPFLAG-5) = .TRUE.
GFLAG(JPPFLAG-4) = .TRUE.

! GRID POINT INDEX (sea only) (as used inside the wave model):
GFLAG(JPPFLAG-3) = .TRUE.
GFLAG(JPPFLAG-2) = .TRUE.

! GLOBAL POINT INDEX (see only) (as used for non parallel binary input/output)
! It is different than the local grid point index if LL1D is false
GFLAG(JPPFLAG-1) = .TRUE.
GFLAG(JPPFLAG) = .TRUE.


! Set output parameter mapping
CALL MPCRTBL

ALLOCATE ( BOUT(NPROMA_WAM, MAX(NIPRMOUT,1), NCHNK) )



! SETUP GLOBAL INDICES 1D ARRAY + LOCAL MASK
! SINCE WE ONLY HAVE OCEAN POINTS THE LOCAL MASK IS ONE FOR ALL POINTS
! SEE INITNEMOCPL
IGLOBAL(:,:) = -1
K = 0
DO J = 1, NGY
  DO I = 1, NLONRGG(NGY-J+1)
    K = K + 1
    IGLOBAL(I,NGY-J+1) = K
  ENDDO
ENDDO
!
DO IJ = IJS, IJL
   IX          = BLK2GLO%IXLG(IJ)
   JSN         = BLK2GLO%KXLT(IJ)
   NLOCMSK(IJ) = 1
   NGLOIND(IJ) = IGLOBAL(IX,JSN)
ENDDO

! Defining the ouput fields:
! -------------------------


IK1 = 0
IK2 = 0
DO ICHNK = 1, NCHNK

   KIJS = 1
   IJSB = IJFROMCHNK(KIJS,ICHNK)
   KIJL = KIJL4CHNK(ICHNK)
   IJLB = IJFROMCHNK(KIJL,ICHNK) 

  IR = 0

  ! MPI RANK:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = REAL(MOD(IRANK,65536), JWRB)

  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = REAL(IRANK/65536, JWRB)


  ! LAND/SEA MASK:
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = REAL(NLOCMSK(IJSB:IJLB), JWRB)


  ! GLOBAL GRID INDEX (sea and land points):
  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = REAL(MOD(NGLOIND(IJSB:IJLB),65536), JWRB)

  IR = IR + 1
  BOUT(KIJS:KIJL, IR, ICHNK) = REAL(NGLOIND(IJSB:IJLB)/65536, JWRB)


  ! LOCAL SEA POINT BLOCK INDEX:
  IR = IR + 1
  DO IPRM = KIJS, KIJL
    IK1 = IK1 + 1
    BOUT(IPRM, IR, ICHNK) = REAL(MOD(IK1,65536), JWRB)
  ENDDO

  IR = IR + 1
  DO IPRM = KIJS, KIJL
    IK2 = IK2 + 1
    BOUT(IPRM, IR, ICHNK) = REAL(IK2/65536, JWRB)
  ENDDO


  ! GLOBAL SEA POINT INDEX:
  IR = IR + 1
  IF ( .NOT. LL1D) THEN
    DO IJGLO = 1, NIBLO
      IJ = IJ2NEWIJ(IJGLO) 
      IF ( IJ >= IJSB .AND. IJ <= IJLB ) THEN 
         IPRM = IPRMFROMIJ(IJ)
         BOUT(IPRM, IR, ICHNK) = REAL(MOD(IJGLO,65536), JWRB)
      ENDIF
    ENDDO
  ELSE
    DO IJ = IJSB, IJLB
      IPRM = IPRMFROMIJ(IJ)
      BOUT(IPRM, IR, ICHNK) = REAL(MOD(IJ,65536), JWRB)
    ENDDO
  ENDIF

  IR = IR + 1
  IF ( .NOT. LL1D) THEN
    DO IJGLO = 1, NIBLO
      IJ = IJ2NEWIJ(IJGLO) 
      IF ( IJ >= IJSB .AND. IJ <= IJLB ) THEN 
         IPRM = IPRMFROMIJ(IJ)
         BOUT(IPRM, IR, ICHNK) = REAL(IJGLO/65536, JWRB)
      ENDIF
    ENDDO
  ELSE
    DO IJ = IJSB, IJLB
      IPRM = IPRMFROMIJ(IJ)
      BOUT(IPRM, IR, ICHNK) = REAL(IJ/65536, JWRB)
    ENDDO
  ENDIF



  IF ( IR /= NIPRMOUT ) THEN
    WRITE(IU06,*) '******************************************'
    WRITE(IU06,*) ' OUTMDLDCP : ABORTING' 
    WRITE(IU06,*) ' IR /= NIPRMOUT ', IR, NIPRMOUT 
    WRITE(IU06,*) '******************************************'
    CALL ABORT1
  ENDIF

ENDDO


! Gather data for output (to IRANK = 1)
CALL OUTGRID(BOUT)



IF(IRANK == 1) THEN
  ! Grib output to file:
  CALL IGRIB_OPEN_FILE(IUOUT, OUTFILE(1:LFILE), 'w')

  ! Prepare grib template
  LLCREATE = .TRUE.
  CALL PRESET_WGRIB_TEMPLATE("I", NGRIB_HANDLE_WAM_I, LLCREATE=LLCREATE, NBITSPERVALUE=NBITSPERVALUE )


  ! keep looping over all posible output varaibles (as in outint)
  ITABLE=140
  ICOUNT=0
  DO IFLAG = JPPFLAG-8, JPPFLAG


    IF ( GFLAG(IFLAG) ) THEN
      ICOUNT = ICOUNT+1

      SELECT CASE (IFLAG)
        CASE (JPPFLAG-8)
          IPARAM=80
          IZLEV=0
        CASE (JPPFLAG-7)
          IPARAM=80
          IZLEV=1
        CASE (JPPFLAG-6)
          IPARAM=81
          IZLEV=0
        CASE (JPPFLAG-5)
          IPARAM=82
          IZLEV=0
        CASE (JPPFLAG-4)
          IPARAM=82
          IZLEV=1
        CASE (JPPFLAG-3)
          IPARAM=83
          IZLEV=0
        CASE (JPPFLAG-2)
          IPARAM=83
          IZLEV=1
        CASE (JPPFLAG-1)
          IPARAM=84
          IZLEV=0
        CASE (JPPFLAG)
          IPARAM=84
          IZLEV=1
      END SELECT

      IK = 0
      IM = 0

      CDATE = '20220220000000' ! set any date and time as the start date of the run might still be unknown
      IFCST = 0
      MARSTYPE_DUM = 'an'

      LGRHDIFS_DUM = .FALSE.
      IDT = 0
      KSTEP = 0
      LPADPOLES_DUM = .FALSE.
      LRSTST0_DUM = .FALSE.

      IGRIB_HANDLE=-99
      CALL IGRIB_CLONE(NGRIB_HANDLE_WAM_I, IGRIB_HANDLE)

      ! grib coding
      CALL WGRIBENCODE( IU06, ITEST, &
 &                      NGX, NGY, &
 &                      GOUT(ICOUNT,:,:),  &
 &                      ITABLE, IPARAM, &
 &                      IZLEV, &
 &                      IK, IM, &
 &                      CDATE, IFCST, MARSTYPE_DUM, &
 &                      PPMISS, PPEPS, PPREC, PPRESOL, PPMIN_RESET, NTENCODE, &
 &                      LGRHDIFS_DUM, &
 &                      IDT, &
 &                      NGRBRESS, LNEWLVTP, LPADPOLES_DUM, &
 &                      SIZE(NLONRGG), NLONRGG, IRGG, &
 &                      AMONOP, AMOSOP, XDELLA, CLDOMAIN, &
 &                      KSTEP, LRSTST0_DUM, &
 &                      ZMISS, &
 &                      IGRIB_HANDLE )

      ! Save grib data to file

      CALL IGRIB_GET_MESSAGE_SIZE(IGRIB_HANDLE, KBYTES)
      ISIZE = (KBYTES+NPRECI-1)/NPRECI
      ALLOCATE(KGRIB_BUFR(ISIZE))

      CALL IGRIB_GET_MESSAGE(IGRIB_HANDLE, KGRIB_BUFR)

      CALL IGRIB_WRITE_BYTES(IUOUT, KGRIB_BUFR, KBYTES)

      CALL IGRIB_RELEASE(IGRIB_HANDLE)

      DEALLOCATE(KGRIB_BUFR)

    ENDIF
  ENDDO

  CALL IGRIB_CLOSE_FILE(IUOUT)
  CALL IGRIB_RELEASE(NGRIB_HANDLE_WAM_I)

  IF (ALLOCATED(GOUT)) DEALLOCATE(GOUT)

ENDIF


DEALLOCATE ( BOUT )

! Restore output configuration: 
FFLAG(:) = FFLAGBAK(:)
GFLAG(:) = GFLAGBAK(:)
NFLAG(:) = NFLAGBAK(:)
LFDB = LFDBBAK
!  reset output field mapping
CALL MPCRTBL

WRITE(IU06,*) ' '

IF (LHOOK) CALL DR_HOOK('OUTMDLDCP',1,ZHOOK_HANDLE)

END SUBROUTINE OUTMDLDCP
