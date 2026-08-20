MODULE CMF_UTILS_MOD
!==========================================================
!* PURPOSE: Shared ulitity functions/subroutines for CaMa-Flood
! (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
!
! Licensed under the Apache License, Version 2.0 (the "License");
!   You may not use this file except in compliance with the License.
!   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software distributed under the License is 
!  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
! See the License for the specific language governing permissions and limitations under the License.
!==========================================================
USE PARKIND1,                ONLY: JPIM,   JPRB, JPRM
USE YOS_CMF_INPUT,           ONLY: LOGNAM, DMIS, RMIS, NX,NY
USE YOS_CMF_MAP,             ONLY: NSEQMAX, NSEQALL
IMPLICIT NONE
CONTAINS
!####################################################################
! map related subroutines & functions
!-- VEC2MAP     : convert 1D vector data -> 2D map data (REAL*4)
!-- VEC2MAPD    : convert 1D vector data -> 2D map data (REAL*8)
!-- MAP2VEC     : convert 2D map data -> 1D vector data (REAL*4)
!-- MAP2VECD    : convert 2D map data -> 1D vector data (REAL*8)
!-- MAP2VECI    : convert 2D map data -> 1D vector data (Integer)
!
! time related subroutines & functions
! -- MIN2DATE  : calculate DATE of KMIN from base time (YYYY0,MM0,DD0)
! -- DATE2MIN  : convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
! -- SPLITDATE : splite date (YYYYMMDD) to (YYYY,MM,DD)
! -- SPLITHOUR : split hour (HHMM) to (HH,MM)
! -- IMDAYS    : function to calculate days in a monty IMDAYS(IYEAR,IMON)
!
! endian conversion
!-- CONV_END    : Convert 2D Array endian (REAL4)
!-- CONV_ENDI   : Convert 2D Array endian (Integer)
!-- ENDIAN4R    : byte swap (REAL*4)
!-- ENDIAN4I    : byte swap (Integer)
!
! file I/O
!-- INQUIRE_FID : inruire unused file FID
!-- NCERROR     : netCDF I/O wrapper
!####################################################################



!####################################################################
! map related subroutines & functions
!-- VEC2MAP     : convert 1D vector data -> 2D map data (REAL*4)
!-- VEC2MAPD    : convert 1D vector data -> 2D map data (REAL*8)
!-- MAP2VEC     : convert 2D map data -> 1D vector data (REAL*4)
!-- MAP2VECD    : convert 2D map data -> 1D vector data (REAL*8)
!-- MAP2VECI    : convert 2D map data -> 1D vector data (Integer)
!####################################################################
SUBROUTINE VEC2MAP(D2VEC,R2MAP)
USE YOS_CMF_MAP,             ONLY: I1SEQX,I1SEQY
IMPLICIT NONE
!* input/output
REAL(KIND=JPRB),INTENT(IN)      :: D2VEC(NSEQMAX,1)
REAL(KIND=JPRM),INTENT(OUT)     :: R2MAP(NX,NY)
!* local variable
!$ SAVE
INTEGER(KIND=JPIM)              ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
!================================================
R2MAP(:,:) = RMIS
!$OMP PARALLEL DO
DO ISEQ=1,NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  R2MAP(IX,IY) = REAL(D2VEC(ISEQ,1),KIND=JPRM)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE VEC2MAP
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE VEC2MAPD(D2VEC,D2MAP)
USE YOS_CMF_MAP,             ONLY: I1SEQX,I1SEQY
IMPLICIT NONE
!* input/output
REAL(KIND=JPRB),INTENT(IN)      :: D2VEC(NSEQMAX,1)
REAL(KIND=JPRB),INTENT(OUT)     :: D2MAP(NX,NY)
!* local variable
!$ SAVE
INTEGER(KIND=JPIM)              ::  IX,IY,ISEQ
!$OMP THREADPRIVATE                (IX,IY)
!================================================
D2MAP(:,:) = DMIS
!$OMP PARALLEL DO
DO ISEQ=1,NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  D2MAP(IX,IY) = D2VEC(ISEQ,1)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE VEC2MAPD
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE MAP2VEC(R2TEMP,D2VAR)
USE YOS_CMF_MAP,             ONLY: I1SEQX,I1SEQY
IMPLICIT NONE
!* input/output
REAL(KIND=JPRM),INTENT(IN)      :: R2TEMP(NX,NY)
REAL(KIND=JPRB),INTENT(OUT)     :: D2VAR(NSEQMAX,1)
!* local variable
!$ SAVE
INTEGER(KIND=JPIM)              :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
DO ISEQ=1,NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  D2VAR(ISEQ,1) = REAL(R2TEMP(IX,IY),KIND=JPRB)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE MAP2VEC
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE MAP2VECD(D2TEMP,D2VAR)
USE YOS_CMF_MAP,             ONLY: I1SEQX,I1SEQY
IMPLICIT NONE
!* input/output
REAL(KIND=JPRB),INTENT(IN)      :: D2TEMP(NX,NY)
REAL(KIND=JPRB),INTENT(OUT)     :: D2VAR(NSEQMAX,1)
!* local variable
!$ SAVE
INTEGER(KIND=JPIM)              :: IX,IY, ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
DO ISEQ=1,NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  D2VAR(ISEQ,1) = D2TEMP(IX,IY)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE MAP2VECD
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE MAP2VECI(I2TEMP,I2VAR)
USE YOS_CMF_MAP,             ONLY: I1SEQX,I1SEQY
IMPLICIT NONE
!* input/output
INTEGER(KIND=JPIM),INTENT(IN)   :: I2TEMP(NX,NY)
INTEGER(KIND=JPIM),INTENT(OUT)  :: I2VAR(NSEQMAX,1)
!* local variable
!$ SAVE
INTEGER(KIND=JPIM)              :: IX,IY,ISEQ
!$OMP THREADPRIVATE               (IX,IY)
!================================================
!$OMP PARALLEL DO
DO ISEQ=1,NSEQALL
  IX=I1SEQX(ISEQ)
  IY=I1SEQY(ISEQ)
  I2VAR(ISEQ,1) = I2TEMP(IX,IY)
ENDDO
!$OMP END PARALLEL DO
END SUBROUTINE MAP2VECI
!####################################################################





!####################################################################
! time related subroutines & functions
! -- MIN2DATE  : calculate DATE of KMIN from base time (YYYY0,MM0,DD0)
! -- DATE2MIN  : convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
! -- SPLITDATE : splite date (YYYYMMDD) to (YYYY,MM,DD)
! -- SPLITHOUR : split hour (HHMM) to (HH,MM)
! -- IMDAYS    : function to calculate days in a monty IMDAYS(IYEAR,IMON)
!==========================================================
SUBROUTINE MIN2DATE(IMIN,YYYYMMDD,HHMM)
!  Return YYYYMMDD and HHMM for IMIN
USE YOS_CMF_TIME,            ONLY: YYYY0, MM0, DD0
IMPLICIT NONE
! local
INTEGER(KIND=JPIM),INTENT(IN)   :: IMIN      !!  input minutes
INTEGER(KIND=JPIM),INTENT(OUT)  :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(OUT)  :: HHMM
INTEGER(KIND=JPIM)              :: YYYY,MM,DD,HH,MI,NDAYS,NDM,ID
INTEGER(KIND=JPIM)              :: D2MIN                   ! minutes in one day
PARAMETER                         (D2MIN=1440)
!================================================
YYYYMMDD = 0
HHMM     = 0

NDAYS = IMIN/D2MIN              !! days  in IMIN : 1440 = (minutes in a day)
MI    = MOD(IMIN,D2MIN)
HH    = INT(MI/60)              !! hours in IMIN
MI    = MOD(MI,60)              !! mins  in IMIN

YYYY  = YYYY0
MM    = MM0
DD    = DD0
NDM   = IMDAYS(YYYY,MM)      !! number of days in a month

! WRITE(LOGNAM,*)  YYYY,MM,DD
DO ID=1,NDAYS
  DD=DD+1
  IF ( DD .GT. NDM ) THEN
    MM=MM+1
    DD=1
    IF ( MM .GT. 12 ) THEN
      MM=1
      YYYY=YYYY+1
    ENDIF
    NDM=IMDAYS(YYYY,MM)
  ENDIF
ENDDO

HHMM     = HH*100+MI
YYYYMMDD = YYYY*10000+MM*100+DD
END SUBROUTINE MIN2DATE
!==========================================================
!+
!+
!+
!==========================================================
FUNCTION DATE2MIN(YYYYMMDD,HHMM)
! convert (YYYYMMDD,HHMM) to KMIN from base time (YYYY0,MM0,DD0)
USE YOS_CMF_TIME,            ONLY: YYYY0
IMPLICIT NONE
INTEGER(KIND=JPIM)              :: DATE2MIN
INTEGER(KIND=JPIM),INTENT(IN)   :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(IN)   :: HHMM
INTEGER(KIND=JPIM)              :: YYYY,MM,DD,HH,MI
INTEGER(KIND=JPIM)              :: IY,IM
INTEGER(KIND=JPIM)              :: D2MIN                   ! minutes in one day
PARAMETER                         (D2MIN=1440)
!================================================
DATE2MIN = 0
CALL SPLITDATE(YYYYMMDD,YYYY,MM,DD)
HH = HHMM/100                          !! hour
MI = HHMM-HH*100                       !! minute
!============================
IF ( YYYY .LT. YYYY0) THEN
  WRITE(LOGNAM,*) 'DATE2MIN: YYYY .LT. YYYY0: Date Problem', YYYY,YYYY0
  STOP
ENDIF
IF ( MM.LT.1 .or. MM .GT. 12 ) THEN
  WRITE(LOGNAM,*) 'DATE2MIN: MM:    Date Problem', YYYYMMDD, HHMM
  STOP
ENDIF
IF ( DD.LT.1 .or. DD .GT. IMDAYS(YYYY,MM)) THEN
  WRITE(LOGNAM,*) 'DATE2MIN: DD:    Date Problem', YYYYMMDD, HHMM
  STOP
ENDIF
IF ( HH.LT.0 .or. HH .GT. 24) THEN
  WRITE(LOGNAM,*) 'DATE2MIN: HH:    Date Problem', YYYYMMDD, HHMM
  STOP
ENDIF
IF ( MI.LT.0 .or. MI .GT. 60) THEN
  WRITE(LOGNAM,*) 'DATE2MIN: MI:    Date Problem', YYYYMMDD, HHMM
  STOP
ENDIF

IY=YYYY0
DO WHILE (IY .LT. YYYY)
  DO IM=1,12
    DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
  ENDDO
  IY=IY+1
ENDDO
IM=1
DO WHILE (IM .LT. MM )
  DATE2MIN=DATE2MIN+IMDAYS(IY,IM)*D2MIN
  IM=IM+1
ENDDO

DATE2MIN = DATE2MIN + (DD-1)*D2MIN
DATE2MIN = DATE2MIN + HH*60 + MI

END FUNCTION DATE2MIN
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SPLITDATE(YYYYMMDD,YYYY,MM,DD)
! sprit YYYYMMDD to (YYYY,MM,DD)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: YYYYMMDD
INTEGER(KIND=JPIM),INTENT(OUT)  :: YYYY,MM,DD
!================================================
YYYY =  YYYYMMDD/10000
MM   = (YYYYMMDD - YYYY*10000) / 100
DD   =  YYYYMMDD -(YYYY*10000+MM*100)
END SUBROUTINE SPLITDATE
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE SPLITHOUR(HHMM,HH,MI)
! sprit YYYYMMDD to (YYYY,MM,DD)
IMPLICIT NONE
INTEGER(KIND=JPIM),INTENT(IN)   :: HHMM
INTEGER(KIND=JPIM),INTENT(OUT)  :: HH,MI
!================================================
HH=INT(HHMM/100)
MI=INT(HHMM-HH*100)
END SUBROUTINE SPLITHOUR
!==========================================================
!+
!+
!+
!==========================================================
FUNCTION IMDAYS(IYEAR,IMON)
!! days in month
USE YOS_CMF_INPUT,           ONLY: LLEAPYR
IMPLICIT NONE
INTEGER(KIND=JPIM)              :: IMDAYS
INTEGER(KIND=JPIM),INTENT(IN)   :: IYEAR
INTEGER(KIND=JPIM),INTENT(IN)   :: IMON
INTEGER(KIND=JPIM)              :: ND(12)
DATA ND /31,28,31,30,31,30,31,31,30,31,30,31/
!================================================
IMDAYS=ND(IMON)
IF ( IMON == 2 .and. LLEAPYR ) THEN
  IF ( MOD(IYEAR,400) == 0 .OR. (MOD(IYEAR,100) .NE. 0 .AND. MOD(IYEAR,4) .EQ. 0 )) IMDAYS=29
ENDIF
END FUNCTION IMDAYS
!==========================================================




!####################################################################
! endian conversion
!-- CONV_END    : Convert 2D Array endian (REAL4)
!-- CONV_ENDI   : Convert 2D Array endian (Integer)
!-- ENDIAN4R    : byte swap (REAL*4)
!-- ENDIAN4I    : byte swap (Integer)
!####################################################################
SUBROUTINE CONV_END(R2TEMP,NX,NY)
!-- Convert 2D Array endian (REAL4)
IMPLICIT NONE
!* input/output
INTEGER(KIND=JPIM),INTENT(IN)   :: NX,NY
REAL(KIND=JPRM),INTENT(INOUT)   :: R2TEMP(NX,NY)
!* local variables
INTEGER(KIND=JPIM)              :: IY,IX
!================================================
DO IY=1, NY
  DO IX=1, NY
    CALL ENDIAN4R(R2TEMP(IX,IY))
  END DO
END DO
END SUBROUTINE CONV_END
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE CONV_ENDI(I2TEMP,NX,NY)
!-- Convert 2D Array endian (INTEGER)
IMPLICIT NONE
!+ input/output
INTEGER(KIND=JPIM),INTENT(IN)     :: NX,NY
INTEGER(KIND=JPIM),INTENT(INOUT)  :: I2TEMP(NX,NY)
!* local variables
INTEGER(KIND=JPIM)                :: IY,IX
!================================================
DO IY=1, NY
  DO IX=1, NY
    CALL ENDIAN4I(I2TEMP(IX,IY))
  END DO
END DO
END SUBROUTINE CONV_ENDI
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE ENDIAN4R( realIn )
!! Byte Swap 
!
! Adpated from: http://www.cgd.ucar.edu/cas/software/endian.html
!           FILE: SUBR_native_4byte_real.f90
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003 
!  LAST MODIFIED: 29 April 2003
IMPLICIT NONE
!* input/output
REAL(KIND=JPRM), INTENT(INOUT)  :: realIn
!* Local variables (generic 32 bit INTEGER spaces):
INTEGER                         :: i_element
INTEGER                         :: i_element_br
!================================================
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
i_element = TRANSFER( realIn, 0 )
! Reverse order of 4 bytes in 32 bit INTEGER space:
CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
CALL MVBITS( i_element,  0, 8, i_element_br, 24 )

! Transfer reversed order bytes to 32 bit REAL space (realOut):
realIn = TRANSFER( i_element_br, 0.0 )
END SUBROUTINE ENDIAN4R
!==========================================================
!+
!+
!+
!==========================================================
SUBROUTINE ENDIAN4I(IntIn)
!! Byte Swap 
IMPLICIT NONE
!* input/output
INTEGER(KIND=JPIM), INTENT(INOUT)    :: IntIn
! Local variables
INTEGER                              :: i_element
INTEGER                              :: i_element_br
!================================================
! Transfer 32 bits of realIn to generic 32 bit INTEGER space:
i_element = TRANSFER( IntIn, 0 )
! Reverse order of 4 bytes in 32 bit INTEGER space:
CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
CALL MVBITS( i_element,  0, 8, i_element_br, 24 )

  intIn = i_element_br
END SUBROUTINE ENDIAN4I
!####################################################################





!####################################################################
! file I/O
!-- INQUIRE_FID : inruire unused file FID
!-- NCERROR     : netCDF I/O wrapper
!####################################################################
FUNCTION INQUIRE_FID() RESULT(FID)
IMPLICIT NONE
!* input/output
INTEGER :: FID ! FILE ID
!* local variable
LOGICAL :: I_OPENED ! FILE ID IS ALREADY USED OR NOT?
!================================================
DO FID = 10, 999
  INQUIRE(FID,OPENED=I_OPENED)
  IF ( .NOT. I_OPENED ) RETURN
ENDDO
END FUNCTION INQUIRE_FID
!==========================================================
!+
!+
!+
!==========================================================
#ifdef UseCDF_CMF
SUBROUTINE NCERROR(STATUS,STRING)
!! NETCDF error handling 
USE NETCDF
IMPLICIT NONE
INTEGER,INTENT(IN)                     :: STATUS
CHARACTER(LEN=*),INTENT(IN),OPTIONAL   :: STRING
!================================================
IF ( STATUS /= 0 ) THEN
  WRITE(LOGNAM,*)  TRIM(NF90_STRERROR(STATUS))
  IF( PRESENT(STRING) ) WRITE(LOGNAM,*) TRIM(STRING)
  WRITE(LOGNAM,*) 'PROGRAM STOP ! '
  STOP 10
ENDIF
END SUBROUTINE NCERROR
#endif
!####################################################################

END MODULE CMF_UTILS_MOD
