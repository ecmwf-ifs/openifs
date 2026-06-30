SUBROUTINE UPDTIM1S(KSTEP,PTDT,PTSTEP)


USE YOMLUN1S , ONLY : NULOUT
USE YOMLOG1S , ONLY : CFFORC
USE YOMRIP   , ONLY : NSSSSS   ,NSTADD   ,NSTASS   ,RTIMST   , &
     &            RSTATI   ,RTIMTR   ,RHGMT    ,REQTIM   ,RSOVR    ,&
     &            RDEASO   ,RDECLI   ,RWSOVR   ,RIP0     ,RCODEC   ,&
     &            RSIDEC   ,RCOVSR   ,RSIVSR   ,RTDT
USE YOMCST   , ONLY : RPI      ,RDAY     ,REA      ,RI0, REPSM
USE YOERIP   , ONLY : RIP0M    ,RCODECM  ,RSIDECM  ,RCOVSRM  ,&
     &            RSIVSRM
USE YOMDYN1S , ONLY : NSTEP    ,TSTEP
USE YOMCC1S  , ONLY : VCALB ,VCLAIL   ,VCLAIH,VCFWET

USE YOMGPD1S , ONLY : VFALBF   ,&
     &                VFALUVP,VFALUVD,VFALNIP,VFALNID , &
     &                VFALUVI,VFALUVV,VFALUVG, &
     &                VFALNII,VFALNIV,VFALNIG, &
     &                VFLAIL,VFLAIH ,VFFWET,VFTVL,VFTVH

USE YOMDPHY  , ONLY : NPOI
USE YOEPHY, ONLY: LECTESSEL, LECLIM10D
USE YOMCT01S , ONLY : NSTART

#ifdef DOC
! (C) Copyright 1995- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!**** *UPDTIM1S* - UPDATE TIME OF THE ONE COLUMN SURFACE MODEL 

!     Purpose.
!     --------
!     UPDATE TIME OF THE ONE COLUMN SURFACE MODEL 

!**   Interface.
!     ----------
!        *CALL* *UPDTIM1S(KSTEP,PTDT,PTSTEP)

!        Explicit arguments :
!        --------------------
!        KSTEP : TIME STEP INDEX
!        PTDT  : TIME STEP LEAPFROG
!        PTSTEP: TIME STEP

!        Implicit arguments :
!        --------------------
!        YOMRIP

!     Method.
!     -------
!        See documentation

!     Externals.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the 
!        one column surface model 

!     Author.
!     -------
!        Jean-Francois Mahfouf and Pedro Viterbo  *ECMWF*

!     Modifications.
!     --------------
!        Original : 95-03-22

!     ------------------------------------------------------------------
#endif

USE PARKIND1  ,ONLY : JPIM     ,JPRB , JPRD
USE YOMHOOK   ,ONLY : LHOOK    ,DR_HOOK, JPHOOK

IMPLICIT NONE

!* Arguments
INTEGER(KIND=JPIM) :: KSTEP
REAL(KIND=JPRB) :: PTDT,PTSTEP

!* Local variables
INTEGER(KIND=JPIM) :: ITIME,IZT,IPR,NRADFR,ISTADD,ISTASS,IYMD,IHM,IDD,&
     &      IMM,IYYYY,IHH,ISS,IMT1,IMT2,IYT1,IYT2,JL,IMT11,IMT12,IDD1,IDD2
REAL(KIND=JPRD) :: ZTETA,ZSTATI,ZHGMT,ZDEASOM,ZDECLIM,ZEQTIMM,ZSOVRM,&
     &      ZWSOVRM,ZJUL,ZTIMTR,ZT1,ZT2,ZT,ZWEI1,ZWEI2
!          ,ZRVCOV(0:20) !original CTESSEL (0:7)

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "fctast.h"
#include "fcttim.h"
#include "surf_inq.h"

IF (LHOOK) CALL DR_HOOK('UPDTIM1S',0,ZHOOK_HANDLE)
 

ITIME=NINT(PTSTEP)
IZT=NINT(REAL(PTSTEP,KIND=JPRD)*(REAL(KSTEP,KIND=JPRD)+0.5_JPRD))
!CBH      IZT=NINT(PTSTEP*(REAL(KSTEP)))
RSTATI=REAL(IZT,KIND=JPRD)
NSTADD=IZT/NINT(RDAY)
NSTASS=MOD(IZT,NINT(RDAY))
RTIMTR=RTIMST+RSTATI
IPR=0
IF(IPR.EQ.1)THEN
  WRITE(UNIT=NULOUT,FMT='(1X,'' TIME OF THE MODEL '',E20.14,&
     & '' TIME SINCE START '',E20.14)') RTIMTR,RSTATI
ENDIF
RHGMT=REAL(MOD(NINT(RSTATI)+NSSSSS,NINT(RDAY)),KIND=JPRD)

ZTETA=RTETA(RTIMTR)
RDEASO=RRS(ZTETA)
RDECLI=RDS(ZTETA)
REQTIM=RET(ZTETA)
RSOVR =REQTIM+RHGMT
RWSOVR=RSOVR*2._JPRD*REAL(RPI/RDAY,KIND=JPRD)
RIP0=RI0*REA*REA/(RDEASO*RDEASO)

RCODEC=COS(RDECLI)
RSIDEC=SIN(RDECLI)

RCOVSR=COS(RWSOVR)
RSIVSR=SIN(RWSOVR)

RTDT=PTDT


!          2.   PARAMETERS FOR ECMWF-STYLE INTERMITTENT RADIATION 
!               -------------------------------------------------

NRADFR=1
ITIME=NINT( TSTEP)
IZT=NINT( REAL(TSTEP,KIND=JPRD)*(REAL(KSTEP,KIND=JPRD)+0.5_JPRD))
!CBH      IZT=NINT( TSTEP*(REAL(KSTEP)))
ZSTATI=REAL(IZT,KIND=JPRD)+REAL(0.5_JPRD*NRADFR*ITIME,KIND=JPRD)
ISTADD=IZT/NINT(RDAY)
ISTASS=MOD(IZT,NINT(RDAY))
ZTIMTR=RTIMST+ZSTATI
ZHGMT=REAL(MOD(NINT(ZSTATI)+NSSSSS,NINT(RDAY)),KIND=JPRD)

ZTETA=RTETA(ZTIMTR)
ZDEASOM=RRS(ZTETA)
ZDECLIM=RDS(ZTETA)
ZEQTIMM=RET(ZTETA)
ZSOVRM =ZEQTIMM+ZHGMT                                 
ZWSOVRM=ZSOVRM*2._JPRD*REAL(RPI/RDAY,KIND=JPRD)
RIP0M=RI0*REA*REA/(ZDEASOM*ZDEASOM)

RCODECM=COS(ZDECLIM)
RSIDECM=SIN(ZDECLIM)

RCOVSRM=COS(ZWSOVRM)
RSIVSRM=SIN(ZWSOVRM)


!          2.   MODIFY SEASONALLY VARYING FIELDS 
!               --------------------------------

call dattim(zjul,iymd,ihm)
if (nstep == 0 .OR. ihm == 0000 .or. nstep==nstart) then
!  Update albedo at 00 GMT, or define it for the first time step
  idd=ndd(iymd)
  imm=nmm(iymd)
  iyyyy=nccaa(iymd)
  ihh=ihm/100
  iss=60*ihh+60*mod(ihm,100)


IF (LECLIM10D) THEN !LECLIM10D TRUE
  if (idd >= 5 .and. idd < 15) then
      idd1=5
      idd2=15
      imt1=imm
      imt2=imm
      iyt1=iyyyy
      iyt2=iyyyy
      zt1=RTIME(iyt1,imt1,5,0)
      zt2=RTIME(iyt2,imt2,15,0)
      imt11=(imm-1)*3+2
      imt12=(imm-1)*3+3

  else if (idd >= 15 .and. idd < 25) then
      idd1=15
      idd2=25
      imt1=imm
      imt2=imm
      iyt1=iyyyy
      iyt2=iyyyy
      zt1=RTIME(iyt1,imt1,15,0)
      zt2=RTIME(iyt2,imt2,25,0)
      imt11=(imm-1)*3+3
      imt12=(imm-1)*3+4

  else if (idd < 5 ) then
      idd1=25
      idd2=5
      imt1=1+mod(imm+10,12)
      imt2=imm
      if(imt1 == 12) then
        iyt1=iyyyy-1
      else
        iyt1=iyyyy
      endif 
      iyt2=iyyyy

      zt1=RTIME(iyt1,imt1,25,0)
      zt2=RTIME(iyt2,imt2,5,0)

      imt11=(imm-1)*3+1
      imt12=(imm-1)*3+2

  else if (idd >= 25 ) then
      idd1=25
      idd2=5
      imt1=imm
      imt2=1+mod(imm,12)
      iyt1=iyyyy
      if(imt2 == 1) then
        iyt2=iyt1+1
      else
        iyt2=iyt1
      endif
      zt1=RTIME(iyt1,imt1,25,0)
      zt2=RTIME(iyt2,imt2,5,0)

      imt11=(imm-1)*3+4
      imt12=(imm-1)*3+5


  endif
ELSE ! LECLIM10D FALSE

   if (idd >= 15) then
      imt1=imm
      imt2=1+mod(imm,12)
      iyt1=iyyyy
     if(imt2 == 1) then
      iyt2=iyt1+1
     else
      iyt2=iyt1
     endif
  else
    imt1=1+mod(imm+10,12)
    imt2=imm
    if(imt1 == 12) then
      iyt1=iyyyy-1
    else
      iyt1=iyyyy
    endif
    iyt2=iyyyy
  endif
  zt1=RTIME(iyt1,imt1,15,0)
  zt2=RTIME(iyt2,imt2,15,0)

  imt11=imt1
  imt12=imt2
ENDIF  ! LECLIM10D


! zt=RTIME(iyyyy,imm,idd,iss)
zt=RTIME(iyyyy,imm,idd,0)  !updated assuming we're at 00UTC 
zwei1=(zt2-zt)/(zt2-zt1)
zwei2=1.-zwei1


vfalbf(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALUVP(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALUVD(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALNIP(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALNID(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFLAIL(:)=zwei1*VCLAIL(:,imt11)+zwei2*VCLAIL(:,imt12)
VFLAIH(:)=zwei1*VCLAIH(:,imt11)+zwei2*VCLAIH(:,imt12)
VFFWET(:)=zwei1*VCFWET(:,imt11)+zwei2*VCFWET(:,imt12)

! 6-component MODIS albedo: to use an albedo independent of solar
! zenith angle, set only the isotropic component
VFALUVI(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALUVV(:)=0.0_jprb
VFALUVG(:)=0.0_jprb
VFALNII(:)=zwei1*vcalb(:,imt11)+zwei2*vcalb(:,imt12)
VFALNIV(:)=0.0_jprb
VFALNIG(:)=0.0_jprb

!IF (LECTESSEL) THEN
! update LAIL LAIH BCL BCH 
!  CALL SURF_INQ(PRVCOV=ZRVCOV)
!
!    ! crop biome cover is dependant on LAI (not used in this version, rcov is kept cte according to lookup table)
!  DO JL=1,NPOI 
!    IF (VFTVL(JL)>=6) THEN
!      VFBCL(JL)=1._JPRB-EXP(-0.6_JPRB*VFLAIL(JL))
!   ELSE
!      VFBCL(JL) =ZRVCOV(INT(VFTVL(JL)))
!   ENDIF 
!   !high vegetation
!   VFBCH(JL) = ZRVCOV(INT(VFTVH(JL)))
!  END DO
! ENDIF

endif

IF (LHOOK) CALL DR_HOOK('UPDTIM1S',1,ZHOOK_HANDLE)

RETURN
END SUBROUTINE UPDTIM1S
