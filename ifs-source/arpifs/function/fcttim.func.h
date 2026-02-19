! (C) Copyright 1989- ECMWF.
! (C) Copyright 1989- Meteo-France.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

!     ------------------------------------------------------------------

! - Time functions
!   the descriptions are in the annex 1 of the documentation

! TIME

! NDD   : extraxt dd from ccaammdd
! NMM   : extract mm from ccaammdd
! NAA   : extract aa from ccaammdd
! NCCAA : extract ccaa from ccaammdd
! NAMD  : extract aammdd from ccaammdd
! NCENT : return centuary of ccaammdd
! NYEARC: returns year of the centuary from ccaammdd
! NCONSTRUCT_DATE : returns ccaammdd given centuary,year,month and day
! NCTH  : turn seconds into hours
! RTIME : returns the time of the model (in seconds of course!)
! RINCDAY : length of a day (in seconds)

INTEGER(KIND=JPIM) :: NDD,NMM,NCCAA,NAA,NAMD,NCTH,NZZAA,NZZMM,NCENT,NYEARC,&
&NCONSTRUCT_DATE
REAL(KIND=JPRD) :: RJUDAT,RTIME
REAL(KIND=JPRB) :: RINCDAY
INTEGER(KIND=JPIM) :: KGRDAT,KSEC,KAAAA,KMM,KDD,KSS
INTEGER(KIND=JPIM) :: KCENT,KYEARC,KMONTH,KDAY

NDD(KGRDAT)  =MOD(KGRDAT,100)
NMM(KGRDAT)  =MOD((KGRDAT-NDD(KGRDAT))/100,100)
NCCAA(KGRDAT)=KGRDAT/10000
NAA(KGRDAT)=MOD(NCCAA(KGRDAT),100)
NAMD(KGRDAT)=MOD(KGRDAT,1000000)
NCTH(KSEC)=KSEC/3600
NCENT(KGRDAT)=NCCAA(KGRDAT)/100+MIN(NAA(KGRDAT),1)
NYEARC(KGRDAT)=NAA(KGRDAT)+100*(1-MIN(NAA(KGRDAT),1))
NCONSTRUCT_DATE(KCENT,KYEARC,KMONTH,KDAY)=&
&(KCENT-1)*10**6+KYEARC*10**4+KMONTH*10**2+KDAY

NZZAA(KAAAA,KMM)=KAAAA-( (1-SIGN(1,KMM-3))/2 )
NZZMM(KMM)=KMM+6*(1-SIGN(1,KMM-3))
RJUDAT(KAAAA,KMM,KDD)=1720994.5_JPRD + REAL(&
  &2-NZZAA(KAAAA,KMM)/100 + (NZZAA(KAAAA,KMM)/100)/4 &
&+ INT(365.25_JPRD*REAL(NZZAA(KAAAA,KMM),JPRD))&
&+ INT(30.601_JPRD*REAL(NZZMM(KMM)+1,JPRD))&
&+ KDD,JPRD)
RTIME(KAAAA,KMM,KDD,KSS,RINCDAY)=(RJUDAT(KAAAA,KMM,KDD)-2451545._JPRD)&
    &*RINCDAY+REAL(KSS,JPRD)
!    -------------------------------------------------------------
