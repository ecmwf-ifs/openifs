! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE WSTREAM_STRG(ISTREAM,CSTREAM,NENSFNB,NTOTENS,          &
     &                        MARSFCTYPE, KSTREAM, LASTREAM)

! ----------------------------------------------------------------------

!**** *WSTREAM_STRG* -

!     J. BIDLOT     ECMWF    AUGUST 2000.


!*    PURPOSE.
!     --------
!     TO RETURN A 4 CHARACTER STRING WHICH CORRESPONDS TO THE VALUE
!     ISTREAM WHEN DEALING WITH GRIB DATA AND THE FDB.
!     IF NO MATCH IS FOUND, CSTREAM IS SET TO '****'
!     BASED ON ISTREAM, NENSFNB (the ensemble number)
!                   and NTOTENS (total number of ensemble members)
!     FORECAST TYPE WILL ALSO BE RETURNED. 
!     IF ISTREAM IS A WAVE STREAM THEN KSTREAM IS RETURNED WITH THE
!     CORRESPONDING IFS STREAM AND VICE VERSA.
!     LASTREAM IS SET TO TRUE IF ISTREAM IS AN IFS STREAM.

!**   INTERFACE.
!     ----------
!     *CALL WSTREAM_STRG*(ISTREAM,CSTREAM)

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: ISTREAM, NENSFNB, NTOTENS
      INTEGER(KIND=JWIM), INTENT(OUT) :: KSTREAM

      CHARACTER(LEN=2), INTENT(OUT) :: MARSFCTYPE
      CHARACTER(LEN=4), INTENT(OUT) :: CSTREAM
      LOGICAL, INTENT(OUT) :: LASTREAM

! ----------------------------------------------------------------------

      IF(ISTREAM.EQ.1045) THEN
!       DETERMINISTIC WAVE FORECAST
        CSTREAM = 'wave'
        MARSFCTYPE = 'fc'
        KSTREAM = 1025
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1077) THEN
!       WAVE ENSEMBLE forecast HINDCAST STATISTICS
        CSTREAM = 'wehs'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1040
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1125) THEN
!       WAVE ENSEMBLE forecast HINDCAST STATISTICS
        CSTREAM = 'wehs'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1121
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1079) THEN
!       ENSEMBLE FORECAST WAVE HINDCAST
        CSTREAM = 'enwh'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1033
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1080) THEN
!       DETERMINISTIC WAVE MONTHLY MEANS AND CLIMATOLOGY
        CSTREAM = 'wamo'
        MARSFCTYPE = 'fc'
        KSTREAM = 1071
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1081) THEN
!       ENSEMBLE WAVE FORECAST
        CSTREAM = 'waef'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1035
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1123) THEN
!       WAVE EXTENDED ENSEMBLE  FORECAST
        CSTREAM = 'weef'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1122
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1082) THEN
!       WAVE SEASONAL FORECAST
        CSTREAM = 'wasf'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1090
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1083) THEN
!       WAVE MULTI ANALYSIS
        CSTREAM = 'mawv'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1037
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1084) THEN
!       ENSEMBLE WAVE HINDCAST
        CSTREAM = 'ewhc'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1039
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1124) THEN
!       WAVE EXTENDED ENSEMBLE HINDCAST
        CSTREAM = 'weeh'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1120
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1024) THEN
!       WAVE HINDCAST
        CSTREAM = 'wvhc'
        MARSFCTYPE = 'fc'
        KSTREAM = 1085
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1078) THEN
!       Ensemble Forecast Wave Hindcast Overlap 
        CSTREAM = 'ewho'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1032
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1086) THEN
!       WAVE ENSEMBLE FORECAST OVERLAP 
        CSTREAM = 'weov'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1034
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1027) THEN
!       WAVE SHORT CUT-OFF FORECAST
        CSTREAM = 'scwv'
        MARSFCTYPE = 'fc'
        KSTREAM = 1026
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1029) THEN
!       WAVE DELAYED CUT-OFF FORECAST
        CSTREAM = 'dcwv'
        MARSFCTYPE = 'fc'
        KSTREAM = 1028
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1248) THEN
!       WAVE LONG CUT-OFF FORECAST
        CSTREAM = 'lwwv'
        MARSFCTYPE = 'fc'
        KSTREAM = 1247
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1092) THEN
!       WAVE SEASONAL MONTHLY MEANS
        CSTREAM = 'swmm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1091
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1095) THEN
!       WAVE MONTHLY FORECAST (old)
        CSTREAM = 'wamf'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1093
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1096) THEN
!       WAVE MONTHLY FORECAST MEANS (old)
        CSTREAM = 'wmfm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1094
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1203) THEN
!       MONTHLY FORECAST (new)
        CSTREAM = 'mnfw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1200
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1204) THEN
!       MONTHLY FORECAST HINDCAST
        CSTREAM = 'mfhw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1201
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1205) THEN
!       MONTHLY FORECAST ANOMALIES
        CSTREAM = 'mfaw'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1202
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1209) THEN
!       MONTHLY FORECAST MEANS
        CSTREAM = 'mfwm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1206
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1210) THEN
!       MONTHLY FORECAST HINDCAST MEANS
        CSTREAM = 'mhwm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1207
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1211) THEN
!       MONTHLY FORECAST HINDCAST ANOMALIES
        CSTREAM = 'mawm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1208
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1222) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST
        CSTREAM = 'wams'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1220
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1223) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST MONTHLY MEANS
        CSTREAM = 'mswm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1221
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1025) THEN
!       DETERMINISTIC WAVE FORECAST
        CSTREAM = 'oper'
        MARSFCTYPE = 'fc'
        KSTREAM = 1045
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1071) THEN
!       DETERMINISTIC WAVE MONTHLY MEANS AND CLIMATOLOGY
        CSTREAM = 'moda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1080
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1035) THEN
!       ENSEMBLE WAVE FORECAST
        CSTREAM = 'enfo'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1081
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1122) THEN
!       WAVE EXTENDED ENSEMBLE FORECAST
        CSTREAM = 'weef'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1123
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1090) THEN
!       WAVE SEASONAL FORECAST
        CSTREAM = 'seas'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1082
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1033) THEN
!       Ensemble Forecast Wave Hindcasts 
        CSTREAM = 'enfh'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1079
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1040) THEN
!       ENSEMBLE forecast HINDCAST STATISTICS
        CSTREAM = 'efhs'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1077
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1121) THEN
!       Extended ENSEMBLE forecast HINDCAST STATISTICS
        CSTREAM = 'eehs'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1125
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1037) THEN
!       WAVE MULTI ANALYSIS
        CSTREAM = 'maed'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1083
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1039) THEN
!       ENSEMBLE WAVE HINDCAST
        CSTREAM = 'efhc'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1084
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1120) THEN
!        WAVE ENSEMBLE FORECAST HINDCAST
        CSTREAM = 'eefh'
        MARSFCTYPE = 'fc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1124
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1085) THEN
!       WAVE HINDCAST
        CSTREAM = 'wvhc'
        MARSFCTYPE = 'fc'
        KSTREAM = 1024
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1032) THEN
!       WAVE ENSEMBLE FORECAST OVERLAP 
        CSTREAM = 'efho'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1078
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1034) THEN
!       WAVE ENSEMBLE FORECAST OVERLAP 
        CSTREAM = 'efov'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
        KSTREAM = 1086
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1026) THEN
!       WAVE SHORT CUT-OFF FORECAST
        CSTREAM = 'scda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1027
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1028) THEN
!       WAVE DELAYED CUT-OFF FORECAST
        CSTREAM = 'dcda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1029
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1247) THEN
!       WAVE LONG CUT-OFF FORECAST
        CSTREAM = 'lwda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1248
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1091) THEN
!       WAVE SEASONAL MONTHLY MEANS
        CSTREAM = 'sfmm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1092
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1093) THEN
!       WAVE MONTHLY FORECAST (old)
        CSTREAM = 'mofc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1095
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1094) THEN
!       WAVE MONTHLY FORECAST MEANS (old)
        CSTREAM = 'mofm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1096
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1200) THEN
!       MONTHLY FORECAST (new)
        CSTREAM = 'mnfc'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1203
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1201) THEN
!       MONTHLY FORECAST HINDCAST
        CSTREAM = 'mnfh'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1204
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1202) THEN
!       MONTHLY FORECAST ANOMALIES
        CSTREAM = 'mnfa'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1205
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1206) THEN
!       MONTHLY FORECAST MEANS
        CSTREAM = 'mnfn'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1209
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1207) THEN
!       MONTHLY FORECAST HINDCAST MEANS
        CSTREAM = 'mfwm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1210
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1208) THEN
!       MONTHLY FORECAST HINDCAST ANOMALIES
        CSTREAM = 'mfam'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1211
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1220) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST
        CSTREAM = 'mmsf'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
!       KSTREAM = 1222  multi-model seasonal forecasts still use original seasonal wave stream
        KSTREAM = 1082
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1221) THEN
!       MULTI MODEL WAVE SEASONAL FORECAST MONTHLY MEANS
        CSTREAM = 'msmm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'fc'
        ENDIF
        KSTREAM = 1223
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1087) THEN
!       NEW WAVE STAND ALONE STREAM 
        CSTREAM = 'wavm'
        IF(NTOTENS.EQ.0) THEN
          MARSFCTYPE = 'fc'
        ELSEIF(NENSFNB.EQ.0) THEN
          MARSFCTYPE = 'cf'
        ELSE
          MARSFCTYPE = 'pf'
        ENDIF
!       there is no such thing as IFS in stand alone
        KSTREAM = 1087
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1030) THEN
!       ENSEMBLE DATA ASSIMILATION 
        CSTREAM = 'enda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1088
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1249) THEN
!       ENSEMBLE DATA ASSIMILATION 
        CSTREAM = 'elda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1250
        LASTREAM=.TRUE.
      ELSE IF(ISTREAM.EQ.1088) THEN
!       ENSEMBLE WAVE DATA ASSIMILATION 
        CSTREAM = 'ewda'
        MARSFCTYPE = 'fc'
        KSTREAM = 1030
        LASTREAM=.FALSE.
      ELSE IF(ISTREAM.EQ.1250) THEN
!       ENSEMBLE WAVE DATA ASSIMILATION 
        CSTREAM = 'ewla'
        MARSFCTYPE = 'fc'
        KSTREAM = 1249
        LASTREAM=.FALSE.
      ELSE
        CSTREAM='****'
        MARSFCTYPE = 'fc'
        KSTREAM = 0 
        LASTREAM=.TRUE.
      ENDIF

      END SUBROUTINE WSTREAM_STRG
