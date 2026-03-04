! (C) Copyright 1989- ECMWF.
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! 
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction
! 
! (C) Copyright 1989- Meteo-France.
! 

MODULE YOMCHEV

USE PARCMA, ONLY : JPMXRE1, JPMXDE1

IMPLICIT NONE

SAVE

!*     YOMCHEV - CMA EVENTS (CHARACTER VERSION)

!        D. VASILJEVIC   ECMWF     15/09/94

!     REPORT EVENTS PART 1:

!     NAME            TYPE                   MEANING
!     ----            ----                   -------
!     CR1EVENT( 1)      C     ' 1=NO DATA IN THE REPORT'
!     CR1EVENT( 2)      C     ' 2=ALL DATA REJECTED'
!     CR1EVENT( 3)      C     ' 3=BAD REPORTING PRACTICE'
!     CR1EVENT( 4)      C     ' 4=REJECTED DUE TO RDB FLAG'
!     CR1EVENT( 5)      C     ' 5=ACTIVATED DUE TO RDB FLAG' 
!     CR1EVENT( 6)      C     ' 6=ACTIVATED BY WHITELIST' 
!     CR1EVENT( 7)      C     ' 7=HORIZONTAL POSITION OUT OF RANGE'
!     CR1EVENT( 8)      C     ' 8=VERTICAL POSITION OUT OF RANGE'
!     CR1EVENT( 9)      C     ' 9=TIME OUT OF RANGE'
!     CR1EVENT(10)      C     '10=REDUNDANT REPORT'
!     CR1EVENT(11)      C     '11=REPORT OVER LAND'
!     CR1EVENT(12)      C     '12=REPORT OVER SEA'
!     CR1EVENT(13)      C     '13=MISSING STATION ALTITUDE'
!     CR1EVENT(14)      C     '14=MODEL SUR. TOO FAR FROM STAT. ALT.'
!     CR1EVENT(15)      C     '15=REPORT REJECTED THROUGH THE NAMELIST'
!     CR1EVENT(16)      C     '16=FAILED QUALITY CONTROL'

!     DATA EVENTS PART 1:

!     NAME            TYPE               MEANING
!     ----            ----               -------
!     CH1EVENT( 1)      C    ' 1=MISSING VERTICAL COORDINATE'
!     CH1EVENT( 2)      C    ' 2=MISSING OBSERVED VALUE'
!     CH1EVENT( 3)      C    ' 3=MISSING FIRST GUESS VALUE'
!     CH1EVENT( 4)      C    ' 4=REJECTED DUE TO RDB FLAG'
!     CH1EVENT( 5)      C    ' 5=ACTIVATED DUE TO RDB FLAG'
!     CH1EVENT( 6)      C    ' 6=ACTIVATED BY WHITELIST'
!     CH1EVENT( 7)      C    ' 7=BAD REPORTING PRACTICE'
!     CH1EVENT( 8)      C    ' 8=VERTICAL POSITION OUT OF RANGE'
!     CH1EVENT( 9)      C    ' 9=REFERENCE LEVEL POSITION OUT OF RANGE'
!     CH1EVENT(10)      C    '10=TOO BIG FIRST GUESS DEPARTURE'
!     CH1EVENT(11)      C    '11=TOO BIG DEPARTURE IN ASSIMILATION'
!     CH1EVENT(12)      C    '12=TOO BIG OBSERVATION ERROR'
!     CH1EVENT(13)      C    '13=REDUNDANT DATUM'
!     CH1EVENT(14)      C    '14=REDUNDANT LEVEL'
!     CH1EVENT(15)      C    '15=REPORT OVER LAND'
!     CH1EVENT(16)      C    '16=REPORT OVER SEA'
!     CH1EVENT(17)      C    '17=NOT AN ANALYSIS VARIABLE'
!     CH1EVENT(18)      C    '18=DUPLICATED DATUM/LEVEL'
!     CH1EVENT(19)      C    '19=TOO MANY SURFACE DATA/LEVELS'
!     CH1EVENT(20)      C    '20=MULTI LEVEL CHECK'
!     CH1EVENT(21)      C    '21=LEVEL SELECTION'
!     CH1EVENT(22)      C    '22=VERTICAL CONSISTENCY CHECK'
!     CH1EVENT(23)      C    '23=VERTICAL COORDINATE CHANGED FROM Z TO P'
!     CH1EVENT(24)      C    '24=DATUM REJECTED THROUGH THE NAMELIST'
!     CH1EVENT(25)      C    '25=COMBINED FLAGGING'
!     CH1EVENT(26)      C    '26=DATUM REJECTED DUE TO REJECTED REPORT'
!     CH1EVENT(27)      C    '27=VARIATIONAL QC PERFORMED'

CHARACTER (LEN = 64) ::  CR1EVENT(JPMXRE1)
CHARACTER (LEN = 64) ::  CH1EVENT(JPMXDE1)

!-----------------------------------------------------------------------

END MODULE YOMCHEV
