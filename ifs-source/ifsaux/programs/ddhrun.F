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

      PROGRAM DDHRUN
C
C     MASTER PROGRAM TO READ AND PLOT DDH DIAGNOSTICS
C
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    03/94.                            
C
#include "param.h"
#include "comcha.h"
#include "comddh.h"
C
C     1. READ NAMELIST, SET UP PART OF COMMON BLOCKS 
C
      CALL SUCOM
C
C     2. OPEN FILES, FILL REST OF COMMON BLOCKS, SET UP AUXILIARY ARRAYS
C
      CALL OPONE
C
C     FILLS REST OF COMMON BLOCK
C
      if (.not.larchive) CALL SUCOM2(.true.)
C
C     3. WRITE "ALAN BETTS" FILES
C
      IF (LARCHIVE) THEN
         CALL CTLARCH
C
C     4. PLOT/PRINT ZONAL FILES
C
      ELSE IF (AREATYP.EQ.'zon') THEN
         CALL CTLZON
C
C     5. PLOT/PRINT A VERTICAL PROFILE OF LIMITED AREA DOMAINS
C
      ELSE IF (AREATYP.EQ.'lar'.AND.PLOTTYP.EQ.'profil') THEN
         CALL CTLDOMPF
C
C     6. PLOT/PRINT A TIME SERIES OF LIMITED AREA DOMAINS
C
      ELSE IF (AREATYP.EQ.'lar'.AND.PLOTTYP.EQ.'timser') THEN
         CALL CTLDOMTS
C
      ENDIF
C
      STOP
      END
