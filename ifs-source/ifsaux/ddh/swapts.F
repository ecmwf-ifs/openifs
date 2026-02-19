! (C) Copyright 1989- Meteo-France.

      SUBROUTINE SWAPTS(KT)
C
C     SSWAPS THE FIELDS WHEN USING TIME SERIES
C
C
C     INPUT
C     -----
C
C     VAR        DIMENSION   DESCRIPTION
C
C  INTEGER
C
C     KT                     PREVIOUS TIME STEP
C
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    12/94.                            
C
#include "param.h"
#include "comddh.h"
#include "comcha.h"
#include "comdat.h"
C
      CHARACTER*4 CDU1,CDU2
      CHARACTER*2 CCDUM
C
C     COPY TIME STEP 1 INTO 2
C
      CLDATA(2) (1:NLENCHBL)=CLDATA(1) (1:NLENCHBL)
C
      DO 110 J=1,NLENRBL
 110  DATA(J,2)=DATA(J,1)
C
      DO 112 J=1,NLENIBL
 112  NDATA(J,2)=NDATA(J,1)
C
      SECOND(2)=SECOND(1)
C
      IF (CEXP(1:1).EQ.'1'.OR.CEXP(1:1).EQ.'2'.OR.CEXP(1:1).EQ.'3'
     *.OR.CEXP(1:1).EQ.'4'.OR.CEXP(1:1).EQ.'5'.OR.CEXP(1:1).EQ.'6'
     *.OR.CEXP(1:1).EQ.'7'.OR.CEXP(1:1).EQ.'8'.OR.CEXP(1:1).EQ.'9')
     *   THEN
         READ(CEXP,"(I4)") IEXP
         WRITE(CDU2,"(I4.4)") IEXP
      ELSE
         CDU2=CEXP
      ENDIF
C
C     TEMPORARILY FOOLS OPONE
C
      NSTEP(2)=0
C
      NSTEP(1)=KT+NSTEPI
      WRITE(CDU1,"(I0)") NSTEP(1)
      IF (AREATYP.EQ.'zon') THEN
         CCDUM='ZO'
      ELSE IF (AREATYP.EQ.'lar') THEN
         CCDUM='DL'
      ENDIF
      CLFNAME(1)='DHF'//CCDUM//CDU2//'+'//CDU1
      WRITE(CDU1,"(I0)") NSTEP(2)
      CLFNAME(2)='DHF'//CCDUM//CDU2//'+'//CDU1
      CALL OPONE
C
      NSTEP(2)=KT
C
      RETURN
      END
