! (C) Copyright 1989- Meteo-France.

      SUBROUTINE GETCOM2(KINDEX,KVERT,KFIL,KPTS,PON)
C     
C     GETS A COMBINED FIELD (ACCORDING TO FIELD TYPE) FOR A SPECIFIED
C     TIME
C     
C     
C     INPUT
C     -----
C     
C     VAR        DIMENSION   DESCRIPTION
C
C     INTEGER
C     
C     KINDEX                 INDEX OF THE VARIABLE
C     KPTS                   DIMENSI0N OF PON
C     KVERT                  VERTICAL 'DIMENSION' OF PON
C     KFIL                   1    FIRST (CURRENT) FILE
C                            2    SECOND (PREVIOUS) FILE
C     
C     INPUT/OUTPUT
C     ------------
C     
C     REAL
C     
C     PON        KPTS        ARRAY CONTAINING VALUES TO BE PLOTTED
C     
C     
C     
C     
C     PEDRO VITERBO         E.C.M.W.F.    12/94.                            
C     
#include "param.h"
#include "comddh.h"
#include "comcha.h"
C     
      REAL PON(JPZLE)
C     
      REAL ZCCC(JPZLE),ZPRES(JPZLE),zdelp(jpzle)
      INTEGER JPDLE, JPTIM
      JPTIM=240
      JPDLE=JPLEVP1*JPTIM
C     
C     1. COMBINATION OF FIELDS
C     
      DO 100 J=1,KPTS
        IF (CDCOMB(KINDEX).EQ.'+') THEN
          ZDELTA=1.
        ELSE IF (CDCOMB(KINDEX).EQ.'-') THEN 
          ZDELTA=-1.
        ENDIF
        PON(J)=ZDELTA*PON(J)
 100  CONTINUE
      DO 104 JC=1,NDCOMB(KINDEX)
        IF (AREATYP.EQ.'zon') THEN
          CALL GETZON(KINDEX+JC,KFIL,ZCCC,ZPRES,zdelp,JPDLE,KPTS,KVERT)
        ELSE IF (AREATYP.EQ.'lar') THEN
          CALL GETDOM(KINDEX+JC,KFIL,NDOM,ZCCC,ZPRES,zdelp,JPZLE,KVERT)
        ENDIF
        DO 102 J=1,KPTS
          IF (CDCOMB(KINDEX+JC).EQ.'+') THEN
            ZDELTA=1.
          ELSE IF (CDCOMB(KINDEX+JC).EQ.'-') THEN 
            ZDELTA=-1.
          ENDIF
          PON(J)=PON(J)+ZDELTA*ZCCC(J)
 102    CONTINUE
 104  CONTINUE
C     
      RETURN
      END
      
