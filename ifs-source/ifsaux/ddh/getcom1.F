! (C) Copyright 1989- Meteo-France.

      SUBROUTINE GETCOM1(KSTEPP,KINDEX,KON,KVERT,KPTS,PSTEPC,PSTEPP,
     * PON,PPRES,PDELP)
C     
C     GET AND COMBINE TWO TIME EVENTS
C     
C     
C     INPUT
C     -----
C     
C     VAR        DIMENSION   DESCRIPTION
C     
C     INTEGER
C     
C     KSTEPP                 STEP OF THE PREVIOUS EVENT (INDEX)
C     KINDEX                 INDEX OF THE VARIABLE
C     KON                    DIMENSI0N OF PON
C     
C     REAL
C     
C     PSTEPC                 STEP OF THE CURRENT EVENT (SECONDS)
C     PSTEPP                 STEP OF THE PREVIOUS EVENT (SECONDS)
C     
C     OUTPUT
C     ------
C     
C     INTEGER
C     
C     KPTS                   NUMBER OF POINTS OF PON
C     KVERT                  VERTICAL 'DIMENSION' OF PON 
C     
C     REAL
C
C     PON        KON         ARRAY CONTAINING VALUES TO BE PLOTTED
C     PPRES      KON         ARRAY TO CONTAIN PRESSURE
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
      REAL PON(JPZLE),PPRES(JPZLE),PDELP(JPZLE)
C     
      REAL ZPONB(JPZLE),ZPRESB(JPZLE),ZDELPB(JPZLE)
      CHARACTER*12 CLRECNM
      CHARACTER*2 CPARTYPE
C     
C     1. RETRIEVES AND SCALES THE END FIELD
C     
C     
C     1.1 COMBINED AND NON-COMBINED FIELDS
C     
      IF (AREATYP.EQ.'zon') THEN
        CALL GETZON(KINDEX,1,PON,PPRES,PDELP,JPZLE,KPTS,KVERT)
        IIP=KPTS
      ELSE IF (AREATYP.EQ.'lar') THEN
        CALL GETDOM(KINDEX,1,NDOM,PON,PPRES,PDELP,JPZLE,KVERT)
        IIP=KVERT
      ENDIF
C     
C     1.2 COMBINED FIELDS ONLY
C     
      IF (NDCOMB(KINDEX).NE.0) THEN
        CALL GETCOM2(KINDEX,KVERT,1,IIP,PON)
      ENDIF
C     
C     2. RETRIEVES, SCALES THE BEG FIELD, AND COMBINE WITH THE END FIELD
C     
      CLRECNM=CDVARDDH(KINDEX)
      CPARTYPE=CLRECNM(1:2)
      
      IF (CDPRO(KINDEX).NE.'HIST'.AND.CDPRO(KINDEX).NE.'INIT'.AND.
     & CDPAR(KINDEX).NE.'SVE'.and.CPARTYPE.NE.'SV') THEN
        IF (KSTEPP.NE.0) THEN
C     
C     2.1 COMBINED AND NON-COMBINED FIELDS
C     
          IF (AREATYP.EQ.'zon') THEN
            CALL GETZON(KINDEX,2,ZPONB,ZPRESB,ZDELPB,JPZLE,
     *       KPTS,KVERT)
          ELSE IF (AREATYP.EQ.'lar') THEN
            CALL GETDOM(KINDEX,2,NDOM,ZPONB,ZPRESB,ZDELPB,
     *       JPZLE,KVERT)
          ENDIF
C     
C     2.2 COMBINED FIELDS ONLY
C     
          IF (NDCOMB(KINDEX).NE.0) THEN
            CALL GETCOM2(KINDEX,KVERT,2,IIP,ZPONB)
          ENDIF
C     
C     2.3 GETTING THE FINAL FIELD
C     
          DO 230 J=1,IIP
            ZTE=PSTEPC
            ZTB=PSTEPP
            ZTD=ZTE-ZTB
            IF (CDPRO(KINDEX).EQ.'TTOT') THEN
              IF(LARCHIVE) THEN
                PON(J)=(PON(J)-ZPONB(J))/ZTD
              ELSE
                PON(J)=(PON(J)-ZPONB(J))*86400./ZTD    
              ENDIF
            ELSEIF(CDTYPE(KINDEX).EQ.'D'.OR.CDTYPE(KINDEX).EQ.'d'
     &         .OR.CDTYPE(KINDEX).EQ.'L'.OR.CDTYPE(KINDEX).EQ.'l')
     &         THEN
              ZDP=(PDELP(J)-ZDELPB(J))
              PON(J)=(PON(J)-ZPONB(J))/ZDP
            ELSE
              PON(J)=(ZTE/ZTD)*PON(J)-(ZTB/ZTD)*ZPONB(J)
            ENDIF
 230      CONTINUE
        ELSE
          do  j=1,iip
            ZTD=PSTEPC
            IF (CDPRO(KINDEX).EQ.'TTOT') THEN
              IF(LARCHIVE) THEN
                PON(J)=PON(J)/ZTD
              ELSE
                PON(J)=PON(J)*86400./ZTD    
              ENDIF
            ELSEIF(CDTYPE(KINDEX).EQ.'D'.OR.CDTYPE(KINDEX).EQ.'d'
     &         .OR.CDTYPE(KINDEX).EQ.'L'.OR.CDTYPE(KINDEX).EQ.'l')
     &         THEN
              ZDP=PDELP(J)
              PON(J)=PON(J)/ZDP
            ENDIF
          enddo
        ENDIF
      ENDIF
C     
      RETURN
      END
