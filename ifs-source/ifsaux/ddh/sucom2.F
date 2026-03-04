! (C) Copyright 1989- Meteo-France.

      SUBROUTINE SUCOM2(ldverbose)
C
C     FIND RECORD LENGTH, FILL REST OF COMMON BLOCKS
C
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    03/94.                            
C
#include "param.h"
#include "comddh.h"
#include "comdat.h"
#include "comcha.h"
      logical ldverbose
C
      CHARACTER CLRECNM*12
      INTEGER IDATE(11),IDOCFI(17)
      REAL ZDOCNNN(11)
C
C
C     FIND DATE, TIME
C
      IBLOCK=1
      CLRECNM='DATE'
      CALL GETRECI(CLRECNM,IBLOCK,NBLLEN,
     S CLDATA(1),NLENCHBL,NDATA(1,1),NLENIBL,NDATA(1,1),NLENIBL,
     S 11,IDATE,ILEN)
      NYEA=IDATE(1)
      NMON=IDATE(2)
      NDAY=IDATE(3)
      NTIME=IDATE(4)
      NDATE=NDAY+100*NMON+10000*NYEA
C
C     FIND OTHER RELEVANT INFO
C
      IBLOCK=1
      CLRECNM='DOCFICHIER'
      CALL GETRECI(CLRECNM,IBLOCK,NBLLEN,
     S CLDATA(1),NLENCHBL,NDATA(1,1),NLENIBL,NDATA(1,1),NLENIBL,
     S 17,IDOCFI,ILEN)
      NAREAS=IDOCFI(15)
      NLEV=IDOCFI(6)
      IF (NLEVBOT.EQ.-999) NLEVBOT=NLEV
C
C     FIND GEOGRAPHICAL INFO FOR LIMITED AREA DOMAINS
C
      IF (AREATYP.EQ.'lar') THEN
         CLRECNM='DOCD'
         IBLOCK=NDOM
         CALL GETREC(CLRECNM,IBLOCK,NBLLEN,
     S   CLDATA(1),NLENCHBL,DATA(1,1),NLENRBL,NDATA(1,1),NLENIBL,
     S   11,ZDOCNNN,ILEN)
         NTYPE=NINT(ZDOCNNN(11))
         nvpln=NINT(ZDOCNNN(1))
         ZRD2DG=90./ASIN(1.)
         IF (NTYPE.EQ.3) THEN
            RWES=ZRD2DG*ZDOCNNN(3)
            RNOR=ZRD2DG*ASIN(ZDOCNNN(4))
            REAS=ZRD2DG*ZDOCNNN(7)
            RSOU=ZRD2DG*ASIN(ZDOCNNN(8))
            if (ldverbose)
     *      write(*,"('Area found: ',4(1x,f8.2))") rwes,rnor,reas,rsou
         ELSE IF (NTYPE.EQ.4.or.ntype.eq.1) THEN
            RWES=ZRD2DG*ZDOCNNN(3)
            RNOR=ZRD2DG*ASIN(ZDOCNNN(4))
            if (ldverbose) 
     *      write(*,"('Point found: ',2(1x,f8.2))") rwes,rnor
         else IF (NTYPE.EQ.2) THEN
            RWES=ZRD2DG*ZDOCNNN(3)
            RNOR=ZRD2DG*ASIN(ZDOCNNN(4))
            REAS=ZRD2DG*ZDOCNNN(5)
            RSOU=ZRD2DG*ASIN(ZDOCNNN(6))
            RWES2=ZRD2DG*ZDOCNNN(7)
            RNOR2=ZRD2DG*ASIN(ZDOCNNN(8))
            REAS2=ZRD2DG*ZDOCNNN(9)
            RSOU2=ZRD2DG*ASIN(ZDOCNNN(10))
            if (ldverbose) write(*,"('Area found: ',8(1x,f8.2))") 
     *            rwes,rnor,reas,rsou,rwes2,rnor2,reas2,rsou2
         ENDIF
      ENDIF
C
      RETURN
      END
