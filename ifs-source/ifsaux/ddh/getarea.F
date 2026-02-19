! (C) Copyright 1989- Meteo-France.

C      SUBROUTINE GETAREA(KINDEX,KST,KAREA,PON,PPRES,ZDELP,KPZV,KVERT)
      SUBROUTINE GETAREA(KINDEX,KST,KAREA,PON,ZDELP,KPZV,KVERT)
C
C     FILLS AND SCALES ONE AREA (ONE VERTICAL)
C
C
C     INPUT
C     -----
C
C     VAR        DIMENSION   DESCRIPTION
C
C  INTEGER
C
C     KINDEX                 INDEX OF THE VARIABLE TO BE PLOTTED
C     KST                    1 FOR FIRST STEP (NSTEPB)
C                            2 FOR SECOND STEP (NSTEPE)
C     KPZV                   DIMENSION OF PON
C     KAREA                  AREA IDENTIFIER
C
C
C     OUTPUT
C     ------
C
C  INTEGER
C
C     KVERT                  NUMBER OF VALUES OF PON 
C
C  REAL
C
C     PON        KPZV        ARRAY TO CONTAIN VALUES
C     PPRES      KPZV        ARRAY TO CONTAIN PRESSURE
C
C
C
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    07/94.                            
#include "param.h"
C
      REAL PON(KPZV),ZDELP(KPZV)
C      REAL PON(KPZV),PPRES(KPZV),ZDELP(KPZV)
#include "comddh.h"
#include "comcha.h"
#include "comdat.h"
#include "namus.h"
C
      CHARACTER*12 CLRECNM
      CHARACTER*1 CSCALE,CPARTYPE
      REAL ZQ(JPLEVP1),ZDUM(JPLEVP1),ZDUMA(JPLEVP1)
      integer idum
      ZEPS=1.E-8

C
C     1. FILLS THE FIELD IN THE RIGHT ORDER
C
      CLRECNM=CDVARDDH(KINDEX)
c      write(0,'(1x,a)') clrecnm
      CALL GETREC(CLRECNM,KAREA,NBLLEN,
     S CLDATA(KST),NLENCHBL,DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,
     S JPLEVP1,ZDUM,KVERT)
c       print*,(zdum(j),j=1,kvert)
c       print*,CDTYPE(KINDEX)
c      print*,'getarea',CLRECNM
C
C     2. PRIMARY (FIRST) "SCALING" - APPLYING A PRE-SPECIFED OPERATOR
C
C
C     2.1 TYPE 'I'      INITIAL FIELDS, NON THERMO
C
      IF (CDTYPE(KINDEX).EQ.'I') THEN
         CLRECNM='VPP0'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         DO 210 J=1,KVERT
         PON(J)=(1./ZDELP(J))*ZDUM(J)
 210     CONTINUE
C
C     2.2 TYPE 'H'      HISTORY FIELDS, NON THERMO
C
      ELSE IF (CDTYPE(KINDEX).EQ.'H') THEN
         CLRECNM='VPP1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         DO 220 J=1,KVERT
         PON(J)=(1./ZDELP(J))*ZDUM(J)
 220     CONTINUE
C
C     2.3 TYPE 'i'      INITIAL FIELDS, THERMO
C
      ELSE IF (CDTYPE(KINDEX).EQ.'i') THEN
         CLRECNM='VPP0'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV0'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 230 J=1,KVERT
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZCP=CPD*(1.+VTMPC2*ZQ(J))
         PON(J)=(1./(ZCP*ZDELP(J)))*ZDUM(J)
 230     CONTINUE
C
C     2.4 TYPE 'h'      HISTORY FIELDS, THERMO
C
      ELSE IF (CDTYPE(KINDEX).EQ.'h') THEN
         CLRECNM='VPP1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 240 J=1,KVERT
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZCP=CPD*(1.+VTMPC2*ZQ(J))
         PON(J)=(1./(ZCP*ZDELP(J)))*ZDUM(J)
 240     CONTINUE
C
C     2.5 TYPE 'A'      TIME-ACCUMULATED PRESSURE 
C
      ELSE IF (CDTYPE(KINDEX).EQ.'A') THEN
         DO 250 J=1,KVERT
         PON(J)=(G/SECOND(KST))*ZDUM(J)
 250     CONTINUE
C
C     2.6 TYPE 'T'      NON-DERIVED TENDENCIES, NON THERMO 
C
      ELSE IF (CDTYPE(KINDEX).EQ.'T') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         DO 260 J=1,KVERT
         PON(J)=ZDUM(J)/ZDELP(J)
 260     CONTINUE
C
C     2.7 TYPE 't'     NON-DERIVED TENDENCIES, THERMO  
C
      ELSE IF (CDTYPE(KINDEX).EQ.'t') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 270 J=1,KVERT
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZCP=CPD*(1.+VTMPC2*ZQ(J))
         PON(J)=ZDUM(J)/(ZCP*ZDELP(J))
 270     CONTINUE
C
C     2.8 TYPE 'F'     FLUXES  
C
      ELSE IF (CDTYPE(KINDEX).EQ.'F') THEN
         DO 280 J=1,KVERT
         PON(J)=ZDUM(J)/(SECOND(KST))
 280     CONTINUE
C
C     2.9 TYPE 'D'     TENDENCIES DERIVED FROM FLUXES, NON THERMO  
C
      ELSE IF (CDTYPE(KINDEX).EQ.'D') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         DO 290 J=1,KVERT-1
         PON(J)=-(ZDUM(J+1)-ZDUM(J))
 290     CONTINUE
         KVERT=KVERT-1
C
C     2.10 TYPE 'd'    TENDENCIES DERIVED FROM FLUXES, THERMO
C
      ELSE IF (CDTYPE(KINDEX).EQ.'d') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 2100 J=1,KVERT-1
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZCP=CPD*(1.+VTMPC2*ZQ(J))
         PON(J)=-(ZDUM(J+1)-ZDUM(J))/ZCP
 2100    CONTINUE
         KVERT=KVERT-1
C
C     2.11 TYPE 'L'    T-TENDENCIES DERIVED FROM Q-FLUXES, Water
C
      ELSE IF (CDTYPE(KINDEX).EQ.'L') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 2110 J=1,KVERT-1
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZLDCP=RLV/(CPD*(1.+VTMPC2*ZQ(J)))
         PON(J)=(ZDUM(J+1)-ZDUM(J))*ZLDCP
 2110    CONTINUE
         KVERT=KVERT-1
C
C     2.12 TYPE 'l'   T-TENDENCIES DERIVED FROM Q-FLUXES, ICE
C
      ELSE IF (CDTYPE(KINDEX).EQ.'l') THEN
         CLRECNM='PPP'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZDELP,ILEN)
         CLRECNM='VQV1'
         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
     S               ZQ,ILEN)
         DO 2120 J=1,KVERT-1
         ZQ(J)=(1./ZDELP(J))*ZQ(J)
         ZLDCP=RLS/(CPD*(1.+VTMPC2*ZQ(J)))
         PON(J)=(ZDUM(J+1)-ZDUM(J))*ZLDCP
 2120    CONTINUE
         KVERT=KVERT-1
C
C     2.13 TYPE 'P' PRESSURE THICKNESSES TO BE CONVERTED TO PRESSURE
C         
      ELSE IF (CDTYPE(KINDEX).EQ.'P'.and.cdpro(kindex).ne.'ACCU') THEN
         PTOP=0.
         DO 2130 J=1,KVERT
         PON(J)=PTOP+(G*ZDUM(J)/2.)
         PTOP=PTOP+(G*ZDUM(J))
 2130    CONTINUE
C
C     2.14 TYPE 'S' SURFACE FIELDS
C         
      ELSE IF (CDTYPE(KINDEX).EQ.'S') THEN
         if (llsm.and.zdum(1) > 0.) then
!  Normalize by land fraction in the area, if non-zero
           zdum(2:6)=zdum(2:6)/zdum(1)
         endif
         IF(CDPAR(KINDEX).EQ.'ST1') THEN
            PON(1)=ZDUM(2)/(CPS*RDS(1))
         ELSEIF(CDPAR(KINDEX).EQ.'STD') THEN
            PON(1)=ZDUM(3)/(CPS*(RDS(2)+RDS(3)+RDS(4)))
         ELSEIF(CDPAR(KINDEX).EQ.'SW1') THEN
            PON(1)=ZDUM(4)/(1000.*RDS(1))
         ELSEIF(CDPAR(KINDEX).EQ.'SWD') THEN
            PON(1)=ZDUM(5)/(1000.*(RDS(2)+RDS(3)+RDS(4)))
         ELSEIF(CDPAR(KINDEX).EQ.'SSN') THEN
            PON(1)=ZDUM(6)
         ELSEIF(CDPAR(KINDEX).EQ.'S') THEN
           PTOP=0.
           DO 2135 J=1,KVERT
             PTOP=PTOP+(G*ZDUM(J))
 2135      CONTINUE
            pon(1)=ptop
         ELSE
            PON(1)=-9999.
         ENDIF
         KVERT=1
C
C     2.15 TYPE 's' SURFACE FLUXES
C         
      ELSE IF (CDTYPE(KINDEX).EQ.'s') THEN
         if (llsm) then
           CLRECNM='S V0'
           CALL GETREC(CLRECNM,KAREA,NBLLEN,
     S  CLDATA(KST),NLENCHBL,DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,
     S     JPLEVP1,ZDUMA,KVERT)
           if (zduma(1) > 0.) then
!  Normalize by land fraction in the area, if non-zero
              zdum(1:17)=zdum(1:17)/zduma(1)
           endif
         endif
         IF(CDPRO(KINDEX).EQ.'FSWR') THEN
            PON(1)=ZDUM(1)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FLWR') THEN
            PON(1)=ZDUM(2)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FLHW') THEN
            PON(1)=ZDUM(3)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FLHN') THEN
            PON(1)=ZDUM(4)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FSH') THEN
            PON(1)=ZDUM(5)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FH1D') THEN
            PON(1)=ZDUM(6)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FMLT'.AND.CDPAR(KINDEX).EQ.'ST1') THEN
            PON(1)=ZDUM(7)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FLSR') THEN
            PON(1)=ZDUM(8)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FLSN') THEN
            PON(1)=ZDUM(9)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FCVR') THEN
            PON(1)=ZDUM(10)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FCVN') THEN
            PON(1)=ZDUM(11)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FEW') THEN
            PON(1)=ZDUM(12)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FEN') THEN
            PON(1)=ZDUM(13)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FW1D') THEN
            PON(1)=ZDUM(14)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FMLT'.AND.CDPAR(KINDEX).EQ.'SW1') THEN
            PON(1)=ZDUM(15)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FRO'.AND.CDPAR(KINDEX).EQ.'SW1') THEN
            PON(1)=ZDUM(16)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FRO'.AND.CDPAR(KINDEX).EQ.'SWD') THEN
            PON(1)=ZDUM(17)/SECOND(KST)
         ELSE
            PON(1)=-9999.
         ENDIF
         KVERT=1
C
C     2.16 TYPE 'E' EXTRA SURFACE FIELDS
C         
       ELSE IF (CDTYPE(KINDEX).EQ.'E') THEN
         IF(CDPRO(KINDEX).EQ.'T2M') THEN
           PON(1)=ZDUM(1)
         ELSEIF(CDPRO(KINDEX).EQ.'Q2M') THEN
           PON(1)=ZDUM(2)
         ELSEIF(CDPRO(KINDEX).EQ.'U10M') THEN
           PON(1)=ZDUM(3)
         ELSEIF(CDPRO(KINDEX).EQ.'V10M') THEN
           PON(1)=ZDUM(4)
         ELSEIF(CDPRO(KINDEX).EQ.'TSK') THEN
           PON(1)=ZDUM(5)
         ELSEIF(CDPRO(KINDEX).EQ.'FSR') THEN
           PON(1)=ZDUM(6)
         ELSEIF(CDPRO(KINDEX).EQ.'FSRH') THEN
           PON(1)=ZDUM(7)
         ELSEIF(CDPRO(KINDEX).EQ.'FAL') THEN
           PON(1)=ZDUM(8)
         ELSEIF(CDPRO(KINDEX).EQ.'BLH') THEN
           PON(1)=ZDUM(9)
         ENDIF
         KVERT=1
C
C     2.17 TYPE 'e' EXTRA SURFACE FLUXES
C         
       ELSE IF (CDTYPE(KINDEX).EQ.'e') THEN
         IF(CDPRO(KINDEX).EQ.'FSRD') THEN
           PON(1)=ZDUM(1)/SECOND(KST)
         ELSEIF(CDPRO(KINDEX).EQ.'FTRD') THEN
           PON(1)=ZDUM(2)/SECOND(KST)
         ENDIF
         KVERT=1


C     
       ELSE IF (CDTYPE(KINDEX).EQ.'N') THEN

c first crack parameter name to find out about type of parameter 
c (flux or variable) and possible scaling (none, land- or sea-ice
c fraction)

         CLRECNM=CDVARDDH(KINDEX)
         CPARTYPE=CLRECNM(2:2)
         CSCALE=CLRECNM(3:3)

c now scale according to parameter type and scaling factor
         
         IF(CSCALE == 'L') THEN
           CLRECNM='SVGFSLF'
           idum = 1
           CALL GETREC(CLRECNM,KAREA,NBLLEN,
     S      CLDATA(KST),NLENCHBL,DATA(1,KST),NLENRBL,
     &      NDATA(1,KST),NLENIBL,JPLEVP1,ZDUMA,idum)
           IF(CPARTYPE == 'V') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)/MAX(ZDUMA(1),ZEPS)
             ENDDO
           ELSEIF(CPARTYPE == 'F') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)/(SECOND(KST)*MAX(ZDUMA(1),ZEPS))
             ENDDO
           ELSE
             WRITE(0,*) ' ILLEGAL PARAMETER TYPE ',CPARTYPE
             CALL ABORT
           ENDIF
         ELSEIF(CSCALE == 'I') THEN
           CLRECNM='SVGFSIF'
           idum = 1        
           CALL GETREC(CLRECNM,KAREA,NBLLEN,
     S      CLDATA(KST),NLENCHBL,DATA(1,KST),NLENRBL,
     &      NDATA(1,KST),NLENIBL,JPLEVP1,ZDUMA,idum)
           IF(CPARTYPE == 'V') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)/MAX(ZDUMA(1),ZEPS)
             ENDDO
           ELSEIF(CPARTYPE == 'F') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)/(SECOND(KST)*MAX(ZDUMA(1),ZEPS))
             ENDDO
           ELSE
             WRITE(0,*) ' ILLEGAL PARAMETER TYPE ',CPARTYPE
             CALL ABORT
           ENDIF
         ELSEIF(CSCALE == 'G') THEN
           IF(CPARTYPE == 'V') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)
             ENDDO
           ELSEIF(CPARTYPE == 'F') THEN
             DO J=1,KVERT
               PON(J)=ZDUM(J)/SECOND(KST)
             ENDDO
           ELSE
             WRITE(0,*) ' ILLEGAL PARAMETER TYPE ',CPARTYPE
             CALL ABORT
           ENDIF
         ELSE
           WRITE(0,*) ' ILLEGAL SCALING TYPE ',CSCALE
           CALL ABORT
         ENDIF
           
C No scaling required
      ELSEIF (CDTYPE(KINDEX).EQ.'Z') THEN
             DO J=1,KVERT
c               write(*,*) 'pon = ',pon(j)
               PON(J)=ZDUM(J)
             ENDDO
       ENDIF
     
C
C     3. GET PRESSURE ARRAY
C
C
C     3.1 TENDENCIES (FULL-LEVELS) FIELDS, INITIAL VALUES
C
c      IF (CDTYPE(KINDEX).EQ.'I'.OR.CDTYPE(KINDEX).EQ.'i') THEN
c         CLRECNM='VPP0'
c         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
c     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
c     S               ZDELP,ILEN)
c         PPRES(1)=0.5*ZDELP(1)
c         ZPPP=ZDELP(1)
c         DO 310 J=2,KVERT
c         PPRES(J)=ZPPP+0.5*ZDELP(J)
c         ZPPP=ZPPP+ZDELP(J)
c 310     CONTINUE
C
C     3.2 TENDENCIES (FULL-LEVELS) FIELDS, INITIAL VALUES
C
c      ELSE IF (CDTYPE(KINDEX).EQ.'H'.OR.CDTYPE(KINDEX).EQ.'h'.OR.
c     *         CDTYPE(KINDEX).EQ.'T'.OR.CDTYPE(KINDEX).EQ.'t'.OR.
c     *         CDTYPE(KINDEX).EQ.'D'.OR.CDTYPE(KINDEX).EQ.'d') THEN
c         CLRECNM='VPP1'
c         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
c     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
c     S               ZDELP,ILEN)
c         PPRES(1)=0.5*ZDELP(1)
c         ZPPP=ZDELP(1)
c         DO 320 J=2,KVERT
c         PPRES(J)=ZPPP+0.5*ZDELP(J)
c         ZPPP=ZPPP+ZDELP(J)
c 320     CONTINUE
C
C     3.3 FLUXES
C
c      ELSE IF (CDTYPE(KINDEX).EQ.'F') THEN
c         CLRECNM='VPP1'
c         CALL GETREC(CLRECNM,KAREA,NBLLEN,CLDATA(KST),NLENCHBL,
c     S               DATA(1,KST),NLENRBL,NDATA(1,KST),NLENIBL,JPLEVP1,
c     S               ZDELP,ILEN)
c         PPRES(1)=0.
c         DO 330 J=2,KVERT
c         PPRES(J)=PPRES(J-1)+ZDELP(J)
c 330     CONTINUE
c      ENDIF
C
      RETURN
      END
