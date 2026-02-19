! (C) Copyright 1989- Meteo-France.

      SUBROUTINE SUCOM
C
C     READ INPUT, INITIALISE COMMON BLOCKS
C
C                                                                               
C     PEDRO VITERBO         E.C.M.W.F.    03/94.                            
C
#include "param.h"
#include "comddh.h"
#include "comcha.h"
#include "namus.h"
      CHARACTER CCDUM*2,CCD*1
      CHARACTER CTABLE*256
C
C     1. INITIALISE CONSTANTS
C
      API=2.*ASIN(1.)
      CPD=1005.46
      CPS=2.19E+06
      RDS(1)=0.07
      RDS(2)=0.21
      RDS(3)=0.72
      RDS(4)=1.89
      RLV=2.5008E+06
      RLS=2.8345E+06
      RLF=RLS-RLV
      ZCPV=1869.46
      VTMPC2=ZCPV/CPD-1
      G=9.80665
C
C     2. READ INPUTS, EXCEPT NAMELIST
C
C     EXPERIMENT NAME
      OPEN(UNIT=10,FILE='input',FORM='FORMATTED',STATUS='OLD')
      READ (10,"(A4)") CEXP
      READ (10,"(I10)") IDATE
      WRITE(*,"(' EXPERIMENT ',A4)") CEXP
      WRITE(*,"(' DATE ',I10)") IDATE
      CLOSE(UNIT=10)
C     DESCRIPTION TABLE
C  Description table of DDH
C
C  C     Combination parameter (INTEGER).
C        If it is 0, only the parameter referred on the line will be used.
C        If it is N, the present line and the subsequent N lines will be 
C        combined, according to the first character on the DDH variable
C        table (only similar variables - starting by the same character -
C        can be combined).
C        
C  S     Specifies the operaters to apply to the variable (CHARACTER).
C        I multiply by g/DELP(0)
C        H multiply by g/DELP
C        i multiply by g/(CP(0)*DELP(0))      Note CP(0)=CP(Q(0))
C        h multiply by g/(CP*DELP)            Note CP=CP(Q)
C        A multiply by g/TIME         TIME is the total forecast time
C        T multiply by g/(ADELP*TIME) ADELP is the time accumulated
C                                     pressure
C        t multiply by g/(CP*ADELP*TIME)
C        F divide by TIME
C        D vertical derivative
C        d (1/CP)*(vertical derivative)
C  
C  After the operation described by S is used , a constant offset (Offset)
C  and scaling (Gain) are applied in order to have the fields in the
C  plotting units (Units). Plotted vertical fluxes are positive downwards.
C
C  Quantities will be plotted, by default, with contour interval Rint;
C  Value Rmin (Rmax) will be used for the axis, as axis minimum (maximum)
C  value. Shading is controlled by Lsh; if TRUE, shading will be done
C  between values Rmin and Rmax. If Rmin*Rmax < 0 contours are 
C  positive/solid, negative/dashed. If Llist is TRUE a LEVEL_LIST will
C  be defined in the setup routine.
C
C
C PAR  PROC  DDH Variable   CS  Description                   Units       Gain      Offset  Rint     Rmin     Rmax     LshLliLze
C
C AxxxxA4CCxAA13CCCCCCCCCCxxIAxxA28CCCCCCCCCCCCCCCCCCCCCCCCCxxA8CCCCCCCxE9.2CCCCCxxF7.2CCCxxF7.1CCCxxF7.1CCCxxF7.1CCCxxLxxLxxL
C
      CALL GET_ENVIRONMENT_VARIABLE("DDH_DESCRIPTION_TABLE", CTABLE)
      IF (LEN_TRIM(CTABLE).EQ.0) THEN
       CTABLE = '/home/emos/src/ddh/runddh/description_table.dat'
      ENDIF
      OPEN(UNIT=11,FILE=CTABLE,FORM='FORMATTED',STATUS='OLD')
      READ (11,"(A3)") CDPAR(1)
      READ (11,"(A3)") CDPAR(1)
      DO 210 J=1,JPD
c      DO 210 J=1,18
      READ (11,"(A3,2X,A4,X,A1,A13,1X,I2,A1,2X,A28,2X,A8,2X,E9.2,2X,
     *      F7.2,3(2X,F7.1),3(2X,L1))",END=212) 
     *     CDPAR(J),CDPRO(J),CDCOMB(J),CDVARDDH(J),NDCOMB(J),CDTYPE(J),
     *     CDTEXT(J),CDUNITS(J),GAIN(J),OFFSET(J),
     *     RINT(J),RMIN(J),RMAX(J),LSHADE(J),LLIST(JPD),LZERO(JPD)

C
C     LIMITS ALL PAR VARIABLES, EXCEPT SURFACE, TO ONE CHARACTER
C
      IF (CDPAR(J)(1:1).NE.'S') THEN
            CCD=CDPAR(J)(1:1)
            CDPAR(J)=CCD
      ENDIF




 210  CONTINUE
 212  CONTINUE
      ND=J-1
      CLOSE(UNIT=11)
C
C     DEFINE HERE VARIABLES WITH CONTOUR_LEVEL_LIST (LLIST(J)=.TRUE.)
C     (DEFINE RLIST(J,JL) AND NLIST(J))
C
C
C     3. READ NAMELIST NAMUS (GENERAL CONTROL NAMELIST)
C
C      
C     SET DEFAULTS
C
      CPARDUM='%' 
C '%%%' Wall
      CPRODUM='%%%%'
      NLEVTOP=1
      NLEVBOT=-999
      NLVTOT=0
      NSTEPB=0
      NSTEPI=4
      NUDOM(:)=-999
      NLVPTS(:)=-999
      PARAM(:)=CPARDUM
      PROCESS(:)=CPRODUM
      MYTEXT=' '
      AREATYP='zon'
      PLOTTYP='profil'
      RLATS=-90.
      RLATN=90.
      LPRFI=.FALSE.
      LPLOT=.TRUE.
      LARCHIVE=.FALSE.
      llsm=.false.
      lheader=.true.
C      
C     SET DEFAULTS FOR CONTOUR-STYLE PARAMETERS
C
      CPARDUM='%'
C '%%%' Wall
      CPRODUM='%%%%'
      DO 314 J=1,50
      RUINT(J)=-999.
      RUMIN(J)=-999.
      RUMAX(J)=-999.
      LUSHADE(J)=.FALSE.
      LULIST(J)=.FALSE.
      LUZERO(J)=.FALSE.
      UGAIN(J)=-999.
      UOFFSET(J)=-999.
      CUNITS(J)=CPRODUM//CPRODUM
      CUPRO(J)=CPRODUM
 314  CUPAR(J)=CPARDUM
      DO 316 JL=1,15
      DO 316 J=1,50
 316  RULIST(J,JL)=-999.
      LCOLOUR=.FALSE.
      ORNTATION=' '
      CARCSUI=' '
C
C     READ NAMELIST
C
      OPEN(UNIT=12,FILE='namelist',FORM='FORMATTED',STATUS='OLD')
      READ (12,NAMUS)
C      WRITE(*,NAMUS)
      CLOSE(UNIT=12)
C     INCOMPATIBILITIES
      IF (LARCHIVE) LPLOT=.FALSE.
C     SET UP SOME VARIABLES
      IF (CEXP(1:1).EQ.'1'.OR.CEXP(1:1).EQ.'2'.OR.CEXP(1:1).EQ.'3'.OR.
     *  CEXP(1:1).EQ.'4'.OR.CEXP(1:1).EQ.'5'.OR.CEXP(1:1).EQ.'6'.OR.
     *  CEXP(1:1).EQ.'7'.OR.CEXP(1:1).EQ.'8'.OR.CEXP(1:1).EQ.'9') THEN
         READ(CEXP,"(I4)") IEXP
         WRITE(CDU2,"(I4.4)") IEXP
      ELSE
         CDU2=CEXP
      ENDIF
      NSTEP(1)=NSTEPE
      NSTEP(2)=NSTEPB
      WRITE(CDU1,"(I0)") NSTEP(1)
      IF (AREATYP.EQ.'zon') THEN
         CCDUM='ZO'
      ELSE IF (AREATYP.EQ.'lar') THEN
         CCDUM='DL'
      ENDIF
      CLFNAME(1)='DHF'//CCDUM//CDU2//'+'//CDU1
      IF (NSTEP(2).NE.0) THEN
         WRITE(CDU1,"(I0)") NSTEP(2)
         CLFNAME(2)='DHF'//CCDUM//CDU2//'+'//CDU1
      ENDIF
      IF (ORNTATION.EQ.' ') THEN
         IF (PLOTTYP.EQ.'timser') THEN
            ORNTATION='LANDSCAPE'
         ELSE
            ORNTATION='PORTRAIT'
         ENDIF
      ENDIF
C
C     4. ESTABLISH WORK TO BE DONE
C

      IF (LARCHIVE) THEN
         CALL SUARCH
      ELSE
         NDOM=NUDOM(1)
         NUDTOT=1
         DO 410 J=1,15
         IF (PARAM(J).EQ.CPARDUM) GO TO 412
 410     CONTINUE
 412     IPAR=J-1
         DO 414 J=1,15
         IF (PROCESS(J).EQ.CPRODUM) GO TO 416
 414     CONTINUE
 416     IPRO=J-1
C
         NTPL=0
C  
         DO 422 JPAR=1,IPAR
         DO 422 JPRO=1,IPRO
         DO 418 J=1,ND
         IF (PARAM(JPAR).EQ.CDPAR(J).AND.PROCESS(JPRO).EQ.CDPRO(J)) THEN
            NTPL=NTPL+1
            NINDEX(NTPL)=J
            GOTO 420
         ENDIF
 418     CONTINUE
         WRITE(*,"(' COMBINATION OF PARAM ',A1,' AND PROCESS ',A4,
     *      ' NOT AVAILABLE ')") PARAM(JPAR),PROCESS(JPRO)
 420     CONTINUE
 422     CONTINUE
C
         IF (NTPL.EQ.0) THEN
            WRITE(*,*) ' NO WORK TO BE DONE '
            STOP
         ENDIF
C
         DO 424 J=1,NTPL
         JI=NINDEX(J)
         WRITE(*,"(' PARAM ',A3,' AND PROCESS ',A4,
     *      ' SELECTED FOR PLOTTING')") CDPAR(JI),CDPRO(JI)
 424     CONTINUE
C
         DO 426 J=1,15
         IF (NLVPTS(J).EQ.-999) GO TO 428
 426     CONTINUE
 428     NLVTOT=J-1
C
C     5. MODIFY CONTOUR CHARACTEERISTICS, ACCORDING TO USER' S REQUESTS
C
         DO 510 J=1,50
         IF (CUPAR(J).EQ.CPARDUM) GO TO 512
 510     CONTINUE
 512     IPAR=J-1
C
         IF (IPAR.NE.0) THEN
            DO 520 JPAR=1,IPAR
            DO 518 J=1,ND
            IF (CUPAR(JPAR).EQ.CDPAR(J).AND.CUPRO(JPAR).EQ.CDPRO(J))
     *          THEN
               IF (RUINT(JPAR).NE.-999.) THEN
                  RINT(J)=RUINT(JPAR)
               ENDIF
               IF (RUMIN(JPAR).NE.-999.) THEN
                  RMIN(J)=RUMIN(JPAR)
               ENDIF
               IF (RUMAX(JPAR).NE.-999.) THEN
                  RMAX(J)=RUMAX(JPAR)
               ENDIF
               IF (UGAIN(JPAR).NE.-999.) THEN
                  GAIN(J)=UGAIN(JPAR)
               ENDIF
               IF (UOFFSET(JPAR).NE.-999.) THEN
                  OFFSET(J)=UOFFSET(JPAR)
               ENDIF
               IF (.NOT.LUZERO(JPAR)) THEN
                  LZERO(J)=LUZERO(JPAR)
               ENDIF
               IF (.NOT.LUSHADE(JPAR)) THEN
                  LSHADE(J)=LUSHADE(JPAR)
               ENDIF
               IF (CUNITS(JPAR).NE.'%%%%%%%%') THEN
                  CDUNITS(J)=CUNITS(JPAR)
               ENDIF
               IF (.NOT.LULIST(JPAR)) THEN
                  LLIST(J)=LULIST(JPAR)
                  DO 514 JLIS=1,15
                  IF (RULIST(JPAR,JLIS).EQ.-999.) GO TO 516
! sucom.F:287:0: warning: ‘jplis’ may be used uninitialized in this function [-Wmaybe-uninitialized] ! FIXME
C 514              RLIST(J,JLIS)=RULIST(JPAR,JPLIS)
 514              RLIST(J,JLIS)=RULIST(JPAR,JLIS)
 516              NLIST(J)=JLIS-1
               ENDIF
               WRITE(*,"(' CONTOURING OF PARAM ',A1,' AND PROCESS ',A4,
     *            ' WILL BE MADE WITH USER-SPECIFIED STYLE')") 
     *            PARAM(JPAR),PROCESS(JPRO)
               GOTO 520
            ENDIF
 518        CONTINUE
 520        CONTINUE
         ENDIF
      ENDIF
C
      RETURN
      END
