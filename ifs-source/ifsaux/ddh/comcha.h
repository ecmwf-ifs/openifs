! (C) Copyright 1989- Meteo-France.

C
C     CONTAINS DESCRIPTION TABLE CHARACTERISTICS AND OTHER CHARACTER
C     VARIABLES
C
C     VAR       DIMENSION     CONTENTS
C
C  CHARACTER VARIABLES
C
C     CEXP                    EXPERIMENT ID (zabc)
C     CDPAR     JPD           PARAMETER NAME
C     CDPRO     JPD           PROCESS NAME
C     CDWRI     JPD           TEXT TO IDENTIFY VARIABLE (FOR archive)
C                                     cdpar//cdpro
C     CDCOMB    JPD           COMBINATION OPERATOR (IF IT IS NON-BLANK
C                             SPECIFIES THE OPERATION TO BE DONE FOR 
C                             COMBINED VARIABLES - SEE NDCOMB BELOW)
C     AREATYP                 TYPE OF AREA   'ZON'     ZONAL MEANS
C                                            'DOM'     LIMITED DOMAINS
C                                            'GLO'     GLOBAL MEANS
C     CDVARDDH  JPD           NAME OF DDH VARIABLE
C     CDTYPE    JPD           TYPE OF VARIABLE: SPECIFIES THE OPERATOR
C                             TO BE APPLIED
C                              I multiply by g/DELP(0)
C                              H multiply by g/DELP
C                              i multiply by g/(CP(0)*DELP(0))
C                                   Note CP(0)=CP(Q(0))
C                              h multiply by g/(CP*DELP)
C                                   Note CP=CP(Q)
C                              A multiply by g/TIME
C                                   TIME is the total forecast time
C                              T multiply by g/(ADELP*TIME)
C                                   ADELP is the time accumulated pressure
C                              t multiply by g/(CP*ADELP*TIME)
C                              F divide by TIME
C                              D vertical derivative
C                              d (1/CP)*(vertical derivative)
C     CDTEXT    JPD           TEXT TO IDENTIFY VARIABLE (FOR PLOTS)
C     CDUNITS   JPD           UNITS OF VARIABLE (FOR PLOTS)
C     MYTEXT                  USER TEXT
C     CLFNAMEB                NAME OF BEGINING FILE
C     CLFNAMEE                NAME OF END FILE
C     ORNTATION               'PORTRAIT' OR 'LANDSCAPE'
C
C
      COMMON /COMCHA/CEXP,CDPAR,CDPRO,CDWRI,CDCOMB,AREATYP,CDVARDDH,
     *               CDTYPE,
     *               CDTEXT,CDUNITS,MYTEXT,CLFNAMEB,CLFNAMEE,ORNTATION,
     *                CLFNAME,PLOTTYP,CARCSUI,param,process,cpardum,
     *               cprodum
      CHARACTER*4 CEXP
      CHARACTER CDPAR(JPD)*3,CDPRO(JPD)*4,CDWRI(JPD)*7,CDCOMB(JPD)*1,
     *          AREATYP*3,
     *          CDVARDDH(JPD)*13,CDTYPE(JPD)*1,CDTEXT(JPD)*28,
     *          CDUNITS(JPD)*8,MYTEXT*48,CLFNAMEB*14,CLFNAMEE*14,
     *          ORNTATION*9,CLFNAME(2)*14,PLOTTYP*6,CARCSUI*2
