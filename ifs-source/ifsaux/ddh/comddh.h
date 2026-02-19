! (C) Copyright 1989- Meteo-France.

C
C     CONTAINS ASSORTED INFORMATION ON DDH', GRID, DATE, ...
C
C     VAR       DIMENSION     CONTENTS
C
C  INTEGER
C
C     NDCOMB    JPD           COMBINATION INTEGER PARAMETER
C                              0   THIS LINE IS USED (NO COMBINATION)
C                              N   THE NEXT N LINES ARE USED (ACCORDING
C                                  TO CDCOMB)
C     ND                      ACTUAL SIZE OF DESCRIPTION TABLE
C     NSTEPB                  BEGINING TIME STEP
C     NSTEPE                  END TIME STEP
C     NSTEPI                  INCREMENT STEP
C     NLEVTOP                 TOP LEVEL              (FOR PLOT)
C     NLEVBOT                 BOTTOM LEVEL           (FOR PLOT)
C     NINDEX    JPD           USER CHOSEN PARAMETERS (COMBINATION OF PAR
C                                                     AND PRO)
C     NTPL                    NUMBER OF USER CHOSEN PARAMETERS
C     NDATE                   DATE   (FORMAT yyYYMMDD)
C     NTIME                   TIME   (FORMAT HH)
C     NYEA                    YEAR   (FORMAT yyYY)
C     NMON                    MONTH  (FORMAT MM)
C     NDAY                    DAY    (FORMAT DD)
C     NAREAS                  NUMBER OF LATITUDE BANDS
C     NLEV                    NUMBER OF LEVELS
C     NLIST     JPD           DIMENSION OF CONTOUR_LEVEL_LIST
C     NDOM                    INDEX OF CURRENT LIMITED DOMAIN
C     NUDOM                   INDEX OF LIMITED DOMAIN
C     NUDTOT                  Number OF LIMITED DOMAINs to be processed
C     NVPLN                   INDEX OF vitrual plane
C     NTYPE                   TYPE OF LIMITED DOMAIN (BOX - 3
C                                                     POINT -4)
C     NSTEP     2             NUMBER OF ELAPSED SECONDS (CORRESPONDING TO
C                               CURRENT AND PREVIOUS STEP, RESPECTIVELY)
C     NLVPTS    15            ARRAY WITH LEVELS TO BE PLOTTED AS A TIME
C                               SERIES
C     NLVTOT                  NUMBER OF LEVELS TO BE PLOTTED AS A TIME
C                               SERIES
C
C  REAL
C 
C     API                     PI
C     CPD                     SPECIFIC HEAT OF DRY AIR
C     CPS                     Volumetric Heat capacity of soil
C     RDS                     Depths of the soil layers
C     RLV                     LATENT HEAT OF VAPOURISATION
C     RLS                     LATENT HEAT OF SUBLIMATION
C     RLF                     LATENT HEAT OF FREEZING
C     VTMPC2                  CPV/CPD-1
C     G                       G
C     GAIN      JPD           GAIN FACTOR FOR PLOTTING VARIABLE
C     OFFSET    JPD           OFFSET FACTOR FOR PLOTTING VARIABLE
C     SECOND    2             NUMBER OF ELAPSED SECONDS (CORRESPONDING
C                              TO NSTEP(1) AND NSTEP(2), RESPECTIVELY)
C     RLATS                   LEFT LATITUDE           (FOR PLOT)
C     RLATN                   RIGHT LATITUDE          (FOR PLOT)
C     RINT      JPD           CONTOUR INTERVAL
C     RMIN      JPD           CONTOUR MINIMUM LEVEL, OR MINIMUM VALUE OF
C                              VARIABLE AXIS
C     RMAX      JPD           CONTOUR MINIMUM LEVEL, OR MINIMUM VALUE OF
C                              VARIABLE AXIS
C     RLIST     JPD,15        CONTOUR_LEVEL_LIST
C     RNOR                    NORTHERN BOUNDARY OF DDH LAR DOMAIN 
C     RWES                    WESTERN  BOUNDARY OF DDH LAR DOMAIN 
C     RSOU                    SOUTHERN BOUNDARY OF DDH LAR DOMAIN 
C     REAS                    EASTERN  BOUNDARY OF DDH LAR DOMAIN 
C     RNOR2                   NORTHERN BOUNDARY OF DDH LAR DOMAIN 
C     RWES2                   WESTERN  BOUNDARY OF DDH LAR DOMAIN 
C     RSOU2                   SOUTHERN BOUNDARY OF DDH LAR DOMAIN 
C     REAS2                   EASTERN  BOUNDARY OF DDH LAR DOMAIN 
C
C  LOGICAL
C
C     LCOLOUR                 TRUE FOR COLOUR PLOTS
C     LSHADE    JPD           TRUE FOR SHADING
C     LLIST     JPD           TRUE FOR CONTOUR_LEVEL_LIST
C                             FALSE FOR CONTOUR_INTERVAL
C     LZERO     JPD           TRUE IF ZERO CONTOUR IS TO BE PLOTTED
C     LPRFI                   TRUE IF A FILE NEEDS TO BE CREATED
C     LPLOT                   TRUE IF A PLOT IS NEEDED
C     LARCHIVE                TRUE IF CREATING AN ARCHIVE, ABETTS STYLE
C     LHEADER                 TRUE IF header lines to be printed in the
c                                   archive files
c     LLSM                    TRUE if normalizing the surface fields by
c                                  the percent of land in the area
C
C
      COMMON /COMDDH/NDCOMB,ND,NSTEPB,NSTEPE,NSTEPI,NLEVTOP,NLEVBOT,
     *               NINDEX,NTPL,NDATE,NTIME,NYEA,NMON,NUDOM,NUDTOT,
     *               NDAY,NAREAS,
     *               NLEV,NLIST,NDOM,nvpln,NTYPE,NSTEP,NLVPTS,NLVTOT,
     *               API,CPD,VTMPC2,G,GAIN,OFFSET,SECOND,RLATS,RLATN,
     *               RINT,RMIN,RMAX,RLIST,RNOR,RWES,RSOU,REAS,
     *               RNOR2,RWES2,RSOU2,REAS2,RLV,RLS,
     *               RLF,CPS,RDS,
     *               LCOLOUR,LSHADE,LLIST,LZERO,LPRFI,LPLOT,LARCHIVE,
     *               lheader,llsm
      INTEGER NDCOMB(JPD),NINDEX(JPD),NLIST(JPD),NSTEP(2),NLVPTS(15),
     *        NUDOM(JPVR)
      REAL GAIN(JPD),OFFSET(JPD),SECOND(2),
     *     RINT(JPD),RMIN(JPD),RMAX(JPD),RLIST(JPD,15),RDS(4)
      LOGICAL LCOLOUR,LSHADE(JPD),LLIST(JPD),LZERO(JPD),LPRFI,LPLOT,
     *        LARCHIVE,lheader,llsm
