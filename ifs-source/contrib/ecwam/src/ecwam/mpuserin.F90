! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

      SUBROUTINE MPUSERIN

! ----------------------------------------------------------------------

!**** *MPUSERIN* - ROUTINE TO READ NAMELIST INPUTS.

!     J. BIDLOT      ECMWF          JULY    1996
!     B. HANSEN      ECMWF          APRIL   1997
!     B. HANSEN      ECMWF          JANUARY 1998 NAMELIST INPUT.
!     S. ABDALLA     ECMWF          OCTOBER 2001 LGUST & LADEN ADDED

!*    PURPOSE.
!     --------

!**   INTERFACE.
!     ----------

!      *CALL* *MPUSERIN

!     METHOD.
!     -------

!     SEE USERIN

!     EXTERNALS.
!     ----------

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWALTAS , ONLY : NUMALT   ,IBUFRSAT  ,ALTSDTHRSH,ALTBGTHRSH, &
     &            HSALTCUT, LALTGRDOUT, LALTPAS,                        &
     &            XKAPPA2  ,HSCOEFCOR,HSCONSCOR ,LALTCOR   ,LALTLRGR,   &
     &            LODBRALT ,CSATNAME
      USE YOWCOUP  , ONLY : LWCOU    ,KCOUSTEP  ,LWFLUX ,LWVFLX_SNL,    &
     &            LWCOUAST,                                             &
     &            LWCOUNORMS, LLNORMIFS2WAM,LLNORMWAM2IFS,LLNORMWAMOUT, &
     &            LLNORMWAMOUT_GLOBAL, CNORMWAMOUT_FILE,                &
     &            LWNEMOCOU, LWNEMOCOUSEND, LWNEMOCOURECV,              &
     &            LWNEMOCOUDEBUG, LWNEMOCOUCIC, LWNEMOCOUCIT,           &
     &            LWNEMOCOUCUR,                                         &
     &            LWNEMOCOUSTK,  LWNEMOCOUSTRN, LWNEMOTAUOC, NEMOFRCO,  &
     &            LLCAPCHNK, LLGCBZ0, LLNORMAGAM
      USE YOWCOUT  , ONLY : COUTT    ,COUTS    ,CASS     ,FFLAG    ,    &
     &            FFLAG20  ,GFLAG    ,                                  &
     &            GFLAG20  ,NFLAG    ,                                  &
     &            NFLAGALL ,UFLAG    ,LFDB     ,NOUTT    ,NOUTS    ,    &
     &            NASS     ,JPPFLAG  ,                                  &
     &            IRWDIR   , IRCD    ,IRU10    ,                        &
     &            LRSTPARALW,LRSTPARALR,LRSTINFDAT,                     &
     &            NTRAIN   ,                                            &
     &            IPFGTBL  ,                                            &
     &            LWAMANOUT,                                            &
     &            NWRTOUTWAM,                                           &
     &            LSECONDORDER,                                         &
     &            LWAM_USE_IO_SERV,                                     &
     &            LOUTMDLDCP,                                           &
     &            NGOUT    ,OUTLONG  ,OUTLAT
      USE YOWCPBO  , ONLY : GBOUNC_MAX, IBOUNC ,CBCPREF
      USE YOWCURR  , ONLY : IDELCUR  ,CDATECURA, LLCFLCUROFF
      USE YOWFPBO  , ONLY : IBOUNF
      USE YOWFRED  , ONLY : IFRE1, FR1, XKMSS_CUTOFF 
      USE YOWGRIBHD, ONLY : NGRIB_VERSION, LGRHDIFS,                    &
     &                      LNEWLVTP, LL_GRID_SIMPLE_MATRIX
      USE YOWGRIB_HANDLES , ONLY : NGRIB_HANDLE_IFS, NGRIB_HANDLE_IFS2
      USE YOWGRID  , ONLY : NPROMA_WAM
      USE YOWICE   , ONLY : LICERUN  ,LMASKICE ,LWAMRSETCI, LCIWABR  ,  &
     &            LICETH
      USE YOWMESPAS, ONLY : LFDBIOOUT,LGRIBIN  ,LGRIBOUT ,LNOCDIN
      USE YOWMAP   , ONLY : CLDOMAIN
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,              &
     &            SWAMPWIND,SWAMPWIND2,DTNEWWIND,LTURN90 ,              &
     &            SWAMPCIFR,SWAMPCITH,LWDINTS  ,LL1D     ,LLUNSTR
      USE YOWPHYS  , ONLY : BETAMAX  ,ZALP     ,ALPHA    ,  ALPHAPMAX,  &
     &            TAUWSHELTER, TAILFACTOR, TAILFACTOR_PM

      USE YOWSTAT  , ONLY : CDATEE   ,CDATEF   ,CDATER   ,CDATES   ,    &
     &            IFRELFMAX, DELPRO_LF, IDELPRO  ,IDELT  ,IDELWI   ,    &
     &            IDELWO   ,IDELALT  ,IREST    ,IDELRES  ,IDELINT  ,    &
     &            IDELBC   ,                                            &
     &            ICASE    ,ISHALLO  ,                                  &
     &            IPHYS    ,                                            &
     &            ISNONLIN ,                                            &
     &            IDAMPING ,                                            &
     &            LBIWBK   ,                                            &
     &            IREFRA   ,IPROPAGS ,IASSI    ,                        &
     &            NENSFNB  ,NTOTENS  ,NSYSNB   ,NMETNB   ,CDATEA   ,    &
     &            YCLASS   ,YEXPVER  ,L4VTYPE  ,LFRSTFLD ,LALTAS   ,    &
     &            LSARAS   ,LSARINV  ,ISTREAM  ,NLOCGRB  ,NCONSENSUS,   &
     &            NDWD     ,NMFR     ,NNCEP    ,NUKM     ,IREFDATE ,    &
     &            LGUST    ,LADEN    ,LSUBGRID ,LLSOURCE ,LNSESTART,    &
     &            LSMSSIG_WAM,CMETER ,CEVENT   ,                        &
     &            LRELWIND ,                                            &
     &            IDELWI_LST, IDELWO_LST, CDTW_LST, NDELW_LST
      USE YOWSHAL  , ONLY : NDEPTH   ,DEPTHA   ,DEPTHD    ,TOOSHALLOW
      USE YOWTEST  , ONLY : IU06     ,ITEST    ,ITESTB
      USE YOWTEXT  , ONLY : LRESTARTED,ICPLEN   ,USERID   ,RUNID    ,   &
     &            PATH     ,CPATH    ,CWI
#ifdef WAM_HAVE_UNWAM
      USE YOWUNPOOL, ONLY : LPREPROC, LVECTOR, IVECTOR
      USE UNWAM    , ONLY : LIMPLICIT, JGS_DIFF_SOLVERTHR,              &
     &            SOURCE_IMPL, WAE_SOLVERTHR,                           &
     &            LNONL, BLOCK_GAUSS_SEIDEL,                            &
     &            LLIMT, L_SOLVER_NORM, LCHKCONV
     USE UNSTRUCT_BOUND , ONLY : LBCWA
#endif
      USE YOWWAMI  , ONLY : CBEGDT   ,CENDDT   ,CBPLTDT  ,CEPLTDT  ,    &
     &            CLSPDT   ,CRSTDT   ,IANALPD  ,IFOREPD  ,IDELWIN  ,    &
     &            IASSIM   ,NFCST    ,ISTAT
      USE YOWWIND  , ONLY : CWDFILE  ,LLWSWAVE ,LLWDWAVE ,RWFAC

      USE YOMHOOK  , ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE YOWABORT, ONLY : WAM_ABORT

! ----------------------------------------------------------------------

      IMPLICIT NONE
#include "iwam_get_unit.intfb.h"
#include "wam_u2l1cr.intfb.h"
#include "mpcrtbl.intfb.h"
#include "wposnam.intfb.h"

      INTEGER(KIND=JWIM) :: IU05, LFILE
      INTEGER(KIND=JWIM) :: ISAT, IC, II, ILEN

      REAL(KIND=JWRU) :: R8_DEPTHA, R8_DEPTHD 
      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

      CHARACTER(LEN=14), PARAMETER :: ZERO = ' '
      CHARACTER(LEN=1) :: CLMTSU(4), CLOTSU(8)
      CHARACTER(LEN=2) :: C2
      CHARACTER(LEN=70) :: CLHEADER
      CHARACTER(LEN=:), ALLOCATABLE :: NAMELIST_FILENAME

      LOGICAL :: LLEOF 
      LOGICAL :: LLEXIST

#ifndef WAM_HAVE_UNWAM
      ! dummy read namelist read variables otherwise present in UNWAM module
      REAL(KIND=JWRU) :: WAE_SOLVERTHR
      REAL(KIND=JWRU) :: JGS_DIFF_SOLVERTHR
      LOGICAL         :: LIMPLICIT
      LOGICAL         :: SOURCE_IMPL
      LOGICAL         :: LNONL
      LOGICAL         :: BLOCK_GAUSS_SEIDEL
      LOGICAL         :: LLIMT
      LOGICAL         :: L_SOLVER_NORM
      LOGICAL         :: LCHKCONV
      ! dummy read namelist read variables otherwise present in UNSTRUCT_BOUND module
      LOGICAL         :: LBCWA
      ! dummy read namelist read variables otherwise present in YOWUNPOOL module
      LOGICAL            :: LPREPROC
      LOGICAL            :: LVECTOR
      INTEGER(KIND=JWIM) :: IVECTOR
#endif

! ----------------------------------------------------------------------

      NAMELIST /NALINE/ CLHEADER,                                       &
     &   CLDOMAIN,                                                      &
     &   NANG, IFRE1, FR1, NFRE, NFRE_RED,                              &
     &   CBPLTDT, CEPLTDT, CDATEF,                                      &
     &   IFRELFMAX, DELPRO_LF, IDELPRO, IDELT, IDELWO, IDELWI, CLMTSU,  &
     &   IDELALT, IDELINT, IDELRES,                                     &
     &   IDELCUR, CDATECURA,                                            &
     &   LLCFLCUROFF,                                                   &
     &   CLOTSU, CDATER, CDATES,                                        &
     &   FFLAG,  GFLAG, NFLAG,                                          &
     &   XKMSS_CUTOFF,                                                  &
     &   LFDB, LGRIBIN, LGRIBOUT, LFDBIOOUT,                            &
     &   LRSTPARALW, LRSTPARALR, LRSTINFDAT,                            &
     &   LWAMANOUT,                                                     &
     &   NWRTOUTWAM,                                                    &
     &   LSECONDORDER,                                                  &
     &   ICASE, ISHALLO, ITEST, ITESTB, IREST, IASSI,                   &
     &   IPROPAGS,                                                      &
     &   IREFRA,                                                        &
     &   IPHYS,                                                         &
     &   ISNONLIN,                                                      &
     &   IDAMPING,                                                      &
     &   LBIWBK  ,                                                      &
     &   LMASKICE,                                                      &
     &   LWAMRSETCI,                                                    &
     &   NDEPTH   ,R8_DEPTHA   ,R8_DEPTHD,                              &
     &   IBOUNC, IBOUNF,                                                &
     &   IDELBC, CBCPREF,                                               &
     &   USERID, RUNID,  PATH, YCLASS, YEXPVER, CPATH,                  &
     &   NGRIB_VERSION,                                                 &
     &   NENSFNB, NTOTENS, NSYSNB, NMETNB,                              &
     &   LWCOU, LWCOUAST, LNOCDIN, LODBRALT,                            &
     &   LALTCOR, L4VTYPE, LFRSTFLD, LALTAS, LSARAS, LSARINV, XKAPPA2,  &
     &   IBUFRSAT, CSATNAME,                                            &
     &   SWAMPWIND, SWAMPWIND2, SWAMPCIFR, SWAMPCITH,                   &
     &   DTNEWWIND, LTURN90,                                            &
     &   LALTLRGR, HSCOEFCOR, HSCONSCOR,ALTSDTHRSH,ALTBGTHRSH,HSALTCUT, &
     &   ISTREAM, NLOCGRB, IREFDATE,                                    &
     &   NCONSENSUS, NDWD, NMFR, NNCEP, NUKM,                           &
     &   LGUST, LADEN, LRELWIND, LALTGRDOUT, LSUBGRID, LALTPAS,         &
     &   LLSOURCE,                                                      &
     &   LNSESTART,                                                     &
     &   LLUNSTR, LPREPROC, LVECTOR, IVECTOR,                           &
     &   WAE_SOLVERTHR, JGS_DIFF_SOLVERTHR,                             &
     &   LIMPLICIT,                                                     &
     &   SOURCE_IMPL,                                                   &
     &   LNONL, BLOCK_GAUSS_SEIDEL,                                     &
     &   LLIMT, L_SOLVER_NORM, LCHKCONV,                                &
     &   LBCWA,                                                         &
     &   LSMSSIG_WAM,CMETER,CEVENT,                                     &
     &   LLWSWAVE, LLWDWAVE,                                            &
     &   NPROMA_WAM, LL1D, LGRHDIFS , LNEWLVTP, LL_GRID_SIMPLE_MATRIX,  &
     &   LWCOUNORMS, LLNORMIFS2WAM, LLNORMWAM2IFS, LLNORMWAMOUT,        &
     &   LLNORMWAMOUT_GLOBAL, CNORMWAMOUT_FILE,                         &
     &   LICERUN, LCIWABR, LICETH,                                      &
     &   LWVFLX_SNL,                                                    &
     &   LWNEMOCOU, NEMOFRCO,                                           &
     &   LWNEMOCOUSEND, LWNEMOCOUSTK, LWNEMOCOUSTRN, LWNEMOTAUOC,       &
     &   LWNEMOCOURECV,                                                 &
     &   LWNEMOCOUCIC, LWNEMOCOUCIT, LWNEMOCOUCUR,                      &
     &   LWNEMOCOUDEBUG,                                                &
     &   LLCAPCHNK, LLGCBZ0, LLNORMAGAM,                                &
     &   LWAM_USE_IO_SERV,                                              &
     &   LOUTMDLDCP


      CHARACTER(LEN=14) :: CLOUT
      NAMELIST /NAOT/ CLOUT

      CHARACTER(LEN=14) :: CLSOUT
      NAMELIST /NAOS/ CLSOUT

      CHARACTER(LEN=14) :: CLAOUT
      NAMELIST /NAAT/ CLAOUT

      INTEGER :: IDWI, IDWO
      CHARACTER(LEN=14) :: CLWOUT
      NAMELIST /NAWI/ IDWI, IDWO, CLWOUT

      REAL(KIND=JWRB) :: ZOUTLAT, ZOUTLONG
      NAMELIST /NAOUTP/ ZOUTLAT, ZOUTLONG

!     NAMELIST NALINE : 
!     ===============
!     CLDOMAIN:  CHARACTER DEFINES THE DOMAIN OF THE MODEL 
!     CBPLTDT: USER INPUT START DATE OF RUN.
!     CEPLTDT: USER INPUT END DATE OF RUN.
!     CDATEF: BEGIN DATE OF FORECAST.
!     IFRELFMAX: INDEX FOR THE LAST FREQUENCY THAT WILL BR PROPAGATED WITH TIME STEP DELPRO_LF (see below).
!     DELPRO_LF: PROPAGATION TIME STEP FOR THE ALL WAVES WITH FREQUENCIES <= FR(IFRELFMAX).
!                if 1 <= IFRELFMAX <= NFRE
!                IF COUPLED TO IFS, THE INPUT VALUE OF DELPRO_LF WILL BE
!                THE UPPER BOUND ON THE TIME STEP, IT WILL BE SET SO THAT IT IS A FRACTION OF IDELPRO,
!                WHILE STILL SATSISFYING THE LIMIT IMPOSED BY THE UPPER BOUND.
!     IDELPRO: PROPAGATION TIME STEP FOR THE ALL WAVES WITH FREQUENCIES > FR(IFRELFMAX).
!              IF COUPLED TO IFS, THE INPUT VALUE OF IDELPRO WILL BE
!              THE UPPER BOUND ON THE TIME STEP, WHICH WOULD OTHERWISE
!              BE SET TO THE COUPLING TIME STEP.
!     IDELT: SOURCE TERM INTEGRATION TIME STEP.
!              IF COUPLED TO IFS, THE INPUT VALUE OF IDELT WILL BE
!              THE UPPER BOUND ON THE TIME STEP, WHICH WOULD OTHERWISE
!              BE SET TO THE COUPLING TIME STEP.
!     IDELWO: WIND OUTPUT TIME STEP.
!             (if specific values are not specified by namelist NAWI).
!     IDELWI: WIND INPUT TIME STEP.
!             (if specific values are not specified by namelist NAWI).
!     CLMTSU: STEP UNIT (S : seconds or H : hours).
!     IDELALT: ALTIMETER DATA TIME WINDOW (in seconds).
!     IDELINT: OUTPUT TIME STEP FOR INTEGRATED PARAMETER OF TOTAL SEA
!              (if output times are not specified by namelist NAOT).
!     IDELRES: OUTPUT TIME STEP for RESTART SPECTRA (if output times
!              are not specified by namelist NAOS).
!     IDELCUR: CURRENTS INPUT TIME STEP IN SECONDS. 
!     CDATECURA: FIRST DATE FOR THE CURRENT INPUT.
!     LLCFLCUROFF: IF TRUE PREVENT CFL CRITERIA FAILURE DUE TO SURFACE CURRENTS
!                  IT RESETS THE CURRENT REFRACTION TERMS TO 0 FOR THOSE POINTS
!                  WHERE CFL WAS NOT SATISFIED (it should only be used in operational applications)
!     CLOTSU: STEP UNIT (S : seconds or H : hours) FOR IDELINT,
!             IDELRES, IDELBC.
!     CDATER: ONE SPECIFIC OUTPUT TIME FOR RESTART FILES if output times
!             are not specified by namelist NAOS).
!     CDATES: LAST OUTPUT TIME ALLOWED BESIDE CDATER FOR RESTART SPECTRA
!             (if output times are not specified by namelist NAOS).
!     FFLAG: OUTPUT FLAG FOR OUTPUT TO USER OUTPUT OF EACH OUTPUT TYPE.
!     GFLAG: OUTPUT FLAG FOR OUTPUT TO GRIB OF EACH OUTPUT TYPE.
!     NFLAG: OUTPUT FLAG FOR USER OUTPUT DISPLAY OF FIELD NORM FOR
!            EACH OUTPUT TYPE.
!     XKMSS_CUTOFF: IF DIFFERENT FROM 0., SETS THE MAXIMUM WAVE NUMBER TO BE USED IN
!                   THE CALCULATION OF THE MEAN SQUARE SLOPE.
!                   OTHERWISE, USE XK_GC(NWAV_GC)
!     TYPE OF INTEGRATED PARAMETERS IN FFLAG GFLAG (see OUTINT) :
!     GFLAG and FFLAG SEE *MPCRTBL* for a definition of which parameters they control
!     LFDB: TRUE IF INTEGRATED PARAMETERS ARE WRITTEN OUT TO FDB (see WVWAMINIT1 when coupled to IFS).
!     LGRIBIN : IF TRUE THE WAVE SPECTRA IS INPUT IN GRIB FORMAT, ELSE
!               THE BINARY RESTART FILES ARE USED.
!     LGRIBOUT : IF TRUE THE WAVE SPECTRA IS OUTPUT IN GRIB FORMAT, ELSE
!                THE BINARY RESTART FILES ARE PRODUCED.
!     LFDBIOOUT : IF TRUE THE GRIB SPECTRA OUTPUT IS SENT TO THE FDB (See WVWAMINIT1 when coupled to IFS).
!     NWRTOUTWAM: STRIDE WITH WHICH FDB OUTPUT PE's ARE SELECTED.
!     LRSTPARALW: IF TRUE BINARY RESTART FILES WILL BE WRITTEN in PARALLEL
!                 IF THAT OPTION IS USED.
!     LRSTPARALR: IF TRUE BINARY RESTART FILES WILL BE READ in PARALLEL
!                 IF THAT OPTION IS USED.
!     LRSTINFDAT:  F TRUE WILL WRITE AN ADDITIONAL wamfile WITH DATE/TIME INFO
!     LSECONDORDER : IF TRUE THEN SECOND ORDER CORRECTIOn IS APPLIED TO
!                    OUTPUT INTEGRATED PARAMETERS.
!     ICASE: 1 FOR SPHERICAL COORDINATES ELSE CARTESIAN COORDINATES.
!     ISHALLO: THIS OPTION IS DEPRICATED. IT IS ALWAYS SHALLOW WATER FORMULATION
!     IREFRA: 0 MODEL RUNS WITHOUT REFRACTION.
!             1 MODEL RUNS WITH DEPTH REFRACTION ONLY.
!             2 MODEL RUNS WITH CURRENT REFRACTION ONLY.
!             3 MODEL RUNS WITH DEPTH AND CURRENT REFRACTION.
!     ITEST: TEST OUTPUT LEVEL.
!     ITESTB: MAX BLOCK NUMBER FOR OUTPUT IN BLOCK LOOPS.
!     IREST: 1 FOR THE PRODUCTION OF RESTART FILE(S).
!     IASSI: 1 ASSIMILATION IS DONE IF ANALYSIS RUN.
!     IPHYS:  WAVE PHYSICS PACKAGE (0 or 1)
!     ISNONLIN : 0 FOR OLD SNONLIN, 1 FOR NEW SNONLIN, 2 FOR LATEST BASED ON JANSSEN 2018 (ECMWF TM 813).
!     IDAMPING : 0 NO WAVE DAMPING, 1 WAVE DAMPING ON.
!                ONLY MEANINGFUl FOR IPHYS=0
!     NDEPTH   NUMBER OF ENTRIED IN DEPTH TABLE (if used) 
!     DEPTHA   MINIMUM DEPTH. IT IS ALSO USED IN DEPTH TABLE (if used)
!     DEPTHD   MAXIMUM DEPTH IN TABLES = DEPTHA*DEPTHD**(NDEPTH-1)
!     IBOUNC: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (COARSE GRID).
!     IBOUNF: 1 FOR RUN WHICH INCLUDES BOUNDARY POINTS (FINE GRID).
!     IDELBC: TIMESTEP FOR THE DISPOSAL OF OUTPUT BC FILE(S) INTO PERMANENT
!             FILE(S). IF SET TO ZERO THEN IT WILL BE RESET TO THE LENGTH
!             OF THE RUN (I.E. THE FILE(S) WILL DISPOSED AT THE END OF THE
!             RUN.
!     CBCPREF: PREFIXES OF BC FILES
!     USERID, RUNID, PATH : OUT OF DATE, SHOULD ONLY BE USED WHENEVER
!                           USER WANTS TO WRITE OUTPUT DIRECTLY TO ECFS
!                           WHICH IS HIGHLY NOT RECOMMENDATED ON THE VPP
!     YCLASS: DATA CLASS USED TO CODE DATA IN GRIB.
!     YEXPVER: EXPERIMENT VERSION USED TO CODE DATA IN GRIB.
!     CPATH: PATH FOR OUTPUT TO DISK (USED WHENEVER OUTPUT IS NOT
!            TO THE FIELD DATA BASE (SEE LFDBIOOUT)).
!     NGRIB_VERSION: GRIB VERSION FOR OUTPUT.
!     NTOTENS: TOTAL ENSEMBLE FORECAST MEMBERS (DEFAULT=0).
!     NENSFNB: ENSEMBLE FORECAST NUMBER (DEFAULT=0).
!     NSYSNB : SYSTEM NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     NMETNB : METHOD NUMBER TO BE USED FOR GRIBBING OF SEASONAL DATA.
!              or MONTHLY FORECAST RUNS
!     LWCOU: FALSE FOR UNCOUPLED RUN (see WVWAMINIT1 if coupled to IFS).
!     LWCOUAST: IF TRUE AND LWCOU THEN THE ATMOSPHERIC SURFACE STRESS OVER THE OCEAN WILL BE USED
!               TO COMPUTE THE OCEAN SIDE STRESS
!     NEMO COUPLING FLAGS:
!     LWNEMOCOU: FALSE FOR NO COUPLING TO NEMO RUN.
!     NEMOFRCO: NEMO COUPLING FREQ IN WAM TIME STEPS
!     LWNEMOCOUSEND: IF FALSE THEN NO DATA WILL BE SENT TO NEMO
!     LWNEMOCOUSTK:  IF TRUE SEND SURFACE STOKES DRIFT TO NEMO
!     LWNEMOCOUSTRN: IF TRUE SEND SEA ICE WAVE STRAIN TO NEMO
!     LWNEMOTAUOC IF TRUE THEN THE ATMOSPHERIC STRESS THAT IS PASSED TO
!                 NEMO VIA THE WAVE MODEL WILL ALSO BE MODULATED BY IT.
!
!     LWNEMOCOURECV: IF FALSE THEN NO DATA WILL BE RECEIVED FROM NEMO
!     LWNEMOCOUCIC: IF TRUE THEN SEA ICE CONCENTRATION WILL BE OBTAINED FROM FROM NEMO 
!     LWNEMOCOUCIT: IF TRUE THEN SEA ICE THICKNESS WILL BE OBTAINED FROM FROM NEMO 
!     LWNEMOCOUCUR: IF TRUE THEN SURFACE CURRENTS WILL BE OBTAINED FROM FROM NEMO 
!
!     LWNEMOCOUDEBUG: FALSE IF NO DEBUGGING OUTPUT IN WAM<->NEMO COUPLING
!
!     LLCAPCHNK : CAP CHARNOCK FOR HIGH WINDS.
!     LLGCBZ0 : USE MODEL FOR BACKGROUND ROUGHNESS.
!     LLNORMAGAM : USE THE RENORMALISTION OF THE GROWTH RATE.
!     LWAM_USE_IO_SERV: TRUE IF SPECTRAL AND INTEGRATED PARAMETER OUTPUT SHOULD BE
!              DONE USING IFS IO SERVER
!     LOUTMDLDCP : OUPUT MODEL DECOMPOSITION TO FILE
!
!     LNOCDIN: IF TRUE THEN GRIB INPUT OF A DRAG COEFFICIENT FIELD IS
!              NOT REQUIRED.  
!     LODBRALT: IF TRUE THEN THE ALTIMETER DATA WILL BE READ AND PASSED 
!               THROUGH OBSERVATION DATABASE (ODB)
!     LALTCOR: IF TRUE THEN THE ALTIMETER DATA WILL BE CORRECTED
!              SEE GRFIELD(but this is different than the bias correction scheme)
!              It was implemented when the ERS altimeters were used in operation.
!              We has since then gone for an overall bias correction scheme and so by default this will not be used.. 
!     L4VTYPE: IF TRUE THEN MARS TYPE 4V IS USED INSTEAD OF FG.
!     LFRSTFLD: IF TRUE THEN INITIAL INTEGRATED PARAMETER FIELDS WILL
!               OUTPUT PROVIDED THEY WERE NOT INPUT.
!     XKAPPA2: KAPPA2 PARAMETER USED IN THE FORMULA FOR THE
!              DETERNMINATION OF ALTIMETER WAVE HEIGHT.
!     LWAMANOUT: CONTROLS WHETHER FIELDS WILL BE WRITTEN OUT AT ANALYSIS TIME OR NOT.
!                (USEFUL WHEN WAVE DATA ASSIMILATION IS DONE IN MORE THAN ONE TRAJECTORY.)
!     LALTAS: CONTROLS WHETHER ALTIMETER DATA ARE ASSIMILATED
!     LSARAS: CONTROLS WHETHER SAR DATA ASSIMILATED
!     LSARINV: CONTROLS WHETHER SAR INVERSION IS DONE
!     SWAMPWIND : CONSTANT WIND SPEED USED BY SWAMP CASE.
!     SWAMPWIND2 : SECOND WIND SPEED VALUE TO BE USED (if different than 0)
!                  IN SWAMP CASE. IT WILL BE APPLIED DTNEWWIND HOURS
!                  AFTER STATING TIME.
!     SWAMPCIFR : CONSTANT SEA ICE FRACTION FOR THE 1/2 NORTHERN PART OF SWAMP DOMAIN
!     SWAMPCITH : CONSTANT SEA ICE THICKNESS FOR THE 1/2 NORTHERN PART OF SWAMP DOMAIN
!     DTNEWWIND : TIME IN HOURS AFTER WHICH SWAMPWIND2 WILL BE APPLIED
!     LTURN90 : IF TRUE THE NEW WIND IN SWAMP CASE WILL TURN BY 90 DEGREES.
!     LALTLRGR: IF TRUE THEN THE ALTIMETER DATA WILL BE CORRECTED
!               PRIOR TO THEIR ASSIMILATION USING A LINEAR REGRESSION AS
!               PROVIDED BY HSCOEFCOR AND HSCONSCOR (see below)
!               This has been replaced by the bias correction scheme
!     HSCOEFCOR: COEFFICIENT OF THE CORRECTIVE LINEAR REGRESSION FOR
!                ALTIMETER WAVE HEIGHTS.
!     HSCONSCOR: CONSTANT OF THE CORRECTIVE LINEAR REGRESSION FOR
!                ALTIMETER WAVE HEIGHTS.
!     ALTSDTHRSH:THRESHOLD FOR SUSPICIOUS DATA (SEE GRFIELD).
!     ALTBGTHRSH:THRESHOLD FOR BACKGROUND CHECK (SEE GRFIELD).
!     HSALTCUT: USER INPUT OF THE MINIMUM WAVE HEIGHT ALLOWED IN ALTAS
!              (SEE GRFIELD). 
!     ISTREAM: STREAM NUMBER USED WHEN GRIBBING THE DATA 
!                     IT MUST BE SPECIFIED EXCEPT IF LGRHDIFS IS TRUE ! 
!     LALTGRDOUT: FLAG USED TO CONTROL THE OUTPUT OF GRIDDED ALTIMETER
!                 PRODUCTS OF SPECIFIC INSTRUMENT (TRUE TO OUTPUT) 
!                (SEE GRFIELD). 
!     LALTPAS: FLAG USED TO CONTROL FEED THROUGH ALT DATA PASSIVELY. 
!                (SEE GRFIELD). 
!     LGUST:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO GUSTINESS 
!     LADEN:   FLAG USED TO ACTIVATE COMPUTATIONS RELATED TO AIR DENSITY
!     LRELWIND: IF TRUE THEN RELATIVE WINDS ARE USED WITH RESPECT TO
!               SURFACE CURRENTS.
!     LLWSWAVE: FLAG USE TO ACTIVATE USE OF WAVE PARAMETER WIND SPEED AS
!               INPUT TO CONSTRUCT THE WIND FORCING - UNCOUPLED RUNS ONLY
!     LLWDWAVE: FLAG USE TO ACTIVATE USE OF WAVE PARAMETER WIND DIRECTION
!               AS INPUT TO CONSTRUCT THE WIND FORCING - UNCOUPLED RUNS ONLY.
!     NLOCGRB: LOCAL GRIB TABLE NUMBER
!     IREFDATE: REFERENCE DATE FOR MONTHLY FORECAST HINDCAST.
!     NCONSENSUS: ORIGINE OF THE INITIAL CONDITONS FOR
!                 MULTI-ANALYSIS ENSEMBLE RUNS
!     NDWD: IF 1 DWD ANALYSES ARE USED IN THE MULTI-ANALYSIS ENSEMBLE
!           RUNS
!     NMFR: IF 1 METEO FRANCE ANALYSES ARE USED IN THE MULTI-ANALYSIS
!           ENSEMBLE RUNS
!     NNCEP: IF 1 NCEP ANALYSES ARE USED IN THE MULTI-ANALYSIS ENSEMBLE
!           RUNS
!     NUKM: IF 1 MET OFFICE ANALYSES ARE USED IN THE MULTI-ANALYSIS
!           ENSEMBLE RUNS
!     NPROMA_WAM: NUMBER OF THE GRID POINTS PER LOOP WHEN CUT INTO CHUNKS
!     LL1D: IF TRUE THEN THE DOMAIN IS ONLY DIVIDED INTO LATITUDINAL BANDS
!     LWCOUNORMS : FLAG THAT CONTROLS WHETHER OR NOT GLOBAL NORMS
!                  ARE PRODUCED FOR THE DIFFERENT FIELDS EXCHANGED 
!                  BETWEEN WAM AND THE COUPLED ATMOSPHERIC MODEL.
!                  IT IS FALSE BY DEFAULT.
!     LLNORMIFS2WAM : IF TRUE NORMS FOR FIELDS PASSED FROM IFS TO WAM WILL BE PRODUCED
!     LLNORMWAM2IFS : IF TRUE NORMS FOR FIELDS PASSED FROM WAM TO IFS WILL BE PRODUCED
!     LLNORMWAMOUT : IF TRUE NORMS OF SELECTED OUPTUT FIELDS WILL BE PRODUCED
!     LLNORMWAMOUT_GLOBAL : IF TRUE NORMS FOR SELECTED OUPTUT FIELDS WILL BE GLOBAL (see above)
!     LGRHDIFS : FLAGS CONTROLLING WHETHER OR NOT GRIB HEADER INFORMATION
!                IS COPIED FROM THE ATMOSPHERIC MODEL (ONLY USEFULL IF
!                COUPLED TO THE IFS).
!     LNEWLVTP : FLAG CONTROLLING WHETHER OR NOT THE NEW GRIB LEVTYPE
!                DEFINITIONS ARE USED.
!     LL_GRID_SIMPLE_MATRIX IF TRUE THEN THE 2D SPECTRA WILL USE THE LEGACY grid_simple_matrix
!                           TO ENCODE THE 2D SPECTRA in GRIB1. THIS SHOULD BE PHASED OUT as soon as feasible!
!     LICERUN : FLAG CONTROLLING WHETHER OR NOT SEA ICE FRACTION (OR SST)
!               FIEDS ARE PROVIDED WITH THE WIND FIELDS TO GENERATE THE
!               SEA ICE MASK (TRUE BY DEFAULT). 
!     LCIWABR : FLAG CONTROLLING THE USE OF SEA ICE BOTTOM FRICTION ATTENUATION  
!     LMASKICE  SET TO TRUE IF ICE MASK IS APPLIED
!     LWAMRSETCI SET TO TRUE IF FIELDS THAT ARE EXCHANGED WITH THE ATMOSPHERE AND THE OCEAN
!                ARE RESET TO WHAT WOULD BE USED IF THERE WERE NO WAVE MODELS.

!     LICETH  : FLAG CONTROLLING WHETHER OR NOT SEA ICE THICKNESS FILEDS
!               ARE PROVIDED WITH THE WIND FIELDS (FALSE BY DEFAULT). 
!     LLSOURCE : FLAG CONTROLLING WHETHER OR NOT THE SOURCE TERM CONTRIBUTION 
!                IS COMPUTED.
!     LNSESTART : FLAG CONTROLLING WHETHER OR NOT THE INITIAL SPECTRA ARE
!                 RESET TO NOISE LEVEL.
!     LSMSSIG_WAM : .T. = send signals to ECFLOW (ECMWF supervisor)
!     CMETER :  SMS or ECFLOW meter command (ECMWF supervisor)
!     CEVENT :  SMS or ECFLOW event command (ECMWF supervisor)


!     NAMELIST NAOT : 
!     ===============
!     CLOUT : LIST OF OUTPUT TIME FOR INTEGRATED PARAMETERS.
      
!     NAMELIST NAOS : 
!     ===============
!     CLSOUT : LIST OF OUTPUT TIME FOR SPECTRA. 

!     NAMELIST NAAT : 
!     ===============
!     CLAOUT : LIST OF TIME FOR ALTIMETER DATA ASSIMILATION. 

!     NAMELIST NAWI (only meaningful for uncoupled runs): 
!     ===============
!     IDWI : LIST OF IDELWI
!     IDWO : LIST OF IDELWO
!     CLWOUT : LIST OF TIME UNTIL WHICH IDWI AND IDWO ARE VALID. 

! ----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPUSERIN',0,ZHOOK_HANDLE)

!*    0. SET DEFAULT VALUES FOR THE NAMELIST ELEMENTS.
!        ---------------------------------------------

      CLDOMAIN  = 'g'

      NANG      = 0
      IFRE1     = 3
      FR1       = 4.177248E-02_JWRB
      NFRE      = 0
      NFRE_RED  = 0

      CLMTSU    = 'S'
      CLOTSU    = 'H'
      CLHEADER  =  ZERO
      CBPLTDT   =  ZERO
      CEPLTDT   =  ZERO
      CDATEF    =  ZERO
      IFRELFMAX = 0
      DELPRO_LF =  0.0_JWRB
      IDELPRO   =  0
      IDELT     =  0
      IDELWO    =  0
      IDELWI    =  0
      IDELALT   =  21600 
      IDELINT   =  0

      IDELCUR   =  0
      CDATECURA = ZERO
      LLCFLCUROFF = .TRUE.
      CDATER    = ZERO
      CDATES    = ZERO
      FFLAG(:)  = .FALSE. 
!     note: if IREST>0 then GFLAG(IR) for IR=IRWDIR, IRCD, IRU10 will always be reset
!     true if grib restart files are requested (LGRIBOUT=.true.).
!     see below.
      GFLAG(:)  = .FALSE.
      NFLAG(:)  = .FALSE. 
!     patial initialisation of the output table
      CALL MPCRTBL
      NFLAG(1)  = .TRUE. 
      NFLAG(6)  = .TRUE. 
      NFLAG(IRWDIR)= .TRUE. 
      NFLAG(IRCD)  = .TRUE. 
      NFLAG(IRU10) = .TRUE. 

      XKMSS_CUTOFF = 0.0_JWRB

      LSECONDORDER = .TRUE.

      LFDB      = .TRUE. 
      LGRIBIN   = .TRUE. 
      LGRIBOUT  = .TRUE. 
      LFDBIOOUT = .TRUE. 
      NWRTOUTWAM= 1
      LRSTPARALW= .TRUE.
      LRSTPARALR= .TRUE.
      LRSTINFDAT= .FALSE.
      LODBRALT  = .FALSE.
      ICASE     = 1 
      ISHALLO   = 0   !! depricated 
      IPHYS     = 1
      ISNONLIN  = 1 
      IDAMPING  = 1 
      IPROPAGS  = 0 
      IREFRA    = 0 
      ITEST     = 0 
      ITESTB    = 0 
      IREST     = 0 
      IASSI     = 0 
      NDEPTH    = 74
      R8_DEPTHA = 1.0_JWRU
      R8_DEPTHD = 1.1_JWRU
      IBOUNC    = 0 
      IBOUNF    = 0 
      IDELBC    = 0
      DO II=1,GBOUNC_MAX
        WRITE(C2,'(I2.2)') II
        CBCPREF(II)='B'//C2
      ENDDO
      USERID    = ZERO
      RUNID     = ZERO
      PATH      = ZERO
      CPATH     = ZERO
      YCLASS    = ZERO
      YEXPVER   = ZERO
      NENSFNB   = 0
      NTOTENS   = 0
      NSYSNB    = -1
      NMETNB    = -1

      NOUTT     = 0

      LNOCDIN   = .FALSE. 

      DO ISAT=1,NUMALT
        IBUFRSAT(ISAT) = 0
        CSATNAME(ISAT) = 'SATELLITE NAME MISSING'
        LALTPAS(ISAT)   = .FALSE. 
        LALTCOR(ISAT)   = .FALSE. 
        XKAPPA2(ISAT)   = 0.0700_JWRB
        LALTLRGR(ISAT)  = .FALSE.
        HSCOEFCOR(ISAT) = 1.0_JWRB
        HSCONSCOR(ISAT) = 0.0_JWRB

        ALTBGTHRSH(ISAT) = 1.5_JWRB
!       if no value is provided in the namelist ALTSDTHRSH will
!       be set in grfield.
        ALTSDTHRSH(ISAT) = -1.0_JWRB

!       HSALTCUT is used in combination with the error estimate of
!       the altimeter data to determine the minimum Hs allowed for
!       altimeter data.
        HSALTCUT(ISAT) = 999999._JWRB
        LALTGRDOUT(ISAT)  = .FALSE.
      ENDDO

      L4VTYPE   = .FALSE. 

      LFRSTFLD  = .FALSE.

      LALTAS    = .FALSE. 
      LSARAS    = .FALSE. 
      LSARINV   = .FALSE. 

      LWAMANOUT = .TRUE.

      SWAMPWIND = 18.45_JWRB
      SWAMPWIND2 = 0.0_JWRB
      DTNEWWIND = 0.0_JWRB
      SWAMPCIFR = 0.0_JWRB
      SWAMPCITH = 0.0_JWRB
      LTURN90 = .FALSE.

      ISTREAM   = 0

      NLOCGRB   = 1

      IREFDATE  = 0

      NCONSENSUS= 0
      NDWD      = 0
      NMFR      = 0
      NNCEP     = 0
      NUKM      = 0

      NOUTS     = 0

      NASS      = 0

      LGUST = .FALSE.
      LADEN = .FALSE.

      LLWSWAVE = .FALSE.
      LLWDWAVE = .FALSE.

      LSUBGRID = .TRUE.

      NPROMA_WAM = 1
      NPROMA_WAM = HUGE(NPROMA_WAM)/2

      LL1D = .TRUE.

      LWCOU = .FALSE.

      LWCOUAST = .TRUE.

      LWVFLX_SNL = .TRUE.

      LWNEMOCOU = .FALSE.

      NEMOFRCO = 0

      LWNEMOCOUSEND = .TRUE.

      LWNEMOCOUSTK=.TRUE.

      LWNEMOCOUSTRN=.FALSE.

      LWNEMOTAUOC = .TRUE.

      LWNEMOCOURECV = .FALSE.

      LWNEMOCOUCIC=.FALSE.

      LWNEMOCOUCIT=.FALSE.

      LWNEMOCOUCUR=.FALSE.

      LWNEMOCOUDEBUG = .FALSE.

      LLCAPCHNK = .FALSE.
      LLGCBZ0 = .TRUE.
      LLNORMAGAM = .TRUE.

      LWAM_USE_IO_SERV = .FALSE.

      LOUTMDLDCP = .FALSE.

      LWCOUNORMS = .FALSE.
      LLNORMIFS2WAM = .FALSE.
      LLNORMWAM2IFS = .FALSE.
      LLNORMWAMOUT = .FALSE.
      LLNORMWAMOUT_GLOBAL = .FALSE.
      CNORMWAMOUT_FILE = ''

      LGRHDIFS = .FALSE.

      NGRIB_VERSION = 1

      LNEWLVTP = .FALSE.

      LL_GRID_SIMPLE_MATRIX = .TRUE.

      LBIWBK = .TRUE.

      LICERUN = .TRUE.

      LCIWABR = .FALSE.

      LMASKICE = .FALSE.

      LWAMRSETCI = .TRUE.

      LICETH = .FALSE.

      LLSOURCE = .TRUE.

      LNSESTART = .FALSE.

      LRELWIND = .TRUE.

      NDELW_LST = 0

      LSMSSIG_WAM =.FALSE.

      CMETER='smsmeter'

      CEVENT='smsevent'

      LLUNSTR=.FALSE.
      LVECTOR=.FALSE.
      IVECTOR=1
      LPREPROC=.FALSE.
      JGS_DIFF_SOLVERTHR = 1.E-5_JWRU
      WAE_SOLVERTHR = 1.E-10_JWRU
      LIMPLICIT = .FALSE.
      SOURCE_IMPL = .FALSE.
      LNONL = .FALSE.
      BLOCK_GAUSS_SEIDEL = .TRUE.
      LLIMT = .FALSE.
      L_SOLVER_NORM = .FALSE.
      LCHKCONV = .TRUE.
      LBCWA = .FALSE.

      NGRIB_HANDLE_IFS = -1
      NGRIB_HANDLE_IFS2 = -1


! ----------------------------------------------------------------------

!*    1. READ NAMELIST NALINE.
!        ---------------------

      NAMELIST_FILENAME='wam_namelist'

      LFILE=0
      LLEXIST=.FALSE.
      IF (NAMELIST_FILENAME /= ' ') LFILE=LEN_TRIM(NAMELIST_FILENAME)
      INQUIRE(FILE=NAMELIST_FILENAME(1:LFILE),EXIST=LLEXIST)
      IF (.NOT. LLEXIST) THEN
        CALL WAM_ABORT("INPUT FILE '"//NAMELIST_FILENAME(1:LFILE)//"' IS MISSING !!!",__FILENAME__,__LINE__)
      ENDIF

      WRITE(IU06,*) ' '
      WRITE(IU06,'(3A)') '  WAVE MODEL NAMELIST "', NAMELIST_FILENAME,'" IS BEING READ'

      IU05 =  IWAM_GET_UNIT (IU06, NAMELIST_FILENAME, 's', 'f', 0, 'READ')

      CALL WPOSNAM (IU05, 'NALINE', LLEOF)
      IF (.NOT. LLEOF) THEN
        READ (IU05, NALINE)
      ELSE
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'+ SUBROUTINE MPUSERIN :                    +'
        WRITE(IU06,*)'+ READ NAMELIST FAILED                     +' 
        WRITE(IU06,*)'+ NAMELIST FILENAME: ',NAMELIST_FILENAME
        WRITE(IU06,*)'+ PROGRAM WILL ABORT                       +'
        WRITE(IU06,*)'+                                          +'
        WRITE(IU06,*)'++++++++++++++++++++++++++++++++++++++++++++'
        CALL WAM_ABORT("Could not read namelist '"//NAMELIST_FILENAME//"'",__FILENAME__,__LINE__)
      ENDIF

      IF( NANG <= 0 ) THEN
         CALL WAM_ABORT( "Expected positive value for NANG in '"//NAMELIST_FILENAME//"' ", __FILENAME__, __LINE__ )
         CALL WAM_ABORT( "Did you provide a meaningful namelist in '"//NAMELIST_FILENAME//"' ?", __FILENAME__, __LINE__ )
      ENDIF
      IF( NFRE <= 0 ) CALL WAM_ABORT( "Expected positive value for NFRE", __FILENAME__, __LINE__ )
      IF( NFRE_RED <= 0 ) NFRE_RED = NFRE
      IF( IFRE1 <= 0 ) CALL WAM_ABORT( "Expected positive value for IFRE1",  __FILENAME__, __LINE__ )

      IF (NFRE_RED > NFRE ) THEN
        WRITE (IU06,*) '**********************************************'
        WRITE (IU06,*) '*                                            *'
        WRITE (IU06,*) '*       FATAL ERROR IN SUB. UIPREP           *'
        WRITE (IU06,*) '*       ==========================           *'
        WRITE (IU06,*) '* THE REDUCED NUMBER OF FREQUENCIES NFRE_RED *'
        WRITE (IU06,*) '* IS LARGER THAN THE TOTAL NUMNBER NFRE  !!  *'
        WRITE (IU06,*) '* NFRE_RED = ', NFRE_RED
        WRITE (IU06,*) '* NFRE     = ', NFRE
        WRITE (IU06,*) '**********************************************'
        CALL WAM_ABORT(__FILENAME__,__LINE__)
      ENDIF


      IF (IRANK > 1) THEN
         CNORMWAMOUT_FILE = ''
      ENDIF

      ! when coupled to IFS, the control will come from it via calls to wavemdl
      IF (LWCOU) LSMSSIG_WAM=.FALSE.

      DEPTHA = R8_DEPTHA
      DEPTHD = R8_DEPTHD

      TOOSHALLOW=0.1_JWRB*DEPTHA


!           **** OUTPUT TIME AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPT0: DO
        CALL WPOSNAM (IU05, 'NAOT', LLEOF)
        IF (LLEOF) EXIT SPT0
        READ (IU05, NAOT, END=1900)
        NOUTT = NOUTT+1
      ENDDO SPT0
1900  CONTINUE
      REWIND(IU05)
      IF ( NOUTT > 0) THEN
        IF (.NOT.ALLOCATED(COUTT)) ALLOCATE(COUTT(NOUTT))
      ENDIF
      SPT: DO IC=1,NOUTT
        CALL WPOSNAM (IU05, 'NAOT', LLEOF)
        IF (LLEOF) EXIT SPT
        READ (IU05, NAOT, END=1910)
        COUTT(IC) = CLOUT
      ENDDO SPT
1910  CONTINUE

!           **** SPECTRA OUTPUT TIME AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPS0: DO
        CALL WPOSNAM (IU05, 'NAOS', LLEOF)
        IF (LLEOF) EXIT SPS0
        READ (IU05, NAOS, END=1901)
        NOUTS = NOUTS+1
      ENDDO SPS0
1901  CONTINUE
      REWIND(IU05)
      IF ( NOUTS > 0) THEN
        IF (.NOT.ALLOCATED(COUTS)) ALLOCATE(COUTS(NOUTS)) 
      ENDIF
      SPS: DO IC=1,NOUTS
        CALL WPOSNAM (IU05, 'NAOS', LLEOF)
        IF (LLEOF) EXIT SPS
        READ (IU05, NAOS, END=1911)
        COUTS(IC) = CLSOUT
      ENDDO SPS
1911  CONTINUE

!           **** ASSIMILATION AT SPECIFIED TIMES ****
      REWIND(IU05)
      SPA0: DO
        CALL WPOSNAM (IU05, 'NAAT', LLEOF)
        IF (LLEOF) EXIT SPA0
        READ (IU05, NAAT, END=1902)
        NASS = NASS+1
      ENDDO SPA0
1902  CONTINUE
      REWIND(IU05)
      IF ( NASS > 0) THEN
        IF (.NOT.ALLOCATED(CASS)) ALLOCATE(CASS(NASS))
      ENDIF
      SPA: DO IC=1,NASS
        CALL WPOSNAM (IU05, 'NAAT', LLEOF)
        IF (LLEOF) EXIT SPA
        READ (IU05, NAAT, END=1912)
        CASS(IC) = CLAOUT
      ENDDO SPA
1912  CONTINUE

      IF (.NOT.LWCOU) THEN
!           **** FORCING TIME STEPPING AT SPECIFIED TIMES ****
        REWIND(IU05)
        SPWI: DO
          CALL WPOSNAM (IU05, 'NAWI', LLEOF)
          IF (LLEOF) EXIT SPWI
          READ (IU05, NAWI, END=1903)
          NDELW_LST = NDELW_LST + 1 
        ENDDO SPWI
1903    CONTINUE
        REWIND(IU05)
        IF ( NDELW_LST > 0) THEN
          IF (.NOT.ALLOCATED(IDELWI_LST)) ALLOCATE(IDELWI_LST(NDELW_LST))
          IF (.NOT.ALLOCATED(IDELWO_LST)) ALLOCATE(IDELWO_LST(NDELW_LST))
          IF (.NOT.ALLOCATED(CDTW_LST)) ALLOCATE(CDTW_LST(NDELW_LST))
        ENDIF
        SPW: DO IC=1,NDELW_LST
          CALL WPOSNAM (IU05, 'NAWI', LLEOF)
          IF (LLEOF) EXIT SPW
                  READ (IU05, NAWI, END=1913)
                  IDELWI_LST(IC) = IDWI 
                  IDELWO_LST(IC) = IDWO 
                  CDTW_LST(IC) = CLWOUT
                ENDDO SPW
        1913    CONTINUE
      ENDIF

!           **** OUTPUT POINTS READ FROM NAMELIST NAOUTP ****
      NGOUT = 0
      REWIND(IU05)
      SPOU: DO
        CALL WPOSNAM (IU05, 'NAOUTP', LLEOF)
        IF (LLEOF) EXIT SPOU
        READ (IU05, NAOUTP, END=1904)
        NGOUT = NGOUT+1
      ENDDO SPOU
1904  CONTINUE
      REWIND(IU05)
      IF (NGOUT > 0) THEN
        IF (.NOT.ALLOCATED(OUTLAT)) ALLOCATE(OUTLAT(NGOUT))
        IF (.NOT.ALLOCATED(OUTLONG)) ALLOCATE(OUTLONG(NGOUT))
      ENDIF
      SPOUS: DO IC=1,NGOUT
        CALL WPOSNAM (IU05, 'NAOUTP', LLEOF)
        IF (LLEOF) EXIT SPOUS
        READ (IU05, NAOUTP, END=1914)
        OUTLAT(IC) = ZOUTLAT 
        OUTLONG(IC) = ZOUTLONG 
      ENDDO SPOUS
1914  CONTINUE


      CLOSE(IU05)

!           **** MODEL TIME STEPS ****
      IF (CLMTSU(1) == 'H') IDELPRO = IDELPRO*3600
      IF (CLMTSU(2) == 'H') IDELT   = IDELT*3600
      IF (CLMTSU(3) == 'H') IDELWO  = IDELWO*3600
      IF (CLMTSU(4) == 'H') IDELWI  = IDELWI*3600
!           **** OUTPUT TIME IN FIXED INTERVALS ****
      IF (CLOTSU(1) == 'H') IDELINT = IDELINT*3600
      IF (CLOTSU(7) == 'H') IDELRES = IDELRES*3600
      IF (CLOTSU(8) == 'H') IDELBC  = IDELBC*3600


      CALL WAM_U2L1CR( YCLASS )



      ICPLEN=LEN_TRIM(CPATH)
      IF (ICPLEN > 0) THEN
        IF (CPATH(ICPLEN:ICPLEN).EQ.'/') THEN
          CPATH=CPATH(1:ICPLEN-1)
          ICPLEN=ICPLEN-1
        ENDIF
      ENDIF

!           **** CHECK LENGTH OF YEXPVER AND PUT IT RIGHT JUSTIFIED ****
      ILEN=LEN_TRIM(YEXPVER)
      YEXPVER(5-ILEN:4)=YEXPVER(1:ILEN)
      YEXPVER(1:4-ILEN)='0000'


!     RESET CERTAIN FLAGS:

      IF (IFRELFMAX <= 0) DELPRO_LF = REAL(IDELPRO, JWRB)

!     WE SHOULD RECEIVE DATA FROM NEMO
      IF (LWNEMOCOUCIC .OR. LWNEMOCOUCIT .OR. LWNEMOCOUCUR) LWNEMOCOURECV = .TRUE.


! Here we set LL1D = .TRUE. for the case of LLUNSTR in order to omit the mapping for the parallel strucutured grid
      IF (LLUNSTR) LL1D = .TRUE.
      IF (LWCOU) LSMSSIG_WAM =.FALSE. ! by definition

!     Most of the namelist selection will be written to the logfiles in userin.

!     Some are printed below

      IF (IRANK == 1) THEN
        WRITE(6,*) '==============================================='
        WRITE(6,*) '*** MPUSERIN has read the following settings'
        WRITE(6,*) '*** NANG=',NANG
        WRITE(6,*) '*** NFRE=',NFRE
        WRITE(6,*) '*** NFRE_RED=',NFRE_RED
        WRITE(6,*) '*** CBPLTDT=',CBPLTDT
        WRITE(6,*) '*** CEPLTDT=',CEPLTDT
        WRITE(6,*) '*** LFDB = ',LFDB
        WRITE(6,*) '*** LRSTPARALW = ',LRSTPARALW
        WRITE(6,*) '*** LRSTPARALR = ',LRSTPARALR
        WRITE(6,*) '*** LRSTINFDAT = ',LRSTINFDAT
        WRITE(6,*) '*** LWCOU= ',LWCOU
        WRITE(6,*) '*** LWNEMOCOU= ',LWNEMOCOU
        WRITE(6,*) '*** LWNEMOCOUSEND  = ',LWNEMOCOUSEND
        WRITE(6,*) '*** LWNEMOCOURECV  = ',LWNEMOCOURECV
        WRITE(6,*) '*** LWNEMOCOUDEBUG = ',LWNEMOCOUDEBUG
        WRITE(6,*) '*** LWNEMOCOUCIC   = ',LWNEMOCOUCIC
        WRITE(6,*) '*** LWNEMOCOUCIT   = ',LWNEMOCOUCIT
        WRITE(6,*) '*** LWNEMOCOUCUR   = ',LWNEMOCOUCUR
        WRITE(6,*) '*** LWNEMOCOUSTK   = ',LWNEMOCOUSTK
        WRITE(6,*) '*** LWNEMOCOUSTRN  = ',LWNEMOCOUSTRN
        WRITE(6,*) '*** LWNEMOTAUOC    = ',LWNEMOTAUOC
        WRITE(6,*) '*** LSUBGRID= ',LSUBGRID
        WRITE(6,*) '*** IPROPAGS= ',IPROPAGS
        WRITE(6,*) '*** IREFRA= ',IREFRA
        WRITE(6,*) '*** LLUNSTR= ',LLUNSTR
        WRITE(6,*) '*** LVECTOR= ',LVECTOR
        WRITE(6,*) '*** IVECTOR= ',IVECTOR
        WRITE(6,*) '*** LPREPROC= ',LPREPROC
        WRITE(6,*) '*** LWCOUNORMS= ',LWCOUNORMS
        WRITE(6,*) '*** LLNORMIFS2WAM= ',LLNORMIFS2WAM
        WRITE(6,*) '*** LLNORMWAM2IFS= ',LLNORMWAM2IFS
        WRITE(6,*) '*** LLNORMWAMOUT= ',LLNORMWAMOUT
        WRITE(6,*) '*** LLNORMWAMOUT_GLOBAL= ',LLNORMWAMOUT_GLOBAL
        WRITE(6,*) '*** CNORMWAMOUT_FILE= ',TRIM(CNORMWAMOUT_FILE)
        WRITE(6,*) '*** LSMSSIG_WAM= ',LSMSSIG_WAM
        WRITE(6,*) '*** LWAM_USE_IO_SERV = ',LWAM_USE_IO_SERV
        WRITE(6,'(2A)') ' *** CPATH = ',CPATH(1:ICPLEN)
        WRITE(6,*) '*** LOUTMDLDCP = ',LOUTMDLDCP
        WRITE(6,*) '*** NLOCGRB = ',NLOCGRB
        WRITE(6,*) '*** YEXPVER = ',YEXPVER
        WRITE(6,*) '*** YCLASS = ',YCLASS
        WRITE(6,*) '*** ISTREAM = ',ISTREAM
        IF (NGOUT > 0) THEN
          WRITE (6,*) " OUTPUT POINTS FOR SPECTRA AS DEFINED BY USER INPUT    NO.    LAT.   LONG. "
          DO IC=1,NGOUT
            WRITE (6,'(3X,I5,2F8.2)') IC,OUTLAT(IC),OUTLONG(IC)
          ENDDO
        ENDIF
        WRITE(6,*) '==============================================='
      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPUSERIN',1,ZHOOK_HANDLE)

      END SUBROUTINE MPUSERIN
