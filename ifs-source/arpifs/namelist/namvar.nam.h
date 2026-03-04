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

!     ------------------------------------------------------------------
NAMELIST/NAMVAR/ &
! logical:
            &LTEST,LZOWA,&
            &LJC,LJCDFI,LUSEJCDFI,&
            &LGRASCAL,LWRIBVEC,LWRIBVEC_FULL,LREABVEC,LWRIEVEC,LWRISIGB,LWRISIGA,&
            &LWRISIGF,LEVECCNTL,LEVECGRIB,&
            &L_CHECK_CONVERGENCE,LBGTRUNC,LBGOBS,LBGM,LANOBS,&
            &LAVCGL, LMPCGL, LINITCV, LENSCV, LTOVSCV, LTOVSREP, LSLADREP, LGCV,&
            &LWREINI, LN1CG1, LCONGRAD, L3DFGAT, &
            &LAMV_REASSIGN_PASSIVE, LSKIPMIN, LREPRO4DVAR,LREO3_BCOR,LCH4_BCOR, &
            &LAEOLUSAMD, L_ABS_CONVERGENCE, L_INFO_CONVERGENCE, LAMV_HEIGHT_ADJUST, &
            &L_GUESS_RUNTIME, &
            &L_INFO_CONTENT,LMODERR, LSTATMERR, LVARBC, LPERTMRPAS, &
            &LPROPTL, LCLDSINK, L_CHECK_GRADIENT,LBGPERT,&
            &LBACKGE, LBACKGECV,LBACKGERENORM,LUSEWAVRENORM,LCONSTANTZFCE, &
            &LWRISB_VPROF, &
            &LJBIMPACT, LJCIMPACT, LFCOBS, LFCOBSTEST, LREFINC, LENDA, LMONITOR_FCDEPAR, &
            &LTOY42, LSUSPQLIM, LSPINT, LDIAG_LCT, LFAMEMBERS, &
            &LJBZERO, LJPZERO, LJHZERO, LJLZERO, LJTZERO, &
            &LCHRESINCR, LINC_TOVSCV, LUSE_EDA_SKT, LECV, L_FGNOTBG, L_BGREC, L_FGREC, L_COMPBGDEP, LWRICV, &
! integer:
            &NBGTRUNC, NITER ,NITER_MIN, NSIMU ,NINFRA,NMIMP,NFRREF, NREFTS,&
            &NFRANA,NANATS,NFRGRA,NGRATS,&
            &NUPTRA, MUPTRA, NPCVECS, NBGVECS, NFGFCLEN, NWRIEVEC,&
            &NSAVEEV, NSAVEPC, NITERGCV, N1IMP, NPRECO, NSELECT,&
            &N_INFO_CONTENT_METHOD, N_INFO_CONTENT_SEED,&
            &N_DIAGS_CONVERGENCE, N_DIAGS_EIGENVECS,N_COUPLED_WINDOWS,NDATE_TIME_WINDOW_END,&
            &NITER_GUESS_RUNTIME, NFREQ_GUESS_RUNTIME, NMEM_GUESS_RUNTIME, NUPTRA_RANGE, NFGREC_MIN,&
! real:
            &RCVGE, RDX, ALPHAG, ALPHAV, RXMIN,  &
            &R_NORM_REDUCTION_ABORT_LEVEL, R_MAX_CNUM_PC,&
            &ZEPSNEG, RTOL_CHECK_GRADIENT, FILTERFACTOR, FILTEREXPO, DELTA, &
            &FILTERRESOL, CTOPBGE , CAMTBGE , &
            &RCOEFCO, RCOEFNO2, RCOEFGO3, &
! character:
            &CFNSIGB, CFNSIGA, CFNSIGF, CFNPCV

!     ------------------------------------------------------------------
