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

!*** Note: MAKE SURE ALL VARIABLES ARE IN lowercase !!
!    Please make sure that variable lengths do not exceed 31 chars (a F90-std limit) !
!    If you want to know what is what, please look at odb/cma2odb/initmdb.F90
CHARACTER(LEN=128) :: cactiveretr
logical :: lactiveretr
logical :: allocated_satpred
!***  Observational array lengths
INTEGER(KIND=JPIM) :: nrows_robhdr
INTEGER(KIND=JPIM) :: ncols_robhdr
INTEGER(KIND=JPIM) :: nrows_satpred
INTEGER(KIND=JPIM) :: ncols_satpred
INTEGER(KIND=JPIM) :: nrows_satbody
INTEGER(KIND=JPIM) :: ncols_satbody
INTEGER(KIND=JPIM) :: nrows_sathdr
INTEGER(KIND=JPIM) :: ncols_sathdr
INTEGER(KIND=JPIM) :: nrows_robody
INTEGER(KIND=JPIM) :: ncols_robody
INTEGER(KIND=JPIM) :: nrows_robsu
INTEGER(KIND=JPIM) :: ncols_robsu
INTEGER(KIND=JPIM) :: nrows_robddr
INTEGER(KIND=JPIM) :: ncols_robddr
!***  observational arrays
REAL(KIND=JPRD)   , pointer :: robhdr(:,:)
REAL(KIND=JPRD)   , pointer :: satpred(:,:)
REAL(KIND=JPRD)   , pointer :: satbody(:,:)
REAL(KIND=JPRD)   , pointer :: sathdr(:,:) ! sometimes TARGET
REAL(KIND=JPRD)   , pointer :: robody(:,:)
!***  for obs-setup retrievals obs_boxes, obatabs, mkglobstab, suobarea, ecset, smtov
REAL(KIND=JPRD),    pointer :: robsu(:,:)
!***  link between header & body (the length is normally NROWS_ROBHDR + 1)
INTEGER(KIND=JPIM), pointer :: mlnkh2b(:)
!***  The "JOBS" of the given "JBODY" (the length is normally NROWS_ROBODY)
INTEGER(KIND=JPIM), pointer :: mbodyjobs(:)

!***  ddrs
REAL(KIND=JPRD), pointer :: robddr(:,:)

!*** auxiliary ROBAUX<n>, up to JPMAX_AUX_CASES
INTEGER(KIND=JPIM) :: nrows_robaux(JPMAX_AUX_CASES)
INTEGER(KIND=JPIM) :: ncols_robaux(JPMAX_AUX_CASES)
TYPE(robaux_t) :: iml(JPMAX_AUX_CASES) ! iml == InterMediateLayer

!*** dynamic mdb-pointers

! Aeolus_hdr, scalars
INTEGER(KIND=JPIM) :: mdb_aeo_hdrflag_at_aeolus_hdr !  'aeo_hdrflag@aeolus_hdr'


! Aeolus_hdr, links
! 
INTEGER(KIND=JPIM) :: mlnk_sat2aeolus_hdr(2)
INTEGER(KIND=JPIM) :: mlnk_aeolus_hdr2aeolus_auxmet(2)
INTEGER(KIND=JPIM) :: mlnk_aeolus_hdr2aeolus_l2c(2)

! Aeolus_auxmet, all scalars
INTEGER(KIND=JPIM) :: mdb_lev_at_aeolus_auxmet     !  'lev@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_ptop_at_aeolus_auxmet    !  'ptop@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_pnom_at_aeolus_auxmet    !  'pnom@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_ztop_at_aeolus_auxmet    !  'ztop@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_znom_at_aeolus_auxmet    !  'znom@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_u_at_aeolus_auxmet       !  'u@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_v_at_aeolus_auxmet       !  'v@aeolus_auxmet'
!
INTEGER(KIND=JPIM) :: mdb_t_at_aeolus_auxmet       !  't@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_rh_at_aeolus_auxmet      !  'rh@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_q_at_aeolus_auxmet       !  'q@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_cc_at_aeolus_auxmet      !  'cc@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_clwc_at_aeolus_auxmet    !  'clwc@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_ciwc_at_aeolus_auxmet    !  'ciwc@aeolus_auxmet'
!
INTEGER(KIND=JPIM) :: mdb_error_t_at_aeolus_auxmet       !  'error_t@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_error_rh_at_aeolus_auxmet      !  'error_rh@aeolus_auxmet'
INTEGER(KIND=JPIM) :: mdb_error_p_at_aeolus_auxmet       !  'error_p@aeolus_auxmet'

! Aeolus_l2c, scalars
INTEGER(KIND=JPIM) :: mdb_hlos_ob_err_at_aeolus_l2c   !  'hlos_ob_err@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_hlos_fg_at_aeolus_l2c       !  'hlos_fg@aeolus_l2c'   
INTEGER(KIND=JPIM) :: mdb_u_fg_at_aeolus_l2c          !  'u_fg@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_u_fg_err_at_aeolus_l2c      !  'u_fg_err@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_v_fg_at_aeolus_l2c          !  'v_fg@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_v_fg_err_at_aeolus_l2c      !  'v_fg_err@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_hlos_fg_err_at_aeolus_l2c   !  'hlos_fg_err@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_hlos_an_at_aeolus_l2c       !  'hlos_an@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_hlos_an_err_at_aeolus_l2c   !  'hlos_an_err@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_u_an_at_aeolus_l2c          !  'u_an@aeolus_l2c'
INTEGER(KIND=JPIM) :: mdb_v_an_at_aeolus_l2c          !  'v_an@aeolus_l2c'


! enkf scalars
INTEGER(KIND=JPIM) :: mdb_member_at_enkf(jpmxenkf)   ! 'member@enkf[$NENKF]'
INTEGER(KIND=JPIM) :: mdb_hprior_at_enkf(jpmxenkf)   ! 'hprior@enkf[$NENKF]'

! enda columns
INTEGER(KIND=JPIM) :: mdb_member_at_enda(jpmxenda)   ! 'member@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_report_status_at_enda(jpmxenda)   ! 'report_status@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_report_event1_at_enda(jpmxenda)   ! 'report_event1@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_obsvalue_at_enda(jpmxenda)   ! 'obsvalue@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_anflag_at_enda(jpmxenda)   ! 'datum_anflag@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_status_at_enda(jpmxenda)   ! 'datum_status@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_event1_at_enda(jpmxenda)   ! 'datum_event1@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_biascorr_at_enda(jpmxenda)   ! 'biascorr@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_biascorr_fg_at_enda(jpmxenda)   ! 'biascorr_fg@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_an_depar_at_enda(jpmxenda)   ! 'an_depar@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_fg_depar_at_enda(jpmxenda)   ! 'fg_depar@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_qc_pge_at_enda(jpmxenda)   ! 'qc_pge@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_final_obs_error_at_enda(jpmxenda)   ! 'final_obs_error@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_obs_error_at_enda(jpmxenda)   ! 'obs_error@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_fg_error_at_enda(jpmxenda)   ! 'fg_error@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_skintemp_at_enda(jpmxup+1,jpmxenda)   ! 'skintemp[$NMXUPD]@enda[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_status_hires_at_enda(jpmxenda)   ! 'datum_status_hires@enda[$NENDA]'

! enda columns for surface analysis feedback
INTEGER(KIND=JPIM) :: mdb_fg_depar_at_sfc_fb(jpmxenda)   ! 'fg_depar@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_an_depar_at_sfc_fb(jpmxenda)   ! 'an_depar@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_sfc_event_at_sfc_fb(jpmxenda)   ! 'datum_sfc_event@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_datum_status_at_sfc_fb(jpmxenda)   ! 'datum_status@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_lsm_at_sfc_fb(jpmxenda)   ! 'lsm@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_snow_depth_at_sfc_fb(jpmxenda)   ! 'snow_depth@surfbody_feedback[$NENDA]'
INTEGER(KIND=JPIM) :: mdb_snow_density_at_sfc_fb(jpmxenda)   ! 'snow_density@surfbody_feedback[$NENDA]'


INTEGER(KIND=JPIM) :: mdb_albedo  ! 'albedo@modsurf' 
INTEGER(KIND=JPIM) :: mdb_andate_at_desc  ! 'andate@desc' 
INTEGER(KIND=JPIM) :: mdb_antime_at_desc  ! 'antime@desc' 
INTEGER(KIND=JPIM) :: mdb_inidate_at_desc  ! 'inidate@desc' 
INTEGER(KIND=JPIM) :: mdb_initime_at_desc  ! 'initime@desc' 
INTEGER(KIND=JPIM) :: mdb_expver_at_desc ! 'expver@desc' 
INTEGER(KIND=JPIM) :: mdb_class_at_desc  ! 'class@desc'
INTEGER(KIND=JPIM) :: mdb_stream_at_desc ! 'stream@desc'
INTEGER(KIND=JPIM) :: mdb_type_at_desc   ! 'type@desc'
INTEGER(KIND=JPIM) :: mdbonm ! 'seqno@hdr'
INTEGER(KIND=JPIM) :: mdbotp ! 'obstype@hdr'
!INTEGER_M :: mdboch ! 'obschar@hdr'
INTEGER(KIND=JPIM) :: mdbdat ! 'date@hdr'
INTEGER(KIND=JPIM) :: mdbetm ! 'time@hdr'
INTEGER(KIND=JPIM) :: mdbrfl ! 'report_rdbflag@hdr'
INTEGER(KIND=JPIM) :: mdbrst ! 'report_status@hdr'
INTEGER(KIND=JPIM) :: mdbrev1 ! 'report_event1@hdr'
INTEGER(KIND=JPIM) :: mdbrble ! 'report_blacklist@hdr'
INTEGER(KIND=JPIM) :: mdbbox ! 'sortbox@hdr'
INTEGER(KIND=JPIM) :: mdbstd ! 'sitedep@hdr'
INTEGER(KIND=JPIM) :: mdbsid ! 'statid@hdr'
INTEGER(KIND=JPIM) :: mdblat ! 'lat@hdr'
INTEGER(KIND=JPIM) :: mdblon ! 'lon@hdr'
INTEGER(KIND=JPIM) :: mdbalt ! 'stalt@hdr'
INTEGER(KIND=JPIM) :: mdb_lsm ! 'lsm@modsurf'
INTEGER(KIND=JPIM) :: mdb_seaice ! 'seaice@modsurf'
INTEGER(KIND=JPIM) :: mdb_orography ! 'orography@modsurf'
INTEGER(KIND=JPIM) :: mdb_snow_depth 
INTEGER(KIND=JPIM) :: mdb_snow_density 
INTEGER(KIND=JPIM) :: mdb_soil_moisture
INTEGER(KIND=JPIM) :: mdb_soil_temp
INTEGER(KIND=JPIM) :: mdb_lai
INTEGER(KIND=JPIM) :: mdb_t2m ! 't2m@modsurf'
INTEGER(KIND=JPIM) :: mdb_windspeed10m ! 'windspeed10m@modsurf'
INTEGER(KIND=JPIM) :: mdb_u10m ! 'u10m@modsurf'
INTEGER(KIND=JPIM) :: mdb_v10m ! 'v10m@modsurf'
INTEGER(KIND=JPIM) :: mdb_window_offset ! 'window_offset@hdr'
INTEGER(KIND=JPIM) :: mdb_surface_class ! 'surface_class@modsurf'
INTEGER(KIND=JPIM) :: mdbtla ! 'trlat@hdr'
INTEGER(KIND=JPIM) :: mdbtlo ! 'trlon@hdr'
INTEGER(KIND=JPIM) :: mdbrev2 ! 'report_event2@hdr'
INTEGER(KIND=JPIM) :: mdbsbcmm ! 'comp_method@satob'
INTEGER(KIND=JPIM) :: mdbsbiup ! 'instdata@satob'
INTEGER(KIND=JPIM) :: mdbsbdpt ! 'dataproc@satob'
INTEGER(KIND=JPIM) :: mdb_qi_fc_at_satob ! 'qi_fc@satob'
INTEGER(KIND=JPIM) :: mdb_rff_at_satob ! 'rff@satob'
INTEGER(KIND=JPIM) :: mdb_qi_nofc_at_satob ! 'qi_nofc@satob'
INTEGER(KIND=JPIM) :: mdb_ee_at_satob ! 'ee@satob'
INTEGER(KIND=JPIM) ::   mdb_tb_at_satob ! 'tb@satob' 
INTEGER(KIND=JPIM) ::   mdb_t_at_satob ! 't@satob' 
INTEGER(KIND=JPIM) ::   mdb_shear_at_satob ! 'shear@satob' 
INTEGER(KIND=JPIM) ::   mdb_t200_at_satob  ! 't200@satob' 
INTEGER(KIND=JPIM) ::   mdb_t500_at_satob ! 't500@satob' 
INTEGER(KIND=JPIM) ::   mdb_top_mean_t_at_satob ! 'top_mean_t@satob' 
INTEGER(KIND=JPIM) ::   mdb_top_wv_at_satob ! 'top_wv@satob' 
INTEGER(KIND=JPIM) ::   mdb_dt_by_dp_at_satob ! 'dt_by_dp@satob' 
INTEGER(KIND=JPIM) ::   mdb_p_best_at_satob ! 'p_best@satob' 
INTEGER(KIND=JPIM) ::   mdb_u_best_at_satob ! 'u_best@satob' 
INTEGER(KIND=JPIM) ::   mdb_v_best_at_satob ! 'v_best@satob' 
INTEGER(KIND=JPIM) ::   mdb_p_old_at_satob ! 'p_old@satob' 
INTEGER(KIND=JPIM) ::   mdb_u_old_at_satob ! 'u_old@satob' 
INTEGER(KIND=JPIM) ::   mdb_v_old_at_satob ! 'v_old@satob' 
INTEGER(KIND=JPIM) ::   mdb_height_assignment_method_at_satob ! 'height_assignment_method@satob'
INTEGER(KIND=JPIM) ::   mdb_tracer_correlation_method_at_satob ! 'tracer_correlation_method@satob
INTEGER(KIND=JPIM) ::   mdb_land_sea_at_satob ! 'land_sea@satob
INTEGER(KIND=JPIM) ::   mdb_tracking_error_u_at_satob ! 'tracking_error_u@satob'
INTEGER(KIND=JPIM) ::   mdb_tracking_error_v_at_satob ! 'tracking_error_v@satob'
INTEGER(KIND=JPIM) ::   mdb_h_assignment_error_u_at_satob ! 'h_assignment_error_u@satob'
INTEGER(KIND=JPIM) ::   mdb_h_assignment_error_v_at_satob  ! 'h_assignment_error_v@satob'
INTEGER(KIND=JPIM) ::   mdb_error_in_h_assignment_at_satob  ! 'error_in_h_assignment@satob'
INTEGER(KIND=JPIM) ::   mdb_ct_p_at_satob  ! 'ct_p@satob'
INTEGER(KIND=JPIM) ::   mdb_cb_p_at_satob  ! 'cb_p@satob'
INTEGER(KIND=JPIM) ::   mdb_umod_old_at_satob  ! 'umod_old@satob'
INTEGER(KIND=JPIM) ::   mdb_vmod_old_at_satob  ! 'vmod_old@satob'
INTEGER(KIND=JPIM) :: mdbssia ! 'sensor@hdr'
INTEGER(KIND=JPIM) :: mdbsccno ! 'cellno@scatt'
INTEGER(KIND=JPIM) :: mdbscpfl ! 'prodflag@scatt'
!%INTEGER(KIND=JPIM) :: mdb1derr ! 'rainrate_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1la ! 'lat_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1lo ! 'lon_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1tq ! 'timequal_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1da ! 'date_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1ti ! 'time_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1ts ! 'surftype_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1vs ! 'vertsign_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1bt1 ! 'temp1_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe1bt2 ! 'temp2_1@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2la ! 'lat_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2lo ! 'lon_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2tq ! 'timequal_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2da ! 'date_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2ti ! 'time_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2ts ! 'surftype_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2vs ! 'vertsign_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2bt1 ! 'temp1_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe2bt2 ! 'temp2_2@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3la ! 'lat_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3lo ! 'lon_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3tq ! 'timequal_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3da ! 'date_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3ti ! 'time_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3ts ! 'surftype_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3vs ! 'vertsign_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3bt1 ! 'temp1_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbe3bt2 ! 'temp2_3@ssmi'
!%INTEGER(KIND=JPIM) :: mdbssmst ! 'surftype@ssmi'
!%INTEGER(KIND=JPIM) :: mdb1dstq ! 'qc_surftype@ssmi'
INTEGER(KIND=JPIM) :: mdb1dnit ! 'iterno_conv_1dvar@ssmi'
INTEGER(KIND=JPIM) :: mdb1dnis ! 'simno_conv_1dvar@ssmi'
INTEGER(KIND=JPIM) :: mdb1deps ! 'epsg_1dvar@ssmi'
INTEGER(KIND=JPIM) :: mdb1dfin ! 'failure_1dvar@ssmi'
INTEGER(KIND=JPIM) :: mdb1dfim ! 'minim_status_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdbies ! 'scattering_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdbssisi ! 'scattering_indep@ssmi'
!%INTEGER(KIND=JPIM) :: mdbsserr ! 'rainrate_indep@ssmi'
!%INTEGER(KIND=JPIM) :: mdb1drep ! 'tpw_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdb1drew ! 'wspeed_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdb1drec ! 'clwpath_1dvar@ssmi'
!%INTEGER(KIND=JPIM) :: mdbierep ! 'tpw_indep@ssmi'
!%INTEGER(KIND=JPIM) :: mdbierew ! 'wspeed_indep@ssmi'
!%INTEGER(KIND=JPIM) :: mdbierec ! 'clwpath_indep@ssmi'
INTEGER(KIND=JPIM) :: mdb1bps ! 'surfpress_1@ssmi' !c
INTEGER(KIND=JPIM) :: mdb1bsts ! 'skintemp_1@ssmi' !c
!%INTEGER(KIND=JPIM) :: mdb1b2ts ! 't2m_1@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1b2qs ! 'q2m_1@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1bpws ! 'pwc_1@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1bcls ! 'clw_1@ssmi_slev'
INTEGER(KIND=JPIM) :: mdb1bsus(2) ! 'u10m[2]@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1bsvs(2) ! 'v10m[2]@ssmi' !n
!%INTEGER(KIND=JPIM) :: mdb1a2qs ! 'q2m_2@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1apws ! 'pwc_2@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1acls ! 'clw_2@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdb1asws ! 'wspeed10m_2@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdbsipws ! 'indep_pwc@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdbsicls ! 'indep_clw@ssmi_slev'
!%INTEGER(KIND=JPIM) :: mdbsisws ! 'indep_wspeed10m@ssmi_slev'
INTEGER(KIND=JPIM) :: mdb1dvps ! 'press@ssmi_body' !c
INTEGER(KIND=JPIM) :: mdb1dvbt ! 'temp_1@ssmi_body' !c
INTEGER(KIND=JPIM) :: mdb1dvat ! 'temp_2@ssmi_body' !c
INTEGER(KIND=JPIM) :: mdb1dvbq ! 'q_1@ssmi_body' !c
INTEGER(KIND=JPIM) :: mdb1dvaq ! 'q_2@ssmi_body' !c
!%INTEGER(KIND=JPIM) :: mdb1dvbc ! 'clw_1@ssmi_mlev'
!%INTEGER(KIND=JPIM) :: mdb1dvac ! 'clw_2@ssmi_mlev'
INTEGER(KIND=JPIM) :: mdb_procid_at_index  ! 'procid@index' 
INTEGER(KIND=JPIM) :: mdb_target_at_index  ! 'target@index' 
INTEGER(KIND=JPIM) :: mdb_distribid_at_hdr  ! 'distribid@hdr'
INTEGER(KIND=JPIM) :: mdb_distribtype_at_hdr  ! 'distribtype@hdr'
INTEGER(KIND=JPIM) :: mdb_gp_dist ! 'gp_dist@hdr'
INTEGER(KIND=JPIM) :: mdb_gp_number ! 'gp_number@hdr'
INTEGER(KIND=JPIM) :: mdb_kset_at_index  ! 'kset@index' 
INTEGER(KIND=JPIM) :: mdb_tslot_at_index  ! 'timeslot@index' 
INTEGER(KIND=JPIM) :: mdb_abnob_at_index  ! 'abnob@index' 
INTEGER(KIND=JPIM) :: mdb_mapomm_at_index ! 'mapomm@index'
INTEGER(KIND=JPIM) :: mdb_maptovscv_at_index ! 'maptovscv@index'
INTEGER(KIND=JPIM) :: mdb_thinningkey_at_hdr(jp_numthbox) ! 'thinningkey@hdr' 
INTEGER(KIND=JPIM) :: mdb_thinningtimekey_at_hdr ! 'thinningtimekey@hdr' 
INTEGER(KIND=JPIM) :: mlnk_index2hdr(2) 
INTEGER(KIND=JPIM) :: mlnk_hdr2body(2) 
INTEGER(KIND=JPIM) :: mlnk_hdr2update(2,jpmxup) 
INTEGER(KIND=JPIM) :: mlnk_hdr2ensemble(2) 
INTEGER(KIND=JPIM) :: mlnk_ensemble2enkf(2,jpmxenkf) 
INTEGER(KIND=JPIM) :: mlnk_ensemble2enda(2,jpmxenda) 
INTEGER(KIND=JPIM) :: mlnk_ensemble2surfbody_feedback(2,jpmxenda) 
INTEGER(KIND=JPIM) :: mlnk_desc2fcdiag(2)
INTEGER(KIND=JPIM) :: mlnk_fcdiag2fcdiag_body(2,jpmxfcdiag)
INTEGER(KIND=JPIM) :: mlnk_hdr2errstat(2) 
INTEGER(KIND=JPIM) :: mlnk_hdr2sat(2) 
INTEGER(KIND=JPIM) :: mlnk_hdr2superob(2)
INTEGER(KIND=JPIM) :: mlnk_sat2limb(2) 
INTEGER(KIND=JPIM) :: mlnk_sat2resat(2)
INTEGER(KIND=JPIM) :: mlnk_sat2smos(2)
INTEGER(KIND=JPIM) :: mlnk_sat2ssmi(2) 
INTEGER(KIND=JPIM) :: mlnk_sat2scatt(2) 
INTEGER(KIND=JPIM) :: mlnk_sat2satob(2) 
INTEGER(KIND=JPIM) :: mlnk_ssmi2ssmi_body(2) 
INTEGER(KIND=JPIM) :: mlnk_scatt2scatt_body(2) 
!*AF new tables
INTEGER(KIND=JPIM) :: mlnk_radiance2allsky(2)
INTEGER(KIND=JPIM) :: mlnk_allsky2allsky_body(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2gbrad(2)
INTEGER(KIND=JPIM) :: mlnk_gbrad2gbrad_body(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2raingg(2)
INTEGER(KIND=JPIM) :: mlnk_raingg2raingg_body(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2auxiliary(2)
INTEGER(KIND=JPIM) :: mlnk_auxiliary2auxiliary_body(2)
INTEGER(KIND=JPIM) :: mlnk_radiance2cloud_sink(2)
INTEGER(KIND=JPIM) :: mlnk_radiance2colloc_im_info(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2modsurf(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2surfbody_feedback(2)
INTEGER(KIND=JPIM) :: mlnk_sat2radiance(2)
INTEGER(KIND=JPIM) :: mlnk_radiance2radiance_body(2)
INTEGER(KIND=JPIM) :: mlnk_sat2gnssro(2)
INTEGER(KIND=JPIM) :: mlnk_gnssro2gnssro_body(2)
!*AF END

!*** New items for cloud radar and lidar observations
INTEGER(KIND=JPIM) :: mlnk_sat2clrad(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2clrad(2)
INTEGER(KIND=JPIM) :: mlnk_clrad2clrad_body(2)
INTEGER(KIND=JPIM) :: mlnk_sat2claerlid(2)
INTEGER(KIND=JPIM) :: mlnk_hdr2claerlid(2)
INTEGER(KIND=JPIM) :: mlnk_claerlid2claerlid_body(2)

INTEGER(KIND=JPIM) :: mdb_channel_qc_at_radiance_body ! 'channel_qc@radiance_body'
INTEGER(KIND=JPIM) :: mdb_cold_nedt_at_radiance_body ! 'cold_nedt@radiance_body'
INTEGER(KIND=JPIM) :: mdb_warm_nedt_at_radiance_body ! 'warm_nedt@radiance_body'
INTEGER(KIND=JPIM) :: mdbvnm ! 'varno@body'
INTEGER(KIND=JPIM) :: mdbvco ! 'vertco_type@body'
INTEGER(KIND=JPIM) :: mdbrdfl ! 'datum_rdbflag@body'
INTEGER(KIND=JPIM) :: mdbflg ! 'datum_anflag@body'
INTEGER(KIND=JPIM) :: mdbdsta ! 'datum_status@body'
INTEGER(KIND=JPIM) :: mdbdev1 ! 'datum_event1@body'
INTEGER(KIND=JPIM) :: mdb_sfc_event ! 'datum_sfc_event@body'
INTEGER(KIND=JPIM) :: mdbdble ! 'datum_blacklist@body'
INTEGER(KIND=JPIM) :: mdbesqn ! 'entryno@body'
INTEGER(KIND=JPIM) :: mdbppp ! 'vertco_reference_1@body'
INTEGER(KIND=JPIM) :: mdbprl ! 'vertco_reference_2@body'
INTEGER(KIND=JPIM) :: mdbvar ! 'obsvalue@body'
INTEGER(KIND=JPIM) :: mdbomn ! 'an_depar@body'
INTEGER(KIND=JPIM) :: mdbomf ! 'fg_depar@body'
INTEGER(KIND=JPIM) :: mdbaso ! 'an_sens_obs@body'
INTEGER(KIND=JPIM) :: mdbfso ! 'fc_sens_obs@body'
INTEGER(KIND=JPIM) :: mdbfoe ! 'final_obs_error@errstat'
INTEGER(KIND=JPIM) :: mdboer ! 'obs_error@errstat'
INTEGER(KIND=JPIM) :: mdbrer ! 'repres_error@errstat'
INTEGER(KIND=JPIM) :: mdbper ! 'pers_error@errstat'
INTEGER(KIND=JPIM) :: mdbfge ! 'fg_error@errstat'
INTEGER(KIND=JPIM) :: mdb_eda_spread ! 'eda_spread@errstat'
INTEGER(KIND=JPIM) :: mdb_standard_deviation_at_superob ! 'standard_deviation@superob'
INTEGER(KIND=JPIM) :: mdb_n_obs_at_superob ! 'n_obs@superob'
INTEGER(KIND=JPIM) :: mdb_cloud_fraction_at_superob ! 'cloud_fraction@superob'
INTEGER(KIND=JPIM) :: mdb_sf_max_at_superob ! 'sf_max@superob'

INTEGER(KIND=JPIM) :: mdb_actual_depar ! 'actual_depar@body'
INTEGER(KIND=JPIM) :: mdb_actual_ndbiascorr ! 'actual_ndbiascorr@body'
INTEGER(KIND=JPIM) :: mdb_obscorev_at_errstat(jp_numev) ! 'obs_corr_ev[$NUMEV]@errstat'
INTEGER(KIND=JPIM) :: mdb_obscormask_at_errstat ! 'obs_corr_mask@errstat'
INTEGER(KIND=JPIM) :: mdb_qc_a! 'qc_a@body'
INTEGER(KIND=JPIM) :: mdb_qc_l ! 'qc_l@body'
INTEGER(KIND=JPIM) :: mdb_qc_pge ! 'qc_pge@body'
INTEGER(KIND=JPIM) ::    mdbifc1(jpmxup) ! hires@update[$NMXUPD]
INTEGER(KIND=JPIM) ::    mdbifc2(jpmxup) ! lores@update[$NMXUPD]
INTEGER(KIND=JPIM) ::    mdb_datum_tbflag_hires(jpmxup) ! datum_tbflag_hires@update[$NMXUPD]
INTEGER(KIND=JPIM) ::    mdb_datum_status_hires(jpmxup) ! datum_status_hires@update[$NMXUPD]
INTEGER(KIND=JPIM) :: mdbrbvc ! 'mf_vertco_type@body'
INTEGER(KIND=JPIM) :: mdbrbpio ! 'mf_log_p@body'
INTEGER(KIND=JPIM) :: mdbrboe ! 'mf_stddev@body'
INTEGER(KIND=JPIM) :: mdbdev2 ! 'event2@body'
!-- Due to the faster MAPVARDB, a column name must be mapped only once to a MDB-pointer
INTEGER(KIND=JPIM) :: mdbtorb ! 'biascorr@body'
!!INTEGER(KIND=JPIM) :: mdbtorba ! 'biascorr@body'
!!INTEGER(KIND=JPIM) :: mdbssrb ! 'biascorr@body'
INTEGER(KIND=JPIM) :: mdbs1dvc ! 'radcost@ssmi_body'
!%INTEGER(KIND=JPIM) :: mdbssrde ! 'depar@ssmi_body'
INTEGER(KIND=JPIM) :: mdbssccf ! 'frequency@ssmi_body'
INTEGER(KIND=JPIM) :: mdbsscbw ! 'bandwidth@ssmi_body'
INTEGER(KIND=JPIM) :: mdbssanp ! 'polarisation@ssmi_body'
INTEGER(KIND=JPIM) :: mdbscbaa ! 'azimuth@scatt_body'
INTEGER(KIND=JPIM) :: mdbscbia ! 'incidence@scatt_body'
INTEGER(KIND=JPIM) :: mdbsckp ! 'kp@scatt_body'
INTEGER(KIND=JPIM) :: mdbscres ! 'invresid@scatt_body'
INTEGER(KIND=JPIM) :: mdbscdis ! 'dirskill@scatt_body'
INTEGER(KIND=JPIM) :: mdbsckpq ! 'kp_qf@scatt_body'
INTEGER(KIND=JPIM) :: mdbscs0q ! 'sigma0_qf@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsm  ! 'sigma0_sm@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsms ! 'soilmoist_sd@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsmc ! 'soilmoist_cf@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsmp ! 'soilmoist_pf@scatt_body'
INTEGER(KIND=JPIM) :: mdbsclfr ! 'land_fraction@scatt_body'
INTEGER(KIND=JPIM) :: mdb_likelihood_at_scatt_body ! 'likelihood@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsmw ! 'wetland_fraction@scatt_body'
INTEGER(KIND=JPIM) :: mdbscsmss ! 'sm_sensitivity@scatt_body'


INTEGER(KIND=JPIM) :: mdbscsmt ! 'topo_complex@scatt_body'
INTEGER(KIND=JPIM) :: mdb_satid_at_sat  ! 'satellite_identifier@sat'
INTEGER(KIND=JPIM) :: mdb_satinst_at_sat ! 'satellite_instrument@sat'
INTEGER(KIND=JPIM) :: mdb_gen_centre_at_sat  ! 'gen_centre@sat'
INTEGER(KIND=JPIM) :: mdb_gen_subcentre_at_sat  ! 'gen_subcentre@sat'
INTEGER(KIND=JPIM) :: mdb_datastream_at_sat      ! 'datastream@sat'
INTEGER(KIND=JPIM) :: mdb_channel_at_sat         ! 'channel@sat'
INTEGER(KIND=JPIM) ::   mdb_cldptop1_at_radiance  ! 'cldptop1@radiance' 
INTEGER(KIND=JPIM) ::   mdb_cldne1_at_radiance  ! 'cldne1@radiance' 
INTEGER(KIND=JPIM) ::   mdb_cldptop2_at_radiance  ! 'cldptop2@radiance' 
INTEGER(KIND=JPIM) ::   mdb_cldne2_at_radiance  ! 'cldne2@radiance' 
INTEGER(KIND=JPIM) ::   mdb_cldptop3_at_radiance  ! 'cldptop3@radiance' 
INTEGER(KIND=JPIM) ::   mdb_cldne3_at_radiance  ! 'cldne3@radiance' 
INTEGER(KIND=JPIM) ::   mdb_nobs_averaged  ! 'nobs_averaged@radiance_body' 
INTEGER(KIND=JPIM) ::   mdb_stdev_averaged  ! 'stdev_averaged@radiance_body'
INTEGER(KIND=JPIM) :: mdbresatit ! 'instrument_type@resat'
INTEGER(KIND=JPIM) :: mdbresatpt ! 'product_type@resat'
INTEGER(KIND=JPIM) :: mdbresatla1 ! 'lat_fovcorner_1@resat'
INTEGER(KIND=JPIM) :: mdbresatlo1 ! 'lon_fovcorner_1@resat'
INTEGER(KIND=JPIM) :: mdbresatla2 ! 'lat_fovcorner_2@resat'
INTEGER(KIND=JPIM) :: mdbresatlo2 ! 'lon_fovcorner_2@resat'
INTEGER(KIND=JPIM) :: mdbresatla3 ! 'lat_fovcorner_3@resat'
INTEGER(KIND=JPIM) :: mdbresatlo3 ! 'lon_fovcorner_3@resat'
INTEGER(KIND=JPIM) :: mdbresatla4 ! 'lat_fovcorner_4@resat'
INTEGER(KIND=JPIM) :: mdbresatlo4 ! 'lon_fovcorner_4@resat'
INTEGER(KIND=JPIM) :: mdbresatsoe ! 'solar_elevation@resat'
INTEGER(KIND=JPIM) :: mdbresatfov ! 'scanpos@resat'
INTEGER(KIND=JPIM) :: mdbresatclc ! 'cloud_cover@resat'
INTEGER(KIND=JPIM) :: mdbresatcp ! 'cloud_top_press@resat'
INTEGER(KIND=JPIM) :: mdbresatqr ! 'quality_retrieval@resat'
INTEGER(KIND=JPIM) :: mdbresatnl ! 'number_layers@resat'
INTEGER(KIND=JPIM) :: mdbaersii ! 'snow_ice_indicator@resat'
INTEGER(KIND=JPIM) :: mdbaerstf ! 'surface_type_indicator@resat'
INTEGER(KIND=JPIM) :: mdb_creadate_at_desc  ! 'creadate@desc' 
INTEGER(KIND=JPIM) :: mdb_creatime_at_desc  ! 'creatime@desc' 
INTEGER(KIND=JPIM) :: mdb_creaby_at_desc    ! 'creaby@desc' 
INTEGER(KIND=JPIM) :: mdb_moddate_at_desc   ! 'moddate@desc' 
INTEGER(KIND=JPIM) :: mdb_modtime_at_desc   ! 'modtime@desc' 
INTEGER(KIND=JPIM) :: mdb_modby_at_desc     ! 'modby@desc' 
INTEGER(KIND=JPIM) :: mlnk_desc2hdr(2) 
INTEGER(KIND=JPIM) :: mdb_subtype_at_hdr   ! 'subtype@hdr' 
INTEGER(KIND=JPIM) :: mdb_bufrtype_at_hdr  ! 'bufrtype@hdr' 
INTEGER(KIND=JPIM) :: mdb_groupid_at_hdr   ! 'groupid@hdr'
INTEGER(KIND=JPIM) :: mdb_reportype_at_hdr ! 'reportype@hdr'
INTEGER(KIND=JPIM) :: mdb_numlev_at_hdr    ! 'numlev@hdr' 
INTEGER(KIND=JPIM) :: mdb_numactiveb_at_hdr ! 'numactiveb@hdr'
INTEGER(KIND=JPIM) :: mdb_duplseqno_at_hdr(jp_maxdupl)  ! 'duplseqno@hdr' 
!*** New items added due to direct BUFR2ODB conversion
INTEGER(KIND=JPIM) :: mdb_mpc_at_scatt_body ! 'mpc@scatt_body' 
INTEGER(KIND=JPIM) :: mdb_wvc_qf ! 'wvc_qf@scatt' 
INTEGER(KIND=JPIM) :: mdb_nretr_amb ! 'nretr_amb@scatt' 
INTEGER(KIND=JPIM) :: mdb_subsetno_at_hdr ! 'subsetno@hdr' 
!*** New items added due to Poolmasking
INTEGER(KIND=JPIM) :: mlnk_desc2poolmask(2)
INTEGER(KIND=JPIM) :: mdb_poolno_at_poolmask ! 'poolno@poolmask' 
INTEGER(KIND=JPIM) :: mdb_obstype_at_poolmask ! 'obstype@poolmask' 
INTEGER(KIND=JPIM) :: mdb_codetype_at_poolmask ! 'codetype@poolmask' 
INTEGER(KIND=JPIM) :: mdb_sensor_at_poolmask ! 'sensor@poolmask' 
INTEGER(KIND=JPIM) :: mdb_tslot_at_poolmask ! 'timeslot@poolmask' 
INTEGER(KIND=JPIM) :: mdb_subtype_at_poolmask ! 'subtype@poolmask' 
INTEGER(KIND=JPIM) :: mdb_bufrtype_at_poolmask ! 'bufrtype@poolmask' 
INTEGER(KIND=JPIM) :: mdb_hdr_count_at_poolmask ! 'hdr_count@poolmask' 
INTEGER(KIND=JPIM) :: mdb_body_count_at_poolmask ! 'body_count@poolmask' 
INTEGER(KIND=JPIM) :: mdb_max_bodylen_at_poolmask ! 'max_bodylen@poolmask' 
!*** New items added due to introduction of 'timeslot_index'-table
INTEGER(KIND=JPIM) :: mdb_timeslot_at_timeslot_index ! 'timeslot@timeslot_index' 
INTEGER(KIND=JPIM) :: mdb_modstep_at_timeslot_index ! 'model_timestep@timeslot_index'
INTEGER(KIND=JPIM) :: mdb_enddate_at_timeslot_index ! 'enddate@timeslot_index' 
INTEGER(KIND=JPIM) :: mdb_endtime_at_timeslot_index ! 'endtime@timeslot_index' 
INTEGER(KIND=JPIM) :: mlnk_desc2timeslot_index(2)
INTEGER(KIND=JPIM) :: mlnk_timeslot_index2index(2)
!*** Commonly needed subwords/members of 'obschar.*@hdr'
!INTEGER_M :: mdb_obschar_d_codetype_at_hdr  ! 'obschar.codetype@hdr' ! 'obschar_d_codetype@hdr' 
!INTEGER_M :: mdb_obschar_d_instype_at_hdr   ! 'obschar.instype@hdr' ! 'obschar_d_instype@hdr' 
!INTEGER_M :: mdb_obschar_d_retrtype_at_hdr  ! 'obschar.retrtype@hdr' ! 'obschar_d_retrtype@hdr' 
!INTEGER_M :: mdb_obschar_d_datagroup_at_hdr ! 'obschar.datagroup@hdr' ! 'obschar_d_datagroup@hdr' 
!*** obschar is replaced by 4 new entries
INTEGER(KIND=JPIM) :: mdb_codetype_at_hdr ! 'codetype@hdr'
INTEGER(KIND=JPIM) :: mdb_insttype_at_hdr ! 'instrument_type@hdr'
INTEGER(KIND=JPIM) :: mdb_retrtype_at_hdr ! 'retrtype@hdr'
INTEGER(KIND=JPIM) :: mdb_areatype_at_hdr ! 'areatype@hdr'
!*** Additional @satob entries requested by Niels Bormann  [20/11/2001]
INTEGER(KIND=JPIM) :: mdb_segment_size_x_at_satob ! 'segment_size_x@satob' 
INTEGER(KIND=JPIM) :: mdb_segment_size_y_at_satob ! 'segment_size_y@satob' 
!*** New @satob table item suggested by Niels Bormann [7/11/2002]
INTEGER(KIND=JPIM) :: mdb_chan_freq_at_satob ! 'chan_freq@satob'
!*** New @body & @atovs entries proposed by Christina Koepken [26/11/2002]
INTEGER(KIND=JPIM) :: mdb_cld_fg_depar ! 'cld_fg_depar@radiance_body'
INTEGER(KIND=JPIM) :: mdb_csr_pclear ! 'csr_pclear@radiance_body'
!*** New @desc & @hdr entries pushed in by Sami Saarinen [17/02/2004]
INTEGER(KIND=JPIM) :: mdb_mxup_traj_at_desc ! 'mxup_traj@desc'
INTEGER(KIND=JPIM) :: mdb_numtsl_at_desc ! 'numtsl@desc'
INTEGER(KIND=JPIM) :: mdb_source_at_hdr     ! 'source@hdr'
!*** VarBC-related @body entries by Dick Dee [23/10/2007]
INTEGER(KIND=JPIM) :: mdb_biascorr_fg_at_body  ! 'biascorr_fg@body'
INTEGER(KIND=JPIM) :: mdb_varbc_ix_at_body  ! 'varbc_ix@body'
!*** Table for ASRs
INTEGER(KIND=JPIM) :: mdb_asr_pclear ! 'asr_pclear@radiance' 
INTEGER(KIND=JPIM) :: mdb_asr_pcloudy ! 'asr_pcloudy@radiance' 
INTEGER(KIND=JPIM) :: mdb_asr_pcloudy_low ! 'asr_pcloudy_low@radiance' 
INTEGER(KIND=JPIM) :: mdb_asr_pcloudy_middle ! 'asr_pcloudy_middle@radiance' 
INTEGER(KIND=JPIM) :: mdb_asr_pcloudy_high ! 'asr_pcloudy_high@radiance' 
!    vector of MDB_OBSCORDIAG_AT_ERRSTAT [15/3/11 AF]
INTEGER(KIND=JPIM) :: mdb_obscordiag_at_errstat(jp_numdiag)  ! 'obs_corr_diag@errstat' 
!*** latlon_rad@desc added to flag whether (lat,lon) are in radians or not [23/1/04 SS]
INTEGER(KIND=JPIM) ::  mdb_latlon_rad_at_desc ! 'latlon_rad@desc'
!*** SSM/I reorganization (Peter Bauer/Sami Saarinen) [11/3/04]
!    New entries (!n), Changed entries (!c), Removed entries (!%)
INTEGER(KIND=JPIM) :: mdb_prec_st_at_ssmi(2) ! 'prec_st[2]@ssmi' !n
INTEGER(KIND=JPIM) :: mdb_prec_cv_at_ssmi(2) ! 'prec_cv[2]@ssmi' !n
INTEGER(KIND=JPIM) :: mdb_rain_at_ssmi_body(2) ! 'rain[2]@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_snow_at_ssmi_body(2) ! 'snow[2]@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb1d_cost ! 'cost@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_sfc_rain_3d_fg ! 'sfc_rain_3d_fg@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_sfc_snow_3d_fg ! 'sfc_snow_3d_fg@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_sfc_rain_3d_an ! 'sfc_rain_3d_an@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_sfc_snow_3d_an ! 'sfc_snow_3d_an@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_rwp(4) ! 'rwp@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_swp(4) ! 'swp@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_cwp(4) ! 'cwp@ssmi' !n
INTEGER(KIND=JPIM) :: mdb1d_iwp(4) ! 'iwp@ssmi' !n
INTEGER(KIND=JPIM) :: mdb_rad_obs_at_ssmi_body ! 'rad_obs@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_fg_3d_at_ssmi_body ! 'rad_fg_3d@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_4dan_at_ssmi_body ! 'rad_4dan@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb1bpws ! 'tcwv_obs@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb1drep ! 'tcwv_fg@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_fg_depar_at_ssmi_body ! 'rad_fg_depar@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_an_depar_at_ssmi_body ! 'rad_an_depar@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_obs_err_at_ssmi_body ! 'rad_obs_err@ssmi_body' !n
INTEGER(KIND=JPIM) :: mdb_rad_bias_at_ssmi_body ! 'rad_bias@ssmi_body' !n
!*** checksum@hdr-word added to get over reproducibility problems for good [18/3/04 SS]
!    contains the SUM(obsvalue) WHERE press IS NOT NULL and obsvalue IS NOT NULL
!    just after bufr2odb (in mergeodb; routine odb/cma2odb/update_obsdb.F90)
INTEGER(KIND=JPIM) :: mdb_checksum_at_hdr ! 'checksum@hdr'
!*** New item for SMOS proposed by J. Muñoz Sabater [06/10/2009-06/2010]
INTEGER(KIND=JPIM) :: mdb_tb_flag_smos      ! 'report_tbflag@smos'
INTEGER(KIND=JPIM) :: mdb_tb_ang_smos       ! 'incidence_angle@smos'
INTEGER(KIND=JPIM) :: mdb_tb_far_smos       ! 'faradey_rot_angle@smos'
INTEGER(KIND=JPIM) :: mdb_tb_geo_smos       ! 'pixel_rot_angle@smos'
INTEGER(KIND=JPIM) :: mdb_polarisation_at_smos  ! 'polarisation@smos'
INTEGER(KIND=JPIM) :: mdb_tb_smos           ! 'tbvalue@smos'
INTEGER(KIND=JPIM) :: mdb_snapshot_id_smos  ! 'snapshot_id@smos'
INTEGER(KIND=JPIM) :: mdb_grid_pt_id_smos   ! 'grid_point_id@smos'
INTEGER(KIND=JPIM) :: mdb_ecount_smos       ! 'electron_count@smos'
INTEGER(KIND=JPIM) :: mdb_sun_bt_smos       ! 'sun_bt@smos'
INTEGER(KIND=JPIM) :: mdb_snapshot_acc_smos ! 'snapshot_acc@smos'
INTEGER(KIND=JPIM) :: mdb_rad_acc_pure_smos ! 'rad_acc_pure@smos'
INTEGER(KIND=JPIM) :: mdb_rad_acc_cros_smos ! 'rad_acc_cross@smos'
INTEGER(KIND=JPIM) :: mdb_footpr_ax1_smos   ! 'footprint_axis_1@smos'
INTEGER(KIND=JPIM) :: mdb_footpr_ax2_smos   ! 'footprint_axis_2@smos'
INTEGER(KIND=JPIM) :: mdb_water_frac_smos   ! 'water_fraction@smos'
INTEGER(KIND=JPIM) :: mdb_info_smos         ! 'info@smos'
INTEGER(KIND=JPIM) :: mdb_snapshot_qlt_smos ! 'snapshot_quality@smos'
!*** New item for radio occultation SBH [27/10/2004]
INTEGER(KIND=JPIM) :: mdb_radcurv            ! 'radcurv@gnssro'
INTEGER(KIND=JPIM) :: mdb_undulation        ! 'undulation@gnssro'
INTEGER(KIND=JPIM) :: mdb_anaprop_at_radar_body
INTEGER(KIND=JPIM) :: mdb_antenht_at_radar_station
INTEGER(KIND=JPIM) :: mdb_beamwidth_at_radar_station
INTEGER(KIND=JPIM) :: mdb_distance_at_radar_body
INTEGER(KIND=JPIM) :: mdb_elevation_at_radar_body
INTEGER(KIND=JPIM) :: mdb_failure_1dv_at_radar
INTEGER(KIND=JPIM) :: mdb_flgdyn_at_radar_body
INTEGER(KIND=JPIM) :: mdb_frequency_at_radar_station
INTEGER(KIND=JPIM) :: mdb_ident_at_radar_station
INTEGER(KIND=JPIM) :: mdb_iternoconv_1dv_at_radar
INTEGER(KIND=JPIM) :: mdb_lat_at_radar_station
INTEGER(KIND=JPIM) :: mdb_lon_at_radar_station
INTEGER(KIND=JPIM) :: mdb_polarisation_at_radar_body
INTEGER(KIND=JPIM) :: mdb_azimuth_at_radar_body
INTEGER(KIND=JPIM) :: mdb_press_at_radar_body
INTEGER(KIND=JPIM) :: mdb_q1_at_radar_body
INTEGER(KIND=JPIM) :: mdb_q2_at_radar_body
INTEGER(KIND=JPIM) :: mdb_q_1dv_at_radar_body
INTEGER(KIND=JPIM) :: mdb_reflcost_at_radar_body
INTEGER(KIND=JPIM) :: mdb_stalt_at_radar_station
INTEGER(KIND=JPIM) :: mdb_temp1_at_radar_body
INTEGER(KIND=JPIM) :: mdb_temp2_at_radar_body
INTEGER(KIND=JPIM) :: mdb_temp_1dv_at_radar_body
INTEGER(KIND=JPIM) :: mdb_time_at_radar_body
INTEGER(KIND=JPIM) :: mdb_type_at_radar_station
INTEGER(KIND=JPIM) :: mlnk_sat2radar_station(2)
INTEGER(KIND=JPIM) :: mlnk_sat2radar(2)
INTEGER(KIND=JPIM) :: mlnk_radar2radar_body(2)
!*AF entries for gnssro (MF) to be moved to gnssro tables...
INTEGER(KIND=JPIM) :: mdb_obs_tvalue_at_robody          ! 'obs_tvalue@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_obs_zvalue_at_robody          ! 'obs_zvalue@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_bg_tvalue_at_robody           ! 'bg_tvalue@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_obs_dndz_at_robody            ! 'obs_dndz@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_obs_refractivity_at_robody    ! 'obs_refractivity@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_bg_dndz_at_robody             ! 'bg_dndz@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_bg_refractivity_at_robody     ! 'bg_refractivity@gnssro_body'
INTEGER(KIND=JPIM) :: mdb_bg_layerno_at_robody          ! 'bg_layerno@gnssro_body'
!*** New item for radar
INTEGER(KIND=JPIM) :: mdb_qmod_at_radar(jpmx_radar_niv)     ! 'qmod@radar'
INTEGER(KIND=JPIM) :: mdb_zsimp_at_radar(jpmx_radar_niv)    ! 'zsimp@radar'
!*** New item for radio occultation SBH [27/10/2004]
INTEGER(KIND=JPIM) :: mdb_ntan_at_limb      ! 'ntan@limb'
INTEGER(KIND=JPIM) :: mdb_ztan_at_limb(jpmx_limb_tan) ! 'ztan@limb'
INTEGER(KIND=JPIM) :: mdb_ptan_at_limb(jpmx_limb_tan) ! 'ptan@limb'
INTEGER(KIND=JPIM) :: mdb_thtan_at_limb(jpmx_limb_tan) ! 'thtan@limb'
INTEGER(KIND=JPIM) :: mdb_cloud_index_at_limb(jpmx_limb_tan) ! 'cloud_index@limb'
INTEGER(KIND=JPIM) :: mdb_window_rad_at_limb(jpmx_limb_tan) ! 'window_rad@limb'

!*** New @atovs entry proposed by Thomas Auligne [03/02/2005]
!*** moved to radiance table
INTEGER(KIND=JPIM) :: mdb_cldcover_at_radiance ! 'cldcover@radiance'

!*** New @atovs entries proposed by Andrew Collard [29/09/2008]

INTEGER(KIND=JPIM) :: mdb_avhrr_num_clusters  ! 'avhrr_num_clusters@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_max_cluster   ! 'avhrr_max_cluster@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_mean_ir      ! 'avhrr_mean_ir@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_stddev_ir    ! 'avhrr_stddev_ir@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_mean_vis      ! 'avhrr_mean_vis@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_stddev_vis    ! 'avhrr_stddev_vis@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_cold_cluster  ! 'avhrr_coldest_cluster_ir@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_warm_cluster  ! 'avhrr_warmest_cluster_ir@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_larg_cluster  ! 'avhrr_largest_cluster_ir@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_provider_qc          ! 'provider_qc@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_stddev_ir2    ! 'avhrr_stddev_ir2@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl1      ! 'avhrr_frac_cl1@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl2      ! 'avhrr_frac_cl2@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl3      ! 'avhrr_frac_cl3@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl4      ! 'avhrr_frac_cl4@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl5      ! 'avhrr_frac_cl5@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl6      ! 'avhrr_frac_cl6@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_frac_cl7      ! 'avhrr_frac_cl7@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl1     ! 'avhrr_m_ir1_cl1@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl2     ! 'avhrr_m_ir1_cl2@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl3     ! 'avhrr_m_ir1_cl3@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl4     ! 'avhrr_m_ir1_cl4@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl5     ! 'avhrr_m_ir1_cl5@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl6     ! 'avhrr_m_ir1_cl6@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir1_cl7     ! 'avhrr_m_ir1_cl7@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl1     ! 'avhrr_m_ir2_cl1@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl2     ! 'avhrr_m_ir2_cl2@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl3     ! 'avhrr_m_ir2_cl3@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl4     ! 'avhrr_m_ir2_cl4@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl5     ! 'avhrr_m_ir2_cl5@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl6     ! 'avhrr_m_ir2_cl6@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_m_ir2_cl7     ! 'avhrr_m_ir2_cl7@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_fg_ir1        ! 'avhrr_fg_ir1@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_fg_ir2        ! 'avhrr_fg_ir2@collocated_imager_information'
INTEGER(KIND=JPIM) :: mdb_avhrr_cloud_flag    ! 'avhrr_cloud_flag@collocated_imager_information'


!*** New items for cloud  sink variable from Tony Mc
INTEGER(KIND=JPIM) :: mdb_ctopbg_at_cloud_sink ! 'ctopbg@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_ctoper_at_cloud_sink ! 'ctoper@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_ctopinc_at_cloud_sink ! 'ctopinc@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_ctop_at_cloud_sink(jpmxup) ! 'ctop@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_camtbg_at_cloud_sink ! 'camtbg@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_camter_at_cloud_sink ! 'camter@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_camtinc_at_cloud_sink ! 'camtinc@cloud_sink'
INTEGER(KIND=JPIM) :: mdb_camt_at_cloud_sink(jpmxup) ! 'camt@cloud_sink'

INTEGER(KIND=JPIM) :: mdb_nensemble_at_ensemble   ! 'nensemble@ensemble'
!*** Ensemble Data Assimilation (enda) by Sami Saarinen [26/08/2005]
INTEGER(KIND=JPIM) :: mdb_enda_member_at_desc ! 'enda_member@desc'

!*** Dialing country code for the future speedup of SELECT ... WHERE inside(lat,lon,'country_name') [28/09/2005/SS]

!*** New @body entries proposed by Richard Engelen [07/07/2006]
INTEGER(KIND=JPIM) :: mdb_tausfc ! 'tausfc@radiance_body'
!*** New items added for new table resat_ak by Antje Dethof [20/02/2006]
INTEGER(KIND=JPIM) :: mdb_nak_at_resat_ak ! 'nak@resat_ak'
INTEGER(KIND=JPIM) :: mdb_pak_at_resat_ak(jpmx_ak) ! 'pak@resat_avergaing_kernel'
INTEGER(KIND=JPIM) :: mdb_wak_at_resat_ak(jpmx_ak) ! 'wak@resat_averaging_kernel'
INTEGER(KIND=JPIM) :: mlnk_resat2resat_ak(2)

!*** Yannick : Note: hardcoded #3 ;-(
INTEGER(KIND=JPIM) ::  mdb_obs_diags(3,jpmxup) ! f.ex. obs_diags[7]@update[5] i.e. obs_diags@update_5 :== ROBODY(JBODY, MDB_OBS_DIAGS(7,5))

!*** Changes due to BoM/MetOffice -compatibility (by SS, from their odb_CY32R2.039.tgz around 21-Nov-2007)
INTEGER(KIND=JPIM) :: mdb_dd_best_at_satob ! 'dd_best@satob'
INTEGER(KIND=JPIM) :: mdb_ff_best_at_satob ! 'ff_best@satob'

!*** New entry to hdr table for extended GEMS resat bufr by A. Dethof [05/03/2007]
INTEGER(KIND=JPIM) :: mdb_retrsource_at_resat ! retrsource@resat
!*** New entry to errstat table for extended GEMS resat bufr by A. Dethof [05/03/2007]
INTEGER(KIND=JPIM) :: mdb_obs_ak_error_at_errstat ! obs_ak_error@errstat
!*** New entry to resat table for Sciamachy CH4 correction factor by R. Engelen [11/02/2008]
INTEGER(KIND=JPIM) :: mdb_ch4cor_at_resat ! 'methane_correction@resat'
!*** New @atovs_body entry proposed by Fatima Karbou [16/04/2007]
INTEGER(KIND=JPIM) :: mdb_emis_atlas  ! emis_atlas@radiance_body
INTEGER(KIND=JPIM) :: mdb_emis_atlas_error  ! emis_atlas_error@radiance_body
INTEGER(KIND=JPIM) :: mdb_emis_retr   ! emis_retr@radiance_body
INTEGER(KIND=JPIM) :: mdb_emis_fg     ! emis_fg@radiance_body
INTEGER(KIND=JPIM) :: mdb_emis_an     ! emis_an@radiance_body
INTEGER(KIND=JPIM) :: mdb_emis_rtin   ! emis_rtin@radiance_body
INTEGER(KIND=JPIM) :: mdb_specularity ! specularity@radiance_body
INTEGER(KIND=JPIM) :: mdb_skintemp_retr  ! skintemp_retr@radiance_body
INTEGER(KIND=JPIM) :: mdb_tsfc        ! tsfc@modsurf

!*** Add prior profile to resat averaging kernels by Richard Engelen [22/10/2008]
INTEGER(KIND=JPIM) :: mdb_apak_at_resat_ak(jpmx_ak) ! 'apak@resat_averaging_kernel'
INTEGER(KIND=JPIM) :: mdb_sfc_height_at_resat ! 'surface_height@resat'
!*** New entry to body table for layer number by R. Dragani [21/01/2009]
INTEGER(KIND=JPIM) :: mdb_nlayer_at_body ! nlayer@body
!*** AF [31/10/2008]  New entry in hdr to keep a copy of seqno (subwindows --> Yannick)
INTEGER(KIND=JPIM) :: mdb_subseqno_at_hdr      ! 'subseqno@hdr'

INTEGER(KIND=JPIM) :: mdb_zenith_at_sat          ! zenith@sat
INTEGER(KIND=JPIM) :: mdb_azimuth_at_sat         ! azimuth@sat
INTEGER(KIND=JPIM) :: mdb_solar_zenith_at_sat    ! solar_zenith@sat
INTEGER(KIND=JPIM) :: mdb_solar_azimuth_at_sat   ! solar_azimuth@sat
INTEGER(KIND=JPIM) :: mdb_zenith_by_channel      ! zenith_by_channel@radiance_body

INTEGER(KIND=JPIM) :: mdb_lsm_fov_at_sat         ! lsm_fov@sat

INTEGER(KIND=JPIM) :: mdb_seaice_estimate
INTEGER(KIND=JPIM) :: mdb_fg_rain_rate          ! fg_rain_rate@allsky
INTEGER(KIND=JPIM) :: mdb_fg_snow_rate          ! fg_snow_rate@allsky
INTEGER(KIND=JPIM) :: mdb_fg_tcwv               ! fg_tcwv@allsky
INTEGER(KIND=JPIM) :: mdb_fg_cwp                ! fg_cwp@allsky
INTEGER(KIND=JPIM) :: mdb_fg_iwp                ! fg_iwp@allsky
INTEGER(KIND=JPIM) :: mdb_fg_rwp                ! fg_rwp@allsky
INTEGER(KIND=JPIM) :: mdb_fg_swp                ! fg_swp@allsky
INTEGER(KIND=JPIM) :: mdb_fg_rttov_cld_frac     ! fg_rttov_cld_fraction@allsky
INTEGER(KIND=JPIM) :: mdb_fg_theta700           ! fg_theta700@allsky
INTEGER(KIND=JPIM) :: mdb_fg_thetasfc           ! fg_thetasfc@allsky
INTEGER(KIND=JPIM) :: mdb_fg_uth                ! fg_uth@allsky
INTEGER(KIND=JPIM) :: mdb_fg_conv               ! fg_conv@allsky
INTEGER(KIND=JPIM) :: mdb_fg_pbl                ! fg_pbl@allsky

INTEGER(KIND=JPIM) :: mdb_an_rain_rate          ! an_rain_rate@allsky
INTEGER(KIND=JPIM) :: mdb_an_snow_rate          ! an_snow_rate@allsky
INTEGER(KIND=JPIM) :: mdb_an_tcwv               ! an_tcwv@allsky
INTEGER(KIND=JPIM) :: mdb_an_cwp                ! an_cwp@allsky
INTEGER(KIND=JPIM) :: mdb_an_iwp                ! an_iwp@allsky
INTEGER(KIND=JPIM) :: mdb_an_rwp                ! an_rwp@allsky
INTEGER(KIND=JPIM) :: mdb_an_swp                ! an_swp@allsky
INTEGER(KIND=JPIM) :: mdb_an_rttov_cld_frac     ! an_rttov_cld_fraction@allsky
INTEGER(KIND=JPIM) :: mdb_an_theta700           ! an_theta700@allsky
INTEGER(KIND=JPIM) :: mdb_an_thetasfc           ! an_thetasfc@allsky
INTEGER(KIND=JPIM) :: mdb_an_uth                ! an_uth@allsky
INTEGER(KIND=JPIM) :: mdb_an_conv               ! an_conv@allsky
INTEGER(KIND=JPIM) :: mdb_an_pbl                ! an_pbl@allsky

INTEGER(KIND=JPIM) :: mdb_gnorm_10mwind         ! gnorm_10mwind@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_skintemp        ! gnorm_skintemp@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_temp            ! gnorm_temp@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_q               ! gnorm_q@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_rainflux        ! gnorm_rainflux@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_snowflux        ! gnorm_snowflux@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_clw             ! gnorm_clw@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_ciw             ! gnorm_ciw@allsky
INTEGER(KIND=JPIM) :: mdb_gnorm_cc              ! gnorm_cc@allsky

INTEGER(KIND=JPIM) :: mdb_report_tbcloud        ! report_tbcloud@allsky
INTEGER(KIND=JPIM) :: mdb_tbvalue               ! tbvalue@allsky_body
INTEGER(KIND=JPIM) :: mdb_tbvaluead             ! tbvaluead@allsky_body
INTEGER(KIND=JPIM) :: mdb_tbvaluetl             ! tbvaluetl@allsky_body
INTEGER(KIND=JPIM) :: mdb_datum_tbflag          ! datum_tbflag@allsky_body
INTEGER(KIND=JPIM) :: mdb_cloudprox_type        ! cloudprox_type@allsky_body
INTEGER(KIND=JPIM) :: mdb_ob_cloudprox          ! ob_cloudprox@allsky_body
INTEGER(KIND=JPIM) :: mdb_fg_cloudprox          ! fg_cloudprox@allsky_body
INTEGER(KIND=JPIM) :: mdb_an_cloudprox          ! an_cloudprox@allsky_body

!*** New item for radar rain rate and rain gauges
INTEGER(KIND=JPIM) :: mdb_rrvalue        ! 'rrvalue@gbrad_body' or 'rrvalue@raingg_body'
INTEGER(KIND=JPIM) :: mdb_rrvaluetl      ! 'rrvaluetl@gbrad_body' or 'rrvaluetl@raingg_body'
INTEGER(KIND=JPIM) :: mdb_rrvaluead      ! 'rrvaluetl@gbrad_body' or 'rrvaluetl@raingg_body'
INTEGER(KIND=JPIM) :: mdb_report_rrflag  ! 'report_rrflag@gbrad' or 'report_rrflag@raingg'

!*** New item for rain gauges
INTEGER(KIND=JPIM) :: mdb_wdeff_bcorr    ! 'wdeff_bcorr@body'

!*** New items for cloud radar and lidar observations
INTEGER(KIND=JPIM) :: mdb_report_clreflflag   ! report_clreflflag@clrad
INTEGER(KIND=JPIM) :: mdb_clreflvalue         ! clreflvalue@clrad_body
INTEGER(KIND=JPIM) :: mdb_clreflvaluetl       ! clreflvaluetl@clrad_body
INTEGER(KIND=JPIM) :: mdb_clreflvaluead       ! clreflvaluead@clrad_body
INTEGER(KIND=JPIM) :: mdb_clreflvaluefg       ! clreflvaluefg@clrad_body
INTEGER(KIND=JPIM) :: mdb_datum_clrefl_stflag ! datum_clrefl_stflag@clrad_body
INTEGER(KIND=JPIM) :: mdb_report_clbscflag    ! report_clbscflag@claerlid
INTEGER(KIND=JPIM) :: mdb_clbscvalue          ! clbscvalue@claerlid_body
INTEGER(KIND=JPIM) :: mdb_clbscvaluetl        ! clbscvaluetl@claerlid_body
INTEGER(KIND=JPIM) :: mdb_clbscvaluead        ! clbscvaluead@claerlid_body
INTEGER(KIND=JPIM) :: mdb_clbscvaluefg        ! clbscvaluefg@claerlid_body
INTEGER(KIND=JPIM) :: mdb_datum_clbsc_stflag  ! datum_clbsc_stflag@claerlid_body

INTEGER(KIND=JPIM) :: mdb_fg_presf            ! mdb_fg_presf@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_z                ! mdb_fg_z@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_temp             ! mdb_fg_temp@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_spech            ! mdb_fg_spech@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_cwc              ! mdb_fg_cwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_iwc              ! mdb_fg_iwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_cc               ! mdb_fg_cc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_rwc              ! mdb_fg_rwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_swc              ! mdb_fg_swc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_pfra             ! mdb_fg_pfra@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_crwc             ! mdb_fg_crwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_fg_cswc             ! mdb_fg_cswc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_presf            ! mdb_an_presf@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_temp             ! mdb_an_temp@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_spech            ! mdb_an_spech@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_cwc              ! mdb_an_cwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_iwc              ! mdb_an_iwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_cc               ! mdb_an_cc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_rwc              ! mdb_an_rwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_swc              ! mdb_an_swc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_pfra             ! mdb_an_pfra@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_crwc             ! mdb_an_crwc@clrad_body/claerlid_body
INTEGER(KIND=JPIM) :: mdb_an_cswc             ! mdb_an_cswc@clrad_body/claerlid_body

! Auxiliary table entries
INTEGER(KIND=JPIM) :: mdb_aux_at_auxiliary(jp_numaux) ! report_aux@auxiliary
INTEGER(KIND=JPIM) :: mdb_aux_at_auxiliary_body(jp_numaux) ! datum_aux@auxiliary_body

! new radiance table
INTEGER(KIND=JPIM) :: mdb_scanline_at_radiance           ! scanline@radiance
INTEGER(KIND=JPIM) :: mdb_scanpos_at_radiance            ! scanpos@radiance
INTEGER(KIND=JPIM) :: mdb_orbit_at_radiance              ! orbit@radiance
INTEGER(KIND=JPIM) :: mdb_typesurf_at_radiance           ! typesurf@radiance
INTEGER(KIND=JPIM) :: mdb_corr_version_at_radiance      ! corr_version@radiance
INTEGER(KIND=JPIM) :: mdb_tbcorr_at_body                ! tbcorr@body
INTEGER(KIND=JPIM) :: mdb_skintemper_at_radiance         ! skintemper@radiance
INTEGER(KIND=JPIM) :: mdb_skintemp_at_radiance(jpmxup+1)  ! skintemp@radiance
INTEGER(KIND=JPIM) :: mdb_sink1er_at_radiance
INTEGER(KIND=JPIM) :: mdb_sink1_at_radiance(jpmxup+1)
INTEGER(KIND=JPIM) :: mdb_sink2er_at_radiance
INTEGER(KIND=JPIM) :: mdb_sink2_at_radiance(jpmxup+1)
INTEGER(KIND=JPIM) :: mdb_sink3er_at_radiance
INTEGER(KIND=JPIM) :: mdb_sink3_at_radiance(jpmxup+1)
INTEGER(KIND=JPIM) :: mdb_sink4er_at_radiance
INTEGER(KIND=JPIM) :: mdb_sink4_at_radiance(jpmxup+1)
INTEGER(KIND=JPIM) :: mdb_scatterindex_89_157           ! scatterindex_89_157@radiance
INTEGER(KIND=JPIM) :: mdb_scatterindex_23_89            ! scatterindex_23_89@radiance
INTEGER(KIND=JPIM) :: mdb_scatterindex_23_165           ! scatterindex_23_165@radiance
INTEGER(KIND=JPIM) :: mdb_lwp_obs                       ! lwp_obs@radiance

!*** New @body entries proposed by Tony McNally [29/11/2002]  TO BE CLEANED for AIRS AND IASI ONLY!!!
INTEGER(KIND=JPIM) :: mdb_rank_cld ! 'rank_cld@radiance_body' 
INTEGER(KIND=JPIM) :: mdb_jacobian_peak   ! 'jacobian_peak@body'
INTEGER(KIND=JPIM) :: mdb_jacobian_peakl  ! 'jacobian_peakl@body'
INTEGER(KIND=JPIM) :: mdb_jacobian_hpeak  ! 'jacobian_hpeak@body'
INTEGER(KIND=JPIM) :: mdb_jacobian_hpeakl ! 'jacobian_hpeakl@body'

!*** New @radiance_body entry to replace that @allsky_body (S. Migliorini and A. Geer) [27/08/2014]
INTEGER(KIND=JPIM) :: mdb_tbclear ! 'tbclear@radiance_body'

!*** New entries for forecast_diagnostic table - asked by Gabor Radnoti
INTEGER(KIND=JPIM) :: mdb_max_fcdiag_at_fcdiag                ! 'max_fcdiag@fcdiagnostic'
INTEGER(KIND=JPIM) :: mdb_fc_depar_at_fcdiag_body(jpmxfcdiag)  ! 'fc_depar@fcdiagnostic_body' 
INTEGER(KIND=JPIM) :: mdb_fc_step_at_fcdiag_body(jpmxfcdiag)   ! 'fc_step@fcdiagnostic_body' 

!*** New entries for conv and conv_body tables
INTEGER(KIND=JPIM) :: mdb_anemoht_at_conv                      ! 'anemoht@conv'
INTEGER(KIND=JPIM) :: mdb_baroht_at_conv                       ! 'baroht@conv'
INTEGER(KIND=JPIM) :: mdb_level_at_conv_body                   ! 'level@conv_body'
INTEGER(KIND=JPIM) :: mdb_ppcode_at_conv_body                  ! 'ppcode@conv_body'
INTEGER(KIND=JPIM) :: mdb_datum_qcflag_at_conv_body            ! 'datum_qcflag@conv_body'
INTEGER(KIND=JPIM) :: mdb_sonde_type_at_conv                   ! 'sonde_type@conv'
INTEGER(KIND=JPIM) :: mdb_flight_phase_at_conv                 ! 'flight_phase@conv'
INTEGER(KIND=JPIM) :: mdb_flight_dp_o_dt_at_conv               ! flight_dp_o_dt@conv
INTEGER(KIND=JPIM) :: mdb_station_type_at_conv                 ! 'station_type@conv'
INTEGER(KIND=JPIM) :: mdb_country_at_conv                      ! 'country@conv'
INTEGER(KIND=JPIM) :: mdb_aircraft_type_at_conv                ! 'aircraft_type@conv'
INTEGER(KIND=JPIM) :: mdb_heading_at_conv                      ! 'heading@conv'
INTEGER(KIND=JPIM) :: mlnk_hdr2conv(2)                            
INTEGER(KIND=JPIM) :: mlnk_conv2conv_body(2)

!*** New entry for reanalysis [21/06/2011] - asked by Hans Hersbach
INTEGER(KIND=JPIM) :: mdb_cid_at_conv      ! 'collection_identifier@conv'
INTEGER(KIND=JPIM) :: mdb_uid_at_conv      ! 'unique_identifier@conv'

!*** New entry for reanalysis [25/06/2012] - added by Paul Poli
INTEGER(KIND=JPIM) :: mdb_tsix_at_conv     ! 'timeseries_index@conv'
INTEGER(KIND=JPIM) :: mdb_biasvolatility   ! 'bias_volatility@body'

!*** New @radiance entries for aerosol [09/09/2016]- added by Julie Letertre-Danczak
INTEGER(KIND=JPIM) :: mdb_dust_aod_ir ! 'dust_aod_ir@radiance_body' 

!*** New entry for radiosonde drift by Tomas Kral [27/02/2017]
INTEGER(KIND=JPIM) :: mdb_reportno_at_hdr  ! 'reportno@hdr'
