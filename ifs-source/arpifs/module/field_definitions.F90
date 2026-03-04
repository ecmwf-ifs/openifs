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

! -- Field definitions - prototype AJG 24/10/2012 --
!
! This module contains all the field-specific information, e.g.
! names, dimensions, ID number, attributes, field classes. 
! There are a number of tedious but undemanding steps 
! required to add a new field:
!
!   (F1) Create a pointer in the field_access type
!   (F2) Add a field ID name...
!   (F3) ... and number (in sequence)
!   (F4) Put the field in any appropriate "classes" (field_class.F90)
!   (F5) Fill in other details (e.g. name as a string, dimensions, GRIB code?)
!   (F6) Set up the mapping from a field ID to the named pointer in field_access

! If a new model field is required, this should be the only place 
! that needs to be changed. Hopefully the tedious stuff can be eliminated
! by auto-generating this file in the future (e.g. from a central database of 
! model fields in xml)

! Problem 1: what about multi-dimension GOM fields, i.e. GHG, GRG, AERO, PHYS?
! Ideally these should be replaced by individual fields, which can be grouped in a "class"

! Problem 2: what about GOM-only "fake" fields like ul, vl, tl, ql? Ideally, GOM should do
! vertical interpolation or sub-sampling to a smaller number of model levels, thus removing the
! need for ul, vl, etc. (but it should be possible for GOM to create and use fields for itself
! if necessary).

! Problem 3: 2D and 3D canari model error fields used by GOM

module field_definitions

use field_definitions_base, only: set_fvar, type_fvar, field_access_base, field_metadata_base, &
  & jp_name_max_len, jp_necv_2d_max, jp_necv_3d_max, jp_name_max_len, jp_comments_max_len

use parkind1, only: jpim, jprb
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMCST   , ONLY : RV, RCPV, RCW, RCS
USE YOMLUN   , ONLY : NULOUT
USE YOM_GRIB_CODES

implicit none
public

! This type is to allow "named" and "shaped" user access to fields. (F1)
type, extends(field_access_base) :: field_access

  ! GMV
  real(kind=jprb), pointer :: u(:,:) => null()
  real(kind=jprb), pointer :: v(:,:) => null()
  real(kind=jprb), pointer :: div(:,:) => null()
  real(kind=jprb), pointer :: vor(:,:) => null()
  real(kind=jprb), pointer :: t(:,:) => null()
  real(kind=jprb), pointer :: pd(:,:) => null()
  real(kind=jprb), pointer :: vd(:,:) => null()   
  real(kind=jprb), pointer :: nhx(:,:) => null()
  real(kind=jprb), pointer :: edot(:,:) => null()
  real(kind=jprb), pointer :: vcasrs(:,:) => null() ! ky: obsolete field (deep-layer models pruned)
  real(kind=jprb), pointer :: rdlasr(:,:) => null() ! ky: obsolete field (deep-layer models pruned)
  real(kind=jprb), pointer :: sgrtl(:,:) => null()
  real(kind=jprb), pointer :: sgrtm(:,:) => null()
  real(kind=jprb), pointer :: unl(:,:) => null()
  real(kind=jprb), pointer :: unl_si(:,:) => null()
  real(kind=jprb), pointer :: vnl(:,:) => null()
  real(kind=jprb), pointer :: vnl_si(:,:) => null()
  real(kind=jprb), pointer :: tnl(:,:) => null()
  real(kind=jprb), pointer :: tnl_si(:,:) => null()
  real(kind=jprb), pointer :: spnl(:,:)  => null()
  real(kind=jprb), pointer :: spnl_si(:,:)  => null()
  real(kind=jprb), pointer :: vwvnl(:,:) => null()
  real(kind=jprb), pointer :: vdnl_si(:,:) => null()
  real(kind=jprb), pointer :: pdnl(:,:) => null()
  real(kind=jprb), pointer :: pdnl_si(:,:) => null()
  real(kind=jprb), pointer :: gw(:,:) => null()
  real(kind=jprb), pointer :: curhs(:,:) => null()
  real(kind=jprb), pointer :: cvrhs(:,:) => null()
  real(kind=jprb), pointer :: ctrhs(:,:) => null()
  real(kind=jprb), pointer :: csprhs(:) => null()
  real(kind=jprb), pointer :: cspdrhs(:,:) => null()
  real(kind=jprb), pointer :: csvdrhs(:,:) => null()
  real(kind=jprb), pointer :: cunl(:,:) => null()
  real(kind=jprb), pointer :: cvnl(:,:) => null()
  real(kind=jprb), pointer :: ctnl(:,:) => null()
  real(kind=jprb), pointer :: cspnl(:,:) => null()
  real(kind=jprb), pointer :: cvwvnl(:,:) => null()
  real(kind=jprb), pointer :: cspdnl(:,:) => null()
  real(kind=jprb), pointer :: cpdnl(:,:) => null()
  real(kind=jprb), pointer :: cupt(:,:) => null()
  real(kind=jprb), pointer :: cvpt(:,:) => null()
  real(kind=jprb), pointer :: ctpt(:,:) => null()
  real(kind=jprb), pointer :: cpdpt(:,:) => null()
  real(kind=jprb), pointer :: cvdpt(:,:) => null()
  real(kind=jprb), pointer :: dphi(:,:) => null()
  real(kind=jprb), pointer :: nhxnl(:,:) => null()
  real(kind=jprb), pointer :: cnhxnl(:,:) => null()
  ! GFL
  real(kind=jprb), pointer :: q(:,:) => null()
  real(kind=jprb), pointer :: o3(:,:) => null()
  real(kind=jprb), pointer :: l(:,:) => null()
  real(kind=jprb), pointer :: i(:,:) => null()
  real(kind=jprb), pointer :: a(:,:) => null()
  real(kind=jprb), pointer :: s(:,:) => null()
  real(kind=jprb), pointer :: r(:,:) => null()
  real(kind=jprb), pointer :: g(:,:) => null()
  real(kind=jprb), pointer :: h(:,:) => null()
  real(kind=jprb), pointer :: lconv(:,:) => null()
  real(kind=jprb), pointer :: iconv(:,:) => null()
  real(kind=jprb), pointer :: rconv(:,:) => null()
  real(kind=jprb), pointer :: sconv(:,:) => null()
  real(kind=jprb), pointer :: tke(:,:) => null()
  real(kind=jprb), pointer :: tte(:,:) => null()
  real(kind=jprb), pointer :: efb1(:,:) => null()
  real(kind=jprb), pointer :: efb2(:,:) => null()
  real(kind=jprb), pointer :: efb3(:,:) => null()
  real(kind=jprb), pointer :: mxl(:,:) => null()
  real(kind=jprb), pointer :: src(:,:) => null()
  real(kind=jprb), pointer :: shtur(:,:) => null()
  real(kind=jprb), pointer :: fqtur(:,:) => null()
  real(kind=jprb), pointer :: fstur(:,:) => null()
  real(kind=jprb), pointer :: cvv(:,:) => null()
  real(kind=jprb), pointer :: rkth(:,:) => null()
  real(kind=jprb), pointer :: rktqv(:,:) => null()
  real(kind=jprb), pointer :: rktqc(:,:) => null()
  real(kind=jprb), pointer :: cpf(:,:) => null()
  real(kind=jprb), pointer :: uom(:,:) => null()
  real(kind=jprb), pointer :: ual(:,:) => null()
  real(kind=jprb), pointer :: dom(:,:) => null()
  real(kind=jprb), pointer :: dal(:,:) => null()
  real(kind=jprb), pointer :: uen(:,:) => null()
  real(kind=jprb), pointer :: unebh(:,:) => null()
  real(kind=jprb), pointer :: spf(:,:) => null()
  real(kind=jprb), pointer :: cvgq(:,:) => null()
  real(kind=jprb), pointer :: lrad(:,:) => null()
  real(kind=jprb), pointer :: irad(:,:) => null()
  real(kind=jprb), pointer :: rspec(:,:) => null()
  real(kind=jprb), pointer :: lrch4(:,:) => null()
  real(kind=jprb), pointer :: extra3d(:,:,:) => null()
  real(kind=jprb), pointer :: ezdiag(:,:,:) => null()
  real(kind=jprb), pointer :: unogw(:,:) => null()
  real(kind=jprb), pointer :: vnogw(:,:) => null()
  
  ! pressure
  real(kind=jprb), pointer :: pre(:,:) => null()
  real(kind=jprb), pointer :: pref(:,:) => null()
  real(kind=jprb), pointer :: delp(:,:) => null()
  ! GMVS
  real(kind=jprb), pointer :: sp(:)  => null()
  real(kind=jprb), pointer :: csppt(:)  => null()
  real(kind=jprb), pointer :: dbbc(:)  => null()
  real(kind=jprb), pointer :: gws(:)  => null()
  real(kind=jprb), pointer :: spnl2(:)  => null()
  real(kind=jprb), pointer :: cspnl2(:)  => null()
  real(kind=jprb), pointer :: prehyds(:)  => null()
  real(kind=jprb), pointer :: extcv2d(:,:)  => null()
  real(kind=jprb), pointer :: extcv3d(:,:)  => null()

  ! For LELAM
  real(kind=jprb), pointer :: u_mean(:) => null()
  real(kind=jprb), pointer :: v_mean(:) => null()
  
  ! SLB2
  real(kind=jprb), pointer :: vvel(:,:) => null()
  ! Constants
  real(kind=jprb), pointer :: orog(:) => null()
  real(kind=jprb), pointer :: cori(:) => null()
  real(kind=jprb), pointer :: gnordl(:) => null()
  real(kind=jprb), pointer :: gnordm(:) => null()
  ! Clouds
  real(kind=jprb), pointer :: ccc(:) => null()
  real(kind=jprb), pointer :: lcc(:) => null()
  real(kind=jprb), pointer :: mcc(:) => null()
  real(kind=jprb), pointer :: hcc(:) => null()
  real(kind=jprb), pointer :: tcc(:) => null()
  ! Surface
  ! Group SB
  real(kind=jprb), pointer :: sb_t(:,:) => null()
  real(kind=jprb), pointer :: sb_q(:,:) => null()
  real(kind=jprb), pointer :: sb_tl(:,:) => null()
  ! Group SG
  real(kind=jprb), pointer :: sg_f(:,:) => null()
  real(kind=jprb), pointer :: sg_a(:,:) => null()
  real(kind=jprb), pointer :: sg_r(:,:) => null()
  real(kind=jprb), pointer :: sg_t(:,:) => null()
  real(kind=jprb), pointer :: sg_w(:,:) => null()
  ! Group SL
  real(kind=jprb), pointer :: sl_lict(:) => null()
  real(kind=jprb), pointer :: sl_lmlt(:) => null()
  real(kind=jprb), pointer :: sl_ltlt(:) => null()
  real(kind=jprb), pointer :: sl_lblt(:) => null()
  real(kind=jprb), pointer :: sl_lshf(:) => null()
  real(kind=jprb), pointer :: sl_licd(:) => null()
  real(kind=jprb), pointer :: sl_lmld(:) => null()
  ! Group RR
  real(kind=jprb), pointer :: rr_t(:) => null()
  real(kind=jprb), pointer :: rr_tir(:) => null()
  real(kind=jprb), pointer :: rr_tmw(:) => null()
  real(kind=jprb), pointer :: rr_w(:) => null()
  real(kind=jprb), pointer :: rr_fc(:) => null()
  real(kind=jprb), pointer :: rr_ic(:) => null()
  real(kind=jprb), pointer :: rr_fp1(:) => null()
  ! Group CL
  real(kind=jprb), pointer :: cl_tcls(:) => null() 
  real(kind=jprb), pointer :: cl_hucls(:) => null() 
  ! Group OM=OML: prognostic quantities for ocean mixed layer model (KPP/TKE)
  real(kind=jprb), pointer :: om_to(:,:) => null() 
  real(kind=jprb), pointer :: om_so(:,:) => null() 
  real(kind=jprb), pointer :: om_uo(:,:) => null() 
  real(kind=jprb), pointer :: om_vo(:,:) => null() 
  ! Group EP=EXTRP: extra 3-d prognostic fields
  real(kind=jprb), pointer :: ep_ep(:,:) => null()
  ! Group X2=XTRP2: extra 2-d prognostic fields
  real(kind=jprb), pointer :: x2_x2(:,:) => null() 
  ! Group CI=CANRI: 2-d prognostic fields for CANARI
  real(kind=jprb), pointer :: ci_ci(:,:) => null() 
  ! Group VF
  real(kind=jprb), pointer :: vf_z0f(:) => null() 
  real(kind=jprb), pointer :: vf_albf(:) => null() 
  real(kind=jprb), pointer :: vf_emisf(:) => null() 
  real(kind=jprb), pointer :: vf_getrl(:) => null() 
  real(kind=jprb), pointer :: vf_lsm(:) => null() 
  real(kind=jprb), pointer :: vf_veg(:) => null() 
  real(kind=jprb), pointer :: vf_vrlan(:) => null() 
  real(kind=jprb), pointer :: vf_vrldi(:) => null() 
  real(kind=jprb), pointer :: vf_sig(:) => null() 
  real(kind=jprb), pointer :: vf_albsf(:) => null() 
  real(kind=jprb), pointer :: vf_lan(:) => null() 
  real(kind=jprb), pointer :: vf_sst(:) => null() 
  real(kind=jprb), pointer :: vf_sss(:) => null() 
  real(kind=jprb), pointer :: vf_lz0h(:) => null() 
  real(kind=jprb), pointer :: vf_cvl(:) => null() 
  real(kind=jprb), pointer :: vf_co2typ(:) => null() 
  real(kind=jprb), pointer :: vf_cvh(:) => null()
  real(kind=jprb), pointer :: vf_fwet(:) => null()
  real(kind=jprb), pointer :: vf_cur(:) => null()
  real(kind=jprb), pointer :: vf_tvl(:) => null() 
  real(kind=jprb), pointer :: vf_tvh(:) => null() 
  real(kind=jprb), pointer :: vf_lail(:) => null() 
  real(kind=jprb), pointer :: vf_laih(:) => null() 
  real(kind=jprb), pointer :: vf_soty(:) => null() 
  real(kind=jprb), pointer :: vf_clk(:) => null() 
  real(kind=jprb), pointer :: vf_dl(:) => null() 
  real(kind=jprb), pointer :: vf_ci(:) => null() 
  real(kind=jprb), pointer :: vf_ucur(:) => null() 
  real(kind=jprb), pointer :: vf_vcur(:) => null()
  real(kind=jprb), pointer :: vf_z0rlf(:) => null()
  real(kind=jprb), pointer :: vf_cgpp(:) => null()
  real(kind=jprb), pointer :: vf_crec(:) => null()
  real(kind=jprb), pointer :: vf_sdfor(:) => null()
  real(kind=jprb), pointer :: vf_aluvp(:) => null() 
  real(kind=jprb), pointer :: vf_aluvd(:) => null() 
  real(kind=jprb), pointer :: vf_alnip(:) => null() 
  real(kind=jprb), pointer :: vf_alnid(:) => null() 
  real(kind=jprb), pointer :: vf_aluvi(:) => null() 
  real(kind=jprb), pointer :: vf_aluvv(:) => null() 
  real(kind=jprb), pointer :: vf_aluvg(:) => null() 
  real(kind=jprb), pointer :: vf_alnii(:) => null() 
  real(kind=jprb), pointer :: vf_alniv(:) => null() 
  real(kind=jprb), pointer :: vf_alnig(:) => null() 
  real(kind=jprb), pointer :: vf_fp1(:) => null()     ! surface orography in the 2nd part of FULLPOS-927
  real(kind=jprb), pointer :: vf_so2dd(:) => null()   ! sulphate dry dep velocity
  real(kind=jprb), pointer :: vf_dmso(:) => null()    ! oceanic DMS
  real(kind=jprb), pointer :: vf_urbf(:) => null()   ! SOA from CO
  real(kind=jprb), pointer :: vf_fca1(:) => null()    ! fraction of calcite over dust 1st bin
  real(kind=jprb), pointer :: vf_fca2(:) => null()    ! fraction of calcite over dust 2st bin
  real(kind=jprb), pointer :: vf_aerdep(:) => null()  ! dust emission potential 
  real(kind=jprb), pointer :: vf_aerlts(:) => null()  ! dust lifting threshold speed 
  real(kind=jprb), pointer :: vf_aerscc(:) => null()  ! dust soil clay content
  real(kind=jprb), pointer :: vf_dsf(:) => null()     ! dust soil clay content
  real(kind=jprb), pointer :: vf_dsz(:) => null()     ! Dust size dist modulate
  real(kind=jprb), pointer :: vf_chemflxo(:,:) => null() !  total chemistry flux (emissions + deposition)
  real(kind=jprb), pointer :: vf_chemwdflx(:,:) => null() !  wet deposition chemistry flux 
  real(kind=jprb), pointer :: vf_chemddflx(:,:) => null() !  dry deposition chemistry flux 
  real(kind=jprb), pointer :: vf_chemdv(:,:) => null()  ! chemistry deposition velocity
  real(kind=jprb), pointer :: vf_nudm(:) => null()    ! nudging mask
  real(kind=jprb), pointer :: vf_emis2d(:,:) => null()  ! 2D emission fields for composition
  real(kind=jprb), pointer :: vf_emis2daux(:,:) => null()  ! 2D emission auxiliary fields for composition
  ! Group VP=VCLIP: deep soil diagnostic fields
  real(kind=jprb), pointer :: vp_tpc(:) => null()
  real(kind=jprb), pointer :: vp_wpc(:) => null()
  ! Group VV=VCLIV: vegetation diagnostic fields
  real(kind=jprb), pointer :: vv_arg(:) => null() 
  real(kind=jprb), pointer :: vv_sab(:) => null() 
  real(kind=jprb), pointer :: vv_d2(:) => null() 
  real(kind=jprb), pointer :: vv_iveg(:) => null() 
  real(kind=jprb), pointer :: vv_rsmin(:) => null() 
  real(kind=jprb), pointer :: vv_lai(:) => null() 
  real(kind=jprb), pointer :: vv_hv(:) => null() 
  real(kind=jprb), pointer :: vv_z0h(:) => null() 
  real(kind=jprb), pointer :: vv_als(:) => null() 
  real(kind=jprb), pointer :: vv_alv(:) => null() 
  ! Group VN=VCLIN: cloudiness diagnostic predictors:
  real(kind=jprb), pointer :: vn_top(:) => null() 
  real(kind=jprb), pointer :: vn_bas(:) => null() 
  real(kind=jprb), pointer :: vn_acpr(:) => null() 
  real(kind=jprb), pointer :: vn_accpr(:) => null() 
  real(kind=jprb), pointer :: vn_accpr5(:) => null() 
  ! Group VH=VCLIH: convective cloud diagnostic fields
   real(kind=jprb), pointer :: vh_tcch(:) => null()
   real(kind=jprb), pointer :: vh_scch(:) => null()
   real(kind=jprb), pointer :: vh_bcch(:) => null()
   real(kind=jprb), pointer :: vh_pblh(:) => null()
   real(kind=jprb), pointer :: vh_spsh(:) => null()
   real(kind=jprb), pointer :: vh_qsh (:) => null()
  ! Group  VA=VCLIA: aerosol diagnostic fields
   real(kind=jprb), pointer :: va_sea(:) => null()
   real(kind=jprb), pointer :: va_lan(:) => null()
   real(kind=jprb), pointer :: va_soo(:) => null()
   real(kind=jprb), pointer :: va_des(:) => null()
   real(kind=jprb), pointer :: va_sul(:) => null()
   real(kind=jprb), pointer :: va_vol(:) => null()
  ! Group  VG=VCLIG: ice-coupler diagnostic fields
   real(kind=jprb), pointer :: vg_icfr(:) => null()
   real(kind=jprb), pointer :: vg_soup(:) => null()
   real(kind=jprb), pointer :: vg_irup(:) => null()
   real(kind=jprb), pointer :: vg_chss(:) => null()
   real(kind=jprb), pointer :: vg_evap(:) => null()
   real(kind=jprb), pointer :: vg_taux(:) => null()
   real(kind=jprb), pointer :: vg_tauy(:) => null()
  ! Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields
   real(kind=jprb), pointer :: vc_a(:) => null()
   real(kind=jprb), pointer :: vc_b(:) => null()
   real(kind=jprb), pointer :: vc_c(:) => null()
  ! Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP)
   real(kind=jprb), pointer :: v2_ocdep(:) => null()
   real(kind=jprb), pointer :: v2_ustrc(:) => null()
   real(kind=jprb), pointer :: v2_vstrc(:) => null()
  ! * Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
   real(kind=jprb), pointer :: v3_difm(:,:)  => null()
   real(kind=jprb), pointer :: v3_dift(:,:)  => null()
   real(kind=jprb), pointer :: v3_difs(:,:)  => null()
   real(kind=jprb), pointer :: v3_advt(:,:)  => null()
   real(kind=jprb), pointer :: v3_advs(:,:)  => null()
   real(kind=jprb), pointer :: v3_tri0(:,:)  => null()
   real(kind=jprb), pointer :: v3_tri1(:,:)  => null()
   real(kind=jprb), pointer :: v3_swdk(:,:)  => null()
   real(kind=jprb), pointer :: v3_zo(:,:)  => null()
   real(kind=jprb), pointer :: v3_ho(:,:)  => null()
   real(kind=jprb), pointer :: v3_do(:,:)  => null()
   real(kind=jprb), pointer :: v3_ho_inv(:,:)  => null()
   real(kind=jprb), pointer :: v3_uoc(:,:)  => null()
   real(kind=jprb), pointer :: v3_voc(:,:)  => null()
   real(kind=jprb), pointer :: v3_otke(:,:)  => null()
   ! Group VD
   real(kind=jprb), pointer :: vd_lsp(:) => null() 
   real(kind=jprb), pointer :: vd_cp(:) => null() 
   real(kind=jprb), pointer :: vd_sf(:) => null() 
   real(kind=jprb), pointer :: vd_fzra(:) => null() 
   real(kind=jprb), pointer :: vd_bld(:) => null() 
   real(kind=jprb), pointer :: vd_sshf(:) => null() 
   real(kind=jprb), pointer :: vd_slhf(:) => null() 
   real(kind=jprb), pointer :: vd_nee(:) => null() 
   real(kind=jprb), pointer :: vd_gpp(:) => null() 
   real(kind=jprb), pointer :: vd_rec(:) => null() 
   real(kind=jprb), pointer :: vd_msl(:) => null() 
   real(kind=jprb), pointer :: vd_sp(:) => null() 
   real(kind=jprb), pointer :: vd_tcc(:) => null() 
   real(kind=jprb), pointer :: vd_10u(:) => null() 
   real(kind=jprb), pointer :: vd_10v(:) => null() 
   real(kind=jprb), pointer :: vd_2t(:) => null() 
   real(kind=jprb), pointer :: vd_2d(:) => null() 
   real(kind=jprb), pointer :: vd_2sh(:) => null() 
   real(kind=jprb), pointer :: vd_ssr(:) => null() 
   real(kind=jprb), pointer :: vd_str(:) => null() 
   real(kind=jprb), pointer :: vd_tsr(:) => null() 
   real(kind=jprb), pointer :: vd_ttr(:) => null() 
   real(kind=jprb), pointer :: vd_ewss(:) => null() 
   real(kind=jprb), pointer :: vd_nsss(:) => null() 
   real(kind=jprb), pointer :: vd_e(:) => null() 
   real(kind=jprb), pointer :: vd_pev(:) => null() 
   real(kind=jprb), pointer :: vd_ccc(:) => null() 
   real(kind=jprb), pointer :: vd_lcc(:) => null() 
   real(kind=jprb), pointer :: vd_mcc(:) => null() 
   real(kind=jprb), pointer :: vd_hcc(:) => null() 
   real(kind=jprb), pointer :: vd_lgws(:) => null() 
   real(kind=jprb), pointer :: vd_mgws(:) => null() 
   real(kind=jprb), pointer :: vd_gwd(:) => null() 
   real(kind=jprb), pointer :: vd_mx2t(:) => null() 
   real(kind=jprb), pointer :: vd_mn2t(:) => null() 
   real(kind=jprb), pointer :: vd_mx2t6(:) => null() 
   real(kind=jprb), pointer :: vd_mn2t6(:) => null() 
   real(kind=jprb), pointer :: vd_ro(:) => null() 
   real(kind=jprb), pointer :: vd_sro(:) => null() 
   real(kind=jprb), pointer :: vd_ssro(:) => null() 
   real(kind=jprb), pointer :: vd_alb(:) => null() 
   real(kind=jprb), pointer :: vd_iewss(:) => null() 
   real(kind=jprb), pointer :: vd_insss(:) => null() 
   real(kind=jprb), pointer :: vd_isshf(:) => null() 
   real(kind=jprb), pointer :: vd_ie(:) => null() 
   real(kind=jprb), pointer :: vd_inee(:) => null() 
   real(kind=jprb), pointer :: vd_igpp(:) => null() 
   real(kind=jprb), pointer :: vd_irec(:) => null()
   real(kind=jprb), pointer :: vd_ich4(:) => null()
   real(kind=jprb), pointer :: vd_csf(:) => null() 
   real(kind=jprb), pointer :: vd_lssf(:) => null() 
   real(kind=jprb), pointer :: vd_mxtpr(:) => null() 
   real(kind=jprb), pointer :: vd_mntpr(:) => null() 
   real(kind=jprb), pointer :: vd_mxtpr6(:) => null() 
   real(kind=jprb), pointer :: vd_mntpr6(:) => null() 
   real(kind=jprb), pointer :: vd_tpr(:) => null() 
   real(kind=jprb), pointer :: vd_lsrr(:) => null() 
   real(kind=jprb), pointer :: vd_crr(:) => null() 
   real(kind=jprb), pointer :: vd_lssfr(:) => null() 
   real(kind=jprb), pointer :: vd_csfr(:) => null() 
   real(kind=jprb), pointer :: vd_ptype(:) => null() 
   real(kind=jprb), pointer :: vd_ilspf(:) => null() 
   real(kind=jprb), pointer :: vd_z0f(:) => null() 
   real(kind=jprb), pointer :: vd_lz0h(:) => null() 
   real(kind=jprb), pointer :: vd_tcw(:) => null() 
   real(kind=jprb), pointer :: vd_tcwv(:) => null() 
   real(kind=jprb), pointer :: vd_tclw(:) => null() 
   real(kind=jprb), pointer :: vd_tciw(:) => null() 
   real(kind=jprb), pointer :: vd_tcrw(:) => null() 
   real(kind=jprb), pointer :: vd_tcsw(:) => null() 
   real(kind=jprb), pointer :: vd_tcslw(:) => null() 
   real(kind=jprb), pointer :: vd_ssrd(:) => null() 
   real(kind=jprb), pointer :: vd_strd(:) => null() 
   real(kind=jprb), pointer :: vd_ssrdc(:) => null() 
   real(kind=jprb), pointer :: vd_strdc(:) => null() 
   real(kind=jprb), pointer :: vd_blh(:) => null() 
   real(kind=jprb), pointer :: vd_sund(:) => null() 
   real(kind=jprb), pointer :: vd_spar(:) => null() 
   real(kind=jprb), pointer :: vd_suvb(:) => null() 
   real(kind=jprb), pointer :: vd_sfdir(:) => null() 
   real(kind=jprb), pointer :: vd_scdir(:) => null() 
   real(kind=jprb), pointer :: vd_sdsrp(:) => null() 
   real(kind=jprb), pointer :: vd_cape(:) => null() 
   real(kind=jprb), pointer :: vd_capes(:) => null() 
   real(kind=jprb), pointer :: vd_mucape(:) => null() 
   real(kind=jprb), pointer :: vd_pdepl(:) => null() 
   real(kind=jprb), pointer :: vd_mlcape50(:) => null() 
   real(kind=jprb), pointer :: vd_mlcape100(:) => null() 
   real(kind=jprb), pointer :: vd_mlcin50(:) => null() 
   real(kind=jprb), pointer :: vd_mlcin100(:) => null() 
   real(kind=jprb), pointer :: vd_tropotp(:) => null() 
   real(kind=jprb), pointer :: vd_tsrc(:) => null() 
   real(kind=jprb), pointer :: vd_ttrc(:) => null() 
   real(kind=jprb), pointer :: vd_ssrc(:) => null() 
   real(kind=jprb), pointer :: vd_strc(:) => null() 
   real(kind=jprb), pointer :: vd_es(:) => null() 
   real(kind=jprb), pointer :: vd_smlt(:) => null() 
   real(kind=jprb), pointer :: vd_10fg(:) => null() 
   real(kind=jprb), pointer :: vd_10fg6(:) => null() 
   real(kind=jprb), pointer :: vd_10fgcv(:) => null() 
   real(kind=jprb), pointer :: vd_i10fg(:) => null() 
   real(kind=jprb), pointer :: vd_lspf(:) => null() 
   real(kind=jprb), pointer :: vd_tco3(:) => null() 
   real(kind=jprb), pointer :: vd_vimd(:) => null() 
   real(kind=jprb), pointer :: vd_sparc(:) => null() 
   real(kind=jprb), pointer :: vd_stinc(:) => null() 
   real(kind=jprb), pointer :: vd_cbase(:) => null() 
   real(kind=jprb), pointer :: vd_0degl(:) => null() 
   real(kind=jprb), pointer :: vd_m10degl(:) => null() 
   real(kind=jprb), pointer :: vd_visih(:) => null() 
   real(kind=jprb), pointer :: vd_cin(:) => null() 
   real(kind=jprb), pointer :: vd_kindex(:) => null() 
   real(kind=jprb), pointer :: vd_ttindex(:) => null() 
   real(kind=jprb), pointer :: vd_cbasea(:) => null() 
   real(kind=jprb), pointer :: vd_ctopc(:) => null() 
   real(kind=jprb), pointer :: vd_ztwetb0(:) => null() 
   real(kind=jprb), pointer :: vd_ztwetb1(:) => null() 
   real(kind=jprb), pointer :: vd_tcghg(:) => null() 
   real(kind=jprb), pointer :: vd_tcchem(:) => null() 
   real(kind=jprb), pointer :: vd_aerodiag(:,:) => null() 
   real(kind=jprb), pointer :: vd_aero_wvl_diag(:,:) => null() 
   real(kind=jprb), pointer :: vd_100u(:) => null() 
   real(kind=jprb), pointer :: vd_100v(:) => null() 
   real(kind=jprb), pointer :: vd_zust(:) => null() 
   real(kind=jprb), pointer :: vd_10nu(:) => null() 
   real(kind=jprb), pointer :: vd_10nv(:) => null() 
   real(kind=jprb), pointer :: vd_dndzn(:) => null() 
   real(kind=jprb), pointer :: vd_dndza(:) => null() 
   real(kind=jprb), pointer :: vd_dctb(:) => null() 
   real(kind=jprb), pointer :: vd_tplb(:) => null() 
   real(kind=jprb), pointer :: vd_tplt(:) => null() 
   real(kind=jprb), pointer :: vd_odss(:) => null() 
   real(kind=jprb), pointer :: vd_oddu(:) => null() 
   real(kind=jprb), pointer :: vd_odom(:) => null() 
   real(kind=jprb), pointer :: vd_odbc(:) => null() 
   real(kind=jprb), pointer :: vd_odsu(:) => null() 
   real(kind=jprb), pointer :: vd_odni(:) => null() 
   real(kind=jprb), pointer :: vd_odam(:) => null() 
   real(kind=jprb), pointer :: vd_odsoa(:) => null() 
   real(kind=jprb), pointer :: vd_odvfa(:) => null() 
   real(kind=jprb), pointer :: vd_odvsu(:) => null() 
   real(kind=jprb), pointer :: vd_odtoacc(:) => null() 
   real(kind=jprb), pointer :: vd_aepm1(:) => null() 
   real(kind=jprb), pointer :: vd_aepm25(:) => null() 
   real(kind=jprb), pointer :: vd_aepm10(:) => null() 
   real(kind=jprb), pointer :: vd_uvbed(:) => null() 
   real(kind=jprb), pointer :: vd_uvbedcs(:) => null() 
   real(kind=jprb), pointer :: vd_litoti(:) => null() 
   real(kind=jprb), pointer :: vd_licgi(:) => null() 
   real(kind=jprb), pointer :: vd_litota6(:) => null() 
   real(kind=jprb), pointer :: vd_licga6(:) => null()
   real(kind=jprb), pointer :: vd_ptypeocc6(:,:) => null()
   real(kind=jprb), pointer :: vd_200u(:) => null() 
   real(kind=jprb), pointer :: vd_200v(:) => null() 

   real(kind=jprb), pointer :: vd_sdsl(:) => null() 

! * Group SM=SATSIM: (ECMWF) simulated satellite images:
   real(kind=jprb), pointer :: sm_clbt(:,:) => null()
   real(kind=jprb), pointer :: sm_csbt(:,:) => null()
   
  ! Group WS
  real(kind=jprb), pointer :: ws_char(:) => null() 
  real(kind=jprb), pointer :: ws_charhq(:) => null() 
  real(kind=jprb), pointer :: ws_ustokes(:) => null() 
  real(kind=jprb), pointer :: ws_vstokes(:) => null() 
  real(kind=jprb), pointer :: ws_tauocx(:) => null() 
  real(kind=jprb), pointer :: ws_tauocy(:) => null() 
  real(kind=jprb), pointer :: ws_phioc(:) => null() 
  real(kind=jprb), pointer :: ws_wsemean(:) => null() 
  real(kind=jprb), pointer :: ws_wsfmean(:) => null()
  !  Group WW 
  real(kind=jprb), pointer :: ww_u10n(:) => null() 
  real(kind=jprb), pointer :: ww_v10n(:) => null() 
  real(kind=jprb), pointer :: ww_rho(:) => null() 
  real(kind=jprb), pointer :: ww_zil(:) => null() 
  real(kind=jprb), pointer :: ww_cif(:) => null() 
  real(kind=jprb), pointer :: ww_clk(:) => null() 
  real(kind=jprb), pointer :: ww_ucurw(:) => null() 
  real(kind=jprb), pointer :: ww_vcurw(:) => null() 

  ! Group VX
  real(kind=jprb), pointer :: vx_oro(:) => null() 
  real(kind=jprb), pointer :: vx_tsc(:) => null() 
  real(kind=jprb), pointer :: vx_pws(:) => null() 
  real(kind=jprb), pointer :: vx_pwp(:) => null() 
  real(kind=jprb), pointer :: vx_sno(:) => null() 
  real(kind=jprb), pointer :: vx_tpc(:) => null() 
  real(kind=jprb), pointer :: vx_sab(:) => null() 
  real(kind=jprb), pointer :: vx_xd2(:) => null() 
  real(kind=jprb), pointer :: vx_lsm(:) => null() 
  real(kind=jprb), pointer :: vx_iveg(:) => null() 
  real(kind=jprb), pointer :: vx_arg(:) => null() 
  real(kind=jprb), pointer :: vx_rsmin(:) => null() 
  real(kind=jprb), pointer :: vx_lai(:) => null() 
  real(kind=jprb), pointer :: vx_veg(:) => null() 

  ! * Group VK=VCLIK: Convective cloud pseudo-historic fields:
   real(kind=jprb), pointer :: vk_udgro(:) => null()

  class(field_access),pointer :: dm => null()
  class(field_access),pointer :: dl => null()

  contains

  procedure :: field_map_storage => main_field_map_storage

end type

type, extends(field_metadata_base) :: field_metadata

  contains

  procedure :: field_set_metadata => main_field_set_metadata
  procedure :: field_get_clevtype => main_field_get_clevtype

end type

! We will always work internally with lists of generic fields, identified 
! by their field ID. 
integer(kind=jpim),parameter :: jpnumfids=4200
! (F2)
type type_field_id
!!$  integer(kind=jpim) :: u, v, div,vor,t, pd, vd, nhx, edot,vcasrs,rdlasr,&
!!$   & sgrtl,sgrtm,unl,unl_si,vnl,vnl_si,tnl,tnl_si,spnl,spnl_si,vwvnl,&
!!$   & vdnl_si,pdnl,pdnl_si,gw,curhs,cvrhs,ctrhs,csprhs,cspdrhs,csvdrhs,&
!!$   & cunl,cvnl,ctnl,cspnl,cvwvnl,cspdnl,cupt,cvpt,ctpt,cpdpt,cvdpt,dphi,nhxnl,cnhxnl,&
!!$    & q, o3, l, i, a, ghg, grg, aero, phys, s, r, g, h, lconv,iconv,rconv,sconv,&
!!$    & tke,tte,efb1,efb2,efb3,mxl,src,shtur,fqtur,fstur,cvv,cpf,spf,cvgq,&
!!$    & uom,ual,dom,dal,uen,unebh,rkth,rktqv,rktqc,&
!!$    & lrad,irad,rspec,pre, pref, delp, extra3d,ezdiag,&
!!$    & sp, csppt,dbbc,gws,spnl2,cspnl2,prehyds,u_mean,v_mean,orog, cori, gnordl, gnordm, ccc, lcc, mcc, hcc, tcc, &
!!$    & sb_q, sg_f, sg_a, sg_r, rr_t, rr_w, rr_ic, cl_tcls, &
!!$    & cl_hucls, x2_prwa, x2_prsn, vf_z0f, vf_albf, vf_emisf, vf_getrl, vf_lsm, vf_veg, vf_cvl, vf_co2typ, &
!!$    & vf_cvh, vf_tvl, vf_tvh, vf_lail, vf_laih, vf_soty, vf_ci, vf_ucur, vf_vcur, &
!!$    & vf_aluvp, vf_aluvd, vf_alnip, vf_alnid, &
!!$    & vf_aluvi, vf_aluvv, vf_aluvg, vf_alnii, vf_alniv, vf_alnig, &
!!$    & vv_arg, vv_sab, vv_hv, vv_z0h, vn_top, vn_bas, vn_acpr, vn_accpr, &
!!$    & ws_char, vd_10nu, vd_10nv, vd_upd, vd_z0f, vx_oro, vx_tsc, vx_sno, est, esn, et2, eh2, &
!!$    & ev1, ez, eh, vvel
! From GMV
integer(kind=jpim) :: u=1
integer(kind=jpim) :: v=2
integer(kind=jpim) :: div=3
integer(kind=jpim) :: vor=4
integer(kind=jpim) :: t=5
integer(kind=jpim) :: pd=6
integer(kind=jpim) :: vd=7
integer(kind=jpim) :: nhx=8
integer(kind=jpim) :: edot=9
! ky: GMV fields VCASRS and RDLASR (containing a/r) have been pruned, but how to renumber fields?
integer(kind=jpim) :: vcasrs=10 ! ky: obsolete field (deep-layer models pruned)
integer(kind=jpim) :: rdlasr=11 ! ky: obsolete field (deep-layer models pruned)
! -----------------------------------------------------------------------------------------------
integer(kind=jpim) :: sgrtl=12
integer(kind=jpim) :: sgrtm=13
integer(kind=jpim) :: unl=14
integer(kind=jpim) :: unl_si=15
integer(kind=jpim) :: vnl=16
integer(kind=jpim) :: vnl_si=17
integer(kind=jpim) :: tnl=18
integer(kind=jpim) :: tnl_si=19
integer(kind=jpim) :: spnl=20
integer(kind=jpim) :: spnl_si=21
integer(kind=jpim) :: vwvnl=22
integer(kind=jpim) :: vdnl_si=23
integer(kind=jpim) :: pdnl=24
integer(kind=jpim) :: pdnl_si=25
integer(kind=jpim) :: gw=26
integer(kind=jpim) :: curhs=27
integer(kind=jpim) :: cvrhs=28
integer(kind=jpim) :: ctrhs=29
integer(kind=jpim) :: csprhs=30
integer(kind=jpim) :: cspdrhs=31
integer(kind=jpim) :: csvdrhs=32
integer(kind=jpim) :: cunl=33
integer(kind=jpim) :: cvnl=34
integer(kind=jpim) :: ctnl=35
integer(kind=jpim) :: cspnl=36
integer(kind=jpim) :: cvwvnl=37
integer(kind=jpim) :: cspdnl=38
integer(kind=jpim) :: cupt=39
integer(kind=jpim) :: cvpt=40
integer(kind=jpim) :: ctpt=41
integer(kind=jpim) :: cpdpt=42
integer(kind=jpim) :: cvdpt=43
integer(kind=jpim) :: dphi=44
integer(kind=jpim) :: nhxnl=45
integer(kind=jpim) :: cnhxnl=46
! From gfl
integer(kind=jpim) :: q=47
integer(kind=jpim) :: o3=48 
integer(kind=jpim) :: l=49 
integer(kind=jpim) :: i=50 
integer(kind=jpim) :: a=51 
!integer(kind=jpim) :: ghg=52   ! }
!integer(kind=jpim) :: chem=53  ! } Split into 4xxx
!integer(kind=jpim) :: aero=54  ! }
integer(kind=jpim) :: lrch4=52
integer(kind=jpim) :: phys=55 
integer(kind=jpim) :: s=56
integer(kind=jpim) :: r=57
integer(kind=jpim) :: g=58 
integer(kind=jpim) :: h=59 
integer(kind=jpim) :: lconv=60
integer(kind=jpim) :: iconv=61
integer(kind=jpim) :: rconv=62
integer(kind=jpim) :: sconv=63
integer(kind=jpim) :: tke=64
integer(kind=jpim) :: tte=65
integer(kind=jpim) :: efb1=66
integer(kind=jpim) :: efb2=67
integer(kind=jpim) :: efb3=68
integer(kind=jpim) :: mxl=69
integer(kind=jpim) :: src=70
integer(kind=jpim) :: shtur=71
integer(kind=jpim) :: fqtur=72
integer(kind=jpim) :: fstur=73
integer(kind=jpim) :: cvv=74
integer(kind=jpim) :: cpf=75
integer(kind=jpim) :: spf=76
integer(kind=jpim) :: cvgq=77
integer(kind=jpim) :: uom=78
integer(kind=jpim) :: ual=79
integer(kind=jpim) :: dom=80
integer(kind=jpim) :: dal=81
integer(kind=jpim) :: uen=82
integer(kind=jpim) :: unebh=83
integer(kind=jpim) :: rkth=84
integer(kind=jpim) :: rktqv=85
integer(kind=jpim) :: rktqc=86
integer(kind=jpim) :: lrad=87
integer(kind=jpim) :: irad=88
integer(kind=jpim) :: rspec=89
integer(kind=jpim) :: pre=90
integer(kind=jpim) :: pref=91 
integer(kind=jpim) :: delp=92 
integer(kind=jpim) :: extra3d=93
integer(kind=jpim) :: ezdiag=94
integer(kind=jpim) :: sp=95 
integer(kind=jpim) :: csppt=96
integer(kind=jpim) :: dbbc=97
integer(kind=jpim) :: gws=98
integer(kind=jpim) :: spnl2=99
integer(kind=jpim) :: cspnl2=100
integer(kind=jpim) :: u_mean=101
integer(kind=jpim) :: v_mean=102
integer(kind=jpim) :: orog=103
integer(kind=jpim) :: cori=104
integer(kind=jpim) :: gnordl=105 
integer(kind=jpim) :: gnordm=106
integer(kind=jpim) :: ccc=107 
integer(kind=jpim) :: lcc=108 
integer(kind=jpim) :: mcc=109 
integer(kind=jpim) :: hcc=110 
integer(kind=jpim) :: tcc=111
integer(kind=jpim) :: est=159
integer(kind=jpim) :: esn=160 
integer(kind=jpim) :: et2=161 
integer(kind=jpim) :: eh2=162 
integer(kind=jpim) :: ev1=163 
integer(kind=jpim) :: ez=164
integer(kind=jpim) :: eh=165 
integer(kind=jpim) :: vvel=166
integer(kind=jpim) :: unogw=167
integer(kind=jpim) :: vnogw=168
integer(kind=jpim) :: extcv2d=170
integer(kind=jpim) :: extcv3d=170+jp_necv_2d_max

! ky: sorry, I don't know the rule of numbering this
integer(kind=jpim) :: prehyds=999

  ! Surface (see surface_fields_mix.F90)
  ! Group SB
integer(kind=jpim) :: sb_q=1101
integer(kind=jpim) :: sb_t=1102
integer(kind=jpim) :: sb_tl=1103
  ! Group SG
integer(kind=jpim) :: sg_f=1201 
integer(kind=jpim) :: sg_a=1202
integer(kind=jpim) :: sg_r=1203
integer(kind=jpim) :: sg_t=1204
integer(kind=jpim) :: sg_w=1205
! Group SL
integer(kind=jpim) :: sl_lict=1301
integer(kind=jpim) :: sl_lmlt=1302
integer(kind=jpim) :: sl_ltlt=1303
integer(kind=jpim) :: sl_lblt=1304
integer(kind=jpim) :: sl_lshf=1305
integer(kind=jpim) :: sl_licd=1306
integer(kind=jpim) :: sl_lmld=1307
! Group RR
integer(kind=jpim) :: rr_t=1401
integer(kind=jpim) :: rr_w=1402
integer(kind=jpim) :: rr_fc=1403 
integer(kind=jpim) :: rr_ic=1404
integer(kind=jpim) :: rr_fp1=1405
integer(kind=jpim) :: rr_tir=1406
integer(kind=jpim) :: rr_tmw=1407
!Group CL
integer(kind=jpim) :: cl_tcls=1501
integer(kind=jpim) :: cl_hucls=1502 
!Group OM
integer(kind=jpim) :: om_to=1601
integer(kind=jpim) :: om_so=1602
integer(kind=jpim) :: om_uo=1603
integer(kind=jpim) :: om_vo=1604
! Group EP
integer(kind=jpim) :: ep_ep=1701
! Group X2
integer(kind=jpim) :: x2_x2=1801
! Group CI
integer(kind=jpim) :: ci_ci=1901
! Group VF
integer(kind=jpim) :: vf_z0f=2001 
integer(kind=jpim) :: vf_albf=2002 
integer(kind=jpim) :: vf_emisf=2003
integer(kind=jpim) :: vf_getrl=2004 
integer(kind=jpim) :: vf_lsm=2005
integer(kind=jpim) :: vf_veg=2006
integer(kind=jpim) :: vf_vrlan=2007 
integer(kind=jpim) :: vf_vrldi=2008 
integer(kind=jpim) :: vf_sig=2009 
integer(kind=jpim) :: vf_albsf=2010
integer(kind=jpim) :: vf_lan=2011
integer(kind=jpim) :: vf_sst=2012 
integer(kind=jpim) :: vf_sss=2013 
integer(kind=jpim) :: vf_lz0h=2014 
integer(kind=jpim) :: vf_cvl=2015
integer(kind=jpim) :: vf_cvh=2016
integer(kind=jpim) :: vf_tvl=2017 
integer(kind=jpim) :: vf_tvh=2018 
integer(kind=jpim) :: vf_lail=2019
integer(kind=jpim) :: vf_laih=2020 
integer(kind=jpim) :: vf_soty=2021
integer(kind=jpim) :: vf_clk=2022
integer(kind=jpim) :: vf_dl=2023
integer(kind=jpim) :: vf_ci=2024 
integer(kind=jpim) :: vf_ucur=2025 
integer(kind=jpim) :: vf_vcur=2026
integer(kind=jpim) :: vf_z0rlf=2027
integer(kind=jpim) :: vf_cgpp=2032
integer(kind=jpim) :: vf_crec=2033
integer(kind=jpim) :: vf_sdfor=2036
integer(kind=jpim) :: vf_aluvp=2037 
integer(kind=jpim) :: vf_aluvd=2038 
integer(kind=jpim) :: vf_alnip=2039
integer(kind=jpim) :: vf_alnid=2040
integer(kind=jpim) :: vf_fp1=2041     ! surface orography in the 2nd part of FULLPOS-927
integer(kind=jpim) :: vf_so2dd=2051   ! sulphate dry dep velocity
integer(kind=jpim) :: vf_dmso=2056    ! oceanic dms
integer(kind=jpim) :: vf_urbf=2058   ! soa from co
integer(kind=jpim) :: vf_fca1=2060    ! fraction of calcite over dust 1st bin
integer(kind=jpim) :: vf_fca2=2061    ! fraction of calcite over dust 2st bin
integer(kind=jpim) :: vf_aerdep=2062  ! dust emission potential 
integer(kind=jpim) :: vf_aerlts=2063  ! dust lifting threshold speed 
integer(kind=jpim) :: vf_aerscc=2064  ! dust soil clay content
integer(kind=jpim) :: vf_dsf=2065     ! dust source function
integer(kind=jpim) :: vf_chemflxo=2067 ! total chemistry flux (emissions + deposition) 
integer(kind=jpim) :: vf_chemdv=2068  ! chemistry deposition velocity
integer(kind=jpim) :: vf_nudm=2069    ! nudging mask

! MODIS 6-component albedo coefficients
integer(kind=jpim) :: vf_aluvi=2070
integer(kind=jpim) :: vf_aluvv=2071
integer(kind=jpim) :: vf_aluvg=2072
integer(kind=jpim) :: vf_alnii=2073
integer(kind=jpim) :: vf_alniv=2074
integer(kind=jpim) :: vf_alnig=2075

integer(kind=jpim) :: vf_cur=2076 ! Urban

integer(kind=jpim) :: vf_emis2d=2077 ! 2D emission fields for composition
integer(kind=jpim) :: vf_emis2daux=2078 ! 2D emission auxiliary fields for composition

integer(kind=jpim) :: vf_dsz=2079     ! dust size variation

! Wet and dry deposition fluxes  
integer(kind=jpim) :: vf_chemwdflx=2080 ! wet deposition chemistry flux  
integer(kind=jpim) :: vf_chemddflx=2081 ! dry deposition chemistry flux  

integer(kind=jpim) :: vf_co2typ=2082 ! C3/C4 CO2 photosynthesis type for low vegetation

integer(kind=jpim) :: vf_fwet=2083 ! Wetland

! * Group VP=VCLIP: deep soil diagnostic fields
integer(kind=jpim) :: vp_tpc=2101 ! climatological deep layer temperature
integer(kind=jpim) :: vp_wpc=2102 ! climatological deep layer moisture
! * Group VV=VCLIV: vegetation diagnostic fields:
integer(kind=jpim) :: vv_arg=2201
integer(kind=jpim) :: vv_sab=2202
integer(kind=jpim) :: vv_d2=2203
integer(kind=jpim) :: vv_iveg=2204
integer(kind=jpim) :: vv_rsmin=2205
integer(kind=jpim) :: vv_lai=2206
integer(kind=jpim) :: vv_hv=2207 
integer(kind=jpim) :: vv_z0h=2208 
integer(kind=jpim) :: vv_als=2209
integer(kind=jpim) :: vv_alv=2210
! * Group VN=VCLIN: cloudiness diagnostic predictors:
integer(kind=jpim) :: vn_top=2301
integer(kind=jpim) :: vn_bas=2302
integer(kind=jpim) :: vn_acpr=2303
integer(kind=jpim) :: vn_accpr=2304
integer(kind=jpim) :: vn_accpr5=2305
! * Group VH=VCLIH: convective cloud diagnostic fields:
integer(kind=jpim) :: vh_tcch=2401
integer(kind=jpim) :: vh_scch=2402
integer(kind=jpim) :: vh_bcch=2403
integer(kind=jpim) :: vh_pblh=2404
integer(kind=jpim) :: vh_spsh=2405
integer(kind=jpim) :: vh_qsh =2406
! Group VA=VCLIA: aerosol diagnostic fields:
integer(kind=jpim) :: va_sea=2501
integer(kind=jpim) :: va_lan=2502
integer(kind=jpim) :: va_soo=2503
integer(kind=jpim) :: va_des=2504
integer(kind=jpim) :: va_sul=2505
integer(kind=jpim) :: va_vol=2506
! Group  VG=VCLIG: ice-coupler diagnostic fields
integer(kind=jpim) :: vg_icfr=2601
integer(kind=jpim) :: vg_soup=2602
integer(kind=jpim) :: vg_irup=2603
integer(kind=jpim) :: vg_chss=2604
integer(kind=jpim) :: vg_evap=2605
integer(kind=jpim) :: vg_taux=2606
integer(kind=jpim) :: vg_tauy=2607
!  Group VC=VO3ABC: A,B and C (Climatological ozone profiles) diagnostic fields
integer(kind=jpim) :: vc_a=2701
integer(kind=jpim) :: vc_b=2702
integer(kind=jpim) :: vc_c=2703
! Group V2=VDIAGO2: 2-D climatological/diagnostic fields for an ocean mixed layer model (KPP)
integer(kind=jpim) :: v2_ocdep=2801
integer(kind=jpim) :: v2_ustrc=2802
integer(kind=jpim) :: v2_vstrc=2803
! Group V3=VDIAGO3: 3-D climatological/diagnostic fields for an ocean mixed layer model (KPP):
integer(kind=jpim) :: v3_difm=2901
integer(kind=jpim) :: v3_dift=2902
integer(kind=jpim) :: v3_difs=2903
integer(kind=jpim) :: v3_advt=2904
integer(kind=jpim) :: v3_advs=2905
integer(kind=jpim) :: v3_tri0=2906
integer(kind=jpim) :: v3_tri1=2907
integer(kind=jpim) :: v3_swdk=2908
integer(kind=jpim) :: v3_zo=2909
integer(kind=jpim) :: v3_ho=2910
integer(kind=jpim) :: v3_do=2911
integer(kind=jpim) :: v3_ho_inv=2912
integer(kind=jpim) :: v3_uoc=2913
integer(kind=jpim) :: v3_voc=2914
integer(kind=jpim) :: v3_otke=2915

! * Group VD=VDIAG: (ECMWF) diagnostic fields:
integer(kind=jpim) :: vd_lsp=3001
integer(kind=jpim) :: vd_cp=3002 
integer(kind=jpim) :: vd_sf=3003 
integer(kind=jpim) :: vd_fzra=3004 
integer(kind=jpim) :: vd_bld=3005 
integer(kind=jpim) :: vd_sshf=3006 
integer(kind=jpim) :: vd_slhf=3007 
integer(kind=jpim) :: vd_nee=3008 
integer(kind=jpim) :: vd_gpp=3009 
integer(kind=jpim) :: vd_rec=3010 
integer(kind=jpim) :: vd_msl=3011 
integer(kind=jpim) :: vd_sp=3012 
integer(kind=jpim) :: vd_tcc=3013 
integer(kind=jpim) :: vd_10u=3014 
integer(kind=jpim) :: vd_10v=3015 
integer(kind=jpim) :: vd_2t=3016 
integer(kind=jpim) :: vd_2d=3017 
integer(kind=jpim) :: vd_ssr=3018 
integer(kind=jpim) :: vd_str=3019 
integer(kind=jpim) :: vd_tsr=3020 
integer(kind=jpim) :: vd_ttr=3021
integer(kind=jpim) :: vd_ewss=3022 
integer(kind=jpim) :: vd_nsss=3023 
integer(kind=jpim) :: vd_e=3024 
integer(kind=jpim) :: vd_pev=3025 
integer(kind=jpim) :: vd_ccc=3026 
integer(kind=jpim) :: vd_lcc=3027 
integer(kind=jpim) :: vd_mcc=3028 
integer(kind=jpim) :: vd_hcc=3029 
integer(kind=jpim) :: vd_lgws=3030 
integer(kind=jpim) :: vd_mgws=3031 
integer(kind=jpim) :: vd_gwd=3032 
integer(kind=jpim) :: vd_mx2t=3033 
integer(kind=jpim) :: vd_mn2t=3034 
integer(kind=jpim) :: vd_mx2t6=3035 
integer(kind=jpim) :: vd_mn2t6=3036 
integer(kind=jpim) :: vd_ro=3037 
integer(kind=jpim) :: vd_sro=3038 
integer(kind=jpim) :: vd_ssro=3039 
integer(kind=jpim) :: vd_alb=3040 
integer(kind=jpim) :: vd_iewss=3041 
integer(kind=jpim) :: vd_insss=3042 
integer(kind=jpim) :: vd_isshf=3043 
integer(kind=jpim) :: vd_ie=3044 
integer(kind=jpim) :: vd_inee=3045 
integer(kind=jpim) :: vd_igpp=3046 
integer(kind=jpim) :: vd_irec=3047 
integer(kind=jpim) :: vd_csf=3048 
integer(kind=jpim) :: vd_lssf=3049 
integer(kind=jpim) :: vd_mxtpr=3050 
integer(kind=jpim) :: vd_mntpr=3051 
integer(kind=jpim) :: vd_mxtpr6=3052 
integer(kind=jpim) :: vd_mntpr6=3053 
integer(kind=jpim) :: vd_tpr=3054 
integer(kind=jpim) :: vd_lsrr=3055 
integer(kind=jpim) :: vd_crr=3056
integer(kind=jpim) :: vd_lssfr=3057 
integer(kind=jpim) :: vd_csfr=3058
integer(kind=jpim) :: vd_ptype=3059 
integer(kind=jpim) :: vd_ilspf=3060 
integer(kind=jpim) :: vd_z0f=3061 
integer(kind=jpim) :: vd_lz0h=3062 
integer(kind=jpim) :: vd_tcw=3063 
integer(kind=jpim) :: vd_tcwv=3064 
integer(kind=jpim) :: vd_tclw=3065 
integer(kind=jpim) :: vd_tciw=3066 
integer(kind=jpim) :: vd_tcrw=3067 
integer(kind=jpim) :: vd_tcsw=3068 
integer(kind=jpim) :: vd_tcslw=3069 
integer(kind=jpim) :: vd_ssrd=3070
integer(kind=jpim) :: vd_strd=3071 
integer(kind=jpim) :: vd_ssrdc=3072 
integer(kind=jpim) :: vd_strdc=3073 
integer(kind=jpim) :: vd_blh=3074 
integer(kind=jpim) :: vd_sund=3075 
integer(kind=jpim) :: vd_spar=3076 
integer(kind=jpim) :: vd_suvb=3077 
integer(kind=jpim) :: vd_sfdir=3078 
integer(kind=jpim) :: vd_scdir=3079 
integer(kind=jpim) :: vd_sdsrp=3080 
integer(kind=jpim) :: vd_cape=3081
integer(kind=jpim) :: vd_capes=3082
integer(kind=jpim) :: vd_tsrc=3083
integer(kind=jpim) :: vd_ttrc=3084 
integer(kind=jpim) :: vd_ssrc=3085 
integer(kind=jpim) :: vd_strc=3086 
integer(kind=jpim) :: vd_es=3087 
integer(kind=jpim) :: vd_smlt=3088 
integer(kind=jpim) :: vd_10fg=3089 
integer(kind=jpim) :: vd_10fg6=3090 
integer(kind=jpim) :: vd_10fgcv=3091 
integer(kind=jpim) :: vd_i10fg=3092
integer(kind=jpim) :: vd_lspf=3093
integer(kind=jpim) :: vd_tco3=3094 
integer(kind=jpim) :: vd_vimd=3095 
integer(kind=jpim) :: vd_sparc=3096 
integer(kind=jpim) :: vd_stinc=3097 
integer(kind=jpim) :: vd_cbase=3098 
integer(kind=jpim) :: vd_0degl=3099 
integer(kind=jpim) :: vd_visih=3100 
integer(kind=jpim) :: vd_cin=3101 
integer(kind=jpim) :: vd_kindex=3102 
integer(kind=jpim) :: vd_ttindex=3103 
integer(kind=jpim) :: vd_cbasea=3104 
integer(kind=jpim) :: vd_ctopc=3105 
integer(kind=jpim) :: vd_ztwetb0=3106 
integer(kind=jpim) :: vd_ztwetb1=3107 
integer(kind=jpim) :: vd_tcghg=3108
integer(kind=jpim) :: vd_tcchem=3109
integer(kind=jpim) :: vd_aerodiag=3110
integer(kind=jpim) :: vd_aero_wvl_diag=3111
integer(kind=jpim) :: vd_100u=3112 
integer(kind=jpim) :: vd_100v=3113 
integer(kind=jpim) :: vd_zust=3114 
integer(kind=jpim) :: vd_10nu=3115 
integer(kind=jpim) :: vd_10nv=3116 
integer(kind=jpim) :: vd_dndzn=3117 
integer(kind=jpim) :: vd_dndza=3118 
integer(kind=jpim) :: vd_dctb=3119
integer(kind=jpim) :: vd_tplb=3120 
integer(kind=jpim) :: vd_tplt=3121 
integer(kind=jpim) :: vd_odss=3122 
integer(kind=jpim) :: vd_oddu=3123 
integer(kind=jpim) :: vd_odom=3124 
integer(kind=jpim) :: vd_odbc=3125 
integer(kind=jpim) :: vd_odsu=3126
integer(kind=jpim) :: vd_odni=3127
integer(kind=jpim) :: vd_odam=3128 
integer(kind=jpim) :: vd_odsoa=3129
integer(kind=jpim) :: vd_odvfa=3130 
integer(kind=jpim) :: vd_odvsu=3131 
integer(kind=jpim) :: vd_aepm1=3132 
integer(kind=jpim) :: vd_aepm25=3133 
integer(kind=jpim) :: vd_aepm10=3134 
integer(kind=jpim) :: vd_uvbed=3135
integer(kind=jpim) :: vd_uvbedcs=3136 
integer(kind=jpim) :: vd_litoti=3137
integer(kind=jpim) :: vd_licgi=3138
integer(kind=jpim) :: vd_litota6=3139
integer(kind=jpim) :: vd_licga6=3140
integer(kind=jpim) :: vd_200u=3141 
integer(kind=jpim) :: vd_200v=3142 
integer(kind=jpim) :: vd_2sh=3143
integer(kind=jpim) :: vd_odtoacc=3144
integer(kind=jpim) :: vd_m10degl=3145
integer(kind=jpim) :: vd_mucape=3146
integer(kind=jpim) :: vd_pdepl=3147
integer(kind=jpim) :: vd_mlcape50=3148
integer(kind=jpim) :: vd_mlcape100=3149
integer(kind=jpim) :: vd_mlcin50=3150
integer(kind=jpim) :: vd_mlcin100=3151
integer(kind=jpim) :: vd_tropotp=3152
integer(kind=jpim) :: vd_sdsl=3153
integer(kind=jpim) :: vd_ich4=3154
integer(kind=jpim) :: vd_ptypeocc6=3155

! Group SM=SATSIM: (ECMWF) simulated satellite images
integer(kind=jpim) :: sm_clbt=3201
integer(kind=jpim) :: sm_csbt=3202

! Group WS=WAVES: surface prognostic quantities over sea (used by IFS)
integer(kind=jpim) :: ws_char=3301
integer(kind=jpim) :: ws_charhq=3302
integer(kind=jpim) :: ws_ustokes=3303
integer(kind=jpim) :: ws_vstokes=3304
integer(kind=jpim) :: ws_tauocx=3305
integer(kind=jpim) :: ws_tauocy=3306
integer(kind=jpim) :: ws_phioc=3307
integer(kind=jpim) :: ws_wsemean=3308
integer(kind=jpim) :: ws_wsfmean=3309

! * Group VX=VCLIX: auxilary climatological diagnostic fields:
integer(kind=jpim) :: vx_oro=3401
integer(kind=jpim) :: vx_tsc=3402
integer(kind=jpim) :: vx_pws=3403
integer(kind=jpim) :: vx_pwp=3404
integer(kind=jpim) :: vx_sno=3405
integer(kind=jpim) :: vx_tpc=3406
integer(kind=jpim) :: vx_sab=3407
integer(kind=jpim) :: vx_xd2=3408
integer(kind=jpim) :: vx_lsm=3409
integer(kind=jpim) :: vx_iveg=3410
integer(kind=jpim) :: vx_arg=3411
integer(kind=jpim) :: vx_rsmin=3412
integer(kind=jpim) :: vx_lai=3413
integer(kind=jpim) :: vx_veg=3414
! * Group VK=VCLIK: Convective cloud pseudo-historic fields:
integer(kind=jpim) :: vk_udgro=3501

! * Dynamically-assigned COMPO 3D fields
!integer(kind=jpim) :: compo_3d_first = 3601
!integer(kind=jpim) :: compo_3d_last = 3800
  !  Group WW 

integer(kind=jpim) :: ww_u10n=3801 
integer(kind=jpim) :: ww_v10n=3802 
integer(kind=jpim) :: ww_rho= 3803
integer(kind=jpim) :: ww_zil= 3804
integer(kind=jpim) :: ww_cif= 3805
integer(kind=jpim) :: ww_clk= 3806
integer(kind=jpim) :: ww_ucurw=3807 
integer(kind=jpim) :: ww_vcurw= 3808

end type type_field_id

! (F3)
type type_dynamic_names
  character(len=JP_NAME_MAX_LEN) :: cname=''
  character(len=JP_COMMENTS_MAX_LEN) :: clongname=''
  integer(kind=jpim) :: igribcode
  integer(kind=jpim) :: ifid
end type type_dynamic_names

! Dynamic namespace (currently for COMPO)
integer(kind=jpim),parameter :: compo_3d_first = 3601
integer(kind=jpim),parameter :: compo_3d_last = 3800
integer(kind=jpim)           :: ndyn_names=0
type(type_dynamic_names) ::  dyn_nam(compo_3d_last-compo_3d_first+1)

type(type_field_id), parameter :: fid=type_field_id()
type(field_metadata) :: main_field_metadata

! 2nd dimension types
integer(kind=jpim), parameter :: d2none      = 0 ! no second dimension
integer(kind=jpim), parameter :: d2full      = 1 ! full model levels
integer(kind=jpim), parameter :: d2half      = 2 ! half model levels
integer(kind=jpim), parameter :: d2soil      = 3 ! soil model levels
integer(kind=jpim), parameter :: d2snow      = 4 ! snow model levels
integer(kind=jpim), parameter :: d2om        = 5 ! Ocean model "levels"
integer(kind=jpim), parameter :: d2chemflxo  = 6
integer(kind=jpim), parameter :: d2chemwdflx = 7
integer(kind=jpim), parameter :: d2chemddflx = 8
integer(kind=jpim), parameter :: d2chemdv    = 9
integer(kind=jpim), parameter :: d2emis2d    = 10
integer(kind=jpim), parameter :: d2emis2daux = 11

integer(kind=jpim), parameter :: nleveltypes_main = 11

! 3rd dimension info - AJGDB where does this come from when it's used for real?

integer(kind=jpim),parameter :: ndim3types_main=2
integer(kind=jpim), parameter :: d3extra3d = 1
integer(kind=jpim), parameter :: d3ezdiag  = 2

!integer(kind=jpim), parameter :: nlev_3d = 2

#ifndef FIELD_MOD_TEST
#include "abor1.intfb.h"
#endif

contains

! ---------------------------------------------------------
! (Essentially a constructor) Initialises names, dimensions 
! and other attributes
! ---------------------------------------------------------
subroutine main_field_set_metadata(this,metadata,kleveltypes,kdim3types)
class(field_metadata),    intent(in)    :: this
type(type_fvar), ALLOCATABLE, intent(inout) :: metadata(:)
integer(kind=jpim),       intent(  out) :: kleveltypes
integer(kind=jpim),       intent(  out) :: kdim3types
REAL(KIND=JPHOOK)  :: ZHOOK_HANDLE
integer(kind=jpim) :: jfid, iecv, jname
character(len=jp_name_max_len) :: clname

IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:MAIN_FIELD_SET_METADATA',0,ZHOOK_HANDLE)
kleveltypes = nleveltypes_main
kdim3types  = ndim3types_main
! (F5) Increment the total number of fields here
allocate(metadata(jpnumfids))

! Main field details and attributes (F6)

! GMV
call set_fvar(metadata,fid%u,'u',2,d2full,'Horizontal wind component',1,KGRIB=NGRBU)
call set_fvar(metadata,fid%v,'v',2,d2full,'Horizontal wind component',1,KGRIB=NGRBV)
call set_fvar(metadata,fid%div,'div',2,d2full,'Divergence',1,KGRIB=NGRBD)
call set_fvar(metadata,fid%vor,'vor',2,d2full,'Vorticity',1,KGRIB=NGRBVO)
call set_fvar(metadata,fid%t,'t',2,d2full,'Temperature',1,KGRIB=NGRBT)
call set_fvar(metadata,fid%pd,'pd',2,d2full,'Pressure Departure (NH only)',1)
call set_fvar(metadata,fid%vd,'vd',2,d2full,'Vertical Divergence (NH only)',1)
call set_fvar(metadata,fid%nhx,'nhx',2,d2full,'NHX term (NH only)',1)
call set_fvar(metadata,fid%edot,'edot',2,d2full,'?',1)
call set_fvar(metadata,fid%vcasrs,'vcasrs',2,d2full,'?',1) ! ky: obsolete field (deep-layer models pruned)
call set_fvar(metadata,fid%rdlasr,'rdlasr',2,d2full,'?',1) ! ky: obsolete field (deep-layer models pruned)
call set_fvar(metadata,fid%sgrtl,'sgrtl',2,d2full,'?',1)
call set_fvar(metadata,fid%sgrtm,'sgrtm',2,d2full,'?',1)
call set_fvar(metadata,fid%unl,'unl',2,d2full,'?',1)
call set_fvar(metadata,fid%unl_si,'unl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%vnl,'vnl',2,d2full,'?',1)
call set_fvar(metadata,fid%vnl_si,'vnl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%tnl,'tnl',2,d2full,'?',1)
call set_fvar(metadata,fid%tnl_si,'tnl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%spnl,'spnl',2,d2full,'?',1)
call set_fvar(metadata,fid%spnl_si,'spnl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%vwvnl,'vwvnl',2,d2full,'?',1)
call set_fvar(metadata,fid%vdnl_si,'vdnl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%pdnl,'pdnl',2,d2full,'?',1)
call set_fvar(metadata,fid%pdnl_si,'pdnl_si',2,d2full,'?',1)
call set_fvar(metadata,fid%gw,'gw',2,d2full,'?',1)
call set_fvar(metadata,fid%curhs,'curhs',2,d2full,'?',1)
call set_fvar(metadata,fid%cvrhs,'cvrhs',2,d2full,'?',1)
call set_fvar(metadata,fid%ctrhs,'ctrhs',2,d2full,'?',1)
call set_fvar(metadata,fid%cspdrhs,'cspdrhs',2,d2full,'?',1)
call set_fvar(metadata,fid%csvdrhs,'csvdrhs',2,d2full,'?',1)
call set_fvar(metadata,fid%cunl,'cunl',2,d2full,'?',1)
call set_fvar(metadata,fid%cvnl,'cvnl',2,d2full,'?',1)
call set_fvar(metadata,fid%ctnl,'ctnl',2,d2full,'?',1)
call set_fvar(metadata,fid%cspnl,'cspnl',2,d2full,'?',1)
call set_fvar(metadata,fid%cvwvnl,'cvwvnl',2,d2full,'?',1)
call set_fvar(metadata,fid%cspdnl,'cspdnl',2,d2full,'?',1)
call set_fvar(metadata,fid%cupt,'cupt',2,d2full,'?',1)
call set_fvar(metadata,fid%cvpt,'cvpt',2,d2full,'?',1)
call set_fvar(metadata,fid%ctpt,'ctpt',2,d2full,'?',1)
call set_fvar(metadata,fid%cpdpt,'cpdpt',2,d2full,'?',1)
call set_fvar(metadata,fid%cvdpt,'cvdpt',2,d2full,'?',1)
call set_fvar(metadata,fid%dphi,'dphi',2,d2full,'?',1)
call set_fvar(metadata,fid%nhxnl,'nhxnl',2,d2full,'?',1)
call set_fvar(metadata,fid%cnhxnl,'cnhxnl',2,d2full,'?',1)


! GFL
call set_fvar(metadata,fid%q,   'q',    2,d2full,'Specific humidity',3,PR=RV,PRCP=RCPV,KGRIB=NGRBQ)
call set_fvar(metadata,fid%o3,  'o3',   2,d2full,'Ozone',3,KGRIB=NGRBO3)
call set_fvar(metadata,fid%l,   'l',    2,d2full,'Cloud liquid water content',3,PR=0.0_JPRB,PRCP=RCW,KGRIB=NGRBCLWC)
call set_fvar(metadata,fid%i,   'i',    2,d2full,'Cloud ice content',3,PR=0.0_JPRB,PRCP=RCS,KGRIB=NGRBCIWC)
call set_fvar(metadata,fid%a,   'a',    2,d2full,'Cloud fraction',3,KGRIB=NGRBCC)
!call set_fvar(metadata,fid%ghg, 'ghg',  2,d2full,'Greenhouse gases',3)         ! NFTBC
!call set_fvar(metadata,fid%grg, 'grg',  2,d2full,'Reactive gases',3)           ! NFTBC
!call set_fvar(metadata,fid%aero,'aero', 2,d2full,'Aerosols',3)                 ! NFTBC
call set_fvar(metadata,fid%lrch4,   'lrch4',    2,d2full,'CH4 loss rate',3,KGRIB=210071)
call set_fvar(metadata,fid%phys,'phys ',2,d2full,'Diagnostic moist physics',3) ! NFTBC
call set_fvar(metadata,fid%s,   's',    2,d2full,'Snow water content',3,PR=0.0_JPRB,PRCP=RCS,KGRIB=NGRBCSWC)
call set_fvar(metadata,fid%r,   'r',    2,d2full,'Rain water content',3,PR=0.0_JPRB,PRCP=RCW,KGRIB=NGRBCRWC)
call set_fvar(metadata,fid%g,   'g',    2,d2full,'Graupel water content',3,PR=0.0_JPRB,PRCP=RCS)
call set_fvar(metadata,fid%h,   'h',    2,d2full,'Hail water content',3,PR=0.0_JPRB,PRCP=RCS)
call set_fvar(metadata,fid%lconv,   'lconv',    2,d2full,'Liquid water (CONV. PART)',3,PR=0.0_JPRB,PRCP=RCW)
call set_fvar(metadata,fid%iconv,   'iconv',    2,d2full,'Ice    water (CONV. PART)',3,PR=0.0_JPRB,PRCP=RCS)
call set_fvar(metadata,fid%rconv,   'rconv',    2,d2full,'Rain         (CONV. PART)',3,PR=0.0_JPRB,PRCP=RCW)
call set_fvar(metadata,fid%sconv,   'sconv',    2,d2full,'Snow         (CONV. PART)',3,PR=0.0_JPRB,PRCP=RCS)
call set_fvar(metadata,fid%tke,   'tke',    2,d2full,'Total kinetic energy',3)
call set_fvar(metadata,fid%tte,   'tte',    2,d2full,'Turbulent Total Energy',3)
call set_fvar(metadata,fid%efb1,   'efb1',    2,d2full,'First variable EFB scheme',3)
call set_fvar(metadata,fid%efb2,   'efb2',    2,d2full,'Second variable EFB scheme',3)
call set_fvar(metadata,fid%efb3,   'efb3',    2,d2full,'Third variable EFB scheme',3)
call set_fvar(metadata,fid%mxl,   'mxl',    2,d2full,'Prognostic mixing length',3)
call set_fvar(metadata,fid%src,   'src',    2,d2full,'Second-order flux for AROME s"rc"/2Sigma_s2 multiplied by Lambda_3',3)
call set_fvar(metadata,fid%shtur,   'shtur',    2,d2full,'Shear source term for turbulence',3)
call set_fvar(metadata,fid%fqtur,   'fqtur',    2,d2full,'Flux form source term for turbulence -moisture',3)
call set_fvar(metadata,fid%fstur,   'fstur',    2,d2full,'Flux form source term for turbulence -enthalpy',3)
call set_fvar(metadata,fid%cvv,   'cvv',    2,d2full,'Convective Vertical Velocity',3)
call set_fvar(metadata,fid%rkth,   'rkth',    2,d2full,'Rasch-Kristjansson H tendency',3)
call set_fvar(metadata,fid%rktqv,   'rktqv',    2,d2full,'Rasch-Kristjansson Qv tendency',3)
call set_fvar(metadata,fid%rktqc,   'rktqc',    2,d2full,'Rasch-Kristjansson Qc tendency',3)
call set_fvar(metadata,fid%cpf,   'cpf',    2,d2full,'Convective Precipitation Flux',3)
call set_fvar(metadata,fid%spf,   'spf',    2,d2full,' Stratiform Precipitation Flux',3)
call set_fvar(metadata,fid%cvgq,   'cvgq',    2,d2full,'Moisture Convergence for french physics',3)
call set_fvar(metadata,fid%uom,   'uom',    2,d2full,'Updraught vert velocity ',3)
call set_fvar(metadata,fid%ual,   'ual',    2,d2full,'Updraught mesh fraction ',3)
call set_fvar(metadata,fid%dom,   'dom',    2,d2full,'Downdraught vert velocity ',3)
call set_fvar(metadata,fid%dal,   'dal',    2,d2full,'Downdraught mesh fraction ',3)
call set_fvar(metadata,fid%uen,   'uen',    2,d2full,'Updraught entrainment ',3)
call set_fvar(metadata,fid%unebh,   'unebh',    2,d2full,'pseudo-historic convective ',3)
call set_fvar(metadata,fid%lrad,   'lrad',    2,d2full,'Radiative cloud Liquid water',3)
call set_fvar(metadata,fid%irad,   'irad',    2,d2full,'Radiative cloud Ice water',3)
call set_fvar(metadata,fid%rspec,   'rspec',    2,d2full,'Specific gas constant',3)
call set_fvar(metadata,fid%extra3d,  'extra3d',   3,d2full,'Extra 3D fields',3,kdim3_type=d3extra3d)
call set_fvar(metadata,fid%ezdiag,   'ezdiag',    3,d2full,'Easy diagnostics?',3,kdim3_type=d3ezdiag)
call set_fvar(metadata,fid%unogw,   'unogw',    2,d2full,'U gw-momentum tendencies',3,KGRIB=228134)
call set_fvar(metadata,fid%vnogw,   'vnogw',    2,d2full,'V gw-momentum tendencies',3,KGRIB=228136)
! Pressure
call set_fvar(metadata,fid%pre, 'pre', 2,d2half,'Half level pressure',0)
call set_fvar(metadata,fid%pref,'pref',2,d2full,'Full level pressure',0)
call set_fvar(metadata,fid%delp,'delp',2,d2full,'Pressure difference across layers',0)
! GMVS
call set_fvar(metadata,fid%sp,'sp',1,d2none,'Surface pressure (GMVS)',2,KGRIB=NGRBLNSP)
call set_fvar(metadata,fid%csppt,'csppt',1,d2none,'?',2)
call set_fvar(metadata,fid%dbbc,'dbbc',1,d2none,'?',2)
call set_fvar(metadata,fid%gws,'gws',1,d2none,'?',2)
call set_fvar(metadata,fid%csprhs,'csprhs',1,d2none,'?',2)
call set_fvar(metadata,fid%cspnl2,'cspnl2',1,d2none,'?',2)
call set_fvar(metadata,fid%spnl2,'spnl2',1,d2none,'?',2)
call set_fvar(metadata,fid%prehyds,'prehyds',1,d2none,'?',2)
do iecv=1, jp_necv_2d_max
  jfid=fid%extcv2d+iecv-1
  write (clname, fmt='("extcv2d_var_ ",I3.3)') jfid
  call set_fvar(metadata,jfid,clname,1,d2none,'Extended control variable 2D',0)
enddo
do iecv=1, jp_necv_3d_max
  jfid=fid%extcv3d+iecv-1
  write (clname, fmt='("extcv3d_var_ ",I3.3)') jfid
  call set_fvar(metadata,jfid,clname,2,d2full,'Extended control variable 3D',0)
enddo

call set_fvar(metadata,fid%u_mean,'u_mean',0,d2full,'?',13)
call set_fvar(metadata,fid%v_mean,'v_mean',0,d2full,'?',13)
! Constants
call set_fvar(metadata,fid%orog,  'orog',  1,d2none,'Orography',0) 
call set_fvar(metadata,fid%cori,  'cori',  1,d2none,'Coriolis parameter',0) 
call set_fvar(metadata,fid%gnordl,'gnordl',1,d2none,'Zonal component of vector directed towards geographic North',0)
call set_fvar(metadata,fid%gnordm,'gnordm',1,d2none,'Meridional component of vector directed towards geographic North',0) 
! Cloud
call set_fvar(metadata,fid%ccc,   'ccc',   1,d2none,'Convective cloud cover',0) 
call set_fvar(metadata,fid%lcc,   'lcc',   1,d2none,'Low cloud cover',0) 
call set_fvar(metadata,fid%mcc,   'mcc',   1,d2none,'Middle cloud cover',0) 
call set_fvar(metadata,fid%hcc,   'hcc',   1,d2none,'High cloud cover',0) 
call set_fvar(metadata,fid%tcc,   'tcc',   1,d2none,'Total cloud cover',0) 
! Surface
! Group SB
call set_fvar(metadata,fid%sb_t,'sb_t',2,d2soil,'temperature (group SB=SOILB)',11)
call set_fvar(metadata,fid%sb_q,'sb_q',2,d2soil,'Liquid water content (group SB=SOILB)',11)
call set_fvar(metadata,fid%sb_tl,'sb_tl',2,d2soil,'ice water content (for MF)(group SB=SOILB)',11)
! Group SG
call set_fvar(metadata,fid%sg_f,'sg_f',2,d2snow,'Content of surface snow (group SG=SNOWG)',12,kgrib=NGRBSD)
call set_fvar(metadata,fid%sg_a,'sg_a',2,d2snow,'Snow albedo (group SG=SNOWG)',12,kgrib=NGRBASN)
call set_fvar(metadata,fid%sg_r,'sg_r',2,d2snow,'Snow density (group SG=SNOWG)',12,kgrib=NGRBRSN)
call set_fvar(metadata,fid%sg_t,'sg_t',2,d2snow,'Total albedo (diagnostic for MF for LVGSN) (group SG=SNOWG)',12,kgrib=NGRBTSN)
call set_fvar(metadata,fid%sg_w,'sg_w',2,d2snow,'Snow liquid water content (group SG=SNOWG)',12,kgrib=NGRBWSN)
!Group SL
call set_fvar(metadata,fid%sl_lict,'sl_lict',1,d2none,'lake ice temperature',13,KGRIB=NGRBLICT)
call set_fvar(metadata,fid%sl_lmlt,'sl_lmlt',1,d2none,'lake mixed-layer temperature',13,KGRIB=NGRBLMLT)
call set_fvar(metadata,fid%sl_ltlt,'sl_ltlt',1,d2none,'lake total layer temperature',13,KGRIB=NGRBLTLT)
call set_fvar(metadata,fid%sl_lblt,'sl_lblt',1,d2none,'lake bottom layer temperature',13,KGRIB=NGRBLBLT)
call set_fvar(metadata,fid%sl_lshf,'sl_lshf',1,d2none,'lake shape factor',13,KGRIB=NGRBLSHF)
call set_fvar(metadata,fid%sl_licd,'sl_licd',1,d2none,'lake ice depth',13,KGRIB=NGRBLICD)
call set_fvar(metadata,fid%sl_lmld,'sl_lmld',1,d2none,'lake mixed-layer depth',13,KGRIB=NGRBLMLD)

! Group RR
call set_fvar(metadata,fid%rr_t, 'rr_t', 1,d2none,'Skin temperature (Ts) (group RR=RESVR)',14,KGRIB=NGRBSKT)
call set_fvar(metadata,fid%rr_w, 'rr_w', 1,d2none,'Superficial reservoir water content (Ws) (group RR=RESVR)',14,KGRIB=NGRBSRC)
call set_fvar(metadata,fid%rr_fc,'rr_fc',1,d2none,'skin water content (Wl) at MF',14)
call set_fvar(metadata,fid%rr_ic,'rr_ic',1,d2none,'Superficial reservoir ice (group RR=RESVR)',14)
call set_fvar(metadata,fid%rr_fp1,'rr_fp1',1,d2none,'interpolated Ts for 2nd part of 927-FULLPOS',14)
call set_fvar(metadata,fid%rr_tir, 'rr_tir', 1,d2none,'Skin temperature (Ts) (IR)',14)
call set_fvar(metadata,fid%rr_tmw, 'rr_tmv', 1,d2none,'Skin temperature (Ts) (MW)',14)
! Group CL
call set_fvar(metadata,fid%cl_tcls, 'cl_tcls', 1,d2none,'2m temperature (group CL=CLS)',15)
call set_fvar(metadata,fid%cl_hucls,'cl_hucls',1,d2none,'2m humidity (group CL=CLS)',15)
! Group OM=OML
call set_fvar(metadata,fid%om_to, 'om_to', 2,d2om,'OML temperature',16)
call set_fvar(metadata,fid%om_so, 'om_so', 2,d2om,'OML salinity',16)
call set_fvar(metadata,fid%om_uo, 'om_uo', 2,d2om,'OML U velocity',16)
call set_fvar(metadata,fid%om_vo, 'om_vo', 2,d2om,'OML V velocity',16)
! Group EP
call set_fvar(metadata,fid%ep_ep,'ep_ep',1,d2none,'',17)
! Group X2
call set_fvar(metadata,fid%x2_x2,'x2_x2',1,d2none,'',18)
! Group CI
call set_fvar(metadata,fid%ci_ci,'ci_ci',1,d2none,'',19)
! Group VF
call set_fvar(metadata,fid%vf_z0f,  'vf_z0f',  1,d2none,'Gravity times surface roughness length (group VF=VARSF)',20,KGRIB=NGRBSR)
call set_fvar(metadata,fid%vf_albf, 'vf_albf', 1,d2none,'Surface shortwave albedo (group VF=VARSF)',20,KGRIB=NGRBAL)
call set_fvar(metadata,fid%vf_emisf,'vf_emisf',1,d2none,'Surface longwave emissivity (group VF=VARSF)',20)
call set_fvar(metadata,fid%vf_getrl,'vf_getrl',1,d2none,'Standard devaition of orography (group VF=VARSF)',20,KGRIB=NGRBSDOR)
call set_fvar(metadata,fid%vf_lsm,  'vf_lsm',  1,d2none,'Land-sea mask (group VF=VARSF)',20,KGRIB=NGRBLSM)
call set_fvar(metadata,fid%vf_veg,  'vf_veg',  1,d2none,'Vegetation cover (group VF=VARSF)',20)
call set_fvar(metadata,fid%vf_vrlan,  'vf_vvrlan',  1,d2none,'anisotropy of the sub-grid scale orography',20,KGRIB=NGRBISOR)
call set_fvar(metadata,fid%vf_vrldi,  'vf_vvrldi',  1,d2none,'angle of the direction of orography with the x axis',20,KGRIB=NGRBANOR)
call set_fvar(metadata,fid%vf_sig,  'vf_sig',  1,d2none,'characteristic orographic slope',20,KGRIB=NGRBSLOR)
call set_fvar(metadata,fid%vf_albsf,  'vf_albsf',  1,d2none,'soil shortwave albedo',20)
call set_fvar(metadata,fid%vf_lan,  'vf_lan',  1,d2none,'fraction of land',20)
call set_fvar(metadata,fid%vf_sst,  'vf_sst',  1,d2none,'(open) sea surface temperature',20,KGRIB=NGRBSSTK)
call set_fvar(metadata,fid%vf_sss,  'vf_sss',  1,d2none,'sea surface salinity',20,KGRIB=NGRBSSS)
call set_fvar(metadata,fid%vf_lz0h,  'vf_lz0h',  1,d2none,'logarithm of roughness length for heat',20,KGRIB=NGRBLSRH)
call set_fvar(metadata,fid%vf_cvl,  'vf_cvl',  1,d2none,'Low vegetation cover (group VF=VARSF)',20,KGRIB=NGRBCVL)
call set_fvar(metadata,fid%vf_co2typ,'vf_co2typ', 1,d2none,'CO2 photosynthesis type for low vegetation cover (group VF=VARSF)',20,KGRIB=NGRBCO2TYP)
call set_fvar(metadata,fid%vf_cvh,  'vf_cvh',  1,d2none,'High vegetation cover (group VF=VARSF)',20,KGRIB=NGRBCVH)
call set_fvar(metadata,fid%vf_fwet, 'vf_fwet',  1,d2none,'Wetland fraction (group VF=VARSF)',20,KGRIB=NGRBFWET)
call set_fvar(metadata,fid%vf_cur,  'vf_cur',  1,d2none,'Urban cover (group VF=VARSF)',20,KGRIB=NGRBCUR)
call set_fvar(metadata,fid%vf_tvl,  'vf_tvl',  1,d2none,'Low vegetation type (group VF=VARSF)',20,KGRIB=NGRBTVL)
call set_fvar(metadata,fid%vf_tvh,  'vf_tvh',  1,d2none,'High vegetation type (group VF=VARSF)',20,KGRIB=NGRBTVH)
call set_fvar(metadata,fid%vf_lail, 'vf_lail', 1,d2none,'Low vegetation LAI (group VF=VARSF)',20,KGRIB=NGRBLAIL)
call set_fvar(metadata,fid%vf_laih, 'vf_laih', 1,d2none,'High vegetation LAI (group VF=VARSF)',20,KGRIB=NGRBLAIH)
call set_fvar(metadata,fid%vf_soty, 'vf_soty', 1,d2none,'Soil type (group VF=VARSF)',20,KGRIB=NGRBSLT)
call set_fvar(metadata,fid%vf_clk,  'vf_clk',  1,d2none,'lake cover',20,KGRIB=NGRBCL)
call set_fvar(metadata,fid%vf_dl,   'vf_dl',  1,d2none, 'lake depth',20,KGRIB=NGRBDL)
call set_fvar(metadata,fid%vf_ci,   'vf_ci',   1,d2none,'Sea ice fraction (group VF=VARSF)',20,KGRIB=NGRBCI)
call set_fvar(metadata,fid%vf_ucur, 'vf_ucur', 1,d2none,'U-component of the ocean current (group VF=VARSF)',20,KGRIB=NGRBUCUR)
call set_fvar(metadata,fid%vf_vcur, 'vf_vcur', 1,d2none,'V-component of the ocean current (group VF=VARSF)',20,KGRIB=NGRBVCUR)
call set_fvar(metadata,fid%vf_z0rlf,'vf_z0rlf',1,d2none,'gravity * vegetation roughness length',20)
call set_fvar(metadata,fid%vf_cgpp,'vf_cgpp',1,d2none,'GPP bias correction factor',20)
call set_fvar(metadata,fid%vf_crec,'vf_crec',1,d2none,'REC bias correction factor',20)
call set_fvar(metadata,fid%vf_sdfor,'vf_sdfor',1,d2none,'SD filtered orography',20,KGRIB=NGRBSDFOR)
call set_fvar(metadata,fid%vf_aluvp,'vf_aluvp',1,d2none,'MODIS-derived parallel albedo for UV/vis radiation (group VF=VARSF)',20,KGRIB=NGRBALUVP)
call set_fvar(metadata,fid%vf_aluvd,'vf_aluvd',1,d2none,'MODIS-derived diffuse albedo for UV/vis radiation (group VF=VARSF)',20,KGRIB=NGRBALUVD)
call set_fvar(metadata,fid%vf_alnip,'vf_alnip',1,d2none,'MODIS-derived parallel albedo for near-IR radiation (group VF=VARSF)',20,KGRIB=NGRBALNIP)
call set_fvar(metadata,fid%vf_alnid,'vf_alnid',1,d2none,'MODIS-derived diffuse albedo for near-IR radiation  (group VF=VARSF)',20,KGRIB=NGRBALNID)
call set_fvar(metadata,fid%vf_aluvi,'vf_aluvi',1,d2none,'MODIS-derived isotropic albedo coeff for UV/vis radiation (group VF=VARSF)',20,KGRIB=NGRBALUVI)
call set_fvar(metadata,fid%vf_aluvv,'vf_aluvv',1,d2none,'MODIS-derived volumetric albedo coeff for UV/vis radiation (group VF=VARSF)',20,KGRIB=NGRBALUVV)
call set_fvar(metadata,fid%vf_aluvg,'vf_aluvg',1,d2none,'MODIS-derived geometric albedo coeff for UV/vis radiation (group VF=VARSF)',20,KGRIB=NGRBALUVG)
call set_fvar(metadata,fid%vf_alnii,'vf_alnii',1,d2none,'MODIS-derived isotropic albedo coeff for near-IR radiation (group VF=VARSF)',20,KGRIB=NGRBALNII)
call set_fvar(metadata,fid%vf_alniv,'vf_alniv',1,d2none,'MODIS-derived volumetric albedo coeff for near-IR radiation (group VF=VARSF)',20,KGRIB=NGRBALNIV)
call set_fvar(metadata,fid%vf_alnig,'vf_alnig',1,d2none,'MODIS-derived geometric albedo coeff for near-IR radiation (group VF=VARSF)',20,KGRIB=NGRBALNIG)
call set_fvar(metadata,fid%vf_fp1 ,    'vf_fp1',1,d2none,' surface orography in the 2nd part of fullpos-927',20)
call set_fvar(metadata,fid%vf_so2dd ,  'vf_so2dd',1,d2none,' sulphate dry dep velocity',20)
call set_fvar(metadata,fid%vf_dmso ,   'vf_dmso',1,d2none,' oceanic DMS',20)
call set_fvar(metadata,fid%vf_urbf ,  'vf_urbf',1,d2none,' Urban fraction',20)
call set_fvar(metadata,fid%vf_fca1,    'vf_fca1',1,d2none,' fraction of calcite over dust 1st bin',20)
call set_fvar(metadata,fid%vf_fca2,    'vf_fca2',1,d2none,' fraction of calcite over dust 2st bin',20)
call set_fvar(metadata,fid%vf_aerdep,  'vf_aerdep',1,d2none,' dust emission potential',20) 
call set_fvar(metadata,fid%vf_aerlts,  'vf_aerlts',1,d2none,' dust lifting threshold speed ',20)
call set_fvar(metadata,fid%vf_aerscc,  'vf_aerscc',1,d2none,' dust soil clay content',20)
call set_fvar(metadata,fid%vf_dsf,  'vf_dsf',1,d2none,' dust source function',20)
call set_fvar(metadata,fid%vf_dsz,  'vf_dsz',1,d2none,' dust size distribution modulation',20)
call set_fvar(metadata,fid%vf_chemflxo, 'vf_chemflxo',2,d2chemflxo,' total chemistry flux',20)
call set_fvar(metadata,fid%vf_chemwdflx, 'vf_chemwdflx',2,d2chemwdflx,' wet deposition chemistry flux',20)
call set_fvar(metadata,fid%vf_chemddflx, 'vf_chemddflx',2,d2chemddflx,' dry deposition chemistry flux',20)
call set_fvar(metadata,fid%vf_chemdv,  'vf_chemdv',2,d2chemdv,' chemistry deposition velocity',20)
call set_fvar(metadata,fid%vf_nudm,    'vf_nudm ',1,d2none,' nudging mask',20)
call set_fvar(metadata,fid%vf_emis2d,  'vf_emis2d',2,d2emis2d,' 2D emissions for composition',20)
call set_fvar(metadata,fid%vf_emis2daux, 'vf_emis2d',2,d2emis2daux,' 2D emission aux fields for composition',20)
! Group VP
call set_fvar(metadata,fid%vp_tpc,   'vp_tpc',1,d2none,'climatological deep layer temperature',21)
call set_fvar(metadata,fid%vp_wpc,   'vp_wpc',1,d2none,'climatological deep layer moisture',21)

! Group VV=VCLIV: vegetation diagnostic fields
call set_fvar(metadata,fid%vv_arg,  'vv_arg',1,d2none,'Silt percentage within soil (group VV=VCLIV)',22)
call set_fvar(metadata,fid%vv_sab,  'vv_sab',1,d2none,'Percentage of sand within the soil (group VV=VCLIV)',22)
call set_fvar(metadata,fid%vv_d2,   'vv_d2',1,d2none,'soil depth',22)
call set_fvar(metadata,fid%vv_iveg, 'vv_iveg',1,d2none,'type of vegetation',22)
call set_fvar(metadata,fid%vv_rsmin,'vv_rsmin',1,d2none,'stomatal minimum resistance',22)
call set_fvar(metadata,fid%vv_lai,  'vv_lai',1,d2none,'leaf area index',22)
call set_fvar(metadata,fid%vv_hv,   'vv_hv', 1,d2none,'Resistance to evapotranspiration (group VV=VCLIV)',22)
call set_fvar(metadata,fid%vv_z0h,  'vv_z0h',1,d2none,'Gravity times roughness length for heat (group VV=VCLIV)',22)
call set_fvar(metadata,fid%vv_als,  'vv_als',1,d2none,'albedo of bare ground',22)
call set_fvar(metadata,fid%vv_alv,  'vv_alv',1,d2none,'albedo of vegetation',22)
! Group VN
call set_fvar(metadata,fid%vn_top,'vn_top',1,d2none,'index of convective cloud top (group VN=VCLIN)',23)
call set_fvar(metadata,fid%vn_bas,'vn_bas',1,d2none,'index of convective cloud base (group VN=VCLIN)',23)
call set_fvar(metadata,fid%vn_acpr,'vn_acpr',1,d2none,'averaged cumulative precipitation (group VN=VCLIN)',23)
call set_fvar(metadata,fid%vn_accpr,'vn_accpr',1,d2none,'accumulated total precipitation for assimilation (group VN=VCLIN)',23)
call set_fvar(metadata,fid%vn_accpr5,'vn_accpr5',1,d2none,'accumulated total precipitation trajectory for assimilation (group VN=VCLIN)',23)
! Group VH
call set_fvar(metadata,fid%vh_tcch,'vh_tcch',1,d2none,'total convective cloudiness',24)
call set_fvar(metadata,fid%vh_scch,'vh_scch',1,d2none,'convective cloud summit',24)
call set_fvar(metadata,fid%vh_bcch,'vh_bcch',1,d2none,'convective cloud base',24)
call set_fvar(metadata,fid%vh_pblh,'vh_pblh',1,d2none,'PBL height',24)
call set_fvar(metadata,fid%vh_spsh,'vh_spsh',1,d2none,'variable for prognostic convection scheme (ALARO)',24)
call set_fvar(metadata,fid%vh_qsh,'vh_qsh',1,d2none,'surface moisture historic variable (used by TOUCANS)',24)
! Group VK
call set_fvar(metadata,fid%vk_udgro,'vk_udgro',1,d2none,'ud top position (accsu)',35)
! Group VA
call set_fvar(metadata,fid%va_sea,'va_sea',1,d2none,'aerosol: sea',25)
call set_fvar(metadata,fid%va_lan,'va_lan',1,d2none,'aerosol: land',25)
call set_fvar(metadata,fid%va_soo,'va_soo',1,d2none,'aerosol: soot',25)
call set_fvar(metadata,fid%va_des,'va_des',1,d2none,'aerosol: desert',25)
call set_fvar(metadata,fid%va_sul,'va_sul',1,d2none,'aerosol: sulfate',25)
call set_fvar(metadata,fid%va_vol,'va_vol',1,d2none,'aerosol: volcano',25)
! Group VG
call set_fvar(metadata,fid%vg_icfr,'vg_icfr',1,d2none,'sea-ice fraction',26)
call set_fvar(metadata,fid%vg_soup,'vg_soup',1,d2none,'upward solar flux over sea-ice',26)
call set_fvar(metadata,fid%vg_irup,'vg_irup',1,d2none,'upward IR flux over sea-ice',26)
call set_fvar(metadata,fid%vg_chss,'vg_chss',1,d2none,'sensible heat over sea-ice',26)
call set_fvar(metadata,fid%vg_evap,'vg_evap',1,d2none,'evaporation over sea-ice',26)
call set_fvar(metadata,fid%vg_taux,'vg_taux',1,d2none,'U-component of stress over sea-ice',26)
call set_fvar(metadata,fid%vg_tauy,'vg_tauy',1,d2none,'V-component of stress over sea-ice',26)
! Group VC
call set_fvar(metadata,fid%vc_a,'vc_a',1,d2none,'A climatological ozone profile',27)
call set_fvar(metadata,fid%vc_b,'vc_b',1,d2none,'B climatological ozone profile',27)
call set_fvar(metadata,fid%vc_c,'vc_c',1,d2none,'C climatological ozone profile',27)
! Group V2
call set_fvar(metadata,fid%v2_ocdep,'v2_ocdep',1,d2none,'bottom layer depth',28)
call set_fvar(metadata,fid%v2_ustrc,'v2_ustrc',1,d2none,'taux clim.',28)
call set_fvar(metadata,fid%v2_vstrc,'v2_vstrc',1,d2none,'tauy clim.',28)
! Group V3 
call set_fvar(metadata,fid%v3_difm,'v3_difm',1,d2none,'viscosity',29)
call set_fvar(metadata,fid%v3_dift,'v3_dift',1,d2none,'diff. coef. of temp',29)
call set_fvar(metadata,fid%v3_difs,'v3_difs',1,d2none,'diff. coef. of salinity',29)
call set_fvar(metadata,fid%v3_advt,'v3_advt',1,d2none,'correction term for temp.',29)
call set_fvar(metadata,fid%v3_advs,'v3_advs',1,d2none,'correction term for sal.',29)
call set_fvar(metadata,fid%v3_tri0,'v3_tri0',1,d2none,'coef. for solving matrix.',29)
call set_fvar(metadata,fid%v3_tri1,'v3_tri1',1,d2none,'coef. for solving matrix.',29)
call set_fvar(metadata,fid%v3_swdk,'v3_swdk',1,d2none,'radiation term',29)
call set_fvar(metadata,fid%v3_zo,'v3_zo',1,d2none,'depth of layer',29)
call set_fvar(metadata,fid%v3_ho,'v3_ho',1,d2none,'depth of interface layer',29)
call set_fvar(metadata,fid%v3_do,'v3_do',1,d2none,' layer thickness',29)
call set_fvar(metadata,fid%v3_ho_inv,'v3_ho_inv',1,d2none,'1 / YHO',29)
call set_fvar(metadata,fid%v3_uoc,'v3_uoc',1,d2none,'U velocity clim.',29)
call set_fvar(metadata,fid%v3_voc,'v3_voc',1,d2none,'V velocity clim.',29)
call set_fvar(metadata,fid%v3_otke,'v3_otke',1,d2none,'ocean turb. kin. energy',29)

! Group VD
call set_fvar(metadata,fid%vd_lsp,'vd_lsp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_cp,'vd_cp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sf,'vd_sf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_fzra,'vd_fzra',1,d2none,'',30)
call set_fvar(metadata,fid%vd_bld,'vd_bld',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sshf,'vd_sshf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_slhf,'vd_slhf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_nee,'vd_nee',1,d2none,'',30)
call set_fvar(metadata,fid%vd_gpp,'vd_gpp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_rec,'vd_rec',1,d2none,'',30)
call set_fvar(metadata,fid%vd_msl,'vd_msl',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sp,'vd_sp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcc,'vd_tcc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10u,'vd_10u',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10v,'vd_10v',1,d2none,'',30)
call set_fvar(metadata,fid%vd_2t,'vd_2t',1,d2none,'',30)
call set_fvar(metadata,fid%vd_2d,'vd_2d',1,d2none,'',30)
call set_fvar(metadata,fid%vd_2sh,'vd_2sh',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ssr,'vd_ssr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_str,'vd_str',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tsr,'vd_tsr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ttr,'vd_ttr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ewss,'vd_ewss',1,d2none,'',30)
call set_fvar(metadata,fid%vd_nsss,'vd_nsss',1,d2none,'',30)
call set_fvar(metadata,fid%vd_e,'vd_e',1,d2none,'',30)
call set_fvar(metadata,fid%vd_pev,'vd_pev',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ccc,'vd_ccc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lcc,'vd_lcc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mcc,'vd_mcc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_hcc,'vd_hcc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lgws,'vd_lgws',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mgws,'vd_mgws',1,d2none,'',30)
call set_fvar(metadata,fid%vd_gwd,'vd_gwd',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mx2t,'vd_mx2t',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mn2t,'vd_mn2t',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mx2t6,'vd_mx2t6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mn2t6,'vd_mn2t6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ro,'vd_ro',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sro,'vd_sro',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ssro,'vd_ssro',1,d2none,'',30)
call set_fvar(metadata,fid%vd_alb,'vd_alb',1,d2none,'',30)
call set_fvar(metadata,fid%vd_iewss,'vd_iewss',1,d2none,'',30)
call set_fvar(metadata,fid%vd_insss,'vd_insss',1,d2none,'',30)
call set_fvar(metadata,fid%vd_isshf,'vd_isshf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ie,'vd_ie',1,d2none,'',30)
call set_fvar(metadata,fid%vd_inee,'vd_inee',1,d2none,'',30)
call set_fvar(metadata,fid%vd_igpp,'vd_igpp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_irec,'vd_irec',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ich4,'vd_ich4',1,d2none,'',30)
call set_fvar(metadata,fid%vd_csf,'vd_csf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lssf,'vd_lssf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mxtpr,'vd_mxtpr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mntpr,'vd_mntpr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mxtpr6,'vd_mxtpr6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mntpr6,'vd_mntpr6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tpr,'vd_tpr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lsrr,'vd_lsrr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_crr,'vd_crr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lssfr,'vd_lssfr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_csfr,'vd_csfr',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ptype,'vd_ptype',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ilspf,'vd_ilspf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_z0f,'vd_z0f',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lz0h,'vd_lz0h',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcw,'vd_tcw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcwv,'vd_tcwv',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tclw,'vd_tclw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tciw,'vd_tciw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcrw,'vd_tcrw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcsw,'vd_tcsw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcslw,'vd_tcslw',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ssrd,'vd_ssrd',1,d2none,'',30)
call set_fvar(metadata,fid%vd_strd,'vd_strd',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ssrdc,'vd_ssrdc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_strdc,'vd_strdc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_blh,'vd_blh',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sund,'vd_sund',1,d2none,'',30)
call set_fvar(metadata,fid%vd_spar,'vd_spar',1,d2none,'',30)
call set_fvar(metadata,fid%vd_suvb,'vd_suvb',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sfdir,'vd_sfdir',1,d2none,'',30)
call set_fvar(metadata,fid%vd_scdir,'vd_scdir',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sdsrp,'vd_sdsrp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_cape,'vd_cape',1,d2none,'',30)
call set_fvar(metadata,fid%vd_capes,'vd_capes',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mucape,'vd_mucape',1,d2none,'',30)
call set_fvar(metadata,fid%vd_pdepl,'vd_pdepl',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mlcape50,'vd_mlcape50',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mlcape100,'vd_mlcape100',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mlcin50,'vd_mlcin50',1,d2none,'',30)
call set_fvar(metadata,fid%vd_mlcin100,'vd_mlcin100',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tropotp,'vd_tropotp',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tsrc,'vd_tsrc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ttrc,'vd_ttrc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ssrc,'vd_ssrc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_strc,'vd_strc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_es,'vd_es',1,d2none,'',30)
call set_fvar(metadata,fid%vd_smlt,'vd_smlt',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10fg,'vd_10fg',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10fg6,'vd_10fg6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10fgcv,'vd_10fgcv',1,d2none,'',30)
call set_fvar(metadata,fid%vd_i10fg,'vd_i10fg',1,d2none,'',30)
call set_fvar(metadata,fid%vd_lspf,'vd_lspf',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tco3,'vd_tco3',1,d2none,'',30)
call set_fvar(metadata,fid%vd_vimd,'vd_vimd',1,d2none,'',30)
call set_fvar(metadata,fid%vd_sparc,'vd_sparc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_stinc,'vd_stinc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_cbase,'vd_cbase',1,d2none,'',30)
call set_fvar(metadata,fid%vd_0degl,'vd_0degl',1,d2none,'',30)
call set_fvar(metadata,fid%vd_m10degl,'vd_m10degl',1,d2none,'',30)
call set_fvar(metadata,fid%vd_visih,'vd_visih',1,d2none,'',30)
call set_fvar(metadata,fid%vd_cin,'vd_cin',1,d2none,'',30)
call set_fvar(metadata,fid%vd_kindex,'vd_kindex',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ttindex,'vd_ttindex',1,d2none,'',30)
call set_fvar(metadata,fid%vd_cbasea,'vd_cbasea',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ctopc,'vd_ctopc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ztwetb0,'vd_ztwetb0',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ztwetb1,'vd_ztwetb1',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcghg,'vd_tcghg',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tcchem,'vd_tcchem',1,d2none,'',30)
call set_fvar(metadata,fid%vd_aerodiag,'vd_aerodiag',1,d2none,'',30)
call set_fvar(metadata,fid%vd_aero_wvl_diag,'vd_aero_vwl_diag',1,d2none,'',30)
call set_fvar(metadata,fid%vd_100u,'vd_100u',1,d2none,'',30)
call set_fvar(metadata,fid%vd_100v,'vd_100v',1,d2none,'',30)
call set_fvar(metadata,fid%vd_zust,'vd_zust',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10nu,'vd_10nu',1,d2none,'',30)
call set_fvar(metadata,fid%vd_10nv,'vd_10nv',1,d2none,'',30)
call set_fvar(metadata,fid%vd_dndzn,'vd_dndzn',1,d2none,'',30)
call set_fvar(metadata,fid%vd_dndza,'vd_dndza',1,d2none,'',30)
call set_fvar(metadata,fid%vd_dctb,'vd_dctb',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tplb,'vd_tplb',1,d2none,'',30)
call set_fvar(metadata,fid%vd_tplt,'vd_tplt',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odss,'vd_odss',1,d2none,'',30)
call set_fvar(metadata,fid%vd_oddu,'vd_oddu',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odom,'vd_odom',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odbc,'vd_odbc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odsu,'vd_odsu',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odni,'vd_odni',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odam,'vd_odam',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odsoa,'vd_odsoa',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odvfa,'vd_odvfa',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odvsu,'vd_odvsu',1,d2none,'',30)
call set_fvar(metadata,fid%vd_odtoacc,'vd_odtoacc',1,d2none,'',30)
call set_fvar(metadata,fid%vd_aepm1,'vd_aepm1',1,d2none,'',30)
call set_fvar(metadata,fid%vd_aepm25,'vd_aepm25',1,d2none,'',30)
call set_fvar(metadata,fid%vd_aepm10,'vd_aepm10',1,d2none,'',30)
call set_fvar(metadata,fid%vd_uvbed,'vd_uvbed',1,d2none,'',30)
call set_fvar(metadata,fid%vd_uvbedcs,'vd_uvbedcs',1,d2none,'',30)
call set_fvar(metadata,fid%vd_litoti,'vd_litoti',1,d2none,'',30)
call set_fvar(metadata,fid%vd_licgi,'vd_licgi',1,d2none,'',30)
call set_fvar(metadata,fid%vd_litota6,'vd_litota6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_licga6,'vd_licga6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_ptypeocc6,'vd_ptypeocc6',1,d2none,'',30)
call set_fvar(metadata,fid%vd_200u,'vd_200u',1,d2none,'',30)
call set_fvar(metadata,fid%vd_200v,'vd_200v',1,d2none,'',30)

call set_fvar(metadata,fid%vd_sdsl,'vd_sdsl',1,d2none,'',30)

! Group SM
call set_fvar(metadata,fid%sm_clbt,'sm_clbt',1,d2none,'Cloudy brightness temperature',32)
call set_fvar(metadata,fid%sm_csbt,'sm_csbt',1,d2none,'Clear-sky brightness temperature',32)


! Group WS
call set_fvar(metadata,fid%ws_char,'ws_char',1,d2none,'Charnock parameter as modified by the wave model ',33,KGRIB=NGRBCHNK)
call set_fvar(metadata,fid%ws_charhq,'ws_charhq',1,d2none,'Charnock for heat and moisture from the wave model ',33)
call set_fvar(metadata,fid%ws_ustokes,'ws_ustokes',1,d2none,'U-component of the surface Stokes drift',33)
call set_fvar(metadata,fid%ws_vstokes,'ws_vstokes',1,d2none,'V-component of the surface Stokes drift',33)
call set_fvar(metadata,fid%ws_tauocx,'ws_tauocx',1,d2none,'U-component of the Momentum flux to ocean',33)
call set_fvar(metadata,fid%ws_tauocy,'ws_tauocy',1,d2none,'V-component of the Momentum flux to ocean',33)
call set_fvar(metadata,fid%ws_phioc,'ws_phioc',1,d2none,'Energy flux to ocean',33)
call set_fvar(metadata,fid%ws_wsemean,'ws_wsemean',1,d2none,'Windsea variance',33)
call set_fvar(metadata,fid%ws_wsfmean,'ws_wsfmean',1,d2none,'Windsea mean frequency',33)

! Group WW
call set_fvar(metadata,fid%ww_u10n,'ww_u10n',1,d2none,'10m neutral wind U-component passed to the wave model (WAM) ',38)
call set_fvar(metadata,fid%ww_v10n,'ww_v10n',1,d2none,'10m neutral wind V-component passed to the wave model (WAM)',38)
call set_fvar(metadata,fid%ww_rho,'ww_rho',1,d2none,'surface density passed to the wave model (WAM)',38)
call set_fvar(metadata,fid%ww_zil,'ww_zil',1,d2none,'ZI/L passed to the wave model (used for gustiness in WAM)',38)
call set_fvar(metadata,fid%ww_cif,'ww_cif',1,d2none,'Sea ice fraction passed to the wave model (WAM)',38)
call set_fvar(metadata,fid%ww_clk,'ww_clk',1,d2none,'Lake cover passed to the wave model (WAM)',38)
call set_fvar(metadata,fid%ww_ucurw,'ww_ucurw',1,d2none,'Ocean current    U-component passed to the wave model (WAM)',38)
call set_fvar(metadata,fid%ww_vcurw,'ww_vcurw',1,d2none,'Ocean current    V-component passed to the wave model (WAM)',38)

! Group VX
call set_fvar(metadata,fid%vx_oro,'vx_oro',1,d2none,'Climatological surface geopotential (group VX=VCLIX)',34)
call set_fvar(metadata,fid%vx_tsc,'vx_tsc',1,d2none,'climatological surface temperature',34)
call set_fvar(metadata,fid%vx_pws,'vx_pws',1,d2none,'climatological surface max. prop. moisture',34)
call set_fvar(metadata,fid%vx_pwp,'vx_pwp',1,d2none,'climatological deep soil max. prop. moisture',34)
call set_fvar(metadata,fid%vx_sno,'vx_sno',1,d2none,'Climatological snow cover (group VX=VCLIX)',34)
call set_fvar(metadata,fid%vx_tpc,'vx_tpc',1,d2none,'climatological deep soil temperature',34)
call set_fvar(metadata,fid%vx_sab,'vx_sab',1,d2none,'climatologic percentage of sand within the soil',34)
call set_fvar(metadata,fid%vx_xd2,'vx_xd2',1,d2none,'climatologic soil depth',34)
call set_fvar(metadata,fid%vx_lsm,'vx_lsm',1,d2none,'climatologic land sea mask',34)
call set_fvar(metadata,fid%vx_iveg,'vx_iveg',1,d2none,'climatologic type of vegetation',34)
call set_fvar(metadata,fid%vx_arg,'vx_arg',1,d2none,'silt percentage within soil',34)
call set_fvar(metadata,fid%vx_rsmin,'vx_rsmin',1,d2none,'climatologic stomatal minimum resistance',34)
call set_fvar(metadata,fid%vx_lai,'vx_lai',1,d2none,'leaf area index',34)
call set_fvar(metadata,fid%vx_veg,'vx_veg',1,d2none,'vegetation cover',34)


! Dynamically-assigned COMPO fields in GFL
do jname=1,ndyn_names !                    compo_3d_first,compo_3d_first+ndyn_names-1
  jfid = jname+compo_3d_first-1
  write(nulout,*) 'dynamic set_fvar ',jfid,dyn_nam(jname)%cname
  call set_fvar(metadata,jfid,dyn_nam(jname)%cname,2,d2full,dyn_nam(jname)%clongname,3,kgrib=dyn_nam(jname)%igribcode)
enddo
do jfid=compo_3d_first+ndyn_names, compo_3d_last
  write (clname, fmt='("compo_3d_",I5.5)') jfid-compo_3d_first+1
  call set_fvar(metadata,jfid,clname,2,d2full,'Atomospheric composition 3D field',3)
enddo


! Canari 2D and 3D model errors
call set_fvar(metadata,fid%est,'est',1,d2none,'2D model errors statistics',0)
call set_fvar(metadata,fid%esn,'esn',1,d2none,'2D model errors statistics',0)
call set_fvar(metadata,fid%et2,'et2',1,d2none,'2D model errors statistics',0)
call set_fvar(metadata,fid%eh2,'eh2',1,d2none,'2D model errors statistics',0)
call set_fvar(metadata,fid%ev1,'ev1',1,d2none,'2D model errors statistics',0)
call set_fvar(metadata,fid%ez, 'ez', 2,d2full,'3D model errors statistics',0)
call set_fvar(metadata,fid%eh, 'eh', 2,d2full,'3D model errors statistics',0)
! SLB2 fields (whatever they are - NFTBC AJG)
call set_fvar(metadata,fid%vvel,'vvel',2,d2full,'Vertical velocity',0)

IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:MAIN_FIELD_SET_METADATA',1,ZHOOK_HANDLE)
end subroutine main_field_set_metadata

function main_field_get_clevtype(this, kleveltype) result(clevtype)
class(field_metadata), intent(in) :: this
integer(kind=jpim), intent(in) :: kleveltype
character(len=3) :: clevtype

select case (kleveltype)
  case(d2full, d2half)
    clevtype = "ML"
  case(d2none, d2soil, d2snow, d2om, d2chemflxo, d2chemwdflx, d2chemddflx, d2chemdv, d2emis2d, d2emis2daux)
    clevtype = "SFC"
  case default
    call abor1("field_definitions:main_field_get_clevtype: unknown level type")
end select

end function main_field_get_clevtype

! -------------------------------------------------
!
! Perhaps clunky. This kind of technique was not
! used for the GOMs, because it could be quite slow 
! in that context. It should be more appropriate 
! here when the number of calls should not be too 
! large.
! -------------------------------------------------
subroutine main_field_map_storage(this, kid, storage_1d, storage_2d, storage_3d, ld_nullify)

class(field_access),               intent(inout) :: this
integer(kind=jpim),                intent(in)    :: kid
real(kind=jprb), optional, target, intent(in)    :: storage_1d(:), storage_2d(:,:),storage_3d(:,:,:)
logical, optional,                 intent(in)    :: ld_nullify 

logical :: ll_nullify
real(kind=jprb), pointer :: z1d(:), z2d(:,:),z3d(:,:,:)
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:MAIN_FIELD_MAP_STORAGE',0,ZHOOK_HANDLE)
ll_nullify = .false.
if (present(ld_nullify)) ll_nullify = ld_nullify

if(ll_nullify) then
  z1d=>null()
  z2d=>null()
  z3d=>null()
else
  if(present(storage_1d)) z1d => storage_1d
  if(present(storage_2d)) z2d => storage_2d
  if(present(storage_3d)) z3d => storage_3d
endif

! Field name mapping (F7)
select case(kid)

! GMV
case(fid%u)
  this%u => z2d 
case(fid%v)
  this%v => z2d 
case(fid%div)
  this%div => z2d
case(fid%vor)
  this%vor => z2d
case(fid%t)
  this%t => z2d
case(fid%pd)
  this%pd => z2d 
case(fid%vd)
  this%vd => z2d 
case(fid%nhx)
  this%nhx => z2d
case(fid%edot)
  this%edot => z2d
case(fid%vcasrs)      ! ky: obsolete field 'vcasrs' (deep-layer models pruned)
  this%vcasrs => z2d
case(fid%rdlasr)      ! ky: obsolete field 'rdlasr' (deep-layer models pruned)
  this%rdlasr => z2d
case(fid%sgrtl)
  this%sgrtl => z2d
case(fid%sgrtm)
  this%sgrtm => z2d
case(fid%unl)
  this%unl => z2d
case(fid%unl_si)
  this%unl_si => z2d
case(fid%vnl)
  this%vnl => z2d
case(fid%vnl_si)
  this%vnl_si => z2d
case(fid%tnl)
  this%tnl => z2d
case(fid%tnl_si)
  this%tnl_si => z2d
case(fid%spnl)
  this%spnl => z2d
case(fid%spnl_si)
  this%spnl_si => z2d
case(fid%vwvnl)
  this%vwvnl => z2d
case(fid%vdnl_si)
  this%vdnl_si => z2d
case(fid%pdnl)
  this%pdnl => z2d
case(fid%pdnl_si)
  this%pdnl_si => z2d
case(fid%gw)
  this%gw => z2d
case(fid%curhs)
  this%curhs => z2d
case(fid%cvrhs)
  this%cvrhs => z2d
case(fid%ctrhs)
  this%ctrhs => z2d
case(fid%cspdrhs)
  this%cspdrhs => z2d
case(fid%csvdrhs)
  this%csvdrhs => z2d
case(fid%cunl)
  this%cunl => z2d
case(fid%cvnl)
  this%cvnl => z2d
case(fid%ctnl)
  this%ctnl => z2d
case(fid%cspnl)
  this%cspnl => z2d
case(fid%cvwvnl)
  this%cvwvnl => z2d
case(fid%cspdnl)
  this%cspdnl => z2d
case(fid%cupt)
  this%cupt => z2d
case(fid%cvpt)
  this%cvpt => z2d
case(fid%ctpt)
  this%ctpt=> z2d
case(fid%cpdpt)
  this%cpdpt => z2d
case(fid%cvdpt)
  this%cvdpt => z2d
case(fid%dphi)
  this%dphi => z2d
case(fid%nhxnl)
  this%nhxnl => z2d
case(fid%cnhxnl)
  this%cnhxnl => z2d

! GFL      
case(fid%q)
  this%q  => z2d 
case(fid%o3)
  this%o3 => z2d 
case(fid%l)
  this%l  => z2d
case(fid%i)
  this%i  => z2d
case(fid%a)
  this%a  => z2d
case(fid%s)
  this%s  => z2d
case(fid%r)
  this%r  => z2d
case(fid%g)
  this%g  => z2d
case(fid%h)
  this%h  => z2d
case(fid%lrch4)
  this%lrch4 => z2d
case(fid%lconv)
  this%lconv  => z2d
case(fid%iconv)
  this%iconv  => z2d
case(fid%rconv)
  this%rconv  => z2d
case(fid%sconv)
  this%sconv  => z2d
case(fid%tke)
  this%tke  => z2d
case(fid%tte)
  this%tte  => z2d
case(fid%efb1)
  this%efb1  => z2d
case(fid%efb2)
  this%efb2  => z2d
case(fid%efb3)
  this%efb3  => z2d
case(fid%mxl)
  this%mxl  => z2d
case(fid%src)
  this%src  => z2d
case(fid%shtur)
  this%shtur  => z2d
case(fid%fqtur)
  this%fqtur  => z2d
case(fid%fstur)
  this%fstur  => z2d
case(fid%cvv)
  this%cvv  => z2d
case(fid%rkth)
  this%rkth  => z2d
case(fid%rktqv)
  this%rktqv  => z2d
case(fid%rktqc)
  this%rktqc  => z2d
case(fid%cpf)
  this%cpf  => z2d
case(fid%spf)
  this%spf  => z2d
case(fid%cvgq)
  this%cvgq  => z2d
case(fid%uom)
  this%uom  => z2d
case(fid%ual)
  this%ual  => z2d
case(fid%dom)
  this%dom  => z2d
case(fid%dal)
  this%dal  => z2d
case(fid%uen)
  this%uen  => z2d
case(fid%unebh)
  this%unebh  => z2d
case(fid%lrad)
  this%lrad  => z2d
case(fid%irad)
  this%irad  => z2d
case(fid%rspec)
  this%rspec  => z2d
case(fid%extra3d)
  this%extra3d  => z3d
case(fid%ezdiag)
  this%ezdiag  => z3d
case(fid%unogw)
  this%unogw  => z2d
case(fid%vnogw)
  this%vnogw  => z2d

! Pressure
case(fid%pre)
  this%pre  => z2d
case(fid%pref)
  this%pref  => z2d
case(fid%delp)
  this%delp  => z2d

! GMVS
case(fid%sp)
  this%sp  => z1d
case(fid%csppt)
  this%csppt  => z1d
case(fid%dbbc)
  this%dbbc  => z1d
case(fid%gws)
  this%gws  => z1d
case(fid%csprhs)
  this%csprhs => z1d
case(fid%spnl2)
  this%spnl2 => z1d
case(fid%cspnl2)
  this%cspnl2 => z1d
case(fid%prehyds)
  this%prehyds  => z1d

! mean winds for lelem
case(fid%u_mean)
  this%u_mean => z1d
case(fid%v_mean)
  this%v_mean => z1d

! Constants
case(fid%orog)
  this%orog  => z1d
case(fid%cori)
  this%cori  => z1d
case(fid%gnordl)
  this%gnordl  => z1d
case(fid%gnordm)
  this%gnordm  => z1d

!Cloud
case(fid%ccc)
  this%ccc  => z1d
case(fid%lcc)
  this%lcc  => z1d
case(fid%mcc)
  this%mcc  => z1d
case(fid%hcc)
  this%hcc  => z1d
case(fid%tcc)
  this%tcc  => z1d

! Surface

! Group SB
case(fid%sb_t)
  this%sb_t  => z2d
case(fid%sb_q)
  this%sb_q  => z2d
case(fid%sb_tl)
  this%sb_tl  => z2d

! Group SG
case(fid%sg_f)
  this%sg_f  => z2d
case(fid%sg_a)
  this%sg_a  => z2d
case(fid%sg_r)
  this%sg_r  => z2d
case(fid%sg_t)
  this%sg_t  => z2d
case(fid%sg_w)
  this%sg_w  => z2d

! Group SL
case(fid%sl_lict)
  this%sl_lict  => z1d
case(fid%sl_lmlt)
  this%sl_lmlt  => z1d
case(fid%sl_ltlt)
  this%sl_ltlt  => z1d
case(fid%sl_lblt)
  this%sl_lblt  => z1d
case(fid%sl_lshf)
  this%sl_lshf  => z1d
case(fid%sl_licd)
  this%sl_licd  => z1d
case(fid%sl_lmld)
  this%sl_lmld  => z1d

! Group RR
case(fid%rr_t)
  this%rr_t  => z1d
case(fid%rr_w)
  this%rr_w  => z1d
case(fid%rr_fc)
  this%rr_fc => z1d
case(fid%rr_ic)
  this%rr_ic => z1d
case(fid%rr_fp1)
  this%rr_fp1 => z1d

! Group CL
case(fid%cl_tcls)
  this%cl_tcls  => z1d
case(fid%cl_hucls)
  this%cl_hucls  => z1d

! Group OM
case(fid%om_to)
  this%om_to  => z2d
case(fid%om_so)
  this%om_so  => z2d
case(fid%om_uo)
  this%om_uo  => z2d
case(fid%om_vo)
  this%om_vo  => z2d

! Group EP
case(fid%ep_ep)
  this%ep_ep  => z2d

! Group X2
case(fid%x2_x2)
  this%x2_x2  => z2d

! Group CI
case(fid%ci_ci)
  this%ci_ci  => z2d

! Group VF
case(fid%vf_z0f)
  this%vf_z0f  => z1d
case(fid%vf_albf)
  this%vf_albf  => z1d
case(fid%vf_emisf)
  this%vf_emisf  => z1d
case(fid%vf_getrl)
  this%vf_getrl  => z1d
case(fid%vf_lsm)
  this%vf_lsm  => z1d
case(fid%vf_veg)
  this%vf_veg  => z1d
case(fid%vf_vrlan)
  this%vf_vrlan  => z1d
case(fid%vf_vrldi)
  this%vf_vrldi  => z1d
case(fid%vf_sig)
  this%vf_sig  => z1d
case(fid%vf_albsf)
  this%vf_albsf  => z1d
case(fid%vf_lan)
  this%vf_lan  => z1d
case(fid%vf_sst)
  this%vf_sst  => z1d
case(fid%vf_sss)
  this%vf_sss => z1d
case(fid%vf_lz0h)
  this%vf_lz0h  => z1d
case(fid%vf_cvl)
  this%vf_cvl  => z1d
case(fid%vf_co2typ)
  this%vf_co2typ  => z1d
case(fid%vf_cvh)
  this%vf_cvh  => z1d
case(fid%vf_fwet)
  this%vf_fwet  => z1d
case(fid%vf_cur)
  this%vf_cur  => z1d
case(fid%vf_tvl)
  this%vf_tvl  => z1d
case(fid%vf_tvh)
  this%vf_tvh  => z1d
case(fid%vf_lail)
  this%vf_lail  => z1d
case(fid%vf_laih)
  this%vf_laih  => z1d
case(fid%vf_soty)
  this%vf_soty  => z1d
case(fid%vf_clk)
  this%vf_clk  => z1d
case(fid%vf_dl)
  this%vf_dl  => z1d
case(fid%vf_ci)
  this%vf_ci  => z1d
case(fid%vf_ucur)
  this%vf_ucur  => z1d
case(fid%vf_vcur)
  this%vf_vcur  => z1d
case(fid%vf_z0rlf)
  this%vf_z0rlf  => z1d
case(fid%vf_cgpp)
  this%vf_cgpp  => z1d
case(fid%vf_crec)
  this%vf_crec  => z1d
case(fid%vf_sdfor)
  this%vf_sdfor  => z1d
case(fid%vf_aluvp)
  this%vf_aluvp  => z1d
case(fid%vf_aluvd)
  this%vf_aluvd  => z1d
case(fid%vf_alnip)
  this%vf_alnip  => z1d
case(fid%vf_alnid)
  this%vf_alnid  => z1d
case(fid%vf_aluvi)
  this%vf_aluvi  => z1d
case(fid%vf_aluvv)
  this%vf_aluvv  => z1d
case(fid%vf_aluvg)
  this%vf_aluvg  => z1d
case(fid%vf_alnii)
  this%vf_alnii  => z1d
case(fid%vf_alniv)
  this%vf_alniv  => z1d
case(fid%vf_alnig)
  this%vf_alnig  => z1d
case(fid%vf_fp1)
  this%vf_fp1  => z1d
case(fid%vf_so2dd)
  this%vf_so2dd  => z1d
case(fid%vf_dmso)
  this%vf_dmso  => z1d
case(fid%vf_urbf)
  this%vf_urbf  => z1d
case(fid%vf_fca1)
  this%vf_fca1  => z1d
case(fid%vf_fca2)
  this%vf_fca2  => z1d
case(fid%vf_aerdep)
  this%vf_aerdep  => z1d
case(fid%vf_aerlts)
  this%vf_aerlts  => z1d
case(fid%vf_aerscc)
  this%vf_aerscc  => z1d
case(fid%vf_dsf)
  this%vf_dsf  => z1d
case(fid%vf_dsz)
  this%vf_dsz  => z1d
case(fid%vf_chemflxo)
  this%vf_chemflxo  => z2d
case(fid%vf_chemwdflx)
  this%vf_chemwdflx  => z2d
case(fid%vf_chemddflx)
  this%vf_chemddflx  => z2d
case(fid%vf_chemdv)
  this%vf_chemdv  => z2d
case(fid%vf_nudm)
  this%vf_nudm  => z1d
case(fid%vf_emis2d)
  this%vf_emis2d  => z2d
case(fid%vf_emis2daux)
  this%vf_emis2daux  => z2d

! Group VP
case(fid%vp_tpc)
  this%vp_tpc  => z1d
case(fid%vp_wpc)
  this%vp_wpc  => z1d


! Group VV
case(fid%vv_arg)
  this%vv_arg  => z1d
case(fid%vv_sab)
  this%vv_sab  => z1d
case(fid%vv_d2)
  this%vv_d2  => z1d
case(fid%vv_iveg)
  this%vv_iveg  => z1d
case(fid%vv_rsmin)
  this%vv_rsmin  => z1d
case(fid%vv_lai)
  this%vv_lai  => z1d
case(fid%vv_hv)
  this%vv_hv  => z1d
case(fid%vv_z0h)
  this%vv_z0h  => z1d
case(fid%vv_als)
  this%vv_als  => z1d
case(fid%vv_alv)
  this%vv_alv  => z1d

! Group VN
case(fid%vn_top)
  this%vn_top  => z1d
case(fid%vn_bas)
  this%vn_bas  => z1d
case(fid%vn_acpr)
  this%vn_acpr  => z1d
case(fid%vn_accpr)
  this%vn_accpr  => z1d
case(fid%vn_accpr5)
  this%vn_accpr5  => z1d
! Group VH
case(fid%vh_tcch)
  this%vh_tcch  => z1d
case(fid%vh_scch)
  this%vh_scch  => z1d
case(fid%vh_bcch)
  this%vh_bcch  => z1d
case(fid%vh_pblh)
  this%vh_pblh  => z1d
case(fid%vh_spsh)
  this%vh_spsh  => z1d
case(fid%vh_qsh)
  this%vh_qsh  => z1d
! Group VK
case(fid%vk_udgro)
  this%vk_udgro  => z1d
! Group VA
case(fid%va_sea)
  this%va_sea  => z1d
case(fid%va_lan)
  this%va_lan  => z1d
case(fid%va_soo)
  this%va_soo  => z1d
case(fid%va_des)
  this%va_des  => z1d
case(fid%va_sul)
  this%va_sul  => z1d
case(fid%va_vol)
  this%va_vol  => z1d

! Group VG
case(fid%vg_icfr)
  this%vg_icfr  => z1d
case(fid%vg_soup)
  this%vg_soup  => z1d
case(fid%vg_irup)
  this%vg_irup  => z1d
case(fid%vg_chss)
  this%vg_chss  => z1d
case(fid%vg_evap)
  this%vg_evap  => z1d
case(fid%vg_taux)
  this%vg_taux  => z1d
case(fid%vg_tauy)
  this%vg_tauy  => z1d

! Group VC
case(fid%vc_a)
  this%vc_a  => z1d
case(fid%vc_b)
  this%vc_b  => z1d
case(fid%vc_c)
  this%vc_c  => z1d

! Group V2
case(fid%v2_ocdep)
  this%v2_ocdep  => z1d
case(fid%v2_ustrc)
  this%v2_ustrc  => z1d
case(fid%v2_vstrc)
  this%v2_vstrc  => z1d

! Group V3
case(fid%v3_difm)
  this%v3_difm  => z2d
case(fid%v3_dift)
  this%v3_dift  => z2d
case(fid%v3_difs)
  this%v3_difs  => z2d
case(fid%v3_advt)
  this%v3_advt  => z2d
case(fid%v3_advs)
  this%v3_advs  => z2d
case(fid%v3_tri0)
  this%v3_tri0  => z2d
case(fid%v3_tri1)
  this%v3_tri1  => z2d
case(fid%v3_swdk)
  this%v3_swdk  => z2d
case(fid%v3_zo)
  this%v3_zo  => z2d
case(fid%v3_ho)
  this%v3_ho  => z2d
case(fid%v3_do)
  this%v3_do  => z2d
case(fid%v3_ho_inv)
  this%v3_ho_inv  => z2d
case(fid%v3_uoc)
  this%v3_uoc  => z2d
case(fid%v3_voc)
  this%v3_voc  => z2d
case(fid%v3_otke)
  this%v3_otke  => z2d

! Group WS
case(fid%ws_char)
  this%ws_char  => z1d
case(fid%ws_charhq)
  this%ws_charhq  => z1d
case(fid%ws_ustokes)
  this%ws_ustokes  => z1d
case(fid%ws_vstokes)
  this%ws_vstokes  => z1d
case(fid%ws_tauocx)
  this%ws_tauocx  => z1d
case(fid%ws_tauocy)
  this%ws_tauocy  => z1d
case(fid%ws_phioc)
  this%ws_phioc  => z1d
case(fid%ws_wsemean)
  this%ws_wsemean  => z1d
case(fid%ws_wsfmean)
  this%ws_wsfmean  => z1d

! Group WW
case(fid%ww_u10n)
  this%ww_u10n  => z1d
case(fid%ww_v10n)
  this%ww_v10n  => z1d
case(fid%ww_rho)
  this%ww_rho  => z1d
case(fid%ww_zil)
  this%ww_zil  => z1d
case(fid%ww_cif)
  this%ww_cif  => z1d
case(fid%ww_clk)
  this%ww_clk  => z1d
case(fid%ww_ucurw)
  this%ww_ucurw  => z1d
case(fid%ww_vcurw)
  this%ww_vcurw  => z1d

! Group VD
case(fid%vd_lsp)
  this%vd_lsp  => z1d
case(fid%vd_cp)
  this%vd_cp  => z1d
case(fid%vd_sf)
  this%vd_sf  => z1d
case(fid%vd_fzra)
  this%vd_fzra  => z1d
case(fid%vd_bld)
  this%vd_bld  => z1d
case(fid%vd_sshf)
  this%vd_sshf  => z1d
case(fid%vd_slhf)
  this%vd_slhf  => z1d
case(fid%vd_nee)
  this%vd_nee  => z1d
case(fid%vd_gpp)
  this%vd_gpp  => z1d
case(fid%vd_rec)
  this%vd_rec  => z1d
case(fid%vd_msl)
  this%vd_msl  => z1d
case(fid%vd_sp)
  this%vd_sp  => z1d
case(fid%vd_tcc)
  this%vd_tcc  => z1d
case(fid%vd_10u)
  this%vd_10u  => z1d
case(fid%vd_10v)
  this%vd_10v  => z1d
case(fid%vd_2t)
  this%vd_2t  => z1d
case(fid%vd_2d)
  this%vd_2d  => z1d
case(fid%vd_2sh)
  this%vd_2sh => z1d
case(fid%vd_ssr)
  this%vd_ssr  => z1d
case(fid%vd_str)
  this%vd_str  => z1d
case(fid%vd_tsr)
  this%vd_tsr  => z1d
case(fid%vd_ttr)
  this%vd_ttr  => z1d
case(fid%vd_ewss)
  this%vd_ewss  => z1d
case(fid%vd_nsss)
  this%vd_nsss  => z1d
case(fid%vd_e)
  this%vd_e  => z1d
case(fid%vd_pev)
  this%vd_pev  => z1d
case(fid%vd_ccc)
  this%vd_ccc  => z1d
case(fid%vd_lcc)
  this%vd_lcc  => z1d
case(fid%vd_mcc)
  this%vd_mcc  => z1d
case(fid%vd_hcc)
  this%vd_hcc  => z1d
case(fid%vd_lgws)
  this%vd_lgws  => z1d
case(fid%vd_mgws)
  this%vd_mgws  => z1d
case(fid%vd_gwd)
  this%vd_gwd  => z1d
case(fid%vd_mx2t)
  this%vd_mx2t  => z1d
case(fid%vd_mn2t)
  this%vd_mn2t  => z1d
case(fid%vd_mx2t6)
  this%vd_mx2t6  => z1d
case(fid%vd_mn2t6)
  this%vd_mn2t6  => z1d
case(fid%vd_ro)
  this%vd_ro  => z1d
case(fid%vd_sro)
  this%vd_sro  => z1d
case(fid%vd_ssro)
  this%vd_ssro  => z1d
case(fid%vd_alb)
  this%vd_alb  => z1d
case(fid%vd_iewss)
  this%vd_iewss  => z1d
case(fid%vd_insss)
  this%vd_insss  => z1d
case(fid%vd_isshf)
  this%vd_isshf  => z1d
case(fid%vd_ie)
  this%vd_ie  => z1d
case(fid%vd_inee)
  this%vd_inee  => z1d
case(fid%vd_igpp)
  this%vd_igpp  => z1d
case(fid%vd_irec)
  this%vd_irec  => z1d
case(fid%vd_ich4)
  this%vd_ich4  => z1d
case(fid%vd_csf)
  this%vd_csf  => z1d
case(fid%vd_lssf)
  this%vd_lssf  => z1d
case(fid%vd_mxtpr)
  this%vd_mxtpr  => z1d
case(fid%vd_mntpr)
  this%vd_mntpr  => z1d
case(fid%vd_mxtpr6)
  this%vd_mxtpr6  => z1d
case(fid%vd_mntpr6)
  this%vd_mntpr6  => z1d
case(fid%vd_tpr)
  this%vd_tpr  => z1d
case(fid%vd_lsrr)
  this%vd_lsrr  => z1d
case(fid%vd_crr)
  this%vd_crr  => z1d
case(fid%vd_lssfr)
  this%vd_lssfr  => z1d
case(fid%vd_csfr)
  this%vd_csfr  => z1d
case(fid%vd_ptype)
  this%vd_ptype  => z1d
case(fid%vd_ilspf)
  this%vd_ilspf  => z1d
case(fid%vd_z0f)
  this%vd_z0f  => z1d
case(fid%vd_lz0h)
  this%vd_lz0h  => z1d
case(fid%vd_tcw)
  this%vd_tcw  => z1d
case(fid%vd_tcwv)
  this%vd_tcwv  => z1d
case(fid%vd_tclw)
  this%vd_tclw  => z1d
case(fid%vd_tciw)
  this%vd_tciw  => z1d
case(fid%vd_tcrw)
  this%vd_tcrw  => z1d
case(fid%vd_tcsw)
  this%vd_tcsw  => z1d
case(fid%vd_tcslw)
  this%vd_tcslw  => z1d
case(fid%vd_ssrd)
  this%vd_ssrd  => z1d
case(fid%vd_strd)
  this%vd_strd  => z1d
case(fid%vd_ssrdc)
  this%vd_ssrdc  => z1d
case(fid%vd_strdc)
  this%vd_strdc  => z1d
case(fid%vd_blh)
  this%vd_blh  => z1d
case(fid%vd_sund)
  this%vd_sund  => z1d
case(fid%vd_spar)
  this%vd_spar  => z1d
case(fid%vd_suvb)
  this%vd_suvb  => z1d
case(fid%vd_sfdir)
  this%vd_sfdir  => z1d
case(fid%vd_scdir)
  this%vd_scdir  => z1d
case(fid%vd_sdsrp)
  this%vd_sdsrp  => z1d
case(fid%vd_cape)
  this%vd_cape  => z1d
case(fid%vd_capes)
  this%vd_capes  => z1d
case(fid%vd_mucape)
  this%vd_mucape  => z1d
case(fid%vd_pdepl)
  this%vd_pdepl  => z1d
case(fid%vd_mlcape50)
  this%vd_mlcape50  => z1d
case(fid%vd_mlcape100)
  this%vd_mlcape100  => z1d
case(fid%vd_mlcin50)
  this%vd_mlcin50  => z1d
case(fid%vd_mlcin100)
  this%vd_mlcin100  => z1d
case(fid%vd_tropotp)
  this%vd_tropotp  => z1d
case(fid%vd_tsrc)
  this%vd_tsrc  => z1d
case(fid%vd_ttrc)
  this%vd_ttrc  => z1d
case(fid%vd_ssrc)
  this%vd_ssrc  => z1d
case(fid%vd_strc)
  this%vd_strc  => z1d
case(fid%vd_es)
  this%vd_es  => z1d
case(fid%vd_smlt)
  this%vd_smlt  => z1d
case(fid%vd_10fg)
  this%vd_10fg  => z1d
case(fid%vd_10fg6)
  this%vd_10fg6  => z1d
case(fid%vd_10fgcv)
  this%vd_10fgcv  => z1d
case(fid%vd_i10fg)
  this%vd_i10fg  => z1d
case(fid%vd_lspf)
  this%vd_lspf  => z1d
case(fid%vd_tco3)
  this%vd_tco3  => z1d
case(fid%vd_vimd)
  this%vd_vimd  => z1d
case(fid%vd_sparc)
  this%vd_sparc  => z1d
case(fid%vd_stinc)
  this%vd_stinc  => z1d
case(fid%vd_cbase)
  this%vd_cbase  => z1d
case(fid%vd_0degl)
  this%vd_0degl  => z1d
case(fid%vd_m10degl)
  this%vd_m10degl  => z1d
case(fid%vd_visih)
  this%vd_visih  => z1d
case(fid%vd_cin)
  this%vd_cin  => z1d
case(fid%vd_kindex)
  this%vd_kindex  => z1d
case(fid%vd_ttindex)
  this%vd_ttindex  => z1d
case(fid%vd_cbasea)
  this%vd_cbasea  => z1d
case(fid%vd_ctopc)
  this%vd_ctopc  => z1d
case(fid%vd_ztwetb0)
  this%vd_ztwetb0  => z1d
case(fid%vd_ztwetb1)
  this%vd_ztwetb1  => z1d
case(fid%vd_tcghg)
  this%vd_tcghg  => z1d
case(fid%vd_tcchem)
  this%vd_tcchem  => z1d
case(fid%vd_aerodiag)
  this%vd_aerodiag  => z2d
case(fid%vd_aero_wvl_diag)
  this%vd_aero_wvl_diag  => z2d
case(fid%vd_100u)
  this%vd_100u  => z1d
case(fid%vd_100v)
  this%vd_100v  => z1d
case(fid%vd_zust)
  this%vd_zust  => z1d
case(fid%vd_10nu)
  this%vd_10nu  => z1d
case(fid%vd_10nv)
  this%vd_10nv  => z1d
case(fid%vd_dndzn)
  this%vd_dndzn  => z1d
case(fid%vd_dndza)
  this%vd_dndza  => z1d
case(fid%vd_dctb)
  this%vd_dctb  => z1d
case(fid%vd_tplb)
  this%vd_tplb  => z1d
case(fid%vd_tplt)
  this%vd_tplt  => z1d
case(fid%vd_odss)
  this%vd_odss  => z1d
case(fid%vd_oddu)
  this%vd_oddu  => z1d
case(fid%vd_odom)
  this%vd_odom  => z1d
case(fid%vd_odbc)
  this%vd_odbc  => z1d
case(fid%vd_odsu)
  this%vd_odsu  => z1d
case(fid%vd_odsoa)
  this%vd_odsoa  => z1d
case(fid%vd_odni)
  this%vd_odni  => z1d
case(fid%vd_odam)
  this%vd_odam  => z1d
case(fid%vd_odvfa)
  this%vd_odvfa  => z1d
case(fid%vd_odtoacc)
  this%vd_odtoacc  => z1d
case(fid%vd_odvsu)
  this%vd_odvsu  => z1d
case(fid%vd_aepm1)
  this%vd_aepm1  => z1d
case(fid%vd_aepm25)
  this%vd_aepm25  => z1d
case(fid%vd_aepm10)
  this%vd_aepm10  => z1d
case(fid%vd_uvbed)
  this%vd_uvbed  => z1d
case(fid%vd_uvbedcs)
  this%vd_uvbedcs  => z1d
case(fid%vd_litoti)
  this%vd_litoti  => z1d
case(fid%vd_licgi)
  this%vd_licgi  => z1d
case(fid%vd_litota6)
  this%vd_litota6 => z1d
case(fid%vd_licga6)
  this%vd_licga6  => z1d
case(fid%vd_ptypeocc6)
  this%vd_ptypeocc6  => z2d
case(fid%vd_200u)
  this%vd_200u  => z1d
case(fid%vd_200v)
  this%vd_200v  => z1d
case(fid%vd_sdsl)
  this%vd_sdsl  => z1d

! Group SM
case(fid%sm_clbt)
  this%sm_clbt  => z2d
case(fid%sm_csbt)
  this%sm_csbt  => z2d

! Group VX
case(fid%vx_oro)
  this%vx_oro   => z1d
case(fid%vx_tsc)
  this%vx_tsc   => z1d
case(fid%vx_pws)
  this%vx_pws   => z1d
case(fid%vx_pwp)
  this%vx_pwp   => z1d
case(fid%vx_sno)
  this%vx_sno   => z1d
case(fid%vx_tpc)
  this%vx_tpc   => z1d
case(fid%vx_sab)
  this%vx_sab   => z1d
case(fid%vx_xd2)
  this%vx_xd2   => z1d
case(fid%vx_lsm)
  this%vx_lsm   => z1d
case(fid%vx_iveg)
  this%vx_iveg  => z1d
case(fid%vx_arg)
  this%vx_arg   => z1d
case(fid%vx_rsmin)
  this%vx_rsmin => z1d
case(fid%vx_lai)
  this%vx_lai   => z1d
case(fid%vx_veg)
  this%vx_veg   => z1d

! SLB2
case(fid%vvel)
  this%vvel   => z2d

case default

  call abor1('Missing variable in pointer mapping case statement')

end select
IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:MAIN_FIELD_MAP_STORAGE',1,ZHOOK_HANDLE)

end subroutine main_field_map_storage

subroutine read_dynamic_namespace
character(len=:),allocatable    :: clf_name
logical :: ll_exists
integer(kind=jpim) :: IUTMP=77,II,io_status,igrib,jj
character(len=JP_NAME_MAX_LEN) :: clname
character(len=JP_COMMENTS_MAX_LEN) :: clongname
REAL(KIND=JPHOOK)    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:READ_DYNAMIC_NAMESPACE',0,ZHOOK_HANDLE)

if(ndyn_names > 0) call abor1('field_definitions:read_dynamic_namespace - Can only be called once')
clf_name='dynamic_namespace'
inquire(file=clf_name,exist=ll_exists)
if(ll_exists) then
  OPEN(IUTMP,FILE=clf_name,FORM='FORMATTED',ACTION='READ',STATUS='OLD')
  ii = 0
  write(nulout,'(A,11X,A,26X,2A)') ' NAME',' COMMENT ','GRIBCODE',' FID'
  write(nulout,'(A,11X,A,26X,2A)') ' ----',' ------- ','--------',' ---'
  do
    read(iutmp,*,iostat=io_status)clname, clongname,igrib 
    if(io_status < 0) exit
    do jj=1,ii
      if(clname == dyn_nam(jj)%cname) then
        write(nulout,*) 'field_definitions:read_dynamic_namespace - duplicate name',clname
        call abor1('field_definitions:read_dynamic_namespace - duplicate name')
      endif
    enddo
    ii = ii+1
    if(compo_3d_first+ii-1 > compo_3d_last) then
      call abor1('field_definitions:read_dynamic_namespace - too many dynamic names')
    endif
    dyn_nam(ii)%ifid=compo_3d_first+ii-1
    dyn_nam(ii)%cname=clname
    dyn_nam(ii)%clongname=clongname
    dyn_nam(ii)%igribcode=igrib
    write(nulout,'(1X,A,A,I8,1X,I6)') dyn_nam(ii)%cname,dyn_nam(ii)%clongname(1:32),dyn_nam(ii)%igribcode,&
     & dyn_nam(ii)%ifid
  enddo
  ndyn_names = ii
  write(nulout,*) 'FIELD_DEFINITIONS:READ_DYNAMIC_NAMESPACE - NAMES READ=',ndyn_names
else
  write(nulout,*) 'FIELD_DEFINITIONS:READ_DYNAMIC_NAMESPACE - No namspace file (',clf_name,') found'
endif

IF (LHOOK) CALL DR_HOOK('FIELD_DEFINITIONS:READ_DYNAMIC_NAMESPACE',1,ZHOOK_HANDLE)
end subroutine read_dynamic_namespace
end module field_definitions

