
!#######################################################################
! lycom_COMMON
!#######################################################################

module lycom_common
contains


! lycom_NAMELIST
!-----------------------------------------------------------------------

subroutine lycom_namelist ()
use lycom_par
implicit none

integer :: stat

! define namelist parameters

namelist /lycompar/ year0,        &
                    cyear0,       &
                    tsindata0,    &
                    accts0,       &
                    tpos0,        &
                    lastyear,     &
                    runperiod,    &
                    tsl,          & ! [s]
                    yearout1,     &
                    yearoutX,     &
                    outint,       &
                    nSites,       &
                    p_nspec,      &
                    frac_s_init,  &
                    specout,      &
                    interCan,     &
                    lrestart,     &
                    NOHONO,       &
                    BSCtypes,     &
                    noVeg,        &
                    inDirect

! read namelist

open( knamelist, file=snamelist, status='old', iostat=stat )

if ( stat .eq. 0 ) then
  read( knamelist, lycompar )
  close( knamelist )
else
  write(*,*) "ERROR reading namelist, using default values"
endif

! check if tsl is le 1 h

if ( tsl .gt. 3600 ) then
  write(*,*) "ERROR: time step may not be longer than 1 hour"
  stop
endif

return
end subroutine lycom_namelist


! lycom_ALLOC
!-----------------------------------------------------------------------

subroutine lycom_alloc ()
use lycom_par
implicit none

! allocate variables

allocate( naccu(nCPts), &
naccu_d(nCPts), &
naccu_n(nCPts), &
naccu_m(nCPts), &
!
vec_o(p_nspec,p_nspecpar), &
o_gS0(p_nspec), &
o_albedo2(p_nspec), &
o_theta_max(p_nspec), &
o_spec_area(p_nspec), &
o_sat_actF(p_nspec), &
o_sat_X(p_nspec), &
o_vcmax_M(p_nspec), &
Run_tot(nCPts,p_ntiles,p_nspec), &
o_spec_Rubisco(p_nspec), &
o_ratio_Resp_Rub(p_nspec), &
o_RubConc(p_nspec), &
Mrt_tot(nCPts,p_nspec), &
SS(p_nspec), &
Als(nCPts,p_nspec), &
Acs(nCPts,p_nspec), &
o_resp_main(p_nspec), &
o_A_leaf(p_nspec), &
htree(p_nspec), &
Lai_new(nCPts,p_nspec), &
Lai_cum(nCPts), &
csum(nCPts), &
CO2_sink(nCPts,p_nspec), &
CO2_pre(nCPts,p_nspec), &
!CO2_sink(nCPts,p_ntiles,p_nspec), &
!CO2_pre(nCPts,p_ntiles,p_nspec), &
Rtres_l(p_nspec), &
as_fH2Ol_td(p_nspec), &
as_gpp(p_nspec), &
as_npp(p_nspec), &
Rtb(nCPts,p_nspec), &
Total_mortality_root(nCPts,p_nspec), &
RtoM(nCPts,p_nspec), &
nfd(p_nspec), &
M0(nCPts,p_nspec), &
Mrt(p_nspec), &
MD0(p_nspec), &
fracrain1(nCPts,p_ntiles), &
Qin(nCPts,p_ntiles,nsoil), &
Qin1(nCPts,p_ntiles,p_nspec,nsoil), &
fH2Ol_ci1(nCPts,p_ntiles), &
fH2Ol_ts1(nCPts,p_ntiles), &
fH2Ol_ux1(nCPts,p_ntiles), &
fracrain(nCPts,p_ntiles), &
fH2Ol_ci(nCPts,p_ntiles), &
fH2Ol_ts(nCPts,p_ntiles), &
fH2Ol_ux(nCPts,p_ntiles), &
o_w_leaf(p_nspec), &
o_X(p_nspec), &
o_G_area(p_nspec), &
o_fracTransm(p_nspec), &
o_turnover(p_nspec), &
o_vomax_M(p_nspec), &
o_ToptP(p_nspec), &
o_Q10_resp(p_nspec), &
o_Eact_Kc(p_nspec), &
o_Eact_Ko(p_nspec), &
o_Eact_Vm(p_nspec), &
o_Eact_Jm(p_nspec), &
N_leaf(p_nspec),&
o_CCM(p_nspec), &
o_DCO2(p_nspec), &
o_DCO2B(p_nspec), &
o_rs(p_nspec), &
o_zt(p_nspec), &
o_prs(p_nspec), &
o_satHph(p_nspec), &
fH2Olg_ga(p_nspec), &
o_LAI(p_nspec), &
mon_cond(nCPts,p_nspec), &
M_cond(p_nspec), &
counttimer(p_nspec), &
res_por(p_nspec),&
Etr_max(nCPts), &
Etr_avg(nCPts), &
!
xT_a(nCPts), &
fRADs_ad(nCPts), &
fRADl_ad(nCPts), &
fH2Ol_ad(nCPts), &
fH2Os_ad(nCPts), &
rH2Og_RH(nCPts), &
fAIR_s(nCPts), &
areaLEAF(nCPts,p_stepLAI), &
areaSTEM(nCPts,p_stepLAI), &
areaLEAF_month(nCPts), &
areaSTEM_month(nCPts), &
fracLEAF(nCPts), &
fracSTEM(nCPts), &
areaSOIL(nCPts), &
Acano(nCPts,p_ntiles), &
biome(nCPts), &
ETdavg(nCPts), &
ETrmax(nCPts),&
frac_noland(nCPts), &
frac_tile(nCPts,p_ntiles), &
albLsf(nCPts,p_ntiles), &
tauD(nCPts,p_ntiles), &
dn_kH2Og(nCPts,p_ntiles), &
rmaxH2Ol_b(nCPts,p_ntiles), &
rH2Os_g(nCPts), &
rH2Ol_0(nCPts,p_ntiles), &
rH2Ol_g1(nCPts,p_ntiles), &
rH2Ol_g2(nCPts,p_ntiles), &
xT_g0(nCPts,p_ntiles), &
!
CSOIL(nCPts), &
kSOIL(nCPts), &
dewmax(nCPts), &
!
rmaxH2Ol_t(p_nspec), &
sat(p_nspec), &
xT_s(p_nspec), &
CO2_p(p_nspec), &
jmax(p_nspec), &
vcmax(p_nspec), &
act(p_nspec), &
expansion(p_nspec), &
wgtspec(p_nspec), &
klife(nCPts,p_ntiles,p_nspec), &
rH2Ol_t(nCPts,p_ntiles,p_nspec), &
rH2Ol_b(nCPts,p_ntiles,p_nspec), &
gpp0(nCPts,p_nspec), &
!gpp(p_nspec), &
!gpp0(nCPts,p_ntiles,p_nspec), &
gpp(nCPts,p_ntiles,p_nspec), &
areaTH_s(nCPts,p_ntiles,p_nspec), &
area_s(nCPts,p_ntiles,p_nspec), &
xT_g(nCPts,p_ntiles,p_nspec), &
netgrowth(nCPts,p_ntiles,p_nspec), &
Q_base1(nCPts,p_ntiles,p_nspec), &
Q_base(nCPts,p_ntiles), &
!
fCO2gc_L(p_nspec), &
fCO2gc_W(p_nspec), &
fCO2gc(nCPts,p_nspec), &
fCO2nc(nCPts,p_nspec), &
fCcg_M(p_nspec), &
fCcg_G(p_nspec), &
fCcb(p_nspec), &
fCbo(p_nspec), &
fH2Ol_xd(p_nspec), &
fH2Olg_xu(p_nspec), &
fH2Ol_bx(p_nspec), &
fH2Ol_bd(p_nspec), &
fQ_tg(p_nspec), &
fQ_ta_L(p_nspec), &
fQ_ta_S(p_nspec), &
fRAD_H(p_nspec), &
W_c1(nCPts,p_ntiles,p_nspec,nsoil), &
W_c0(nCPts,p_ntiles,nsoil), &
S_wc1(p_nspec,nsoil), &
S2(p_nspec), &
Layer_con(p_nspec), &
!
fH2Ol_ug2(nCPts,p_ntiles), &
QR(nCPts,p_ntiles,nsoil), &
QR1(nCPts,p_ntiles,p_nspec,nsoil), &
Wx(nCPts,p_ntiles), &                                   
Wx1(nCPts,p_ntiles,p_nspec), &                                
Runoff(nCPts,p_ntiles), &                                
Runoff1(nCPts,p_ntiles,p_nspec), &
Rtres(nCPts,p_nspec), &
Rspec(p_nspec), &
!npp0(nCPts,p_ntiles,p_nspec), &
npp0(nCPts,p_nspec), &
bsum(nCPts,p_nspec), &
Bl(nCPts,p_nspec), &                              
Br(nCPts,p_nspec,nsoil), &                              
!                             
!Bl(nCPts,p_ntiles,p_nspec), &                              
!Br(nCPts,p_ntiles,p_nspec,nsoil), &                              
!
ah_areaTH_s(p_nhabM), &
ah_area_s(p_nhabM), &
!ah_fH2Ol_lsat(p_nhabM), &
ah_fH2Ol_runoff_l(p_nhabM), &
ah_fCc_gpp(p_nhabM), &
ah_fCc_npp(p_nhabM), &
ah_fH2Ol_bsat(p_nhabM), & 
ah_rH2Ol(p_nhabM), &
ah_rmaxH2Ol(p_nhabM), &
ah_act(p_nhabM), &
ah_fH2Ol_xd(p_nhabM), &
ah_fH2Ogl_ux(p_nhabM), &
ah_fH2Olg_xu(p_nhabM), &
ah_Ts(p_nhabM), &
ah_Tg(p_nhabM), &
ah_H(p_nhabM), &
ah_G(p_nhabM), &
ah_E(p_nhabM), &
ah_C(p_nhabM), &
ah_EB(p_nhabM), &
ah_rCO2d(p_nhabM), &
ah_sCO2d(p_nhabM), &
ah_Lai(p_nhabM), &
ah_rCb(p_nhabM), &
ah_fCO2gc(p_nhabM), &
ah_fCcg(p_nhabM), &
ah_fCcb(p_nhabM), &
ah_fCcb_l(p_nhabM), &
ah_fCcb_c(p_nhabM), &
ah_fCbo(p_nhabM), &
!
av_areaTH_s(p_nvertM), &
av_area_s(p_nvertM), &
!av_fH2Ol_lsat(p_nvertM), &
av_fH2Ol_runoff_l(p_nvertM), &
av_fCc_gpp(p_nvertM), &
av_fCc_npp(p_nvertM), &
av_fH2Ol_bsat(p_nvertM), &
av_rH2Ol(p_nvertM), &
av_rmaxH2Ol(p_nvertM), &
av_act(p_nvertM), &
av_fH2Ol_ux(p_nvertM), &
av_fH2Ol_xd(p_nvertM), &
av_fH2Ogl_ux(p_nvertM), &
av_fH2Olg_xu(p_nvertM), &
av_fH2Ol_bx(p_nvertM), &
av_fCO2gc(p_nvertM), &
av_Ts(p_nvertM), &
av_H(p_nvertM), &
av_E(p_nvertM), &
av_C(p_nvertM), &
av_EB(p_nvertM), &
av_rCO2d(p_nvertM), &
av_rCb(p_nvertM), &
av_fCcg(p_nvertM), &
av_fCcb(p_nvertM), &
av_fCcb_l(p_nvertM), &
av_fCcb_c(p_nvertM), &
av_fCbo(p_nvertM), &
!
at_fH2Ol_ux(p_ntiles,1), &
at_fH2Ol_xd(p_ntiles,1), &
!at_fH2Ol_lsat(p_ntiles,1), &
at_fH2Ol_runoff_l(p_ntiles,1), &
at_fCc_gpp(p_ntiles,1), &
at_fCc_npp(p_ntiles,1), &
at_fH2Ol_bsat(p_ntiles,1), &
at_areaTH_s(p_ntiles,1), &
at_area_s(p_ntiles,1), &
at_rH2Ol(p_ntiles,1), &
at_rmaxH2Ol(p_ntiles,1), &
at_act(p_ntiles,1), &
at_rCO2d(p_ntiles,1), &
at_sCO2d(p_ntiles), &
at_Lai(p_ntiles), &
at_rCb(p_ntiles,1), &
at_fCO2gc(p_ntiles,1), &
at_fCcg(p_ntiles,1), &
at_fCcb(p_ntiles,1), &
at_fCcb_l(p_ntiles,1), &
at_fCcb_c(p_ntiles,1), &
at_fCbo(p_ntiles,1), &
at_fH2Ogl_ux(p_ntiles,1), &
at_fH2Olg_xu(p_ntiles,1), &
at_fH2Ol_bx(p_ntiles,1), &
at_fH2Ol_bd(p_ntiles), &
at_fH2Ol_ug2(p_ntiles), &
at_Ts(p_ntiles,1), &
at_H(p_ntiles,1), &
at_E(p_ntiles,1), &
at_C(p_ntiles,1), &
at_EB(p_ntiles,1), &
at_Tg(p_ntiles), &
at_G(p_ntiles), &
at_rH2Ol_g1(p_ntiles), &
at_rH2Ol_g2(p_ntiles), &
at_fH2Ol_ug(p_ntiles), &
at_fH2Ol_go(p_ntiles), &
at_fH2Ol_gb(p_ntiles), &
at_fH2Olg_ga(p_ntiles), &
!
atN_fH2Ol_xd(nCPts,p_ntiles,p_nvertM), &
atN_fH2Ol_bd(nCPts,p_ntiles), &
!
ag_fH2Ol_ux(nCPts,p_ntiles), &
ag_fH2Ol_xd(nCPts,p_ntiles), &
ag_areaTH_s(nCPts,p_ntiles), &
ag_area_s(nCPts,p_ntiles), &
ag_fH2Ol_bsat(nCPts,p_ntiles), &
!ag_fH2Ol_lsat(nCPts,p_ntiles), &
ag_fH2Ol_runoff_l(nCPts,p_ntiles), &
ag_fCc_gpp(nCPts,p_ntiles), &
ag_fCc_npp(nCPts,p_ntiles), &
ag_rH2Ol(nCPts,p_ntiles), &
ag_rmaxH2Ol(nCPts,p_ntiles), &
ag_act(nCPts,p_ntiles), &
ag_rCO2d(nCPts,p_ntiles), &
ag_sCO2d(nCPts), &
ag_Lai(nCPts), &
ag_rCb(nCPts,p_ntiles), &
ag_fCO2gc(nCPts,p_ntiles), &
ag_fCcg(nCPts,p_ntiles), &
ag_fCcb(nCPts,p_ntiles), &
ag_fCcb_l(nCPts,p_ntiles), &
ag_fCcb_c(nCPts,p_ntiles), &
ag_fCbo(nCPts,p_ntiles), &
ag_fH2Ogl_ux(nCPts,p_ntiles), &
ag_fH2Olg_xu(nCPts,p_ntiles), &
ag_fH2Ol_bx(nCPts,p_ntiles), &
ag_fH2Ol_bd(nCPts), &
ag_fH2Ol_ug2(nCPts), &
ag_Ts(nCPts,p_ntiles), &
ag_H(nCPts,p_ntiles), &
ag_E(nCPts,p_ntiles), &
ag_C(nCPts,p_ntiles), &
ag_EB(nCPts,p_ntiles), &
ag_Tg(nCPts), &
ag_G(nCPts), &
ag_rH2Ol_g1(nCPts), &
ag_rH2Ol_g2(nCPts), &
ag_fH2Ol_ug(nCPts), &
ag_fH2Ol_go(nCPts), &
ag_fH2Ol_gb(nCPts), &
ag_fH2Olg_ga(nCPts), &
ag_rH2Os_g(nCPts), &
ag_fH2Osl_g(nCPts), &
ag_fH2Os_ad(nCPts), &
ag_xT_a(nCPts), &
ag_fH2Ol_ad(nCPts), &
ag_fRADs(nCPts), &
!
count_spec(nCPts), &
count_spec_h(p_nhabA,nCPts), &
vec_o_avg(p_nspecpar,p_nhabA,nCPts) )

return
end subroutine lycom_alloc


! lycom_DEALLOC
!-----------------------------------------------------------------------

subroutine lycom_dealloc ()
use lycom_par
implicit none

! deallocate variables

deallocate( naccu, &
naccu_d, &
naccu_n, &
naccu_m, &
!
vec_o, &
o_albedo2, &
o_theta_max, &
o_spec_area, &
o_sat_actF, &
o_resp_main, &
o_gS0, &
o_G_area, &
o_w_leaf, &
o_A_leaf, &
htree, &
Lai_new, &
csum, &
Lai_cum, &
CO2_sink, &
CO2_pre, &
as_fH2Ol_td, &
as_gpp, &
as_npp, &
Run_tot, &
Rtb, &
Rtres_l, &
RtoM, &
Total_mortality_root, &
nfd, &
M0, &
MD0, &
Mrt, &
fracrain1, &
fH2Olg_ga, &
Rspec, &
fH2Ol_ci1, &
fH2Ol_ts1, &
fH2Ol_ux1, &
fracrain, &
fH2Ol_ci, &
fH2Ol_ts, &
fH2Ol_ux, &
o_X, &
o_sat_X, &
o_fracTransm, &
o_vcmax_M, &
o_spec_Rubisco, &
N_leaf, &
o_ratio_Resp_Rub, &
o_RubConc, &
o_turnover, &
o_vomax_M, &
o_ToptP, &
o_Q10_resp, &
o_Eact_Kc, &
o_Eact_Ko, &
o_Eact_Vm, &
o_Eact_Jm, &
o_CCM, &
o_DCO2, &
o_DCO2B, &
o_rs, &
o_zt, &
o_prs, &
o_satHph, &
o_LAI, &
mon_cond, &
M_cond, &
res_por, &
counttimer,&
QR, &
QR1, &
Wx, &                                   
Wx1, & 
Mrt_tot, & 
SS, &
Als, &
Acs, &                              
Runoff, &                                
Runoff1, &
Rtres, &
npp0, &                              ! NPP [mol C / (m2 T * s)]
Bl, &                              ! leaf biomass [mol C / (m2 T * s)]
bsum, &
Br, &                              ! root biomass [mol C / (m2 T * s)]
!
xT_a, &
fRADs_ad, &
fRADl_ad, &
fH2Ol_ad, &
fH2Os_ad, &
rH2Og_RH, &
fAIR_s, &
areaLEAF, &
areaSTEM, &
areaLEAF_month, &
areaSTEM_month, &
fracLEAF, &
fracSTEM, &
areaSOIL, &
Etr_max, &
Etr_avg, &
Acano, &
biome, &
ETdavg, &
ETrmax, &
frac_noland, &
frac_tile, &
albLsf, &
tauD, &
dn_kH2Og, &
rmaxH2Ol_b, &
rH2Os_g, &
rH2Ol_0, &
rH2Ol_g1, &
rH2Ol_g2, &
xT_g0, &
Q_base, &
Q_base1, &
!
CSOIL, &
kSOIL, &
dewmax, &
!
rmaxH2Ol_t, &
sat, &
xT_s, &
CO2_p, &
jmax, &
vcmax, &
act, &
netgrowth, &
expansion, &
wgtspec, &
W_c1, &
W_c0, &
!
klife, &
rH2Ol_t, &
rH2Ol_b, &
gpp0, &
gpp, &
areaTH_s, &
area_s, &
xT_g, &
S_wc1, &
Layer_con, &
S2, &
!
fCO2gc_L, &
fCO2gc_W, &
fCO2gc, &
fCO2nc, &
fCcg_M, &
fCcg_G, &
fCcb, &
fCbo, &
fH2Ol_xd, &
fH2Olg_xu, &
fH2Ol_bx, &
fH2Ol_bd, &
fQ_tg, &
fQ_ta_L, &
fQ_ta_S, &
fRAD_H, &
!
fH2Ol_ug2, &
Qin, &
Qin1, &

!
ah_areaTH_s, &
ah_area_s, &
!ah_fH2Ol_lsat, &
ah_fH2Ol_runoff_l, &
ah_fCc_gpp, &
ah_fCc_npp, &
ah_fH2Ol_bsat, &
ah_rH2Ol, &
ah_rmaxH2Ol, &
ah_act, &
ah_fH2Ol_xd, &
ah_fH2Ogl_ux, &
ah_fH2Olg_xu, &
ah_Ts, &
ah_Tg, &
ah_H, &
ah_G, &
ah_E, &
ah_C, &
ah_EB, &
ah_rCO2d, &
ah_sCO2d, &
ah_Lai, &
ah_rCb, &
ah_fCO2gc, &
ah_fCcg, &
ah_fCcb, &
ah_fCcb_l, &
ah_fCcb_c, &
ah_fCbo, &
!
av_areaTH_s, &
av_area_s, &
av_rH2Ol, &
!av_fH2Ol_lsat, &
av_fH2Ol_runoff_l, &
av_fCc_gpp, &
av_fCc_npp, &
av_fH2Ol_bsat, &
av_rmaxH2Ol, &
av_act, &
av_fH2Ol_ux, &
av_fH2Ol_xd, &
av_fH2Ogl_ux, &
av_fH2Olg_xu, &
av_fH2Ol_bx, &
av_Ts, &
av_H, &
av_E, &
av_C, &
av_EB, &
av_rCO2d, &
av_rCb, &
av_fCcg, &
av_fCcb, &
av_fCcb_l, &
av_fCcb_c, &
av_fCbo, &
av_fCO2gc, &
!
at_fH2Ol_ux, &
at_fH2Ol_xd, &
at_areaTH_s, &
at_area_s, &
!at_fH2Ol_lsat, &
at_fH2Ol_runoff_l, &
at_fCc_gpp, &
at_fCc_npp, &
at_fH2Ol_bsat, &
at_rH2Ol, &
at_rmaxH2Ol, &
at_act, &
at_rCO2d, &
at_sCO2d, &
at_Lai, &
at_rCb, &
at_fCO2gc, &
at_fCcg, &
at_fCcb, &
at_fCcb_l, &
at_fCcb_c, &
at_fCbo, &
at_fH2Ogl_ux, &
at_fH2Olg_xu, &
at_fH2Ol_bx, &
at_fH2Ol_bd, &
at_Ts, &
at_H, &
at_E, &
at_C, &
at_EB, &
at_Tg, &
at_G, &
at_rH2Ol_g1, &
at_rH2Ol_g2, &
at_fH2Ol_ug, &
at_fH2Ol_go, &
at_fH2Ol_gb, &
at_fH2Olg_ga, &
at_fH2Ol_ug2, &
!
atN_fH2Ol_xd, &
atN_fH2Ol_bd, &
!
ag_fH2Ol_ux, &
ag_fH2Ol_xd, &
ag_areaTH_s, &
ag_area_s, &
!ag_fH2Ol_lsat, &
ag_fH2Ol_runoff_l, &
ag_fCc_gpp, &
ag_fCc_npp, &
ag_fH2Ol_bsat, &
ag_rH2Ol, &
ag_rmaxH2Ol, &
ag_act, &
ag_rCO2d, &
ag_sCO2d, &
ag_Lai, &
ag_rCb, &
ag_fCO2gc, &
ag_fCcg, &
ag_fCcb, &
ag_fCcb_l, &
ag_fCcb_c, &
ag_fCbo, &
ag_fH2Ogl_ux, &
ag_fH2Olg_xu, &
ag_fH2Ol_bx, &
ag_fH2Ol_bd, &
ag_Ts, &
ag_H, &
ag_E, &
ag_C, &
ag_EB, &
ag_Tg, &
ag_G, &
ag_rH2Ol_g1, &
ag_rH2Ol_g2, &
ag_fH2Ol_ug, &
ag_fH2Ol_go, &
ag_fH2Ol_gb, &
ag_fH2Olg_ga, &
ag_fH2Ol_ug2, &
ag_rH2Os_g, &
ag_fH2Osl_g, &
ag_fH2Os_ad, &
ag_xT_a, &
ag_fH2Ol_ad, &
ag_fRADs, &
!
count_spec, &
count_spec_h, &
vec_o_avg )

return
end subroutine lycom_dealloc


! lycom_READSPEC -- READ STRATEGY PARAMETERS
!-----------------------------------------------------------------------

subroutine lycom_readSpec ()
use lycom_par
implicit none

integer :: stat,j

character(9) :: fd0

write( fd0 ,'(I9)' ) p_nspecpar

! read species parameter file with random numbers [0, 1]

open( kspecpar, file=sspecpar, status='old', action='read', iostat=stat )

if ( stat .ne. 0 ) then
  write(*,*) "ERROR opening species parameter file"
  stop
endif

do j = 1,p_nspec
  read( kspecpar, '('//fd0//'F6.3)' ) vec_o(j,:)
enddo

close( kspecpar )
write (*,*)"P_nspec after reaSpec in common is", p_nspec

return
end subroutine lycom_readSpec


! lycom_READMON -- READ MONTHLY INPUT
!-----------------------------------------------------------------------

subroutine lycom_readMon(m)
use lycom_par
implicit none

integer :: m

! assign monthly canopy properties
areaLEAF_month(:) = areaLEAF(:,m)
areaSTEM_month(:) = areaSTEM(:,m)

return
end subroutine lycom_readMon


! lycom_RESET
!-----------------------------------------------------------------------

subroutine lycom_reset ()
use lycom_par
implicit none

  naccu(:) = 0

  naccu_d(:) = 0
  naccu_n(:) = 0

  ag_fH2Ol_ux(:,:)              = 0.0
  ag_fH2Ol_xd(:,:)              = 0.0

  ag_area_s(:,:)                = 0.0
  !ag_fH2Ol_lsat(:,:)            = 0.0
  ag_fH2Ol_runoff_l(:,:)        = 0.0
  ag_fCc_gpp(:,:)               = 0.0
  ag_fCc_npp(:,:)               = 0.0
  ag_fH2Ol_bsat(:,:)            = 0.0
  ag_rH2Ol(:,:)                 = 0.0
  ag_rmaxH2Ol(:,:)              = 0.0
  ag_act(:,:)                   = 0.0
  ag_rCO2d(:,:)                 = 0.0
  ag_sCO2d(:)                 = 0.0
  ag_Lai(:)                  = 0.0
  ag_rCb(:,:)                   = 0.0
  ag_fCO2gc(:,:)                = 0.0
  ag_fCcg(:,:)                  = 0.0
  ag_fCcb(:,:)                  = 0.0
  ag_fCcb_l(:,:)                = 0.0
  ag_fCcb_c(:,:)                = 0.0
  ag_fCbo(:,:)                  = 0.0

  ag_fH2Ogl_ux(:,:)             = 0.0
  ag_fH2Olg_xu(:,:)             = 0.0
  ag_fH2Ol_bx(:,:)              = 0.0

  ag_fH2Ol_bd(:)                = 0.0
  ag_fH2Ol_ug2(:)               = 0.0

  ag_Ts(:,:)                    = 0.0
  ag_H(:,:)                     = 0.0
  ag_E(:,:)                     = 0.0
  ag_C(:,:)                     = 0.0
  ag_EB(:,:)                    = 0.0

  ag_Tg(:)                      = 0.0
  ag_G(:)                       = 0.0
  
  ag_rH2Ol_g1(:)                = 0.0
  ag_rH2Ol_g2(:)                = 0.0
  
  ag_fH2Ol_ug(:)                = 0.0
  ag_fH2Ol_go(:)                = 0.0
  ag_fH2Ol_gb(:)                = 0.0
  ag_fH2Olg_ga(:)               = 0.0
  
  ag_rH2Os_g(:)                 = 0.0
  ag_fH2Osl_g(:)                = 0.0
  
  ag_fH2Os_ad(:)                = 0.0
  ag_xT_a(:)                    = 0.0
  ag_fH2Ol_ad(:)                = 0.0
  ag_fRADs(:)                   = 0.0

return
end subroutine lycom_reset

! lycom_AV_OUTPUT
!-----------------------------------------------------------------------

subroutine lycom_av_output(nCPts2)
use lycom_par
! use lycom_opt
implicit none

integer :: zero_count, i
integer :: nCPts2
zero_count = 0
do i = 1, nCPts2
  if (naccu(i) .le. 0) then
    zero_count = zero_count + 1
  endif
enddo
if (zero_count .gt. 0) then
  write(*,*) "WARNING: Rank", rank, " has", zero_count, " points with naccu <= 0"
  write(*,*) "  Total points:", nCPts2
  write(*,*) "  First few naccu values:", naccu(1:min(10,nCPts2))
endif






do i = 1,nCPts2

  if (naccu(i) .le. 0) then
    write(*,*) "WARNING: naccu(i) is zero or negative at i=", i, " naccu=", naccu(i)
    write(*,*) "Skipping averaging for this point"
    ag_fH2Ol_ux(i,:)              = -9999.0  !ag_fH2Ol_ux(i,:) / real(naccu(i))
    ag_fH2Ol_xd(i,:)              = -9999.0  !ag_fH2Ol_xd(i,:) / real(naccu(i))

    ag_area_s(i,:)                =  -9999.0  !ag_area_s(i,:) / real(naccu(i))
    !ag_fH2Ol_lsat(i,:)            = -9999.0   !ag_fH2Ol_lsat(i,:) / real(naccu(i))
    ag_fH2Ol_runoff_l(i,:)        = -9999.0   !ag_fH2Ol_runoff_l(i,:)/real(naccu(i))
    ag_fCc_gpp(i,:)               = -9999.0   !ag_fCc_gpp(i,:)/real(naccu(i))
    ag_fCc_npp(i,:)               = -9999.0   !ag_fCc_npp(i,:)/real(naccu(i))
    ag_fH2Ol_bsat(i,:)            = -9999.0   !ag_fH2Ol_bsat(i,:) / real(naccu(i))
  !ag_rH2Ol(i,:)                 = ag_rH2Ol(i,:)    / real(naccu(i))
  !ag_rmaxH2Ol(i,:)              = ag_rmaxH2Ol(i,:) / real(naccu(i))
  !ag_act(i,:)                   = ag_act(i,:)      / real(naccu(i))
    ag_rCO2d(i,:)                 = -9999.0  !ag_rCO2d(i,:)    / real(naccu(i))
    ag_sCO2d(i)                 = -9999.0    !ag_sCO2d(i)    / real(naccu(i))
    ag_Lai(i)                   =-9999.0  !ag_Lai(i) /real(naccu(i))
    ag_rCb(i,:)                   = -9999.0  ! ag_rCb(i,:)      / real(naccu(i))
    ag_fCO2gc(i,:)                = -9999.0  !ag_fCO2gc(i,:)   / real(naccu(i))
    ag_fCcg(i,:)                  = -9999.0   !ag_fCcg(i,:)     / real(naccu(i))
    ag_fCcb(i,:)                  = -9999.0  !ag_fCcb(i,:)     / real(naccu(i))
    ag_fCcb_l(i,:)                = -9999.0  !ag_fCcb_l(i,:)   / real(naccu(i))
    ag_fCcb_c(i,:)                = -9999.0  !ag_fCcb_c(i,:)   / real(naccu(i))
    ag_fCbo(i,:)                  = -9999.0  !ag_fCbo(i,:)     / real(naccu(i))

    ag_fH2Ogl_ux(i,:)             = -9999.0  !ag_fH2Ogl_ux(i,:)/ real(naccu(i))
    ag_fH2Olg_xu(i,:)             = -9999.0  !ag_fH2Olg_xu(i,:)/ real(naccu(i))
    ag_fH2Ol_bx(i,:)              = -9999.0  !ag_fH2Ol_bx(i,:) / real(naccu(i))

    ag_fH2Ol_bd(i)                = -9999.0  !ag_fH2Ol_bd(i) / real(naccu(i))
    ag_fH2Ol_ug2(i)               = -9999.0  !ag_fH2Ol_ug2(i) / real(naccu(i))

    ag_Ts(i,:)                    = -9999.0  !ag_Ts(i,:)       / real(naccu(i))      ! _d
    ag_H(i,:)                     = -9999.0  !ag_H(i,:)        / real(naccu(i))
    ag_E(i,:)                     = -9999.0  !ag_E(i,:)        / real(naccu(i))
    ag_C(i,:)                     = -9999.0  !ag_C(i,:)        / real(naccu(i))
    ag_EB(i,:)                    = -9999.0  !ag_EB(i,:)       / real(naccu(i))

    ag_Tg(i)                      = -9999.0  !ag_Tg(i)        / real(naccu(i))
    ag_G(i)                       = -9999.0  !ag_G(i)         / real(naccu(i))

    ag_rH2Ol_g1(i)                = -9999.0  !ag_rH2Ol_g1(i)  / real(naccu(i))
    ag_rH2Ol_g2(i)                = -9999.0  !ag_rH2Ol_g2(i)  / real(naccu(i))

    ag_fH2Ol_ug(i)                = -9999.0  !ag_fH2Ol_ug(i)  / real(naccu(i))
    ag_fH2Ol_go(i)                = -9999.0  !ag_fH2Ol_go(i)  / real(naccu(i))
    ag_fH2Ol_gb(i)                = -9999.0  !ag_fH2Ol_gb(i)  / real(naccu(i))
    ag_fH2Olg_ga(i)               = -9999.0  !ag_fH2Olg_ga(i) / real(naccu(i))

    ag_rH2Os_g(i)                 = -9999.0  !ag_rH2Os_g(i)   / real(naccu(i))
    ag_fH2Osl_g(i)                = -9999.0  !ag_fH2Osl_g(i)  / real(naccu(i))

    ag_fH2Os_ad(i)                = -9999.0  !ag_fH2Os_ad(i)  / real(naccu(i))
    ag_xT_a(i)                    = -9999.0  !ag_xT_a(i)      / real(naccu(i))
    ag_fH2Ol_ad(i)                = -9999.0  !ag_fH2Ol_ad(i)  / real(naccu(i))
    ag_fRADs(i)                   = -9999.0  !ag_fRADs(i)     / real(naccu(i))





    cycle  ! Skip this iteration
  endif
  ag_fH2Ol_ux(i,:)              = ag_fH2Ol_ux(i,:) / real(naccu(i))
  ag_fH2Ol_xd(i,:)              = ag_fH2Ol_xd(i,:) / real(naccu(i))

  ag_area_s(i,:)                = ag_area_s(i,:) / real(naccu(i))
  !ag_fH2Ol_lsat(i,:)            = ag_fH2Ol_lsat(i,:) / real(naccu(i))
  
  if (naccu(i) > 0) then
    ag_fH2Ol_runoff_l(i,:) = ag_fH2Ol_runoff_l(i,:)/real(naccu(i), kind(1.0d0))
  else
    ag_fH2Ol_runoff_l(i,:) = -9999.0   ! or 0.0d0, whichever your convention is
  end if

!ag_fH2Ol_runoff_l(i,:)        = ag_fH2Ol_runoff_l(i,:)/real(naccu(i))
  ag_fCc_gpp(i,:)               = ag_fCc_gpp(i,:)/real(naccu(i))
  ag_fCc_npp(i,:)               = ag_fCc_npp(i,:)/real(naccu(i))
  ag_fH2Ol_bsat(i,:)            = ag_fH2Ol_bsat(i,:) / real(naccu(i))
  !ag_rH2Ol(i,:)                 = ag_rH2Ol(i,:)    / real(naccu(i))
  !ag_rmaxH2Ol(i,:)              = ag_rmaxH2Ol(i,:) / real(naccu(i))
  !ag_act(i,:)                   = ag_act(i,:)      / real(naccu(i))
  ag_rCO2d(i,:)                 = ag_rCO2d(i,:)    / real(naccu(i))
  ag_sCO2d(i)                 = ag_sCO2d(i)    / real(naccu(i))
  ag_Lai(i)                   =ag_Lai(i) /real(naccu(i))
  ag_rCb(i,:)                   = ag_rCb(i,:)      / real(naccu(i)) 
  ag_fCO2gc(i,:)                = ag_fCO2gc(i,:)   / real(naccu(i))
  ag_fCcg(i,:)                  = ag_fCcg(i,:)     / real(naccu(i))
  ag_fCcb(i,:)                  = ag_fCcb(i,:)     / real(naccu(i))
  ag_fCcb_l(i,:)                = ag_fCcb_l(i,:)   / real(naccu(i))
  ag_fCcb_c(i,:)                = ag_fCcb_c(i,:)   / real(naccu(i))
  ag_fCbo(i,:)                  = ag_fCbo(i,:)     / real(naccu(i))

  ag_fH2Ogl_ux(i,:)             = ag_fH2Ogl_ux(i,:)/ real(naccu(i))
  ag_fH2Olg_xu(i,:)             = ag_fH2Olg_xu(i,:)/ real(naccu(i))
  ag_fH2Ol_bx(i,:)              = ag_fH2Ol_bx(i,:) / real(naccu(i))

  ag_fH2Ol_bd(i)                = ag_fH2Ol_bd(i) / real(naccu(i))
  ag_fH2Ol_ug2(i)               = ag_fH2Ol_ug2(i) / real(naccu(i))

  ag_Ts(i,:)                    = ag_Ts(i,:)       / real(naccu(i))      ! _d
  ag_H(i,:)                     = ag_H(i,:)        / real(naccu(i)) 
  ag_E(i,:)                     = ag_E(i,:)        / real(naccu(i)) 
  ag_C(i,:)                     = ag_C(i,:)        / real(naccu(i)) 
  ag_EB(i,:)                    = ag_EB(i,:)       / real(naccu(i)) 

  ag_Tg(i)                      = ag_Tg(i)        / real(naccu(i)) 
  ag_G(i)                       = ag_G(i)         / real(naccu(i)) 

  ag_rH2Ol_g1(i)                = ag_rH2Ol_g1(i)  / real(naccu(i))
  ag_rH2Ol_g2(i)                = ag_rH2Ol_g2(i)  / real(naccu(i))

  ag_fH2Ol_ug(i)                = ag_fH2Ol_ug(i)  / real(naccu(i))
  ag_fH2Ol_go(i)                = ag_fH2Ol_go(i)  / real(naccu(i))
  ag_fH2Ol_gb(i)                = ag_fH2Ol_gb(i)  / real(naccu(i))
  ag_fH2Olg_ga(i)               = ag_fH2Olg_ga(i) / real(naccu(i))

  ag_rH2Os_g(i)                 = ag_rH2Os_g(i)   / real(naccu(i)) 
  ag_fH2Osl_g(i)                = ag_fH2Osl_g(i)  / real(naccu(i)) 

  ag_fH2Os_ad(i)                = ag_fH2Os_ad(i)  / real(naccu(i)) 
  ag_xT_a(i)                    = ag_xT_a(i)      / real(naccu(i)) 
  ag_fH2Ol_ad(i)                = ag_fH2Ol_ad(i)  / real(naccu(i)) 
  ag_fRADs(i)                   = ag_fRADs(i)     / real(naccu(i)) 
enddo

!  if (BSCtypes) call lycom_avBSC(i)
!  if (BSCtypes) call lycom_avgcBSC(i,i2,j)

return
end subroutine lycom_av_output

!***********************************************************************
! lycom_AV_SPECIES
!***********************************************************************

subroutine lycom_av_species(nCPts2)
use lycom_par
implicit none

integer :: i,k,l
integer :: nCPts2
real    :: Asum

! distribution of species parameters

vec_o_avg(:,:,:) = 0.0
count_spec_h(:,:) = 0.0
count_spec(:) = 0.0

do i = 1,nCPts2

  count_spec(i) = real(sum(klife(i,:,:)))

  do k = 1,p_nhabA ! total number of tile/level/habitat combinations

    select case( k )
      case( 1 ) ! forest floor
        Asum                    = sum(area_s(i,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,:)))

        if (Asum .gt. p_critD) then
          do l = 1,p_nspecpar
            vec_o_avg(l,k,i)    = sum( vec_o(:,l) * area_s(i,1,:) / Asum )
          end do
        else
          vec_o_avg(:,k,i)      = 0.5
        endif
      case( 2 ) ! forest canopy stems
        Asum                    = sum(area_s(i,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,:)))

        if (Asum .gt. p_critD) then
          do l = 1,p_nspecpar
            vec_o_avg(l,k,i)    = sum( vec_o(:,l) * area_s(i,1,:) / Asum )
          end do
        else
          vec_o_avg(:,k,i)      = 0.5
        endif
      case( 3 ) ! forest canopy leaves
        Asum                    = sum(area_s(i,1,:))
        count_spec_h(k,i)       = real(sum(klife(i,1,:)))

        if (Asum .gt. p_critD) then
          do l = 1,p_nspecpar
            vec_o_avg(l,k,i)    = sum( vec_o(:,l) * area_s(i,1,:) / Asum )
          end do
        else
          vec_o_avg(:,k,i)      = 0.5
        endif
      case( 4 ) ! bare / grass
        Asum                    = sum(area_s(i,2,:))
        count_spec_h(k,i)       = real(sum(klife(i,2,:)))

        if (Asum .gt. p_critD) then
          do l = 1,p_nspecpar
            vec_o_avg(l,k,i)    = sum( vec_o(:,l) * area_s(i,2,:) / Asum )
          end do
        else
          vec_o_avg(:,k,i)      = 0.5
        endif
      case default
        Asum                    = 0.0
        count_spec_h(k,i)       = 0.0
        vec_o_avg(l,k,i)        = 0.5
    end select
  end do
enddo ! loop over all grid cells

return
end subroutine lycom_av_species


end module lycom_common

