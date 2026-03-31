module lycom_nova
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!HOW TO UNDERTAKE THE PREPROCESING!!!!!!!!!!!!!!!!!!
!***********************************************************************
! NOVA_INIT
!***********************************************************************

!NEW: V8 -------------------------------------------------------
!
! update Hydrophobicity; Ref.Respiration; kB
!
! introduce LAI of NVV = F(zT)
!
! introduce Tmin for activity
!
! introduce liquid water uptake from below (soil & bark)
!
! rm bug in snow scheme
!
! new concept for water relations: 2-pool model + surface resistance! ........ + DCO2!!
!
! weighting by height or equal weights as alternative to netgrowth
!
! fix monthly cover update
!
! simplify albedo scheme
!
! account for rainfall when calculating water available for ETpot
!
!---------------------------------------------------------------

subroutine nova_init (nCPts3)
use lycom_par
use lycom_opt
implicit none


integer :: i,t,v,h,j,m,l
integer :: nCPts3
!write (*,*) "LYcophyte_init: is starting"
! initialise fields

!klife(:,:,:)                = 0   !!check if it is required!!    !CHANGED coordinate               ! switch between dead and alive
!rH2Ol_t(:,:,:,:,:)              = 0.0 !not required                                  ! thallus water content [m3 H2O / m2 T]
!rH2Ol_b(:,:,:,:,:)              = 0.0                                   ! thallus water content [m3 H2O / m2 T]

netgrowth(:,:,:)            = 0.0  !!Need to check!!             !CHANGED coordinate                             ! net growth [1 / ts]
!!!!!!AREATH_S!!!!! is not req most likely
area_s(:,:,:)             = 0.0  !!check if this is mandatory                                 ! surface cover [ m2 T / m2 V ]

!gpp0(i,:,:)                 = 0.0                                   ! initial GPP
!npp0(i,:,:)                 = 0.0                                   ! initial NPP


!the random parameter is diff from lycom

! initialise species parameters
!write (*,*) "The total number of species is ", p_nspec
!write (*,*) "The total number of speciesvariables is ", p_nspecpar

do j = 1,p_nspec ! loop over all species

!!!!!!!!!!!!***********IS albedo2 =palbveg????!!!!!!!!  

!!!!!!ALL varibales in this sector need to be defined

  o_albedo2(j)                  = 0.0!vec_o(j,4)*(p_alb_h-p_alb_l)+p_alb_l  ! pure NVV albedo []       *NEW GLOBAL VARIABLE
                                                                        
  !o_zt(j)                       = p_zt_l * exp(vec_o(j,2) &             ! thallus height (photosynthesising) [m]
   !                             * log(p_zt_h / p_zt_l))                 

  !o_prs(j)                      = vec_o(j,3) * (p_prs_h - p_prs_l) &    ! total thallus porosity when dry []
    !                            + p_prs_l

  !o_LAI(j)                      = (vec_o(j,2)*(p_LAInvv_h-p_LAInvv_l) & ! LAI of NVV []                                                         ! NEEDS EXP. VALIDATION I
     !                           + p_LAInvv_l)

  !o_rs(j)                       = (1.0-vec_o(j,3)) &                    ! thallus resistance to H2O [s / m]                                     ! NEEDS EXP. VALIDATION II
      !                          * (p_rs_h - p_rs_l) + p_rs_l

  !fracAir                       = vec_o(j,4) * (p_fracA_h -p_fracA_l) & ! fraction of air at saturation []                                              ! HOW TO DETERMINE IN LAB ?
       !                         + p_fracA_l

  o_spec_area(j)                = vec_o(j,5)*(13-2)+2!(0.500-0.075)+0.075   ! specific area [m2 T / kg C]                                           ! NEEDS EXP. VALIDATION III
                                
                                                                        
  !o_theta_max(j)                = o_prs(j)*(1.0-fracAir)*o_zt(j) &      ! water storage capacity [kg H2O / kg C]                                ! NEEDS EXP. VALIDATION III
         !                       * c_rhoH2Ol *o_spec_area(j) * 0.5 !calib

  !o_sat_X(j)                    = vec_o(j,5)                            ! water saturation at which potential becomes negative []               ! NEEDS EXP. VALIDATION II

  !o_sat_actF(j)                 = (1.0 -vec_o(j,5)) *(p_sat_actF_h &    ! water saturation needed for full activity []                          ! NEEDS EXP. VALIDATION IV
  !                              - p_sat_actF_l) + p_sat_actF_l

  !o_DCO2(j)                     = p_kCO2g_satl * exp(vec_o(j,4) &       ! DCO2 min [mol / (m2 T * s)]                                           ! NEEDS EXP. VALIDATION IV
  !                              * log(p_kCO2g_sath / p_kCO2g_satl))     

  !o_DCO2B(j)                    = vec_o(j,5) *(p_kCO2gB_h-p_kCO2gB_l) & ! DCO2 slope [ ]                                                        ! NEEDS EXP. VALIDATION IV
  !                              + p_kCO2gB_l

                                                                        
  !o_satHph(j)                   = vec_o(j,9)*(p_satHph_h-p_satHph_l) &  ! saturation below which hydrophobicity occurs []
          !                      + p_satHph_l
! write (*,*) "Te total number of species is ", p_nspec 
!  write (*,*) "The j is ", j
!  write (*,*) "The random number is from libry_specpar", vec_o(j,8)                                                                      
  o_vcmax_M(j)                  = p_vcmaxM_l * exp(vec_o(j,8) &         ! carboxylation rate of Rubisco (molar vcmax) [1 / s] p_vcmaxM_l=0.0139,p_vcmaxM_h=26.8   *all global var even the limits 
                                * log(p_vcmaxM_h / p_vcmaxM_l))         
                                                                        
  o_vomax_M(j)                  = vec_o(j,9) * (p_vomaxM_h &            ! oxygenation rate of Rubisco (molar vomax) [1 / s] p_vomaxM_l=0.391,p_vomaxM_h=2.5  *all global var even the limits
                                - p_vomaxM_l) + p_vomaxM_l              
                                                                        
  !o_resp_main(j)                = p_Rref_l * exp(vec_o(j,8) &           ! reference maintenance Respiration  [mol / (m2 T * s)] 
   !                             * log(p_Rref_h / p_Rref_l)) &           
    !                            / o_spec_area(j)

  o_RubConc(j)			=(Rub_h - Rub_l) * vec_o(j,10)+ Rub_l ! Rub_h=10e-6  Rub_l=8.5e-6                    *all global var even the limits

  o_ratio_Resp_Rub(j) 		= (resp_rub_h-resp_rub_l) * vec_o(j,12)+ resp_rub_l !resp_rub_h=0.1,resp_rub_l=0.01   *all global var even the limits

  o_gS0(j) 			= (gS0_h-gS0_l)*vec_o(j,2)+gS0_l !gS0_h=0.350,gS0_l=0.200                           *all global var even the limits
                                                                        
  !!!!o_spec_Rubisco(j)             = o_RubConc(j) * o_ratio_Resp_Rub(j)                 ! specific Rubisco content [mol / m2 T]                                 ! NEEDS EXP. VALIDATION V

  o_resp_main(j)		= o_RubConc(j) * o_ratio_Resp_Rub(j)                 ! specific Rubisco content [mol / m2 T] *all global var even the limits     ! NEEDS EXP. VALIDATION V
                                                                        
  !o_turnover(j)                 = p_turnover_l * exp(vec_o(j,8) &       ! turnover [1 / yr]                                                     ! NEEDS EXP. VALIDATION VI 




  o_X(j)=(0.99-0.9)*vec_o(j,7) + 0.9                                                                                 !*all global var even the limits

  o_G_area(j)=(((50-30)*vec_o(j,14)) + 30) / 10000                                                                     !*all global var even the limits

  o_w_leaf(j)= (1.0E-3-1.0E-4)*vec_o(j,14) + 1.0E-4                                                                   !*all global var even the limits

  o_A_leaf(j)=((0.15-0.01)*vec_o(j,14)+0.01)/1000000                                                               !*all global var even the limits

                                !* log(p_turnover_h / p_turnover_l)) 
  o_fracTransm(j) 		= (Trans_h-Trans_l)*vec_o(j,15)+Trans_l !Trans_l=0.60,Trans_h=0.75                 !*all global var even the limits                  
                                                                        
  o_ToptP(j)                    = vec_o(j,11) * (p_ToptP_max &          ! optimum temperature for photosynthesis [K] p_ToptP_max=323,p_ToptP_min=273 *all global var even the limits
                                - p_ToptP_min) + p_ToptP_min            
                                                                        
  o_Q10_resp(j)                 = vec_o(j,1) * (p_Q10R_h &             ! Q10 value of respiration [] p_Q10R_h=2.3,p_Q10R_l=1.5 *all global var even the limits
                                - p_Q10R_l) + p_Q10R_l                  

  o_Eact_Kc(j)                  = vec_o(j,3) * (p_EaKc_h - p_EaKc_l) & ! enzyme activation energy of Kc [J / mol] ! vec_o(j,12) p_EaKc_h=120000,p_EaKc_l=50000 *all global var even the limits
                                + p_EaKc_l                              
                                                                        
  o_Eact_Ko(j)                  = vec_o(j,3) * (p_EaKo_h - p_EaKo_l) & ! enzyme activation energy of Ko [J / mol] ! vec_o(j,13)   p_EaKo_h=50000,p_EaKo_l=10000 *all global var even the limits
                                + p_EaKo_l                              

  o_Eact_Vm(j)                  = vec_o(j,3) * (p_EaVm_h - p_EaVm_l) & ! enzyme activation energy of Vcmax [J / mol] ! vec_o(j,12)  p_EaVm_h=110000,p_EaVm_l=40000 *all global var even the limits 
                                + p_EaVm_l                              
                                                                        
  o_Eact_Jm(j)                  = vec_o(j,3) * (p_EaJm_h - p_EaJm_l) & ! enzyme activation energy of Vcmax [J / mol] ! vec_o(j,12) p_EaJm_h=80000,p_EaJm_l=30000   *all global var even the limits
                                + p_EaJm_l                              
                                                                        
                                                                        

enddo
do i = 1,nCPts3
  csum(i) = 0.0
  Lai_cum(i) =0.0
      
  do t = 1,p_ntiles

    do j = 1,p_nspec

      Rspec(j) =0.0
      fCcg_M(j)   = 0.0
      Rtres(i,j)=0.0
      fQ_tg(j)=0.0
      xT_g(i,t,j)=0.0
      counttimer(j)=1
      Layer_con(j)=0.5
      S2(j)=0.5
      Mrt_tot(i,j)=0.0
      as_fH2Ol_td(j)=0.0
      MD0(j)=0.0
      M0(i,j)=0.0
      RtoM(i,j)=0.0
      fCcb(j)=0.0
      Mrt(j)=0.0
      nfd(j)=0.0
      Als(i,j)=0.0
      Acs(i,j)=0.0
      fCbo(j)=0.0
      mon_cond(i,j)=0.0
      fH2Ol_ux1(i,t)=0.0
      klife(i,t,j)=0.0
      fH2Ol_xd(j)=0.0
      CO2_sink(i,j)=0.0
      CO2_pre(i,j)= 0.0
      netgrowth(i,t,j)            = 0.0  !!Need to check!!             !CHANGED coordinate                             ! net growth [1 / ts]
      area_s(i,t,j)= 0.0 
      gpp0(i,j)     = 0.0            ! GPP in [mol C / (m2 T * s)]
      npp0(i,j) =0.0
      Bl(i,j)=0.0
      bsum(i,j)=0.0
      xT_g(i,t,j)               = 288.0
      Wx1(i,t,j)=0.5*por*0.65
      Lai_new(i,j)=0.0
      Run_tot(i,t,j)=0.0
      Total_mortality_root(i,j)=0.0
      do l = 1,nsoil
        Qin1(i,t,j,l) =0.0
        Br(i,j,l)=0.0
        W_c1(i,t,j,l)=0.5*por*0.03
      enddo
      
    enddo

  enddo
enddo

frac_s_crit                   = frac_s_init/real(p_nspec) *fracratiocrit
!The alive constraints need to be changed!!!!! LOOK INTO THIS AFTER COMPLETING THE CODES!!SHOULD CONTAIN WITH RESPECT TO BIOMASS!!
!write (*,*) "LYcophyte_init: is ending"

return
end subroutine nova_init

!***********************************************************************
! lycom_STEP
!***********************************************************************

subroutine nova_step (i,t,v,h)
use lycom_par
use lycom_opt
implicit none
! temporary variables

integer ::i,t,v,h,j,l,m
integer ::kcccc
integer :: simhour
real    :: kH2Og
!real    :: csum
real    :: fracRADs2 !, albedo
real    :: RH_red, satb
real    :: kH2Og_sM
real    :: kCO2g_t
real    :: sO2, sCO2, P
real    :: ETpot, ETact, ETpot_can
real    :: gamma2
real    :: grh, cdg, crh
real    :: xT_s_wet, xT_s_dry
real    :: fracL
real    :: dew, wetfrac
real    :: waterUp0, waterUp, waterUp_b, Overflow
real    :: Rnet_v, ETpot_v, rootuptk, transpiration, soilevap, fH2Ol_tb_f1, w_rain_canopy
real    :: vcmaxTo, vcmax25
real    :: HcpO,Hcp
real    :: Resp, IR
!real    :: D1l, D2l, al, bl, cl, discl, xl, convLiq
real    :: D1c, D2c, bc, cc, discc, xc, K0,Jo, Ix
real    :: Bin,Bout
real    :: ngsum, wgtsum, expsum, hsum
real    :: retreat, disturbance
real    :: csum_area
real    :: rain
real    :: gSleaf, gSleaf2
real    :: dRAD, fRAD_Hw, fRAD_Hd,fracRADs_0,fracRADl_c
real    :: Rnet
real    :: a, b, d, xl, Al, Al_r, Ac, Ac_r, K
real    :: hmon
real    :: QRpot
real, parameter  :: wlim = 1.0e-4


kcccc=0
simhour = lastyear*365*24 
ETact=0.0
rCO2g_a                         = 360.0                                 ! [ppm]
rO2g_a                          = 210000.0                              ! [ppm]
!csum=0.0
Al=0.0 
Al_r =0.0
Ac =0.0
Ac_r=0.0
fH2Ol_tb_f1= 0.0
csum(i)    = 0.0
Lai_cum(i) = 0.0
!write (*,*) "LYcophyte_nova_crucial calculátions in lycophytes: is starting"

! Boundary layer conductance
                          
kH2Og                           = max( p_vonKarman * p_vonKarman &      ! [m / s]    !local                                                           REF: Allen,1998!
                                * max(p_critD,fAIR_s(i)) &
                                / dn_kH2Og(i,t), &
                                  0.004 )


rain=fH2Ol_ad(i) !(ESTONIA/SA/INDi)


!COMMENTED fH2Ol_tb_f1 =fH2Ol_tbf1
!commented !if (t .eq. 1) then  !forest

  !write (*,*) "Forest tile"
fracrain1(i,t)                 =2.0/p_LAImax
  !write (*,*) "The frac Rain", fracrain1(i,t)
  !write (*,*) "The Rain", rain*p_dt
fH2Ol_ci1(i,t)                 =max(0.0, rain*fracrain1(i,t)*0.35 * p_dt)!mul 0.65        ! water input into canopy [m3 H2O / (m2 C * s)]     !!!!!! THIS CANOPY WATER--may be utilised for direct evaporation
!write (*,*) "canopy water", fH2Ol_ci1(i,t)
w_rain_canopy= rain* p_dt * (1.0 - fracrain1(i,t)*0.35) !Water not entrapped by canopy
  !write (*,*) "The water required in the top soil is ", fH2Ol_tb_f1
  !write (*,*) "The rain water after canopy interception ", w_rain_canopy
!COMMENTED TOPSOILfH2Ol_ts1(i,t)    =max(0.0, min(w_rain_canopy, fH2Ol_tb_f1))  !!!0.0 !topsoilwater
  
  !write (*,*) "Top soil water", fH2Ol_ts1(i,t)

fH2Ol_ux1(i,t)                  =max(0.0, w_rain_canopy) ! water input into soil as throughfall after the loss in the topsoil [ m3 H2O / m2 G / s ]



if (Wmax .le. p_critD .or. Wmax .ne. Wmax) then
  write(*,*) "FATAL: invalid Wmax in nova_step"
  write(*,*) "rank=", rank, " i=", i, " t=", t, " Wmax=", Wmax
endif

if (Wxmax .le. p_critD .or. Wxmax .ne. Wxmax) then
  write(*,*) "FATAL: invalid Wxmax in nova_step"
  write(*,*) "rank=", rank, " i=", i, " t=", t, " Wxmax=", Wxmax
endif

if (Qp0 .ne. Qp0 .or. Qb0 .ne. Qb0 .or. p_dt .le. 0.0 .or. p_dt .ne. p_dt) then
  write(*,*) "FATAL: invalid hydraulic/time parameter in nova_step"
  write(*,*) "rank=", rank, " i=", i, " t=", t
  write(*,*) "Qp0=", Qp0, " Qb0=", Qb0, " p_dt=", p_dt
endif

!-----------------------------------------------------------------------
! Start loop over all species in a grid cell - I -
!-----------------------------------------------------------------------
!!!!!!!!!!!**********SEE to that I can access the Preprocessing here


if (v .eq. 2) then ! canopy
  lground                       = 0.0  !at canopy(not req for my case)
else
  lground                       = 1.0   !at groudlevel (soil)
endif

! No-vegetation run: disable all lycophyte processes
do j = 1, p_nspec
  klife(i,t,j)      = 0.0
  fCO2gc(i,j)       = 0.0
  fCO2nc(i,j)       = 0.0
  gpp0(i,j)         = 0.0
  npp0(i,j)         = 0.0
  netgrowth(i,t,j)  = 0.0
  area_s(i,t,j)     = 0.0
  Runoff1(i,t,j)    = 0.0
  Run_tot(i,t,j)    = 0.0
  Mrt_tot(i,j)      = 0.0
  CO2_sink(i,j)     = 0.0
  CO2_pre(i,j)      = 0.0
  fH2Ol_xd(j)       = 0.0
  fCcg_M(j)         = 0.0
  fCcb(j)           = 0.0
  fCbo(j)           = 0.0
  Rspec(j)          = 0.0
  Layer_con(j)      = 0.0
  S2(j)             = 0.0
enddo

csum(i)    = 0.0
Lai_cum(i) = 0.0
as_lai_s   = 0.0
as_area_s  = 0.0

if (writeout) then
  as_rCO2d          = 0.0
  as_sCO2d          = 0.0
  as_rCb            = 0.0
  as_Lai            = 0.0
  as_fCO2gc         = 0.0
  as_fCcg           = 0.0
  as_fCcb           = 0.0
  as_fCcb_l         = 0.0
  as_fCcb_c         = 0.0
  as_fCbo           = 0.0
  as_fH2Ol_lsat     = 0.0
  as_fH2Ol_bsat     = 0.0
  as_fH2Ol_runoff_l = 0.0
  as_fCc_npp        = 0.0
  as_fCc_gpp        = 0.0
  as_rH2Ol_t        = 0.0
  as_rmaxH2Ol_t     = 0.0
  as_Ts             = 0.0
  as_Tg             = 0.0
  as_H              = 0.0
  as_G              = 0.0
  as_E              = 0.0
  as_C              = 0.0
  as_EB             = 0.0
endif
fH2Ol_ci1(i,t) = 0.0
fH2Ol_ux1(i,t) = 0.0
return
end subroutine nova_step

end module lycom_nova
