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

do j = 1,p_nspec ! loop over all species

!!!!!!!!!!!!***********IS albedo2 =palbveg????!!!!!!!!  

!!!!!!ALL varibales in this sector need to be defined

  o_albedo2(j)                  = vec_o(j,4)*(p_alb_h-p_alb_l)+p_alb_l  ! pure NVV albedo []       *NEW GLOBAL VARIABLE
                                                                        
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

  o_G_area(j)=(((30-10)*vec_o(j,14)) + 10) / 10000                                                                     !*all global var even the limits

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
  
  do t = 1,p_ntiles
    do j = 1,p_nspec
      Rtres(j)=0.0
      fQ_tg(i)=0.0
      xT_g(i,t,j)=0.0
      counttimer(j)=0
      Layer_con(j)=0.5
      S2(j)=0.5
      fCcg_M(j)=0.0
      Mrt_tot(j)=0.0
      as_fH2Ol_td(j)=0.0
      MD0(j)=0.0
      M0(j)=0.0
      RtoM(j)=0.0
      fCcb(j)=0.0
      Mrt(j)=0.0
      nfd(j)=0.0
      Als(j)=0.0
      Acs(j)=0.0
      fCbo(j)=0.0
      mon_cond(j)=0.0
      fH2Ol_ux1(i,t)=0.0
      !if (biome(i) .eq. 16) then
      !  klife(i,t,j) = 0.0
      !  area_s(i,t,j) = 0.0
      !else
      klife(i,t,j) = 0.0
      area_s(i,t,j)= frac_s_init / real(p_nspec)

      !endif

      fH2Ol_xd(j)=0.0
      CO2_sink(j)=0.0
      CO2_pre(j)= 0.0
      netgrowth(i,t,j)            = 0.0  !!Need to check!!             !CHANGED coordinate                             ! net growth [1 / ts]
      !area_s(i,t,j)= frac_s_init / real(p_nspec)
      gpp0(j)     = 0.001 / 1.0 /c_MCo2 /p_dt             ! GPP in [mol C / (m2 T * s)]
     ! gpp(j)     = 0.001 / 1.0 /c_MCo2 /p_dt
      npp0(j) =0.0
      Bl(j)=20
      bsum(j)=40.0
      xT_g(i,t,j)               = 288.0
      Wx1(i,t,j)=0.5*por*0.65
      Lai_new(j)=2.0
      Run_tot(i,t,j)=0.0
      Total_mortality_root(j)=0.0
      do l = 1,nsoil
        Br(j,l)=4
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

integer ::i,t,v,h,j,l
integer ::kcccc
real    :: kH2Og
real    :: csum
!real    :: Lai_cum
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
!real    :: rain
real    :: gSleaf, gSleaf2
real    :: dRAD, fRAD_Hw, fRAD_Hd,fracRADs_0,fracRADl_c
real    :: Rnet
real    :: a, b, d, xl, Al, Al_r, Ac, Ac_r, K
real    :: hmon

kcccc=0
ETact=0.0
rCO2g_a                         = 400.0                                 ! [ppm]
rO2g_a                          = 210000.0                              ! [ppm]
csum = 0.0
!Lai_cum = 0.0
Al=0.0 
Al_r =0.0
Ac =0.0
Ac_r=0.0
fH2Ol_tb_f1= 0.0


!write (*,*) "LYcophyte_nova_crucial calculátions in lycophytes: is starting"

! Boundary layer conductance
                          
kH2Og                           = max( p_vonKarman * p_vonKarman &      ! [m / s]    !local                                                           REF: Allen,1998!
                                * max(p_critD,fAIR_s(i)) &
                                / dn_kH2Og(i,t), &
                                  0.004 )


rain=fH2Ol_ad(i)!*0.001 !(ESTONIA/SA/IND)
fH2Ol_tb_f1 =fH2Ol_tbf1
!commented !if (t .eq. 1) then  !forest

  !write (*,*) "Forest tile"
fracrain1(i,t)                 =areaLEAF_month(i)/p_LAImax
  !write (*,*) "The frac Rain", fracrain1(i,t)
!write (*,*) "The Rain", rain*p_dt
fH2Ol_ci1(i,t)                 =max(0.0, rain*fracrain1(i,t)*0.35 * p_dt)!mul 0.65        ! water input into canopy [m3 H2O / (m2 C * s)]     !!!!!! THIS CANOPY WATER--may be utilised for direct evaporation
  !write (*,*) "canopy water", fH2Ol_ci1(i,t)
w_rain_canopy= rain* p_dt*(1.0 - fracrain1(i,t)*0.35)
  !write (*,*) "The water required in the top soil is ", fH2Ol_tb_f1
  !write (*,*) "The rain water after canopy interception ", w_rain_canopy
fH2Ol_ts1(i,t)    =max(0.0, min(w_rain_canopy, fH2Ol_tb_f1))  !!!0.0
  
  !write (*,*) "Top soil water", fH2Ol_ts1(i,t)

fH2Ol_ux1(i,t)                  =max(0.0, w_rain_canopy - fH2Ol_ts1(i,t)) ! water input into soil as throughfall after the loss in the topsoil [ m3 H2O / m2 G / s ]
  !write (*,*) "water into ground", fH2Ol_ux1(i,t)
if (fH2Ol_ts1(i,t) - p_rmaxH2Ol_g1 .le. 0.0) then 

  fH2Ol_tb_f1                =   0 !if filled
else
      
  fH2Ol_tb_f1                =  p_rmaxH2Ol_g1 - fH2Ol_ts1(i,t)  !if not filled

endif                      
  
 !commented section 
!else
!  !write (*,*) "Non forest tile"
!  fH2Ol_ci1(i,t)                 =0
!  !write (*,*) "canopy water", fH2Ol_ci1(i,t)
!  fH2Ol_ts1(i,t)                 =0.0!commentedmax(0.0 ,min(rain*p_dt, fH2Ol_tb_f1))       !0.05*rain    water on top of the soil('may be 5% of the rainfall')      !!!!!! THIS Water on top of the soil--may be utilised for direct evaporation

!  !write (*,*) "Top soil water", fH2Ol_ts1(i,t)
!  fH2Ol_ux1(i,t)                 =max(0.0,rain*p_dt - fH2Ol_ts1(i,t))    !water into soil
!  !write (*,*) "water into ground", fH2Ol_ux1(i,t)
!  if (fH2Ol_ts1(i,t) - p_rmaxH2Ol_g1 .le. 0.0) then 
!    fH2Ol_tb_f1                =0 !if filled
!  else
!      
!    fH2Ol_tb_f1                = (p_rmaxH2Ol_g1-fH2Ol_ts1(i,t))  !if not filled
!
!  endif                      
!  
!  !fH2Ol_ux1                 =rain-fH2Ol_ts1     !water into soil   'This becomes a bit more complicated in case of bareground and snow or frozen ground'
!endif




!-----------------------------------------------------------------------
! Start loop over all species in a grid cell - I -
!-----------------------------------------------------------------------
!!!!!!!!!!!**********SEE to that I can access the Preprocessing here


if (v .eq. 2) then ! canopy
  lground                       = 0.0  !at canopy(not req for my case)
else
  lground                       = 1.0   !at groudlevel (soil)
endif


do j = 1,p_nspec
  
  !write (*,*) "klife is", klife(i,t,j)
  wetfrac =0.0
  ! Check if lycophyte is alive

  if (Bl(j) .ge. 0.001 .and. sum(Br(j,:)) .ge. 0.001) then                         !Klife(i,t,j) .eq. 1) then
  
    Klife(i,t,j) = 1.0
   ! write (*,*) "The plant species still alive", j
! *** WATER and ENERGY BALANCE ***

    ! Light absorption by the canopy(Aslai -Surrounding LAI)
    !!!!!!!****INITIALISE ALL the variables below
    !!!!THIS ASLAI needs to be initialised!!!!
    fracRADs_0                    = 1.0 - exp(-p_beer_s * Lai_new(j))         !*all global var even the limits
    

    fracRADl_c                    = 1.0 - exp(-p_beer_l * Lai_new(j))         !*all global var even the limits


    fracRADs 			  =(1- fracRADs_0)* (1.0 - exp(-p_beer_s*Lai_new(j))) *(1.0-o_albedo2(j))    !*all global var even the limits

    fracRADl			  = (1-fracRADl_c)*exp(-p_beer_l * Lai_new(j))        !*all global var even the limits

! Maximum lichen water storage capacity: I initialised it to 0

     !!!!!******************************************************Upto here

! Saturation water vapour pressure
    Ta3                           =xT_a(i)*xT_a(i)*xT_a(i)           !local
    Ta4                           =Ta3*xT_a(i)                       !local
    zT_a                          =xT_a(i)- c_TH2Osl                 !local

    esatAIR 			  = p_esatAIR3 *exp(p_esatAIR1*zT_a / (p_esatAIR2 + zT_a))   !local
   
    desatdT 			  = exp(p_esatAIR1*zT_a/(p_esatAIR2 &     !local slope of saturation vapour pressure vs temperature relationship []
             		                   + zT_a)) &
             		                   * (p_esatAIR1 *p_esatAIR2 *p_esatAIR3 &
             		                   / ((p_esatAIR2 + zT_a) &
             		                   * (p_esatAIR2 + zT_a)))





    

    Rnet 			  = fRADs_ad(i)*(1.0-o_albedo2(j)) + p_eps*fRADl_ad(i) - p_eps*c_sigma* xT_a(i)**4  !local

    if (zT_a .gt. 0.0) then
      !Evapotranspiration[m/s]
      ETdata = 1.3 * Rnet * desatdT / ( desatdT + c_gamma) / c_HH2Olg / c_rhoH2Ol                                 !*global
      !write (*,*)"The Rnet", Rnet
      !write (*,*)"The desatdt", desatdT
      !write (*,*)"The c_gamma", c_gamma
      !write (*,*)"The h2olg", c_HH2Olg
      !write (*,*)"The c_rhoh2ol", c_rhoH2Ol
      !write (*,*) "The timestep ETD ", ETdata
    else

      ETdata = 0.0                                                                                       
    endif
    
    
    !ETratio = max(p_critD, ETdata /ETdavgdata(i)) !NEW COMMENTED-BELOW VALUE USED
    ETratio =  max(p_critD, ETdata /4.2976E-8)                                                                   !global
    !write (*,*) "The ETDavg", ETdavg
    !write (*,*) "The timestep ETD ", ETdata
    !gSleaf = min(o_gS0(j) , o_gS0(j)/(ETratio**(log(p_gS1) / log(ETrmaxdata(i)))))  !NEW COMMENTED-BELOW VALUE USEDi
    gSleaf =  min(o_gS0(j) , o_gS0(j)/(ETratio**(log(p_gS1) / log(6.8957))))                                            !global all
    
    gamma2 = c_gamma * (1.0 + kH2Og / gSleaf)                                                                     !global

    gSleaf2 = o_gS0(j) * Layer_con(j)


    gS = max(p_critD, min(gSleaf, gSleaf2))  

    !write (*,*) "The stomatal conductivity is ", gS   
                                                                      !global

    if (dsnow .ge. p_H2Os_crit) then             ! snow layer too thick                            
       !!!!!cHECK THIS OUT!!!!!!!!!!
      ETpot                     = 0.0                                                                             !global
      fH2Olg_xu(j)              = 0.0  !evaporation from thallus not req                                          !global
      fH2Ol_bx(j)               = 0.0  !bark water uptake.. not req                                               !global
      fH2Osl_g                  = min(3.22 * max(0.0,xT_a(i)-c_TH2Osl)& ! snow melt [m/s]                                                       REF: Bergstrom,1992
                                / c_day / 1000.0, &
                                  rH2Os_g(i) / p_dt + fH2Os_ad(i)/1000)

      rH2Os_g(i)                = max(0.0, rH2Os_g(i) &
                                + fH2Os_ad(i)/1000 * p_dt &
                                - fH2Osl_g * p_dt &
                                - rH2Os_g(i) * p_H2Os_loss*p_dt)
      fCO2gc_L(j)              = 0.0
      fCO2gc_W(j)              = 0.0
      wetfrac                  = 0.0
      dsnow                    = rH2Os_g(i) * c_rhoH2Ol / p_rhoH2Os    
      fH2Ol_bd(j)               = 0.0  !bark wter overflow                                                        !global
      fH2Ol_xd(j)               = fH2Ol_xd(j)+ 0.0   !runoff.......'what is this?--check out now'
      xT_s_dry                  = min(c_TH2Osl - 0.1, xT_a(i))                                                    !global
      xT_s(j)                   = xT_s_dry                                                                        !local
      fRAD_H(j)                 = 0.0
      fCbo(j)                   = 0.0                                                                            !local
      fQ_ta_L(j)                = 0.0                                                                             !local
      fQ_ta_S(j)                = 0.0
      Layer_con(j)              = sum(W_c1(i,t,j,:))/(Wmax*5)
      mon_cond(j)               = mon_cond(j) + Layer_con(j)
      counttimer(j)                =counttimer(j)+1
      S2(j)                     = max(0.0, Wx1(i,t,j) / Wxmax) 
      nfd(j)                    = nfd(j) + 1.0
      Rspec(j)                  = 0.0!o_resp_main(j) &                      ! [mol / (m2 T * s)]                                                    REF: Kruse,2010
                                      !      * o_Q10_resp(j) &
                                      !      **((xT_s(j) - o_ToptP(j)) / 10.0)

      fCO2nc(j)                   = 0.0
      fCO2gc(j)                   =0.0
      npp0(j)       =npp0(j)+ fCO2nc(j)
      gpp0(j)=  fCO2gc(j)
      fH2Olg_ga(j)     = 0.0
      fCcb(j)          =  0.0 !Same question 
      fQ_tg(j)                    = 0.0!kSOIL(i) * (xT_s(j) )& ! Ground heat flux [W / m2 T]
                                !/ p_dz_SOIL * lground

      xT_g(i,t,j)               = xT_g(i,t,j) +0.0!&
                                !+ fQ_tg(j) /CSOIL(i) /p_dz_SOIL *p_dt 
!!!!!More stuffs need to be added here                                                                         !local
    else
    !write(*,*) "The part when snow is less than critical" 
    
      grh		 = gamma2 / kH2Og  
                                                                       !local
      cdg 		 = c_CAIR *(desatdT+gamma2)                                                               !local
                             
      xT_s_wet 		 = ( grh * (fRADs_ad(i)*fracRADs +fracRADl *p_eps * fRADl_ad(i) &                         !local
              		+ ((1.0-fracRADl) + 3.0) *p_eps *c_sigma *Ta4 +kSOIL(i)/p_dz_SOIL*xT_g(i,t,j)) &
              		+ xT_a(i)*cdg -c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) ) /( grh*(4.0*p_eps*c_sigma*Ta3 &
              		+ kSOIL(i) /p_dz_SOIL ) + cdg )
      
      crh 		 = 1.0/(c_CAIR*kH2Og)                                                                     !local
                               
      xT_s_dry 		 = ( xT_a(i) +crh*(fRADs_ad(i)*fracRADs +fracRADl *p_eps *fRADl_ad(i) &                   !local
                	+ ((1.0-fracRADl)+ 3.0) *p_eps *c_sigma *Ta4 +kSOIL(i)/p_dz_SOIL*xT_g(i,t,j)) ) &
                	/ ( 1.0 + crh*(4.0*p_eps*c_sigma*Ta3 +kSOIL(i)/p_dz_SOIL) )
      
      
! Net radiation

      dRAD 		 = fRADs_ad(i)*fracRADs +fracRADl *p_eps *fRADl_ad(i) &                                   !local
           		+ ((1.0-fracRADl) + 3.0)*p_eps*c_sigma*Ta4

      fRAD_Hd 		= dRAD -4.0*p_eps*c_sigma*Ta3*xT_s_dry - kSOIL(i)*(xT_s_dry - xT_g(i,t,j)) /p_dz_SOIL   !local
    
      if (xT_s_dry .le. c_TH2Osl-5) then ! frozen surface
        !write (*,*) "Entered part with temp less than -5"
        fRAD_Hw                 = 0.0
        ETpot                = 0.0
        fH2Olg_xu(j)            = 0.0
        fH2Ol_bx(j)             = 0.0
        fH2Ol_bd(j)             = 0.0
        wetfrac                 = 0.0
        fH2Ol_xd(j)             =fH2Ol_xd(j) +0.0
        xT_s(j)                 = xT_s_dry        
        fH2Osl_g	     = 0.0                                                                                        !*global

        rH2Os_g(i)           = max(0.0, rH2Os_g(i) + fH2Os_ad(i)/1000 * p_dt- fH2Osl_g * p_dt - rH2Os_g(i) * p_H2Os_loss*p_dt )  ! %snow layer [m3 H2O / m2 G] *global

        dsnow                = rH2Os_g(i) *c_rhoH2Ol/p_rhoH2Os !  %snow depth   *global
         
        fCbo(j)=  Mrt_tot(j)/c_year /c_MCo2 &         ! litterfall of photobiont [mol / (m2 T * s)]
                                / o_spec_area(j) ! * act(j)

!        fCcb(j)= gpp0(i,t,j) - (fCcg_M(j)*p_dt* c_MCo2)! check for
                         !Respiration under such low temp 
  !      Rspec(j)         = 0.0!o_resp_main(j) &                      ! [mol / (m2 T * s)]                                                    REF: Kruse,2010  !NEW COMMENTED
                               !  * o_Q10_resp(j) &
                                !   **((xT_s(j) - o_ToptP(j)) / 10.0)

        fCcg_M(j) = Rspec(j) !NEW COMMENTED


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!AADDEEDD please check
        do l =1,nsoil

          if (l .eq. 1) then

            Qin1(i,t,j,l) = 0.0 !!check if the coordinate system is right
        
          else
        
            Qin1(i,t,j,l)= Q_Per1 + Q_Oflow1
        
          endif

          
          W_c1(i,t,j,l)=W_c1(i,t,j,l)+Qin1(i,t,j,l)
          
          S_wc1(j,l) =  max(0.0,min(1.0,W_c1(i,t,j,l) / Wmax))                            !!! Check if W designation is right and must have 0.5 as initial value !!Layer relative water content
        
          Q_Per1 = max(0.0, min(Qp0 * S_wc1(j,l)*p_dt, W_c1(i,t,j,l)))            !Percolation from the layer
        

          W_c1(i,t,j,l) = W_c1(i,t,j,l) -  Q_Per1     !Water remaining in the layer
        
          Q_Oflow1  = max(0.0, W_c1(i,t,j,l)-Wmax)             !Overflow from the layer
          
          W_c1(i,t,j,l) = W_c1(i,t,j,l) - Q_Oflow1                !Water remaining in the layer
        
          W_con_l=W_con_l + W_c1(i,t,j,l)/Wmax
          
          ETact=ETact + 0.0
          
          nfd(j) = nfd(j) + 1.0 !NEW UNCOMMENTED
          
          if (S_wc1(j,l) .lt. 0.01) MD0(j) = MD0(j)+1
          
        
        enddo
        !nfd(j)= nfd(j) +1.0 !NEW COMMENTED
        Layer_con(j)=  sum(W_c1(i,t,j,:))/(Wmax*5)
        mon_cond(j)=mon_cond(j) + Layer_con(j)
        counttimer(j)=counttimer(j)+1
        W_con_l =0.0
      
        Rnet_v                = fRADs_ad(i)*0.85 +p_eps*fRADl_ad(i) &
                                - p_eps*c_sigma*Ta4
               
        
        fH2Ol_tb_f1 = p_rmaxH2Ol_g1-fH2Ol_ts1(i,t)

      !Bulk soil (bucket)

        Wx1(i,t,j)= Wx1(i,t,j) + Q_Oflow1 + Q_Per1
        
        
        
        Overflow= max(0.0, Wx1(i,t,j)-Wxmax)
        
        Wx1(i,t,j) = Wx1(i,t,j)-Overflow
        
        S2(j) = max(0.0,min(1.0,Wx1(i,t,j) / Wxmax))
        !write (*,*)"Saturationlevel Bucket", S2(j)
        
        Q_base1(i,t,j) = min(Qb0 * S2(j)*p_dt, Wx1(i,t,j) )      !Baseflow from bucket
     
        Wx1(i,t,j) = Wx1(i,t,j) - Q_base1(i,t,j) !Remaining water after baseflow loss
        
        S2(j) = max(0.0,Wx1(i,t,j) / Wxmax)  ! Relative water content in the bucket
      
        Runoff1(i,t,j)=Q_base1(i,t,j) + Overflow !!
        
!        fH2Olg_ga(j) =0.0  !fH2Olg_xu1 + transpiration + soilevap + ETact !NEW
!        COMMENTED


        fH2Ol_xd(j)=fH2Ol_xd(j)+max(0.0,Runoff1(i,t,j))

      else
      
        !write (*,*) "Enters part when the temp is not less than -5 degree C"
        fCcg_M(j)                 = Rspec(j)                         ! maintenance resp. [mol / (m2 T * s)]

        fCcb(j)                   = gpp0(j) - (fCcg_M(j)*p_dt* c_MCo2)          ! NPP [mol / (m2 T * s)]
!        write(*,*) "fCcb value is ", fCcb(j)
        if (fCcb(j) .gt. 0.0) fCcb(j) = fCcb(j) * Ngrow ! growth efficiency

        fCbo(j)                   =  Mrt_tot(j)/c_year /c_MCo2 &         ! litterfall of photobiont [mol / (m2 T * s)]
                                / o_spec_area(j) ! * act(j)

        
        fH2Olg_xu(j)            = 0.0
        fH2Ol_bx(j)             = 0.0
        fH2Ol_bd(j)             = 0.0
        fH2Ol_xd(j)             = 0.0 !NEW ADDED  
        fH2Osl_g                  = min(3.22 * max(0.0,xT_a(i)-c_TH2Osl)& ! snow melt [m/s]                                                       REF: Bergstrom,1992
                                / c_day / 1000.0, &
                                  rH2Os_g(i) / p_dt + fH2Os_ad(i)/1000)
        wetfrac                =0.0
        !write (*,*) 'wetfrac', wetfrac
        rH2Os_g(i)                      = max(0.0, rH2Os_g(i) &
                                + fH2Os_ad(i)/1000 * p_dt &
                                - fH2Osl_g * p_dt &
                                - rH2Os_g(i) * p_H2Os_loss*p_dt)


        dsnow                           = rH2Os_g(i) * c_rhoH2Ol / p_rhoH2Os
        
        fRAD_Hw = dRAD - 4.0*p_eps*c_sigma*Ta3*xT_s_wet -kSOIL(i) *(xT_s_wet - xT_g(i,t,j)) /p_dz_SOIL      !!CHEK HERE!!wet tem is used here but in lycom uses dry temp
      
        ! Potential Evapotranspiration Monteith 1981 m/s
        ETpot = ( fRAD_Hw * desatdT +c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) *kH2Og ) &
               / (desatdT+gamma2) / c_HH2Olg / c_rhoH2Ol
 	
! water uptake from below
        if (ETpot .le. 0.0) then

          wetfrac             = 1.0
        else
          wetfrac             =  ETact / ETpot/3600  
          
        endif
     

      

      !fH2Ol_tbf will increase with evaporation


        fH2Ol_ux1(i,t)= fH2Ol_ux1(i,t) + fH2Osl_g *p_dt  !add snow melt
        
        !write (*,*) "Water going into soil rain + snowmelt", fH2Ol_ux1(i,t)
        
        ETpot_can                    = ( fRAD_Hw * desatdT &                ! potential evaporation [m3 H2O / (m2 C * s)]                           REF: Monteith,1965
                                + c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) *kH2Og ) &
                                / (desatdT+c_gamma) / c_HH2Olg / c_rhoH2Ol
        
        fH2Olg_xu1            = max(0.0,min(fH2Ol_ci1(i,t),ETpot_can*p_dt))

        fH2Ol_xd1             = max(0.0, fH2Ol_ci1(i,t) - fH2Olg_xu1)

        !write (*,*) "water into soil from canopy after evaporation", fH2Ol_xd1
        !write (*,*) "Water evaporation from canopy", (ETpot_can*p_dt)
        !write (*,*) "Water present in canopy ", fH2Ol_ci1(i,t)
        fH2Ol_gwl1    =max(0.0, fH2Ol_ux1(i,t) + fH2Ol_xd1) !new variable

        !write (*,*) "Water going into the soil", fH2Ol_gwl1 
    !rH2Ol_t(i,t,v,h,j)      = max(0.0, rH2Ol_t(i,t,v,h,j) &         ! water reservoir [m3 H2O / m2 T]
     !                           + waterUp * p_dt &                      ! rainfall/throughfall + dew
      !                          + fH2Ol_bx(j) * p_dt &                  ! from below
       !                         - max(0.0, fH2Olg_xu(j)) * p_dt)

    
        do l =1,nsoil

          if (l .eq. 1) then

            Qin1(i,t,j,l) = max(0.0, fH2Ol_gwl1) !!check if the coordinate system is right
        
          else
        
            Qin1(i,t,j,l)= Q_Per1 + Q_Oflow1
        
          endif

          !write (*,*) "The water input in lycom layer", Qin1(i,t,j,l)

          W_c1(i,t,j,l)=W_c1(i,t,j,l)+Qin1(i,t,j,l)
        
          S_wc1(j,l) = max(0.0, min(1.0,W_c1(i,t,j,l) / Wmax))                            !!! Check if W designation is right and must have 0.5 as initial value !!Layer relative water content
        
          Q_Per1 = max(0.0,min(Qp0 * S_wc1(j,l)*p_dt, W_c1(i,t,j,l)))            !Percolation from the layer
        

          W_c1(i,t,j,l) = W_c1(i,t,j,l) -  Q_Per1     !Water remaining in the layer
        
          Q_Oflow1  = max(0.0, W_c1(i,t,j,l)-Wmax)             !Overflow from the layer
          
          W_c1(i,t,j,l) = W_c1(i,t,j,l) - Q_Oflow1                !Water remaining in the layer
        
          QR1(i,t,j,l) = min(ETpot, W_c1(i,t,j,l))   !Water lost pertaining to the evapotranpiration from the layer (do we need to put )
        
          W_c1(i,t,j,l) = W_c1(i,t,j,l) - QR1(i,t,j,l)  !remaining water in the layer

          ETact=ETact+ QR1(i,t,j,l)
        
        
          if (ETpot .lt. 0.0) then  !checking for the Water requirement of the plant
            ETpot=0.0
          else
            ETpot= ETpot -QR1(i,t,j,l)
          endif
          nfd(j) = nfd(j) + 1.0 !NEW UNCOMMENTED
          if (S_wc1(j,l) .lt. 0.01) MD0(j) = MD0(j)+1
          W_con_l= W_con_l+W_c1(i,t,j,l)/ Wmax  
        enddo
!        nfd(j)  =nfd(j) +1.0 !NEW COMMENTED

        SS(j)=sum(W_c1(i,t,j,:))/(Wmax*5)
       ! write (*,*) "The total sum of layer water content", W_con_l
        Layer_con(j)= sum(W_c1(i,t,j,:))/(Wmax*5)
        W_con_l =0.0
        mon_cond(j)=mon_cond(j) + Layer_con(j)
       ! write (*,*) "The layer water saturation is", Layer_con(j)
        counttimer(j)=counttimer(j)+1
      
        Rnet_v                        = fRADs_ad(i)*0.85 +p_eps*fRADl_ad(i) &
                                - p_eps*c_sigma*Ta4
               
        ETpot_v                     = 1.4 *Rnet_v*desatdT/(desatdT+c_gamma) & ! [m3 H2O / (m2 G * s)]
                                / c_HH2Olg/c_rhoH2Ol

        soilevap                   =min(fH2Ol_ts1(i,t), (max(0.0,ETpot_v)*p_dt))

        fH2Ol_ts1(i,t)=fH2Ol_ts1(i,t)-soilevap

        fH2Ol_tb_f1 = p_rmaxH2Ol_g1-fH2Ol_ts1(i,t)
        
        Wx1(i,t,j)= Wx1(i,t,j) + Q_Per1 + Q_Oflow1
        
        Overflow= max(0.0, Wx1(i,t,j)-Wxmax)
        
        Wx1(i,t,j) = Wx1(i,t,j)-Overflow
        
      !Bulk soil (bucket)
        rootuptk                    = min( Wx1(i,t,j)/p_dt, p_kH2Ol_sv & ! [m3 H2O / (m2 G * s)]         p_kH2Ol_sv = 5.0e-8
                                * (Wx1(i,t,j)/Wxmax)**2 )                  

        transpiration            = min( max(0.0,ETpot_v), rootuptk ) *p_dt    ! [m3 H2O / (m2 G * s)]
        
        Wx1(i,t,j)               = Wx1(i,t,j) - transpiration    ! [m3 H2O / m2 G]



        fH2Olg_ga(j) =fH2Olg_xu1 + transpiration + soilevap + ETact

        ETact=0.0
   
        S2(j) = max(0.0, Wx1(i,t,j) / Wxmax)  ! Relative water content in the bucket

        Q_base1(i,t,j) = min(Qb0 * S2(j)*p_dt, Wx1(i,t,j) )      !Baseflow from bucket
     
        Wx1(i,t,j) = max(0.0,Wx1(i,t,j) - Q_base1(i,t,j)) !Remaining water after baseflow loss
        
        
      
        Runoff1(i,t,j)=Q_base1(i,t,j) + Overflow !!THis needs to declared universally available
        
        S2(j) = max(0.0, Wx1(i,t,j) / Wxmax)

        fH2Ol_xd(j)=fH2Ol_xd(j)+max(0.0,Runoff1(i,t,j)) 
        
        
      endif 
    !write (*,*) 'wetfrac', wetfrac
      xT_s(j)                   = wetfrac * xT_s_wet + (1.0- wetfrac)*xT_s_dry

    !write (*,*) 'wetfrac', wetfrac
      fRAD_H(j)                 = wetfrac*fRAD_Hw + (1.0-wetfrac)*fRAD_Hd

    !fQ_ta_L(j)                = fH2Olg_xu(j) * c_HH2Olg * c_rhoH2Ol   ! Latent heat flux [W / m2 T]
      fQ_ta_L(j)                = fH2Olg_ga(j) * c_HH2Olg * c_rhoH2Ol   ! Latent heat flux [W / m2 T]     !!!!!!!check etact 


      ! Sensible heat
      fQ_ta_S(j)                = c_CAIR * (xT_s_wet - xT_a(i)) &       ! Sensible heat flux [W / m2 T]
                                * kH2Og * wetfrac &
                                + c_CAIR * (xT_s_dry - xT_a(i)) &
                                * kH2Og * (1.0 - wetfrac)
!!!!!!THE SENSIBLE AND GROUND HEAT NEEDS CORRECTION!!!!!
! Ground heat
      fQ_tg(j)                    = kSOIL(i) * (xT_s(j) - xT_g(i,t,j))& ! Ground heat flux [W / m2 T]
                                / p_dz_SOIL * lground

    ! Heat balance
      xT_g(i,t,j)               = xT_g(i,t,j) &
                                + fQ_tg(j) /CSOIL(i) /p_dz_SOIL *p_dt   ! balance of ground heat reservoir [K]
    !write (*,*) "The productivity calculation starts here"
    !write (*,*) "i value", i
      if (gS .gt. p_critD) then

      !write (*,*) "The productivity condition is entered"
       

        KcfT  = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
              * o_Eact_Kc(j) &
              / (o_ToptP(j) * c_Rgas * xT_s(j)))
      !write (*,*) "The temp in KcfT", xT_s(j)   
      !write (*,*) "The KcfT", KcfT                 
        KofT                       = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Ko(j) &
                                / (o_ToptP(j) * c_Rgas * xT_s(j)))
      !write (*,*) "The KofT", KofT
        VmfT                        = exp((xT_s(j) - o_ToptP(j)) &          ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Vm(j) &
                                / (o_ToptP(j) * c_Rgas * xT_s(j)))
      !write (*,*) "The VmfT", VmfT                          
        JmfT                        = exp((xT_s(j) - 298.15) &              ! temperature response of Michaelis-Menten-Constant []                  REF: Medlyn,2002
                                * o_Eact_Jm(j) &
                                / (xT_s(j) * c_Rgas * 298.15))
      !write (*,*) "The JmfT", JmfT
        vcmaxTo                     = o_vcmax_M(j) * o_RubConc(j)      ! vcmax [mol / (m2 T * s)] *** at T opt ***      !!!!!!SHOULD WE USE O_Rubconc(j) instead!!!! NEED TO CLARIFY- DIFFERENT IN lycom AND LYCOM
      !write (*,*) "The vcmaxTo", vcmaxTo                          
        vcmax(j)                    = vcmaxTo * VmfT                        ! vcmax at current Ts
      !write (*,*) "The vcmax", vcmax(j)                          
        vcmax25                     = vcmaxTo*exp((298.15 - o_ToptP(j)) &   ! vcmax at standard temperature
                                * o_Eact_Ko(j) &
                                / (o_ToptP(j) * c_Rgas * 298.15))
      !write (*,*) "The vcmax25", vcmax25
        jmax(j)                     = vcmax25 * 2.1 * JmfT                  ! jmax [mol / (m2 T * s)] at current Ts                                 REF: Wullschleger,1993
      !write (*,*) "The jmax", jmax(j)                          
      !commented recently !KcM                         = p_KcM1 * o_vcmax_M(j)**p_KcM2         ! Michaelis-Menten-Constant for CO2 [muM]                               REF: Savir,2009
      !!write (*,*) "The kcM", KcM                          
      !commented recently ! KoM                         = o_vomax_M(j) / (p_KoM1 &              ! Michaelis-Menten-Constant for O2 [muM]                                REF: Savir,2009
                                         !* (o_vcmax_M(j) / KcM)**p_KoM2)
      !!write (*,*) "The koM", KoM                          
        Kc                          = p_KcM1 * o_vcmax_M(j)**p_KcM2 * 0.001 * KcfT      !KcM * 0.001 * KcfT                    ! temperature correction [mol / m3]
      !write (*,*) "The kc", Kc 
      !write (*,*) "The vcmM is",o_vcmax_M(j)                         
        Ko                          = o_vomax_M(j) / (p_KoM1*(o_vcmax_M(j)/(p_KcM1*o_vcmax_M(j)**p_KcM2))**p_KoM2)*0.001*KofT    !KoM * 0.001 * KofT                    ! temperature correction [mol / m3]
      !write (*,*) "The ko", Ko

!!!!!!!!!!!!!!!

        Ix     			= cpar * fRADs_ad(i) * fracRADs * o_fracTransm(j)
      !write (*,*) "Ix", Ix
        Jo     			= jmax(j)*Ix*o_X(j) / (2.1*jmax(j) + Ix)
      !write (*,*) "Jo", Jo
        Rspec(j)                       = o_resp_main(j) &                      ! [mol / (m2 T * s)]                                                    REF: Kruse,2010
                                            * o_Q10_resp(j) &
                                            **((xT_s(j) - o_ToptP(j)) / 10.0)

      !write (*,*) "Respiration", Rspec(j)
      





        sO2   			= 0.00126*exp(1700.0*(1.0/xT_s(j) - 1.0/298.15))
      !write (*,*) "The sO2", sO2
        sCO2  			= 0.0334*exp(2400.0*(1.0/xT_s(j) - 1.0/298.15))
      !write (*,*) "The sCO2", sCO2
        P     			= 0.5*rO2g_a/1.0E6*sO2*1000.0*o_vomax_M(j)/ o_vcmax_M(j) * Kc/Ko
      !write (*,*) "The P", P


! %%%% light-limited rate
     
        a     = -4.5 * gS / (sCO2 * 1000)
      !write (*,*) "The a", a
     
        b     = gS * (4.5 * rCO2g_a /1.0e6 - 10.5 * P / (sCO2 *1000)) - Jo + 4.5 * Rspec(j)
      !write (*,*) "The b", b
        d     = P * (10.5 * gS * rCO2g_a /1.0e6 + Jo + 10.5 * Rspec(j) )
      !write (*,*) "The d", d
        xl    = (-b - sqrt(b**2 - 4.0*a * d)) / (2.0 * a)
      !write (*,*) "The xl", xl
        Al    = (Jo * (xl - P) / (4.5* xl + 10.5*P) -Rspec(j))

            
      !write (*,*) "The productivity LIght calculation is", Al
!%      Al    = (J * (xl-P) / (4.5*xl + 10.5*P) -R)*p_dt;
      
        Al_r  = (Jo * (xl - P) / (4.5* xl + 10.5*P))
      !write (*,*) "The productivity LIght GPP calculation is", Al_r
        fCO2gc_L(j)=Al

!      %%%% CO2-limited rate
     
        a     = -gS / (sCO2 * 1000)
      !write (*,*) "The a", a
        K     = Kc * (1.0 + rO2g_a /1.0E6 *sO2*1000.0 / Ko)
      !write (*,*) "The k", K
        b     = gS * (rCO2g_a /1.0e6 - K / (sCO2 *1000)) - vcmax(j) + Rspec(j)
      !write (*,*) "The b", b
        d     = (gS * rCO2g_a /1.0e6 + Rspec(j))*K + P*vcmax(j)
      !write (*,*) "The d", d
        xc    = (-b - sqrt(b**2 - 4.0*a * d)) / (2.0 * a)
      !write (*,*) "The xc", xc
        Ac    = ((vcmax(j) * (xc-P) / (xc + Kc*(1.0 + rO2g_a /1.0e6 * sO2 *1000  / Ko)) -Rspec(j)))
      !write (*,*) "The productivity Carbon calculation is", Ac
     
!%     Ac    = (vcm * (xc-P) / (xc + Kc*(1.0 + O2 /1.0e6 * sO2 *1000  / Ko)) -R);
        Ac_r  = (vcmax(j) * (xc-P) / (xc + Kc*(1.0 + rO2g_a /1.0e6 * sO2 *1000  / Ko)))
      !write (*,*) "The productivity carbon calculation GPP is", Ac_r
        fCO2gc_W(j)=Ac
      else
       ! Acs=Acs+1.0 commented out
        Ac_r= 0.0
        Al_r = 0.0
        Al=0.0
        Ac=0.0 
        !Rspec(j)                       = o_resp_main(j) &                      ! [mol / (m2 T * s)]                                                    REF: Kruse,2010
         !                                   * o_Q10_resp(j) &
          !                                  **((xT_s(j) - o_ToptP(j)) / 10.0)
          !                                  !NEW COMMENTED

       fCO2gc_L(j) = 0.0
       fQ_tg(j)                    = kSOIL(i) * (xT_s(j) - xT_g(i,t,j))& ! Ground heat flux [W / m2 T]
                                / p_dz_SOIL * lground


       xT_g(i,t,j)               = xT_g(i,t,j) &
                                + fQ_tg(j) /CSOIL(i) /p_dz_SOIL *p_dt 
       fCO2gc_W(j) = 0.0
      endif
     
    !write (*,*) "i value", i

      if (Al .lt. Ac) then
      
      
        fCO2gc(j)             = Al_r * c_MCo2 * Ngrow*p_dt
        fCO2nc(j)             =Al * c_MCo2 *  Ngrow*p_dt
        
        Als(j)                   =Als(j) +1.0
      else
      
        fCO2gc(j)             =Ac_r * c_MCo2* Ngrow*p_dt
        fCO2nc(j) 	    = Ac * c_MCo2* Ngrow*p_dt
        Acs(j)                   =Acs(j) +1.0
      
        if (Al .eq. Ac) then
          Acs(j)=Acs(j)+0.5!(Because 1 is already added before Acs=Acs-1+0.5)
          Als(j)=Als(j)+0.5
        endif
      endif
      bsum(j)                        =(Bl(j) +sum(Br(j,:)))

      M0(j)                          = M0(j) + kmrt0* bsum(j)

      
    
      gpp0(j)=  fCO2gc(j)
      npp0(j)= npp0(j)+ fCO2nc(j)
      if (biome(i) .eq. 16) then
        gpp0(j) = 0.0
        npp0(j) = 0.0
      endif      

!write (*,*) "The productivity end is reached and GPP is", fCO2gc(j)
      !write (*,*) "The productivity cumulative NPP is", npp0(i,t,j)
    endif
    fCcb(j)                   = gpp0(j) - (fCcg_M(j)*p_dt* c_MCo2)          ! NPP [mol / (m2 T * s)]
    !write(*,*) "fCcb value is ", fCcb(j)
!    if (fCcb(j) .gt. 0.0) fCcb(j) = fCcb(j) * Ngrow ! growth efficiency
   ! if (fCO2nc(j) .gt. 0.0) fCO2nc(j) = fCO2nc(j) * Ngrow

!    netgrowth(i,t,j)= netgrowth(i,t,j)+(fCO2nc(j)*p_dt)*o_spec_area(j)*c_MCo2
    netgrowth(i,t,j)   =  netgrowth(i,t,j)  &                ! [1 / month]
                          + (fCcb(j)*p_dt - fCbo(j)*p_dt) &
                          * o_spec_area(j) * c_MCo2
                          
    csum = csum + area_s(i,t,j)
    
    Lai_cum = Lai_cum + Lai_new(j)!!!!
    if (biome(i) .eq. 16) then
      fCcb(j) = 0.0
      netgrowth(i,t,j) = netgrowth(i,t,j) +0.0
      klife(i,t,j) = 0.0
    endif
  else
    !write (*,*) "The Klife is changed to 0. THe species dies"
    klife(i,t,j)=0.0

  endif  !klife!
  !write (*,*) "Species number ", j
  !write (*,*) "i value", i
  !write (*,*) "t", t
  
enddo
!-----------------------------------------------------------------------
! Execute once per month
!-----------------------------------------------------------------------
if (day .eq. dpm .and. ts .eq. tspd .and. Klife(i,t,j) .ne. 0.0) then

  do j = 1,p_nspec
    !bsum(j)                        =(Bl(j) +sum(Br(j,:)))

    !M0(j)                          = M0(j) + kMrt0* bsum(j)*dpm*tspd
  
    if (Als(j)+Acs(j) .le. 0.0) then
    
      fracL=0.5
    else
      fracL = Als(j) /(Als(j)+Acs(j))

    endif
    !write (*,*)"Als is", Als(j)
    !write (*,*)"Als is", Acs(j)
  
    !if (SS(j) .le. 0.15) then
    !  fracL=0.3
    !endif  
    if (fracL .eq. 1.0) then
      fracL=0.8
    endif
    if (fracL .le. 0.0) then
      fracL=0.2
    endif
    
  
    !write(*,*) "Monthly package is entered"
!!!!! Disturbance and retreat meaning
!This might be done in the next part

    !bsum(j)                        =(Bl(j) +sum(Br(j,:)))
    !M0(j)                          = M0(j) + kMrt(j)0* bsum(j)
   ! write (*,*) "The Mortality check"
   ! write (*,*) "Bsum here is", bsum(j)
   ! write (*,*) "M0 is", M0(j)
   ! write (*,*) "MD0 is ", MD0(j)
   ! write (*,*) "nfd is", nfd(j)
    Mrt(j)  = min(bsum(j),max(M0(j) *p_dt, bsum(j) * MD0(j)/max(p_critD,nfd(j))))
    !Mrt(j)  = min(bsum(j),max(M0(j), bsum(j) * MD0(j)/max(1.0,nfd(j))))
    !if (j .eq. 9) then
   ! write (*,*) "First step Mortal", Mrt(j)
      !write (*,*) "Biomass Sum", bsum(j)
!    write (*,*) "M0(j) is", M0(j)
!    write (*,*) "MD0(j) is", MD0(j)
!    write (*,*) "NFD is", nfd(j)
      !write (*,*) "The species number is ", j
    !write (*,*) "The productivity cumulative NPP is", npp0(j)
      !write (*,*) "The Tile_fraction is :", frac_tile(i,t)
    !endif
    do l =1,nsoil
      if (l .eq. 1) then
        Bin =npp0(j) * (1-fracL)
      else
        Bin = Bout
      endif
      Br(j,l) =Br(j,l) + fracR0*Bin - Mrt(j)* Br(j,l)/max(p_critD, bsum(j))
      RtoM(j) =RtoM(j) + max(0.0, Br(j,l) - (Brmax*1000) )
      Total_mortality_root(j)=Total_mortality_root(j)+Mrt(j)*Br(j,l)/max(p_critD,bsum(j))+RtoM(j)
      Bout=(1-fracR0)*Bin
  !!write (*,*) "Frac root layer wise", fracL
    enddo
    !write (*,*) "RtoM is", RtoM(J)
    Mrt(j) = Mrt(j) +RtoM(j)
    
    Rtb(j)	=sum(Br(j,:))    
    if (Rtb(j) .eq. 0.0) then
      Rtres(j)=Rtres(j)+0
    else
      do l =1,nsoil

        Rtres_l(j) = Rspec(j) * o_Q10_resp(j)**((xT_s(j) - o_ToptP(j)) / 10.0 )*c_MCo2*p_dt*Br(j,l)/Rtb(j)  !*c_MCo2
        Br(j,l) = Br(j,l)- Rtres_l(j)
        Rtres(j)=Rtres(j)+Rtres_l(j)
        !!write (*,*) "Root respiration amount", Rtres(j)
        !!write (*,*) "Root in layer", Br(j,l)
      enddo
    endif
    
    
    Bl(j)= Bl(j)+ npp0(j)*fracL- Mrt(j)*Bl(j)/max(p_critD, bsum(j))
    Rtb(j) = sum(Br(j,:))
    !write (*,*) "Fraction Leaf", fracL
    !write (*,*) "Mortal", Mrt(j)
    
    !write (*,*) "The leaf biomass", Bl(j)
    !write (*,*) "Root biomass", sum(Br(j,:))
    !write (*,*) "Root Mortality Total monthly", Total_mortality_root(j)
    !write (*,*) "Root respiration Total monthly", Rtres(j)
      !write (*,*) "Total Mortality for the timestep", Mrt_tot(j)
     
    Mrt_tot(j) = Mrt(j)*Bl(j)/max(p_critD, bsum(j))+Rtres(j)+Total_mortality_root(j)     !total dead matter in soil
   ! if (j .eq. 40) then
     ! write (*,*) "Fraction Leaf", fracL
      !write (*,*) "Mortal", Mrt(j)
    
      !write (*,*) "The leaf biomass", Bl(j)
      !write (*,*) "Root biomass", sum(Br(j,:))
      !write (*,*) "Root Mortality Total monthly", Total_mortality_root(j)
      !write (*,*) "Root respiration Total monthly", Rtres(j)
      !write (*,*) "Total Mortality for the timestep", Mrt_tot(j)
      
   ! endif
   ! write (*,*), "MON_CON IS", mon_cond(j)
   ! write (*,*), "Timer", counttimer(j)
    if (klife(i,t,j) .eq. 1 )then

      M_cond(j)=mon_cond(j)/counttimer(j)
      res_por(j)=por*(1-M_cond(j))
      if (res_por(j) .le. 0) then
        res_por(j)=0.1
      endif
    endif
    hmon=24*30.5
    if (klife(i,t,j) .eq. 1.0) then
       
     ! write(*,*) "The respor is ", res_por(j) 
     ! write (*,*) "Total Mortality for the timestep", Mrt_tot(j)
   ! write (*,*) "The saturation of the layers for the timestep is",Layer_con(j)
   ! write (*,*) "The remaining porosity is", res_por(j)
     ! write (*,*) "The remaining CO2 in soil from previous month", CO2_pre(j)
     ! write (*,*) "Atmospheric CO2",rCO2g_a
     ! write (*,*) "DPM IS", dpm
     ! write (*,*) "TSPD", tspd
   ! write (*,*) "Root Mortality Total monthly", Total_mortality_root(j)
   ! write (*,*) "Root respiration Total monthly", Rtres(j)
       
   ! write (*,*) "The Current co2 amount in soil", CO2_sink(j)

      CO2_sink(j)=max(rCO2g_a, &
                (((Mrt_tot(j)*1000)*0.15/(0.66*res_por(j)*0.05*dpm*tspd)) &
                *22.71108/44)+((CO2_pre(j)*1.94)*0.15/(0.66*res_por(j)*0.05*dpm*tspd)*22.71108/44) + rCO2g_a)
      
    else
      CO2_sink(j)= max(rCO2g_a, CO2_sink(j)+0.0 )
    endif
   ! write (*,*) "Total Mortality for the timestep", Mrt_tot(j)
   ! write (*,*) "The saturation of the layers for the timestep is",Layer_con(j)
   ! write (*,*) "The remaining porosity is", res_por(j)
   ! write (*,*) "The remaining CO2 in soil from previous month", CO2_pre(j)
   ! write (*,*) "Atmospheric CO2",rCO2g_a
   ! write (*,*) "The leaf biomass", Bl(j)
   ! write (*,*) "Root biomass", sum(Br(j,:))
   ! write (*,*) "Root Mortality Total monthly", Total_mortality_root(j)
   ! write (*,*) "Root respiration Total monthly", Rtres(j)
       
   ! write (*,*) "The Current co2 amount in soil", CO2_sink(j)
    
      
    
    CO2_pre(j)=max(rCO2g_a,  CO2_sink(j)) 
    !if (CO2_pre(j) .ge. 100000) then
    !write (*,*) "Location", i
    !write (*,*) "Species", j
    !write (*,*) "The Co2 soil sink is", CO2_pre(j)
    !endif
    
    Rtres(j)=0.0
    
    Als(j)=0.0
    Acs(j)=0.0
    
    N_leaf(j)=0.7*Bl(j)/o_w_leaf(j)
    htree(j)=0.3*Bl(j)/(p_rhoBl*1.5)/(pi*rstem*rstem)        !height of tree, density of tree 1.5 times that of leaf
    Lai_new(j)=o_A_leaf(j)*N_leaf(j)/o_G_area(j)
    kcccc=kcccc + 1
    !write (*,*) "Counting", kcccc 
    !write (*,*) "THe species number ending here is:", kcccc
    !if (Lai_new(j) .le. 0.0) Lai_new(j)= 0.0
    !endif
    nfd(j)=0.0
    MD0(j)=0.0
    M0(j)=0.0
    Mrt(j)=0.0
    RtoM(j)=0.0
    Total_mortality_root(j)=0.0
    Run_tot(i,t,j)=Run_tot(i,t,j) + max(0.0,fH2Ol_xd(j))
    mon_cond(j)=0.0
    counttimer(j)=0
  enddo

  ngsum=0.0
  !do j=1,p_nspec
  !  write (*,*) "Species number:", j
  !  write (*,*) "Lai of tree", Lai_new(j)
  !  !write (*,*) "Runoff for the species_ last timestep", fH2Ol_xd(j)
  !  !write (*,*) "Runoff for the species_ last timestep", Run_tot(i,t,j)
  !  write (*,*) "Layer water content", Layer_con(j)
  !  write (*,*) "Bucket water content", S2(j)
  !  !if (j .eq. 40) then
  !  !write (*,*) "Carbon Dioxide in the soil due to the species ", CO2_sink(j)
  !  !endif
  !  
  !  if (klife(i,t,j) .eq. 1.0) ngsum=ngsum +max(0.0, netgrowth(i,t,j))
  
!    if (CO2_pre(j) .ge. 50000) then
!      write (*,*) "Species number exceeding CO2 sink is", j
!      write (*,*) "The Co2 soil sink is", CO2_pre(j)
!    endif
  !enddo
! Start loop over all species - IIb -

!  if (hsum .gt. p_critD) then

  if (ngsum .gt. p_critD) then

    wgtsum                      = 0.0

    do j = 1,p_nspec
      
      if (klife(i,t,j) .eq. 1.0) then

!        wgtspec(j)              = 1.0 &                                 ! equal weights at low cover (no competition)
!                                * ( (o_zt(j)/hsum)**csum )**1           ! competition under low disturbance -> weighting by growth height
!        wgtspec(j)              = 1.0 &                                 ! equal weights at low cover (no competition)
!                                * ( ( max(0.0,netgrowth(i,t,j)) &                                 
!                                / ngsum )**(Lai_new(j)/Lai_cum))**2 

        wgtspec(j)              = 1.0 &                                 ! equal weights at low cover (no competition)
                                * ( ( max(0.0,netgrowth(i,t,j)) &                                 
                                / ngsum )**csum )**2                    ! competition under low disturbance -> weighting by net growth    !May be keep it or use "LAI"

!        wgtspec(j)              = 1.0                                   ! equal weights (neutral model)

        wgtsum                  = wgtsum + wgtspec(j)
      endif
    enddo ! End of loop over all species - IIb -
     
! Start loop over all species - IIc -
 
    do j = 1,p_nspec

      if (klife(i,t,j) .eq. 1.0) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!**********!!!!!!!!!!!!!!***********
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!******
!        expansion(j)       = max(0.0, min( &                       ! [m2 T / m2 V / month]
!                                  netgrowth(i,t,j) *area_s(i,t,j), & ! mass balance constraint
!                                  (1.0 - (Lai_new(j)/Lai_cum)) &                        ! available area
!                                * p_NCbt &                              ! establishment
!                                * wgtspec(j) /wgtsum) )                 ! competition                                          !!must decide


        expansion(j)       = max(0.0, min( &                       ! [m2 T / m2 V / month]
                                  netgrowth(i,t,j) *area_s(i,t,j), & ! mass balance constraint
                                  (1.0 - csum) &                        ! available area
                                * p_NCbt &                              ! establishment
                                * wgtspec(j) /wgtsum) )                 ! competition                                          !!must decide

      else
        expansion(j)       = 0.0
      endif
    enddo
  else
    expansion(:)           =0.0
  endif
  expsum                        = sum(expansion(:))
  if (expsum .gt. p_critD) then
    !must decide
    expansion(:)                = expansion(:) * min(1.0, &       !!!!must decide
                                  (1.0-csum) / expsum)
!   expansion(:)                = expansion(:) * min(1.0, &       !!!!must decide
!                                 (1.0-(Lai_new(j)/Lai_cum)) / expsum)

  else
    expansion(:)                = 0.0
  endif

    
!-----------------------------------------------------------------------
! Start loop over all species - III -
!-----------------------------------------------------------------------


  do j = 1,p_nspec


    if (klife(i,t,j) .eq. 1.0) then

! *** Cover change II ***

! *** Cover change II ***

      !retreat                   = min( max(0.0, - netgrowth(i,t,j)) &   ! [m2 T / m2 V / month]
       !                         * area_s(i,t,j), &
        !                          area_s(i,t,j) )
                            
!      disturbance               = 1.0/tauD(i,t,v,h) * area_s(i,t,j)  ! [m2 T / m2 V / month]  This one may be used later
                              
      area_s(i,t,j)       = max(0.0, area_s(i,t,j) &           ! [m2 T / m2 V]
                                + expansion(j))! &
                                !- retreat )                  !- disturbance &  retreat removed 
                                
                            
      netgrowth(i,t,j)      = 0.0
      
      if (area_s(i,t,j) .le. frac_s_crit) then

        klife(i,t,j)        = 0.0
        area_s(i,t,j)     = 0.0

      endif 

! Set lichen to dead if cover is too low

    endif ! check for survival

  enddo ! End of loop over all species - III -
endif ! End of execute once per month

!write (*,*) "LYcophyte_nova_crucial calculátions in lycophytes: is close to end"


! *** Average strategies ***

as_lai_s                        = Lai_cum !new var

as_area_s                       = csum ! total area of all strategies

!if (as_area_s .gt. max_area_grid)then
!  as_area_s=  0.4
!else
!  as_area_s= csum/max_area_grid
!endif
if (writeout) then

  ! initialise accumulated variables with zero
  as_rCO2d                      = 0.0
  as_sCO2d                      = 0.0
  as_rCb                        = 0.0
  as_Lai                        = 0.0
  as_fCO2gc                     = 0.0
  as_fCcg                       = 0.0
  as_fCcb                       = 0.0
  as_fCcb_l                     = 0.0
  as_fCcb_c                     = 0.0
  as_fCbo                       = 0.0
  as_fH2Ogl_ut                  = 0.0
  as_fH2Olg_tu                  = 0.0
  
  as_fH2Ol_lsat                 = 0.0 !new var
  as_fH2Ol_bsat                 = 0.0 !new var
  as_fH2Ol_runoff_l             = 0.0
  as_fCc_npp             = 0.0
  as_fCc_gpp             = 0.0
  as_Ts                         = 0.0
  as_Tg                         = 0.0
  as_H                          = 0.0
  as_G                          = 0.0
  as_E                          = 0.0
  as_C                          = 0.0
  as_EB                         = 0.0
endif
! loop over all strategies
do j = 1,p_nspec

  as_npp(j)=npp0(j)
  as_gpp(j)=gpp0(j)
    

  ! Check if lichen is alive
  if (day .eq. dpm .and. ts .eq. tspd) then
    
    if (klife(i,t,j) .eq. 1.0) then
      
      as_fH2Ol_td(j)                 = Run_tot(i,t,j)!fH2Ol_xd(j) !as_fH2Ol_td + fH2Ol_xd(j) *1000.0*c_year *area_s(i,t,j)
    !as_fH2Ol_td(j)                 = as_fH2Ol_td + fH2Ol_xd(j) *area_s(i,t,j)
    !as_fH2Ol_runoff_l= as_fH2Ol_td(j)
    else
      as_fH2Ol_td(j)                 =0.0!!newly added
    !as_fH2Ol_runoff_l=0.0

    endif
    fH2Ol_xd(j)=0.0
    npp0(j)=0.0
    gpp0(j)=0.0
 
  endif
enddo

!write (*,*) "LYcophyte_nova_crucial calculátions in lycophytes: is ending and the file writing will be started"

!-----------------------------------------------------------------------
! Start loop over all species - V -
!-----------------------------------------------------------------------

if (writeout) then
 write (*,*) "enters writeout section"
  tim =tim +1 
  do j = 1,p_nspec
   
    ! Check if lichen is alive
   
    if (klife(i,t,j) .eq. 1.0) then
      if(tim == tstep_total)then

        write (*,*) "The area sum is", sum(area_s(i,t,:))
        write (*,*) "Species number:", j
        write (*,*) "Lai of tree", Lai_new(j)
    !write (*,*) "Runoff for the species_ last timestep", fH2Ol_xd(j)
    !write (*,*) "Runoff for the species_ last timestep", Run_tot(i,t,j)
        write (*,*) "Layer water content", Layer_con(j)
        write (*,*) "Bucket water content", S2(j)
        write (*,*) "The remaining CO2 in soil from previous month", CO2_pre(j)
      endif

      ! intensive variables are weighted by cover   !STILL I USE AREA BUT LATER MAY BE CHANGE TO lai
      !if (csum .le. 0.0) csum = 1.0
      !if (Lai_cum .le. 0.0) Lai_cum = 1.0  
      cweight                   = area_s(i,t,j) / csum

      lweight                   = Lai_new(j) / Lai_cum

      as_rCO2d                  = as_rCO2d + Mrt_tot(j) *area_s(i,t,j)

      as_sCO2d                  = as_sCO2d + CO2_pre(j)*cweight!area_s(i,t,j)
      as_Lai                    = as_Lai + Lai_new(j) * cweight!area_s(i,t,j)
      as_rCb                    = as_rCb + o_spec_area(j) * area_s(i,t,j)

      as_fH2Ol_lsat             = as_fH2Ol_lsat  + Layer_con(j) *area_s(i,t,j) !new var         **********WEIGHTING SCHEME***THINK ABOUT IT
 
      as_fH2Ol_bsat             = as_fH2Ol_bsat  + S2(j) *area_s(i,t,j) !new var                **********WEIGHTING SCHEME***THINK ABOUT IT
      
      as_fH2Ol_runoff_l         = as_fH2Ol_runoff_l + as_fH2Ol_td(j)*area_s(i,t,j)
      
      as_fCc_npp                = as_fCc_npp + fCO2nc(j)*area_s(i,t,j)
      
      as_fCc_gpp                = as_fCc_gpp + as_gpp(j) *area_s(i,t,j)
                                    
      as_fCO2gc                 = as_fCO2gc + fCO2gc(j) *area_s(i,t,j)
                                
      as_fCcg                   = as_fCcg + (Rspec(j)) *c_MCo2*p_dt *area_s(i,t,j)
                                
      as_fCcb                   = as_fCcb + fCcb(j) *area_s(i,t,j)

      as_fCcb_l                 = as_fCcb_l + (fCO2gc_L(j)-fCcg_M(j))  *area_s(i,t,j)
                                
      as_fCcb_c                 = as_fCcb_c + (fCO2gc_W(j)-fCcg_M(j))  *area_s(i,t,j)

      as_fCbo                   = as_fCbo + fCbo(j) *area_s(i,t,j)

      as_fH2Ogl_ut              = as_fH2Ogl_ut + min(0.0,fH2Olg_ga(j))  *area_s(i,t,j)
                                
      as_fH2Olg_tu              = as_fH2Olg_tu + max(0.0,fH2Olg_ga(j))  *area_s(i,t,j)

      as_Ts                     = as_Ts + xT_s(j) *cweight
   
      as_Tg                     = as_Tg + xT_g(i,t,j) *lground *cweight
   
      as_H                      = as_H + fRAD_H(j) *area_s(i,t,j)
   
      as_G                      = as_G + fQ_tg(j) *area_s(i,t,j)
   
      as_E                      = as_E + fQ_ta_L(j) *area_s(i,t,j)
   
      as_C                      = as_C + fQ_ta_S(j) *area_s(i,t,j)
   
      as_EB                     = as_EB + (fRAD_H(j)-fQ_ta_L(j)-fQ_ta_S(j)) *area_s(i,t,j)
      
   
      ! average BSC-related properties
   
      !     if (i2 .eq. 1 .and. BSCtypes) call lycom_accBSC(i,i2,j)
   
      endif ! check for survival
   
  enddo ! End of loop over all species - VI -
  Lai_cum = 0.0
  csum = 0.0
  write(*,*) "writing complete"
endif ! write output?

return
end subroutine nova_step

end module lycom_nova

