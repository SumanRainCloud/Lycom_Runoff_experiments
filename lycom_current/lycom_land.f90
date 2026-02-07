module lycom_land
contains



subroutine first_land()


use lycom_par
use lycom_opt

implicit none

integer :: i,i2,j,m,t,x,y
integer :: nCPts3             !!number of sites 
real    :: taC
real    :: Rnet
integer    :: nyrsim
integer    :: hmon, hyear
integer    :: hdata
real    :: ETpot0
real    :: desatdT3
real    :: ETd
integer, dimension(7) :: stat3
integer :: stat
real, dimension(24*31*12*40) :: p, q, r
! open climate forcing files


!write (*,*) "Preprocessing started!"

!nyrsim=40
!hmon=24*30.5
!hyear=hmon*12
!hdata=(hyear*nyrsim)

!desatdT3=0.0
!ETd=0
!j=0
!if (rank .eq. 1) then
! opening the file for reading
!open (1, file = 'srad', status = 'old')
!open (2, file = 'lrad', status = 'old')
!open (3, file = 'tair', status = 'old')
!do i2 = 1,hdata  
!  read(1,*) p(i2)
!end do 
   

!do i2 = 1,hdata  
!  read(2,*) q(i2)
!end do 
   

!do i2 = 1,hdata  
!  read(3,*) r(i2)
!end do 
   

   
   

!do y =1,hdata

!  Rnet = p(y)*(1.0-p_albVEG) + p_eps*q(y) - p_eps*c_sigma*r(y)*r(y)*r(y)*r(y)

!  taC = r(y) - c_TH2Osl

!  if (taC .gt. 0.0) then

!    desatdT3 = exp(p_esatAIR1 *taC /(p_esatAIR2 + taC)) *p_esatAIR3 *p_esatAIR1 * p_esatAIR2 /(p_esatAIR2 +taC)**2

!    ETpot0 = 1.3 * Rnet * desatdT3 / ( desatdT3 + c_gamma) / c_HH2Olg / c_rhoH2Ol
!  else

!    ETpot0 = 0.0
!  endif

!  ETd = ETd + ETpot0
!enddo

!ETdavg2 = ETd / float(hdata)
!write (*,*) "The ETd:preprocessing ", ETd
!write (*,*) "The hdata:preprocessing ", hdata
!write (*,*) "The ETdavg:preprocessing ", ETdavg

!do y = 1,hdata

!  Rnet = p(y)*(1.0-p_albVEG) + p_eps*q(y) - p_eps*c_sigma*r(y)*r(y)*r(y)*r(y)

!  taC = r(y) - c_TH2Osl

!  if (taC .gt. 0.0) then

!    desatdT3 = exp(p_esatAIR1 *taC /(p_esatAIR2 + taC)) *p_esatAIR3 *p_esatAIR1 * p_esatAIR2 /(p_esatAIR2 +taC)**2

!    ETpot0 = 1.3 * Rnet * desatdT3 / ( desatdT3 + c_gamma) / c_HH2Olg / c_rhoH2Ol
!  else

!    ETpot0 = 0.0
!  endif
!  ETratio = ETpot0 / ETdavg(i)

!  ETrmax2 = max(ETrmax(i), ETratio)
  
!enddo
!nCPts3=1
!do y = 1,nCPts3
!  Etr_avg(y)=ETdavg2
!  Etr_max(y)=ETrmax2
!enddo
!write (*,*) "The ETrmax:preprocessing ", ETrmax
!write (*,*) "The ETdavg:preprocessing ", ETdavg
!write (*,*) "first_land:preprocessing module is finished"
!else
!endif
  
  
!close(1)
!close(2)
!close(3)

return
end subroutine first_land





!***********************************************************************
! LAND_INIT
!***********************************************************************

subroutine land_init (nCPts3)
use lycom_par
use lycom_opt

implicit none

integer :: i,t,l, nCPts3  
! Initialize soil parameters FIRST
  ! Soil porosity (example value)
        ! Maximum water content per layer [m]
        ! Maximum bulk soil water content [m]
       ! Percolation rate parameter
       ! Baseflow rate parameter
          ! Number of soil layers



! Check for valid initialization





! initialise atmospheric conditions

rCO2g_a                         = 360.0                                 ! [ppm]
rO2g_a                          = 210000.0                              ! [ppm]


! initialise fields

rH2Os_g(:)                      = 0.0                                   ! snow layer [m3 H2O / m2 G]

rH2Ol_0(:,:)                = 0.0                                   ! water on leaves [m3 H2O / m2 C]

rH2Ol_g1(:,:)                   = 0.0                                   ! water in the top soil [m3 H2O / m2 C]



!!!!!watch out for the variable below                                
rH2Ol_g2(:,:)                   = 0.0                                   ! water in the bulk soil [m3 H2O / m2 C]
!!!!!Watch out for the variable below too
do i = 1, nCPts3
  do t = 1, p_ntiles
    Wx(i,t)    = 0.5 * por * 0.65
    if (Wx(i,t) .le. 0.0 .or. Wx(i,t) .ne. Wx(i,t)) then
      write(*,*) "WARNING: Invalid Wx at rank=", rank, " i=", i, " t=", t
      Wx(i,t) = 0.15  ! Fallback value
    endif

    xT_g0(i,t) = 288.0

    do l = 1,nsoil
      W_c0(i,t,l) = 0.5 * por * 0.03
      if (W_c0(i,t,l) .le. 0.0 .or. W_c0(i,t,l) .ne. W_c0(i,t,l)) then
        write(*,*) "WARNING: Invalid W_c0 at rank=", rank, " i=", i, " t=", t, " l=", l
        W_c0(i,t,l) = 0.015  ! Fallback value
      endif
    enddo
  enddo
enddo

S_2        = 0.0
layer_con0 = 0.0
bucket_con0 = 0.0
W_con_ll   = 0.0

Q_Oflow    = 0.0
Q_Oflow2   = 0.0
Q_per      = 0.0
Q_Oflow_b  = 0.0
Q_Oflow_b2 = 0.0

fH2Ol_gwl  = 0.0
fH2Ol_xd_land = 0.0
Runoff_land   = 0.0


!atN_fH2Ol_xd(:,:,:)             = 0.0                                   ! water overflow from layers [m/s]    or may be used as ! water input from upper layers [m/s]

!atN_fH2Ol_bd(:,:)               = 0.0                                   ! overflow from bark [m/s]

!fH2Ol_ug2(:,:)                  = 0.0                                   ! throughfall [m/s]

! initialise Soil & air properties

do i=1,nCPts3

    

  ! heat conduction
  if ( (biome(i) .le. 6.0) .or. (biome(i) .eq. 12.0) ) then

    CSOIL(i)                    = p_CSOIL_for
    kSOIL(i)                    = p_kSOIL_for
  
    dewmax(i)                   = -p_dewmax_for / 1000.0 / 365.25 / 12.0 / 3600.0                       ! transform value mm/yr to m/s, assuming 12h night
                                
  else
    CSOIL(i)                    = p_CSOIL_des
    kSOIL(i)                    = p_kSOIL_des
  
    dewmax(i)                   = -p_dewmax_des / 1000.0 / 365.25 / 12.0 / 3600.0
                                
  endif
enddo

return
end subroutine land_init








!***********************************************************************
! LAND_STEP G
!***********************************************************************
subroutine land_stepG(i)
use lycom_par
use lycom_opt
implicit none
integer :: i

! day and time step counter
if (fRADs_ad(i) .gt. 0.0) then
  kday                          = 1
  naccu_d(i)                    = naccu_d(i) + 1
else
  kday                          = 0
  naccu_n(i)                    = naccu_n(i) + 1
endif

naccu(i)                        = naccu(i) + 1





! saturation vapor pressure

Ta3                             = xT_a(i) * xT_a(i) * xT_a(i)

Ta4                             = Ta3 * xT_a(i)

zT_a                            = xT_a(i) - c_TH2Osl                    ! air temperature [deg C]
                 
esatAIR                         = p_esatAIR3 *exp(p_esatAIR1*zT_a / (p_esatAIR2 + zT_a))     ! saturation vapour pressure [kg * m / (s2 * m2) = Pa]                  REF: Allen,1998
                                

desatdT                         = exp(p_esatAIR1*zT_a/(p_esatAIR2 + zT_a)) * (p_esatAIR1 *p_esatAIR2 *p_esatAIR3 / ((p_esatAIR2 + zT_a)* (p_esatAIR2 + zT_a)))     ! slope of saturation vapour pressure vs temperature relationship []
                                
                                 
                                

! Snow balance

fH2Osl_g                        = min(3.22 * max(0.0,xT_a(i)-c_TH2Osl)/ c_day / 1000.0,rH2Os_g(i) / p_dt + fH2Os_ad(i)/1000) ! snow melt [m/s]                                                       REF: Bergstrom,1992
                                 
                                  

rH2Os_g(i)                      = max(0.0, rH2Os_g(i) + fH2Os_ad(i)/1000 * p_dt - fH2Osl_g * p_dt- rH2Os_g(i) * p_H2Os_loss * p_dt)
                                 
                                


dsnow                           = rH2Os_g(i) * c_rhoH2Ol / p_rhoH2Os    ! snow depth [m]

return
end subroutine land_stepG

!***********************************************************************
! LAND_STEP V
!***********************************************************************


subroutine land_stepV (i,t,v)
use lycom_par
use lycom_opt
implicit none

integer :: i,t,v,m,n,o
real :: fracRADs_0
! initial ground temperature
real    :: rain

! Switches

if (v .eq. 2) then ! canopy
  lground                       = 0.0  !at canopy(not req for my case)
else
  lground                       = 1.0   !at groudlevel (soil)
endif


! Radiation

Aslai                           =  0.0  !areaLEAF_month(i)   !for lycom I have this as surrounding forest Lai Needs to be input according to region ***Need to make a map and read it from there

!*PLEASE CECK THIS PART*!

if (Aslai .le. p_critD) Aslai = 0.0
!!Be it forest lycophytes grow only on ground    THats the reason for the commented next section
!if (t .eq. 1) then ! forest   

!  !if (v .eq. 2) then ! canopy

!  fracRADs_0                  = 1.0 - exp(-p_beer_s * Aslai)          ! fraction of total absorbed short-wave radiation in canopy             REF: Bonan,2008

!  fracRADs                    = exp(-p_beer_s * Aslai )      ! fraction of average absorbed short-wave radiation in canopy           REF: Bonan,2008

!  fracRADl                    = 1.0 - exp(-p_beer_l * Aslai)          ! fraction of total absorbed long-wave radiation in canopy              REF: Bonan,2008

!!not req here but where lycophytes grow  fracRADs 			  =  (1-fracRADs_0) * (1.0 - exp(-p_beer_s*Lai_new)) * (1.0-p_albVEG);    !fraction of total absorbed long-wave radiation by lycophytes

!!not req here but where lycophytes grow  fracRADl 			  =  (1-fracRADl_c) * (1.0 - exp(-p_beer_l * Lai_new));                   !fraction of total absorbed long-wave radiation by lycophytes              REF: Bonan,2008


    
!else ! no forest

fracRADs                      = exp(-p_beer_s * Aslai)                ! fraction of total absorbed short-wave radiation on ground             REF: Bonan,2008

fracRADl                      = exp(-p_beer_l * Aslai)                ! fraction of total absorbed long-wave radiation on ground              REF: Bonan,2008
!endif


fH2Ol_tb_f   =fH2Ol_tbf  !!!!needs to be initialisedglobal variable!!!!!!!!!!!


!water fluxes

!the if statement is changed sine there is no need to have two part such as canopy and ground separately

rain =fH2Ol_ad(i)!1.5E-06 !fH2Ol_ad(i)!*0.001 !(ESTONIA/SA/IND)
!write (*,*) "LYCOMLAND prints here"
!if (rain .gt. 0.0) then
!  write (*,*) "rain from Lycom land", rain, "for timestep", i 
!endif
if (t .eq. 1) then  !forest

  
  fracrain(i,t)                 = 2.0/p_LAImax

  !fH2Ol_ci(i,t)                 = max(0.0,rain*fracrain(i,t)/max(1.0,Acano(i,t))*0.35 *p_dt)!mul 0.65        ! water input into canopy [m3 H2O / (m2 C * s)]     !!!!!! THIS CANOPY WATER--may be utilised for direct evaporation
  fH2Ol_ci(i,t)                 = max(0.0,rain*fracrain(i,t)*0.35 * p_dt)
  fH2Ol_ts(i,t)                =  max(0.0,min(rain * p_dt-fH2Ol_ci(i,t),fH2Ol_tb_f))!!!0.0 
  
  

  fH2Ol_ux(i,t)                  = max(0.0,rain * p_dt- fH2Ol_ci(i,t) - fH2Ol_ts(i,t)) ! water input into soil as throughfall after the loss in the topsoil [ m3 H2O / m2 G / s ]

  if (fH2Ol_ts(i,t)-p_rmaxH2Ol_g1 .eq. 0.0) then 
    fH2Ol_tb_f                =0.0 !if filled
  else
      
    fH2Ol_tb_f                = (p_rmaxH2Ol_g1-fH2Ol_ts(i,t))  !if not filled

  endif                      
  
  
else
  
  fH2Ol_ci(i,t)                 =0
  fH2Ol_ts(i,t)                 =min(rain * p_dt,fH2Ol_tb_f)!!!0.0!commented       !0.05*rain    water on top of the soil('may be 5% of the rainfall')      !!!!!! THIS Water on top of the soil--may be utilised for direct evaporation


  fH2Ol_ux(i,t)                 =rain * p_dt- fH2Ol_ts(i,t)    !water into soil

  if (fH2Ol_ts(i,t) -p_rmaxH2Ol_g1 .le. 0.0) then 
    fH2Ol_tb_f                =0.0 !if filled
  else
      
    fH2Ol_tb_f                = (p_rmaxH2Ol_g1-fH2Ol_ts(i,t))  !if not filled

  endif                      
  
  fH2Ol_ux                 =rain-fH2Ol_ts(i,t)  !commented   !water into soil   'This becomes a bit more complicated in case of bareground and snow or frozen ground'
endif


!fH2Ol_tbf will increase with evaporation


fH2Ol_ux(i,t)= fH2Ol_ux(i,t) + fH2Osl_g * p_dt  !add snow melt

!do m=1,i
!  do n=1,t
!    Wx(m,n) =  0.5*por*0.65
!    xT_g0(m,n)                = 288.0 !xT_a(i)                               ! bare ground temperature [K] 	*changed accorsing to lycom
!    do o=1,nsoil
!      W_c0(m,n,o)=0.5*por*0.03
!    enddo
!  enddo
!enddo

return 
end subroutine land_stepV

!***********************************************************************
! LAND_STEP     !HERE need to incorporate the soil part as well! interactions
!***********************************************************************

subroutine land_step (i,t,v,h)
use lycom_par
use lycom_opt
implicit none

!integer :: i,t,v,h,j,m,l
integer :: i,t,v,h,j,m,l

real    :: kH2Og0, Rnet
real    :: grh0, cdg0, xT_s_wet0, xT_s_dry0, crh0
real    :: wetfrac_0
real    :: dRAD, fRAD_Hw0, fRAD_Hd0, ETpot0
real    :: ETpot_v, soilevap, Evap


!debug
!real    :: wso1, wso2


!'fH2Ol_ci' & 'fH2Ol_ts' canopy and top soil water both goes into surface runoff
if (t .le. 2) then !only forest and grassland

  kH2Og0                 =p_vonKarman * p_vonKarman * max(p_critD,fAIR_s(i))/ dn_kH2Og(i,t)           ! [m / s]                                                               REF: Allen,1998
                                 
                                
  

  if ( dsnow .ge. p_H2Os_crit) then ! snow layer at ground too thick
  !in this case 'fH2Ol_ci' & 'fH2Ol_ts' canopy and top soil water both goes into surface runoff

    
    ETpot0                      =0.0
    fH2Olg_xu0                  =0.0  !upward water flux from vascular(evaporation from canopy in my case)
    fH2Ol_xd0                   =0.0  !!!or 0.0 becomes complicated   !'new' !Downward water flux from vascular (excess water after evaporation in my case)
    !THE above one gets complicated in the frozen state
    fH2Ol_go                    =fH2Ol_ci(i,t)         !This is not req (+fH2Ol_ts) top soil water remains same   !surface runoff
    fH2Ol_bd0                   =0.0                 !stem outflow from bark (not utilised in my case)
    xT_s_dry0                   = min(c_TH2Osl - 0.1, xT_a(i))
    xT_s0                       = xT_s_dry0
    fRAD_H0                     = 0.0
    fQ_ta_L0                    = 0.0
    fQ_ta_S0                    = 0.0

  else

    grh0                        = c_gamma / kH2Og0

    cdg0                        = c_CAIR*(desatdT+c_gamma)
                               
    xT_s_wet0                   = ( grh0 * (fRADs_ad(i)*fracRADs &      ! surface temperature [K]                                               REF: Monteith,1981
                                * (1.0-albLsf(i,t)) &
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL*xT_g0(i,t) * lground ) &
                                + xT_a(i)*cdg0 &
                                - c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) ) &
                                / ( grh0*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) + cdg0 )

    
    crh0                        = 1.0/(c_CAIR*kH2Og0)
                               
    xT_s_dry0                   = ( xT_a(i) &
                                + crh0*(fracRADs*(1.0-albLsf(i,t))*fRADs_ad(i) &
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) &
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL*xT_g0(i,t) * lground) ) &
                                / ( 1.0 + crh0*(4.0*p_eps*c_sigma*Ta3 * (2.0-lground) &
                                + kSOIL(i)/p_dz_SOIL * lground) )

    dRAD                        = fRADs_ad(i)*fracRADs*(1.0-albLsf(i,t))&  ! downwelling radiation [W / m2 T]
                                + fracRADl * p_eps * fRADl_ad(i) &
                                + fracRADl * p_eps * c_sigma*Ta4 * (1.0-lground) & ! lrad from below
                                + (1.0-fracRADl)*p_eps*c_sigma*Ta4 * lground &
                                + 3.0*p_eps*c_sigma*Ta4 * (2.0-lground)

    fRAD_Hd0                    = dRAD &                                ! net radiation (dry) [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_dry0 * (2.0-lground) &
                                - kSOIL(i) *(xT_s_dry0 - xT_g0(i,t)) /p_dz_SOIL * lground

    if (xT_s_dry0 .lt. c_TH2Osl) then ! frozen surface

      fRAD_Hw0                  = 0.0
      ETpot0                    = 0.0
      fH2Olg_xu0                = 0.0
      fH2Ol_xd0                 = 0.0
      fH2Ol_go                  =fH2Ol_ci(i,t) !+ fH2Ol_ts(i,t)
      !fH2Ol_bd0                 = 0.0
      wetfrac_0                 = 0.0
    else


      fRAD_Hw0                  = dRAD &                                ! net radiation (wet) [W / m2 T]
                                - 4.0*p_eps*c_sigma*Ta3*xT_s_wet0 * (2.0-lground) &
                                - kSOIL(i) *(xT_s_wet0 - xT_g0(i,t)) /p_dz_SOIL * lground

      ETpot0                    = ( fRAD_Hw0 * desatdT &                ! potential evaporation [m3 H2O / (m2 C * s)]                           REF: Monteith,1965
                                + c_CAIR*(esatAIR-rH2Og_RH(i)*esatAIR) *kH2Og0 ) &
                                / (desatdT+c_gamma) / c_HH2Olg / c_rhoH2Ol

      ! evaporation
      if (ETpot0 .lt. 0.0) then
        fH2Ol_bd0                 = 0.0
        fH2Ol_xd0            = fH2Ol_ci(i,t)  
          
        fH2Olg_xu0            = 0.0
        wetfrac_0             = 1.0
      else
        fH2Ol_bd0                 = 0.0
        fH2Olg_xu0            = min(fH2Ol_ci(i,t),ETpot0 * p_dt)
        fH2Ol_xd0             = (fH2Ol_ci(i,t) - fH2Olg_xu0)
        wetfrac_0             = 0.0
      endif    
      !Water balance
      fH2Ol_go              =0.0


      fH2Ol_gwl    =max(0.0, fH2Ol_ux(i,t) + fH2Ol_xd0) !new variable
      
      

      
    
    endif
    
    xT_s0                      = wetfrac_0 * xT_s_wet0 + (1.0 - wetfrac_0)* xT_s_dry0

    fRAD_H0                    = wetfrac_0 * fRAD_Hw0 + (1.0- wetfrac_0)* fRAD_Hd0
    ! Latent heat
    fQ_ta_L0                    = fH2Olg_xu0 * c_HH2Olg * c_rhoH2Ol     ! Latent heat flux [W / m2 T]

   ! Sensible heat
    fQ_ta_S0                    = c_CAIR * (xT_s_wet0 - xT_a(i)) &      ! Sensible heat flux [W / m2 T]
                                * kH2Og0 * wetfrac_0 &
                                + c_CAIR * (xT_s_dry0 - xT_a(i)) &
                                * kH2Og0 * (1.0 - wetfrac_0)
  endif
  
else !wetland or rock

  ETpot0                        = 0.0
  fH2Olg_xu0                    = 0.0
  fH2Ol_xd0                     = 0.0
  fH2Ol_bd0                     = 0.0
  xT_s0                         = xT_a(i)
  fRAD_H0                       = 0.0




endif
! Potential evaporation

Rnet                        = fRADs_ad(i)*0.85 + p_eps*fRADl_ad(i) - p_eps*c_sigma*Ta4
               
ETpot_v                     = 1.4 *Rnet*desatdT/(desatdT+c_gamma) & ! [m3 H2O / (m2 G * s)]
                                / c_HH2Olg/c_rhoH2Ol

soilevap=  min(fH2Ol_ts(i,t), max(0.0,ETpot_v * p_dt))

fH2Ol_ts(i,t)=fH2Ol_ts(i,t)-soilevap

fH2Ol_tb_f =p_rmaxH2Ol_g1-fH2Ol_ts(i,t)

fH2Olg_ga_1= fH2Olg_xu0 +soilevap   !evaporation from ground and canopy

return
end subroutine land_step

!I want to do this separately since I donot want the water required by the lycophytes to be used up for the transpiration of the surrounding vasular vegetation This should be done after the water use up by lycophytes

!***********************************************************************
! LAND_STEP_VTranspiration     !HERE need to incorporate the soil part as well! interactions
!***********************************************************************
subroutine land_stepvTrans (i,t,v,h)
use lycom_par
use lycom_opt
implicit none

integer :: i,t,v,h,j,m,l
real    :: ETpot_v, Rnet, Evap
real    :: percolation, rootuptk, trans, Lay_upt

! soil water balance
Lay_upt = 0.0

if (xT_s0 .lt. c_TH2Osl) then ! inactivity below zero degrees

  fH2Ol_gb = max(p_critD, min(Qb0 * S_2 * p_dt, Wx(i,t)))
  
  Runoff(i,t) = fH2Ol_gb
  Wx(i,t) = Wx(i,t) - Runoff(i,t)
  
  fH2Ol_xd_land = Runoff(i,t)
  Runoff_land = fH2Ol_xd_land
  
  ! FIX: Add safety check
  if (Wxmax .gt. p_critD) then
    bucket_con0 = Wx(i,t) / Wxmax
  else
    bucket_con0 = 0.0
  endif
  
else
  
  Rnet = fRADs_ad(i)*0.85 + p_eps*fRADl_ad(i) - p_eps*c_sigma*Ta4
  
  ETpot_v = 1.4 * Rnet * desatdT / (desatdT + c_gamma) / c_HH2Olg / c_rhoH2Ol
  
  if (Wx(i,t) .gt. p_critD) then
    rootuptk = min( Wx(i,t)/ p_dt, p_kH2Ol_sv * (Wx(i,t)/p_rmaxH2Ol_g2)**2 )
  else
    rootuptk = 0.0
  endif
  
  trans = min( max(0.0, ETpot_v), rootuptk ) * p_dt
  
  ! Upper soil layer
  do l = 1, nsoil
    
    if (l .eq. 1) then
      Qin(i,t,l) = max( 0.0, fH2Ol_gwl )
    else
      Qin(i,t,l) = Q_Per + Q_Oflow + Q_Oflow2
    endif
    
    W_c0(i,t,l) = W_c0(i,t,l) + Qin(i,t,l)
    Q_Oflow = max(0.0, W_c0(i,t,l) - Wmax)
    
    W_c0(i,t,l) = W_c0(i,t,l) - Q_Oflow
    S_wc0 = max(0.0, min(1.0, W_c0(i,t,l) / Wmax))
    
    Q_per = max(0.0, min(Qp0 * S_wc0 * p_dt, W_c0(i,t,l)))
    
    W_c0(i,t,l) = W_c0(i,t,l) - Q_per
    
    Q_Oflow2 = max(0.0, W_c0(i,t,l) - Wmax)
    
    W_c0(i,t,l) = W_c0(i,t,l) - Q_Oflow2
    
    QR(i,t,l) = min(trans, W_c0(i,t,l))
    
    W_c0(i,t,l) = W_c0(i,t,l) - QR(i,t,l)
    
    Lay_upt = Lay_upt + QR(i,t,l)
    
    ! FIX: Add safety check
    if (Wmax .gt. p_critD) then
      W_con_ll = W_con_ll + W_c0(i,t,l) / Wmax
    endif
    
    if (trans .lt. 0.0) then
      trans = 0.0
    else
      trans = trans - QR(i,t,l)
    endif
    
  enddo
  
  layer_con0 = max(0.0, min(1.0, W_con_ll / 5))
  
  W_con_ll = 0.0
  
  ! bulk soil below
  Wx(i,t) = Wx(i,t) + Q_Oflow + Q_per + Q_Oflow2
  Q_Oflow_b = max(0.0, Wx(i,t) - Wxmax)
  Wx(i,t) = Wx(i,t) - Q_Oflow_b
  
  Evap = min(Wx(i,t), max(0.0, trans))
  
  Wx(i,t) = Wx(i,t) - Evap
  
  S_2 = max(0.0, min(1.0, Wx(i,t) / Wxmax))
  
  Q_Oflow_b2 = max(0.0, Wx(i,t) - Wxmax)
  
  Wx(i,t) = Wx(i,t) - Q_Oflow_b2
  
  Q_base(i,t) = min(Qb0 * S_2 * p_dt, Wx(i,t))
  
  Wx(i,t) = Wx(i,t) - Q_base(i,t)
  
  fH2Olg_ga1 = fH2Olg_ga_1 + Lay_upt + Evap
  
  Lay_upt = 0.0
  
  Runoff(i,t) = Q_base(i,t) + Q_Oflow_b + Q_Oflow_b2
  
  fH2Ol_xd_land = Runoff(i,t)
  
  Runoff_land = fH2Ol_xd_land
  
  ! FIX: Add safety check
  if (Wxmax .gt. p_critD) then
    bucket_con0 = Wx(i,t) / Wxmax
  else
    bucket_con0 = 0.0
  endif
  
endif ! surface not frozen

! Ground heat
fQ_tg0 = kSOIL(i) * (xT_s0 - xT_g0(i,t)) / p_dz_SOIL * lground

! Heat balance
! FIX: Add safety check
if (CSOIL(i) .gt. p_critD .and. p_dz_SOIL .gt. p_critD) then
  xT_g0(i,t) = xT_g0(i,t) + fQ_tg0 / CSOIL(i) / p_dz_SOIL * p_dt
endif

a0_Ts = xT_s0

return
end subroutine land_stepvTrans
!**************************************************************************
subroutine check_initialization(nCPts3)
use lycom_par



implicit none

integer :: nCPts3, i, t, l

write(*,*) "=== INITIALIZATION CHECK ==="
write(*,*) "por =", por
write(*,*) "Wmax =", Wmax
write(*,*) "Wxmax =", Wxmax
write(*,*) "nsoil =", nsoil

do i = 1, min(3, nCPts3)
  do t = 1, p_ntiles
    if (Wx(i,t) .le. 0.0) then
      write(*,*) "WARNING: Wx not initialized at i=", i, " t=", t
    endif
    if (xT_g0(i,t) .le. 0.0) then
      write(*,*) "WARNING: xT_g0 not initialized at i=", i, " t=", t
    endif
    do l = 1, nsoil
      if (W_c0(i,t,l) .le. 0.0) then
        write(*,*) "WARNING: W_c0 not initialized at i=", i, " t=", t, " l=", l
      endif
    enddo
  enddo
enddo
write(*,*) "=== END CHECK ==="

end subroutine check_initialization

end module lycom_land






