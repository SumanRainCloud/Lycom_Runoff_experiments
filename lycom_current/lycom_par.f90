
module lycom_par

!***********************************************************************
! run control
!***********************************************************************

integer                 :: year                                         ! year in simulation
integer                 :: month                                        ! month in year
integer                 :: dpm                                          ! days per month
integer                 :: day                                          ! day in simulation
integer                 :: ts                                           ! time step in day
integer                 :: tspd                                         ! time steps per day
integer                 :: accts                                        ! time step accumulated over simulation

integer, dimension(12)   :: months30                                     ! months with 30 days

integer, allocatable    :: naccu(:)                                     ! accumulation counter
integer, allocatable    :: naccu_d(:)                                   ! accumulation counter day 
integer, allocatable    :: naccu_n(:)                                   ! accumulation counter night
integer, allocatable    :: naccu_m(:)                                   ! accumulation counter months
!integer                 :: p_nspec =100
logical                 :: leapyear                                     ! switch for leap years
logical                 :: writeout                                     ! switch for writing output
logical                 :: out1                                         ! switch for opening output files
!integer                 :: pnspec = 100
real                    :: p_dt = 3600                                        ! time step length [s]

integer                 :: nCPts                                        ! number of (potential grid points / sites)
!integer                 :: tstep_total=(cyear - cyear0) *30.5 * 24 * 12 !263520 !30.5(days) * 24(hour)* 12(month)*30(year)
integer                 :: Total_timstep = 86400
! namelist parameters - with default values

integer                 :: cyear0               = 2010                  ! first calendar year
integer                 :: cyear                = 2020                  ! calendar year
integer                 :: lastyear             = 10                   ! last year of simulation
integer                 :: tsl                  = 3600                  ! time step length [s]
integer                 :: yearout1             = 1                   ! first year for writing output
integer                 :: yearoutX             = 2                   ! last year for writing output
integer                 :: outint               = 1                     ! output interval [time steps]
integer                 :: nSites               = 1                     ! number of local sites
integer                 :: p_nspec              = 100                             ! initial number of species
real                    :: frac_s_init          = 0.5                   ! initial surface coverage
logical                 :: specout              = .true.                                                                        ! DEACTIVATED - remove soon!             ! write fluxes for each species
logical                 :: interCan             = .false.               ! canopy species affect water at soil surface
logical                 :: NOHONO               = .false.               ! calculate NO/HONO emissions
logical                 :: BSCtypes             = .false.               ! calculate cover of BSC types
logical                 :: noVeg                = .false.               ! run without vegetation
integer                 :: inDirect             = 0                     ! direct climate data input
integer                 :: tstep_total=86400!263520!(cyear - cyear0) *30.5 * 24 * 12 !263520 !30.5(days) * 24(hour)* 12(month)*30(year)
integer                 :: tim  = 0
logical                 :: lrestart = .false.
!integer                 :: simhour = lastyear *30.5*24 *12
!***********************************************************************
! MPI
!***********************************************************************

logical                 :: para                                         ! switch for parallel simulation

integer                 :: rank                                         ! rank of processor
integer                 :: runperiod                                         ! rank of processor
integer                 :: tpos0                                         ! rank of processor
integer                 :: tsindata0                                         ! rank of processor
integer                 :: accts0                                         ! rank of processor
integer                 :: year0                                         ! rank of processor
integer                 :: numproc                                      ! number of processors
integer                 :: nx, ny                                       ! lon, lat
integer                 :: nt                                           ! time
integer                 :: nland                                        ! number of land points
integer                 :: pp                                           ! land points per processor
integer                 :: ppnp                                         ! total land points + extra points
integer                 :: tsindata                                     ! time step for global input data
integer                 :: tpos                                         ! time posiition for output file
integer                 :: t_varID                                      ! time variable for output file
integer                 :: outID                                        ! ID for output file
integer                 :: outIDS                                       ! ID for specific output file
integer, dimension(3)   :: dimIDs3                                      ! dimension IDs for output file

integer, allocatable    :: ppvec(:)                                     ! list of points per processor**
integer, allocatable    :: indvec(:,:)                                  ! list of land point indices**
integer, allocatable    :: indvec1(:,:)                                  ! list of land point indices**
!integer, dimension(86400,144,143) :: ncIIDD
integer, dimension(7)   :: ncID0                                        ! file index
integer, dimension(200) :: outvarID                                     ! output variable index
integer, dimension(200) :: outvarIDS                                    ! species output variable index
integer, dimension(100) :: outvarID2                                    ! output variable index (optional)

real, allocatable       :: lsdata(:,:)                                  ! land mask**
real, allocatable       :: xpos(:)                                      ! lon coordinates**
real, allocatable       :: ypos(:)                                      ! lat coordinates**

real, allocatable       :: indata(:,:)                                  ! global input data**

real, allocatable       :: R_net(:,:) 
real, allocatable       :: ta_C(:,:)
real, allocatable       :: desat_dT3(:,:)
real, allocatable       :: ET_pot0(:,:)
real, allocatable       :: E_Tdavg(:)
real, allocatable       :: E_Trmax(:)
real, allocatable       :: E_Tratio(:,:)
real, allocatable       :: indataL(:)                                   ! global input data (list)**
real, allocatable       :: tair_for_avg(:,:)
real, allocatable       :: srad_for_avg(:,:)
real, allocatable       :: lrad_for_avg(:,:)
real                    :: ETrmax2 = 0.0
real                    :: ETdavg2 = 0.0
!real                    :: ETdavg = 0.0
real, allocatable       :: indata7(:,:)                                 ! global input data (list,7 files)**
real, allocatable       :: ETrmaxdata(:)                          ! global input data (list)
real, allocatable       :: ETdavgdata(:)                             ! global input data (list)
real, allocatable       :: indata_etmapdata(:)                          ! global input data (list)
real, allocatable       :: indata_etravg(:)                             ! global input data (list)
real, allocatable       :: bcdata(:,:,:)                                ! boundary condition data**
real, allocatable       :: bcdataL(:,:)                                 ! boundary condition data (list)**
real, allocatable       :: bcdata2(:,:)                                 ! boundary condition data**
real, allocatable       :: bcdata2L(:)                                  ! boundary condition data (list)**
real, allocatable       :: laidata(:,:)                                 ! LAI
real, allocatable       :: saidata(:,:)                                 ! SAI
real, allocatable       :: ssadata(:)                                   ! Soil Surface Area
real, allocatable       :: biomedata(:)                                 ! biomes

character (len=*), parameter :: landmask        = "landsea.nc"          ! land mask file

character (len=*), parameter :: varName         = "var1"                ! name of the input variables
character (len=*), parameter :: tairfileG       = "tair.nc"             ! air temperature [K]
character (len=*), parameter :: rhumfileG       = "rhum.nc"             ! relative humidity []
character (len=*), parameter :: windfileG       = "wind.nc"             ! wind speed [m/s]
character (len=*), parameter :: rainfileG       = "rain.nc"             ! rainfall [m/s]
character (len=*), parameter :: snowfileG       = "snow.nc"             ! snowfall [m/s]
character (len=*), parameter :: sradfileG       = "srad.nc"             ! shortwave radiation [W/m2]
character (len=*), parameter :: lradfileG       = "lrad.nc"             ! longwave radiation [W/m2]

character (len=*), parameter :: laifileG        = 'LAI.nc'              ! LAI
character (len=*), parameter :: saifileG        = 'SAI.nc'              ! SAI
character (len=*), parameter :: ssafileG        = 'SSA.nc'              ! soil surface area
character (len=*), parameter :: biomefileG      = 'biome.nc'            ! Biomes
character (len=*), parameter :: ETrmaxfileG      = 'ETrmax.nc'            ! Evapotranspiration
character (len=*), parameter :: ETdavgfileG      = 'ETdavg.nc'            ! Evapotranspiration




!***********************************************************************
! Natural constants
!***********************************************************************

real, parameter         :: c_sigma              = 5.67E-8               ! Stefan-Boltzmann constant [W / (m2 * K4)]
real, parameter         :: c_HH2Olg             = 2.45E6                ! enthalpy of vaporisation [J / kg]
real, parameter         :: c_CAIR               = 1297.0                ! heat capacity of air [J / (m3 * K)]
real, parameter         :: c_rhoH2Ol            = 1000.0                ! density of liquid water [kg / m3]
real, parameter         :: c_rhoAir             = 1.3                   ! density of air [kg Air / m3]
real, parameter         :: c_rhoOrg             = 700.0 ! 1000.0                ! density of cell wall tissue [kg B / m3]
real, parameter         :: c_gamma              = 65.0                  ! psychrometric constant [kg * m / (s2 * m2 * K) = Pa / K]
real, parameter         :: c_MH2O               = 0.018                 ! molar weight of water [kg / mol]
real, parameter         :: c_MCo2               = 44
real, parameter         :: c_MC                 = 0.012                 ! molar weight of carbon [kg / mol]
real, parameter         :: c_Rgas               = 8.3145                ! universal gas constant  [J / (mol * K)]
real, parameter         :: c_TH2Osl             = 273.0                 ! melting temperature of water [K]
real, parameter         :: c_year               = 31536000.0            ! seconds per year [s]
real, parameter         :: c_month              = 2628000.0             ! seconds per month [s]
real, parameter         :: c_day                = 86400.0               ! seconds per day [s]
real, parameter         :: c_hour               = 3600.0                ! seconds per hour [s]

                                                                        
!***********************************************************************
! Environmental parameters
!***********************************************************************

! LAI and SAI files

integer, parameter      :: p_stepLAI            = 12                    ! temporal resolution of canopy properties
                                                                       
! Solubility

real, parameter         :: p_sCO2gl             = 0.0334                ! solubility of CO2 [mol/l] (= S [g/l] / M [g/mol])                             REF: vonCaemmerer,2000
real, parameter         :: p_sO2gl              = 0.00126               ! solubility of O2 [mol/l] (= S [g/l] / M [g/mol])                              REF: vonCaemmerer,2000

! Radiation

real                    :: p_eps                = 0.97                  ! average emissivity of land surface []                                         REF: Brutsaert,1982
real, parameter         :: p_beer_s             = 0.5                   ! extinction coefficient for short-wave radiation []                            REF: Bonan,2008
real, parameter         :: p_beer_l             = 0.95                  ! extinction coefficient for long-wave radiation []                             REF: Kustas,2000
real, parameter         :: p_albVEG             = 0.1                   ! albedo of bare canopy / gras []                                               REF: -
real, parameter         :: p_albDES             = 0.3                   ! albedo of desert soil []                                                      REF: -

! Vapour pressure

real, parameter         :: p_esatAIR1           = 17.27                 ! parameter for saturation vapour pressure []                                   REF: Allen,1998
real, parameter         :: p_esatAIR2           = 237.3                 ! parameter for saturation vapour pressure [deg C]                              REF: Allen,1998
real, parameter         :: p_esatAIR3           = 610.8                 ! parameter for saturation vapour pressure [Pa]                                 REF: Allen,1998

! Aerodynamic resistance to heat transfer

real, parameter         :: p_vonKarman          = 0.41                  ! von Karman constant [ ]                                                       REF: Stull,2003
real, parameter         :: p_mheight            = 10.0                  ! measurement height for wind speed [m]                                         REF: Stull,2003
real, parameter         :: p_roughlen_leaves    = 0.5                   ! roughness length of leaves [m]                                                REF: Stull,2003
real, parameter         :: p_roughlen_stems     = 0.5                   ! roughness length of stems [m]                                                 REF: Stull,2003
real, parameter         :: p_roughlen_grd       = 0.25 !ALMERIA   0.05  ! roughness length of soil surface [m]                                          REF: Stull,2003
real, parameter         :: p_roughlen_floor     = 0.05                  ! roughness length of forest floor [m]                                          REF: Stull,2003
real, parameter         :: p_ratio_mh           = 0.1                   ! ratio between roughness length of momentum and humidity [ ]                   REF: Allen,1998
real, parameter         :: p_ratio_dr           = 2.0/3.0/0.123         ! ratio between displacement height and roughness length [ ]                    REF: Allen,1998

! Ground heat flux

real, parameter         :: p_dz_SOIL            = 0.15                  ! damping depth of the soil for a diurnal cycle [m]                             REF: Bonan,2008 (page 134)
!real, parameter         :: p_dz_AIR             = 0.72                  ! damping depth of air for a diurnal cycle [m]                                  REF: sqrt(86400 * D / Pi); Therm.Diff.(D) = 1.9e-5 m2/s (Wikipedia)
!real, parameter         :: p_dz_WOOD            = 0.05                  ! damping depth of the soil for a diurnal cycle [m]                             REF: sqrt(86400 * D / Pi); Therm.Diff.(D) = 8.2e-8 m2/s (Wikipedia)
real, parameter         :: p_CSOIL_des          = 1.5E6                 ! heat capacity of desert soil [J / m3 / K]                                     REF: Ekici,2014
real, parameter         :: p_CSOIL_for          = 2.5E6                 ! heat capacity of forest soil [J / m3 / K]                                     REF: Ekici,2014
real, parameter         :: p_kSOIL_des          = 0.9 !ALMERIA   0.3                   ! thermal conductivity of desert soil [W / m / K]                               REF: Ekici,2014
real, parameter         :: p_kSOIL_for          = 0.5                   ! thermal conductivity of forest soil [W / m / K]                               REF: Ekici,2014
!real, parameter         :: p_kAIR               = 0.026                 ! thermal conductivity of air [W / m / K]                                       REF: Wikipedia
!real, parameter         :: p_kWOOD              = 0.2                   ! thermal conductivity of wood [W / m / K]                                      REF: Wikipedia

! Dew

real, parameter         :: p_dewmax_des         = 20.0 !40.0                  ! maximum dew in dry regions [mm / yr]                                          REF: Vuollekoski,2015
real, parameter         :: p_dewmax_for         = 30.0 !60.0                  ! maximum dew in moist regions [mm / yr]                                        REF: Vuollekoski,2015 

! Snow layer

real, parameter         :: p_kH2Os              = 0.15                  ! thermal conductivity of average snow [W / m / K]                              REF: Ekici,2014
real, parameter         :: p_rhoH2Os            = 250.0                 ! density of average snow [kg / m3]                                             REF: Ekici,2014
real, parameter         :: p_H2Os_crit          = 0.1                   ! critical snow depth [m]                                                       REF: Kappen,1995
real, parameter         :: p_H2Os_loss          = 3.2E-10               ! glacier loss [1 / s]                                                          REF: F.Gans

! Interception

real, parameter         :: p_rmaxH2Ol_cl0       = 0.0002                ! water storage on leaves [m]                                             REF: Miralles,2010
real, parameter         :: p_rmaxH2Ol_cs0       = 0.00005               ! water storage on stems [m]                                             REF: Miralles,2010
real, parameter         :: p_rmaxH2Ol_s0        = 0.0002                ! water storage at soil surface [m]                                             REF: Miralles,2010
real, parameter         :: p_fracIcpt           = 0.5 !0.7                  ! shape of interception efficiency curve [ ]

! Disturbance

real, parameter         :: p_tauLEAF_trop       = 17.0                  ! leaf lifespan of broadleaf evergreen trees [months]                           REF: Porada,2013
real, parameter         :: p_tauLEAF_bor        = 72.0                  ! leaf lifespan of needleleaf evergreen trees [months]                          REF: Porada,2013
real, parameter         :: p_tauLEAF_med        = 27.0                  ! leaf lifespan of mediterranean trees [months]                                 REF: Porada,2013
!real, parameter         :: p_tauLEAF_temp       = 6.0                   ! leaf lifespan of temperate trees [months]                                     REF: Porada,2013
real, parameter         :: p_tauG_tropwet     = 100.0                   ! turnover time of tropical rainforest [years]                                  REF: Porada,2013
real, parameter         :: p_tauG_tropdry     = 32.0                    ! turnover time of tropical dry forest [years]                                  REF: Porada,2013
real, parameter         :: p_tauG_tropground  = 100.0                   ! turnover time of tropical rainforest [years]                                  REF: Porada,2013
real, parameter         :: p_tauG_temp        = 100.0                   ! turnover time of temperate forest [years]                                     REF: Porada,2013
real, parameter         :: p_tauG_bor         = 100.0                   ! turnover time of boreal forest [years]                                        REF: Porada,2013
real, parameter         :: p_tauG_med         = 30.0 !50.0                    ! turnover time of mediterranean vegetation [years]                             REF: Porada,2013
real, parameter         :: p_tauG_sav         = 30.0                    ! turnover time of savanna [years]                                              REF: Porada,2013
real, parameter         :: p_tauG_gra         = 30.0                    ! turnover time of prairie & montane grassland [years]                          REF: Porada,2013
real, parameter         :: p_tauG_tun         = 50.0                    ! turnover time of tundra [years]                                               REF: Porada,2013
real, parameter         :: p_tauG_des         = 100.0                   ! turnover time of desert [years]                                               REF: Porada,2013
real, parameter         :: p_tauG_pal         = 100.0                   ! turnover time of Ordovician landscape [years]                                 REF: -

! Area

real, parameter         :: p_LAImax             = 5.7                   ! Maximum LAI in data set [ ]                                                   REF: Bonan,2002
real, parameter         :: p_LAIcorr            = 4.0                   ! Correction factor for forest cover [ ]
real, parameter         :: p_xcano              = 0.7                   ! weight for average posiition in canopy [ ]                                    REF: -

! Soil properties
real, parameter         :: p_rmaxH2Ol_g1        = 0.0005*0.35            ! topsoil water storage capacity [ m3 / m2 G ] -> 1.5cm depth, 35% porosity
real, parameter         :: p_rmaxH2Ol_g2        = 1.5*0.35              ! bulk soil water storage capacity [ m3 / m2 G ] -> 1.5m (rooting) depth, 35% porosity
real, parameter         :: p_kH2Ol_sv           = 5.0e-8                ! soil-plant-atmosphere conductivity [ m3 / (m2 G * s) ]
real, parameter         :: p_kH2Ol_sb           = 3.0e-8                ! soil-groundwater conductivity [ m3 / (m2 G * s) ]
real, parameter         :: por                  = 0.45
real, parameter         :: Wmax                 = 0.45*0.03
real, parameter         :: Wxmax                = 0.45*0.65
real, parameter         :: Qp0                  = 0.5e-6!7
real, parameter         :: Qb0                  = 1.0e-7!8

! Surface properties

!real, parameter         :: p_rmaxH2Ol_bg        = 0.005                 ! soil water accessible to NVV [ m ]
real, parameter         :: p_rmaxH2Ol_bs        = 0.0014                ! bark water accessible to NVV [ m ]

real, parameter         :: p_ksatBark           = 1.5e-6                ! saturated hydraulic conductivity [ m/s ]
real, parameter         :: p_fracBark           = 0.5                   ! fraction stemflow / layer overflow [ ]

real, parameter         :: p_fracWup0L          = 0.4                   ! water uptake efficiency bare leaves [ ]
real, parameter         :: p_fracWup0S          = 0.1                   ! water uptake efficiency bare stems [ ]
real, parameter         :: p_fracWup            = 0.6                   ! water uptake efficiency NVV [ ]


!Lycophyte properties
real, parameter         :: rstem                 =0.02
real, parameter         :: pi                    =2.314
!real, parameter         :: max_area_grid         =

!***********************************************************************
! Lichen and bryophyte parameters
!***********************************************************************

! Light absorption:

real                    :: p_alb_h              = 0.4                   ! maximum lichen albedo []                                                      REF: Porada,2013(Thesis)
real                    :: p_alb_l              = 0.05                  ! minimum lichen albedo []                                                      REF: Porada,2013(Thesis)
real                    :: p_LAInvv_h           = 20.0                  ! maximum LAI of NVV []                                                         REF: -
real                    :: p_LAInvv_l           = 1.0                   ! minimum LAI of NVV []                                                         REF: -
real                    :: p_extL               = 0.5                   ! parameter for light extinction []                                             REF: -

! Lichen height:

real                    :: p_zt_l               = 0.001                 ! minimum thallus height (photosynthesising) [m]                                REF: Porada,2013(Thesis)
real                    :: p_zt_h               = 0.15                  ! maximum thallus height (photosynthesising) [m]                                REF: Porada,2013(Thesis)

! Lichen porosity:

real                    :: p_prs_l              = 0.1                   ! minimum porosity []                                                           REF: Porada,2015(JSBACH)
real                    :: p_prs_h              = 0.99                  ! maximum porosity []                                                           REF: Porada,2015(JSBACH)
real                    :: p_fracA_l            = 0.05                  ! minimum air fraction at saturation []                                                           REF: Porada,2015(JSBACH)
real                    :: p_fracA_h            = 0.8                   ! maximum air fraction at saturation []                                                           REF: Porada,2015(JSBACH)

! Lichen density:

!real                    :: p_rhoCb              = 120.0                 ! density of biomass [kg B / m3]                                                REF: Porada,2015(JSBACH)

! Hydrophobicity

real                    :: p_redHph             = 0.05                  ! fraction of thallus reservoir which can be filled per hour
real                    :: p_satHph_l           = 0.0                   ! minimum saturation below which hydrophobicity occurs []
real                    :: p_satHph_h           = 0.3                   ! maximum saturation below which hydrophobicity occurs []

! Lichen activity:

real                    :: p_sat_act0           = 0.03 !0.08            ! saturation at minimum water content for activation []                         REF: NashIII,1996
real                    :: p_sat_actF_l         = 0.1                   ! minimum saturation for full activity (1.75-3.75) []                  REF: NashIII,1996
real                    :: p_sat_actF_h         = 0.5                   ! maximum saturation for full activity (1.75-3.75) []                  REF: NashIII,1996
real                    :: p_sat_act0CCM         = 0.15                 ! saturation at minimum CCM water content for activation (2.1-2.5) []           REF: NashIII,1996
!real                    :: p_sat_act2CCM_l      = 0.5                   ! minimum saturation at full CCM cell water content (7.5-15) []                 REF: NashIII,1996
!real                    :: p_sat_act2CCM_h      = 1.0                   ! maximum saturation at full CCM cell water content (7.5-15) []                 REF: NashIII,1996

!real                    :: p_act_time_h         = 48.0*3600.0           ! maximum activation time [s]                                                   REF: -
!real                    :: p_act_time_l         = 0.5*3600.0            ! minimum activation time [s]                                                   REF: -
real                    :: p_satX_l             = 0.01                  ! minimum saturation for full activity []                                       REF: NashIII,1996
real                    :: p_satX_h             = 0.95                  ! minimum saturation for full activity []                                       REF: NashIII,1996

! Lichen CO2 diffusion:

real                    :: p_kCO2g_max          = 0.15                  ! maximum thallus conductivity for CO2 [mol / (m2 T * s)]                       REF: Cowan,1992; Williams,1998
real                    :: p_kCO2g_satl         = 2.0e-3 ! 5.7e-4                ! minimum DCO2 at saturation [mol / (m2 T * s)]                                 REF: Cowan,1992; Williams,1998
real                    :: p_kCO2g_sath         = 0.01                  ! maximum DCO2 at saturation [mol / (m2 T * s)]                                 REF: Cowan,1992; Williams,1998
real                    :: p_kCO2gB_h           = 12.0                  ! shape parameter for curve                                                     REF: Cowan,1992; Williams,1998
real                    :: p_kCO2gB_l           = 2.0                   ! shape parameter for curve                                                     REF: Cowan,1992; Williams,1998
!real                    :: p_kCO2g_1            = 1.0                   ! shape parameter for curve                                                     REF: Cowan,1992; Williams,1998
!real                    :: p_kCO2g_2            = 2.0 !10.0                  ! shape parameter for curve                                                     REF: Cowan,1992; Williams,1998

! Lichen H2O diffusion:

real                    :: p_rs_h               = 50.0                  ! maximum thallus resistance to H2O [s / m]                                     REF: Monteith,1981
real                    :: p_rs_l               = 0.0                   ! minimum thallus resistance to H2O [s / m]                                     REF: Monteith,1981

! Photosynthesis model:

real                    :: p_conv_PAR           = 2.0699E-6             ! conversion of W/(m2) to [mol / (m2 * s)]                                      REF: Ting,1987
real                    :: p_conv_alpha         = 0.5*0.83*0.84         ! conversion of 2 quanta light to 1 electron (* quantum yield & absorption)     REF: Skillman,2008; Schulze,1995(Book)
real                    :: p_vcmaxM_h           = 26.8                  ! maximum value of molar vcmax [1 / s], min: 1.2(5.8) max: 13.4(13.4)           REF: Porada,2013
real                    :: p_vcmaxM_l           = 0.0139                   ! minimum value of molar vcmax [1 / s]                                          REF: Porada,2013
real                    :: p_vomaxM_h           = 2.5                   ! maximum value of molar vomax [1 / s], min: 0.5(0.8) max: 2.1(1.6)             REF: Porada,2013
real                    :: p_vomaxM_l           = 0.391                   ! minimum value of molar vomax [1 / s]                                          REF: Porada,2013
real                    :: Trans_h           = 0.75                   ! Light transmission                                          REF: Porada,2013
real                    :: Trans_l           = 0.60                   ! Transmission                                          REF: Porada,2013

real                    :: Rub_h              = 10.0E-6               ! maximum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013
real                    :: Rub_l              = 8.5E-6                ! minimum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013
real                    :: resp_rub_h         = 0.1               ! maximum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013
real                    :: resp_rub_l         = 0.01                ! minimum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013
real                    :: gS0_h              = 0.35               ! maximum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013
real                    :: gS0_l              = 0.20                ! minimum value of specific Rubisco content [mol / m2 T]                        REF: Porada,2013



real                    :: p_RR                 = 5.0 !3.3                        ! factor to convert respiration to Rubisco [s]                                  REF: Porada,2013
real                    :: p_KcM1               = 1.32                  ! parameter for KcM as a function of vcmax []                                   REF: Savir,2009
real                    :: p_KcM2               = 2.03                  ! parameter for KcM as a function of vcmax []                                   REF: Savir,2009
real                    :: p_KoM1               = 0.0057                ! parameter for KoM as a function of vcmax and vomax []                         REF: Savir,2009
real                    :: p_KoM2               = 0.51                  ! parameter for KoM as a function of vcmax and vomax []                         REF: Savir,2009

! Kattge 2007: Ea(Vcmax): 40-110, Ea(Jmax): 30-80
! Medlyn 2002: Ea(Kc): 60-120, Ea(Ko): 10-50
real                    :: p_EaKc_h             = 120000.0 !130000.0              ! maximum enzyme activation energy of Kc [J / mol]                              REF: Porada,2013
real                    :: p_EaKc_l             = 50000.0 !30000.0               ! minimum enzyme activation energy of Kc [J / mol]                              REF: Porada,2013
real                    :: p_EaKo_h             = 50000.0 !55000.0               ! maximum enzyme activation energy of Ko [J / mol]                              REF: Porada,2013
real                    :: p_EaKo_l             = 10000.0 !5000.0                ! minimum enzyme activation energy of Ko [J / mol]                              REF: Porada,2013

real                    :: p_EaVm_h             = 110000.0              ! maximum enzyme activation energy of Vcmax [J / mol]
real                    :: p_EaVm_l             = 40000.0               ! minimum enzyme activation energy of Vcmax [J / mol]
real                    :: p_EaJm_h             = 80000.0               ! maximum enzyme activation energy of Jmax [J / mol]
real                    :: p_EaJm_l             = 30000.0               ! minimum enzyme activation energy of Jmax [J / mol]

real                    :: p_ToptP_max          = 323.0                 ! maximum optimum temperature photosynthesis [K]                                REF: NashIII,1996 (page 206)
real                    :: p_ToptP_min          = 273.0                 ! minimum optimum temperature photosynthesis [K]                                REF: NashIII,1996 (page 206)
!real                    :: p_Tref_PS            = 298.0                 ! reference temperature of photosynthesis [K]                                   REF: June,2004
!real                    :: p_omega              = 18.0                  ! shape parameter for temperature response [K]                                  REF: June,2004
                                                                                                   
! Carbon concentration mechanism                                        
                                                                               
real                    :: p_conv_CCM           = 2.0 / 3.0             ! instead of 2 quanta light, 3 are required for fixation of 1/4 CO2             REF: Fridlyand,1996
real                    :: p_NCCM               = 40.0                  ! efficiency of CCM                                                             REF: Porada,2013
real                    :: p_FCCM               = 0.5                   ! default fraction of species with CCM []
                                                                               
! Respiration:                                                           
                                                                               
real                    :: p_Q10R_h             = 2.3                   ! maximum Q10 value of respiration []                                           REF: Porada,2013
real                    :: p_Q10R_l             = 1.5                   ! minimum Q10 value of respiration []                                           REF: Porada,2013
!real                    :: p_Tref_R             = 283.0                 ! reference temperature of respiration [K]
real                    :: p_Rref_h             = 25.0E-6               ! maximum reference Respiration [mol / (kg C * s)]                              REF: Porada,2013
real                    :: p_Rref_l             = 1.0E-6                ! minimum reference Respiration [mol / (kg C * s)]                              REF: Porada,2013
real                    :: p_NCcb               = 1.0 !!!onlyNCB 0.75                  ! efficiency of sugar to biomass conversion []                                  REF: Cannell,2000
real                    :: cpar                 =2.07E-6
real                    :: p_gS1                =1000
real                    :: fracR0=0.25
real                    :: Ngrow =0.75
! Growth:                                                           
                                                                               
!real                    :: p_Q10G_h             = 4.5                   ! maximum Q10 value of growth []                                                REF: Porada,2016
!real                    :: p_Q10G_l             = 1.5                   ! minimum Q10 value of growth []                                                REF: Porada,2016
!real                    :: p_ToptG_max          = 303.0                 ! maximum optimum temperature growth [K]                                        REF: Porada,2016
!real                    :: p_ToptG_min          = 278.0                 ! minimum optimum temperature growth [K]                                        REF: Porada,2016
!real                    :: p_Tref_G             = 263.0                 ! reference temperature of initial growth [K]                                   REF: Porada,2016
!real                    :: p_Tref_Gx            = 313.0                 ! reference temperature of growth limit [K]                                     REF: Porada,2016
!real                    :: p_Tref_Gz            = 333.0                 ! reference temperature of cell death [K]                                       REF: Porada,2016

! Litterfall                                                             
                                                                               
real                    :: p_turnover_h         = 1.0                   ! maximum turnover rate [ 1 / yr ]                                              REF: Porada,2013
real                    :: p_turnover_l         = 0.01                  ! minimum turnover rate [ 1 / yr ]                                              REF: Porada,2013
                                                                               
! Competition                                                        

real                    :: p_NCbt               = 0.01 !0.85                  ! efficiency of growth to new cover conversion due to dispersal []              REF: -
!real                    :: p_areaTH_s0          = 0.01                  ! efficiency of growth to new cover conversion due to dispersal []              REF: -

! Lichen critical conditions
real, parameter         ::fH2Ol_tbf            =0.5 * 0.0005*0.35
real, parameter         ::fH2Ol_tbf1           =0.5 * 0.0005*0.35
real,allocatable        :: as_fH2Ol_td(:)
real,allocatable        :: as_gpp(:)
real,allocatable        :: as_npp(:)
real, parameter         :: B_area              = 0.100 !(55 g/m2) literature leaf density per sq. meter
real, parameter         :: p_rhoBl             = 350.0  !!!!!NEEDS PROPER NUMBER
real, parameter         :: Brmax0              = 200.0
real, parameter         :: Brmax               = 200.0*0.03
real, parameter         :: kmrt0               = 1.0E-7
real, parameter         :: fracratiocrit       = 0.0001                   ! critical fraction of fracTH_s_init 0.1
real                    :: monty               =0.0
real,allocatable        :: nfd(:)
real,allocatable        :: MD0(:) 
real,allocatable        :: M0(:,:)
real,allocatable        :: Mrt(:)
real,allocatable        :: Total_mortality_root(:,:)
real,allocatable        :: RtoM(:,:)
real,allocatable        :: Rtb(:,:)
real,allocatable        :: Run_tot(:,:,:)
real,allocatable        :: Rtres_l(:)
real                    :: W_con_l=0.0
real                    :: W_con_ll=0.0
real                    :: val_soilco2max= 5000000
!***********************************************************************
! Species parameters
!***********************************************************************

integer, parameter      :: p_nspecpar           = 16                    ! number of species parameters

real, allocatable       :: vec_o(:,:)                                   ! vector of species parameters

real, allocatable       :: o_albedo2(:)                                 ! pure NVV albedo []
real, allocatable       :: o_theta_max(:)                               ! water storage capacity [kg H2O / kg C]
real, allocatable       :: o_spec_area(:)                               ! specific area [m2 T / kg C]
real, allocatable       :: o_sat_actF(:)                                ! water saturation necessary for full activity []
real, allocatable       :: o_sat_X(:)                                   ! water saturation when potential becomes negative []
real, allocatable       :: o_gS0(:)                                   ! water saturation when potential becomes negative []
real, allocatable       :: o_resp_main(:)                                   ! Respiration
real, allocatable       :: o_X(:)                                   
real, allocatable       :: o_G_area(:)                                   ! groundarea occupied by the plant []
real, allocatable       :: o_w_leaf(:)                                   ! dry weight of leaf []
real, allocatable       :: o_A_leaf(:)                                   ! Area of Leaf  []
real, allocatable       :: o_fracTransm(:)                                   ! water saturation when potential becomes negative []
real       :: ETdata=0.0
real       :: ETratio=0.0
real, allocatable       :: Lai_cum(:)
!real       :: ETdavg
real       :: gS
integer                 :: counts = 0
integer,allocatable     :: counttimer(:)
real, allocatable       :: o_vcmax_M(:)                                 ! carboxylation rate of Rubisco (molar vcmax) [1 / s]
real, allocatable       :: o_RubConc(:)                                 ! Rubisco content [mol / m2 T]
real, allocatable       :: o_ratio_Resp_Rub(:)                          ! 
real, allocatable       :: o_spec_Rubisco(:)                            ! specific Rubisco content [mol / m2 T]
real, allocatable       :: o_turnover(:)                                ! turnover [mol / (m2 T * s)]
real, allocatable       :: o_vomax_M(:)                                 ! oxygenation rate of Rubisco (molar vomax) [1 / s]
real, allocatable       :: o_ToptP(:)                                   ! optimum temperature of photosynthesis [K]
real, allocatable       :: o_Q10_resp(:)                                ! Q10 value of respiration []
real, allocatable       :: o_Eact_Kc(:)                                 ! enzyme activation energy of Kc [J / mol]
real, allocatable       :: o_Eact_Ko(:)                                 ! enzyme activation energy of Ko [J / mol]
real, allocatable       :: o_Eact_Vm(:)                                 ! enzyme activation energy of Vcmax [J / mol]
real, allocatable       :: o_Eact_Jm(:)                                 ! enzyme activation energy of Jmax [J / mol]
real, allocatable       :: o_CCM(:)                                     ! fraction of electrons allocated to CCM []
real, allocatable       :: o_DCO2(:)                                    ! shape parameter for CO2 diffusivity curve []
real, allocatable       :: o_DCO2B(:)                                   ! shape parameter for CO2 diffusivity curve []
real, allocatable       :: o_rs(:)                                      ! thallus resistance to H2O diffusion [s / m]
real, allocatable       :: o_zt(:)                                      ! thallus height (photosynthesising) [m]
real, allocatable       :: o_prs(:)                                     ! thallus porosity []
!real, allocatable       :: o_act_time(:)                                ! activation time [s]
real, allocatable       :: o_satHph(:)                                  ! hydrophobicity of thallus []
real, allocatable       :: o_LAI(:)                                     ! LAI of non-vascular vegetation []
real, allocatable       :: mon_cond(:,:)
real, allocatable       :: M_cond(:)
real, allocatable       :: res_por(:)

!***********************************************************************
! Miscellaneous parameters
!***********************************************************************

real, parameter         :: p_critD              = 1.0E-12               ! critical value for any denominator

integer, parameter      :: p_ntiles             = 1!4                     ! number of tiles (1=forest, 2=bare/grass, 3=wetland, 4=rocks)
integer, parameter      :: p_nvertM             = 1                     ! max number of levels (1=ground, 2=canopy)
integer, parameter      :: p_nhabM              = 1                     ! max number of habitats (1=stems, 2=leaves)
integer, parameter      :: p_nhabA              = 1!4                     ! actual number of habitats (sum of stems, leaves, forest floor, bare/grass)
integer, parameter      :: nsoil=5
integer, dimension(p_ntiles) :: p_nvert                                 ! number of vertical levels per tile(n=4)
integer, dimension(p_ntiles,p_nvertM) :: p_nhab                         ! number of habitats per level(nmax=2)

!***********************************************************************
! Boundary conditions
!***********************************************************************

! Climate input data

real                    :: rCO2g_a                                      ! atmospheric CO2 [ppm]
real                    :: rO2g_a                                       ! atmospheric O2 [ppm]
real, allocatable       :: xT_a(:)                                      ! air temperature [K]
real, allocatable       :: fRADs_ad(:)                                  ! shortwave radiation [W / m2]
real, allocatable       :: fRADl_ad(:)                                  ! longwave radiation [W / m2]
real, allocatable       :: fH2Ol_ad(:)                                  ! precipitation [m3 / (m2 * s)]
real, allocatable       :: fH2Os_ad(:)                                  ! snowfall [m3 / (m2 * s)]
real, allocatable       :: rH2Og_RH(:)                                  ! relative humidity []
real, allocatable       :: fAIR_s(:)                                    ! surface wind [m / s]

! Ecosystem properties

real, allocatable       :: areaLEAF(:,:)                                ! maximum thallus area per m2 leaf surface []
real, allocatable       :: areaSTEM(:,:)                                ! maximum thallus area per m2 stem surface []
real, allocatable       :: areaSOIL(:)                                  ! maximum thallus area per m2 soil surface []
real, allocatable       :: areaLEAF_month(:)                            ! monthly maximum thallus area per m2 leaf surface []
real, allocatable       :: areaSTEM_month(:)                            ! monthly maximum thallus area per m2 stem surface []
real, allocatable       :: biome(:)                                     ! biome number []
real, allocatable       :: N_leaf(:)
real, allocatable       :: htree(:)
real, allocatable       :: Lai_new(:,:)
real                    :: Lai_ini=2.0
real, allocatable       :: CO2_sink(:,:)
real, allocatable       :: CO2_pre(:,:)
real, allocatable       :: ETrmax(:)
real, allocatable       :: ETdavg(:)
real, allocatable       :: Etr_avg(:)
real, allocatable       :: Etr_max(:)
! Derived BC

real, allocatable       :: tauD(:,:)                                ! disturbance interval [years]
real, allocatable       :: dn_kH2Og(:,:)                            ! denominator for roughness length of X surface [m]
real, allocatable       :: rmaxH2Ol_b(:,:)                          ! denominator for roughness length of X surface [m]
real                    :: lground                                      ! ground or canopy ?


!..............................................

real, allocatable       :: fracSTEM(:)                                  ! ...
real, allocatable       :: fracLEAF(:)                                  ! ...

real, allocatable       :: frac_noland(:)                               ! fraction of non-vegetated surface

real, allocatable       :: frac_tile(:,:)                               ! fraction of non-vegetated surface
real, allocatable       :: albLsf(:,:)                                  ! land surface albedo


!***********************************************************************
! Initial conditions
!***********************************************************************

real                    :: fracTH_s_crit                                ! critical area fraction [ ]
real                    :: frac_s_crit                                ! critical area fraction [ ]

!***********************************************************************
! Environmental state variables
!***********************************************************************

real, allocatable       :: rH2Os_g(:)                                   ! snow reservoir [m3 water / m2 ground]
real, allocatable       :: rH2Ol_0(:,:)                             ! water on leaves [m3 H2O / m2 ground]

real, allocatable       :: rH2Ol_g1(:,:)                                ! water in top soil [m3 H2O / m2 ground]
real, allocatable       :: rH2Ol_g2(:,:)                                ! water in bulk soil [m3 H2O / m2 ground]
real, allocatable       :: xT_g0(:,:)                                 ! bare soil temperature [K]

!***********************************************************************
! Environmental variables
!***********************************************************************

integer                 :: kday                                         ! day or night ?

! Available area

real                    :: Aslai                                        ! SLAI []
real, allocatable       :: Acano(:,:)                                   ! vegetation area in canopy [m2]

! Radiation
real                    :: fracRADs                                     ! fraction of absorbed short-wave radiation (average) []
real                    :: fracRADl                                     ! fraction of absorbed long-wave radiation []

! Vapour pressure

real                    :: Ta3                                          ! air temperature power3 [K]
real                    :: Ta4                                          ! air temperature power4 [K]
real                    :: zT_a                                         ! air temperature in [deg C]
real                    :: esatAIR                                      ! saturation vapour pressure of air [Pa]
real                    :: esatAIR2                                     ! corrected saturation vapour pressure of air [Pa]
real                    :: desatdT                                      ! slope of saturation vapour pressure vs temperature relationship []
real                    :: desatdT2                                     ! corrected slope of saturation vapour pressure vs temperature relationship []

! Aerodynamic resistance to heat transfer

!real, allocatable       :: kH2Og(:,:)                                   ! boundary layer conductance [m / s]
real                    :: roughlen                                     ! roughness length for momentum [m]
real                    :: roughlen_h                                   ! roughness length for humidity [m]
real                    :: dheight                                      ! displacement height of surface [m]
! Surface temperature

real                    :: xT_s0

! Ground heat flux

real, allocatable       :: CSOIL(:)                                     ! heat capacity of soil [J / m3 / K]
real, allocatable       :: kSOIL(:)                                     ! thermal conductivity of soil [W / m / K]

! Dew

real, allocatable       :: dewmax(:)                                    ! maximum rate of dew [m / s]

! Snow layer

real                    :: dsnow                                        ! snow depth  [m3 snow / m2 ground]
real                    :: fH2Osl_g                                     ! snowmelt [m3 / (m2 * s)]
real                    :: fH2Ol_tb_f

! Water fluxes


real,  allocatable      :: fracrain1(:,:)                                     ! fraction of rain absorbed in habitat []
real,  allocatable      :: fH2Ol_ci1(:,:)                                     ! fraction of rain absorbed in habitat []
real,  allocatable      :: fH2Ol_ts1(:,:)                                     ! fraction of rain absorbed in habitat []
real,  allocatable      :: fH2Ol_ux1(:,:)                                     ! water flux from upper layer to current layer [m3 H2O / (m2 ? * s)]
real,  allocatable      :: fracrain(:,:)  
real,  allocatable      :: fH2Ol_ci(:,:)                                     ! fraction of rain absorbed in habitat []
real,  allocatable      :: fH2Ol_ts(:,:)                                     ! fraction of rain absorbed in habitat []
real,  allocatable      :: fH2Ol_ux(:,:)                                     ! water flux from upper layer to current layer [m3 H2O / (m2 ? * s)]
real,  allocatable      :: fH2Ol_ug                                     ! water flux from ground layer into soil [m3 H2O / (m2 ? * s)]
real, allocatable       :: fH2Ol_ug2(:,:)                               ! throughfall (without canopy overflow) [m3 H2O / (m2 ? * s)]
real      :: fH2Ol_xd0=0.0                                    ! downward water flux from vasc layer [m3 H2O / (m2 ? * s)]
real      :: fH2Ol_xd1=0.0
real     :: fH2Olg_xu0=0.0                                   ! upward water flux from vasc layer [m3 H2O / (m2 ? * s)]
real      :: fH2Olg_xu1=0.0
real      :: fH2Ol_ub=0.0                                     ! stemflow (input into bark reservoir) [m3 H2O / (m2 ? * s)]
real      :: fH2Ol_bd0=0.0                                    ! stemflow (output from bark reservoir) [m3 H2O / (m2 ? * s)]
real      :: fH2Ol_gwl=0.0
real      :: fH2Ol_gwl1=0.0
real,  allocatable      :: Layer_con(:)
! Soil hydrology
real,  allocatable      :: Qin(:,:,:)
real,  allocatable      :: Qin1(:,:,:,:)
real      :: Q_Per=0.0
real      :: Q_Per1=0.0

real      :: Q_Oflow=0.0
real      :: Q_Oflow1=0.0
real      :: Q_Oflow2=0.0
real      :: Q_Oflow_b2=0.0
real      :: Q_Oflow_b=0.0
!real      :: rain =0.0!-1.5E-08
real,  allocatable      :: Q_base1(:,:,:)
real,  allocatable      :: Q_base(:,:)
real,  allocatable      :: W_c1(:,:,:,:)                                !Layer soil
real,  allocatable      :: W_c0(:,:,:)
real,  allocatable      :: S_wc1(:,:)
real                    :: S_wc0
real,  allocatable      :: S2(:)
real                    :: S_2=0.0 
real,  allocatable      :: QR(:,:,:)
real,  allocatable      :: QR1(:,:,:,:)

real,  allocatable      :: Wx(:,:)                                   !Bulk soil water content lycomland
real,  allocatable      :: Wx1(:,:,:)                                !Bulk soil water content lycom
real,  allocatable      :: Runoff(:,:)                               !Runoff from 
real,  allocatable      :: Runoff1(:,:,:)
real                    :: Runoff_land=0.0

real                    :: fH2Ol_go=0.0                                     ! Surface runoff [ m3 H2O / (m2 G * s) ]
real                    :: fH2Ol_gb=0.0                                     ! Baseflow [ m3 H2O / (m2 G * s) ]
real,allocatable        :: fH2Olg_ga(:)                                    ! soil evaporation and transpiration [ m3 H2O / (m2 G * s) ]
real                    :: fH2Olg_ga1=0.0
real                    :: fH2Olg_ga_1=0.0
! Energy balance

real                    :: fRAD_H0                                      ! net radiation [W / m2]
real                    :: fQ_tg0                                       ! Ground heat flux [W / m2]
real                    :: fQ_ta_L0                                     ! Latent heat flux [W / m2]
real                    :: fQ_ta_S0                                     ! Sensible heat flux [W / m2]
!***********************************************************************
! Lichen variables
!***********************************************************************

! Lichen Water content:

real, allocatable       :: rmaxH2Ol_t(:)                                ! maximum thallus water content [m] 
real, allocatable       :: sat(:)                                       ! thallus water saturation []
real                    :: sat_act                                      ! saturation at activation []

real                    :: fracMP                                       ! fraction of micropores []
real                    :: fracAir                                      ! fraction of air at saturation []

! Lichen water potential:

real                    :: xH2Ol =0.0                                       ! lichen water potential [J / kg]
real                    :: xH2Ol_b =0.0                                     ! lichen water potential below [J / kg]
real                    :: layer_con0=0.0
real                    :: bucket_con0=0.0
! Lichen surface:

real, allocatable       :: xT_s(:)                                      ! surface temperature [K]
!real                    :: fracRADs_0                                   ! fraction of absorbed short-wave radiation (integral) []
!real                    :: fracRADl_0                                   ! fraction of absorbed long-wave radiation (integral) []
! Lichen CO2 diffusion

real, allocatable       :: CO2_p(:)                                     ! CO2 concentration in pore space [ppm]
real, allocatable       :: Rtres(:,:)
real                    :: kCO2g=0.0                                        ! thallus conductivity for CO2 [1 / s]

! Photosynthesis model:
real                    :: KcfT                                         ! temperature response of Michaelis-Menten-Constant Kc []
real                    :: KofT                                         ! temperature response of Michaelis-Menten-Constant Ko []
!real                    :: vjfT                                         ! temperature response of vcmax and jmax []
real                    :: VmfT                                         ! temperature response of vcmax []
real                    :: JmfT                                         ! temperature response of jmax []
real                    :: KcM                                          ! Michaelis-Menten-Constant for CO2 [muM]
real                    :: KoM                                          ! Michaelis-Menten-Constant for O2 [muM]
real                    :: Kc                                           ! temperature corrected Michaelis-Menten-Constant for CO2 [mol / m3]
real                    :: Ko                                           ! temperature corrected Michaelis-Menten-Constant for O2 [mol / m3]
!real                    :: Gs                                           ! gammastar [mol / m3]
real, allocatable       :: jmax(:)                                      ! jmax [mol / (kg C * s)]
real, allocatable       :: vcmax(:)                                     ! vcmax [mol / (kg C * s)]
real, allocatable       :: Mrt_tot(:,:)
real, allocatable       :: SS(:)
real, allocatable       :: Als(:,:)
real, allocatable       :: Acs(:,:)
real, allocatable       :: csum(:)
real                    :: VA                                           ! carboxylation capacity [mol /(m2 * s)]
real                    :: Je                                           ! electron flux  [mol / (m2 * s)]

! Respiration:

real, allocatable       :: Rspec(:)                                        ! respiration per biomass [mol / (kg * s)]

! Lichen activity:

real, allocatable       :: act(:)                                       ! activity of lichen (depends on water status) [] 

! growth

real, allocatable       :: expansion(:)                                 ! cover expansion [ ]
real, allocatable       :: wgtspec(:)                                   ! species weight [ ]

! cover

real                    :: cweight                                      ! cover weight
real                    :: lweight
real                    :: csumL                                        ! cover sum of light BSC
real                    :: csumD                                        ! cover sum of dark BSC
real                    :: csumC                                        ! cover sum of chlorolichen BSC
real                    :: csumM                                        ! cover sum of moss BSC

!***********************************************************************
! Lichen state variables
!***********************************************************************

real, allocatable    :: klife(:,:,:)                             ! switch to check if lichen is alive

real, allocatable       :: rH2Ol_t(:,:,:)                           ! thallus water content [m3 H2O / m2 T]    !NOT REQ.
real, allocatable       :: rH2Ol_b(:,:,:)                           ! below thallus water content [m3 H2O / m2 T] !NOT REQ.
!real, allocatable       :: act_time(:,:,:)                              ! active time [s]
real, allocatable       :: gpp0(:,:)                              ! GPP [mol C / (m2 T * s)]
real, allocatable       :: gpp(:,:,:)                              ! GPP [mol C / (m2 T * s)]
real, allocatable       :: npp0(:,:)                              ! NPP [mol C / (m2 T * s)]
!real, allocatable       :: Bl(:,:,:)                              ! leaf biomass [mol C / (m2 T * s)]
real, allocatable       :: bsum(:,:)
!real, allocatable       :: Br(:,:,:,:)                              ! root biomass [mol C / (m2 T * s)]
real, allocatable       :: Br(:,:,:)                              ! root biomass [mol C / (m2 T * s)]
real, allocatable       :: Bl(:,:)                              ! leaf biomass [mol C / (m2 T * s)]



real, allocatable       :: areaTH_s(:,:,:)                          ! lichen area per available area (cover) [m2 T / m2 V] 
real, allocatable       :: area_s(:,:,:)                          ! lichen area per available area (cover) [m2 T / m2 V] 

real, allocatable       :: xT_g(:,:,:)                                ! soil temperature [K]

real, allocatable       :: netgrowth(:,:,:)                         ! net growth [kg C / m2 G]

!***********************************************************************
! Lichen fluxes
!***********************************************************************

! Carbon fluxes

real, allocatable       :: fCO2gc_L(:)                                  ! light-limited CO2 assimilation rate [mol / (m2 T * s)]
real, allocatable       :: fCO2gc_W(:)                                  ! water-limited CO2 assimilation rate [mol / (m2 T * s)]
real, allocatable       :: fCO2gc(:,:)                                    ! NPP [mol / (m2 T * s)]
real, allocatable       :: fCO2nc(:,:)                                    ! GPP [mol / (m2 T * s)]
real, allocatable       :: fCcg_M(:)                                    ! maintenance respiration [mol / (m2 T * s)]
real, allocatable       :: fCcg_G(:)                                    ! growth respiration [mol / (m2 T * s)]
real, allocatable       :: fCcb(:)                                      ! Net Primary Productivity [mol / (m2 T * s)]
real, allocatable       :: fCbo(:)                                      ! litterfall [mol / (m2 T * s)]

! Water fluxes

!real, allocatable       :: fH2Ol_xd1(:,:,:)                                  ! runoff [m3 / (m2 T * s)]
real, allocatable       :: fH2Ol_xd(:)                                  !runoff from lycom
real                    :: fH2Ol_xd_land =0.0                                !runoff from lycomland
real, allocatable       :: fH2Olg_xu(:)                                 ! evaporation from the thallus [m3 / (m2 T * s)]
real, allocatable       :: fH2Ol_bx(:)                                  ! bark water uptake [m3 / (m2 T * s)]
real, allocatable       :: fH2Ol_bd(:)                                  ! bark water overflow [m3 / (m2 T * s)]

! Heat fluxes

real, allocatable       :: fQ_tg(:)                                     ! Ground heat flux [W / m2]
real, allocatable       :: fQ_ta_L(:)                                   ! Latent heat flux [W / m2]
real, allocatable       :: fQ_ta_S(:)                                   ! Sensible heat flux [W / m2]
real, allocatable       :: fRAD_H(:)                                    ! net radiation [W / m2]

!***********************************************************************
! Averaged variables
!***********************************************************************

! Species-averaged variables

real                    :: as_rCO2d=0.0
real                    :: as_sCO2d=0.0
real                    :: as_Lai =0.0
real                    :: as_rCb=0.0
real                    :: as_rH2Ol_t=0.0
real                    :: as_rmaxH2Ol_t=0.0
real                    :: as_areaTH_s=0.0
real                    :: as_area_s=0.0
real                    :: as_lai_s=0.0
real                    :: as_act=0.0
real                    :: as_fCO2gc=0.0
real                    :: as_fCcg=0.0
real                    :: as_fCcb=0.0
real                    :: as_fCcb_l=0.0
real                    :: as_fCcb_c=0.0
real                    :: as_fCbo=0.0
real                    :: as_fH2Ogl_ut=0.0
real                    :: as_fH2Olg_tu=0.0
real                    :: as_fH2Ol_lsat=0.0                   
real                    :: as_fH2Ol_bsat=0.0                   
real                    :: as_fH2Ol_runoff_l=0.0
real                    :: as_fCc_npp=0.0
real                    :: as_fCc_gpp=0.0
real                    :: as_fH2Ol_bx=0.0
real                    :: as_fH2Ol_bd=0.0
real                    :: as_Ts=0.0
real                    :: as_Tg=0.0
real                    :: as_H=0.0
real                    :: as_G=0.0
real                    :: as_E=0.0
real                    :: as_C=0.0
real                    :: as_EB=0.0
real                    :: as_ET=0.0

! land surface
real                    :: a0_Ts

! habitat-averaged variables

real, allocatable       :: ah_areaTH_s(:)
real, allocatable       :: ah_area_s(:)
real, allocatable       :: ah_rH2Ol(:)
real, allocatable       :: ah_rmaxH2Ol(:)
real, allocatable       :: ah_act(:)
real, allocatable       :: ah_fH2Ol_xd(:)
real, allocatable       :: ah_fH2Ol_runoff_l(:)
real, allocatable       :: ah_fCc_gpp(:)
real, allocatable       :: ah_fCc_npp(:)
real, allocatable       :: ah_fH2Ol_lsat(:)              
real, allocatable       :: ah_fH2Ol_bsat(:)            
real, allocatable       :: ah_fH2Ogl_ux(:)
real, allocatable       :: ah_fH2Olg_xu(:)
real, allocatable       :: ah_fH2Ol_bx
real, allocatable       :: ah_fH2Ol_bd

real, allocatable       :: ah_Ts(:)
real, allocatable       :: ah_Tg(:)
real, allocatable       :: ah_H(:)
real, allocatable       :: ah_G(:)
real, allocatable       :: ah_E(:)
real, allocatable       :: ah_C(:)
real, allocatable       :: ah_EB(:)

real, allocatable       :: ah_rCO2d(:)
real, allocatable       :: ah_sCO2d(:)
real, allocatable       :: ah_Lai(:)
real, allocatable       :: ah_rCb(:)
real, allocatable       :: ah_fCO2gc(:)
real, allocatable       :: ah_fCcg(:)
real, allocatable       :: ah_fCcb(:)
real, allocatable       :: ah_fCcb_l(:)
real, allocatable       :: ah_fCcb_c(:)
real, allocatable       :: ah_fCbo(:)

! level-averaged variables

real, allocatable       :: av_areaTH_s(:)
real, allocatable       :: av_area_s(:)
real, allocatable       :: av_rH2Ol(:)
real, allocatable       :: av_rmaxH2Ol(:)
real, allocatable       :: av_act(:)
real, allocatable       :: av_fH2Ol_ux(:)
real, allocatable       :: av_fH2Ol_lsat(:)
real, allocatable       :: av_fH2Ol_runoff_l(:)
real, allocatable       :: av_fCc_npp(:)
real, allocatable       :: av_fCc_gpp(:)
real, allocatable       :: av_fH2Ol_bsat(:)
real, allocatable       :: av_fH2Ol_xd(:)
real, allocatable       :: av_fH2Ogl_ux(:)
real, allocatable       :: av_fH2Olg_xu(:)
real, allocatable       :: av_fH2Ol_bx(:)

real                    :: av_fH2Ol_bd
real                    :: av_fH2Ol_ug2

real, allocatable       :: av_Ts(:)
real, allocatable       :: av_H(:)
real, allocatable       :: av_E(:)
real, allocatable       :: av_C(:)
real, allocatable       :: av_EB(:)

real, allocatable       :: av_rCO2d(:)
real                   :: av_sCO2d
real                   :: av_Lai
real, allocatable       :: av_rCb(:)
real, allocatable       :: av_fCO2gc(:)
real, allocatable       :: av_fCcg(:)
real, allocatable       :: av_fCcb(:)
real, allocatable       :: av_fCcb_l(:)
real, allocatable       :: av_fCcb_c(:)
real, allocatable       :: av_fCbo(:)

real                    :: av_Tg
real                    :: av_G
real                    :: av_rH2Ol_g1
real                    :: av_rH2Ol_g2
real                    :: av_fH2Ol_ug
real                    :: av_fH2Ol_go
real                    :: av_fH2Ol_gb
real                    :: av_fH2Olg_ga

! tile-averaged variables

real, allocatable       :: at_fH2Ol_ux(:,:)
real, allocatable       :: at_fH2Ol_xd(:,:)
real, allocatable       :: at_fH2Ol_lsat(:,:)
real, allocatable       :: at_fH2Ol_runoff_l(:,:)
real, allocatable       :: at_fCc_gpp(:,:)
real, allocatable       :: at_fCc_npp(:,:)
real, allocatable       :: at_fH2Ol_bsat(:,:)
real, allocatable       :: at_areaTH_s(:,:)
real, allocatable       :: at_area_s(:,:)
real, allocatable       :: at_rH2Ol(:,:)
real, allocatable       :: at_rmaxH2Ol(:,:)
real, allocatable       :: at_act(:,:)
real, allocatable       :: at_rCO2d(:,:)
real, allocatable       :: at_sCO2d(:)
real, allocatable       :: at_Lai(:)
real, allocatable       :: at_rCb(:,:)
real, allocatable       :: at_fCO2gc(:,:)
real, allocatable       :: at_fCcg(:,:)
real, allocatable       :: at_fCcb(:,:)
real, allocatable       :: at_fCcb_l(:,:)
real, allocatable       :: at_fCcb_c(:,:)
real, allocatable       :: at_fCbo(:,:)

real, allocatable       :: at_fH2Ogl_ux(:,:)
real, allocatable       :: at_fH2Olg_xu(:,:)
real, allocatable       :: at_fH2Ol_bx(:,:)

real, allocatable       :: at_fH2Ol_bd(:)
real, allocatable       :: at_fH2Ol_ug2(:)

real, allocatable       :: at_Ts(:,:)
real, allocatable       :: at_H(:,:)
real, allocatable       :: at_E(:,:)
real, allocatable       :: at_C(:,:)
real, allocatable       :: at_EB(:,:)

real, allocatable       :: at_Tg(:)
real, allocatable       :: at_G(:)

real, allocatable       :: at_rH2Ol_g1(:)
real, allocatable       :: at_rH2Ol_g2(:)

real, allocatable       :: at_fH2Ol_ug(:)
real, allocatable       :: at_fH2Ol_go(:)
real, allocatable       :: at_fH2Ol_gb(:)
real, allocatable       :: at_fH2Olg_ga(:)

! tile-averaged variables for next time step

real, allocatable       :: atN_fH2Ol_xd(:,:,:)
real, allocatable       :: atN_fH2Ol_bd(:,:)

!***********************************************************************
! Output variables
!***********************************************************************

! grid cell-averaged variables

real, allocatable       :: ag_fH2Ol_ux(:,:)
real, allocatable       :: ag_fH2Ol_xd(:,:)
real, allocatable       :: ag_area_s(:,:)
real, allocatable       :: ag_areaTH_s(:,:)
real, allocatable       :: ag_rH2Ol(:,:)
real, allocatable       :: ag_rmaxH2Ol(:,:)
real, allocatable       :: ag_act(:,:)
real, allocatable       :: ag_rCO2d(:,:)
real, allocatable       :: ag_sCO2d(:)
real, allocatable       :: ag_Lai(:)
real, allocatable       :: ag_rCb(:,:)
real, allocatable       :: ag_fCO2gc(:,:)
real, allocatable       :: ag_fCcg(:,:)
real, allocatable       :: ag_fCcb(:,:)
real, allocatable       :: ag_fCcb_l(:,:)
real, allocatable       :: ag_fCcb_c(:,:)
real, allocatable       :: ag_fCbo(:,:)
real, allocatable       :: ag_fH2Ol_lsat(:,:)
real, allocatable       :: ag_fH2Ol_runoff_l(:,:)
real, allocatable       :: ag_fCc_gpp(:,:)
real, allocatable       :: ag_fCc_npp(:,:)
real, allocatable       :: ag_fH2Ol_bsat(:,:)
real, allocatable       :: ag_fH2Ogl_ux(:,:)
real, allocatable       :: ag_fH2Olg_xu(:,:)
real, allocatable       :: ag_fH2Ol_bx(:,:)
real, allocatable       :: ag_fH2Ol_bd(:)
real, allocatable       :: ag_Ts(:,:)
real, allocatable       :: ag_H(:,:)
real, allocatable       :: ag_E(:,:)
real, allocatable       :: ag_C(:,:)
real, allocatable       :: ag_EB(:,:)

real, allocatable       :: ag_Tg(:)
real, allocatable       :: ag_G(:)

real, allocatable       :: ag_rH2Ol_g1(:)
real, allocatable       :: ag_rH2Ol_g2(:)

real, allocatable       :: ag_fH2Ol_ug(:)
real, allocatable       :: ag_fH2Ol_go(:)
real, allocatable       :: ag_fH2Ol_gb(:)
real, allocatable       :: ag_fH2Olg_ga(:)
real, allocatable       :: ag_fH2Ol_ug2(:)

real, allocatable       :: ag_rH2Os_g(:)
real, allocatable       :: ag_fH2Osl_g(:)

real, allocatable       :: ag_fH2Os_ad(:)
real, allocatable       :: ag_xT_a(:)
real, allocatable       :: ag_fH2Ol_ad(:)
real, allocatable       :: ag_fRADs(:)

! Species parameter variables

real, allocatable       :: count_spec(:)
real, allocatable       :: count_spec_h(:,:)
real, allocatable       :: vec_o_avg(:,:,:)

!***********************************************************************
! Input files
!***********************************************************************

! file names

character(len=*),parameter :: snamelist         = 'lycom_namelist'      ! namelist file
character(len=*),parameter :: sspecpar          = 'lycom_specpar'       ! species parameter file

character(len=*),parameter :: slaifile          = 'LAI'                 ! LAI file
character(len=*),parameter :: ssaifile          = 'SAI'                 ! SAI file 
character(len=*),parameter :: sssafile          = 'SSA'                 ! soil surface area file 
character(len=*),parameter :: sbiomefile        = 'biome'               ! Biome file 
character(len=*),parameter :: sETrmaxfile       = 'ETrmax'
character(len=*),parameter :: sETdavgfile       = 'ETdavg'

character(len=*),parameter :: sstatus           = 'lycom_status'        ! simulation status file

character(len=*),parameter :: ssradfile         = 'srad'                ! shortwave radiation
character(len=*),parameter :: slradfile         = 'lrad'                ! longwave radiation
character(len=*),parameter :: srainfile         = 'rain'                ! rainfall
character(len=*),parameter :: sssnowfile         = 'snow'                ! snow
character(len=*),parameter :: ssrhumfile         = 'rhum'                ! relative humidity
character(len=*),parameter :: stairfile         = 'tair'                ! air temperature
character(len=*),parameter :: swindfile         = 'wind'                ! wind

! file units

integer, parameter      :: knamelist            = 100                   ! namelist file 
integer, parameter      :: kspecpar             = 101                   ! species parameter file 

integer, parameter      :: klaifile             = 110                   ! LAI file
integer, parameter      :: ksaifile             = 111                   ! SAI file 
integer, parameter      :: kssafile             = 112                   ! soil surface area file 
integer, parameter      :: kbiomefile           = 113                   ! Biome file 

integer, parameter      :: kstatus              = 199                   ! simulation status file

integer, parameter      :: ksradfile            = 200                   ! shortwave radiation
integer, parameter      :: klradfile            = 201                   ! longwave radiation
integer, parameter      :: krainfile            = 202                   ! rainfall
integer, parameter      :: ksnowfile            = 203                   ! snow
integer, parameter      :: krhumfile            = 204                   ! relative humidity
integer, parameter      :: ktairfile            = 205                   ! air temperature
integer, parameter      :: kwindfile            = 206                   ! wind

!***********************************************************************
! Output files
!***********************************************************************

! global file names

character(len=*),parameter :: sfile_outputG     = 'lycom_output.nc'
character(len=*),parameter :: sfile_outputGS    = 'lycom_outputS.nc'

! local file names

character(len=*),parameter :: sfile_count_spec  = 'Number_of_strategies'

character(len=*),parameter :: sfile_rCO2d       = 'CO2_pore_space'
character(len=*),parameter :: sfile_sCO2d       = 'Soil CO2 content'
character(len=*),parameter :: sfile_Lai         = 'Leaf Area Index'
character(len=*),parameter :: sfile_rCb         = 'Biomass'
character(len=*),parameter :: sfile_fH2Ol_lsat  = 'soil Layer relative saturation' 
character(len=*),parameter :: sfile_fH2Ol_bsat  = 'soil bucket relative saturation'
character(len=*),parameter :: sfile_fH2Ol_runoff_l= 'Actual Runoff'
character(len=*),parameter :: sfile_fCc_gpp = 'Lycophyte GPP'
character(len=*),parameter :: sfile_fCc_npp = 'Lycophyte NPP'
character(len=*),parameter :: sfile_rH2Ol       = 'Water_saturation'
character(len=*),parameter :: sfile_rH2Ol_g1    = 'TopSoil_water_saturation'    !G
character(len=*),parameter :: sfile_rH2Ol_g2    = 'DeepSoil_water_saturation'   !G
character(len=*),parameter :: sfile_rmaxH2Ol    = 'Water_storage_capacity'
character(len=*),parameter :: sfile_rH2Os_g     = 'Snow_cover'                  !G

character(len=*),parameter :: sfile_area        = 'Cover_fraction'

character(len=*),parameter :: sfile_act         = "Active_time_fraction"

character(len=*),parameter :: sfile_fCO2gc      = 'GPP'
character(len=*),parameter :: sfile_fCcg        = 'Respiration'
character(len=*),parameter :: sfile_fCcb        = 'NPP'
character(len=*),parameter :: sfile_fCcb_l      = 'NPPl'
character(len=*),parameter :: sfile_fCcb_c      = 'NPPc'
character(len=*),parameter :: sfile_fCbo        = 'Biomass_loss'

character(len=*),parameter :: sfile_fH2Ol_ux    = 'Water_input'
character(len=*),parameter :: sfile_fH2Ol_xd    = 'Layer_overflow'
character(len=*),parameter :: sfile_fH2Ogl_ux   = 'Dew_input'
character(len=*),parameter :: sfile_fH2Olg_xu   = 'Evaporation'
character(len=*),parameter :: sfile_fH2Ol_bx    = 'Layer_uptake'
character(len=*),parameter :: sfile_fH2Ol_bd    = 'Stemflow'
character(len=*),parameter :: sfile_fH2Ol_ug    = 'Soilwater_input'             !G
character(len=*),parameter :: sfile_fH2Ol_ug2   = 'Throughfall'                 !G
character(len=*),parameter :: sfile_fH2Ol_go    = 'Surface_runoff'              !G
character(len=*),parameter :: sfile_fH2Ol_gb    = 'Baseflow'                    !G
character(len=*),parameter :: sfile_fH2Olg_ga   = 'Evapotranspiration'          !G
character(len=*),parameter :: sfile_fH2Osl_g    = 'Snowmelt'                    !G
character(len=*),parameter :: sfile_fH2Os_ad    = 'Snowfall'                    !I
character(len=*),parameter :: sfile_fH2Ol_ad    = 'Rainfall'                    !I

character(len=*),parameter :: sfile_Ts          = 'Surface_temperature'
character(len=*),parameter :: sfile_xT_a        = 'Air_temperature'             !I
character(len=*),parameter :: sfile_Tg          = 'Soil_temperature'            !G
character(len=*),parameter :: sfile_H           = 'Net_enthalpy_flux'
character(len=*),parameter :: sfile_G           = 'Ground_heat_flux'            !G
character(len=*),parameter :: sfile_E           = 'Latent_heat_flux'
character(len=*),parameter :: sfile_C           = 'Sensible_heat_flux'
character(len=*),parameter :: sfile_EB          = 'Energy_balance'
character(len=*),parameter :: sfile_fRADs       = 'Shortwave_radiation'         !I

character(len=*),parameter :: sfile_survspec    = 'Surviving_strategies'
                                                   
character(len=*),parameter :: sfile_Albedo      = 'Albedo'
character(len=*),parameter :: sfile_Height      = 'Height'
character(len=*),parameter :: sfile_Porosity    = 'Porosity'
character(len=*),parameter :: sfile_FracAir     = 'FracAir'
character(len=*),parameter :: sfile_DCO2B       = 'DCO2B'
character(len=*),parameter :: sfile_MolVcmax    = 'MolVcmax'
character(len=*),parameter :: sfile_MolVomax    = 'MolVomax'
character(len=*),parameter :: sfile_PhsynCap    = 'PhsynCap'
character(len=*),parameter :: sfile_Hydrophob   = 'Hydrophob'
character(len=*),parameter :: sfile_Topt        = 'Topt'
character(len=*),parameter :: sfile_Q10value    = 'Q10value'
character(len=*),parameter :: sfile_CCM         = 'CCM'

character(len=*),parameter :: sfile_canopy      = 'Canopy_'
character(len=*),parameter :: sfile_ground      = 'Ground_'

character(len=*),parameter :: sfile_forLeaves   = 'Leaves_'
character(len=*),parameter :: sfile_forStems    = 'Stems_'
character(len=*),parameter :: sfile_forFloor    = 'Floor_'
character(len=*),parameter :: sfile_grass       = 'Grass_'

character(len=*),parameter :: sfile_weightS     = 'WeightS'

! file units

integer, parameter      :: kfile_rCO2d          = 100 
integer, parameter      :: kfile_sCO2d          = 101
integer, parameter      :: kfile_Lai            = 102                                                      
integer, parameter      :: kfile_rCb            = 200 
                                                      
integer, parameter      :: kfile_rH2Ol          = 300 
integer, parameter      :: kfile_rH2Ol_g1       = 301   !G
integer, parameter      :: kfile_rH2Ol_g2       = 302   !G
integer, parameter      :: kfile_rmaxH2Ol       = 310     
integer, parameter      :: kfile_rH2Os_g        = 320   !G
                                                          
integer, parameter      :: kfile_area           = 400     
                                                          
integer, parameter      :: kfile_act            = 500     
                                                          
integer, parameter      :: kfile_fCO2gc         = 600     
integer, parameter      :: kfile_fCcg           = 610     
integer, parameter      :: kfile_fCcb           = 620     
integer, parameter      :: kfile_fCcb_l         = 621     
integer, parameter      :: kfile_fCcb_c         = 622     
integer, parameter      :: kfile_fCbo           = 630     
                                                          
integer, parameter      :: kfile_fH2Ol_ux       = 700     
integer, parameter      :: kfile_fH2Ol_ug       = 701   !G
integer, parameter      :: kfile_fH2Ol_ug2      = 702   !G
integer, parameter      :: kfile_fH2Ol_xd       = 710     
integer, parameter      :: kfile_fH2Ol_go       = 711   !G
integer, parameter      :: kfile_fH2Ol_gb       = 712   !G
integer, parameter      :: kfile_fH2Ogl_ux      = 720     
integer, parameter      :: kfile_fH2Olg_xu      = 730     
integer, parameter      :: kfile_fH2Olg_ga      = 731   !G
integer, parameter      :: kfile_fH2Osl_g       = 740   !G
integer, parameter      :: kfile_fH2Ol_ad       = 750   !I
integer, parameter      :: kfile_fH2Os_ad       = 751   !I
integer, parameter      :: kfile_fH2Ol_bx       = 760
integer, parameter      :: kfile_fH2Ol_lsat     = 761
integer, parameter      :: kfile_fH2Ol_bsat     = 762
integer, parameter      :: kfile_fH2Ol_runoff_l      = 763
integer, parameter      :: kfile_fCc_gpp      = 764
integer, parameter      :: kfile_fCc_npp      = 765
integer, parameter      :: kfile_fH2Ol_bd       = 770
                                                          
integer, parameter      :: kfile_Ts             = 800     
integer, parameter      :: kfile_xT_a           = 801   !I
integer, parameter      :: kfile_Tg             = 820   !G
integer, parameter      :: kfile_H              = 830     
integer, parameter      :: kfile_G              = 840   !G
integer, parameter      :: kfile_E              = 850     
integer, parameter      :: kfile_C              = 860     
integer, parameter      :: kfile_EB             = 870     
integer, parameter      :: kfile_fRADs          = 880   !I

integer, parameter      :: kfile_weight         = 7000

integer, parameter      :: kfile_count_spec     = 8000
integer, parameter      :: kfile_survspec       = 8999

! output directories

character(len=*),parameter :: fluxdir           = 'output_fluxes/'
character(len=*),parameter :: specdir           = 'output_strategies/'


!###############################################################


! OPTIONAL PARAMETERS AND VARIABLES
!---------------------------------------------------------------

! Species parameters

real, allocatable       :: o_LC(:)                                      ! light crust []
real, allocatable       :: o_DC(:)                                      ! dark crust []
real, allocatable       :: o_CC(:)                                      ! chlorolichen crust []
real, allocatable       :: o_MC(:)                                      ! moss crust []

! NO/HONO Emissions

real, allocatable       :: NO_HONO_fSat(:,:)                            ! NO / HONO emissions (OPTIONAL) [ ng / m2 / s ]
real, allocatable       :: fNO_N(:)                                     ! NO Emissions [ng / m2 T / s ]
real, allocatable       :: fHONO_N(:)                                   ! HONO Emissions [ng / m2 T / s ]

! Accumulation variables

real, allocatable       :: a_areaTHLC_g(:,:,:)
real, allocatable       :: a_areaTHDC_g(:,:,:)
real, allocatable       :: a_areaTHCC_g(:,:,:)
real, allocatable       :: a_areaTHMC_g(:,:,:)
real, allocatable       :: a_rH2OlLC(:,:,:)
real, allocatable       :: a_rH2OlDC(:,:,:)
real, allocatable       :: a_rH2OlCC(:,:,:)
real, allocatable       :: a_rH2OlMC(:,:,:)
real, allocatable       :: a_TsLC(:,:,:)
real, allocatable       :: a_TsDC(:,:,:)
real, allocatable       :: a_TsCC(:,:,:)
real, allocatable       :: a_TsMC(:,:,:)
real, allocatable       :: a_fCcbLC(:,:,:)
real, allocatable       :: a_fCcbDC(:,:,:)
real, allocatable       :: a_fCcbCC(:,:,:)
real, allocatable       :: a_fCcbMC(:,:,:)

real, allocatable       :: a_fNO_N(:,:,:)
real, allocatable       :: a_fHONO_N(:,:,:)

! Averaged accumulation variables

real, allocatable       :: a_MareaTHLC_g_S(:)
real, allocatable       :: a_MareaTHDC_g_S(:)
real, allocatable       :: a_MareaTHCC_g_S(:)
real, allocatable       :: a_MareaTHMC_g_S(:)
real, allocatable       :: a_MfNO_N(:)
real, allocatable       :: a_MfHONO_N(:)
real, allocatable       :: a_MrH2OlLC_S(:)
real, allocatable       :: a_MrH2OlDC_S(:)
real, allocatable       :: a_MrH2OlCC_S(:)
real, allocatable       :: a_MrH2OlMC_S(:)
real, allocatable       :: a_MTsLC_S(:)
real, allocatable       :: a_MTsDC_S(:)
real, allocatable       :: a_MTsCC_S(:)
real, allocatable       :: a_MTsMC_S(:)
real, allocatable       :: a_MfCcbLC_S(:)
real, allocatable       :: a_MfCcbDC_S(:)
real, allocatable       :: a_MfCcbCC_S(:)
real, allocatable       :: a_MfCcbMC_S(:)

! input file names

character(len=*),parameter :: sNOHONO           = 'NO_HONO_fSatBins'    ! NO / HONO emissions file (OPTIONAL)

! local file names

character(len=*),parameter :: sfile_MareaLCgS   = "Light_crust_coverage"
character(len=*),parameter :: sfile_MareaDCgS   = "Dark_crust_coverage"
character(len=*),parameter :: sfile_MareaCCgS   = "Chlorolichen_crust_coverage"
character(len=*),parameter :: sfile_MareaMCgS   = "Moss_crust_coverage"
character(len=*),parameter :: sfile_MfNO_N      = "NO-N_Emissions"
character(len=*),parameter :: sfile_MfHONO_N    = "HONO-N_Emissions"
character(len=*),parameter :: sfile_MrH2OlLC_S  = "Light_crust_moisture"
character(len=*),parameter :: sfile_MrH2OlDC_S  = "Dark_crust_moisture"
character(len=*),parameter :: sfile_MrH2OlCC_S  = "Chlorolichen_crust_moisture"
character(len=*),parameter :: sfile_MrH2OlMC_S  = "Moss_crust_moisture"
character(len=*),parameter :: sfile_MTsLC_S     = "Light_crust_temperature"
character(len=*),parameter :: sfile_MTsDC_S     = "Dark_crust_temperature"
character(len=*),parameter :: sfile_MTsCC_S     = "Chlorolichen_crust_temperature"
character(len=*),parameter :: sfile_MTsMC_S     = "Moss_crust_temperature"
character(len=*),parameter :: sfile_MfCcbLC_S   = "Light_crust_NPP"
character(len=*),parameter :: sfile_MfCcbDC_S   = "Dark_crust_NPP"
character(len=*),parameter :: sfile_MfCcbCC_S   = "Chlorolichen_crust_NPP"
character(len=*),parameter :: sfile_MfCcbMC_S   = "Moss_crust_NPP"

! file units

integer, parameter      :: kNOHONO              = 900                   ! NO / HONO emissions file (OPTIONAL)

integer, parameter      :: kfile_MareaLCgS      = 9001
integer, parameter      :: kfile_MareaDCgS      = 9002
integer, parameter      :: kfile_MareaCCgS      = 9003
integer, parameter      :: kfile_MareaMCgS      = 9004
integer, parameter      :: kfile_MrH2OlLC_S     = 9011
integer, parameter      :: kfile_MrH2OlDC_S     = 9012
integer, parameter      :: kfile_MrH2OlCC_S     = 9013
integer, parameter      :: kfile_MrH2OlMC_S     = 9014
integer, parameter      :: kfile_MTsLC_S        = 9021
integer, parameter      :: kfile_MTsDC_S        = 9022
integer, parameter      :: kfile_MTsCC_S        = 9023
integer, parameter      :: kfile_MTsMC_S        = 9024
integer, parameter      :: kfile_MfCcbLC_S      = 9031
integer, parameter      :: kfile_MfCcbDC_S      = 9032
integer, parameter      :: kfile_MfCcbCC_S      = 9033
integer, parameter      :: kfile_MfCcbMC_S      = 9034

integer, parameter      :: kfile_MfNO_N         = 9100
integer, parameter      :: kfile_MfHONO_N       = 9101


!contains
!
!
!! CHECK STATUS
!!---------------------------------------------------------------
!
!subroutine check(status)
!use netcdf
!
!  integer, intent(in) :: status
!                     
!  if(status /= nf90_noerr) then
!    print *, trim(nf90_strerror(status))
!    stop 2
!  end if
!
!end subroutine check

end module lycom_par

