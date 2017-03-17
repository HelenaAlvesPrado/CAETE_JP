module global_pars
  implicit none
  integer,parameter :: l1 = kind(.TRUE.)
  integer,parameter :: i1 = selected_int_kind(2)
  integer,parameter :: i4 = kind(0)
  integer,parameter :: r4 = kind(0.0)
  integer,parameter :: r8 = kind(0.0D0)
  integer,parameter :: rbig = selected_real_kind(16,300)

  real(kind=r4),parameter :: H = 1.0                         ! soil layer thickness (meters)
  real(kind=r4),parameter :: DIFFU = 1.036800e14 ! soil thermal diffusivity (m2/mes)
  real(kind=r4),parameter :: TAU = 4.822530864197531e-15  ! e-folding times (months) 
  real(kind=r4),parameter :: rcmax = 2000.0
  real(kind=r4),parameter :: rcmin = 110.0
  real(kind=r4),parameter :: ca = 363.0 ! ppmv - atm[CO2]
  real(kind=r4),parameter :: wmax = 500.0
  
  real(kind=r8),parameter :: csru = 0.5_r8
  real(kind=r8),parameter :: alfm = 1.391_r8
  real(kind=r8),parameter :: gm = 3.26_r8  !(*86400 transform s/mm to dia/mm)

  real(kind=r8),parameter :: ncl = 1.0_r8/29.0_r8              !(gN/KgC) from lpj3 
  real(kind=r8),parameter :: ncf = 1.0_r8/29.0_r8              !(gN/KgC)
  real(kind=r8),parameter :: ncs = 1.0_r8/330.0_r8             !(gN/KgC)

  integer(kind=i4) :: ndmonth(12)       !Number of months
  data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month

  integer(kind=i4),parameter :: npls = 20
  integer(kind=i4),parameter :: ntimes = 12
  integer(kind=i4),parameter :: ntraits = 8

  integer(kind=i4),parameter :: nx = 720, ny = 360
  
end module global_pars

module photo_par
  use global_pars, only : r8
  implicit none
  ! setting some parameter values
  real(kind=r8), parameter ::       &
       a   = 0.8300_r8    ,&          !Photosynthesis co-limitation coefficient
       a2  = 0.930_r8     ,&          !Photosynthesis co-limitation coefficient
       p3  = 21200.0_r8   ,&          !Atmospheric oxygen concentration (Pa)
       p4  = 0.080_r8     ,&          !Quantum efficiency (mol electrons/Ein)
       p5  = 0.150_r8     ,&          !Light scattering rate
       p6  = 2.0_r8       ,&          !Parameter for jl
       p7  = 0.50_r8      ,&          !Ratio of light limited photosynthesis to Rubisco carboxylation
       p8  = 5200.0_r8    ,&          !Photo-respiration compensation point
       p9  = 0.570_r8     ,&          !Photosynthesis co-limitation coefficient
       p10 = 0.100_r8     ,&          !Q10 function
       p11 = 25.0_r8      ,&          !Q10 function reference temperature (oC)
       p12 = 30.0_r8      ,&          !Michaelis-Menten constant for CO2 (Pa)
       p13 = 2.100_r8     ,&          !Michaelis-Menten constant for CO2
       p14 = 30000.0_r8   ,&          !Michaelis-Menten constant for O2 (Pa)
       p15 = 1.20_r8      ,&          !Michaelis-Menten constant for O2
       p19 = 0.90_r8      ,&          !Maximum ratio of internal to external CO2
       p20 = 0.10_r8      ,&          !Critical humidity deficit (kg/kg)
       p22 = 2.0_r8       ,&          !Rubisco carboxylation rate
       p23 = 0.30_r8      ,&          !Rubisco carboxylation rate
       p24 = 36.0_r8      ,&          !Rubisco carboxylation rate (oC)
       p25 = 0.000008_r8  ,&          !Maximum gross photosynthesis rate (molCO2/m2/s)
       p26 = 0.50_r8      ,&          !light extinction coefficient for IPAR/sun (0.5/sen90)
       p27 = 1.50_r8      ,&          !light extinction coefficient for IPAR/shade (0.5/sen20)
       p28 = 0.500_r8     ,&          !Soil moisture at field capacity
       p29 = 0.205_r8     ,&          !Soil moisture at wilting point
       p30 = 0.015_r8     ,&          !Ratio of respiration to Rubisco carboxylation rates
       p31 = 3.850_r8     ,&            !Whole plant to leaf respiration ratio
       p32 = 2.00_r8      ,&
       p33 = 0.10_r8      ,&
       p34 = 0.30_r8      ,&
       p35 = 0.05_r8      ,&
       p36 = 0.25_r8
    
end module photo_par
