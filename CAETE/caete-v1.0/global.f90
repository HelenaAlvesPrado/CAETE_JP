module global_pars
  implicit none
  integer,parameter :: l1 = kind(.TRUE.)
  integer,parameter :: i4 = kind(0)
  integer,parameter :: r4 = kind(0.0)
  integer,parameter :: r8 = kind(0.0D0)
  integer,parameter :: rbig = selected_real_kind(16,300)

  real(kind=r4),parameter :: H = 1.0                         ! soil layer thickness (meters)
  real(kind=r4),parameter :: DIFFU = 4.e7 * (30.0 * 86400.0) ! soil thermal diffusivity (m2/mes)
  real(kind=r4),parameter :: TAU = (H ** 2) / (2.0 * DIFFU)  ! e-folding times (months) 
  real(kind=r4),parameter :: rcmax = 2000.0
  real(kind=r4),parameter :: rcmin = 110.0
  real(kind=r4),parameter :: ca = 363.0 ! ppmv - atm[CO2]
  real(kind=r4),parameter :: wmax = 500.0
  
  real(kind=r8),parameter :: csru = 0.5
  real(kind=r8),parameter :: alfm = 1.391
  real(kind=r8),parameter :: gm = 3.26  !(*86400 transform s/mm to dia/mm)

  real(kind=r8),parameter :: ncl = 1./29.              !(gN/KgC) from lpj3 
  real(kind=r8),parameter :: ncf = 1./29.              !(gN/KgC)
  real(kind=r8),parameter :: ncs = 1./330.             !(gN/KgC)

  integer(kind=i4) :: ndmonth(12)       !Number of months
  data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/ !Number of days for each month

  integer(kind=i4),parameter :: npls = 20
  integer(kind=i4),parameter :: ntimes = 12
  integer(kind=i4),parameter :: ntraits = 8
  
end module global_pars

module photo_par
  use global_pars, only : r8
  implicit none
  ! setting some parameter values
  real(kind=r8), parameter ::       &
       a   = 0.8300    ,&          !Photosynthesis co-limitation coefficient
       a2  = 0.930     ,&          !Photosynthesis co-limitation coefficient
       p3  = 21200.0   ,&          !Atmospheric oxygen concentration (Pa)
       p4  = 0.080     ,&          !Quantum efficiency (mol electrons/Ein)
       p5  = 0.150     ,&          !Light scattering rate
       p6  = 2.0       ,&          !Parameter for jl
       p7  = 0.50      ,&          !Ratio of light limited photosynthesis to Rubisco carboxylation
       p8  = 5200.0    ,&          !Photo-respiration compensation point
       p9  = 0.570     ,&          !Photosynthesis co-limitation coefficient
       p10 = 0.100     ,&          !Q10 function
       p11 = 25.0      ,&          !Q10 function reference temperature (oC)
       p12 = 30.0      ,&          !Michaelis-Menten constant for CO2 (Pa)
       p13 = 2.100     ,&          !Michaelis-Menten constant for CO2
       p14 = 30000.0   ,&          !Michaelis-Menten constant for O2 (Pa)
       p15 = 1.20      ,&          !Michaelis-Menten constant for O2
       p19 = 0.90      ,&          !Maximum ratio of internal to external CO2
       p20 = 0.10      ,&          !Critical humidity deficit (kg/kg)
       p22 = 2.0       ,&          !Rubisco carboxylation rate
       p23 = 0.30      ,&          !Rubisco carboxylation rate
       p24 = 36.0      ,&          !Rubisco carboxylation rate (oC)
       p25 = 0.000008  ,&          !Maximum gross photosynthesis rate (molCO2/m2/s)
       p26 = 0.50      ,&          !light extinction coefficient for IPAR/sun (0.5/sen90)
       p27 = 1.50      ,&          !light extinction coefficient for IPAR/shade (0.5/sen20)
       p28 = 0.500     ,&          !Soil moisture at field capacity
       p29 = 0.205     ,&          !Soil moisture at wilting point
       p30 = 0.015     ,&          !Ratio of respiration to Rubisco carboxylation rates
       p31 = 3.850     ,&            !Whole plant to leaf respiration ratio
       p32 = 2.00      ,&
       p33 = 0.10      ,&
       p34 = 0.30      ,&
       p35 = 0.05      ,&
       p36 = 0.25
    
end module photo_par
