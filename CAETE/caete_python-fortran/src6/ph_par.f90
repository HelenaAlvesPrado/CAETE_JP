module photo_par
  ! setting some parameter values
  real(kind=r8), parameter::       &
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
       p31 = 3.850                 !Whole plant to leaf respiration ratio
end module photo_par