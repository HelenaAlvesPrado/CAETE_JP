c234567
c=======================================================================
c Subroutine carbon1 calculates photosynthesis, plant respiration and
c net primary productivity. Subroutine carbon2 calculates microbial
c respiration (yet to be improved).
c
c Code written by David Lapola
c Last update: Aug/2007
c
c A "bug" has been found on Sep/2007: when the water soil
c is between 0.205 and 0.265, NPP unrealistically drops to a level
c below those when wsoil is lesser than 0.205 (because of f5).
c
c=======================================================================
c
      subroutine carbon1 (temp,p0,w,wmax,ca,ipar,i6,tsoil,emax, rc2,       !input
     &                     ca2,cf2,cb2,cl1,ca1,cf1,cb1,
     &					   ocp_tbe2,ocp_tbd2,ocp_herb2,
     &                     beta_leaf,beta_awood,beta_bwood,beta_froot,
     &                     ph,ar,nppa,laia,f5,rm,rml,rmf,rms,
     &                     rg,rgl,rgf,rgs)		            !output
c
c i/o variables
      real temp		!mean monthly temperature (oC)
      real p0		!mean surface pressure (hPa)
      real wa,w,wmax	!soil moisture (dimensionless)
      real ca		!Atmospheric CO2 pressure (Pa)
      real ipar		!incident photos. active radiation
c
      real ph		!Canopy gross photosyntheis (kgC/m2/yr)tivity
      real laia		!Leaf area index (m2 leaf/ m2 area)
      real ar		!Autotrophic (plant) respiration (kgC/m2/yr)
      real nppa		!Net primary productivity (kgC/m2/yr)
c
c internal variables
      real vm,mgama,es,rmax,r,ci,a,b,c,a2,b2,c2,
     &     delta,delta2,f1,f1a,f2,f3,f4,f5,k16,k17,rl,rp
      real f10	!CO2 sensitivity

c Rubisco, light and transport limited photosynthesis rate respectively
      real jc,jl,je,
     &     jp,jp1,jp2,j1,j2	!auxiliars

c Variables to calculate real NPP 
         integer i6	
         real tleaf(3)!leaf turnover time (yr)
	     data tleaf /1.0,0.5,1.0/ !leaf turnover time for the 3 PFTs
         real sla !specific leaf area (m2/kg)
	     real cl2 !leaf compartment's carbon content (kgC/m2)
	     real ca2 !aboveground woody compartment's carbon content (kgC/m2)
	     real cf2 !fineroots compartment's carbon content (kgC/m2)
         real cb2 !belowground woody compartment's carbon content (kgC/m2)
         real csa  !sapwood compartment´s carbon content (5% of woody tissues) (kgC/m2)
         real rm !maintenance respiration (kgC/m2/yr)
	     real rml !leaf maintenance respiration (kgC/m2/yr)
         real rmf !fine roots maintenance respiration (kgC/m2/yr)
         real rms !sapwood maintenance respiration (kgC/m2/yr)
	     real rg !growth respiration (kgC/m2/yr)
	     real ncl !leaf N:C ratio (gN/gC)
	     real ncf !fine roots N:C ratio (gN/gC)
	     real ncs !sapwood N:C ratio(gN/gC)
	     real tsoil !soil temperature (ºC)
         real cl1
	     real ca1
         real cf1
	     real cb1
	     real rgl
         real rgf
         real rgs		 
	     real csai 
         real pt      !taxa potencial de fornecimento para transpiração (mm/dia)
         real csru !Specific root water uptake (0.5 mm/gC/dia; based in Pavlick et al (2013))
         real ad   !atmospheric demand for transpiration (mm/dia;based in Gerten et al. 2004)
         real emax !potential evapotranspiration (mm/dia)
         real alfm !maximum Priestley-Taylor coefficient (based in Gerten et al. 2004; used to calculate ad)
         real gm   !scaling conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)
         real gc   !Canopy conductance (mm/dia)(based in Gerten et al. 2004; used to calculate ad)
	     real beta_leaf,beta_awood,beta_bwood,beta_froot 
	     real ocp_tbe2,ocp_tbd2,ocp_herb2
       		 
		 
c
c
c
c=======================================================================
c Photosynthesis =======================================================
c
c Rubisco maximum carboxylaton rate
c (vm ; molCO2/m^2s) [Eq. 12]
c
      vm = (0.00004*2.0**(0.1*(temp-25.0)))/
     &           (1.0+exp(0.3*(temp-36.0)))
c
c Photo-respiration compensation point
c (mgama ; Pa) [Eq. 8]
c
      mgama = 21200.0/(5200.0*(0.57**(0.1*(temp-25.0))))
c
c Michaelis-Menton CO2 constant
c (f2 ; Pa) [Eq. 9]
c
      f2 = 30*(2.1**(0.1*(temp-25.0)))
c
c Michaelis-Menton O2 constant
c (f3 ; Pa) [Eq. 10]
c
      f3 = 30000.0*(1.2**(0.1*(temp-25.0)))
c
c Saturation partial pressure of water vapour at temperature 'temp' 
c (es ; Pa) [from WBM, to be used in Eq. 14]
c
      call tetens (temp,es)
c
c Saturated mixing ratio
c (rmax ; kg/kg) [Eq. 14]
c
      rmax = 0.622*(es/(p0-es))
c
c Moisture deficit at leaf level
c (r ; kg/kg) [Eq. 13]
c Considering that relative humidity is constant (0.685)
c
      r = (-0.315)*rmax
c
c Internal leaf CO2 partial pressure
c (ci ; Pa) [Eq. 11]
c
c      k16 = 0.875	!Generic value for forest types
c      k17 = 0.08	!Generic value for forest types
      k16 = 0.9	!Value for C3 grass types
      k17 = 0.1	!Value for C3 grass types
c
      ci = k16*(1-(r/k17))*(ca-mgama)+mgama
c
c Empirical function for reduction of NPP sensitivity to CO2
c based in simulations for pre-Industrial period (dimensionless)
      f10 = 1.23/(1+exp(-1*(ca-24.11)/8.93))
c
c Rubisco carboxilation limited photosynthesis rate
c (jc ; molCO2/m^2s) [Eq. 4]
c
      jc = vm*((ci-mgama)/(ci+(f2*f10*(1+(21200.0/f3)))))
c
c Light limited photosynthesis rate
c (jl ; molCO2/m^2s) [Eq. 5]
c
      jl = 0.08*(1.0-0.15)*ipar*((ci-mgama)/(ci+(2.0*mgama)))
c
c Transport limited photosynthesis rate
c (je ; molCO2/m^2s) [Eq. 6]
c
      je = 0.5*vm
c
c jp (minimum between jc and jl)
c [Eq. 3]
c
      a = 0.83
      b = (-1)*(jc+jl)
      c = jc*jl
      delta = (b**2)-4.0*a*c
c
      jp1=(-b-(sqrt(delta)))/(2.0*a)
      jp2=(-b+(sqrt(delta)))/(2.0*a)
	 jp= amin1(jp1,jp2)
c
c Leaf level gross photosynthesis (minimum between jc, jl and je)
c [f1, Eq. 2]
c
      a2 = 0.93
      b2 = (-1)*(jp+je)
      c2 = jp*je
      delta2 = (b2**2)-4.0*a2*c2
c
      j1=(-b2-(sqrt(delta2)))/(2.0*a2)
      j2=(-b2+(sqrt(delta2)))/(2.0*a2)
         f1a = amin1(j1,j2)
c
c Soil water
      wa = w/wmax
	  
c Soil water disponibility proportional to the PFT biomass occupation in a grid cell(ocp) 
	     if (i6.eq.1) then
		    wa=wa*ocp_tbe2
		    else if (i6.eq.2) then
		    wa=wa*ocp_tbd2
		    else if (i6.eq.3) then
		    wa=wa*ocp_herb2
			   endif
c
c Water stress response modifier (dimensionless)
c [f5 ; Eq. 21]
c 
      if (wa.gt.0.5) f5 = 1.0		!Not too lower in e.g. Amazonian dry season
      if ((wa.ge.0.205).and.(wa.le.0.5))
     &   f5 = (wa-0.205)/(0.5-0.205)
      if (wa.lt.0.205) f5 = wa		!Below wilting point f5 accompains wa (then Sahara is well represented)
	  
	    ! if (i2.eq.1) then  !first loop for potential NPP  
        !    f5 = f5	   
	    ! else if (i2.eq.2) then !second loop for real NPP 
!		       csru = 0.5 
!		     pt = csru*(cf2*1000.)*wa !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
!			   alfm = 1.391
!			   gm = 3.26*86400           !(*86400 transform mm/s to mm/dia)
!			   gc = (1/rc2)*86400000   !*86400000 transfor m/s to mm/dia)
!			   d =(emax*alfm)/(1+gm/gc)                 !(based in Gerten et al. 2004)
			   !d = emax
!			   f5 = 1-(exp(-1*(pt/d)))        !(Based in Pavlick et al. 2013)
!			    endif 
				
				 csru = 0.5 
		     pt = csru*(cf1*1000.)*wa !(based in Pavlick et al. 2013; *1000. converts kgC/m2 to gC/m2)
			   alfm = 1.391
			   gm = 3.26*86400           !(*86400 transform mm/s to mm/dia)
			   gc = (1/rc2)*86400000   !*86400000 transfor m/s to mm/dia)
			   d =(emax*alfm)/(1+gm/gc)                 !(based in Gerten et al. 2004)
			   !d = emax
			   f5 = 1-(exp(-1*(pt/d)))        !(Based in Pavlick et al. 2013)
			 
             
c
c ...f1
c
c Photosysthesis minimum and maximum temperature
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         f1 = f1a*f5
      else
         f1 = 0.0	!Temperature above/below photosynthesis windown 
      endif
	    
c
c Leaf area index (m2 leaf/m2 area)
c [lai ; Eq. 15]
c
c      laia  = 0.2*exp(2.5*(f1/0.000008))
      !laia  = 0.25*exp(2.5*(f1/0.000008)) !adjusted after using observed ipar
       !if (i2.eq.1) then
	!laia = laia
	!   print*, 'i1 laia',laia
	 !    else if (i2.eq.2) then
		!laia = (cl2*365*sla)
		   ! print*, 'i2 laia',laia
		
      !endif
	 sla = (0.030*1000.)*((365/(((tleaf(i6))*365)/12))**(-0.46)) !based on Pavlick et al (2013) - (*1000) transform to m2/kg and (/30.4) transform to months
  	 laia = (cl1*365*sla)
      
c 
c LAI with direct incident sun
c [sunlai; Eq. 16]      
      sunlai = (1.0-(exp(-0.5*laia)))/0.5
	 
c
c LAI below sunlai
c [shadelai, Eq. 17]
      shadelai = laia - sunlai
	  
c
c Scaling-up to canopy level (dimensionless)
c [f4 ; Eq. 18]
c
      f4 = (1.0-(exp(-0.5*laia)))/0.5	!sun 90° in the whole canopy, to be used for respiration
	    
c Sun/Shade approach to canopy scaling (based in de Pury & Farquhar 1997) 
c [f4sun, f4shade ; Eqs. 19, 20]
c
      f4sun = (1.0-(exp(-0.5*sunlai)))/0.5	!sun 90°
	     
      f4shade = (1.0-(exp(-1.5*shadelai)))/1.5	!sun ~20°
	      
c
c Canopy gross photosynthesis (kgC/m^2/yr)
c [ph ; Eq. 1]
c (0.012 converts molCO2 to kgC)
c [31557600 converts seconds to year (with 365.25 days)]
c
      ph = 0.012*31557600.0*f1*f4sun*f4shade
	 
	  
	   
c
c=======================================================================
c Plant (autotrophic) respiration ======================================
c 
      !if (i2.eq.1) then  !first loop for potential NPP      
c Leaf respiration (kgC/m2/yr)
c [rl ; Eq. 23]
c
      !rl = 0.012*31557600.0*0.015*vm*f4*f5
c
c Non-leaf parts respiration (kgC/m2/yr)
c [rp ; Eq. 24]
c
      !rp = 3.85*rl
c
c Autotrophic (plant) respiration (kgC/m2/yr)
c [ar ; Eq. 22]
      !ar = rl+rp
	  
	    
c
c	 	   
	   !  else if (i2.eq.2) then !second loop for real NPP 
c Maintenance respiration (kgC/m2/yr) (based in Ryan 1991)
!         csa = 0.05*(ca2+cb2) !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
!		    ncl = 0.034 !(gN/gC)
!		    ncf = 0.034 !(gN/gC)
!		    ncs = 0.003 !(gN/gC)
!		    rml = (ncl*cl2)*27*(exp(0.07*temp))
!		    rmf = (ncf*cf2)*27*(exp(0.07*tsoil))
!		    rms = (ncs*csa)*27*(exp(0.07*temp))
				
    		 csa = 0.05*(ca1+cb1) !sapwood carbon content (kgC/m2). 5% of woody tissues (Pavlick, 2013)
		    ncl = 0.034 !(gN/gC)
		    ncf = 0.034 !(gN/gC)
		    ncs = 0.003 !(gN/gC)
		    rml = (ncl*cl1)*27*(exp(0.07*temp))
		    rmf = (ncf*cf1)*27*(exp(0.07*tsoil))
		    rms = (ncs*csa)*27*(exp(0.07*temp))
			
		    rm = rml + rmf + rms
			
c Growth respiration (KgC/m2/yr)(based in Ryan 1991; Sitch et al. 2003; Levis et al. 2004)		 
!         csai = 0.05*((ca2-ca1)+(cb2-cb1))
!		    rgl = (0.25*((cl2-cl1)*365))
!		    rgf = (0.25*((cf2-cf1)*365))
!		    rgs = (0.25*(csai)*365)
			
			csai = 0.05*(beta_awood+beta_bwood)
		    rgl = (0.25*((beta_leaf)*365))
		    rgf = (0.25*((beta_froot)*365))
		    rgs = (0.25*(csai)*365)
		  
			rg = rgl + rgf + rgs

		    if (rg.lt.0) then
			   rg = 0
	        endif 		 
	
		 
c Autotrophic (plant) respiration (kgC/m2/yr)		 

		     ar = rm + rg

	   !  endif 

c
c Respiration minimum and maximum temperature
      if ((temp.ge.-10.0).and.(temp.le.50.0)) then
         ar = ar
      else
         ar = 0.0	!Temperature above/below respiration windown 
      endif
	  

c=======================================================================
c Productivity =========================================================
c
c Net primary productivity (kgC/m2/yr)
c [npp ; Eq. 25]
c

	     nppa = ph-ar

      
	     if (nppa.lt.0.0) nppa = 0.0 !maybe this is incorrect, but demands retuning of every biome limits
	   
         
c
      return
      end

c=======================================================================
c=======================================================================
c Microbial (heterotrophic) respiration ================================
c234567
      subroutine carbon2 (tsoil,f5,evap,laia, !input
     &                    cl,cs,hr)		 !output
c
c i/o variables
      real tsoil	!mean monthly soil temperature (oC)
      real f5		!stress response to soil moisture (dimensionless)
      real evap		!actual evapotranspiration (mm/day)
      real laia		!Leaf area index (m2 leaf/ m2 area)
c
      real cl		!Litter carbon (kgC/m2)
      real cs		!Soil carbon (kgC/m2)
      real hr		!Heterotrophic (microbial) respiration (kgC/m2/yr)
c
c internal variables
      real lf,f6,f7
c
c========================================================================
c Litter decayment function 
c (controlled by annual evapotranspiration)
c [f6 ; Eq. 29]
c
      f6 = 1.16*10**(-1.4553+0.0014175*(evap*365.0))
c
c Soil carbon storage function
c (controlled by temperature)
c [f7 ; Eq.31]
c
      f7 = 2.0**(0.1*(tsoil-25.0))
c
c Litterfall (kgC/m2)
c [lf ; Eq. 27]
c
      lf = 0.1*laia
c
c Litter carbon (kgC/m2)
c [cl ; Eq. 28]
c
      cl = lf/f6
c
c Soil carbon (kgC/m2)
c [cs ; Eq. 30]
c
      cs = ((0.3*cl)/(0.05*f7))*f5
c
c Heterotrophic (Microbial) respiration (kgC/m2/yr)
c [hr ; Eq. 26]
c
c Respiration minimum and maximum temperature
      if ((tsoil.ge.-10.0).and.(tsoil.le.50.0)) then
         hr = 0.25*(cl*(f6**2)+(cs*f5*evap*(f7**2))) !Litter and Soil respectively

      else
         hr = 0.0	!Temperature above/below respiration windown 
      endif
c
      return
      end
c
c=======================================================================
