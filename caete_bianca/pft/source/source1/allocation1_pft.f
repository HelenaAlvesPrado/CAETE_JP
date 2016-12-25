c234567
c=====================================================================
c
c subroutine allocation calculates the daily carbon content of each
c compartment
c
c code written by Bianca Rius & David Lapola (27.Ago.2015)
c
c=====================================================================

      subroutine allocation (npp,cl1,ca1,cf1,i6,
     &                       alc_leaf2,alc_froot2,alc_awood2,
     &                       cl2,ca2,cf2)   !output
c
c
!     variables
	     
	     integer i6  !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
         real npp     !potential npp (KgC/m2/yr)
	     real npp_aux !auxiliary variable to calculate potential npp in KgC/m2/day
         real cl1     !previous day carbon content on leaf compartment (KgC/m2)
         real cl2     !final carbon content on leaf compartment (KgC/m2)
         real ca1	  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
         real ca2     !final carbon content on aboveground woody biomass compartment (KgC/m2)
         real cf1	  !previous day carbon content on fine roots compartment (KgC/m2)
	     real cf2     !final carbon content on fine roots compartment (KgC/m2)
        
         
		 
		        
		 
         real aleaf(3) 		 !npp percentage alocated to leaf compartment
	     data aleaf /0.35,0.35,0.45/
         real aawood (3) !npp percentage alocated to aboveground woody biomass compartment
         data aawood /0.40,0.40,0.001/
   	     real afroot(3)  !npp percentage alocated to fine roots compartment
	     data afroot /0.25,0.25,0.55/ 
         real tleaf(3)   !turnover time of the leaf compartment (yr)
         data tleaf /1.0,0.5,1.0/ 
         real tawood (3)  !turnover time of the aboveground woody biomass compartment (yr)
	     data tawood /30.0,30.0,30.0/
  	     real tfroot(3)  !turnover time of the fine roots compartment
	     data tfroot /1.0,1.0,1.0/

	     
           
c
c=====================================================================
c
c
! Carbon content of each compartment(KgC/m2)
c
c
c initialization
               
		         npp_aux = npp/365.0 	  !transform (KgC/m2/yr) in (KgC/m2/day)
	  

		          cl2 = ((aleaf(i6)*npp_aux) - 
     &				  (cl1/((tleaf(i6))*365)))+ cl1
	              cf2 = ((afroot(i6)*npp_aux) - 
     &				  (cf1/((tfroot(i6))*365)))+ cf1
	              ca2 = ((aawood(i6)*npp_aux) - 
     &				  (ca1/((tawood(i6))*365)))+ ca1
	             
				
           
             return
      end
c
c
c=====================================================================
