c234567
c=====================================================================
c
c subroutine allocation calculates the daily carbon content of each
c compartment
c
c code written by Bianca Rius & David Lapola (27.Ago.2015)
c
c=====================================================================

      subroutine allocation (npp,cl1,ca1,cf1,cb1,cs1,cr1,co1,i6,   !input
     &                       cl2,ca2,cf2,cb2,cs2,cr2,co2)   !output
c
c
!     variables
	     integer i6   !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
         real npp     !potential npp (KgC/m2/yr)
	     real npp_aux !auxiliary variable to calculate potential npp in KgC/m2/day
         real cl1     !previous day carbon content on leaf compartment (KgC/m2)
         real cl2     !final carbon content on leaf compartment (KgC/m2)
         real ca1	  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
         real ca2     !final carbon content on aboveground woody biomass compartment (KgC/m2)
         real cf1	  !previous day carbon content on fine roots compartment (KgC/m2)
	     real cf2     !final carbon content on fine roots compartment (KgC/m2)
         real cb1     !previous day carbon content on belowground woody compartment (KgC/m2)
         real cb2     !final carbon content on belowground woody compartment (KgC/m2)
         real cs1	  !previous day carbon content on storage compartment (KgC/m2)
         real cs2     !final carbon content on storage compartment (KgC/m2)
         real cr1     !previous day carbon content on reproduction compartment (KgC/m2)
         real cr2     !final carbon content on reproduction compartment (KgC/m2)
         real co1     !previous day carbon content on other compartment (KgC/m2)
         real co2     !final carbon content on other compartment (KgC/m2)
         
		 
         real aleaf(3) 		 !npp percentage alocated to leaf compartment
	     data aleaf /0.23,0.23,0.32/
         real aawood (3) !npp percentage alocated to aboveground woody biomass compartment
         data aawood /0.22,0.22,0.001/
   	     real afroot(3)  !npp percentage alocated to fine roots compartment
	     data afroot /0.15,0.15,0.42/ 
         real abwood(3)  !npp percentage alocated to belowground woody biomass compartment
	     data abwood /0.10,0.10,0.001/
!         real asto(3)    !npp percentage alocated to storage compartment
!         data asto /0.10,0.10,0.10/
!         real arep(3)    !npp percentage alocated to reproduction compartment
!	     data arep /0.15,0.15,0.10/
!         real aother(3)  !npp percentage alocated to other compartment
!	     data aother /0.05,0.05,0.06/ 
         real tleaf(3)   !turnover time of the leaf compartment (yr)
         data tleaf /1.0,0.5,1.0/ 
         real tawood (3)  !turnover time of the aboveground woody biomass compartment (yr)
	     data tawood /20.0,20.0,20.0/
  	     real tfroot(3)  !turnover time of the fine roots compartment
	     data tfroot /1.0,1.0,1.0/
         real tbwood (3)		 !turnover time of the belowground woody biomass compartment
	     data tbwood /30.0,30.0,30.0/
!        real tsto  (3)    !turnover time of the storage compartmentturn
!	     data tsto /5.0,5.0,5.0/ 
!         real trep (3)   !turnover time of the reproduction compartment
!	     data trep /0.25,0.25,0.25/ 
!         real tother (3) !turnover time of the other compartment
!	     data tother /0.12,0.12,0.12/
         real sla     !specific leaf area (m2/kg)
		   
		   !!!!!!!!!!!ATENÇÃO AO SLA!!!!!!!!!!!!!
	     		 
c
c
c=====================================================================
c
c
! Carbon content of each compartment(KgC/m2)
c
c
c initialization
               
		         npp_aux = npp/365.0 	  !transform (KgC/m2/yr) in (KgC/m2/day)
	  

			 cl2 = (((aleaf(i6))*npp_aux) - (cl1/((tleaf(i6))*365)))+ cl1
			 ca2 = (((aawood(i6))*npp_aux) - (ca1/((tawood(i6))*365))) + ca1
             cf2 = (((afroot(i6))*npp_aux)-(cf1/((tfroot(i6))*365)))+cf1
			 cb2 = (((abwood(i6))*npp_aux)- (cb1/((tbwood(i6))*365))) + cb1
!			 cs2 = (((asto(i6))*npp_aux) - (cs1/((tsto(i6))*365))) + cs1
!			 cr2 = (((arep(i6))*npp_aux) - (cr1/((trep(i6))*365))) + cr1
!			 co2 = (((aother(i6))*npp_aux)- (co1/((tother(i6))*365))) + co1
			 

             return
      end
c
c
c=====================================================================
