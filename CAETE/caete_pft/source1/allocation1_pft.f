c234567
c=====================================================================
c
c subroutine allocation calculates the daily carbon content of each
c compartment
c
c code written by Bianca Rius & David Lapola (27.Ago.2015)
c
c=====================================================================
      
      subroutine allocation (i6,npp,cl1,ca1,cf1,
     &     cl2,ca2,cf2)         !output
c     
c     
!     variables
      
      integer i6                !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
      real npp                  !potential npp (KgC/m2/yr)
c      real npp_aux              !auxiliary variable to calculate potential npp in KgC/m2/day
      real cl1                  !previous day carbon content on leaf compartment (KgC/m2)
      real cl2                  !final carbon content on leaf compartment (KgC/m2)
      real ca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real ca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
      real cf1                  !previous day carbon content on fine roots compartment (KgC/m2)
      real cf2                  !final carbon content on fine roots compartment (KgC/m2)
      double precision ca2_64, cl2_64, cf2_64, npp_aux
      double precision aleaf64
      double precision afroot64
      double precision aawood64 
      double precision tleaf64
      double precision tfroot64
      double precision tawood64
      
      
      
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.40, 0.25, 0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.30, 0.40, 0.0/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.30, 0.25, 0.55/ 
      real tleaf(3)             !turnover time of the leaf compartment (yr)
      data tleaf /2.0, 0.7, 1.0/ 
      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
      data tawood /30.0, 3.0, 0.0/
      real tfroot(3)            !turnover time of the fine roots compartment
      data tfroot /3.0, 2.0, 1.0/
      
      
      
c     
c=====================================================================
c     
c     
!     Carbon content of each compartment(KgC/m2)
c     
c     
c     initialization

         
      npp_aux = real(npp,8)/365.0 !transform (KgC/m2/yr) in (KgC/m2/day)
      call critical_value2(npp_aux)
      
      
      aleaf64 = real(aleaf(i6),8)
      afroot64 = real(aleaf(i6),8)
      aawood64 = real(aawood(i6),8)
      
      tleaf64 =  real(tleaf(i6),8)
      tfroot64 = real(tfroot(i6),8)
      tawood64 = real(tawood(i6),8)
      
      cl2_64 = ((aleaf64*npp_aux) - (real(cl1,8)/((tleaf64)*365.)))+
     &    real(cl1,8)
      cf2_64 = ((afroot64*npp_aux) - (real(cl1,8)/((tfroot64)*365.)))
     &    +real(cl1,8)
      if(aawood(i6) .gt. 0.0) then
         ca2_64 = ((aawood64*npp_aux) - (real(cl1,8)/((tawood64)
     &       *365.)))+ real(cl1,8)
      else
         ca2_64 = 0.0
      endif
      
      cl2 = real(cl2_64,4)
      call critical_value(cl2)
      if(cl2 .lt. 0.0) cl2 = 0.0
      cf2 = real(cf2_64,4)
      call critical_value(cf2)
      if(cf2 .lt. 0.0) cf2 = 0.0
      ca2 = real(ca2_64,4)
      call critical_value(ca2)
      if(ca2 .lt. 0.0) ca2 = 0.0
      
      return
      end
c     
c     
c=====================================================================
