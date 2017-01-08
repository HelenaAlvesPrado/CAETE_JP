c=====================================================================
c     
c     subroutine allocation calculates the daily carbon content of each
c     compartment
c     
c     code written by Bianca Rius & David Lapola (27.Ago.2015)
c     
c=====================================================================
      
      subroutine allocation (pft, npp ,scl1,sca1,scf1,
     &    scl2,sca2,scf2)          !output
c     
c     
!     variables
      
      integer pft   
      real npp                  !potential npp (KgC/m2/yr)
      real npp_aux              !auxiliary variable to calculate potential npp in KgC/m2/day
      real scl1                  !previous day carbon content on leaf compartment (KgC/m2)
      real scl2                  !final carbon content on leaf compartment (KgC/m2)
      real sca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
      real sca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
      real scf1                  !previous day carbon content on fine roots compartment (KgC/m2)
      real scf2                  !final carbon content on fine roots compartment (KgC/m2)
c      real cb1                  !previous day carbon content on belowground woody compartment (KgC/m2)
c      real cb2                  !final carbon content on belowground woody compartment (KgC/m2)
c      real cs1                  !previous day carbon content on storage compartment (KgC/m2)
c      real cs2                  !final carbon content on storage compartment (KgC/m2)
c     real cr1                  !previous day carbon content on
c     reproduction compartment (KgC/m2)
c     real cr2                  !final carbon content on reproduction
c     compartment (KgC/m2)
c     real co1                  !previous day carbon content on other
c     compartment (KgC/m2)
c     real co2                  !final carbon content on other
c     compartment (KgC/m2)
      
      
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.40, 0.25, 0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.30, 0.40, 0.0/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.30, 0.25, 0.55/ 
c      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
c      data abwood /0.10,0.10,0.001/
c      real asto(3)              !npp percentage alocated to storage compartment
c      data asto /0.10,0.10,0.10/
c      real arep(3)              !npp percentage alocated to reproduction compartment
c      data arep /0.15,0.15,0.10/
c      real aother(3)            !npp percentage alocated to other compartment
c      data aother /0.05,0.05,0.06/ 
      real tleaf(3)             !turnover time of the leaf compartment (yr)
      data tleaf /2.0, 0.7, 1.0/ 
      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
      data tawood /30.0, 3.0, 0.0/
      real tfroot(3)            !turnover time of the fine roots compartment
      data tfroot /3.0, 2.0, 1.0/
c      real tbwood (3)           !turnover time of the belowground woody biomass compartment
c      data tbwood /40.0,40.0,40.0/
c      real tsto  (3)            !turnover time of the storage compartmentturn
c      data tsto /5.0,5.0,5.0/ 
c      real trep (3)             !turnover time of the reproduction compartment
c      data trep /0.25,0.25,0.25/ 
c      real tother (3)           !turnover time of the other compartment
c      data tother /0.12,0.12,0.12/
!       real sla                  !specific leaf area (m2/kg)
      
!     !!!!!ATENÇÃO AO SLA!!!!!!!!!!!!!
      
c     ================================================================
c     
!     Carbon content of each compartment(KgC/m2)
c     
c     
c     initialization
      if((scl1 .lt. 0.0000001) .or. (scf1 .lt. 0.0000001)) then
         IF(NPP .lt. 0.0000001) THEN
            scl2 = 0.0
            scf2 = 0.0
            sca2 = 0.0 
            goto 10
         ENDIF
      endif   
      npp_aux = npp/365.0       !transform (KgC/m2/yr) in (KgC/m2/day)
c      call critical_value(npp_aux)
      scl2 = scl1 + (aleaf(pft) * npp_aux) -(scl1 /(tleaf(pft)*365.0))
         
      scf2 = scf1 +(afroot(pft) * npp_aux)-(scf1 /(tfroot(pft)*365.0))
      if(aawood(pft) .gt. 0.0) then
         sca2 = sca1 +(aawood(pft)*npp_aux)-(sca1/(tawood(pft)*365.0))
      else
         sca2 = 0.0
      endif

      
c      call critical_value(scl2)
c      call critical_value(scf2)
c      call critical_value(sca2)

      if(scl2 .lt. 0.0) scl2 = 0.0
      if(scf2 .lt. 0.0) scf2 = 0.0
      if(sca2 .lt. 0.0) sca2 = 0.0
      
C     cb2 = (((abwood(pft))*npp_aux)- (cb1/((tbwood(pft))*365))) + cb1
C      cs2 = (((asto(pft))*npp_aux) - (cs1/((tsto(pft))*365))) + cs1
C      cr2 = (((arep(pft))*npp_aux) - (cr1/((trep(pft))*365))) + cr1
C      co2 = (((aother(pft))*npp_aux)- (co1/((tother(pft))*365))) + co1
      
c      if(cl2 .gt. 0) print*, cl2, cf2, ca2, 'carbon final'
 10   continue
      return
      end
c     
      
c=====================================================================
c
c subroutine allocation calculates the daily carbon content of each
c compartment
c
c code written by Bianca Rius & David Lapola (27.Ago.2015)
c
c=====================================================================
c      
c      subroutine allocation (i6,npp,cl1,ca1,cf1,
c     &     cl2,ca2,cf2)         !output
cc     
cc     
c!     variables
c      
c      integer i6                !index for calculate the 3 PFTs (i6.eq.1 - TBE; i6.eq.2 - TBD; i6.eq.3 - HERB)
c      real npp                  !potential npp (KgC/m2/yr)
cc      real npp_aux              !auxiliary variable to calculate potential npp in KgC/m2/day
c      real cl1                  !previous day carbon content on leaf compartment (KgC/m2)
c      real cl2                  !final carbon content on leaf compartment (KgC/m2)
c      real ca1                  !previous day carbon content on aboveground woody biomass compartment(KgC/m2)
c      real ca2                  !final carbon content on aboveground woody biomass compartment (KgC/m2)
c      real cf1                  !previous day carbon content on fine roots compartment (KgC/m2)
c      real cf2                  !final carbon content on fine roots compartment (KgC/m2)
c      double precision ca2_64, cl2_64, cf2_64, npp_aux
c      double precision aleaf64
c      double precision afroot64
c      double precision aawood64 
c      double precision tleaf64
c      double precision tfroot64
c      double precision tawood64
c      
c      
c      
c      real aleaf(3)             !npp percentage alocated to leaf compartment
c      data aleaf /0.40, 0.25, 0.45/
c      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
c      data aawood /0.30, 0.40, 0.0/
c      real afroot(3)            !npp percentage alocated to fine roots compartment
c      data afroot /0.30, 0.25, 0.55/ 
c      real tleaf(3)             !turnover time of the leaf compartment (yr)
c      data tleaf /2.0, 0.7, 1.0/ 
c      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
c      data tawood /30.0, 3.0, 0.0/
c      real tfroot(3)            !turnover time of the fine roots compartment
c      data tfroot /3.0, 2.0, 1.0/
c      
c      
c      
cc     
cc=====================================================================
cc     
cc     
c!     Carbon content of each compartment(KgC/m2)
cc     
cc     
cc     initialization
c
c         
c      npp_aux = real(npp,8)/365.0 !transform (KgC/m2/yr) in (KgC/m2/day)
c      call critical_value2(npp_aux)
c      
c      
c      aleaf64 = real(aleaf(i6),8)
c      afroot64 = real(aleaf(i6),8)
c      aawood64 = real(aawood(i6),8)
c      
c      tleaf64 =  real(tleaf(i6),8)
c      tfroot64 = real(tfroot(i6),8)
c      tawood64 = real(tawood(i6),8)
c      
c      cl2_64 = ((aleaf64*npp_aux) - (real(cl1,8)/((tleaf64)*365.)))+
c     &    real(cl1,8)
c      cf2_64 = ((afroot64*npp_aux) - (real(cl1,8)/((tfroot64)*365.)))
c     &    +real(cl1,8)
c      if(aawood(i6) .gt. 0.0) then
c         ca2_64 = ((aawood64*npp_aux) - (real(cl1,8)/((tawood64)
c     &       *365.)))+ real(cl1,8)
c      else
c         ca2_64 = 0.0
c      endif
c      
c      cl2 = real(cl2_64,4)
c      call critical_value(cl2)
c      if(cl2 .lt. 0.0) cl2 = 0.0
c      cf2 = real(cf2_64,4)
c      call critical_value(cf2)
c      if(cf2 .lt. 0.0) cf2 = 0.0
c      ca2 = real(ca2_64,4)
c      call critical_value(ca2)
c      if(ca2 .lt. 0.0) ca2 = 0.0
c      
c      return
c      end
cc     
cc     
cc=====================================================================
