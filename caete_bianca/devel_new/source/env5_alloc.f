c234567
      program env
c     
c=======================================================================
c     CPTEC-PVM2
c     Adapted from env4.f. 
c     Code written by David Lapola
c     Last update: Aug/2007
c     
c     Compile with: g95 (or gfortran) env5.f wbm4.f carbon2.f -o a.exe (a.exe is the executable)
c     Execute with ./a.exe
c     Then run g95 (gfortran) pvm5.f -o pvm5.exe
c     Execute with ./pvm5.exe
c     
c=======================================================================
c     
c     Parameters and variables
c     ------------------------
c     
      implicit none
      
      integer, parameter :: nx=192, ny=96, npft=3
      real, parameter :: no_data = -9999.0
      integer i, j, k, p
      
      real lsmk(nx,ny),p0(nx,ny), npp_pot(nx,ny)
      real prec(nx,ny,12)
      real pr(nx,ny,12)
      real t(nx,ny,12)
      real temp(nx,ny,12)
      real par(nx,ny,12),ipar(nx,ny,12)
      real ca, npp_sca
      
C     INPUTS PARA WBM/ OUTPUTS DE SPINUP
      real  cleafini(nx,ny,npft)
      real cfrootini(nx,ny,npft)
      real cawoodini(nx,ny,npft)
      real cbwoodini(nx,ny,npft)
      real   cstoini(nx,ny,npft)
      real cotherini(nx,ny,npft)
      real   crepini(nx,ny,npft)
      real, dimension(npft) :: aux1, aux2, aux3, aux4
      real, dimension(npft) :: aux5, aux6, aux7
      
      
C     OUTPUTS PARA WBM
      real emaxm(nx,ny,12)
      real tsoil(nx,ny,12)
      
c     agora as variaveis para pfts
      real photo_pft(nx,ny,12,npft) !Monthly photosynthesis   (kgC/m2)
      real aresp_pft(nx,ny,12,npft) !Monthly autotrophic res  (kgC/m2)
      real   npp_pft(nx,ny,12,npft) !Monthly net primary produ (kgC/m2)
      
      real   lai_pft(nx,ny,12,npft) !Monthly leaf area index
      real  clit_pft(nx,ny,12,npft) !Monthly litter carbon
      real csoil_pft(nx,ny,12,npft) !Monthly soil carbon
      real hresp_pft(nx,ny,12,npft) !Monthly het resp          (kgC/m2)
      real   rcm_pft(nx,ny,12,npft) !stomatal (or canopy?) resistence s m-1
      
      
c     VARIAVEIS HIDROLOGICAS IMPORTANTES   
      real runom_pft(nx,ny,12,npft) !Runoff
      real evapm_pft(nx,ny,12,npft) !Actual evapotranspiration        
      real wsoil_pft(nx,ny,12,npft) !Soil moisture (mm)
      
      
c     carbon in plant tissues
      real  cleaf_pft(nx,ny,12,npft)
      real cawood_pft(nx,ny,12,npft)
      real cbwood_pft(nx,ny,12,npft)
      real cfroot_pft(nx,ny,12,npft)
      real   csto_pft(nx,ny,12,npft)
      real   crep_pft(nx,ny,12,npft)
      real cother_pft(nx,ny,12,npft)
      
C     maintenance and growth respiration    
      real monrml(nx,ny,12,npft)
      real monrmf(nx,ny,12,npft)
      real monrms(nx,ny,12,npft)
      real  monrm(nx,ny,12,npft)
      real monrgl(nx,ny,12,npft)
      real monrgf(nx,ny,12,npft)
      real monrgs(nx,ny,12,npft)
      real  monrg(nx,ny,12,npft)

c --------------------------------E N D--------------------------------
      

c Open files
c ----------


      open( 9,file='../inputs/psup.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(10,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(11,file='../inputs/prec.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(12,file='../inputs/temp.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(13,file='../inputs/ipar.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      open(20,file='../inputs/nppot.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
      
!     Read data
c     ---------
      
      read (9,rec=1) p0         !mb
      read (10,rec=1) lsmk      !land = 1; ocean = 0
      read (20,rec=1) npp_pot
      call read12 (11,pr)       !mm/month
      call read12 (12,t)        !oC
      call read12 (13,ipar)     !w/m2
      
c     Close files
c     ---------------
      
      close (9)
      close (10)
      close (11)
      close (12)
      close (13)
      close (20)
      
      
c     --- Atmospheric CO2 pressure (Pa) ! ppmv / Pa ------------------------
c     1 Pa CO2 = 9.901 ppmv CO2 (Adams et al. 2004)
!     ca= 18.18 !Pa (=180 ppmv; Last Glacial Maximum)
!     ca= 28.28 !Pa (=280 ppmv; Pre-Industrial Rev.)
!     ca= 35.35 !Pa (=350 ppmv; 1961-1990)
      ca= 350/9.901             !Pa (=350 ppmv; 1961-1990)
!     ca= 54.03 !Pa (=535 ppmv; SRES-B1 2080's)
!     ca= 73.73 !Pa (=730 ppmv; SRES-A2 2080's)
!     ca= ((73.73-35.35)*0.5)+35.35 !Pa half effect!!!
      
c     
c     
      do k=1,12
         do i=1,nx
            do j=1,ny
               
               
c     --- Photosynthetically active radiation reaching canopy --------------
c     (ipar ; Ein/m2/s) [Eq. 7] observed data from ISLSCP2
               par(i,j,k) = ipar(i,j,k)/(2.18e5) !converting to Ein/m2/s
               temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
               prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
               if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
               
            enddo
         enddo
      enddo
c     
c     
c     ------------------------------------------------------
c     Spinup to calculate initial carbon content for the compartments
c     -----------------------------------------------------
c     
      print*, 'spinup run ...'
      do i=1,nx
         do j=1,ny
            
            if (nint(lsmk(i,j)) .ne. 0) then
               
               npp_sca = npp_pot(i,j)
               
               do p = 1,npft
                  aux1(p) = 0.0
                  aux2(p) = 0.0
                  aux3(p) = 0.0
                  aux4(p) = 0.0
                  aux5(p) = 0.0
                  aux6(p) = 0.0
                  aux7(p) = 0.0
               enddo
               
               call spinup(npp_sca, aux1, aux2, aux3, aux4, aux5, aux6,
     $              aux7)
               
               do p=1,npft                 
                  cleafini(i,j,p)  = aux1(p)
                  cfrootini(i,j,p) = aux2(p)
                  cawoodini(i,j,p) = aux3(p)
                  cbwoodini(i,j,p) = aux4(p)
                  cstoini(i,j,p)   = aux5(p)
                  cotherini(i,j,p) = aux6(p)
                  crepini(i,j,p)   = aux7(p)
               enddo
            endif
         enddo
         if(mod(i,10) .eq. 0) print*, (real(i)/real(nx))*100.0, '%'
      enddo 
      print*, 'spinup end.'
c     -------------------------------------------------------
c     Calculate environmental variables (water balance model)
c     -------------------------------------------------------
      
      
      call wbm (prec,temp,lsmk,p0,ca,par,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini,
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft, ! out
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,
     &     monrml,monrmf,monrms,monrm,
     &     monrgl,monrgf,monrgs,monrg,
     &     cleaf_pft,cfroot_pft,cawood_pft,cbwood_pft, !outputs
     &     csto_pft, crep_pft, cother_pft)
      
c     
c     Program end
c     -----------
c     
      stop
      end program env
c     
c=======================================================================
      
      subroutine read12(nunit,var)
c     auxiliar reading routine
      parameter(nx=192,ny=96)
      integer nunit
      real var(nx,ny,12)
      real aux(nx,ny)
      do k=1,12
         read(nunit,rec=k) aux
         do i=1,nx
            do j=1,ny
               var(i,j,k) = aux(i,j)      
            enddo
         enddo
      enddo
      return
      end
      
