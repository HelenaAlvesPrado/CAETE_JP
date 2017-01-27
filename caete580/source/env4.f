c234567
!     
      program env
!     
      implicit none
!     
!     =======================================================================
!     CAETE-DV
!     Adapted from env4.f.
!     Code written by David Lapola and Helena Alves do Prado e Bianca Rius
c     &  jpdarela  

!     Last update: jan/2017
!     
!     Compile with: g95 (or gfortran) -g -Wall -mcmodel=medium env4.f wbm4.f prod4.f 
!     Execute with ./a.out
!     =======================================================================



!     Parameters and variables
!     ------------------------
!     
      integer i,j,k,p
      integer, parameter :: nx=720,ny=360,q=580
      real,parameter :: no_data = -9999.0
!     
!     Inputs
!     ------
!     
      real ca                   !CO2 concentration (Pa)
      real lsmk(nx,ny)          !Land=1/Ocean=0
      real p0(nx,ny,12)         !Atmospheric pressure (mb)
      real ps(nx,ny,12)         ! auxiliar to read atm pressure information
      real prec(nx,ny,12)       !Precipitation (mm/month)
      real pr(nx,ny,12)         !Auxiliar_precipitation (mm/month)
      real temp(nx,ny,12)       !Temperature (oC)
      real t(nx,ny,12)          !Auxiliar_temperature (oC)
      real par(nx,ny,12)        !Incident photosynthetic active radiation (Ein/m2/s)
      real ipar(nx,ny,12)       !Auxiliar_incident photosynthetic active radiation (w/m2)
      real rhs(nx,ny,12)        !Relative humidity
      real rhaux(nx,ny,12)      !RHS auxiliar
      real, dimension(nx,ny,q) :: cleafin = no_data
     $                          ,cawoodin = no_data
     $                          ,cfrootin = no_data


!     
!     Outputs
!     -------

c     Vou declarar aqui os outputas para a wbm, note que estas varia-
c     eis recebem os mesmos nomes das variaveis declaradas na definicao da wbm
c     porem, elas nao sao as mesmas variveis-- todas as variaveis da wbm sao
c     secretas para o env5.f.
      
      
c     -------------------------O U T P U T S--P A R A--W B M ------------
      
c     primeiro algumas variaveis que nao sao dependentes de pfts

      real emaxm(nx,ny,12)
      real tsoil(nx,ny,12)
      
      

c      real photo_pft(nx,ny,12) !Monthly photosynthesis   (kgC/m2)
c      real aresp_pft(nx,ny,12) !Monthly autotrophic res  (kgC/m2)
c      real npp_pft(nx,ny,12)  !Monthly net primary produ (kgC/m2)
c      
c      real lai_pft(nx,ny,12)  !Monthly leaf area index
c      real clit_pft(nx,ny,12) !Monthly litter carbon
c      real csoil_pft(nx,ny,12) !Monthly soil carbon
c      real hresp_pft(nx,ny,12) !Monthly het resp          (kgC/m2)
c      real rcm_pft(nx,ny,12) 
c      
cc     VARIAVEIS HIDROLOGICAS IMPORTANTES   
c      real runom_pft(nx,ny,12) !Runoff
c      real evapm_pft(nx,ny,12) !Actual evapotranspiration        
c      real wsoil_pft(nx,ny,12) !Soil moisture (mm)
c      
c     variables related to carbon allocation and autothrophic respiration (bianca)
      real, dimension(nx,ny) :: cleaf_pft,cawood_pft,cfroot_pft


C      WARNING - NEW VARIABLES ---
      
      real gridcell_ocp(nx,ny,q) !  final grid cell occupation for each pft (percentage of area)
      
c     variaveis do spinup
      real, dimension(q) :: aux1, aux2, aux3
      real npp_pot(nx,ny,12)
      real, dimension(nx,ny) :: aux_npp
      real npp_sca  

c     -----FIM DA DEFINICAO DE VARIAVEIS PARA RODAR O MODELO--
C
c     variaveis anuais -
      real, dimension(nx,ny) :: ave_ph = 0.0
      real, dimension(nx,ny) :: ave_ar = 0.0
      real, dimension(nx,ny) :: ave_npp = 0.0
      real, dimension(nx,ny) :: ave_lai = 0.0
      real, dimension(nx,ny) :: ave_clit = 0.0
      real, dimension(nx,ny) :: ave_cs = 0.0
      real, dimension(nx,ny) :: ave_hr = 0.0
      real, dimension(nx,ny) :: ave_rc = 0.0
      real, dimension(nx,ny) :: ave_runom = 0.0
      real, dimension(nx,ny) :: ave_evap = 0.0
      real, dimension(nx,ny) :: ave_wsoil = 0.0    
      
C     THESE WILL RECEIVE MEANS BETWEEN q PFTs and for each pft (ex. ph to mean; ph1 to pft 1)
      real, dimension(nx,ny,12) :: ph
      real, dimension(nx,ny,12) :: ar
      real, dimension(nx,ny,12) :: npp
      real, dimension(nx,ny,12) :: lai
      real, dimension(nx,ny,12) :: clit
      real, dimension(nx,ny,12) :: csoil
      real, dimension(nx,ny,12) :: hr
      real, dimension(nx,ny,12) :: rcm
      real, dimension(nx,ny,12) :: evaptr
      real, dimension(nx,ny,12) :: wsoil
      real, dimension(nx,ny,12) :: runom
      
C     NEW OUTPUTS (AUTOTRF RESPIRATION, ALLOCATION)
      

      real, dimension(nx,ny,12) :: rm
      real, dimension(nx,ny,12) :: rg 
      
C     -------END DECLARATION----------------------------------------
      
!     Open INPUT files
!     ==========
      open( 9,file='../inputs/lsmk.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(10,file='../inputs/ps.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(11,file='../inputs/pr.bin',status='old',
     &    form='unformatted',access='direct',recl=4*nx*ny)
      
      open(12,file='../inputs/tas.bin',status='old',
     &    form='unformatted',access='direct',recl=4*nx*ny)
      
      open(13,file='../inputs/rsds.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(14,file='../inputs/hurs.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

      open(26,file='../inputs/npp.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)

c      open(27,file='../spinup/clini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(28,file='../spinup/cfini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
c
c      open(29,file='../spinup/cwini.bin',status='old',
c     &    form='unformatted',access='direct',recl=4*nx*ny)
!     Read data
!     =========
      
      read (9,rec=1) lsmk
c     read (26,rec=1) aux_npp 
      
      call readx(10,ps,12)
      call readx(11,pr,12)
      call readx(12,t,12)
      call readx(13,ipar,12)
      call readx(14,rhaux,12)
      call readx(26,npp_pot,12)
c     call readx(29,cleafin,q)
c     call readx(28,cfrootin,q)
c     call readx(29,cawoodin,q)
      
!     Close files
!     ===========
!     
       close( 9)
       close(10)
       close(11)
       close(12)
       close(13)
       close(14)
       close(26)
c       close(27)
c       close(28)
c       close(29)
       

c     Calculating annual npp
       do i =1,nx
          do j=1,ny
             if(nint(lsmk(i,j)) .ne. 0) then 
                aux_npp(i,j) = 0.0
                do k = 1,12
                   aux_npp(i,j) = aux_npp(i,j) + (npp_pot(i,j,k)/12.) 
                enddo
             else
                aux_npp(i,j) = no_data
             endif
          enddo
       enddo
       
!     calling spinup
       print*, 'running spinup'
       do i=1,nx
          do j=1,ny
             if (nint(lsmk(i,j)) .ne. 0) then
                npp_sca = aux_npp(i,j)
                
                do p=1,q   
                   aux1(p) = 0.0
                   aux2(p) = 0.0
                   aux3(p) = 0.0
                   gridcell_ocp(i,j,p) = 0.0
                enddo
                
                call spinup(npp_sca, aux1, aux2, aux3)
                
                do p=1,q   
                  cleafin(i,j,p)  = aux1(p)
                  cfrootin(i,j,p) = aux2(p)
                  cawoodin(i,j,p) = aux3(p)
               enddo
            else
               do p = 1,q
                  cleafin(i,j,p)  = no_data
                  cfrootin(i,j,p) = no_data
                  cawoodin(i,j,p) = no_data
               enddo
            endif
         enddo
!     if(mod(nx,10) .eq. 0)print*, (real(i)/real(nx))*100.0, '%'
      enddo

!     estes arquivos ficaram muito grandes ... pensar em uma estrategia melhor
      
c      call nan2ndt(cleafin, q) 
c      open(10,file='../spinup/clini.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call savex(10, cleafin, q)
c      
c      call nan2ndt(cfrootin, q)
c      open(10,file='../spinup/cfini.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call savex(10, cfrootin, q)
c      
c      call nan2ndt(cawoodin, q)
c      open(10,file='../spinup/cwini.bin',
c     &    status='unknown',form='unformatted',
c     &    access='direct',recl=4*nx*ny)
c      call savex(10, cawoodin, q)
c!     ===========
      
      
      
      do i=1,nx
         do j=1,ny            
!     this block set ocean grid cells to no_data 
            if(nint(lsmk(i,j)) .eq. 0) then
               cleaf_pft(i,j)  = no_data
               cfroot_pft(i,j) = no_data
               cawood_pft(i,j) = no_data
               do p = 1,q
                  gridcell_ocp(i,j,p) = no_data
               enddo
            endif
!     ------------------------------------------
            
            do k=1,12
               rhs(i,j,k) =  rhaux(i,j,k) / 100. !(rhaux(60,41,10) / 100.0) !Humidade relativa de manaus em outubro
               par(i,j,k) = ipar(i,j,k)/2.18E5 !Converting to Ein/m2/s
               temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
               p0(i,j,k) = ps(i,j,k) * 0.01 ! transforamando de pascal pra mbar (kPa)
               prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
c     if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0  
            enddo
         enddo
      enddo
      
!     Atmospheric CO2 pressure (Pa) !Ppmv / Pa
      ca= 363/9.901             !Pa (=363 ppmv; 1981-2010)
      
      
!     =======================================
!     Calculate environmental variables (wbm)
!     =======================================
      
      call wbm (prec,temp,lsmk,p0,ca,par,rhs,cleafin,cawoodin,
     &    cfrootin, emaxm, tsoil, ph,ar,npp,lai,clit,csoil,hr,
     &    rcm,runom,evaptr,wsoil,rm,rg,cleaf_pft,cawood_pft,
     &    cfroot_pft,gridcell_ocp)   
      
      
!     SAVE RESULTS TO FILES ! arquivos gigantes > 700Mbytes
      open(10,file='../outputs/gridcell_ocp.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      call savex(10, gridcell_ocp, q)


      open(10,file='../outputs/cleaf.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      write(10,rec=1) cleaf_pft
      close(10)

      open(10,file='../outputs/cawood.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      write(10, rec=1) cawood_pft
      close(10)

      open(10,file='../outputs/cfroot.bin',
     &    status='unknown',form='unformatted',
     &    access='direct',recl=4*nx*ny)
      write(10, rec=1) cfroot_pft
      close(10)

      
      open(10,file='../outputs/ph.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, ph)
      
      open(10,file='../outputs/ar.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, ar)

      open(10,file='../outputs/npp.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, npp)

      open(10,file='../outputs/clit.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, clit)

      open(10,file='../outputs/csoil.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, csoil)

      open(10,file='../outputs/hr.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, hr)

      open(10,file='../outputs/rcm.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rcm)

      open(10,file='../outputs/runom.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, runom)

      open(10,file='../outputs/evaptr.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, evaptr)

      open(10,file='../outputs/wsoil.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, wsoil)

      open(10,file='../outputs/rm.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rm)

      open(10,file='../outputs/rg.bin',
     &     status='unknown',form='unformatted',
     &     access='direct',recl=4*nx*ny)
      call save_file12(10, rg)

      do i = 1,nx
         do j = 1,ny
            if(nint(lsmk(i,j)) .ne. 0) then
               do k = 1,12
                  ave_ph(i,j) = ave_ph(i,j) + ph(i,j,k)/12.
                  ave_ar(i,j) = ave_ar(i,j) + ar(i,j,k)/12.
                  ave_npp(i,j) = ave_npp(i,j) + npp(i,j,k)/12.
                  ave_lai(i,j) = ave_lai(i,j) + lai(i,j,k)/12.
                  ave_clit(i,j) = ave_clit(i,j) + clit(i,j,k)/12.
                  ave_cs(i,j) = ave_cs(i,j) + csoil(i,j,k)/12.
                  ave_hr(i,j) = ave_hr(i,j) + hr(i,j,k)/12.
                  ave_rc(i,j) = ave_rc(i,j) + rcm(i,j,k)/12.
                  ave_runom(i,j) = ave_runom(i,j) + runom(i,j,k)/12.
                  ave_evap(i,j) = ave_evap(i,j) + evaptr(i,j,k)/12.
                  ave_wsoil(i,j) = ave_wsoil(i,j) + wsoil(i,j,k)/12.
               enddo
            else
               ave_ph(i,j) = no_data
               ave_ar(i,j) = no_data
               ave_npp(i,j) = no_data
               ave_lai(i,j) = no_data
               ave_clit(i,j) = no_data
               ave_cs(i,j) = no_data
               ave_hr(i,j) = no_data
               ave_rc(i,j) = no_data
               ave_runom(i,j) = no_data
               ave_evap(i,j) = no_data
               ave_wsoil(i,j) = no_data
            endif
         enddo
      enddo

      open(50,file='../outputs/ambientais.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(50,rec=1) ave_npp
      write(50,rec=2) ave_rc
      write(50,rec=3) ave_ar
      write(50,rec=4) ave_lai
      write(50,rec=5) ave_clit
      write(50,rec=6) ave_cs
      write(50,rec=7) ave_hr
      write(50,rec=8) ave_runom
      write(50,rec=9) ave_evap
      write(50,rec=10) ave_wsoil
      write(50,rec=11) ave_ph
      close(50)

      
      stop
      end program env

!     ==================================================
      subroutine pft_par(par, dt) 
      implicit none
      integer, parameter :: vars = 580 
      integer :: par            ! parameter number 
      real, dimension(vars) :: dt
      
!     dt1 = aleaf
!     dt2 = aawood
!     dt3 = afroot
!     dt4 = tleaf
!     dt5 = tawood
!     dt6 = tfroot
!     dt7 = g1
!     dt8 = p21
!     DT9 = JMAX

      open(233,file='../inputs/pls.bin',status='old',
     &    form='unformatted',access='direct',recl=4*580)

      if(par .gt. 0 .and. par .lt. 10) then
         read(233,rec=par) dt
      else
         print*, 'search failed'
      endif
      close(233)
      return
      end subroutine pft_par
      
c     ==================================================
      subroutine spinup(nppot,
     &     cleafini,cfrootini,cawoodini)
c     &     cbwoodini,cstoini,cotherini,crepini) 
      IMPLICIT NONE

      integer, parameter :: nt=7000
      integer, parameter :: npfts=580
      
c     inputs
      integer i6, kk, k
      
      real :: nppot
      real :: sensitivity

c     outputs
      real :: cleafini(npfts)
      real :: cawoodini(npfts)
      real :: cfrootini(npfts)
      real :: ocp(npfts)
      logical :: inutil(npfts)
      real*8 cleafi_aux(nt)
      real*8 cfrooti_aux(nt)
      real*8 cawoodi_aux(nt)

    
      real aleaf(npfts)             !npp percentage alocated to leaf compartment
      real aawood (npfts)           !npp percentage alocated to aboveground woody biomass compartment
      real afroot(npfts)            !npp percentage alocated to fine roots compartmentc 
      real tleaf(npfts)             !turnover time of the leaf compartment (yr)
      real tawood (npfts)           !turnover time of the aboveground woody biomass compartment (yr)
      real tfroot(npfts)            !turnover time of the fine roots compartment

!     dt1 = aleaf
!     dt2 = aawood
!     dt3 = afroot
!     dt4 = tleaf
!     dt5 = tawood
!     dt6 = tfroot
!     dt7 = g1
!     dt8 = p21
!     DT9 = JMAX

      call pft_par(1, aleaf)
      call pft_par(2, aawood)
      call pft_par(3, afroot)
      call pft_par(4, tleaf)
      call pft_par(5, tawood)
      call pft_par(6, tfroot)

 
      sensitivity = 1.10
      if(nppot .le. 0.0) goto 200
      do i6=1,npfts
         do k=1,nt
            if (k.eq.1) then
               cleafi_aux (k) =  aleaf(i6)*(nppot)
               cawoodi_aux(k) = aawood(i6)*(nppot)
               cfrooti_aux(k) = afroot(i6)*(nppot)

            else
               if(aawood(i6) .gt. 0.0) then
                  cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &                (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
                  cawoodi_aux(k) = ((aawood(i6)*(nppot))-
     &                (cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k
     &                -1)
                  cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &                (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k
     &                -1)
               else
                  cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &                (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
                  cawoodi_aux(k) = 0.0
                  cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &                (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k
     &                -1)
               endif
                  

               
               kk =  nint(k*0.66)
               if(cawoodi_aux(kk) .gt. 0.0) then
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $                .and.(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)
     $                .and.(cawoodi_aux(k)/cawoodi_aux(kk).lt.
     $                sensitivity)) then
                     
                     cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                     cfrootini(i6) = real(cfrooti_aux(k),4)
                     cawoodini(i6) = real(cawoodi_aux(k),4)
                     exit
                  ENDIF
               else
                  if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $                .and.(cleafi_aux(k)
     &                /cleafi_aux(kk).lt.sensitivity)) then
                     
                     cleafini(i6) = real(cleafi_aux(k),4) ! carbon content (kg m-2)
                     cfrootini(i6) = real(cfrooti_aux(k),4)
                     cawoodini(i6) = 0.0
                     exit
                  endif   
               endif
            endif
         enddo                  !nt
      enddo                     ! npfts 
      call pft_area_frac(cleafini, cfrootini, cawoodini, ocp, inutil)
      do i6 = 1,npfts
         cleafini(i6) = cleafini(i6) * ocp(i6)
         cleafini(i6) = cleafini(i6) * ocp(i6)
         cleafini(i6) = cleafini(i6) * ocp(i6)
      enddo
 200  continue
      return
      end subroutine spinup
!     ================================

      subroutine readx(nunit,var,x)
!     auxiliar reading routine
      integer,parameter :: nx=720,ny=360
      integer nunit,x
      real var(nx,ny,x)
      real aux(nx,ny)
      do k=1,x
         read(nunit,rec=k) aux
         do i=1,nx
            do j=1,ny
               var(i,j,k) = aux(i,j)
            enddo
         enddo
      enddo
      return
      end
!     ================================

      subroutine save_file12(nunit, var)
      parameter(nx = 720, ny = 360, nt = 12)
      integer i, j, k
      real var(nx,ny,nt)
      real waux(nx,ny)
      do k=1,nt
         do i=1,nx
            do j=1,ny
               waux(i,j) = var(i,j,k)
            enddo
         enddo
         write(nunit,rec=k) waux
      enddo
      close(nunit)
      return
      end subroutine save_file12
!     ================================

      subroutine savex(nunit, var, x)
      integer, parameter :: nx=720,ny=360
      integer i, j, k, x
      real var(nx,ny,x)
      real waux(nx,ny)
      do k=1,x
         do i=1,nx
            do j=1,ny
               waux(i,j) = var(i,j,k)
            enddo
         enddo
         write(nunit,rec=k) waux
      enddo
      close(nunit)
      return
      end subroutine savex

      
!     ===============================
      subroutine nan2ndt(var, nl)

      integer nl
      integer, parameter :: nx=720, ny=360
      integer i,j
      real, parameter :: no_data = -9999.0
      real var(nx,ny,nl)
      
      do i = 1,nx
         do j = 1,ny
            do k = 1,nl
               
               if(var(i,j,k) .gt. 1E12)then
                  var(i,j,k) = 0.0
c                  print*, 'var .gt. 1e32',i,j,k
               endif
               
               if(var(i,j,k) .lt. 1E-12) then
                  var(i,j,k) = 0.0
c                  print*, 'var .lt. 1e-32',i,j,k
               endif
               
               if(isnan(var(i,j,k))) then
                  var(i,j,k) = no_data
c                  print*, 'NaN found', '---place',i,j,k
               endif
               
            enddo
         enddo
      enddo  
      
      return
      end subroutine nan2ndt
