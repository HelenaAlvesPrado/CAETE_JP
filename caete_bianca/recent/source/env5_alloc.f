c234567
      program env
c
c=======================================================================
c CPTEC-PVM2
c Adapted from env4.f. 
c Code written by David Lapola
c Last update: Aug/2007
c
c Compile with: g95 (or gfortran) env5.f wbm4.f carbon2.f -o a.exe (a.exe is the executable)
c Execute with ./a.exe
c Then run g95 (gfortran) pvm5.f -o pvm5.exe
c Execute with ./pvm5.exe
c
c=======================================================================
c
c Parameters and variables
c ------------------------
c
c     implicit none
      
      integer cell_id(6003), var_id, ni
      integer, parameter :: nx=192
      integer, parameter :: ny=96
      integer, parameter :: npls=3, npft=3 ! nao precisa de 2 
      real waux(nx,ny)
      integer i6                ! loop index 3 PFTs (i6.eq.1-TBE; i6.eq.2-TBD; i6.eq.3-HERB)

c     spinup
      integer supindex 
      
!     inputs for wbm!
      
      real lsmk(nx,ny)
      real p0(nx,ny)
      real pr(nx,ny,12)
      real t(nx,ny,12)
      real prec(nx,ny,12)
      real temp(nx,ny,12)
      real par(nx,ny,12)
      real ipar(nx,ny,12)
      real ca

!     wbm outputs!
      
      real tsoil(nx,ny,12)
      real emaxm(nx,ny,12)

!     estas serao por pft (nx,ny,12,npft)
      real, dimension(nx,ny,12,npft) :: wsoil_pft = 0.0
      real, dimension(nx,ny,12,npft) :: evapm_pft = 0.0
      real, dimension(nx,ny,12,npft) :: rcm_pft = 0.0
      real, dimension(nx,ny,12,npft) :: runom_pft = 0.0


c     carbon cycle
      real, dimension(nx,ny,12,npft) ::  photo_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  aresp_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  npp_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  lai_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  clit_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  csoil_pft = 0.0
      real, dimension(nx,ny,12,npft) ::  hresp_pft = 0.0

c     maintenance resp PRECISO SABER O QUE SO ESTAS VARIAVEIS 
      real, dimension(nx,ny,12,npft) ::  monrml = 0.0
      real, dimension(nx,ny,12,npft) ::  monrmf = 0.0
      real, dimension(nx,ny,12,npft) ::  monrms = 0.0 
      real, dimension(nx,ny,12,npft) ::  monrm = 0.0
      real, dimension(nx,ny,12,npft) ::  monrgl = 0.0
      real, dimension(nx,ny,12,npft) ::  monrgf = 0.0
      real, dimension(nx,ny,12,npft) ::  monrgs = 0.0
      real, dimension(nx,ny,12,npft) ::  monrg = 0.0
      
c     carbon allocation in plant tissues
      real, dimension(nx,ny,npft) :: cleafini = 0.0
      real, dimension(nx,ny,npft) :: cfrootini = 0.0
      real, dimension(nx,ny,npft) :: cawoodini = 0.0
      real, dimension(nx,ny,npft) :: cbwoodini = 0.0
      real, dimension(nx,ny,npft) :: cstoini = 0.0
      real, dimension(nx,ny,npft) :: cotherini = 0.0
      real, dimension(nx,ny,npft) :: crepini = 0.0
       
c     carbon pools in vegetation
      real, dimension(nx,ny,12,npft) :: cleaf = 0.0
      real, dimension(nx,ny,12,npft) :: cawood  = 0.0
      real, dimension(nx,ny,12,npft) :: cfroot = 0.0
      real, dimension(nx,ny,12,npft) :: cbwood = 0.0
      real, dimension(nx,ny,12,npft) :: csto = 0.0
      real, dimension(nx,ny,12,npft) :: crep = 0.0
      real, dimension(nx,ny,12,npft) :: cother = 0.0


c     Open files
c     ----------


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
      open(20,file='../outputs/nppot.bin',status='old',
     &     form='unformatted',access='direct',recl=4*nx*ny)
      
! Read data
c ---------

      read (9,rec=1) p0     !mb
      read (10,rec=1) lsmk      !land = 1; ocean = 0
      read (20,rec=1) nppot
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


      
      do i=1,nx
         do j=1,ny
            do k=1,12
               
               par(i,j,k) = ipar(i,j,k)/(2.18e5) !converting to Ein/m2/s
               temp(i,j,k) = t(i,j,k) !+ant(i,j,k) !uncomment to use future anomalies
               prec(i,j,k) = pr(i,j,k) !+anpr(i,j,k) !+pr(i,j,k)*0.2 !uncomment to use future anomalies
               if (prec(i,j,k).lt.0.0) prec (i,j,k) = 0.0
               
            enddo
         enddo
      enddo
      
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
c ------------------------------------------------------
c  Spinup to calculate initial carbon content for the compartments
c -----------------------------------------------------
c
      
      call spinup (nppot,lsmk,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini)	 
      
      
c     -------------------------------------------------------
c     Calculate environmental variables (water balance model)
c     -------------------------------------------------------
c     
      call wbm (prec,temp,lsmk,p0,ca,par,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini,
     &     emaxm, tsoil, photo_pft,aresp_pft,npp_pft,lai_pft,
     &     clit_pft,csoil_pft, hresp_pft,rcm_pft,runom_pft,
     &     evapm_pft,wsoil_pft,
     &     monrml,monrmf,monrms,monrm,
     &     monrgl,monrgf,monrgs,monrg,
     &     cleaf,cawood,cfroot,cbwood,csto,crep,cother)

      
c Program end
c -----------
c
      stop
      end
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
	  
c========================================================================

      subroutine spinup (nppot,lsmk,i6,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini)

      
      integer, parameter :: nx=192,ny=96,nt=12
      integer, parameter :: npft=3,mal=80,maf=80,maw=80

c     inputs

      integer i6, kk
      real nppot(nx,ny),lsmk(nx,ny)

c     outputs
      real cawoodini(nx,ny,npft),cbwoodini(nx,ny,npft)
      real cfrootini(nx,ny,npft),cstoini(nx,ny,npft)
      real cotherini(nx,ny,npft),crepini(nx,ny,npft)
      real cleafini(nx,ny, npft)

c     internal vars

      real cleafi_aux(nx,ny,nt),cawoodi_aux(nx,ny,nt),
     &     cfrooti_aux(nx,ny,nt),cbwoodi_aux(nx,ny,nt),
     &     cstoi_aux(nx,ny,nt),cotheri_aux(nx,ny,nt),crepi_aux(nx,ny,nt)

C    AUX VARS PARA GUARDAR OS VALORES DO MES ANTERIOR

      real cleafi_aux1
      real cawoodi_aux1
      real cbwoodi_aux1
      real cfrooti_aux1
      real cstoi_aux1
      real cotheri_aux1
      real crepi_aux1

      real alc_leaf(mal,maf,maw),alc_froot(mal,maf,maw),
     &	   alc_awood(mal,maf,maw)
      
      real sensitivity,sensitivity2
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.35,0.35,0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.40,0.40,0.000000001/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.25,0.25,0.55/ 
      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
      data abwood /0.10,0.10,0.000000001/
      real asto(3)              !npp percentage alocated to storage compartment
      data asto /0.10,0.10,0.10/
      real arep(3)              !npp percentage alocated to reproduction compartment
      data arep /0.15,0.15,0.10/
      real aother(3)            !npp percentage alocated to other compartment
      data aother /0.05,0.05,0.06/ 
      real tleaf(3)             !turnover time of the leaf compartment (yr)
      data tleaf /1.0,0.5,1.0/ 
      real tawood (3)           !turnover time of the aboveground woody biomass compartment (yr)
      data tawood /30.0,30.0,30.0/
      real tfroot(3)            !turnover time of the fine roots compartment
      data tfroot /1.0,1.0,1.0/
      real tbwood (3)           !turnover time of the belowground woody biomass compartment
      data tbwood /40.0,40.0,40.0/
      real tsto  (3)            !turnover time of the storage compartmentturn
      data tsto /5.0,5.0,5.0/ 
      real trep (3)             !turnover time of the reproduction compartment
      data trep /0.25,0.25,0.25/ 
      real tother (3)           !turnover time of the other compartment
      data tother /0.12,0.12,0.12/
      sensitivity = 1.10
      sensitivity2 = 1.40
      do i=1,nx
         do j=1,ny
            if (nint(lsmk(i,j)).eq.1) then
               do i6=1,npft

                  n = 0
                  
 10               continue

                  n = n + 1
                  k = mod(n,12)
                  if (k .eq. 0) k =12
                  
                  if (k.eq.1) then
                     cleafi_aux(i,j,k) = aleaf(i6)*(nppot(i,j))
                     cawoodi_aux(i,j,k) = aawood(i6)*(nppot(i,j))
                     cfrooti_aux(i,j,k) = afroot(i6)*(nppot(i,j))
                     cbwoodi_aux(i,j,k) = abwood(i6)*(nppot(i,j))
                     cstoi_aux(i,j,k) = asto(i6)*(nppot(i,j))
                     cotheri_aux(i,j,k) = aother(i6)*(nppot(i,j)/365)
                     crepi_aux(i,j,k) = arep(i6)*(nppot(i,j)/365)

                     
                  else
                     ! guardando os valores de k-1 
                     cleafi_aux1 = cleafi_aux(i,j,k-1)
                     cawoodi_aux1  = cawoodi_aux(i,j,k-1)
                     cbwoodi_aux1  = cbwoodi_aux(i,j,k-1)
                     cfrooti_aux1  = cfrooti_aux(i,j,k-1)
                     cstoi_aux1  = cstoi_aux(i,j,k-1)
                     cotheri_aux1 = cotheri_aux(i,j,k-1)
                     crepi_aux1  = crepi_aux(i,j,k-1)

                     
                     cleafi_aux(i,j,k) = ((aleaf(i6)*(nppot(i,j)))
     $                    -(cleafi_aux1/(tleaf(i6)))) + cleafi_aux1

                     cawoodi_aux(i,j,k) = ((aawood(i6)*(nppot(i,j)))
     $                    -(cawoodi_aux1/(tawood(i6))))+ cawoodi_aux1
                     
                     cfrooti_aux(i,j,k) = ((afroot(i6)*(nppot(i,j)))
     $                    -(cfrooti_aux1/(tfroot(i6)))) + cfrooti_aux1
                     
                     cbwoodi_aux(i,j,k) = ((abwood(i6)*(nppot(i,j)))
     $                    -(cbwoodi_aux1/(tbwood(i6)))) + cbwoodi_aux1

                     cstoi_aux(i,j,k) = ((asto(i6)*(nppot(i,j)))
     $                    -(cstoi_aux1/(tsto(i6)))) + cstoi_aux1
                     
                     cotheri_aux(i,j,k) = ((aother(i6)*(nppot(i,j)))-
     $                    (cotheri_aux1/(tother(i6)*365)))+cotheri_aux1
                     
                     crepi_aux(i,j,k) = ((arep(i6)*(nppot(i,j)))-
     $                    (crepi_aux1/(trep(i6)*365))) + crepi_aux1
                     
                     kk =  int(k*0.66)
                     
                     if((cfrooti_aux(i,j,k)/cfrooti_aux(i,j
     $                    ,kk).lt.sensitivity).and.(cleafi_aux(i,j,k)
     $                    /cleafi_aux(i,j
     $                    ,kk).lt.sensitivity).and.(cawoodi_aux(i,j
     $                    ,k)/cawoodi_aux(i,j
     $                    ,kk).lt.sensitivity2).and.(cbwoodi_aux(i,j
     $                    ,k)/cbwoodi_aux(i,j
     $                    ,kk).lt.sensitivity2).and.(cstoi_aux(i,j,k)
     $                    /cstoi_aux(i,j
     $                    ,kk).lt.sensitivity).and.(cotheri_aux(i,j
     $                    ,k)/cotheri_aux(i,j
     $                    ,kk).lt.sensitivity).and.(crepi_aux(i,j,k)
     $                    /crepi_aux(i,j,kk).lt.sensitivity))   then
                        
                        cbwoodini(i,j,i6) = cbwoodi_aux(i,j,k)
                        cstoini(i,j,i6) = cstoi_aux(i,j,k)
                        cotherini(i,j,i6) = cotheri_aux(i,j,k)
                        crepini(i,j,i6) = crepi_aux(i,j,k)
                        cleafini(i,j,i6) = cleafi_aux(i,j,k)
                        cawoodini(i,j,i6) = cawoodi_aux(i,j,k)
                        cfrootini(i,j,i6) = cfrooti_aux(i,j,k)
                        goto 100
                     else
                        goto 10
                     endif
                  endif
 100              continue
               enddo
            endif
         enddo
      enddo	
      
      return
      end				   
      
      
	 
	 
	 
	 
