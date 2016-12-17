      subroutine spinup (nppot,lsmk,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini)

      
      integer, parameter :: nx=192,ny=96,nt=100
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

c      real alc_leaf(mal,maf,maw),alc_froot(mal,maf,maw),
c     &	   alc_awood(mal,maf,maw)
      
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
            if (nint(lsmk(i,j)).ne.0) then
               do i6=1,npft

                  n = 0
                  
 10               continue

                  n = n + 1
                  k = mod(n,100)
                  if (k .eq. 0) k =100
                  
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
      end subroutine spinup
c     
      
      subroutine spinup2 (nppot,lsmk,
     &     cleafini_pft,cawoodini_pft,cfrootini_pft)

      integer i6
      parameter(nx=192,ny=96,nt=2500,npls=3,mal=80,maf=80,maw=80)
      real nppot(nx,ny),cleafini(nx,ny),lsmk(nx,ny),cawoodini(nx,ny),
     &     cfrootini(nx,ny),cbwoodini(nx,ny),cstoini(nx,ny),
     &     cotherini(nx,ny),crepini(nx,ny)
      real cleafini_tbe(nx,ny),cawoodini_tbe(nx,ny),
     &     cfrootini_tbe(nx,ny),cbwoodini_tbe(nx,ny),
     &     cstoini_tbe(nx,ny),cotherini_tbe(nx,ny),
     &     crepini_tbe(nx,ny),
     &     cleafini_tbd(nx,ny),cawoodini_tbd(nx,ny),
     &     cfrootini_tbd(nx,ny),cbwoodini_tbd(nx,ny),
     &     cstoini_tbd(nx,ny),cotherini_tbd(nx,ny),
     &     crepini_tbd(nx,ny),
     &     cleafini_herb(nx,ny),cawoodini_herb(nx,ny),
     &     cfrootini_herb(nx,ny),cbwoodini_herb(nx,ny),
     &     cstoini_herb(nx,ny),cotherini_herb(nx,ny),
     &     crepini_herb(nx,ny)
      real cleafini_pft(nx,ny,npls),cawoodini_pft(nx,ny,npls),
     &	   cfrootini_pft(nx,ny,npls)
      real cleafi_aux(nx,ny,nt),cawoodi_aux(nx,ny,nt),
     &     cfrooti_aux(nx,ny,nt),cbwoodi_aux(nx,ny,nt),
     &     cstoi_aux(nx,ny,nt),cotheri_aux(nx,ny,nt),crepi_aux(nx,ny,nt)
      real alc_leaf(mal,maf,maw),alc_froot(mal,maf,maw),
     &	   alc_awood(mal,maf,maw)
      
      real sensitivity,sensitivity2
      integer n
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
            if (lsmk(i,j).eq.1) then
               do i6=1,npls
                  do k=1,nt
                     n = n+1
		     if (k.eq.1) then
                        cleafi_aux(i,j,k) = aleaf(i6)*(nppot(i,j))
                        cawoodi_aux(i,j,k) = aawood(i6)*(nppot(i,j))
                        cfrooti_aux(i,j,k) = afroot(i6)*(nppot(i,j))
                        cbwoodi_aux(i,j,k) = abwood(i6)*(nppot(i,j))
                        cstoi_aux(i,j,k) = asto(i6)*(nppot(i,j))
                        cotheri_aux(i,j,k) = aother(i6)*(nppot(i,j)/365)
                        crepi_aux(i,j,k) = arep(i6)*(nppot(i,j)/365)
                     else
                        cleafi_aux(i,j,k) = ((aleaf(i6)*(nppot(i,j)))-
     &                       (cleafi_aux(i,j,k-1)/(tleaf(i6)))) +
     $                       cleafi_aux(i,j,k-1)
                        cawoodi_aux(i,j,k) = ((aawood(i6)*(nppot(i,j)))-
     &                       (cawoodi_aux(i,j,k-1)/(tawood(i6)))) +
     $                       cawoodi_aux(i,j,k-1)
                        cfrooti_aux(i,j,k) = ((afroot(i6)*(nppot(i,j)))-
     &                       (cfrooti_aux(i,j,k-1)/(tfroot(i6)))) +
     $                       cfrooti_aux(i,j,k-1)
                        cbwoodi_aux(i,j,k) = ((abwood(i6)*(nppot(i,j)))-
     &                       (cbwoodi_aux(i,j,k-1)/(tbwood(i6)))) +
     $                       cbwoodi_aux(i,j,k-1)
                        cstoi_aux(i,j,k) = ((asto(i6)*(nppot(i,j)))-
     &                       (cstoi_aux(i,j,k-1)/(tsto(i6)))) +
     $                       cstoi_aux(i,j,k-1)
                        cotheri_aux(i,j,k) = ((aother(i6)*(nppot(i,j)))-
     &                       (cotheri_aux(i,j,k-1)/(tother(i6)*365))) +
     $                       cotheri_aux(i,j,k-1)
                        crepi_aux(i,j,k) = ((arep(i6)*(nppot(i,j)))-
     &                       (crepi_aux(i,j,k-1)/(trep(i6)*365))) +
     $                       crepi_aux(i,j,k-1)
                        kk =  int(k*0.66)
                        if((cfrooti_aux(i,j,k)/cfrooti_aux(i,j
     $                       ,kk).lt.sensitivity).and.(cleafi_aux(i,j,k)
     $                       /cleafi_aux(i,j
     $                       ,kk).lt.sensitivity).and.(cawoodi_aux(i,j
     $                       ,k)/cawoodi_aux(i,j
     $                       ,kk).lt.sensitivity2).and.(cbwoodi_aux(i,j
     $                       ,k)/cbwoodi_aux(i,j
     $                       ,kk).lt.sensitivity2).and.(cstoi_aux(i,j,k)
     $                       /cstoi_aux(i,j
     $                       ,kk).lt.sensitivity).and.(cotheri_aux(i,j
     $                       ,k)/cotheri_aux(i,j
     $                       ,kk).lt.sensitivity).and.(crepi_aux(i,j,k)
     $                       /crepi_aux(i,j,kk).lt.sensitivity))   then
                           
                           cleafini_pft(i,j,i6)=cleafi_aux(i,j,k)
                           cawoodini_pft(i,j,i6)=cawoodi_aux(i,j,k)
                           cfrootini_pft(i,j,i6)=cfrooti_aux(i,j,k)
                           
                           
                           
                           exit
                        endif
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo	      
      
      return
      end subroutine spinup2			   
      
      
