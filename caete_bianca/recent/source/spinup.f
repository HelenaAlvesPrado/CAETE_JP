!     NEW SUBROUTINE SPINUP
      
      subroutine spinup(nppot,
     &     cleafini,cfrootini,cawoodini,
     &     cbwoodini,cstoini,cotherini,crepini) 

      IMPLICIT NONE

      integer, parameter :: nt=5000
      integer, parameter :: npft=3
c     inputs
      integer i6, kk, k
      
      real :: nppot
      real :: sensitivity,sensitivity2

      
c     outputs
      real :: cleafini(npft)
      real :: cawoodini(npft)
      real :: cbwoodini(npft)
      real :: cfrootini(npft)
      real :: cstoini(npft)
      real :: cotherini(npft)
      real :: crepini(npft)

c     internal vars

      real cleafi_aux(nt)
      real cfrooti_aux(nt)
      real cawoodi_aux(nt)
      real cbwoodi_aux(nt)
      real cstoi_aux(nt)
      real cotheri_aux(nt)
      real crepi_aux(nt)

      
    
      real aleaf(3)             !npp percentage alocated to leaf compartment
      data aleaf /0.35,0.35,0.45/
      real aawood (3)           !npp percentage alocated to aboveground woody biomass compartment
      data aawood /0.40,0.40,0.001/
      real afroot(3)            !npp percentage alocated to fine roots compartment
      data afroot /0.25,0.25,0.55/ 
      real abwood(3)            !npp percentage alocated to belowground woody biomass compartment
      data abwood /0.10,0.10,0.001/
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

      do i6=1,npft
         do k=1,nt
            if (k.eq.1) then
               cleafi_aux(k) = aleaf(i6)*(nppot)
               cawoodi_aux(k) = aawood(i6)*(nppot)
               cfrooti_aux(k) = afroot(i6)*(nppot)
               cbwoodi_aux(k) = abwood(i6)*(nppot)
               cstoi_aux(k) = asto(i6)*(nppot)
               cotheri_aux(k) = aother(i6)*(nppot/365)
               crepi_aux(k) = arep(i6)*(nppot/365)
            else
               cleafi_aux(k) = ((aleaf(i6)*(nppot))-
     &              (cleafi_aux(k-1)/(tleaf(i6)))) + cleafi_aux(k-1)
               cawoodi_aux(k) = ((aawood(i6)*(nppot))-
     &              (cawoodi_aux(k-1)/(tawood(i6)))) + cawoodi_aux(k-1)
               cfrooti_aux(k) = ((afroot(i6)*(nppot))-
     &              (cfrooti_aux(k-1)/(tfroot(i6)))) + cfrooti_aux(k-1)
               cbwoodi_aux(k) = ((abwood(i6)*(nppot))-
     &              (cbwoodi_aux(k-1)/(tbwood(i6)))) + cbwoodi_aux(k-1)
               cstoi_aux(k) = ((asto(i6)*(nppot))-
     &              (cstoi_aux(k-1)/(tsto(i6)))) + cstoi_aux(k-1)
               cotheri_aux(k) = ((aother(i6)*(nppot))-
     &              (cotheri_aux(k-1)/(tother(i6)*365))) + cotheri_aux(k
     $              -1)
               crepi_aux(k) = ((arep(i6)*(nppot))-
     &              (crepi_aux(k-1)/(trep(i6)*365))) + crepi_aux(k-1)
               
               kk =  int(k*0.66)
               if((cfrooti_aux(k)/cfrooti_aux(kk).lt.sensitivity)
     $              .and.(cleafi_aux(k)/cleafi_aux(kk).lt.sensitivity)
     $              .and.(cawoodi_aux(k)/cawoodi_aux(kk).lt.
     $              sensitivity2).and.(cbwoodi_aux(k)/cbwoodi_aux(kk)
     $              .lt.sensitivity2).and.(cstoi_aux(k)
     $              /cstoi_aux(kk).lt.sensitivity).and.(cotheri_aux(k)
     $              /cotheri_aux(kk).lt.sensitivity).and.(crepi_aux(k)
     $              /crepi_aux(kk).lt.sensitivity))   then
                  
                  
                  cleafini(i6) = cleafi_aux(k) ! carbon content (kg m-2) 
                  cawoodini(i6) = cawoodi_aux(k)
                  cfrootini(i6) = cfrooti_aux(k)
                  cbwoodini(i6) = cbwoodi_aux(k)
                  cstoini(i6) = cstoi_aux(k)
                  cotherini(i6) = cotheri_aux(k)
                  crepini(i6) = crepi_aux(k)
                  exit
                  
               endif
            endif
         enddo
      enddo
      
      
      return
      end subroutine spinup				   
      
      
