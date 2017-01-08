c     23456
      
      program test
      
      integer, parameter ::  p = 1
      
      real, dimension(7) :: attr
      
      call pft_par(p, attr)
      print*, attr
      
      stop
      
      end program test
      
      
      subroutine pft_par(par, dt)
      
!     input
      integer, parameter :: vars = 7 
      integer :: par            ! parameter number 
      real, dimension(vars) :: dt,dt1,dt2,dt3,dt4,dt5,dt6,dt7,dt8
      
      
!     1 = g1
!     2 = p21 
!     3 = aleaf
!     4 = aawood
!     5 = afroot
!     6 = tleaf
!     7 = tawood
!     8 = tfroot
      
!     PFTS
      
!     1 = tropical evergreen tree
!     2 = tropical deciduous-forest-tree
!     3 = tropical deciduous-savana-tree
!     4 = tropical herb
!     5 = tropical grass
!     6 = temperate tree
!     7 = temperate herb
      
!     
      if(par .eq. 1 ) then      ! g1
!     PFT         1       2       3       4       5       6       7 
         data dt1 /3.04,   2.67,   2.0,    3.0,    1.4,    2.05,   1.95/
         dt(:) = dt1(:)
      else if(par .eq. 2) then  ! p21
         data dt2 /3.2e-5, 3.1e-5, 3.0e-5, 3.3e-5, 2.8e-5, 5.0e-5, 4.0e
     &       -5/
         dt(:) = dt2(:)
      else if(par .eq. 3) then
         data dt3 /0.30,   0.38,   0.36,   0.45,   0.55,   0.35,   0.55/
         dt(:) = dt3(:)
      else if(par .eq. 4) then  ! awood
         data dt4 /0.28,   0.30,   0.27,   0.10,   0.0,    0.40,   0.12/
         dt(:) = dt4(:)
      else if(par .eq. 5) then  ! afroot
         data dt5 /0.42,   0.32,   0.37,   0.55,   0.45,   0.25,   0.23/
         dt(:) = dt5(:)
      else if(par .eq. 6) then  ! tleaf
         data dt6 /7.0,    2.0,    1.0,    2.0,    1.0,    1.0,    2.0 /
         dt(:) = dt6(:)
      else if(par .eq. 7) then  ! tawood
         data dt7 /35.0,   30.0,   25.0,   2.0,    0.0,    32.0,   2.5/
         dt(:) = dt7(:)
      else if(par .eq. 8) then  ! tfroot
         data dt8 /4.0,    3.5,    3.8,    1.5,    1.0,    3.8,    2.2/ 
         dt(:) = dt8(:)
      else
         print*, "your search failed"
      endif
      
      return
      end subroutine pft_par
      
      
