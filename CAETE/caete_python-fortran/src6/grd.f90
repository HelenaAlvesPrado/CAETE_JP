module grdcell
  use global_pars
  implicit none
  private

!!$  interface
!!$     function photosynthesis_rate(vm,temp,p0,ipar,ll) result(f1a)
!!$       !returns instantaneous photosynthesis rate at leaf level (molCO2/m2/s)
!!$       
!!$       use global_pars
!!$       use photo_par
!!$       
!!$       implicit none
!!$       
!!$       real(kind=r4),intent(in) :: vm
!!$       real(kind=r4),intent(in) :: temp
!!$       real(kind=r4),intent(in) :: p0
!!$       real(kind=r4),intent(in) :: ipar
!!$       logical(kind=l1),intent(in) :: ll
!!$       
!!$       real(kind=r8) :: f2,f3            !Michaelis-Menten CO2/O2 constant (Pa)
!!$       real(kind=r8) :: mgama           !Photo-respiration compensation point (Pa)
!!$       real(kind=r8) :: rmax, r
!!$       real(kind=r8) :: ci
!!$       real(kind=r8) :: jp1
!!$       real(kind=r8) :: jp2
!!$       real(kind=r8) :: jp
!!$       real(kind=r8) :: jc
!!$       real(kind=r8) :: jl
!!$       real(kind=r8) :: je
!!$       real(kind=r8) :: b,c,c2,b2,es,j1,j2
!!$       real(kind=r8) :: delta, delta2,aux_ipar
!!$       real(kind=r8) :: f1a
!!$       real(kind=r4) :: tetens
!!$     end function photosynthesis_rate
!!$     ! external procedures declaration
!!$  end interface
  
contains
  type, public :: gridcell
     sequence
     ! state vars
     
     logical(kind=l1) :: in_data = .false.
     logical(kind=l1) :: completed = .false.
     ! these must be filled to run
     integer(kind=i4) :: xpos
     integer(kind=i4) :: ypos
     integer(kind=i4) :: count
     
     real(kind=r4),dimension(ntimes) :: temp
     real(kind=r4),dimension(ntimes) :: p0
     real(kind=r4),dimension(ntimes) :: ipar
     real(kind=r4) :: ca

     ! these are variables that store primary rates
     real(kind=r8),dimension(npls,ntimes) :: photo_rate ! (micromolCO2/m2/s) [leaf_level]

   contains
     function input_temp(this, xpos, ypos) result(this, temp)
       
       function input_pres(this, xpos, ypos) result()
         function input_ipar(this, xpos, ypos) result()
           
  end type gridcell
  
end module grdcell
