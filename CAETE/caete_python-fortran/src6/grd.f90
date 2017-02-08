module grdcell
  use global_pars
  implicit none
  
contains
  type gridcell
     sequence
     ! state vars
     logical(kind=l1) :: in_data = .false.
     logical(kind=l1) :: completed = .false.
     ! these must be filled to run
     real(kind=r4),dimension(ntimes) :: temp
     real(kind=r4),dimension(ntimes) :: p0
     real(kind=r4),dimension(ntimes) :: w
     real(kind=r4),dimension(ntimes) :: rh
     real(kind=r4),dimension(ntimes) :: ipar
     real(kind=r4) :: ca

     ! these are variables that store primary rates
     real(kind=r8),dimension(npls,ntimes) :: photo_rate ! (micromolCO2/m2/s) [leaf_level]
  end type gridcell
  
end module grdcell
