module grdcell
  use global_pars
  implicit none
  private

!       month,w1,g1,s1,ts,temp,prec,p0,ae,ca,ipar,rh&
!     &,cl1_pft,ca1_pft,cf1_pft,w2,g2,s2,smavg,ruavg,evavg,epavg&
!     &,phavg,aravg,nppavg,laiavg,clavg,csavg,hravg,rcavg,rmavg,rgavg&
!     &,cleafavg_pft,cawoodavg_pft,cfrootavg_pft,ocpavg
  
contains
  type, public :: gridcell
     sequence

     
     logical(kind=l1) :: in_data = .false.
     logical(kind=l1) :: completed = .false.
     integer(kind=i4) :: xpos
     integer(kind=i4) :: ypos
     integer(kind=i4) :: count

     
     ! model inputs
     real(kind=r4),dimension(ntimes) :: temp
     real(kind=r4),dimension(ntimes) :: p0
     real(kind=r4),dimension(ntimes) :: ipar
     real(kind=r4),dimension(ntimes) :: rh
     real(kind=r4),dimension(ntimes) :: prec
     real(kind=r4),dimension(ntimes) :: ipar
     real(kind=r4) :: ca

     ! 
     real(kind=r8),dimension(npls,ntimes) :: photo_rate ! (micromolCO2/m2/s) [leaf_level]

  end type gridcell
  
end module grdcell



module model
  use global_pars
  implicit none
  private

!  
contains

  subroutine wbm()
    !---
  end subroutine wbm


  
  subroutine budget
    !=---
  end subroutine budget


  
end module model
