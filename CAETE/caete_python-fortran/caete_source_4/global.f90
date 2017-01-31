module global_pars
  implicit none
  integer,parameter :: i4 = kind(0)
  integer,parameter :: r4 = kind(0.0)
  integer,parameter :: r8 = kind(0.0D0)
  integer,parameter :: rbig = selected_real_kind(16,300)
  integer(kind=i4),parameter :: npls = 12
  integer(kind=i4),parameter :: ntimes = 12
end module global_pars
