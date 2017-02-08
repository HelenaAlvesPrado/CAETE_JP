module global_pars
  implicit none
  integer,parameter :: i4 = kind(0)
  integer,parameter :: r4 = kind(0.0)
  integer,parameter :: r8 = kind(0.0D0)
  integer,parameter :: rbig = selected_real_kind(16,300)
  integer(kind=i4),parameter :: npls = 12
  integer(kind=i4),parameter :: ntimes = 12
  integer(kind=i4),parameter :: nx = 720
  integer(kind=i4),parameter :: ny = 360
  
contains
  
  subroutine readx(nunit,var,x)
    
    integer(kind=i4), intent(in) ::  nunit,x
    real(kind=r4),intent(inout) :: var(nx,ny,x)
    real(kind=r4) :: aux(nx,ny)
    integer(kind=i4) :: i,j,k
    do k=1,x
       read(nunit,rec=k) aux
       do i=1,nx
          do j=1,ny
             var(i,j,k) = aux(i,j)
          enddo
       enddo
    enddo
    return
  end subroutine readx
  !     ================================
  
  subroutine save_file12(nunit, var)
    integer(kind=i4), intent(in) ::  nunit
    real(kind=r4),intent(in) :: var(nx,ny,ntimes)
    real(kind=r4) :: waux(nx,ny)
    integer(kind=i4) :: i,j,k
    do k=1,ntimes
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
  
  subroutine savex(nunit, var, x)
    integer(kind=i4), intent(in) ::  nunit,x
    real(kind=r4),intent(in) :: var(nx,ny,x)
    real(kind=r4) :: waux(nx,ny)
    integer(kind=i4) :: i,j,k
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
  
end module global_pars
