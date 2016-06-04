program read_bin

implicit none
integer :: j
integer, parameter :: nx=720, ny=360
real*4 arr(nx,ny)    
    
    
    open (unit=21,file='pr.bin',status='unknown',&
          form='unformatted',access='direct',recl=nx*ny*4)
    open (unit=25,file='pr_jan.flt',status='unknown',&
          form='unformatted',access='direct',recl=nx*ny*4)
    read(21,rec = 1) arr
    close(21)
    write(25,rec = 1) arr
    close(25)
    !do j = 1, ny
    !    print*, arr(:,j)
    !end do
    
end program
     