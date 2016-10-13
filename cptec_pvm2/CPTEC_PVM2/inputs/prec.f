c234567
      program dry_prec
c
      integer nx, ny
      parameter(nx=192,ny=96)
      real prec(nx,ny,12),pr(nx,ny,12),anpr(nx,ny,12)
      real ds(nx,ny,16),avep(nx,ny,16)
      real aux1(nx,ny),aux2(nx,ny)
      character*3 gcm1(2)
      data gcm1 /'INM','MRI'/
      character*4 gcm2(5)
      data gcm2 /'BCCR','CNRM','GISS','IPSL','NCAR'/
      character*5 gcm3(3)
      data gcm3 /'CGCM3','CSIRO','MIROC'/
      character*6 gcm4(4)
      data gcm4 /'ECHAM5','GFDL20','GFDL21','HADCM3'/
      integer ndmonth(12) !number of days for each month
      data ndmonth /31,28,31,30,31,30,31,31,30,31,30,31/
      character*2 scen
c
      do kk=1,16
         do i=1,nx
         do j=1,ny
      ds(i,j,kk)=0.0
      avep(i,j,kk)=0.0
         enddo
         enddo
      enddo
c
      scen = 'A2' !'B1'
c 
      open(11,file='prec.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (11,pr)  !mm/month
      close (11)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      if(pr(i,j,k).lt.100.0) ds(i,j,1)=ds(i,j,1)+1
      avep(i,j,1)=avep(i,j,1)+(pr(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
      kk=1
      do ii=1,2
      kk=kk+1
      write(*,*) 'GCM:',gcm1(ii),', scenario SRES',scen
      open(12,file='../anomalias/SRES'//scen//'/prec/'//
     & gcm1(ii)//'_'//scen//'_2070_2099_anom_pr.bin',status='old',
     & form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (12,anpr)  !mm/mo
      close (12)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      prec(i,j,k) = (pr(i,j,k)+anpr(i,j,k)) !+pr(i,j,k)*0.2
      if(prec(i,j,k).lt.100.0) ds(i,j,kk)=ds(i,j,kk)+1
      avep(i,j,kk)=avep(i,j,kk)+(prec(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
      enddo !GCM loop
c
      do ii=1,5
      kk=kk+1
      write(*,*) 'GCM:',gcm2(ii),', scenario SRES',scen
      open(12,file='../anomalias/SRES'//scen//'/prec/'//
     & gcm2(ii)//'_'//scen//'_2070_2099_anom_pr.bin',status='old',
     & form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (12,anpr)  !mm/mo
      close (12)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      prec(i,j,k) = (pr(i,j,k)+anpr(i,j,k)) !+pr(i,j,k)*0.2
      if(prec(i,j,k).lt.100.0) ds(i,j,kk)=ds(i,j,kk)+1
      avep(i,j,kk)=avep(i,j,kk)+(prec(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
      enddo !GCM loop
c
      do ii=1,3
      kk=kk+1
      write(*,*) 'GCM:',gcm3(ii),', scenario SRES',scen
      open(12,file='../anomalias/SRES'//scen//'/prec/'//
     & gcm3(ii)//'_'//scen//'_2070_2099_anom_pr.bin',status='old',
     & form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (12,anpr)  !mm/mo
      close (12)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      prec(i,j,k) = (pr(i,j,k)+anpr(i,j,k)) !+pr(i,j,k)*0.2
      if(prec(i,j,k).lt.100.0) ds(i,j,kk)=ds(i,j,kk)+1
      avep(i,j,kk)=avep(i,j,kk)+(prec(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
      enddo !GCM loop
c
      do ii=1,4
      kk=kk+1
      write(*,*) 'GCM:',gcm4(ii),', scenario SRES',scen
      open(12,file='../anomalias/SRES'//scen//'/prec/'//
     & gcm4(ii)//'_'//scen//'_2070_2099_anom_pr.bin',status='old',
     & form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (12,anpr)  !mm/mo
      close (12)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      prec(i,j,k) = (pr(i,j,k)+anpr(i,j,k)) !+pr(i,j,k)*0.2
      if(prec(i,j,k).lt.100.0) ds(i,j,kk)=ds(i,j,kk)+1
      avep(i,j,kk)=avep(i,j,kk)+(prec(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
      enddo !GCM loop
c
c HadCM3 with increased precipitation
      open(12,file='../anomalias/SRES'//scen//'/prec/'//
     & 'HADCM3_'//scen//'_2070_2099_anom_pr.bin',status='old',
     & form='unformatted',access='direct',recl=4*nx*ny)
      call read12 (12,anpr)  !mm/mo
      close (12)
c
      do k=1,12
         do i=1,nx
         do j=1,ny
      prec(i,j,k) = (pr(i,j,k)+anpr(i,j,k))+pr(i,j,k)*0.2
      if(prec(i,j,k).lt.100.0) ds(i,j,16)=ds(i,j,16)+1
      avep(i,j,16)=avep(i,j,16)+(prec(i,j,k)/12/ndmonth(k))
         enddo
         enddo
      enddo
c
c open output files
      open(20,file='dry_season_'//scen//'.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      open(21,file='ave_prec'//scen//'.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      do k=1,16
         do i=1,nx
         do j=1,ny
            aux1(i,j) = ds(i,j,k)
            aux2(i,j) = avep(i,j,k)
         enddo
         enddo
         write(20,rec=k) aux1
         write(21,rec=k) aux2
      enddo
      close(20)
      close(21)
c
      stop
      end
c
c=======================================================================

      subroutine read12 (nunit,var)
c auxiliar reading routine
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

