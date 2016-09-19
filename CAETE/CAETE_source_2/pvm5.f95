c234567
      program pvm5
c
c=====================================================================
c
c Adapted from pvm4.f
c Code written by David Lapola & Marcos Oyama
c Last updates
c MDO: 24 Aug 2007.
c DML: 27 Aug 2007.
c
c=====================================================================
c
c i/o variables
c d2wm : month of the transition between wet and dry season
c
      parameter(nx=192,ny=96)
      real*4 tmin(nx,ny),npp(nx,ny),snpp(nx,ny),hresp(nx,ny),
     &       csoil(nx,ny),photo(nx,ny),aresp(nx,ny),ave_wsoil(nx,ny),
     &       ave_evap(nx,ny),ave_rc(nx,ny)
      real*4 lsmk(nx,ny),bpot(nx,ny)
      real*4 wsoil2(nx,ny,12),waux2(12)
      real*4 ipar(nx,ny,12),avgpar(nx,ny)
      real ratio
      real*4 d2wm(nx,ny)
      real*4 u(nx,ny,12),anu(nx,ny,12),uwind(nx,ny,12),ud2w(nx,ny)
      real*4 auxbpot(nx,ny),auxnpp(nx,ny),auxsnpp(nx,ny),
     & auxhresp(nx,ny),auxcsoil(nx,ny),auxphoto(nx,ny),auxaresp(nx,ny),
     & auxwsoil(nx,ny),auxevap(nx,ny),auxrc(nx,ny)
!      real*4 auxi(nx)
c
c set fire on/off
      mfire = 1 !0=off ; 1=on
c
c read environmental variables
      open(13,file='../outputs/ambientais4.bin',
     &        status='old',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      read(13,rec=1) tmin   !oC
      read(13,rec=2) npp    !kgC/m2/y
      read(13,rec=3) snpp   !dimensionless
      read(13,rec=4) hresp  !kgC/m2/y
      read(13,rec=5) csoil  !kgC/m2/y
      read(13,rec=6) photo  !kgC/m2/y
      read(13,rec=7) aresp  !kgC/m2/y
      read(13,rec=8) ave_wsoil  !mm
      read(13,rec=9) ave_evap  !mm/d
      read(13,rec=10) ave_rc  !s/m

      close(13)

***** MDO BEGIN
c     Read IPAR file and calculate annual average IPAR (avgpar).
      open(15,file='../inputs/ipar.bin',
     &        status='old',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      call read12(15,ipar)
      close(15)

      do i=1,nx
      do j=1,ny
         avgpar(i,j) = 0.0
         do k=1,12
            avgpar(i,j) = avgpar(i,j) + ipar(i,j,k)/12.0
         enddo
      enddo
      enddo
***** MDO END
c
c read soil moisture (W, degree of saturation)
      open(14,file='../outputs/soilm.bin',
     &        status='old',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      call read12(14,wsoil2)
      close(14)
c
c read uwind (850 hPa)
      open(15,file='../inputs/uwnd.bin',
     &        status='old',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      call read12(15,u)
      close(15)
c
! read uwind anomalies (see line 71)
!      open(18,file='../ipcc_ar4/anomalias/SRESA2/uwind/'//
!     &               'HADCM3_A2_2070_2099_anom_u850.bin',status='old',
!     &        form='unformatted',access='direct',recl=4*nx*ny)
!      call read12(18,anu)
!      close(18)
c
c read lsmk
      open(16,file='../inputs/lsmk.bin',status='old',
     &        form='unformatted',access='direct',recl=4*nx*ny)
      read(16,rec=1) lsmk
c
c dry to wet season transition
      do i=1,nx
      do j=1,ny
         d2wm(i,j) = -999.99
         ud2w(i,j) = -999.99
ccc	 npp(i,j) = meanpp(i,j)*0.58	!test for reduction of NPP (Lapola's dissertation)
         if (wsoil2(i,j,1).ge.0.0) then
            do k=1,12
               waux2(k) = wsoil2(i,j,k)
	       uwind(i,j,k) = u(i,j,k)!+anu(i,j,k)
            enddo
            call findm (waux2, moaux2) !find month
            d2wm(i,j) = real(moaux2)
            mr = moaux2 + 1 !month right after d2w transition month
!!!            mr = moaux2 !d2w month
            if (mr.eq.13) mr = 1
            ud2w(i,j) = uwind(i,j,mr) !dry to wet U-wind
         endif
      enddo
      enddo
c
c for all grid points
      do i=1,nx
      do j=1,nx
      nbio = int(lsmk(i,j))
c
c ocean grid points
      if (nbio.eq.0) then
      bpot(i,j) = 0.
      goto 1000
      endif
c
c land grid points
      bpot(i,j) = -1.

***** MDO BEGIN
* Rules for biome allocation.

c desert and semi-desert
      if (tmin(i,j).ge.-15.0) then
         if (tmin(i,j).le.10.0) then !extratropical region
            if (ratio(npp(i,j),avgpar(i,j)).le.0.05) then
               bpot(i,j) = 11.0 !desert
            else
               if (ratio(npp(i,j),avgpar(i,j)).le.0.28) bpot(i,j) = 9.0 !semi-desert
            endif
         else !tropical region
            if (ratio(npp(i,j),avgpar(i,j)).le.0.005) then
               bpot(i,j) = 11.0 !desert
            else
               if (ratio(npp(i,j),avgpar(i,j)).le.0.15) bpot(i,j) = 9.0 !semi-desert
            endif
         endif
      endif

      if (bpot(i,j).eq.-1.0) then !remaining land grid points...

c ice, tundra, larch
      if (tmin(i,j).le.-27.0) then
         if (ratio(npp(i,j),avgpar(i,j)).eq.0.0) then
            bpot(i,j) = 20.0 !ice
         else
            if (ratio(npp(i,j),avgpar(i,j)).le.0.3) then
               bpot(i,j) = 10.0 !tundra
            else
               bpot(i,j) = 5.0 !larch
            endif
         endif
      endif

c remaining tundra, boreal forest, grasslands
      if (tmin(i,j).gt.-27.0 .and. tmin(i,j).le.-9.0) then
         if (ratio(npp(i,j),avgpar(i,j)).le.0.08) then
            bpot(i,j) = 10.0 !remaining tundra
         else
            if (avgpar(i,j).le.58.0) then
               bpot(i,j) = 4.0 !boreal forest
            else
               bpot(i,j) = 7.0 !grasslands
            endif
         endif
      endif

c mixed forest, grasslands
      if (tmin(i,j).gt.-9.0 .and. tmin(i,j).le.-6.0) then
         if (ratio(npp(i,j),avgpar(i,j)).le.0.5) then
c         if (ratio(npp(i,j),avgpar(i,j)).le.0.6) then
            bpot(i,j) = 7.0 !grasslands
         else
            bpot(i,j) = 3.0 !mixed forest
         endif
      endif

c temperate forest, grasslands
      if (tmin(i,j).gt.-6.0 .and. tmin(i,j).le.10.0) then
         if (ratio(npp(i,j),avgpar(i,j)).le.0.7) then
            bpot(i,j) = 7.0 !grasslands
         else
            bpot(i,j) = 2.0 !temperate forest
         endif
      endif

c tropical biomes
      if (tmin(i,j).gt.10.0) then
         if (snpp(i,j).gt.0.1) then
            if (ratio(npp(i,j),avgpar(i,j)).le.0.38) then
               bpot(i,j) = 8.0 !caatinga
            else
               bpot(i,j) = 6.0 !savanna
            endif
         else
            if (ratio(npp(i,j),avgpar(i,j)).le.0.72) then
                bpot(i,j) = 13.0 !seasonal forest
            else
                bpot(i,j) = 1.0 !ombrophyllous forest
            endif
         endif
      endif

      endif
*
***** MDO END
*
      if (mfire.eq.1) then !fire parameterization is on
***** Natural fire parameterization from U-wind (to correct India)
      ulim = 1.5 !obtained after the best results with observed lightning
         if (int(bpot(i,j)).eq.6.) then
            if (ud2w(i,j).le.ulim) then !fire!!!
            bpot(i,j) = 6.0
	    else !no fire!!!
            bpot(i,j) = 13.0
	    endif
         endif
      else !fire parameterization is off
      continue
      endif
c
      if (bpot(i,j).eq.-1.0) then
         print *, 'Error: no biome grid cell!'
         stop
      endif
c
 1000 continue
c
      enddo
      enddo

! loop inserted to have maps atarting at longitude -180 and with
! inverted latitude (to read maps in ArcGIS [11/12/2013]
      do i=1,nx
         do j=1,ny
         if (i.ge.nx/2) then
         auxbpot(i,j) = bpot((i+1)-(nx/2),(ny+1)-j)
         auxnpp(i,j) = npp((i+1)-(nx/2),(ny+1)-j)
         auxsnpp(i,j) = snpp((i+1)-(nx/2),(ny+1)-j)
         auxhresp(i,j) = hresp((i+1)-(nx/2),(ny+1)-j)
         auxcsoil(i,j) = csoil((i+1)-(nx/2),(ny+1)-j)
         auxphoto(i,j) = photo((i+1)-(nx/2),(ny+1)-j)
         auxaresp(i,j) = aresp((i+1)-(nx/2),(ny+1)-j)
         auxwsoil(i,j) = ave_wsoil((i+1)-(nx/2),(ny+1)-j)
         auxevap(i,j) = ave_evap((i+1)-(nx/2),(ny+1)-j)
         auxrc(i,j) = ave_rc((i+1)-(nx/2),(ny+1)-j)
         else
         auxbpot(i,j) = bpot((i+1)+(nx/2),(ny+1)-j)
         auxnpp(i,j) = npp((i+1)+(nx/2),(ny+1)-j)
         auxsnpp(i,j) = snpp((i+1)+(nx/2),(ny+1)-j)
         auxhresp(i,j) = hresp((i+1)+(nx/2),(ny+1)-j)
         auxcsoil(i,j) = csoil((i+1)+(nx/2),(ny+1)-j)
         auxphoto(i,j) = photo((i+1)+(nx/2),(ny+1)-j)
         auxaresp(i,j) = aresp((i+1)+(nx/2),(ny+1)-j)
         auxwsoil(i,j) = ave_wsoil((i+1)+(nx/2),(ny+1)-j)
         auxevap(i,j) = ave_evap((i+1)+(nx/2),(ny+1)-j)
         auxrc(i,j) = ave_rc((i+1)+(nx/2),(ny+1)-j)
         endif
         enddo
      enddo
c
      if (mfire.eq.0) then
      open(17,file='../outputs/bpot4off.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      else
      open(17,file='../outputs/bpot4on.flt',
c      open(17,file='../ipcc_env_vars2/bpot4on_'//
c     &    'HADCM3_A2_co2onlyRc.bin',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      endif
      write(17,rec=1) auxbpot
!      write(17,rec=2) npp
!      write(17,rec=3) snpp
!      write(17,rec=4) hresp
!      write(17,rec=5) csoil
!      write(17,rec=6) d2wm
!      write(17,rec=7) ud2w
!      write(17,rec=8) photo
!      write(17,rec=9) aresp
c
      close(17)
!
! Save variables in different files for the PVM2 workshop [11/12/13]
      open(18,file='../outputs/npp.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(18,rec=1) auxnpp
      close(18)
c
      open(19,file='../outputs/snpp.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(19,rec=1) auxsnpp
      close(19)
c
      open(20,file='../outputs/hresp.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(20,rec=1) auxhresp
      close(20)
c
      open(21,file='../outputs/csoil.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(21,rec=1) auxcsoil
      close(21)
c
      open(22,file='../outputs/photo.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(22,rec=1) auxphoto
      close(22)
c
      open(23,file='../outputs/aresp.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(23,rec=1) auxaresp
      close(23)
c
      open(24,file='../outputs/wsoil.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(24,rec=1) auxwsoil
      close(24)
c
      open(25,file='../outputs/ave_evap.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(25,rec=1) auxevap
      close(25)
c
      open(26,file='../outputs/ave_rc.flt',
     &        status='unknown',form='unformatted',
     &        access='direct',recl=4*nx*ny)
      write(26,rec=1) auxrc
      close(26)
c
      stop
      end
c
c=======================================================================
c
***** MDO BEGIN
* Ratio between NPP and maximum reference NPP (from IPAR).
      real function ratio (anpp, apar)
      real anpp, apar
      real refval
      refval = amax1((apar - 35.0)/45.0,0.0)
      if (refval .eq. 0.0) then
         ratio = 0.0 !independente de PAR
      else
         ratio = anpp/refval
         ratio = amax1(ratio,0.0)
         ratio = amin1(ratio,1.0)
      endif
      return
      end
***** MDO END
c
c=======================================================================
c
      subroutine read12 (nunit,var)
c auxiliar reading routine
      parameter(nx=192,ny=96)
      integer nunit
      real*4 var(nx,ny,12)
      real*4 aux(nx,ny)
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

c=======================================================================

      subroutine findm (w, m)

      real*4 w(12), w2(12), w3(12), dw(12)
      integer m, m3(12)

      wmax = -10.0
      wmin =  10.0
      do k=1,12
         w2(k) = w(k)
         if (w(k).gt.wmax) wmax=w(k)
         if (w(k).lt.wmin) wmin=w(k)
      enddo
c em media, deveria ser dividido por 6; a ideia eh que a restricao
c nao seja tao grande
      wlim = (wmax-wmin)/8.0
c      print *, wlim

      do iter=1,12
         wmin = 10.0
         do k=1,12
            if (w2(k).lt.wmin) then
               wmin = w2(k)
               mmin = k
            endif
         enddo
         w3(iter) = wmin
         m3(iter) = mmin
         if (mmin.ne.12) then
            dw(iter) = w(mmin+1)-w(mmin)
         else
            dw(iter) = w(1) - w(12)
         endif
         w2(mmin) = 100.0
      enddo
 
      iter=0
   10 continue
      iter=iter+1
      if (dw(iter).ge.wlim) then
      m = m3(iter)
      goto 99
      endif
      goto 10
   99 continue
 
      return
      end

c=======================================================================

