program pvm5
    implicit none
! c
! c=====================================================================
! c
! c Adapted from pvm4.f
! c Code written by David Lapola & Marcos Oyama
! c Last updates
! c MDO: 24 Aug 2007.
! c DML: 27 Aug 2007.
! c
! c=====================================================================
! c
! c i/o variables
! c d2wm : month of the transition between wet and dry season
! c
    integer  i, j, k, mfire, moaux2, iter
    real ratio, ulim
    integer, parameter :: nx=720, ny=360, nt=12
    real, parameter :: NO_DATA = -9999.000000
    real, dimension(nt) :: waux2
    real, dimension(nx,ny,nt) :: wsoil2, ipar, u, uwind

    real, dimension(nx,ny) :: tmin, npp, snpp, hresp, csoil, photo,&
                              aresp, ave_wsoil, ave_evap, ave_rc,&
                              lsmk, bpot, avgpar, d2wm, ud2w, auxbpot,&
                              auxnpp, auxsnpp, auxhresp, auxcsoil, auxphoto,&
                              auxaresp, auxwsoil, auxevap, auxrc


! c
! c set fire on/off
       mfire = 1 !0=off ; 1=on
! c
!c read environmental variables
    open(13,file='../outputs/ambientais4.bin', status='old',form='unformatted',&
         access='direct',recl=4*nx*ny)
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
!
! c***** MDO BEGIN
! c     Read IPAR file and calculate annual average IPAR (avgpar).
    open(15,file='../inputs/ipar.bin',&
        status='old',form='unformatted',&
        access='direct',recl=4*nx*ny)

        call read12(15,ipar)
    close(15)

    do i=1,nx
        do j=1,ny
            avgpar(i,j) = 0.0
            do k=1,12
                avgpar(i,j) = avgpar(i,j) + ipar(i,j,k)/12.
            enddo
        enddo
    enddo
!================================================================


    open(14,file='../outputs/soilm.bin',&
        status='old',form='unformatted',&
        access='direct',recl=4*nx*ny)

        call read12(14,wsoil2)
    close(14)

    open(15,file='../inputs/uwnd.bin',&
        status='old',form='unformatted',&
        access='direct',recl=4*nx*ny)

        call read12(15,u)
    close(15)

    open(16,file='../inputs/lmsk.bin',status='old',&
         form='unformatted',access='direct',recl=4*nx*ny)

         read(16,rec=1) lsmk
    close(16)
! c
! c dry to wet season transition
    do i=1,nx
        do j=1,ny
            d2wm(i,j) = NO_DATA
            ud2w(i,j) = NO_DATA
!ccc     npp(i,j) = meanpp(i,j)*0.58    !test for reduction of NPP (Lapola's dissertation)
            if (wsoil2(i,j,1) >= 0.0) then
                do k=1,nt
                    waux2(k) = wsoil2(i,j,k)
                    uwind(i,j,k) = u(i,j,k)
                enddo
                call findm (waux2, moaux2) !find month
                d2wm(i,j) = real(moaux2)
                mr = moaux2 + 1 !month right after d2w transition month
!!!            mr = moaux2 !d2w month
                if (mr == 13) mr = 1
                ud2w(i,j) = uwind(i,j,mr) !dry to wet U-wind
            endif
        enddo
    enddo
! c
! c for all grid points
    do i=1,nx
        do j=1,ny
            nbio = nint(lsmk(i,j))
            ! c
            ! c ocean grid points
            if (nbio == 0) then
                bpot(i,j) = NO_DATA
                goto 1000
            endif
            ! c
            ! c land grid points
            bpot(i,j) = -1.

! ***** MDO BEGIN
! * Rules for biome allocation.
!
! c desert and semi-desert
        if (tmin(i,j) >= -15.0) then
            if (tmin(i,j) <= 10.0) then !extratropical region
                if (ratio(npp(i,j),avgpar(i,j)) <= 0.05) then
                    bpot(i,j) = 11.0 !desert
                else
                    if (ratio(npp(i,j),avgpar(i,j)) <= 0.28) bpot(i,j) = 9.0 !semi-desert
                endif
            else !tropical region
                if (ratio(npp(i,j),avgpar(i,j)) <= 0.005) then
                    bpot(i,j) = 11.0 !desert
                else
                    if (ratio(npp(i,j),avgpar(i,j)) <= 0.15) bpot(i,j) = 9.0 !semi-desert
                endif
            endif
        endif

        if (bpot(i,j) ==  -1.0) then !remaining land grid points...
            !c ice, tundra, larch
            if (tmin(i,j) <= -27.0) then
                if (ratio(npp(i,j),avgpar(i,j)) ==  0.0) then
                bpot(i,j) = 20.0 !ice
            else
                if (ratio(npp(i,j),avgpar(i,j)) <= 0.3) then
                    bpot(i,j) = 10.0 !tundra
                else
                    bpot(i,j) = 5.0 !larch
                endif
            endif
        endif

        !c remaining tundra, boreal forest, grasslands
        if (tmin(i,j) > .-27.0 .and. tmin(i,j) <= -9.0) then
            if (ratio(npp(i,j),avgpar(i,j)) <= 0.08) then
                bpot(i,j) = 10.0 !remaining tundra
             else
                 if (avgpar(i,j) <= 58.0) then
                     bpot(i,j) = 4.0 !boreal forest
                 else
                     bpot(i,j) = 7.0 !grasslands
                 endif
             endif
         endif

         !c mixed forest, grasslands
        if (tmin(i,j) > .-9.0 .and. tmin(i,j) <= -6.0) then
            if (ratio(npp(i,j),avgpar(i,j)) <= 0.5) then
                bpot(i,j) = 7.0 !grasslands
            else
                bpot(i,j) = 3.0 !mixed forest
            endif
        endif

        !c temperate forest, grasslands
        if (tmin(i,j) > .-6.0 .and. tmin(i,j) <= 10.0) then
            if (ratio(npp(i,j),avgpar(i,j)) <= 0.7) then
                bpot(i,j) = 7.0 !grasslands
            else
                bpot(i,j) = 2.0 !temperate forest
            endif
        endif

        !c tropical biomes
        if (tmin(i,j) > .10.0) then
            if (snpp(i,j) > .0.1) then
                if (ratio(npp(i,j),avgpar(i,j)) <= 0.38) then
                    bpot(i,j) = 8.0 !caatinga
                else
                    bpot(i,j) = 6.0 !savanna
                endif
             else
                 if (ratio(npp(i,j),avgpar(i,j)) <= 0.72) then
                     bpot(i,j) = 13.0 !seasonal forest
                 else
                     bpot(i,j) = 1.0 !ombrophyllous forest
                 endif
             endif
         endif
    ! *
    ! ***** MDO END
    ! *
        if (mfire ==  1) then !fire parameterization is on
    !***** Natural fire parameterization from U-wind (to correct India)
            ulim = 1.5 !obtained after the best results with observed lightning
            if (nint(bpot(i,j)) ==  6) then
                if (ud2w(i,j) <= ulim) then !fire!!!
                    bpot(i,j) = 6.0
                else !no fire!!!
                    bpot(i,j) = 13.0
                endif
            endif
        else !fire parameterization is off
            continue
        endif

        if (nint(bpot(i,j)) ==  -1) then
            print *, 'Error: no biome grid cell!'
            stop
        endif

1000 continue
        enddo
    enddo

! whiting outputs
    if (mfire ==  0) then
        open(17,file='../outputs/bpot4off.bin',&
            status='unknown',form='unformatted',&
            access='direct',recl=4*nx*ny)
    else
        open(17,file='../outputs/bpot4on.bin',&
            status='unknown',form='unformatted',&
            access='direct',recl=4*nx*ny)
    endif
    write(17,rec=1) bpot
    write(17,rec=2) npp
    write(17,rec=3) snpp
    write(17,rec=4) hresp
    write(17,rec=5) csoil
    write(17,rec=6) d2wm
    write(17,rec=7) ud2w
    write(17,rec=8) photo
    write(17,rec=9) aresp

    close(17)
    stop "end program"

contains

! * Ratio between NPP and maximum reference NPP (from IPAR).
    real function ratio (anpp, apar)
        real, intent(in) :: anpp, apar
        real, intent(out) :: ratio
        real refval
        refval = amax1((apar - 35.0)/45.0,0.0)
      if (refval == 0.0) then
          ratio = 0.0 !independente de PAR
      else
          ratio = anpp/refval
          ratio = amax1(ratio,0.0)
          ratio = amin1(ratio,1.0)
      endif
      return
  end function ratio

subroutine read12(nunit,var)
      implicit none
  ! c auxiliar reading routine
      integer, parameter :: nx=720,ny=360, nt=12
      integer, intent(in) :: nunit
      real, dimension(nx,ny,12), intent(out) :: var
      integer :: i,j,k
      real, dimension(nx,ny) :: aux

      do k=1,nt
          read(nunit,rec=k) aux
          do i=1,nx
              do j=1,ny
                  var(i,j,k) = aux(i,j)
              enddo
          enddo
       enddo
       return
end subroutine read12
c=======================================================================
subroutine findm (w, m)

    real, dimension(12), intent(in) :: w
    integer, intent(out) :: m
    real, dimension(12) :: w2, w3, dw
    integer, dimenssion(12) :: m3

    real wmax, wmin, wlim
    integer mmin


    wmax = -10.0
    wmin =  10.0
    do k=1,12
        w2(k) = w(k)
        if (w(k) > wmax) wmax=w(k)
        if (w(k) < wmin) wmin=w(k)
    enddo

! c em media, deveria ser dividido por 6; a ideia eh que a restricao
! c nao seja tao grande
    wlim = (wmax-wmin)/8.0
         !print *, wlim

    do iter=1,12
        wmin = 10.0
        do k=1,12
            if (w2(k) < wmin) then
                wmin = w2(k)
                mmin = k
            endif
        enddo
        w3(iter) = wmin
        m3(iter) = mmin
        if (mmin /= 12) then
            dw(iter) = w(mmin+1)-w(mmin)
        else
            dw(iter) = w(1) - w(12)
        endif
        w2(mmin) = 100.0
    enddo

    iter=0
10 continue
    iter=iter+1
    if (dw(iter) >= wlim) then
        m = m3(iter)
        goto 99
    endif
    goto 10
99 continue
    return
end subroutine findm


stop
end
