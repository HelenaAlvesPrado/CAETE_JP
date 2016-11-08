c234567
c======================================================================

      program desempenho

c=====================================================================

      parameter (nmax = 12)

      real*4 bobs(720,360),bmod(720,360)

      real r(nmax)

c=====================================================================

      open(11,file='Lapola_et_al_T062.bin',status='old',
     &        form='unformatted',access='direct',recl=720*360*4)
      open(12,file='bpot4off.bin', status='old',
     &        form='unformatted',access='direct',recl=720*360*4)

      read(11,rec=1) bobs  !observado
      read(12,rec=1) bmod  !modelo

      close(11)
      close(12)

c floresta trop. estacional passa a ser o bioma 12
      do i=1,720
      do j=1,360
      if (int(bobs(i,j)).eq.13) bobs(i,j) = 12.
      if (int(bmod(i,j)).eq.13) bmod(i,j) = 12.
      enddo
      enddo
c      write(*,*) '================= ATENCAO ================='
c      write(*,*) '  bioma 12 eh floresta tropical estacional '
c      write(*,*) '==========================================='
c      write(*,*)


c gelo passa a ser o bioma 13
      do i=1,720
      do j=1,360
      if (int(bobs(i,j)).eq.20) bobs(i,j) = 13.
      if (int(bmod(i,j)).eq.20) bmod(i,j) = 13.
      enddo
      enddo
c      write(*,*) '===== ATENCAO ====='
c      write(*,*) '  bioma 13 eh gelo '
c      write(*,*) '==================='
c      write(*,*)

c=====================================================================
c indice de acerto

c      write(*,*) ' indice de acerto '
      do nbio=1,nmax
      call acerto (nbio,bobs,bmod,raux)
      r(nbio) = raux
      enddo
c      write(*,*)

c=====================================================================

      call kappa (bobs,bmod,akbio,akappa,r,nmax)

c=====================================================================

      stop
      end

c=====================================================================
      subroutine acerto (nbio,bobs,bmod,raux)
      integer nbio
      real*4 bobs(720,360),bmod(720,360)
      real raux
      ntotal = 0
      ncerto = 0
      do i=1,720
      do j=1,360
      if (int(bobs(i,j)).eq.nbio) then
        ntotal = ntotal + 1
        if (int(bobs(i,j)).eq.int(bmod(i,j))) ncerto = ncerto + 1
      endif
      enddo
      enddo
      raux = 100.*real(ncerto)/real(ntotal)
c      write(*,11) nbio,ntotal,ncerto,100.*real(ncerto)/real(ntotal)
c   11 format(i3,2i5,f6.1)
      return
      end
c=====================================================================
      subroutine kappa (bobs,bmod,akbio,akappa,r,nmax)

c kappa statistics

      real*4 bobs(720,360),bmod(720,360)
      real akbio(nmax),akappa,r(nmax)

      integer nfreq(nmax,nmax),ns_obs(nmax),ns_mod(nmax)
      real prob(nmax,nmax),p_obs(nmax),p_mod(nmax)
      real peso(nmax)

c frequencia

      do nobs=1,nmax
      do nmod=1,nmax
        nfreq(nobs,nmod) = 0
      enddo
      enddo

      nsum = 0
      do nobs=1,nmax
        do i=1,720
        do j=1,360
          if (int(bobs(i,j)).eq.nobs) then
            do nmod=1,nmax
              if (int(bmod(i,j)).eq.nmod) then
                nfreq(nobs,nmod) = nfreq(nobs,nmod) + 1
                nsum = nsum + 1
              endif
            enddo
          endif
        enddo
        enddo
      enddo

c      write(*,*) ' matriz de concordancia (linhas=obs;colunas=mod) '
c      write(*,8) (i,i=1,11)
      do nobs=1,nmax
        ns_obs(nobs) = 0 !linha fixa: soma ao longo das colunas
        do nmod=1,nmax
          ns_obs(nobs) = ns_obs(nobs) + nfreq(nobs,nmod)
        enddo
c        write(*,9) nobs,(nfreq(nobs,nmod),nmod=1,nmax),ns_obs(nobs)
      enddo

      do nmod=1,nmax
        ns_mod(nmod) = 0 !coluna fixa: soma ao longo das linhas
        do nobs=1,nmax
          ns_mod(nmod) = ns_mod(nmod) + nfreq(nobs,nmod)
        enddo
      enddo
c      write(*,10) (ns_mod(nmod),nmod=1,nmax),nsum
c      write(*,*)

    8 format(3x,11i5)
    9 format(i2,1x,11i5,1x,'-',1x,i5)
   10 format(3x,11i5,1x,'-',1x,i5)

      do nobs=1,nmax
      do nmod=1,nmax
        prob(nobs,nmod) = real(nfreq(nobs,nmod))/real(nsum)
      enddo
      enddo

      do n=1,nmax
        p_obs(n) = real(ns_obs(n))/real(nsum)
        p_mod(n) = real(ns_mod(n))/real(nsum)
      enddo

      write(*,*) ' indice de acerto (%) e kappa por bioma '
      p_0 = 0.
      p_e = 0.
      do n=1,nmax
        p_0 = p_0 + prob(n,n)
        vezes = p_obs(n)*p_mod(n)
        p_e = p_e + vezes
        akbio(n) = ( prob(n,n) - vezes )/
     &             ( ((p_obs(n)+p_mod(n))/2.) - vezes )
      write(*,11) n,r(n),akbio(n)
      enddo
      write(*,*)

      write(*,*) ' kappa global '
      akappa = (p_0 - p_e)/(1. - p_e)
      write(*,12) ' p0    = ',p_0
      write(*,12) ' kappa = ',akappa
      write(*,*)

   11 format(i2,1x,f6.1,f6.2)
   12 format(a9,f5.2)

c teste (kappa global como uma media ponderada)
c conclusao: OK
      speso = 0.
      akappa2 = 0.
      do n=1,nmax
      peso(n) = ((p_obs(n)+p_mod(n))/2.) - p_obs(n)*p_mod(n)
      speso = speso + peso(n)
      akappa2 = akappa2 + akbio(n)*peso(n)
      enddo
      akappa2 = akappa2/speso
c      write(*,*) akappa,akappa2 

      return
      end
