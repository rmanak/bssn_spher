      subroutine evol_beta_hgd_adv(ctfmp,n_BB,n_beta,np1_BB,np1_beta,x,N
     &x,advc,ht,hx,mus,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 advc
      real*8 ht
      real*8 hx
      real*8 mus
      real*8 myzero
      real*8 ctfmp(Nx)
      real*8 n_BB(Nx)
      real*8 n_beta(Nx)
      real*8 np1_BB(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = np1_beta(i) - 0.4D1 / (advc * ht * np1_beta(i - 1) - 0.1D1 * 
     #advc * ht * np1_beta(i + 1) + 0.4D1 * hx * ctfmp(i)) * ht * hx * c
     #tfmp(i) * ((-0.1D1 * n_beta(i) + np1_beta(i)) / ht - 0.25000000000
     #00000D0 * advc * np1_beta(i) * (-0.1D1 * np1_beta(i - 1) + np1_bet
     #a(i + 1)) / hx / ctfmp(i) - 0.5000000000000000D0 * mus * np1_BB(i)
     # - 0.2500000000000000D0 * advc * n_beta(i) * (-0.1D1 * n_beta(i - 
     #1) + n_beta(i + 1)) / hx / ctfmp(i) - 0.5000000000000000D0 * mus *
     # n_BB(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = -0.1D1 * myzero * x(i)
      res(i)=qb
      end do
      endif
      END
