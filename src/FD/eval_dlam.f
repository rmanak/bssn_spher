      subroutine eval_dlam(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,n
     &_Sx,n_alpha,n_beta,n_divbeta,n_phi,octfmp,x,Nx,hx,myzero,vee,phys_
     &bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 vee
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_Lamx(Nx)
      real*8 n_Sx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 n_phi(Nx)
      real*8 octfmp(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = 0.6666666666666667D0 * n_divbeta(i) * n_Lamx(i) * vee - 0.2D1
     # / n_A(i) * n_alpha(i) * n_Sx(i) + (-0.4D1 / n_B(i) / n_A(i) * n_A
     #thth(i) * n_alpha(i) + 0.4D1 / n_B(i) ** 2 * n_alpha(i) * n_Athth(
     #i)) / ctfm(i) - 0.2D1 / ctfm(i) ** 2 / n_B(i) * n_beta(i) + ((0.50
     #00000000000000D0 / n_A(i) ** 3 * n_Axx(i) * n_alpha(i) * (-0.1D1 *
     # n_A(i - 1) + n_A(i + 1)) - 0.1D1 / n_B(i) ** 2 / n_A(i) * (-0.1D1
     # * n_B(i - 1) + n_B(i + 1)) * n_Athth(i) * n_alpha(i) + 0.66666666
     #66666667D0 / n_A(i) * n_alpha(i) * (n_K(i - 1) - 0.1D1 * n_K(i + 1
     #)) + 0.5000000000000000D0 * n_beta(i) * (-0.1D1 * n_Lamx(i - 1) + 
     #n_Lamx(i + 1)) + 0.1D1 / n_A(i) ** 2 * n_Axx(i) * (n_alpha(i - 1) 
     #- 0.1D1 * n_alpha(i + 1)) - 0.5000000000000000D0 * n_Lamx(i) * (-0
     #.1D1 * n_beta(i - 1) + n_beta(i + 1)) - 0.5000000000000000D0 * (-0
     #.1D1 * n_beta(i - 1) + n_beta(i + 1)) / n_A(i) * octfmp(i) - 0.166
     #6666666666667D0 / n_A(i) * vee * (n_divbeta(i - 1) - 0.1D1 * n_div
     #beta(i + 1)) - 0.6D1 / n_A(i) ** 2 * n_alpha(i) * n_Axx(i) * (n_ph
     #i(i - 1) - 0.1D1 * n_phi(i + 1))) / ctfmp(i) + 0.1D1 / n_B(i) * (-
     #0.1D1 * n_beta(i - 1) + n_beta(i + 1)) / ctfmp(i) / ctfm(i)) / hx 
     #- 0.1D1 / n_A(i) / ctfmp(i) ** 2 * (-0.1D1 * n_beta(i - 1) + 0.2D1
     # * n_beta(i) - 0.1D1 * n_beta(i + 1)) / hx ** 2
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
