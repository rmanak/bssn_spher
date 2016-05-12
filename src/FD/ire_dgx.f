      subroutine ire_dgx(ctfmp,n_A,n_Axx,n_alpha,n_beta,n_divbeta,nm1_A,
     &np1_A,Nx,ht,hx,vee,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 vee
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Axx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 nm1_A(Nx)
      real*8 np1_A(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (nm1_A(i) - 0.1D1 * np1_A(i)) / ht + 
     #0.2D1 * n_Axx(i) * n_alpha(i) + 0.6666666666666667D0 * n_A(i) * ve
     #e * n_divbeta(i) - 0.1D1 * (-0.1D1 * (n_beta(i - 1) - 0.1D1 * n_be
     #ta(i + 1)) / hx * n_A(i) + 0.5000000000000000D0 * n_beta(i) * (-0.
     #1D1 * n_A(i - 1) + n_A(i + 1)) / hx) / ctfmp(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
