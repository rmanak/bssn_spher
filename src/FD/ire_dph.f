      subroutine ire_dph(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,nm1_ph
     &i,np1_phi,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 ctfmp(Nx)
      real*8 n_K(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 n_phi(Nx)
      real*8 nm1_phi(Nx)
      real*8 np1_phi(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_phi(i) + np1_phi(i)) / h
     #t + 0.5000000000000000D0 / ctfmp(i) * (n_phi(i - 1) - 0.1D1 * n_ph
     #i(i + 1)) / hx * n_beta(i) + 0.1666666666666667D0 * n_alpha(i) * n
     #_K(i) - 0.1666666666666667D0 * n_divbeta(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
