      subroutine ire_eqipsi(n_a1,n_alpha,n_beta,n_iPHI,n_iPI,nm1_ipsi,np
     &1_ipsi,Nx,ht,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_iPHI(Nx)
      real*8 n_iPI(Nx)
      real*8 nm1_ipsi(Nx)
      real*8 np1_ipsi(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=1, Nx, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_ipsi(i) + np1_ipsi(i)) /
     # ht - 0.1D1 * n_alpha(i) / n_a1(i) * n_iPI(i) - 0.1D1 * n_beta(i) 
     #* n_iPHI(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
