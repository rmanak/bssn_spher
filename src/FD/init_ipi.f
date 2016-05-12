      subroutine init_ipi(n_a1,n_alpha,n_beta,n_iPHI,n_ipsidot,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_iPHI(Nx)
      real*8 n_ipsidot(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_a1(i) / n_alpha(i) * (n_ipsidot(i) - 0.1D1 * n_beta(i) * n_
     #iPHI(i))
      res(i)=qb
      end do
      END
