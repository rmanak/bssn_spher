      subroutine compute_ts(n_PHI2,n_PI2,n_U,n_a2,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_PHI2(Nx)
      real*8 n_PI2(Nx)
      real*8 n_U(Nx)
      real*8 n_a2(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.5000000000000000D0 * (0.3D1 * n_PI2(i) - 0.1D1 * n_PHI2(i))
     # / n_a2(i) - 0.1500000000000000D1 * n_U(i)
      res(i)=qb
      end do
      END
