      subroutine compute_tt(n_PHI2,n_PI2,n_a2,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_PHI2(Nx)
      real*8 n_PI2(Nx)
      real*8 n_a2(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = (n_PI2(i) - n_PHI2(i)) / n_a2(i)
      res(i)=qb
      end do
      END
