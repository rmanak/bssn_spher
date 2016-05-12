      subroutine compute_b1(n_B,n_psi,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_B(Nx)
      real*8 n_psi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_psi(i) ** 2 * sqrt(n_B(i))
      res(i)=qb
      end do
      END
