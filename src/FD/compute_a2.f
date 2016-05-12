      subroutine compute_a2(n_A,n_psi,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_A(Nx)
      real*8 n_psi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_psi(i) ** 4 * n_A(i)
      res(i)=qb
      end do
      END
