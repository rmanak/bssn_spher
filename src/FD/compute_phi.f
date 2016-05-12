      subroutine compute_phi(n_psi,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_psi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = log(n_psi(i))
      res(i)=qb
      end do
      END
