      subroutine compute_psi(n_phi,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_phi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = exp(n_phi(i))
      res(i)=qb
      end do
      END
