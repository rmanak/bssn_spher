      subroutine compute_c5psi(n_C5s,n_psi,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_C5s(Nx)
      real*8 n_psi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_C5s(i) * n_psi(i) ** 5
      res(i)=qb
      end do
      END
