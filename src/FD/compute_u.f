      subroutine compute_u(n_ipsi,n_rpsi,Nx,mass,res)
      implicit none
      integer i
      integer Nx
      real*8 mass
      real*8 n_ipsi(Nx)
      real*8 n_rpsi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = mass * (n_rpsi(i) ** 2 + n_ipsi(i) ** 2)
      res(i)=qb
      end do
      END
