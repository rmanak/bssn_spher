      subroutine compute_a1(n_a2,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_a2(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = sqrt(n_a2(i))
      res(i)=qb
      end do
      END
