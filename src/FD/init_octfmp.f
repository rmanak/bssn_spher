      subroutine init_octfmp(x,Nx,rv,res)
      implicit none
      integer i
      integer Nx
      real*8 rv
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.1D1 / (exp(x(i)) + 0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) ** 2
     # / rv) ** 2 * (exp(x(i)) + 0.2D1 / (0.1D1 - 0.1D1 * x(i) / rv) ** 
     #3 / rv ** 2)
      res(i)=qb
      end do
      END
