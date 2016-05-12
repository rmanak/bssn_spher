      subroutine init_ctfmp2(x,Nx,lv,rv,res)
      implicit none
      integer i
      integer Nx
      real*8 lv
      real*8 rv
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = exp(x(i)) + (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0
     #.1D1 - 0.1D1 * lv / rv)) ** (rv - 0.1D1 * x(i) + 0.1D1) * (-0.1D1 
     #* log(0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 0.1D1
     # * lv / rv)) + (rv - 0.1D1 * x(i) + 0.1D1) / (0.1D1 - 0.1D1 * x(i)
     # / rv) ** 2 / rv / (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / 
     #(0.1D1 - 0.1D1 * lv / rv)))
      res(i)=qb
      end do
      END
