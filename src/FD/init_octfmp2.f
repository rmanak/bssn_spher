      subroutine init_octfmp2(x,Nx,lv,rv,res)
      implicit none
      integer i
      integer Nx
      real*8 lv
      real*8 rv
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.1D1 / (exp(x(i)) + (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0
     #.1D1 / (0.1D1 - 0.1D1 * lv / rv)) ** (rv - 0.1D1 * x(i) + 0.1D1) *
     # (-0.1D1 * log(0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D
     #1 - 0.1D1 * lv / rv)) + (rv - 0.1D1 * x(i) + 0.1D1) / (0.1D1 - 0.1
     #D1 * x(i) / rv) ** 2 / rv / (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) -
     # 0.1D1 / (0.1D1 - 0.1D1 * lv / rv)))) ** 2 * (exp(x(i)) + (0.1D1 /
     # (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 0.1D1 * lv / rv)) 
     #** (rv - 0.1D1 * x(i) + 0.1D1) * (-0.1D1 * log(0.1D1 / (0.1D1 - 0.
     #1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 0.1D1 * lv / rv)) + (rv - 0.1D
     #1 * x(i) + 0.1D1) / (0.1D1 - 0.1D1 * x(i) / rv) ** 2 / rv / (0.1D1
     # / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 0.1D1 * lv / rv)
     #)) ** 2 + (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 
     #0.1D1 * lv / rv)) ** (rv - 0.1D1 * x(i) + 0.1D1) * (-0.2D1 / (0.1D
     #1 - 0.1D1 * x(i) / rv) ** 2 / rv / (0.1D1 / (0.1D1 - 0.1D1 * x(i) 
     #/ rv) - 0.1D1 / (0.1D1 - 0.1D1 * lv / rv)) + 0.2D1 * (rv - 0.1D1 *
     # x(i) + 0.1D1) / (0.1D1 - 0.1D1 * x(i) / rv) ** 3 / rv ** 2 / (0.1
     #D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D1 / (0.1D1 - 0.1D1 * lv / r
     #v)) - 0.1D1 * (rv - 0.1D1 * x(i) + 0.1D1) / (0.1D1 - 0.1D1 * x(i) 
     #/ rv) ** 4 / rv ** 2 / (0.1D1 / (0.1D1 - 0.1D1 * x(i) / rv) - 0.1D
     #1 / (0.1D1 - 0.1D1 * lv / rv)) ** 2))
      res(i)=qb
      end do
      END
