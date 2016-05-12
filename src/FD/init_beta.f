      subroutine init_beta(x,Nx,myzero,res)
      implicit none
      integer i
      integer Nx
      real*8 myzero
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      END
