      subroutine init_iphi(ctfmp,x,Nx,myzero,res)
      implicit none
      integer i
      integer Nx
      real*8 myzero
      real*8 ctfmp(Nx)
      real*8 x(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.2D1 * myzero * x(i) / ctfmp(i)
      res(i)=qb
      end do
      END
