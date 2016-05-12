      subroutine init_rpsi(ctfm,Nx,amp,deltx,x0,res)
      implicit none
      integer i
      integer Nx
      real*8 amp
      real*8 deltx
      real*8 x0
      real*8 ctfm(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = amp * exp(-0.1D1 * (ctfm(i) - 0.1D1 * x0) ** 2 / deltx ** 2)
      res(i)=qb
      end do
      END
