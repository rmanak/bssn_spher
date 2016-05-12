      subroutine init_rphi(ctfm,ctfmp,Nx,amp,deltx,hx,x0,res)
      implicit none
      integer i
      integer Nx
      real*8 amp
      real*8 deltx
      real*8 hx
      real*8 x0
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = amp * (ctfm(i) - 0.1D1 * x0) / deltx ** 2 * (ctfm(i - 1) - 0.
     #1D1 * ctfm(i + 1)) / hx * exp(-0.1D1 * (ctfm(i) - 0.1D1 * x0) ** 2
     # / deltx ** 2) / ctfmp(i)
      res(i)=qb
      end do
      END
