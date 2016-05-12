      subroutine init_rpi(n_a1,n_alpha,n_beta,n_rPHI,n_rpsidot,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_rPHI(Nx)
      real*8 n_rpsidot(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_a1(i) / n_alpha(i) * (n_rpsidot(i) - 0.1D1 * n_beta(i) * n_
     #rPHI(i))
      res(i)=qb
      end do
      END
