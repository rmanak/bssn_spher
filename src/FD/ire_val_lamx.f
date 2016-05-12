      subroutine ire_val_lamx(ctfm,ctfmp,n_A,n_B,n_Lamx,Nx,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_B(Nx)
      real*8 n_Lamx(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = (0.5000000000000000D0 / n_A(i) / n_B(i) * (n_B(i - 1) - 0.1D1
     # * n_B(i + 1)) / hx + 0.2500000000000000D0 / n_A(i) ** 2 * (-0.1D1
     # * n_A(i - 1) + n_A(i + 1)) / hx) / ctfmp(i) + (0.2D1 / n_B(i) - 0
     #.2D1 / n_A(i)) / ctfm(i) - 0.1D1 * n_Lamx(i)
      res(i)=qb
      end do
      END
