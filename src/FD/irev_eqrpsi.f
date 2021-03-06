      subroutine irev_eqrpsi(n_a1,n_alpha,n_beta,n_rPHI,n_rPI,nm1_rpsi,n
     &p1_rpsi,Nx,ht,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_rPHI(Nx)
      real*8 n_rPI(Nx)
      real*8 nm1_rpsi(Nx)
      real*8 np1_rpsi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) /
     # ht - 0.1D1 * n_alpha(i) / n_a1(i) * n_rPI(i) - 0.1D1 * n_beta(i) 
     #* n_rPHI(i)
      res(i)=qb
      end do
      END
