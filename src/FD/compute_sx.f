      subroutine compute_sx(n_a1,n_iPHI,n_iPI,n_rPHI,n_rPI,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_a1(Nx)
      real*8 n_iPHI(Nx)
      real*8 n_iPI(Nx)
      real*8 n_rPHI(Nx)
      real*8 n_rPI(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = -0.1D1 * (n_rPI(i) * n_rPHI(i) + n_iPI(i) * n_iPHI(i)) / n_a1
     #(i)
      res(i)=qb
      end do
      END
