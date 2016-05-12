      subroutine compute_pi2(n_iPI,n_rPI,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_iPI(Nx)
      real*8 n_rPI(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_rPI(i) ** 2 + n_iPI(i) ** 2
      res(i)=qb
      end do
      END
