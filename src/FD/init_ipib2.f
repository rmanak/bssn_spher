      subroutine init_ipib2(n_b1,n_iPI,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_b1(Nx)
      real*8 n_iPI(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_iPI(i) * n_b1(i) ** 2
      res(i)=qb
      end do
      END
