      subroutine compute_ipi(n_b1,n_iPIb2,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_b1(Nx)
      real*8 n_iPIb2(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_iPIb2(i) / n_b1(i) ** 2
      res(i)=qb
      end do
      END
