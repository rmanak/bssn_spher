      subroutine compute_rpi(n_b1,n_rPIb2,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_b1(Nx)
      real*8 n_rPIb2(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_rPIb2(i) / n_b1(i) ** 2
      res(i)=qb
      end do
      END
