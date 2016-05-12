      subroutine compute_phi2(n_iPHI,n_rPHI,Nx,res)
      implicit none
      integer i
      integer Nx
      real*8 n_iPHI(Nx)
      real*8 n_rPHI(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=1, Nx, 1
      qb = n_rPHI(i) ** 2 + n_iPHI(i) ** 2
      res(i)=qb
      end do
      END
