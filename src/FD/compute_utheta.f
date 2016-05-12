      subroutine compute_utheta(n_Athth,n_B,n_K,x,Nx,myzero,phys_bdy,res
     &)
      implicit none
      integer i
      integer Nx
      real*8 myzero
      real*8 n_Athth(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = n_Athth(i) * n_K(i) - 0.2D1 * n_Athth(i) ** 2 / n_B(i)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
