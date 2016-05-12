      subroutine compute_emph(n_phi,x,Nx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 myzero
      real*8 n_phi(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = myzero * x(i) + 0.1D1
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = 0.1D1 / exp(n_phi(i)) ** 4
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i) + 0.1D1
      res(i)=qb
      end do
      endif
      END
