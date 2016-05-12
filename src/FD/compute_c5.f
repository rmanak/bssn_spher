      subroutine compute_c5(n_A,n_Athth,n_Axx,n_B,n_K,n_rho,x,Nx,myzero,
     &phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 myzero
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_rho(Nx)
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
      qb = -0.8333333333333333D-1 * n_K(i) ** 2 + 0.2500000000000000D0 *
     # n_rho(i) + 0.1250000000000000D0 / n_A(i) ** 2 * n_Axx(i) ** 2 + 0
     #.2500000000000000D0 / n_B(i) ** 2 * n_Athth(i) ** 2
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
