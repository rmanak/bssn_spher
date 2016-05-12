      subroutine compute_db(ctfm,ctfmp,n_A,n_B,n_beta,x,Nx,hx,myzero,phy
     &s_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_B(Nx)
      real*8 n_beta(Nx)
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
      qb = 0.2D1 / ctfm(i) * n_beta(i) + (-0.5000000000000000D0 * n_beta
     #(i - 1) + 0.5000000000000000D0 * n_beta(i + 1) + 0.250000000000000
     #0D0 / n_A(i) * n_beta(i) * (-0.1D1 * n_A(i - 1) + n_A(i + 1)) + 0.
     #5000000000000000D0 / n_B(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) *
     # n_beta(i)) / ctfmp(i) / hx
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
