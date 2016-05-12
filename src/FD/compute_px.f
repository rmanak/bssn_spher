      subroutine compute_px(ctfmp,n_Axx,n_beta,n_divbeta,x,Nx,hx,myzero,
     &vee,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 vee
      real*8 ctfmp(Nx)
      real*8 n_Axx(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
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
      qb = -0.6666666666666667D0 * n_Axx(i) * vee * n_divbeta(i) + (n_Ax
     #x(i) * (-0.1D1 * n_beta(i - 1) + n_beta(i + 1)) + 0.50000000000000
     #00D0 * (-0.1D1 * n_Axx(i - 1) + n_Axx(i + 1)) * n_beta(i)) / ctfmp
     #(i) / hx
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
