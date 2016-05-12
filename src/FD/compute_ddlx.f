      subroutine compute_ddlx(ctfmp,n_alpha,n_metxx,octfmp,x,Nx,hx,myzer
     &o,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfmp(Nx)
      real*8 n_alpha(Nx)
      real*8 n_metxx(Nx)
      real*8 octfmp(Nx)
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
      qb = 0.5000000000000000D0 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 
     #1)) / hx / ctfmp(i) * octfmp(i) + (0.1250000000000000D0 / n_metxx(
     #i) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) * (-0.1D1 * n_metxx
     #(i - 1) + n_metxx(i + 1)) + n_alpha(i - 1) - 0.2D1 * n_alpha(i) + 
     #n_alpha(i + 1)) / ctfmp(i) ** 2 / hx ** 2
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
