      subroutine compute_drpsi(ctfm,ctfmp,n_rpsi,x,Nx,hx,myzero,phys_bdy
     &,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_rpsi(Nx)
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
      qb = myzero * ctfm(i) - 0.5000000000000000D0 / ctfmp(i) * (n_rpsi(
     #i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
