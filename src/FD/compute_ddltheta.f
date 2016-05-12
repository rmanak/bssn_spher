      subroutine compute_ddltheta(ctfm,ctfmp,n_alpha,n_metthetatheta,n_m
     &etxx,x,Nx,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_alpha(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
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
      qb = -0.1250000000000000D0 / ctfmp(i) ** 2 * (n_alpha(i - 1) - 0.1
     #D1 * n_alpha(i + 1)) / hx ** 2 / n_metxx(i) * (-0.1D1 * n_mettheta
     #theta(i - 1) + n_metthetatheta(i + 1)) - 0.5000000000000000D0 * (n
     #_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx / ctfmp(i) / ctfm(i) 
     #/ n_metxx(i) * n_metthetatheta(i)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
