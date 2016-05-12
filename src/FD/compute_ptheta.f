      subroutine compute_ptheta(ctfm,ctfmp,n_Athth,n_beta,n_divbeta,x,Nx
     &,hx,myzero,vee,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 vee
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_Athth(Nx)
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
      qb = 0.5000000000000000D0 / ctfmp(i) * (-0.1D1 * n_Athth(i - 1) + 
     #n_Athth(i + 1)) / hx * n_beta(i) + 0.2D1 / ctfm(i) * n_Athth(i) * 
     #n_beta(i) - 0.6666666666666667D0 * n_Athth(i) * vee * n_divbeta(i)
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
