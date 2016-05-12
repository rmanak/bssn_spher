      subroutine resid_alpha_hgd_adv(ctfmp,n_K,n_alpha,n_beta,np1_K,np1_
     &alpha,np1_beta,x,Nx,advc,ht,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 advc
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 ctfmp(Nx)
      real*8 n_K(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 np1_K(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = (-0.1D1 * n_alpha(i) + np1_alpha(i)) / ht - 0.250000000000000
     #0D0 * advc * np1_beta(i) * (-0.1D1 * np1_alpha(i - 1) + np1_alpha(
     #i + 1)) / hx / ctfmp(i) + np1_alpha(i) * np1_K(i) - 0.250000000000
     #0000D0 * advc * n_beta(i) * (-0.1D1 * n_alpha(i - 1) + n_alpha(i +
     # 1)) / hx / ctfmp(i) + n_alpha(i) * n_K(i)
      res = res + qb**2
      end do
      endif
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_alpha(i) + np1_alpha(i)) / ht - 0.250000000000000
     #0D0 * advc * np1_beta(i) * (-0.1D1 * np1_alpha(i - 1) + np1_alpha(
     #i + 1)) / hx / ctfmp(i) + np1_alpha(i) * np1_K(i) - 0.250000000000
     #0000D0 * advc * n_beta(i) * (-0.1D1 * n_alpha(i - 1) + n_alpha(i +
     # 1)) / hx / ctfmp(i) + n_alpha(i) * n_K(i)
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_alpha(i) + myzero * x(i) - 0.1D1
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
