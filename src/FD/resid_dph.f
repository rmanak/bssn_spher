      subroutine resid_dph(ctfmp,n_K,n_alpha,n_beta,n_divbeta,n_phi,np1_
     &K,np1_alpha,np1_beta,np1_divbeta,np1_phi,x,Nx,ht,hx,myzero,zepsdis
     &,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 zepsdis
      real*8 ctfmp(Nx)
      real*8 n_K(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 n_phi(Nx)
      real*8 np1_K(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_divbeta(Nx)
      real*8 np1_phi(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = -0.5000000000000000D0 * (0.3D1 * np1_phi(i) - 0.4D1 * np1_phi
     #(i + 1) + np1_phi(i + 2)) / hx
      res = res + qb**2
      end do
      endif
      if (phys_bdy(1) .eq. 1) then
      do i=2, 2, 1
      qb = (-0.1D1 * n_phi(i) + np1_phi(i)) / ht + 0.2500000000000000D0 
     #/ ctfmp(i) * (np1_phi(i - 1) - 0.1D1 * np1_phi(i + 1)) / hx * np1_
     #beta(i) + 0.8333333333333333D-1 * np1_alpha(i) * np1_K(i) - 0.8333
     #333333333333D-1 * np1_divbeta(i) + 0.2500000000000000D0 / ctfmp(i)
     # * (n_phi(i - 1) - 0.1D1 * n_phi(i + 1)) / hx * n_beta(i) + 0.8333
     #333333333333D-1 * n_alpha(i) * n_K(i) - 0.8333333333333333D-1 * n_
     #divbeta(i) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * np1_p
     #hi(i) + np1_phi(i + 2) - 0.4D1 * np1_phi(i + 1) - 0.4D1 * np1_phi(
     #i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * (0.7D1 * n_phi(i)
     # + n_phi(i + 2) - 0.4D1 * n_phi(i + 1) - 0.4D1 * n_phi(i - 1))
      res = res + qb**2
      end do
      endif
      do i=3, Nx-2, 1
      qb = (-0.1D1 * n_phi(i) + np1_phi(i)) / ht + 0.2500000000000000D0 
     #/ ctfmp(i) * (np1_phi(i - 1) - 0.1D1 * np1_phi(i + 1)) / hx * np1_
     #beta(i) + 0.8333333333333333D-1 * np1_alpha(i) * np1_K(i) - 0.8333
     #333333333333D-1 * np1_divbeta(i) + 0.2500000000000000D0 / ctfmp(i)
     # * (n_phi(i - 1) - 0.1D1 * n_phi(i + 1)) / hx * n_beta(i) + 0.8333
     #333333333333D-1 * n_alpha(i) * n_K(i) - 0.8333333333333333D-1 * n_
     #divbeta(i) + 0.3125000000000000D-1 * zepsdis / ht * (0.6D1 * np1_p
     #hi(i) + np1_phi(i + 2) + np1_phi(i - 2) - 0.4D1 * np1_phi(i + 1) -
     # 0.4D1 * np1_phi(i - 1)) + 0.3125000000000000D-1 * zepsdis / ht * 
     #(0.6D1 * n_phi(i) + n_phi(i + 2) + n_phi(i - 2) - 0.4D1 * n_phi(i 
     #+ 1) - 0.4D1 * n_phi(i - 1))
      res = res + qb**2
      end do
      do i=Nx-1, Nx-1, 1
      qb = (-0.1D1 * n_phi(i) + np1_phi(i)) / ht + 0.2500000000000000D0 
     #/ ctfmp(i) * (np1_phi(i - 1) - 0.1D1 * np1_phi(i + 1)) / hx * np1_
     #beta(i) + 0.8333333333333333D-1 * np1_alpha(i) * np1_K(i) - 0.8333
     #333333333333D-1 * np1_divbeta(i) + 0.2500000000000000D0 / ctfmp(i)
     # * (n_phi(i - 1) - 0.1D1 * n_phi(i + 1)) / hx * n_beta(i) + 0.8333
     #333333333333D-1 * n_alpha(i) * n_K(i) - 0.8333333333333333D-1 * n_
     #divbeta(i)
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = (-0.1D1 * n_phi(i) + np1_phi(i)) / ht - 0.1D1 * myzero * x(i)
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
