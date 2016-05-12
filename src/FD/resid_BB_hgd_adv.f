      subroutine resid_BB_hgd_adv(ctfmp,n_BB,n_DLamx,n_beta,np1_BB,np1_D
     &Lamx,np1_beta,x,Nx,advc,eta,ht,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 advc
      real*8 eta
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 ctfmp(Nx)
      real*8 n_BB(Nx)
      real*8 n_DLamx(Nx)
      real*8 n_beta(Nx)
      real*8 np1_BB(Nx)
      real*8 np1_DLamx(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = np1_BB(i) + myzero * x(i)
      res = res + qb**2
      end do
      endif
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_BB(i) + np1_BB(i)) / ht - 0.2500000000000000D0 * 
     #advc * np1_beta(i) * (-0.1D1 * np1_BB(i - 1) + np1_BB(i + 1)) / hx
     # / ctfmp(i) + 0.5000000000000000D0 * eta * np1_BB(i) - 0.500000000
     #0000000D0 * np1_DLamx(i) - 0.2500000000000000D0 * advc * n_beta(i)
     # * (-0.1D1 * n_BB(i - 1) + n_BB(i + 1)) / hx / ctfmp(i) + 0.500000
     #0000000000D0 * eta * n_BB(i) - 0.5000000000000000D0 * n_DLamx(i)
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_BB(i) + myzero * x(i)
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
