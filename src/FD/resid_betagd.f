      subroutine resid_betagd(n_Lamx,n_beta,np1_Lamx,np1_beta,x,Nx,eta,h
     &t,mus,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 eta
      real*8 ht
      real*8 mus
      real*8 myzero
      real*8 n_Lamx(Nx)
      real*8 n_beta(Nx)
      real*8 np1_Lamx(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = (-0.1D1 * n_beta(i) + np1_beta(i)) / ht - 0.5000000000000000D
     #0 * mus * np1_Lamx(i) + 0.5000000000000000D0 * eta * np1_beta(i) -
     # 0.5000000000000000D0 * mus * n_Lamx(i) + 0.5000000000000000D0 * e
     #ta * n_beta(i)
      res = res + qb**2
      end do
      endif
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_beta(i) + np1_beta(i)) / ht - 0.5000000000000000D
     #0 * mus * np1_Lamx(i) + 0.5000000000000000D0 * eta * np1_beta(i) -
     # 0.5000000000000000D0 * mus * n_Lamx(i) + 0.5000000000000000D0 * e
     #ta * n_beta(i)
      res = res + qb**2
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = np1_beta(i) + myzero * x(i)
      res = res + qb**2
      end do
      endif
      res = sqrt(res/(1*Nx))
      END
