      subroutine resid_alphadyopl(n_K,n_alpha,np1_K,np1_alpha,x,Nx,ht,my
     &zero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 myzero
      real*8 n_K(Nx)
      real*8 n_alpha(Nx)
      real*8 np1_K(Nx)
      real*8 np1_alpha(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res
      real*8 qb
      res = 0.0D0
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = (-0.1D1 * n_alpha(i) + np1_alpha(i)) / ht + np1_alpha(i) * np
     #1_K(i) + n_alpha(i) * n_K(i)
      res = res + qb**2
      end do
      endif
      do i=2, Nx-1, 1
      qb = (-0.1D1 * n_alpha(i) + np1_alpha(i)) / ht + np1_alpha(i) * np
     #1_K(i) + n_alpha(i) * n_K(i)
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
