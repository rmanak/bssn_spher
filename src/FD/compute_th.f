      subroutine compute_th(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,np1_a1,n
     &p1_alpha,np1_b1,np1_beta,x,Nx,ht,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_b1(Nx)
      real*8 n_beta(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_beta(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = 0.1D1 + myzero * x(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = ((0.5000000000000000D0 / np1_alpha(i) + 0.5000000000000000D0 
     #/ n_alpha(i)) * (-0.1D1 * n_b1(i) + np1_b1(i)) / ht + 0.5000000000
     #000000D0 * (0.1D1 / np1_a1(i) - 0.1D1 * np1_beta(i) / np1_alpha(i)
     #) * (-0.5000000000000000D0 * (np1_b1(i - 1) - 0.1D1 * np1_b1(i + 1
     #)) / hx / ctfmp(i) + np1_b1(i) / ctfm(i)) + 0.5000000000000000D0 *
     # (0.1D1 / n_a1(i) - 0.1D1 * n_beta(i) / n_alpha(i)) * (-0.50000000
     #00000000D0 * (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) / hx / ctfmp(i) +
     # n_b1(i) / ctfm(i))) / (0.5000000000000000D0 * np1_b1(i) + 0.50000
     #00000000000D0 * n_b1(i))
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
