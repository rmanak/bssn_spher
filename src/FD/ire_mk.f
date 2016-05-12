      subroutine ire_mk(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Sx,n_psi,
     &x,Nx,hx,myzero,phys_bdy,res)
      implicit none
      integer i
      integer Nx
      real*8 hx
      real*8 myzero
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_Sx(Nx)
      real*8 n_psi(Nx)
      real*8 x(Nx)
      integer phys_bdy(2)
      real*8 res(Nx)
      real*8 qb
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      qb = -0.1D1 * n_psi(i) ** 6 / n_A(i) * n_Sx(i)
      res(i)=qb
      end do
      endif
      do i=2, Nx-1, 1
      qb = -0.1D1 * n_psi(i) ** 6 / n_A(i) * n_Sx(i) + (0.2D1 * n_psi(i)
     # ** 6 / n_A(i) ** 2 * n_Axx(i) - 0.2D1 * n_psi(i) ** 6 / n_A(i) / 
     #n_B(i) * n_Athth(i)) / ctfm(i) + (-0.5000000000000000D0 * n_psi(i)
     # ** 6 / n_A(i) ** 3 * n_Axx(i) * (-0.1D1 * n_A(i - 1) + n_A(i + 1)
     #) + 0.5000000000000000D0 * n_psi(i) ** 6 / n_A(i) ** 2 * (-0.1D1 *
     # n_Axx(i - 1) + n_Axx(i + 1)) + 0.5000000000000000D0 * n_psi(i) **
     # 6 / n_A(i) ** 2 / n_B(i) * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) * n
     #_Axx(i) - 0.5000000000000000D0 * n_psi(i) ** 6 / n_A(i) / n_B(i) *
     #* 2 * (-0.1D1 * n_B(i - 1) + n_B(i + 1)) * n_Athth(i) + 0.33333333
     #33333333D0 * n_psi(i) ** 6 / n_A(i) * (n_K(i - 1) - 0.1D1 * n_K(i 
     #+ 1)) + 0.3D1 * n_psi(i) ** 5 / n_A(i) ** 2 * (-0.1D1 * n_psi(i - 
     #1) + n_psi(i + 1)) * n_Axx(i)) / ctfmp(i) / hx
      res(i)=qb
      end do
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      qb = myzero * x(i)
      res(i)=qb
      end do
      endif
      END
