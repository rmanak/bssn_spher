      subroutine ire_dtk(ctfmp,n_A,n_Athth,n_Axx,n_B,n_DDLthetatheta,n_D
     &DLxx,n_K,n_TS,n_alpha,n_beta,n_metthetatheta,n_metxx,n_rho,nm1_K,n
     &p1_K,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_DDLthetatheta(Nx)
      real*8 n_DDLxx(Nx)
      real*8 n_K(Nx)
      real*8 n_TS(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_metthetatheta(Nx)
      real*8 n_metxx(Nx)
      real*8 n_rho(Nx)
      real*8 nm1_K(Nx)
      real*8 np1_K(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = -0.5000000000000000D0 * (nm1_K(i) - 0.1D1 * np1_K(i)) / ht - 
     #0.3333333333333333D0 * n_alpha(i) * n_K(i) ** 2 - 0.50000000000000
     #00D0 / ctfmp(i) * (-0.1D1 * n_K(i - 1) + n_K(i + 1)) / hx * n_beta
     #(i) - 0.5000000000000000D0 * n_alpha(i) * n_rho(i) - 0.50000000000
     #00000D0 * n_alpha(i) * n_TS(i) - 0.1D1 / n_A(i) ** 2 * n_Axx(i) **
     # 2 * n_alpha(i) - 0.2D1 / n_B(i) ** 2 * n_Athth(i) ** 2 * n_alpha(
     #i) + 0.2D1 / n_metthetatheta(i) * n_DDLthetatheta(i) + 0.1D1 / n_m
     #etxx(i) * n_DDLxx(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
