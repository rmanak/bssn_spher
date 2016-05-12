      subroutine ire_eqipib2(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_iPIb2
     &,n_ipsi,nm1_iPIb2,np1_iPIb2,octfmp,Nx,ht,hx,mass,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 mass
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_b1(Nx)
      real*8 n_beta(Nx)
      real*8 n_iPIb2(Nx)
      real*8 n_ipsi(Nx)
      real*8 nm1_iPIb2(Nx)
      real*8 np1_iPIb2(Nx)
      real*8 octfmp(Nx)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_iPIb2(i) + np1_iPIb2(i))
     # / ht + n_a1(i) * mass * n_ipsi(i) * n_alpha(i) * n_b1(i) ** 2 - 0
     #.1D1 * (-0.5000000000000000D0 * n_iPIb2(i) * (n_beta(i - 1) - 0.1D
     #1 * n_beta(i + 1)) / hx - 0.5000000000000000D0 * n_beta(i) * (n_iP
     #Ib2(i - 1) - 0.1D1 * n_iPIb2(i + 1)) / hx + 0.5000000000000000D0 /
     # n_a1(i) * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) / hx * octfmp(i
     #) * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) - 0.1D1 * (-0.2500000000
     #000000D0 / n_a1(i) ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) 
     #/ hx ** 2 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) * n_alpha(i) * n_b
     #1(i) ** 2 - 0.2500000000000000D0 / n_a1(i) * (n_ipsi(i - 1) - 0.1D
     #1 * n_ipsi(i + 1)) / hx ** 2 * (-0.1D1 * n_alpha(i - 1) + n_alpha(
     #i + 1)) * n_b1(i) ** 2 - 0.5000000000000000D0 / n_a1(i) * (n_ipsi(
     #i - 1) - 0.1D1 * n_ipsi(i + 1)) / hx ** 2 * (-0.1D1 * n_b1(i - 1) 
     #+ n_b1(i + 1)) * n_alpha(i) * n_b1(i) - 0.1D1 / n_a1(i) * (-0.1D1 
     #* n_ipsi(i - 1) + 0.2D1 * n_ipsi(i) - 0.1D1 * n_ipsi(i + 1)) / hx 
     #** 2 * n_alpha(i) * n_b1(i) ** 2) / ctfmp(i) ** 2 - 0.1D1 * (0.2D1
     # * n_iPIb2(i) * n_beta(i) - 0.1D1 / n_a1(i) * n_alpha(i) * n_b1(i)
     # ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) / hx / ctfmp(i)) /
     # ctfm(i)
      res = res + qb**2
      end do
      res = sqrt(res/(1*Nx))
      END
