      subroutine ire_restx(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_ipsi,n_
     &rpsi,nm1_a1,nm1_b1,nm1_ipsi,nm1_rpsi,np1_a1,np1_b1,np1_ipsi,np1_rp
     &si,octfmp,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_b1(Nx)
      real*8 n_beta(Nx)
      real*8 n_ipsi(Nx)
      real*8 n_rpsi(Nx)
      real*8 nm1_a1(Nx)
      real*8 nm1_b1(Nx)
      real*8 nm1_ipsi(Nx)
      real*8 nm1_rpsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_ipsi(Nx)
      real*8 np1_rpsi(Nx)
      real*8 octfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = (-0.5000000000000000D0 / n_alpha(i) ** 2 / n_b1(i) / n_a1(i) 
     #* (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht * (-0.1D1 * n_b1(i - 1) + 
     #n_b1(i + 1)) / hx + 0.5000000000000000D0 / n_alpha(i) ** 3 / n_b1(
     #i) * (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht * (n_alpha(i - 1) - 0.1
     #D1 * n_alpha(i + 1)) / hx + (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) /
     # hx / n_b1(i) / n_alpha(i) ** 2 * n_beta(i) * octfmp(i) + 0.250000
     #0000000000D0 / n_alpha(i) ** 2 * (-0.1D1 * n_rpsi(i - 1) + n_rpsi(
     #i + 1)) / hx * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / ht - 0.25000
     #00000000000D0 / n_alpha(i) ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(
     #i + 1)) / hx * (-0.1D1 * nm1_ipsi(i) + np1_ipsi(i)) / ht + 0.50000
     #00000000000D0 / n_alpha(i) ** 2 / n_b1(i) * (nm1_b1(i - 1) - 0.1D1
     # * nm1_b1(i + 1) - 0.1D1 * np1_b1(i - 1) + np1_b1(i + 1)) / hx / h
     #t) / ctfmp(i) + (-0.5000000000000000D0 / n_alpha(i) ** 2 / n_b1(i)
     # / n_a1(i) * n_beta(i) * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx 
     #** 2 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) - 0.5000000000000000D0
     # / n_alpha(i) ** 3 / n_b1(i) * n_beta(i) * (n_alpha(i - 1) - 0.1D1
     # * n_alpha(i + 1)) / hx ** 2 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)
     #) - 0.2500000000000000D0 / n_alpha(i) ** 2 * n_beta(i) * (-0.1D1 *
     # n_rpsi(i - 1) + n_rpsi(i + 1)) ** 2 / hx ** 2 - 0.250000000000000
     #0D0 / n_alpha(i) ** 2 * n_beta(i) * (n_ipsi(i - 1) - 0.1D1 * n_ips
     #i(i + 1)) ** 2 / hx ** 2 + 0.2D1 / n_alpha(i) ** 2 / n_b1(i) * n_b
     #eta(i) * (-0.1D1 * n_b1(i - 1) + 0.2D1 * n_b1(i) - 0.1D1 * n_b1(i 
     #+ 1)) / hx ** 2) / ctfmp(i) ** 2 + (-0.1D1 * (-0.1D1 * nm1_a1(i) +
     # np1_a1(i)) / ht / n_alpha(i) ** 2 / n_a1(i) + 0.1D1 / n_alpha(i) 
     #** 2 / n_b1(i) * (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht + (-0.1D1 *
     # (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_alpha(i) ** 2 / n_a1
     #(i) * n_beta(i) - 0.1D1 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)
     #) / hx / n_alpha(i) ** 3 * n_beta(i) - 0.2D1 * (-0.1D1 * n_b1(i - 
     #1) + n_b1(i + 1)) / hx / n_alpha(i) ** 2 / n_b1(i) * n_beta(i)) / 
     #ctfmp(i)) / ctfm(i)
      res(i)=qb
      end do
      END
