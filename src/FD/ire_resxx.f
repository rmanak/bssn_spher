      subroutine ire_resxx(ctfm,ctfmp,n_U,n_a1,n_alpha,n_b1,n_beta,n_ips
     &i,n_rpsi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_ipsi,nm1_rpsi,np1_a1
     &,np1_alpha,np1_b1,np1_beta,np1_ipsi,np1_rpsi,Nx,ht,hx,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_U(Nx)
      real*8 n_a1(Nx)
      real*8 n_alpha(Nx)
      real*8 n_b1(Nx)
      real*8 n_beta(Nx)
      real*8 n_ipsi(Nx)
      real*8 n_rpsi(Nx)
      real*8 nm1_a1(Nx)
      real*8 nm1_alpha(Nx)
      real*8 nm1_b1(Nx)
      real*8 nm1_beta(Nx)
      real*8 nm1_ipsi(Nx)
      real*8 nm1_rpsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_ipsi(Nx)
      real*8 np1_rpsi(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = -0.2500000000000000D0 / n_alpha(i) ** 2 / n_b1(i) ** 2 * (-0.
     #1D1 * nm1_b1(i) + np1_b1(i)) ** 2 / ht ** 2 - 0.1250000000000000D0
     # / n_alpha(i) ** 2 * (-0.1D1 * nm1_ipsi(i) + np1_ipsi(i)) ** 2 / h
     #t ** 2 - 0.1250000000000000D0 / n_alpha(i) ** 2 * (-0.1D1 * nm1_rp
     #si(i) + np1_rpsi(i)) ** 2 / ht ** 2 - 0.5000000000000000D0 * (nm1_
     #alpha(i) - 0.1D1 * np1_alpha(i)) / ht ** 2 / n_b1(i) / n_alpha(i) 
     #** 3 * (-0.1D1 * nm1_b1(i) + np1_b1(i)) + 0.2D1 / n_alpha(i) ** 2 
     #/ n_b1(i) * (-0.1D1 * nm1_b1(i) + 0.2D1 * n_b1(i) - 0.1D1 * np1_b1
     #(i)) / ht ** 2 + 0.5000000000000000D0 * n_U(i) + (0.50000000000000
     #00D0 * (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht / n_alpha(i) ** 2 / n
     #_b1(i) / n_a1(i) * n_beta(i) * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)
     #) / hx + 0.5000000000000000D0 / n_b1(i) / n_alpha(i) ** 2 * n_beta
     #(i) * (nm1_b1(i - 1) - 0.1D1 * nm1_b1(i + 1) - 0.1D1 * np1_b1(i - 
     #1) + np1_b1(i + 1)) / hx / ht + 0.5000000000000000D0 * (-0.1D1 * n
     #_b1(i - 1) + n_b1(i + 1)) / hx / n_alpha(i) ** 2 / n_b1(i) * (-0.1
     #D1 * nm1_beta(i) + np1_beta(i)) / ht + 0.5000000000000000D0 * (nm1
     #_alpha(i) - 0.1D1 * np1_alpha(i)) / ht / n_b1(i) * n_beta(i) / n_a
     #lpha(i) ** 3 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) / hx + 0.50000
     #00000000000D0 * (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht / n_alpha(i)
     # ** 2 * n_beta(i) / n_b1(i) ** 2 * (-0.1D1 * n_b1(i - 1) + n_b1(i 
     #+ 1)) / hx) / ctfmp(i) + (0.5000000000000000D0 * (n_a1(i - 1) - 0.
     #1D1 * n_a1(i + 1)) / hx ** 2 / n_alpha(i) ** 2 / n_b1(i) / n_a1(i)
     # * n_beta(i) ** 2 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) + 0.25000
     #00000000000D0 / n_b1(i) ** 2 / n_a1(i) ** 2 * (-0.1D1 * n_b1(i - 1
     #) + n_b1(i + 1)) ** 2 / hx ** 2 + 0.1250000000000000D0 / n_alpha(i
     #) ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) ** 2 / hx ** 2 * 
     #n_beta(i) ** 2 + 0.1250000000000000D0 / n_alpha(i) ** 2 * (-0.1D1 
     #* n_rpsi(i - 1) + n_rpsi(i + 1)) ** 2 / hx ** 2 * n_beta(i) ** 2 -
     # 0.1250000000000000D0 / n_a1(i) ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_
     #ipsi(i + 1)) ** 2 / hx ** 2 - 0.1250000000000000D0 / n_a1(i) ** 2 
     #* (-0.1D1 * n_rpsi(i - 1) + n_rpsi(i + 1)) ** 2 / hx ** 2 - 0.2500
     #000000000000D0 / n_alpha(i) ** 2 / n_b1(i) ** 2 * (-0.1D1 * n_b1(i
     # - 1) + n_b1(i + 1)) ** 2 / hx ** 2 * n_beta(i) ** 2 - 0.500000000
     #0000000D0 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx ** 2 / 
     #n_alpha(i) / n_b1(i) / n_a1(i) ** 2 * (-0.1D1 * n_b1(i - 1) + n_b1
     #(i + 1)) + 0.5000000000000000D0 * (-0.1D1 * n_b1(i - 1) + n_b1(i +
     # 1)) / hx ** 2 / n_alpha(i) ** 2 / n_b1(i) * n_beta(i) * (n_beta(i
     # - 1) - 0.1D1 * n_beta(i + 1))) / ctfmp(i) ** 2 + ((nm1_alpha(i) -
     # 0.1D1 * np1_alpha(i)) / ht / n_alpha(i) ** 3 * n_beta(i) + 0.1D1 
     #/ n_alpha(i) ** 2 * (-0.1D1 * nm1_beta(i) + np1_beta(i)) / ht + 0.
     #2D1 * (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht / n_alpha(i) ** 2 / n_
     #b1(i) * n_beta(i) + (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht / n_alph
     #a(i) ** 2 / n_a1(i) * n_beta(i) + ((-0.1D1 * n_b1(i - 1) + n_b1(i 
     #+ 1)) / hx / n_b1(i) / n_a1(i) ** 2 + 0.1D1 / n_alpha(i) ** 2 * n_
     #beta(i) * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx - 0.1D1 * (
     #n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx / n_alpha(i) / n_a1(
     #i) ** 2 + (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_a1(i) * n_b
     #eta(i) ** 2 / n_alpha(i) ** 2 - 0.1D1 * (-0.1D1 * n_b1(i - 1) + n_
     #b1(i + 1)) / hx / n_alpha(i) ** 2 / n_b1(i) * n_beta(i) ** 2) / ct
     #fmp(i)) / ctfm(i) + (-0.1D1 * n_beta(i) ** 2 / n_alpha(i) ** 2 + 0
     #.1D1 / n_a1(i) ** 2 - 0.1D1 / n_b1(i) ** 2) / ctfm(i) ** 2
      res(i)=qb
      end do
      END
