      subroutine ire_rpsi_direct(ctfm,ctfmp,n_a1,n_alpha,n_b1,n_beta,n_r
     &psi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_rpsi,np1_a1,np1_alpha,np1
     &_b1,np1_beta,np1_rpsi,octfmp,Nx,ht,hx,mass,res)
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
      real*8 n_rpsi(Nx)
      real*8 nm1_a1(Nx)
      real*8 nm1_alpha(Nx)
      real*8 nm1_b1(Nx)
      real*8 nm1_beta(Nx)
      real*8 nm1_rpsi(Nx)
      real*8 np1_a1(Nx)
      real*8 np1_alpha(Nx)
      real*8 np1_b1(Nx)
      real*8 np1_beta(Nx)
      real*8 np1_rpsi(Nx)
      real*8 octfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = -0.2500000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) 
     #/ ht ** 2 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / n_alpha(i) ** 
     #3 - 0.1D1 * mass * n_rpsi(i) + (-0.1D1 * nm1_rpsi(i) + 0.2D1 * n_r
     #psi(i) - 0.1D1 * np1_rpsi(i)) / ht ** 2 / n_alpha(i) ** 2 + 0.5000
     #000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / ht ** 2 * 
     #(nm1_b1(i) - 0.1D1 * np1_b1(i)) / n_alpha(i) ** 2 / n_b1(i) - 0.25
     #00000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / ht ** 2 
     #* (-0.1D1 * nm1_a1(i) + np1_a1(i)) / n_alpha(i) ** 2 / n_a1(i) + (
     #-0.2500000000000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / 
     #hx * n_beta(i) * (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht / n_alpha(i
     #) ** 2 / n_a1(i) - 0.2500000000000000D0 * (-0.1D1 * nm1_rpsi(i) + 
     #np1_rpsi(i)) / ht * n_beta(i) * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)
     #) / hx / n_alpha(i) ** 2 / n_a1(i) + 0.5000000000000000D0 * (n_rps
     #i(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx * n_beta(i) * (nm1_b1(i) - 
     #0.1D1 * np1_b1(i)) / ht / n_alpha(i) ** 2 / n_b1(i) - 0.5000000000
     #000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / ht * n_beta(i) *
     # (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) / hx / n_alpha(i) ** 2 / n_b1
     #(i) - 0.5000000000000000D0 * octfmp(i) * (n_rpsi(i - 1) - 0.1D1 * 
     #n_rpsi(i + 1)) / hx * n_beta(i) ** 2 / n_alpha(i) ** 2 - 0.2500000
     #000000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx * n_bet
     #a(i) * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / ht / n_alpha(i) ** 
     #3 + 0.2500000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / 
     #ht * n_beta(i) * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx / 
     #n_alpha(i) ** 3 + 0.2500000000000000D0 * (n_rpsi(i - 1) - 0.1D1 * 
     #n_rpsi(i + 1)) / hx * (nm1_beta(i) - 0.1D1 * np1_beta(i)) / ht / n
     #_alpha(i) ** 2 - 0.2500000000000000D0 * (-0.1D1 * nm1_rpsi(i) + np
     #1_rpsi(i)) / ht * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx / n
     #_alpha(i) ** 2 + 0.5000000000000000D0 * octfmp(i) * (n_rpsi(i - 1)
     # - 0.1D1 * n_rpsi(i + 1)) / hx / n_a1(i) ** 2 - 0.5000000000000000
     #D0 * n_beta(i) * (-0.1D1 * nm1_rpsi(i - 1) + nm1_rpsi(i + 1) + np1
     #_rpsi(i - 1) - 0.1D1 * np1_rpsi(i + 1)) / hx / ht / n_alpha(i) ** 
     #2) / ctfmp(i) + ((n_rpsi(i - 1) - 0.2D1 * n_rpsi(i) + n_rpsi(i + 1
     #)) / hx ** 2 / n_a1(i) ** 2 - 0.2500000000000000D0 * (n_rpsi(i - 1
     #) - 0.1D1 * n_rpsi(i + 1)) / hx ** 2 * n_beta(i) ** 2 * (n_a1(i - 
     #1) - 0.1D1 * n_a1(i + 1)) / n_alpha(i) ** 2 / n_a1(i) - 0.50000000
     #00000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx ** 2 * n
     #_beta(i) ** 2 * (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) / n_alpha(i) *
     #* 2 / n_b1(i) + 0.2500000000000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_
     #rpsi(i + 1)) / hx ** 2 * n_beta(i) ** 2 * (n_alpha(i - 1) - 0.1D1 
     #* n_alpha(i + 1)) / n_alpha(i) ** 3 + 0.2500000000000000D0 * (n_rp
     #si(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx ** 2 * (n_alpha(i - 1) - 0
     #.1D1 * n_alpha(i + 1)) / n_alpha(i) / n_a1(i) ** 2 + 0.50000000000
     #00000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx ** 2 * (n_b
     #1(i - 1) - 0.1D1 * n_b1(i + 1)) / n_b1(i) / n_a1(i) ** 2 - 0.50000
     #00000000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx ** 2 
     #* n_beta(i) * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / n_alpha(i)
     # ** 2 - 0.2500000000000000D0 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i +
     # 1)) / hx ** 2 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / n_a1(i) ** 
     #3 - 0.1D1 * n_beta(i) ** 2 * (n_rpsi(i - 1) - 0.2D1 * n_rpsi(i) + 
     #n_rpsi(i + 1)) / hx ** 2 / n_alpha(i) ** 2) / ctfmp(i) ** 2 + ((-0
     #.1D1 * nm1_rpsi(i) + np1_rpsi(i)) / ht * n_beta(i) / n_alpha(i) **
     # 2 + ((n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) / hx * n_beta(i) ** 
     #2 / n_alpha(i) ** 2 - 0.1D1 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 
     #1)) / hx / n_a1(i) ** 2) / ctfmp(i)) / ctfm(i)
      res(i)=qb
      end do
      END
