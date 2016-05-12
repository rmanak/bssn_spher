      subroutine ire_resthth(ctfm,ctfmp,n_U,n_a1,n_alpha,n_b1,n_beta,n_i
     &psi,n_rpsi,nm1_a1,nm1_alpha,nm1_b1,nm1_beta,nm1_ipsi,nm1_rpsi,np1_
     &a1,np1_alpha,np1_b1,np1_beta,np1_ipsi,np1_rpsi,octfmp,Nx,ht,hx,res
     &)
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
      real*8 octfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = -0.1250000000000000D0 / n_alpha(i) ** 2 * (-0.1D1 * nm1_ipsi(
     #i) + np1_ipsi(i)) ** 2 / ht ** 2 - 0.1250000000000000D0 / n_alpha(
     #i) ** 2 * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i)) ** 2 / ht ** 2 - 0.
     #2500000000000000D0 * (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht ** 2 / 
     #n_a1(i) / n_alpha(i) ** 3 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) 
     #- 0.2500000000000000D0 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / h
     #t ** 2 / n_b1(i) / n_alpha(i) ** 3 * (-0.1D1 * nm1_b1(i) + np1_b1(
     #i)) - 0.2500000000000000D0 / n_alpha(i) ** 2 / n_b1(i) / n_a1(i) *
     # (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht ** 2 * (-0.1D1 * nm1_a1(i) 
     #+ np1_a1(i)) + 0.1D1 / n_alpha(i) ** 2 / n_a1(i) * (-0.1D1 * nm1_a
     #1(i) + 0.2D1 * n_a1(i) - 0.1D1 * np1_a1(i)) / ht ** 2 + 0.1D1 / n_
     #alpha(i) ** 2 / n_b1(i) * (-0.1D1 * nm1_b1(i) + 0.2D1 * n_b1(i) - 
     #0.1D1 * np1_b1(i)) / ht ** 2 + 0.5000000000000000D0 * n_U(i) + (-0
     #.2500000000000000D0 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / ht /
     # n_alpha(i) ** 3 * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx + 
     #0.5000000000000000D0 / n_b1(i) / n_alpha(i) ** 2 * n_beta(i) * (nm
     #1_b1(i - 1) - 0.1D1 * nm1_b1(i + 1) - 0.1D1 * np1_b1(i - 1) + np1_
     #b1(i + 1)) / hx / ht + 0.2500000000000000D0 * (-0.1D1 * nm1_a1(i) 
     #+ np1_a1(i)) / ht / n_a1(i) * n_beta(i) / n_alpha(i) ** 3 * (n_alp
     #ha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx - 0.2500000000000000D0 * 
     #(n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_a1(i) * n_beta(i) / n
     #_alpha(i) ** 3 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / ht - 0.50
     #00000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_al
     #pha(i) ** 2 / n_a1(i) * n_beta(i) ** 2 * octfmp(i) - 0.25000000000
     #00000D0 * (-0.1D1 * nm1_b1(i) + np1_b1(i)) / ht / n_alpha(i) ** 2 
     #/ n_b1(i) * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx - 0.50000
     #00000000000D0 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) / hx / n_b1(i
     #) / n_a1(i) ** 2 * octfmp(i) + 0.2500000000000000D0 / n_alpha(i) *
     #* 2 * (nm1_beta(i - 1) - 0.1D1 * nm1_beta(i + 1) - 0.1D1 * np1_bet
     #a(i - 1) + np1_beta(i + 1)) / hx / ht - 0.2500000000000000D0 * (n_
     #a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_alpha(i) ** 2 / n_a1(i) 
     #* (-0.1D1 * nm1_beta(i) + np1_beta(i)) / ht - 0.5000000000000000D0
     # * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx / n_alpha(i) ** 2 
     #* n_beta(i) * octfmp(i) - 0.5000000000000000D0 * (-0.1D1 * nm1_a1(
     #i) + np1_a1(i)) / ht / n_alpha(i) ** 2 / n_a1(i) * (n_beta(i - 1) 
     #- 0.1D1 * n_beta(i + 1)) / hx + 0.5000000000000000D0 * (n_alpha(i 
     #- 1) - 0.1D1 * n_alpha(i + 1)) / hx / n_alpha(i) / n_a1(i) ** 2 * 
     #octfmp(i) + 0.2500000000000000D0 * (-0.1D1 * n_b1(i - 1) + n_b1(i 
     #+ 1)) / hx / n_alpha(i) ** 2 / n_b1(i) * (-0.1D1 * nm1_beta(i) + n
     #p1_beta(i)) / ht + 0.5000000000000000D0 * (-0.1D1 * n_b1(i - 1) + 
     #n_b1(i + 1)) / hx / n_alpha(i) ** 2 / n_b1(i) * octfmp(i) * n_beta
     #(i) ** 2 + 0.2500000000000000D0 * (n_alpha(i - 1) - 0.1D1 * n_alph
     #a(i + 1)) / hx / n_b1(i) * n_beta(i) / n_alpha(i) ** 3 * (-0.1D1 *
     # nm1_b1(i) + np1_b1(i)) / ht - 0.2500000000000000D0 / n_alpha(i) *
     #* 2 * n_beta(i) * (-0.1D1 * nm1_ipsi(i) + np1_ipsi(i)) / ht * (n_i
     #psi(i - 1) - 0.1D1 * n_ipsi(i + 1)) / hx + 0.2500000000000000D0 / 
     #n_alpha(i) ** 2 * n_beta(i) * (-0.1D1 * nm1_rpsi(i) + np1_rpsi(i))
     # / ht * (-0.1D1 * n_rpsi(i - 1) + n_rpsi(i + 1)) / hx - 0.50000000
     #00000000D0 / n_a1(i) / n_alpha(i) ** 2 * n_beta(i) * (-0.1D1 * nm1
     #_a1(i - 1) + nm1_a1(i + 1) + np1_a1(i - 1) - 0.1D1 * np1_a1(i + 1)
     #) / hx / ht + 0.2500000000000000D0 * (nm1_alpha(i) - 0.1D1 * np1_a
     #lpha(i)) / ht / n_b1(i) * n_beta(i) / n_alpha(i) ** 3 * (-0.1D1 * 
     #n_b1(i - 1) + n_b1(i + 1)) / hx + 0.2500000000000000D0 * (-0.1D1 *
     # nm1_a1(i) + np1_a1(i)) / ht / n_alpha(i) ** 2 / n_b1(i) / n_a1(i)
     # * n_beta(i) * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) / hx - 0.25000
     #00000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_alpha
     #(i) ** 2 / n_b1(i) / n_a1(i) * n_beta(i) * (-0.1D1 * nm1_b1(i) + n
     #p1_b1(i)) / ht) / ctfmp(i) + (-0.2500000000000000D0 / n_alpha(i) *
     #* 2 * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) ** 2 / hx ** 2 + 0.2
     #500000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx ** 2 
     #/ n_b1(i) / n_a1(i) ** 3 * (-0.1D1 * n_b1(i - 1) + n_b1(i + 1)) + 
     #0.2500000000000000D0 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) /
     # hx ** 2 * n_beta(i) / n_alpha(i) ** 3 * (n_beta(i - 1) - 0.1D1 * 
     #n_beta(i + 1)) - 0.1D1 / n_alpha(i) ** 2 / n_a1(i) * n_beta(i) ** 
     #2 * (n_a1(i - 1) - 0.2D1 * n_a1(i) + n_a1(i + 1)) / hx ** 2 - 0.25
     #00000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx ** 2 /
     # n_alpha(i) / n_a1(i) ** 3 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i +
     # 1)) - 0.2500000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) 
     #/ hx ** 2 / n_alpha(i) ** 2 / n_b1(i) / n_a1(i) * n_beta(i) ** 2 *
     # (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) - 0.1250000000000000D0 / n_al
     #pha(i) ** 2 * (n_ipsi(i - 1) - 0.1D1 * n_ipsi(i + 1)) ** 2 / hx **
     # 2 * n_beta(i) ** 2 - 0.1250000000000000D0 / n_alpha(i) ** 2 * (n_
     #rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) ** 2 / hx ** 2 * n_beta(i) **
     # 2 - 0.1D1 / n_b1(i) / n_a1(i) ** 2 * (-0.1D1 * n_b1(i - 1) + 0.2D
     #1 * n_b1(i) - 0.1D1 * n_b1(i + 1)) / hx ** 2 + 0.1D1 / n_alpha(i) 
     #/ n_a1(i) ** 2 * (n_alpha(i - 1) - 0.2D1 * n_alpha(i) + n_alpha(i 
     #+ 1)) / hx ** 2 + 0.1D1 / n_alpha(i) ** 2 * n_beta(i) * (-0.1D1 * 
     #n_beta(i - 1) + 0.2D1 * n_beta(i) - 0.1D1 * n_beta(i + 1)) / hx **
     # 2 + 0.1250000000000000D0 / n_a1(i) ** 2 * (n_ipsi(i - 1) - 0.1D1 
     #* n_ipsi(i + 1)) ** 2 / hx ** 2 + 0.1250000000000000D0 / n_a1(i) *
     #* 2 * (n_rpsi(i - 1) - 0.1D1 * n_rpsi(i + 1)) ** 2 / hx ** 2 + 0.1
     #D1 / n_alpha(i) ** 2 / n_b1(i) * (-0.1D1 * n_b1(i - 1) + 0.2D1 * n
     #_b1(i) - 0.1D1 * n_b1(i + 1)) / hx ** 2 * n_beta(i) ** 2 - 0.75000
     #00000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx ** 2 / n_
     #alpha(i) ** 2 / n_a1(i) * n_beta(i) * (n_beta(i - 1) - 0.1D1 * n_b
     #eta(i + 1)) + 0.2500000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i
     # + 1)) / hx ** 2 / n_alpha(i) ** 3 / n_a1(i) * (n_alpha(i - 1) - 0
     #.1D1 * n_alpha(i + 1)) * n_beta(i) ** 2 + 0.2500000000000000D0 * (
     #n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx ** 2 / n_alpha(i) / 
     #n_b1(i) / n_a1(i) ** 2 * (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) + 0.2
     #500000000000000D0 * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx
     # ** 2 / n_b1(i) * n_beta(i) ** 2 / n_alpha(i) ** 3 * (n_b1(i - 1) 
     #- 0.1D1 * n_b1(i + 1)) - 0.5000000000000000D0 * (n_b1(i - 1) - 0.1
     #D1 * n_b1(i + 1)) / hx ** 2 / n_alpha(i) ** 2 / n_b1(i) * n_beta(i
     #) * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1))) / ctfmp(i) ** 2 + (0.
     #5000000000000000D0 * (nm1_alpha(i) - 0.1D1 * np1_alpha(i)) / ht / 
     #n_alpha(i) ** 3 * n_beta(i) - 0.5000000000000000D0 / n_alpha(i) **
     # 2 * (nm1_beta(i) - 0.1D1 * np1_beta(i)) / ht - 0.1D1 * (nm1_b1(i)
     # - 0.1D1 * np1_b1(i)) / ht / n_alpha(i) ** 2 / n_b1(i) * n_beta(i)
     # + 0.5000000000000000D0 * (-0.1D1 * nm1_a1(i) + np1_a1(i)) / ht / 
     #n_alpha(i) ** 2 / n_a1(i) * n_beta(i) + (-0.5000000000000000D0 * (
     #n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx * n_beta(i) ** 2 / n
     #_alpha(i) ** 3 - 0.1D1 * (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) / hx 
     #/ n_b1(i) / n_a1(i) ** 2 + 0.1D1 / n_alpha(i) ** 2 * n_beta(i) * (
     #n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx - 0.5000000000000000D0
     # * (n_alpha(i - 1) - 0.1D1 * n_alpha(i + 1)) / hx / n_alpha(i) / n
     #_a1(i) ** 2 + 0.5000000000000000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i
     # + 1)) / hx / n_a1(i) ** 3 + (n_b1(i - 1) - 0.1D1 * n_b1(i + 1)) /
     # hx / n_alpha(i) ** 2 / n_b1(i) * n_beta(i) ** 2 + 0.5000000000000
     #000D0 * (n_a1(i - 1) - 0.1D1 * n_a1(i + 1)) / hx / n_a1(i) * n_bet
     #a(i) ** 2 / n_alpha(i) ** 2) / ctfmp(i)) / ctfm(i)
      res(i)=qb
      end do
      END
