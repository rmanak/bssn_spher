      subroutine irev_dlamx(ctfm,ctfmp,n_A,n_Athth,n_Axx,n_B,n_K,n_Lamx,
     &n_Sx,n_alpha,n_beta,n_divbeta,n_phi,nm1_Lamx,np1_Lamx,octfmp,Nx,ht
     &,hx,vee,res)
      implicit none
      integer i
      integer Nx
      real*8 ht
      real*8 hx
      real*8 vee
      real*8 ctfm(Nx)
      real*8 ctfmp(Nx)
      real*8 n_A(Nx)
      real*8 n_Athth(Nx)
      real*8 n_Axx(Nx)
      real*8 n_B(Nx)
      real*8 n_K(Nx)
      real*8 n_Lamx(Nx)
      real*8 n_Sx(Nx)
      real*8 n_alpha(Nx)
      real*8 n_beta(Nx)
      real*8 n_divbeta(Nx)
      real*8 n_phi(Nx)
      real*8 nm1_Lamx(Nx)
      real*8 np1_Lamx(Nx)
      real*8 octfmp(Nx)
      real*8 res(Nx)
      real*8 qb
      do i=2, Nx-1, 1
      qb = 0.5000000000000000D0 * (-0.1D1 * nm1_Lamx(i) + np1_Lamx(i)) /
     # ht - 0.6666666666666667D0 * n_divbeta(i) * n_Lamx(i) * vee + 0.2D
     #1 / n_A(i) * n_alpha(i) * n_Sx(i) - 0.1D1 * (0.5000000000000000D0 
     #/ n_A(i) ** 3 * n_Axx(i) * n_alpha(i) * (-0.1D1 * n_A(i - 1) + n_A
     #(i + 1)) / hx + 0.1D1 / n_B(i) ** 2 / n_A(i) * (n_B(i - 1) - 0.1D1
     # * n_B(i + 1)) / hx * n_Athth(i) * n_alpha(i) - 0.6666666666666667
     #D0 / n_A(i) * n_alpha(i) * (-0.1D1 * n_K(i - 1) + n_K(i + 1)) / hx
     # - 0.5000000000000000D0 * n_beta(i) * (n_Lamx(i - 1) - 0.1D1 * n_L
     #amx(i + 1)) / hx - 0.1D1 / n_A(i) ** 2 * n_Axx(i) * (-0.1D1 * n_al
     #pha(i - 1) + n_alpha(i + 1)) / hx + 0.5000000000000000D0 * n_Lamx(
     #i) * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx + 0.500000000000
     #0000D0 * (n_beta(i - 1) - 0.1D1 * n_beta(i + 1)) / hx / n_A(i) * o
     #ctfmp(i) + 0.1666666666666667D0 / n_A(i) * vee * (-0.1D1 * n_divbe
     #ta(i - 1) + n_divbeta(i + 1)) / hx - 0.6D1 / n_A(i) ** 2 * n_alpha
     #(i) * n_Axx(i) * (n_phi(i - 1) - 0.1D1 * n_phi(i + 1)) / hx) / ctf
     #mp(i) + 0.1D1 / n_A(i) / ctfmp(i) ** 2 * (-0.1D1 * n_beta(i - 1) +
     # 0.2D1 * n_beta(i) - 0.1D1 * n_beta(i + 1)) / hx ** 2 - 0.1D1 * (-
     #0.4D1 / n_B(i) / n_A(i) * n_Athth(i) * n_alpha(i) + 0.4D1 / n_B(i)
     # ** 2 * n_alpha(i) * n_Athth(i) - 0.1D1 / n_B(i) * (n_beta(i - 1) 
     #- 0.1D1 * n_beta(i + 1)) / hx / ctfmp(i)) / ctfm(i) + 0.2D1 / ctfm
     #(i) ** 2 / n_B(i) * n_beta(i)
      res(i)=qb
      end do
      END
